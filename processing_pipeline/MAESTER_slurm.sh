#!/bin/bash
# -*- coding: utf-8 -*-

# MAESTER_slurm.sh

# scRNAseq data analysis pipeline (MAESTER) adapted for execution on an HPC cluster using SLURM. Based on the original MAESTER pipeline, but modified to work with SLURM and to include additional steps such as deduplication and variant calling:
# Miller, T.E., Lareau, C.A., Verga, J.A. et al. Mitochondrial variant enrichment from high-throughput single-cell RNA sequencing resolves clonal populations. Nat Biotechnol 40, 1030–1034 (2022). https://doi.org/10.1038/s41587-022-01210-8

# PREPROCESSING (PART 1) - A
## A.0 - FastQC
## A.1 - scRNASeq -> to be done in R
## A.2 - MAESTER
### A.2.1 - Add UMI + barcodes to R2 reads
### A.2.2 - Quality filtering + poly-A trimming with fastp
### A.2.3 - Remove primers
### A.2.4 - Parse read names and compress file
### A.2.4B - Run FastQC again
### A.2.5 - Map the reads
### A.2.6 - Add cell barcode and UMI tags to the BAM
### A.2.7 - Select mitochondrial chromosome reads
# PREPROCESSING (PART 2) - B
## B.1 - Split BAM by barcode and index
## B.2 - Deduplicate
### B.2.fgbio - Deduplicate using fgbio and index
#### B.2.fgbio.1 - Mapped BAM -> Grouped BAM
#### B.2.fgbio.2 - GroupedBam -> Consensus uBam
#### B.2.fgbio.3 - Consensus uBam -> Consensus
#### B.2.fgbio.4 - Consensus Mapped -> Consensus Filtered & Sorted BAM (filter and sort the consensus reads)
### B.2.umitools - Deduplicate using umitools and index
## B.3 - Modify RG tags
## B.4 - Adjust mapping quality if deduplicated with umitools, since some callers interpret STAR's mapping quality as unmapped
# VARIANT CALLING - C
## C.1 - Select only those cells with good coverage
## C.2 - Run recalibration + variant calling with lofreq + variant normalization (only for cells with good coverage)
## C.3 - Merge all VCFs

# Author: UB INCLIVA (fmartinez@incliva.es, smzuniga@incliva.es), based on the previously mentioned MAESTER pipeline by T.E. Miller et al.
# Date: 27/08/2024
########################################################################################################################

usage() {
    echo "Usage: $(basename $0) -n sample_name -c config_file"
    echo "Description: MAESTER pipeline prepared for use with SLURM."
    echo
    echo "Options:"
    echo "  -n string    Text string corresponding to the sample name (e.g. 109_entire_E1)."
    echo "  -c string    Path to the configuration file."
    echo
    exit 1
}

while getopts "n:c:" opt; do
    case ${opt} in
        n) sample=${OPTARG} ;; # sample name passed as an argument (to be used in SLURM)
        c) conf_file=${OPTARG} ;; # configuration file passed as an argument
        :) echo "Option -$OPTARG requires an argument." >&2; exit 1 ;;
        *) usage ;;
    esac
done

if [ -z "${sample}" ] || [ -z "${conf_file}" ];
then
    echo "ERROR: -n and -c are mandatory parameters"
    usage
    exit
fi

# VARIABLES
set -x
source ${conf_file} # load configuration file
set +x # Keep "set -x" during development and debugging

echo -e "\nStarting analysis of sample ${sample}."
echo -e "\nConfiguration file: ${conf_file}"
echo -e "\nAnalysis directory: ${analysis_directory}"

# SCRIPT

# PREPROCESSING (PART 1) - A ############################################################################################

# A.0 - FastQC
mkdir -p ${analysis_directory}/00_raw_QC && cd ${analysis_directory}/00_raw_QC

echo -e "\nCOMMAND - ${fastqc} ${rawdata_directory}/${sample}* --outdir ${analysis_directory}/00_raw_QC"
${fastqc} ${rawdata_directory}/${sample}* --outdir ${analysis_directory}/00_raw_QC

## A.1 - scRNASeq -> previously done in R

## A.2 - MAESTER

### A.2.1 - Add UMI + barcodes to R2 reads

# - R1: 28 nucleotides (16 nucleotides for the cell barcode followed by 12 for the UMI)
# - R2: 256 nucleotides

mkdir -p ${analysis_directory}/01_assembled_fastqs && cd ${analysis_directory}/01_assembled_fastqs

echo -e "\nCOMMAND - ${assemble_fastq} --muestra ${sample} --cell_barcodes ${cell_barcodes} --rawdata ${rawdata_directory} --output_dir ${analysis_directory}/01_assembled_fastqs --CBlength ${cb_legth} --UMIlength ${umi_length}"
${assemble_fastq} --muestra ${sample} --cell_barcodes ${cell_barcodes} --rawdata ${rawdata_directory} --output_dir ${analysis_directory}/01_assembled_fastqs --CBlength ${cb_legth} --UMIlength ${umi_length}

### A.2.2 - Quality filtering + poly-A removal with fastp

# NOTE: Before performing this step, it is recommended to visualize the *fastq* files in IGV. Originally, the 3' end looked very "dirty", which is why this quality filtering is applied with the following parameters.

# In 3 steps, we will trim 50 nts from the 3' end, then remove poly-X with a minimum of 3 nts, and also apply quality filtering. Finally, to remove poly-A, we will treat Ax15 as an adapter.
#   - trim x reads from the 3' end             -> -t 50
#   - trim poly-As with minimum length of 3    -> --trim_poly_x --poly_x_min_len 3
#   - disable length filtering                 -> -L
#   - disable adapter trimming                 -> -A
#   - disable poly-G trimming                  -> -G
#   - disable quality filtering                -> -Q

# We perform fastp in 3 steps.
mkdir -p ${analysis_directory}/02_fastp && cd ${analysis_directory}/02_fastp

# Step 1. Trim 50 nts at the 3' end
echo -e "\nCOMMAND - ${fastp} -t ${trim_3_nts} -L -A -G -Q -i ${analysis_directory}/01_assembled_fastqs/${sample}.fastq.gz -o ${analysis_directory}/02_fastp/${sample}.fastp1.fastq.gz --html fastp1.html --json fastp1.json"
${fastp} -t ${trim_3_nts} -L -A -G -Q -i ${analysis_directory}/01_assembled_fastqs/${sample}.fastq.gz -o ${analysis_directory}/02_fastp/${sample}.fastp1.fastq.gz --html fastp1.html --json fastp1.json

# Step 2. Remove poly-As and filter by quality
echo -e "\nCOMMAND - ${fastp} --trim_poly_x --poly_x_min_len 3 -L -A -G -i ${analysis_directory}/02_fastp/${sample}.fastp1.fastq.gz -o ${analysis_directory}/02_fastp/${sample}.fastp2.fastq.gz --html fastp2.html --json fastp2.json"
${fastp} --trim_poly_x --poly_x_min_len 3 -L -A -G -i ${analysis_directory}/02_fastp/${sample}.fastp1.fastq.gz -o ${analysis_directory}/02_fastp/${sample}.fastp2.fastq.gz --html fastp2.html --json fastp2.json

# Step 3. Remove Ax15 “adapters” (in reality, remove remaining poly-As not trimmed in the previous step). Also remove poly-Gs as suggested by fastqc.
echo -e "\nCOMMAND - ${fastp} -L -Q -a AAAAAAAAAAAAAAA -i ${analysis_directory}/02_fastp/${sample}.fastp2.fastq.gz -o ${analysis_directory}/02_fastp/${sample}.fastp3.fastq.gz --html fastp3.html --json fastp3.json"
${fastp} -L -Q -a AAAAAAAAAAAAAAA -i ${analysis_directory}/02_fastp/${sample}.fastp2.fastq.gz -o ${analysis_directory}/02_fastp/${sample}.fastp3.fastq.gz --html fastp3.html --json fastp3.json

### A.2.3 - Remove primers

# Remove the first 24 bases from the 5' end of the reads corresponding to primers.

# Trim 24bp from 5'
mkdir -p ${analysis_directory}/03_trimmed_fastqs && cd ${analysis_directory}/03_trimmed_fastqs
ln -s ${analysis_directory}/02_fastp/${sample}.fastp3.fastq.gz ${analysis_directory}/03_trimmed_fastqs

echo -e "\nCOMMAND - ${homerTools} trim -5 ${trim_5_nts} ${analysis_directory}/03_trimmed_fastqs/${sample}.fastp3.fastq.gz"
${homerTools} trim -5 ${trim_5_nts} ${analysis_directory}/03_trimmed_fastqs/${sample}.fastp3.fastq.gz

### A.2.4 (I) - Parse read names and compress file

# The fastq file must have this format in the first line of each read: @SRR15598774_1_GCCTGTTTCCGAACGC_TCATTCTGGCTT
# The _1_ refers to the read number and is unique, so you’ll have as many different numbers as reads. Then comes the cell barcode and the UMI.

# In my case, the previous step’s output had read names like:
# @SRR15598774.1 1 length=264_GCCTGTTTCCGAACGC_TCATTCTGGCTT, and I modified them to: @SRR15598774_1_GCCTGTTTCCGAACGC_TCATTCTGGCTT so that the next scripts work.

# Either run this step to compress the file and parse the name:
mkdir -p ${analysis_directory}/04_mod_fastqs && cd ${analysis_directory}/04_mod_fastqs

echo -e "\nCOMMAND - sed \"s/@SRR15598774\\./@SRR15598774_/\" ${analysis_directory}/03_trimmed_fastqs/${sample}.fastp3.fastq.gz.trimmed | sed \"s/ .*length=2..//\" | ${bgzip} > ${analysis_directory}/04_mod_fastqs/${sample}.trimmed.mod.fastq.gz"
sed "s/@SRR15598774\./@SRR15598774_/" ${analysis_directory}/03_trimmed_fastqs/${sample}.fastp3.fastq.gz.trimmed | sed "s/ .*length=2..//" | ${bgzip} > ${analysis_directory}/04_mod_fastqs/${sample}.trimmed.mod.fastq.gz

# Or just compress the file if the reads already have the correct naming:

# path_to_modified_trimmed_fastq <- file.path(output_dir, "MAESTER", paste0(SampleName, ".trimmed.mod.fastq.gz"))
# system(paste(bgzip, "-f", trimmed_file))
# system(paste("mv", paste0(trimmed_file,".gz"), path_to_modified_trimmed_fastq))

# A.2.4 (II) - FastQC (after fastp and trimming)
mkdir -p ${analysis_directory}/04B_processed_QC && cd ${analysis_directory}/04B_processed_QC

echo -e "\nCOMMAND - ${fastqc} ${analysis_directory}/04_mod_fastqs/${sample}.trimmed.mod.fastq.gz --outdir ${analysis_directory}/04B_processed_QC"
${fastqc} ${analysis_directory}/04_mod_fastqs/${sample}.trimmed.mod.fastq.gz --outdir ${analysis_directory}/04B_processed_QC

### A.2.5 - Align reads

# It's important to ensure that the mitochondrial chromosome in the STAR genome reference is named **chrM** and not chrMT. If it's chrMT, replace it later in the BAM.
mkdir -p ${analysis_directory}/05_mapeo && cd ${analysis_directory}/05_mapeo

echo -e "\nCOMMAND - ${STAR} --runThreadN ${cpus} --genomeDir ${path_to_genome} --readFilesIn ${analysis_directory}/04_mod_fastqs/${sample}.trimmed.mod.fastq.gz --readFilesCommand zcat --outFileNamePrefix ${analysis_directory}/05_mapeo/${sample}_"
${STAR} --runThreadN ${cpus} --genomeDir ${path_to_genome} --readFilesIn ${analysis_directory}/04_mod_fastqs/${sample}.trimmed.mod.fastq.gz --readFilesCommand zcat --outFileNamePrefix ${analysis_directory}/05_mapeo/${sample}_

### A.2.6 - Add cell barcode and UMI tags to BAM

# This is done by the Tag_CB_UMI.sh Bash script from MAESTER-2021 Pre-processing directory. Samtools must be in PATH. Script is in /20230621_UB153_scRNASeq_cell_lineage/software/MAESTER-2021/Pre-processing/Tag_CB_UMI.sh. Below is the equivalent code.

# What it does:
# Move cell barcode (CB) and unique molecular identifier (UMI) from read identifier to SAM tags.
# Resulting BAM will contain tags for CB and UMI following 10X convention: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/bam
# These tags are required for maegatk.
mkdir -p ${analysis_directory}/06_modified_bams && cd ${analysis_directory}/06_modified_bams

echo -e "\nCOMMAND - ${samtools} view -h ${analysis_directory}/05_mapeo/${sample}_Aligned.out.sam | \
awk 'BEGIN{FS=\"\t\"; OFS=\"\t\"}{if (substr(\$1,1,1) == \"@\"){print \$0} else {split(\$1, a, \"_\"); \$1=\"\" ; print a[1]\"_\"a[2]\$0\"\tCB:Z:\"a[3]\"-1\tUB:Z:\"a[4] } }' | \
${samtools} view -bh > ${analysis_directory}/06_modified_bams/${sample}_alignment.out.modified.bam"
${samtools} view -h ${analysis_directory}/05_mapeo/${sample}_Aligned.out.sam | \
awk 'BEGIN{FS="\t"; OFS="\t"}{if (substr($1,1,1) == "@"){print $0} else {split($1, a, "_"); $1="" ; print a[1]"_"a[2]$0"\tCB:Z:"a[3]"-1\tUB:Z:"a[4] }}' | \
${samtools} view -bh > ${analysis_directory}/06_modified_bams/${sample}_alignment.out.modified.bam

### A.2.7 - Select mitochondrial reads
mkdir -p ${analysis_directory}/07_mitochondrial_bams && cd ${analysis_directory}/07_mitochondrial_bams

echo -e "\nCOMMAND - ${samtools} sort -@ ${cpus} ${analysis_directory}/06_modified_bams/${sample}_alignment.out.modified.bam > ${analysis_directory}/07_mitochondrial_bams/${sample}_alignment.out.modified.sorted.bam"
${samtools} sort -@ ${cpus} ${analysis_directory}/06_modified_bams/${sample}_alignment.out.modified.bam > ${analysis_directory}/07_mitochondrial_bams/${sample}_alignment.out.modified.sorted.bam
echo -e "\nCOMMAND - ${samtools} index -@ ${cpus} ${analysis_directory}/07_mitochondrial_bams/${sample}_alignment.out.modified.sorted.bam"
${samtools} index -@ ${cpus} ${analysis_directory}/07_mitochondrial_bams/${sample}_alignment.out.modified.sorted.bam
echo -e "\nCOMMAND - ${samtools} view -b -@ ${cpus} ${analysis_directory}/07_mitochondrial_bams/${sample}_alignment.out.modified.sorted.bam chrM > ${analysis_directory}/07_mitochondrial_bams/${sample}_maester.bam"
${samtools} view -b -@ ${cpus} ${analysis_directory}/07_mitochondrial_bams/${sample}_alignment.out.modified.sorted.bam chrM > ${analysis_directory}/07_mitochondrial_bams/${sample}_maester.bam
echo -e "\nCOMMAND - ${samtools} index -@ ${cpus} ${analysis_directory}/07_mitochondrial_bams/${sample}_maester.bam"
${samtools} index -@ ${cpus} ${analysis_directory}/07_mitochondrial_bams/${sample}_maester.bam
echo -e "\nCOMMAND - ${samtools} sort  -@ ${cpus} -t CB ${analysis_directory}/07_mitochondrial_bams/${sample}_maester.bam -o ${analysis_directory}/07_mitochondrial_bams/${sample}_maester_sorted.bam"
${samtools} sort  -@ ${cpus} -t CB ${analysis_directory}/07_mitochondrial_bams/${sample}_maester.bam -o ${analysis_directory}/07_mitochondrial_bams/${sample}_maester_sorted.bam # sort BAM by CB tag

# PREPROCESSING (PART 2) - B ###########################################################################################

mkdir -p ${analysis_directory}/08_splitted_cells/${sample} && cd ${analysis_directory}/08_splitted_cells/${sample}

## B.1 - Split BAM by barcode and index
echo -e "\nCOMMAND - ${split_bam} -i ${analysis_directory}/07_mitochondrial_bams/${sample}_maester_sorted.bam -c ${cell_barcodes}"
${split_bam} -i ${analysis_directory}/07_mitochondrial_bams/${sample}_maester_sorted.bam -c ${cell_barcodes} # warning: sorted BAM is not indexed (problem because sorting is by CB tag, not by position), but it seems to work anyway.

echo -e "\nINFO - Indexing individual BAMs..."
for file in ${analysis_directory}/08_splitted_cells/${sample}/*.bam; do
    echo -e "\nCOMMAND - ${samtools} index -@ ${cpus} ${file}"
    ${samtools} index -@ ${cpus} ${file}
done

mkdir -p ${analysis_directory}/09_dedup_cells/${sample} && cd ${analysis_directory}/09_dedup_cells/${sample}

## B.2 - Deduplication
echo -e "\nINFO - Starting deduplication..."
if [ "${dedup_method}" == "fgbio" ]
then
    echo -e "\nINFO - Deduplicating with fgbio..."
    ### B.2.fgbio - Deduplicate using fgbio and index
    for file in ${analysis_directory}/08_splitted_cells/${sample}/*.bam; do
        cell_name=$(basename "${file}" | cut -d'.' -f1) # base filename without extension
        echo -e "\nINFO - Deduplicating cell ${cell_name}..."

        #### B.2.fgbio.1 - Mapped BAM -> Grouped BAM
        echo -e "\nCOMMAND - ${fgbio} --compression 1 --async-io GroupReadsByUmi --input ${file} --strategy Adjacency -t UB --output ${cell_name}.grouped.bam --family-size-histogram ${cell_name}.tag-family-sizes.txt"
        ${fgbio} --compression 1 --async-io GroupReadsByUmi --input ${file} --strategy Adjacency -t UB --output ${cell_name}.grouped.bam --family-size-histogram ${cell_name}.tag-family-sizes.txt

        #### B.2.fgbio.2 - GroupedBam -> Consensus uBam
        echo -e "\nCOMMAND - ${fgbio} --compression 1 CallMolecularConsensusReads --min-reads 1 --input ${cell_name}.grouped.bam --threads ${cpus} --output ${cell_name}.cons.unmapped.bam"
        ${fgbio} --compression 1 CallMolecularConsensusReads --min-reads 1 --input ${cell_name}.grouped.bam --threads ${cpus} --output ${cell_name}.cons.unmapped.bam

        #### B.2.fgbio.3 - Consensus uBam -> Consensus Mapped BAM (re-align the consensus reads)
        echo -e "\nCOMMAND - ${samtools} fastq ${cell_name}.cons.unmapped.bam | ${bwa_mem2} -t ${cpus} -p -K 150000000 -Y ${genome} - | ${fgbio} --compression 1 --async-io ZipperBams --unmapped ${cell_name}.cons.unmapped.bam --ref ${genome} --tags-to-reverse Consensus --tags-to-revcomp Consensus --output ${cell_name}.cons.mapped.bam"
        ${samtools} fastq ${cell_name}.cons.unmapped.bam | ${bwa_mem2} -t ${cpus} -p -K 150000000 -Y ${genome} - | ${fgbio} --compression 1 --async-io ZipperBams --unmapped ${cell_name}.cons.unmapped.bam --ref ${genome} --tags-to-reverse Consensus --tags-to-revcomp Consensus --output ${cell_name}.cons.mapped.bam

        echo -e "\nCOMMAND - ${samtools} sort --threads ${cpus} -o ${cell_name}.cons.mapped.sorted.bam ${cell_name}.cons.mapped.bam"
        ${samtools} sort --threads ${cpus} -o ${cell_name}.cons.mapped.sorted.bam ${cell_name}.cons.mapped.bam
        echo -e "\nCOMMAND - ${samtools} index -@ ${cpus} ${cell_name}.cons.mapped.sorted.bam"
        ${samtools} index -@ ${cpus} ${cell_name}.cons.mapped.sorted.bam

        #### B.2.fgbio.4 - Consensus Mapped -> Consensus Filtered & Sorted BAM (filter and sort the consensus reads)
        echo -e "\nCOMMAND - ${fgbio} --compression 0 FilterConsensusReads --input ${cell_name}.cons.mapped.bam --output ${cell_name}.cons.filtered.bam --ref ${genome} --min-reads ${dedup_min_reads} --min-base-quality ${dedup_min_bq} --max-base-error-rate ${dedup_max_ber} --max-read-error-rate ${dedup_max_rer} --max-no-call-fraction ${dedup_max_no_cf}"
        ${fgbio} --compression 0 FilterConsensusReads --input ${cell_name}.cons.mapped.bam --output ${cell_name}.cons.filtered.bam --ref ${genome} --min-reads ${dedup_min_reads} --min-base-quality ${dedup_min_bq} --max-base-error-rate ${dedup_max_ber} --max-read-error-rate ${dedup_max_rer} --max-no-call-fraction ${dedup_max_no_cf}

        echo -e "\nCOMMAND - ${samtools} sort --threads ${cpus} -o ${cell_name}.cons.filtered.sorted.bam ${cell_name}.cons.filtered.bam"
        ${samtools} sort --threads ${cpus} -o ${cell_name}.cons.filtered.sorted.bam ${cell_name}.cons.filtered.bam

        # Rename to match the rest of the pipeline
        echo -e "\nCOMMAND - mv ${cell_name}.cons.filtered.sorted.bam ${cell_name}.deduplicated_raw.bam"
        mv ${cell_name}.cons.filtered.sorted.bam ${cell_name}.deduplicated_raw.bam
    done

    # Update BAM headers; this method works with fgbio (using only picard caused errors). Prepares for RG change required by fgbio.
    echo -e "\nINFO - Updating BAM headers..."
    for file in ${analysis_directory}/09_dedup_cells/${sample}/*deduplicated_raw.bam; do
        dedup_name=$(basename "${file}" .deduplicated_raw.bam)
        echo -e "\nINFO - Updating header of ${dedup_name}..."
        # First update @RG to include all fields, including SM (ID must remain as A, as reads have A).
        echo -e "\nCOMMAND - ${samtools} view -H ${dedup_name}.deduplicated_raw.bam | sed 's,^@RG.*,@RG\tID:A\tPU:'${dedup_name}'\tSM:'${dedup_name}'\tLB:lib1\tPL:ILLUMINA,g' | ${samtools} reheader - ${dedup_name}.deduplicated_raw.bam > ${dedup_name}.deduplicated.bam"
        ${samtools} view -H ${dedup_name}.deduplicated_raw.bam | sed 's,^@RG.*,@RG\tID:A\tPU:'${dedup_name}'\tSM:'${dedup_name}'\tLB:lib1\tPL:ILLUMINA,g' | ${samtools} reheader - ${dedup_name}.deduplicated_raw.bam > ${dedup_name}.deduplicated.bam
    done

elif [ "${dedup_method}" == "umitools" ]
then
    ### B.2.umitools - Deduplicate using umitools and index
    echo -e "\nINFO - Deduplicating with umitools..."
    for file in ${analysis_directory}/08_splitted_cells/${sample}/*.bam; do
        bam_name=$(basename "${file}" | cut -d'.' -f1) # base filename without extension
        echo -e "\nINFO - Deduplicating cell ${bam_name}..."

        echo -e "\nCOMMAND - ${umi_tools} dedup -I ${file} --output-stats=${bam_name}_umitools_stats.txt -S ${bam_name}.deduplicated.bam --umi-tag UB --extract-umi-method=tag -L ${bam_name}.log"
        ${umi_tools} dedup -I ${file} --output-stats=${bam_name}_umitools_stats.txt -S ${bam_name}.deduplicated.bam --umi-tag UB --extract-umi-method=tag -L ${bam_name}.log

        echo -e "\nCOMMAND - ${samtools} index -@ ${threads} ${bam_name}.deduplicated.bam"
        ${samtools} index -@ ${threads} ${bam_name}.deduplicated.bam

        # Additionally, extract reads per UMI information
        info=`grep "Mean number of unique UMIs per position" ${bam_name}.log | sed 's/.*per position: //'`

        echo ${bam_name} " " ${info} >> ${analysis_directory}/reads_per_umi.txt
    done

    else
        echo -e "\nERROR - Unrecognized deduplication method. Must choose 'fgbio' or 'umitools'.\n"
        exit 1
fi

## B.3 - Change RG
echo -e "\nINFO - Changing RG..."
for file in ${analysis_directory}/09_dedup_cells/${sample}/*deduplicated.bam; do
    dedup_name=$(basename "${file}" .deduplicated.bam)  # Extracts the name without the .deduplicated.bam extension
    echo -e "\nINFO - Changing RG of ${dedup_name}..."

    echo -e "\nCOMMAND - ${picard} AddOrReplaceReadGroups -I ${file} -O ${analysis_directory}/09_dedup_cells/${sample}/${dedup_name}.dedup.RG.bam -RGLB lib1 -RGPL ILLUMINA -RGPU ${dedup_name} -RGSM ${dedup_name}"
    ${picard} AddOrReplaceReadGroups -I "${file}" -O "${analysis_directory}/09_dedup_cells/${sample}/${dedup_name}.dedup.RG.bam" -RGLB lib1 -RGPL ILLUMINA -RGPU "${dedup_name}" -RGSM "${dedup_name}"
done

mkdir -p ${analysis_directory}/10_preprocessed_bams/${sample} && cd ${analysis_directory}/10_preprocessed_bams/${sample}

## B.4 - Change mapping quality if we deduplicated with umitools; the mapping quality set by STAR is sometimes interpreted by some callers as unmapped
if [ "${dedup_method}" == "umitools" ]
then
    echo -e "\nINFO - Changing mapping quality since we deduplicated with umitools..."
    for file in ${analysis_directory}/09_dedup_cells/${sample}/*.dedup.RG.bam; do
        RG_name=$(basename "${file}" .deduplicated.bam)
        echo -e "\nINFO - Changing mapping quality of ${RG_name}..."

        echo -e "\nCOMMAND - ${samtools} view -H ${file} > header;  ${samtools} view ${file} | awk '{if ($5 == 255) print $0}' | awk '$5=60' | sed 's/ /\t/gi' > body; cat header body | ${samtools} view -bSh > ${analysis_directory}/10_preprocessed_bams/${sample}/${RG_name}.preprocessed.bam"
        ${samtools} view -H ${file} > header;  ${samtools} view ${file} | awk '{if ($5 == 255) print $0}' | awk '$5=60' | sed 's/ /\t/gi' > body; cat header body | ${samtools} view -bSh > ${analysis_directory}/10_preprocessed_bams/${sample}/${RG_name}.preprocessed.bam

        rm header body

        echo -e "\nCOMMAND - ${samtools} index -@ ${cpus} ${analysis_directory}/10_preprocessed_bams/${sample}/${RG_name}.preprocessed.bam"
        ${samtools} index -@ ${cpus} ${analysis_directory}/10_preprocessed_bams/${sample}/${RG_name}.preprocessed.bam
    done

elif [ "${dedup_method}" == "fgbio" ] # If we deduplicated with fgbio, changing mapping quality is not necessary, since the reads were remapped with bwa-mem2.
then
    echo -e "\nINFO - No need to change mapping quality since we deduplicated with fgbio..."
    for file in ${analysis_directory}/09_dedup_cells/${sample}/*.dedup.RG.bam; do
        RG_name=$(basename "${file}" .dedup.RG.bam)
        ln -s ${file} ${analysis_directory}/10_preprocessed_bams/${sample}/${RG_name}.preprocessed.bam
        echo -e "\nINFO - Indexing ${RG_name}..."

        echo -e "\nCOMMAND - ${samtools} index -@ ${cpus} ${analysis_directory}/10_preprocessed_bams/${sample}/${RG_name}.preprocessed.bam"
        ${samtools} index -@ ${cpus} ${analysis_directory}/10_preprocessed_bams/${sample}/${RG_name}.preprocessed.bam
    done

else
    echo -e "\nERROR - Unrecognized deduplication method. Must choose 'fgbio' or 'umitools'.\n"
    exit 1
fi

# VARIANT CALLING - C ##################################################################################################
## C.1 - Select only those cells with good coverage
mkdir -p ${analysis_directory}/11_coverages/${sample} && cd ${analysis_directory}/11_coverages/${sample}

echo -e "\nINFO - Calculating coverages..."
for file in ${analysis_directory}/10_preprocessed_bams/${sample}/*.preprocessed.bam; do
    cell_bc=$(basename "${file}" .preprocessed.bam)
    depth_file="${cell_bc}_depth.tsv"

    # Check if the depth file already exists, since it takes quite some time...
    if [ ! -f "${depth_file}" ]; then
        echo -e "COMMAND -- ${samtools} depth -@ ${cpus} -aa ${file} -b ${intervals} -d 0 > ${depth_file}"
        ${samtools} depth -@ ${cpus} -aa ${file} -b ${intervals} -d 0 > "${depth_file}"
    else
        echo "File ${depth_file} already exists, skipping..."
    fi
done

# Once coverages of all cells are calculated, we need to set the coverage cutoff to consider a cell as having good coverage.
# The cutoff is proposed as median minus one standard deviation.
echo -e "\nCOMMAND - ${calculate_cov_cut} --folder ${analysis_directory}/11_coverages/${sample} --output ${analysis_directory}/11_coverages/${sample}/coverages.tsv --cutoff ${analysis_directory}/11_coverages/${sample}/chrM_coverage_cut.txt"
${calculate_cov_cut} --folder ${analysis_directory}/11_coverages/${sample} --output ${analysis_directory}/11_coverages/${sample}/coverages.tsv --cutoff ${analysis_directory}/11_coverages/${sample}/chrM_coverage_cut.txt # This step creates a tsv with coverage per cell and position and also saves the cutoff value (median - 1 std dev) in ${analysis_directory}/coverages/chrM_coverage_cut.txt

corte_cobertura=$(cat ${analysis_directory}/11_coverages/${sample}/chrM_coverage_cut.txt)
rounded_coverage_cut=$(printf "%.0f" ${corte_cobertura})

echo -e "\nINFO - The coverage cutoff for sample ${sample} is ${rounded_coverage_cut}x."

mkdir -p ${analysis_directory}/11_coverages/${sample}/good_coverage_cells && cd ${analysis_directory}/11_coverages/${sample}/good_coverage_cells

echo -e "\nINFO - Selecting cells with good coverage (greater or equal to the established cutoff)..."

# Initialize the empty list-like string before the loop
failed_cell_bcs=""

for file in ${analysis_directory}/10_preprocessed_bams/${sample}/*.preprocessed.bam; do
    cell_bc=$(basename "${file}" .preprocessed.bam)
    echo -e "\nINFO - Checking coverage of ${cell_bc}..."

    # Here we have an issue because we don't have the capture bed. We extract the % of chrM that is covered

    total_positions=$(wc -l < ${analysis_directory}/11_coverages/${sample}/${cell_bc}_depth.tsv) # The mitochondrial genome has 16569 positions

    covered_positions=$(awk '$3 >= '${cov}' {count ++} END {print count}' "${analysis_directory}/11_coverages/${sample}/${cell_bc}_depth.tsv") # Positions covered above the 'cov' threshold

    coverage_percentage=$(awk -v cpo="${covered_positions}" -v tpo="${total_positions}" 'BEGIN {print (cpo / tpo) * 100}') # Percentage of positions covered above 'cov' relative to total
    rounded_coverage=$(awk -v val="${coverage_percentage}" 'BEGIN { printf("%.0f", val) }') # Needed for the 'if' condition below.

    # Extract bed of covered regions
    echo -e "\nCOMMAND - awk '{if ($3 >= '${cov}') print $1\"\t\"$2-1\"\t\"$2}' ${analysis_directory}/11_coverages/${sample}/${cell_bc}_depth.tsv | ${bedtools} merge > ${analysis_directory}/11_coverages/${sample}/${cell_bc}_covered_regions_${cov}x.bed"
    awk '{if ($3 >= '${cov}') print $1"\t"$2-1"\t"$2}' ${analysis_directory}/11_coverages/${sample}/${cell_bc}_depth.tsv | ${bedtools} merge > ${analysis_directory}/11_coverages/${sample}/${cell_bc}_covered_regions_${cov}x.bed

    # If coverage is greater or equal to the cutoff we decided (remember cutoff was median - 1*std dev), declare the cell as having good coverage and keep working with it
    if [ ${rounded_coverage} -ge ${rounded_coverage_cut} ]; then
        echo -e "${cell_bc} has good coverage: (${rounded_coverage}%)."

        ln -sf ${analysis_directory}/10_preprocessed_bams/${sample}/${cell_bc}.preprocessed.bam .
        ln -sf ${analysis_directory}/10_preprocessed_bams/${sample}/${cell_bc}.preprocessed.bam.bai .
    else
        failed_cell_bcs+="${cell_bc}," # Add the cell to the list of cells with poor coverage
        echo ${failed_cell_bcs} > ${analysis_directory}/11_coverages/${sample}/bad_coverage_cells.txt
    fi
done

# Print the list-like string of $cell_bc that failed the evaluation, removing the last comma
if [ -n "$failed_cell_bcs" ]; then
    # Remove last comma if it exists
    failed_cell_bcs=${failed_cell_bcs%,}
    echo "The following cells do not have good coverage: ${failed_cell_bcs}"
else
    echo "All cells meet the minimum coverage."
fi

mkdir -p ${analysis_directory}/12_variant_calling/${sample} && cd ${analysis_directory}/12_variant_calling/${sample}

## C.2 - Run recalibration + variant calling with lofreq + variant normalization (only for cells with good coverage)
echo -e "\nINFO - Running recalibration + variant calling with lofreq + variant normalization..."
for file in  ${analysis_directory}/11_coverages/${sample}/good_coverage_cells/*.preprocessed.bam; do
    cell_bc=$(basename ${file} | cut -d'.' -f1) # base filename without extension
    echo -e "\nINFO - Processing cell ${cell_bc}..."

    # Calculate recalibration table
    echo -e "\nCOMMAND - ${gatk} --java-options ${java_options} BaseRecalibrator -I ${file} -R ${genome} --intervals ${intervals} --known-sites ${GATK_indels} --known-sites ${dbSNP} --output ${cell_bc}.table"
    ${gatk} --java-options ${java_options} BaseRecalibrator -I ${file} -R ${genome} --intervals ${intervals} --known-sites ${GATK_indels} --known-sites ${dbSNP} --output ${cell_bc}.table

    # Apply recalibration
    echo -e "\nCOMMAND - ${gatk} --java-options ${java_options} ApplyBQSR -R ${genome} -I ${file} --bqsr-recal-file ${cell_bc}.table -O ${cell_bc}.recalibrated.bam"
    ${gatk} --java-options ${java_options} ApplyBQSR -R ${genome} -I ${file} --bqsr-recal-file ${cell_bc}.table -O ${cell_bc}.recalibrated.bam

    # Indel quality with lofreq
    echo -e "\nCOMMAND - ${lofreq} indelqual --dindel -f ${genome} -o ${cell_bc}.recalibrated.indelqual.bam ${cell_bc}.recalibrated.bam"
    ${lofreq} indelqual --dindel -f ${genome} -o ${cell_bc}.recalibrated.indelqual.bam ${cell_bc}.recalibrated.bam

    # Index the recalibrated indelqual file
    echo -e "\nCOMMAND - ${samtools} index -@ ${cpus} ${cell_bc}.recalibrated.indelqual.bam"
    ${samtools} index -@ ${cpus} ${cell_bc}.recalibrated.indelqual.bam

    # Variant calling with lofreq
    echo -e "\nCOMMAND - ${lofreq} call -f ${genome} -l ${intervals} --force-overwrite --no-default-filter --call-indels ${cell_bc}.recalibrated.indelqual.bam ${lofreq_relaxed} -o ${cell_bc}.vcf"
    ${lofreq} call -f ${genome} -l ${intervals} --force-overwrite --no-default-filter --call-indels ${cell_bc}.recalibrated.indelqual.bam ${lofreq_relaxed} -o ${cell_bc}.vcf

    echo -e "\nCOMMAND - ${bgzip} ${cell_bc}.vcf && ${tabix} ${cell_bc}.vcf.gz"
    ${bgzip} ${cell_bc}.vcf && ${tabix} ${cell_bc}.vcf.gz

    # Normalize variants
    echo -e "\nCOMMAND - ${bcftools} norm -m-any ${cell_bc}.vcf.gz | ${bcftools} norm -Oz --check-ref -w -f ${genome} > ${cell_bc}_norm.vcf.gz"
    ${bcftools} norm -m-any ${cell_bc}.vcf.gz | ${bcftools} norm -Oz --check-ref -w -f ${genome} > ${cell_bc}_norm.vcf.gz

    echo -e "\nCOMMAND - ${tabix} ${cell_bc}_norm.vcf.gz"
    ${tabix} ${cell_bc}_norm.vcf.gz

    # Process lofreq vcf since it does not genotype by default
    echo -e "\nCOMMAND - ${python} ${parse_vcf_lofreq} -i ${cell_bc}_norm.vcf.gz -o ${cell_bc}.mod.vcf -he ${het_vcf_lofreq} -ho ${hom_vcf_lofreq}"
    ${python} ${parse_vcf_lofreq} -i ${cell_bc}_norm.vcf.gz -o ${cell_bc}.mod.vcf -he ${het_vcf_lofreq} -ho ${hom_vcf_lofreq}

    echo -e "\nCOMMAND - ${bgzip} ${cell_bc}.mod.vcf && ${tabix} ${cell_bc}.mod.vcf.gz"
    ${bgzip} ${cell_bc}.mod.vcf && ${tabix} ${cell_bc}.mod.vcf.gz
done

## C.3 - Merge all vcfs
echo -e "\nINFO - Merging all vcfs..."

mkdir -p ${analysis_directory}/13_tablas_variantes/${sample} && cd ${analysis_directory}/13_tablas_variantes/${sample}

echo -e "\nCOMMAND - ulimit -n ${max_ficheros_abiertos}"
ulimit -n ${max_ficheros_abiertos} # increase open file limit, otherwise bcftools merge fails

# Intermediate merge: for each letter, merge files starting with that letter
echo -e "\nINFO - Merging intermediate files letter by letter..."

# Define the letters to classify the files
letters=("A" "C" "G" "T")

for letter in "${letters[@]}"; do
    files_to_merge=$(find ${analysis_directory}/12_variant_calling/${sample} -maxdepth 1 -type f -name "${letter}*.mod.vcf.gz" | tr '\n' ' ')

    if [ -n "${files_to_merge}" ]; then
        echo -e "\nCOMMAND - ${bcftools} merge -m none -O z -o ${analysis_directory}/13_tablas_variantes/${sample}/merged_${letter}.vcf.gz ${files_to_merge}"
        ${bcftools} merge -m none -O z -o ${analysis_directory}/13_tablas_variantes/${sample}/merged_${letter}.vcf.gz ${files_to_merge}

        echo -e "\nCOMMAND - ${tabix} ${analysis_directory}/13_tablas_variantes/${sample}/merged_${letter}.vcf.gz"
        ${tabix} ${analysis_directory}/13_tablas_variantes/${sample}/merged_${letter}.vcf.gz
    else
        echo "No VCF files found starting with ${letter}."
    fi
done

# Final merge: combine the 4 files generated in the previous step
files_to_merge_final=$(find ${analysis_directory}/13_tablas_variantes/${sample} -maxdepth 1 -type f -name "merged_*.vcf.gz" | tr '\n' ' ')

if [ -n "${files_to_merge_final}" ]; then
    echo -e "\nCOMMAND - ${bcftools} merge -m none -O z -o ${analysis_directory}/13_tablas_variantes/${sample}/merged.vcf.gz ${files_to_merge_final}"
    ${bcftools} merge -m none -O z -o ${analysis_directory}/13_tablas_variantes/${sample}/merged.vcf.gz ${files_to_merge_final}

    echo -e "\nCOMMAND - ${tabix} ${analysis_directory}/13_tablas_variantes/${sample}/merged.vcf.gz"
    ${tabix} ${analysis_directory}/13_tablas_variantes/${sample}/merged.vcf.gz
else
    echo "No intermediate merged VCF files found."
fi

## C.4 - Create the table of variants contained, by cell, with DP and AF
echo -e "\nCOMMAND - ${create_af_dp_table} --folder ${analysis_directory}/12_variant_calling/${sample} --output_dir ${analysis_directory}/13_tablas_variantes/${sample}"
${create_af_dp_table} --folder ${analysis_directory}/12_variant_calling/${sample} --output_dir ${analysis_directory}/13_tablas_variantes/${sample}

echo -e "\nFinished!\n"
