#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
parse_metrics_slurm.py

Parses an output from the maester pipeline and returns a csv file with the quality metric values. 
This script is created in a relatively hardcoded way for UB375, but it can be generalized for other projects.

This is the modified version that can be run on a cluster with slurm.

This script has many hardcoded variables.

Author: UB INCLIVA (fmartinez@incliva.es)
Date: 05/09/2024
"""
import argparse
import glob
import os
import logging
import subprocess
import pandas as pd
import numpy as np

def count_matching_files(directory, pattern):
    """
    Counts the number of files matching a pattern in a directory.

    INPUT
    - directory (str): path to the directory.
    - pattern (str): pattern of the files you want to count.

    RETURN:
    - len(matching_files) (int): number of files matching the pattern.
    """
    matching_files = glob.glob(os.path.join(directory, pattern))
    print(os.path.join(directory, pattern))

    return len(matching_files)

def count_file_lines(file):
    """
    Counts the number of lines in a file.

    INPUT
    - file (str): path to the file.

    RETURN:
    - counter (int): number of lines in the file.
    """
    with open(file, 'r') as f:
        counter = sum(1 for _ in f)  # counts each line

    return counter

def get_barcodes_by_sample(cellbarcodes_dir, df, column, nombre_muestra):
    """
    Gets the number of cells per sample from the cellbarcodes files and puts these values in a dataframe.

    INPUT
    - cellbarcodes_dir (str): path to the folder with the cellbarcodes files. The filename is hardcoded.
    - df (pd.DataFrame): pandas dataframe that collects the quality metric values.
    - column (str): name of the dataframe column where the number of cellbarcodes will be saved.
    - nombre_muestra (str): sample name.

    RETURN:
    - df (pd.DataFrame): pandas dataframe that collects the quality metric values.
    """
    logging.info(f'Getting the number of high-quality cellbarcodes per sample from the files in "{cellbarcodes_dir}"...\n')

    pattern = f"*_{nombre_muestra}_cellBarcodes.csv" # hardcoded filename
    matching_files = glob.glob(os.path.join(cellbarcodes_dir, pattern)) # returns a list

    # Raise an exception if there is more than one matching file
    if len(matching_files) > 1:
        raise ValueError(f"Error: More than one matching cellbarcodes file found for '{nombre_muestra}'. Files: {matching_files}")

    sample_n_barcodes = count_file_lines(matching_files[0])

    df.loc[nombre_muestra, column] = sample_n_barcodes

    logging.info(f"\n{df}")
    logging.info(f"Dataframe dimensions: {df.shape}\n")

    return df

def count_reads_by_sample_fastq(path, filename_pattern, bioawk, df, column, nombre_muestra):
    """
    Counts the number of reads in the fastq.gz files of the samples.

    INPUT
    - path (str): path to the folder with the fastq.gz files.
    - filename_pattern (str): pattern of the fastq.gz files whose number of reads you want to count.
    - bioawk (str): command to run bioawk.
    - df (pd.DataFrame): pandas dataframe that collects the quality metric values.
    - column (str): name of the dataframe column where the number of reads will be saved.
    - nombre_muestra (str): sample name.

    RETURN:
    - df (pd.DataFrame): pandas dataframe that collects the quality metric values.
    """
    logging.info(f'Counting the number of reads in the files "{filename_pattern}" from "{path}"...\n')
    logging.info(f'The column {column} in the dataframe will be updated...\n')

    os.chdir(path)

    fastq_file = nombre_muestra + filename_pattern
    logging.info(f"{fastq_file}:")

    conteo = subprocess.run(f'{bioawk} {fastq_file} | wc -l', capture_output = True, text = True, shell = True)
    logging.info(f"{conteo.stdout.strip()}")

    df.loc[nombre_muestra, column] = int(conteo.stdout.strip())

    logging.info(f"\n{df}")
    logging.info(f"Dataframe dimensions: {df.shape}\n")

    return df

def count_reads_by_sample_bam(path, filename_pattern, samtools_view, df, column, nombre_muestra):
    """
    Counts the number of reads in the bam files of the samples.

    INPUT
    - path (str): path to the folder with the bam files.
    - filename_pattern (str): pattern of the bam files whose number of reads you want to count.
    - samtools_view (str): command to run samtools_view.
    - df (pd.DataFrame): pandas dataframe that collects the quality metric values.
    - column (str): name of the dataframe column where the number of reads will be saved.
    - nombre_muestra (str): sample name.

    RETURN:
    - df (pd.DataFrame): pandas dataframe that collects the quality metric values.
    """
    logging.info(f'Counting the number of reads in the files "{filename_pattern}" from "{path}"...\n')
    logging.info(f'The column {column} in the dataframe will be updated...\n')

    os.chdir(path)

    alignment_file = nombre_muestra + filename_pattern
    logging.info(f"{alignment_file}:")

    conteo = subprocess.run(f'{samtools_view} {alignment_file}', capture_output = True, text = True, shell = True)
    logging.info(f"{conteo.stdout.strip()}")

    df.loc[nombre_muestra, column] = int(conteo.stdout.strip())

    logging.info(f"\n{df}")
    logging.info(f"Dataframe dimensions: {df.shape}\n")

    return df

def count_cells_by_sample(path, filename_pattern, df, column, nombre_muestra):
    """
    Counts the number of files in the cell folders. The cells generated in 08_splitted_cells/ are counted.

    INPUT
    - path (str): path to the folder with the bam files.
    - filename_pattern (str): pattern of the bam files whose number of reads you want to count.
    - df (pd.DataFrame): pandas dataframe that collects the quality metric values.
    - column (str): name of the dataframe column where the number of reads will be saved.
    - nombre_muestra (str): sample name.

    RETURN:
    - df (pd.DataFrame): pandas dataframe that collects the quality metric values.
    """
    logging.info(f'Counting the number of files in the folders "{path}"...\n')
    logging.info(f'The column {column} in the dataframe will be updated...\n')

    os.chdir(path)

    os.chdir(f"{path}")
    archivos = os.listdir(f"{path}")

    archivos_bam = [archivo for archivo in archivos if archivo.endswith(filename_pattern)]

    cantidad_bam = len(archivos_bam)

    df.loc[nombre_muestra, column] = cantidad_bam

    logging.info(f"\n{df}")
    logging.info(f"Dataframe dimensions: {df.shape}\n")

    return df

def join_and_count_dedup_bams(path, metrics_path, filename_pattern, samtools_merge, samtools_view, df, column, nombre_muestra):
    """
    Merges the bams of the cells with deduplicated reads into a single one and counts its reads.

    INPUT
    - path (str): path to the folder with the bam files.
    - metrics_path (str): path to the folder where the metric files will be saved.
    - filename_pattern (str): pattern of the bam files of the deduplicated cells.
    - samtools_merge (str): command to run samtools_merge.
    - samtools_view (str): command to run samtools_view.
    - df (pd.DataFrame): pandas dataframe that collects the quality metric values.
    - column (str): name of the dataframe column where the number of reads will be saved.
    - nombre_muestra (str): sample name.

    RETURN:
    - df (pd.DataFrame): pandas dataframe that collects the quality metric values.
    """
    logging.info(f'Creating the bam file with the reads of all cells for sample {nombre_muestra}...\n')
    logging.info(f'The column {column} in the dataframe will be updated...\n')

    merged_bam = f"{metrics_path}/{nombre_muestra}.merged.bam"

    if not os.path.exists(merged_bam):
        dedup_preprocessed_files = glob.glob(f"{path}/{filename_pattern}") # these files are the bams of the already deduplicated cells

        bam_list_file = f"{metrics_path}/{nombre_muestra}_bam_list.txt" # temporary file to save the list of BAM files (too many to pass as arguments)

        with open(bam_list_file, 'w') as f: # write the names of the BAM files in the temporary file
            for bam_file in dedup_preprocessed_files:
                f.write(f"{bam_file}\n")

        create_merged_bam = subprocess.run(f'{samtools_merge} -b {bam_list_file} {merged_bam}', capture_output = True, text = True, shell = True) # merges the reads of all bams into one

        # Check the output
        if create_merged_bam.returncode == 0:
            logging.info(f"Merged BAM successfully created: {merged_bam}")

        else:
            logging.info(f"Error merging BAMs: {create_merged_bam.stderr}")

        # Clean up the temporary file after use
        if os.path.exists(bam_list_file):
            os.remove(bam_list_file)

    logging.info(f"Counting the (already deduplicated) reads of bam {merged_bam}...")

    conteo = subprocess.run(f'{samtools_view} {merged_bam}', capture_output = True, text = True, shell = True)
    logging.info(f"{conteo.stdout.strip()}")

    df.loc[nombre_muestra, column] = int(conteo.stdout.strip())

    return df

def get_median_sd_coverage_cut(coverages_table, df, nombre_muestra):
    """
    Calculates the median, standard deviation, and coverage cut-off from the coverage table.

    INPUT
    - coverages_table (str): path to the coverage table.

    RETURN:
    - df (pd.DataFrame): pandas dataframe that collects the quality metric values.
    """
    logging.info(f'Parsing the coverage table {coverages_table}...\n')
    logging.info(f'The columns "median_HQcb_coverage", "sd_HQcb_coverage" and "coverage_cut" will be updated in the dataframe...\n')
    df_c = pd.read_csv(coverages_table, sep = '\t')

    # Calculate the % of bases covered at more than 10x for each cell
    coverage_of_chrM = (df_c > 10).sum() * 100 / len(df_c)

    # Calculate the coverage cut-off based on the median minus the standard deviation
    chrM_coverage_cut = np.median(coverage_of_chrM) - np.std(coverage_of_chrM)

    df.loc[nombre_muestra, 'median_HQcb_coverage'] = int(np.median(coverage_of_chrM))
    df.loc[nombre_muestra, 'sd_HQcb_coverage'] = float(np.std(coverage_of_chrM))
    df.loc[nombre_muestra, 'coverage_cut'] = float(chrM_coverage_cut)

    return df

def main(arguments):
    sample = arguments.sample_name
    metrics_path = arguments.metrics_path

    output_file=  f"{metrics_path}/metrics_{sample}.csv"

    # TODO - update paths here
    base_dir = "/mnt/zonahpc/home/bioinformatica/servicios/20240724_UB375_MAESTER_Asherman"
    rawdata_dir = f"{base_dir}/rawdata"
    analysis_dir = f"{base_dir}/analysis"
    cellbarcodes_dir = f"{base_dir}/nuevos_datasets"
    assembled_fastq_dir = f"{analysis_dir}/01_assembled_fastqs"
    fastp_dir = f"{analysis_dir}/02_fastp"
    trimmed_fastq_dir = f"{analysis_dir}/03_trimmed_fastqs"
    mod_fastq_dir = f"{analysis_dir}/04_mod_fastqs"
    mapeo_dir = f"{analysis_dir}/05_mapeo"
    mitochondiral_bam_dir = f"{analysis_dir}/07_mitochondrial_bams"
    splitted_cells_dir = f"{analysis_dir}/08_splitted_cells"
    preprocessed_bam_dir = f"{analysis_dir}/10_preprocessed_bams"
    coverages_dir = f"{analysis_dir}/11_coverages"
    variants_dir= f'{analysis_dir}/13_tablas_variantes'
    singularity_dir = "/mnt/zonahpc/home/bioinformatica/comun/software/singularity"

    # software used to count reads
    bioawk_path = f"{singularity_dir}/bioawk:1.0--he4a0461_10"
    bioawk = f"singularity exec -B {base_dir},{analysis_dir} {bioawk_path} bioawk -c fastx '{{print $name}}'" # command to run bioawk

    samtools_path = f"{singularity_dir}/samtools:1.19.2--h50ea8bc_1"
    samtools_view = f"singularity exec -B {base_dir},{analysis_dir} {samtools_path} samtools view -c -F 0x900" # excludes reads with SUPPLEMENTARY (0x100) or SECONDARY (0x800) flag, meaning it only counts primary reads and not secondary or supplementary.
    samtools_merge = f"singularity exec -B {base_dir},{analysis_dir} {samtools_path} samtools merge" # command to run samtools merge

    bcftools_path = f"{singularity_dir}/bcftools:1.18--h8b25389_0"
    bcftools_count = f"singularity exec -B {base_dir},{analysis_dir} {bcftools_path} bcftools view -H" # command to run bcftools view -H -c
    # Set up logging
    logging.basicConfig(
        format = '%(asctime)s - %(message)s',
        level = logging.INFO,
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    # pandas dataframe that collects the quality metric values
    columnas = ["muestra", "n_cellbarcodes", "splitted_cells", "rawdata", "assembled_fastq", "fastp1", "fastp2", "fastp3", "trimmed", "mod_fastq", "alineamiento_sam", "maester_bam", "dedup_bam", "median_HQcb_coverage", "sd_HQcb_coverage", "coverage_cut", "good_coverage_HQcb", "bad_coverage_HQcb", "n_variants"]
    df = pd.DataFrame(columns = columnas)
    df['muestra'] = sample
    df.set_index("muestra", inplace = True)

    # n_cellbarcodes
    df = get_barcodes_by_sample(cellbarcodes_dir, df, "n_cellbarcodes", sample)
    # splitted_cells
    df = count_cells_by_sample(f'{splitted_cells_dir}/{sample}', '.bam', df, 'splitted_cells', sample)

    # Count reads in fastq.gz files:
    # rawdata
    df = count_reads_by_sample_fastq(rawdata_dir, "_R1.fastq.gz", bioawk, df, "rawdata", sample)
    # assembled_fastq
    df = count_reads_by_sample_fastq(assembled_fastq_dir, ".fastq.gz", bioawk, df, "assembled_fastq", sample)
    # fastp1
    df = count_reads_by_sample_fastq(fastp_dir, ".fastp1.fastq.gz", bioawk, df, "fastp1", sample)
    # fastp2
    df = count_reads_by_sample_fastq(fastp_dir, ".fastp2.fastq.gz", bioawk, df, "fastp2", sample)
    # fastp3
    df = count_reads_by_sample_fastq(fastp_dir, ".fastp3.fastq.gz", bioawk, df, "fastp3", sample)
    # trimmed
    df = count_reads_by_sample_fastq(trimmed_fastq_dir, ".fastp3.fastq.gz.trimmed", bioawk, df, "trimmed", sample)
    # mod_fastq
    df = count_reads_by_sample_fastq(mod_fastq_dir, ".trimmed.mod.fastq.gz", bioawk, df, "mod_fastq", sample)

    # Count reads in sam/bam files:
    # alineamiento_sam
    df = count_reads_by_sample_bam(mapeo_dir, "_Aligned.out.sam", samtools_view, df, "alineamiento_sam", sample)
    # mod_bam
    df = count_reads_by_sample_bam(mitochondiral_bam_dir, "_maester.bam", samtools_view, df, "maester_bam", sample)

    # Write a first csv
    df.to_csv(output_file)

    # dedup_bam: merge all bams of the splitted cells into one and count the reads it has.
    df = join_and_count_dedup_bams(f'{preprocessed_bam_dir}/{sample}', metrics_path, "*.preprocessed.bam", samtools_merge, samtools_view, df, "dedup_bam", sample)

    # median_HQcb_coverage: parse and process the table 11_coverages/${sample}/coverages.tsv created by inc_calculate_coverage_cut.R
    # sd_HQcb_coverage: parse and process the table 11_coverages/${sample}/coverages.tsv created by inc_calculate_coverage_cut.R
    # coverage_cut: parse and process the table 11_coverages/${sample}/chrM_coverage_cut.txt created by inc_calculate_coverage_cut.R
    df = get_median_sd_coverage_cut(f'{coverages_dir}/{sample}/coverages.tsv', df, sample)

    # good_coverage_HQcb: count files in 11_coverages/${sample}/good_coverage_cells
    # bad_coverage_HQcb: count cells (lines) in file bad_coverage_cells.txt
    df = count_cells_by_sample(f'{coverages_dir}/{sample}/good_coverage_cells/', ".preprocessed.bam", df, "good_coverage_HQcb", sample)

    logging.info(f'Calculating cells with low coverage for sample {sample}...\n')

    all_HQcb = int(count_matching_files(f'{preprocessed_bam_dir}/{sample}', '*.preprocessed.bam'))
    good_coverage_HQcb = int(count_matching_files(f'{coverages_dir}/{sample}/good_coverage_cells', '*.preprocessed.bam'))

    df.loc[sample, "good_coverage_HQcb"] = good_coverage_HQcb
    df.loc[sample, "bad_coverage_HQcb"] = all_HQcb - good_coverage_HQcb

    # n_variants: number of variants in the table 13_tabla_variantes/${sample}/tabnla_AF_variantes.csv
    logging.info(f'Counting variants in {variants_dir}/{sample}/merged.vcf.gz...\n')
    conteo = subprocess.run(f'{bcftools_count} {variants_dir}/{sample}/merged.vcf.gz | wc -l', capture_output = True, text = True, shell = True)
    df.loc[sample, "n_variants"] = conteo.stdout.strip()

    # Write the final csv
    df.to_csv(output_file)
    logging.info(df)

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--sample_name', action = 'store', type = str, dest = 'sample_name', required = True, help = 'Name of the sample to be analyzed.')
parser.add_argument('-m', '--metrics_path', action = 'store', type = str, dest = 'metrics_path', required = True, help = 'Full path where the generated metrics will be saved.')

arguments = parser.parse_args()

if __name__ == "__main__":
    main(arguments)
