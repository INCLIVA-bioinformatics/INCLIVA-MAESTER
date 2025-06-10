# INCLIVA - MAESTER

This pipeline is based on the methodology described in the paper by Miller et al. (2022), Mitochondrial variant enrichment from high-throughput single-cell RNA sequencing resolves clonal populations (*Nat Biotechnol* **40**, 1030–1034). The original code and methods have been optimized to improve performance and efficiency. For more details, see: https://doi.org/10.1038/s41587-022-01210-8

Repository strucuture is as following:

```
maester_pipeline/
├── analyze_variant
│   └── explore_variants_final.r # used to explore and get informative variants obtained through `MAESTER_slurm.sh`. Also, creates several plots.
├── get_HQ_barcodes.R # obtains HQcb (High Quality cell-barcode) from scRNAseq object.
├── parse_metrics_slurm.py # calculates several quality statistics form a `MAESTER_slurm.sh` execution.
├── processing_pipeline
│   ├── inc_assemble_fastq.R # used in `MAESTER_slurm.sh` to assemble the fastq files for high-quality cells from raw data
│   ├── inc_calculate_coverage_cut.R # used in `MAESTER_slurm.sh` to calculate "good quality coverage" threshold in our HQcb.
│   ├── inc_create_af_dp_table.R # used in `MAESTER_slurm.sh` to creates allele frequency and depth table for all variants in all cells.
│   ├── inc_parse_vcf_lofreq.py # used in `MAESTER_slurm.sh`. Corrects `lofreq` vcf format.
│   ├── inc_split_bam.py # used in `MAESTER_slurm.sh`. Splits MAESTER bam in one bam per HQcb.
│   ├── MAESTER_slurm.conf # used along `MAESTER_slurm.sh`. Contains paths, software and parameters.
│   └── MAESTER_slurm.sh # main script for processing and analyizing maester data.
└── README.md
```
