# INCLIVA - MAESTER

This pipeline is partially based on the methodology described in the paper by Miller et al. (2022), Mitochondrial variant enrichment from high-throughput single-cell RNA sequencing resolves clonal populations (*Nat Biotechnol* **40**, 1030–1034). The original code and methods have been optimized to improve performance and efficiency. For more details, see: https://doi.org/10.1038/s41587-022-01210-8

Repository strucuture is as following:

```
maester_pipeline/
├── analyze_variant
│   └── explore_variants_final.r
├── get_HQ_barcodes.R
├── parse_metrics_slurm.py 
├── processing_pipeline
│   ├── inc_assemble_fastq.R
│   ├── inc_calculate_coverage_cut.R
│   ├── inc_create_af_dp_table.R
│   ├── inc_parse_vcf_lofreq.py
│   ├── inc_split_bam.py
│   ├── MAESTER_slurm.conf
│   └── MAESTER_slurm.sh
└── README.md
```

- `analyze_variant/`
  - `explore_variants_final.r`: Used to explore and extract informative variants obtained through `MAESTER_slurm.sh`. Also generates several plots.

- `get_HQ_barcodes.R`: Obtains HQcb (High Quality cell-barcodes) from a scRNA-seq object.

- `parse_metrics_slurm.py`: Calculates several quality statistics from a `MAESTER_slurm.sh` execution.

- `processing_pipeline/`
  - `inc_assemble_fastq.R`: Assembles FASTQ files for high-quality cells from raw data (used in `MAESTER_slurm.sh`).
  - `inc_calculate_coverage_cut.R`: Calculates the "good quality coverage" threshold in HQcb (used in `MAESTER_slurm.sh`).
  - `inc_create_af_dp_table.R`: Creates an allele frequency and depth table for all variants in all cells (used in `MAESTER_slurm.sh`).
  - `inc_parse_vcf_lofreq.py`: Corrects `lofreq` VCF format (used in `MAESTER_slurm.sh`).
  - `inc_split_bam.py`: Splits the MAESTER BAM file into one BAM per HQcb (used in `MAESTER_slurm.sh`).
  - `MAESTER_slurm.conf`: Configuration file with paths, software, and parameters (used with `MAESTER_slurm.sh`).
  - `MAESTER_slurm.sh`: Main script for processing and analyzing MAESTER data.

- `README.md`
