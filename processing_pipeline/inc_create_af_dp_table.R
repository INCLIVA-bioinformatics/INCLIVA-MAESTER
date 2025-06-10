#!/usr/bin/Rscript
# -*- coding: utf-8 -*-

# inc_create_af_dp_table.R

# Starting from a VCF, creates tables with samples as rows and variants as columns.
# Generates 2 tables: one with AF (allele frequency) and another with DP (read depth).
# It is designed for LoFreq VCFs, which should have the INFO field in the following format:
# DP=117;AF=0.264957;SB=0;DP4=83,0,31,0
# Additionally, the sample name must be present in the GT column.

# Author: UB INCLIVA (jfcatala@incliva.es, fmartinez@incliva.es)
# Date: 27/02/2024
#-------------------------------------------------------------------------------
library(tidyverse)
library(vcfR)
library(argparse)

# Define arguments
parser <- ArgumentParser()
parser$add_argument("--folder", help = "Path to the folder", type = "character")
parser$add_argument("--output_dir", help = "Path to the output directory", type = "character")
args <- parser$parse_args()

if (is.null(args$folder)) {
  stop("You must provide the folder path using --folder")
}

if (is.null(args$output)) {
  stop("You must provide the output path using --output_dir")
}

print(paste("Folder path:", args$folder))
print(paste("Output file path:", args$output))

extract <- function(vcf_file){
  vcf <- read.vcfR(vcf_file, verbose = FALSE)
  df <- as_tibble(vcf@fix) %>% separate(col = "INFO", into = c("DP", "AF", "SB", "DP4"), remove = FALSE, sep = ";")
  df <- df %>% mutate(
    DP = as.numeric(str_remove(DP, "DP=")),
    AF = as.numeric(str_remove(AF, "AF=")),
    MUT = str_c(POS, REF, ">", ALT),
    cell = colnames(vcf@gt)[2] # sample name
  ) %>% dplyr::select(MUT, cell, DP, AF)
  return(df)
}

path <- args$folder
out_dir <- args$output

# List all files in a vector
files <- list.files(path = path, pattern = "*.mod.vcf.gz$", full.names = TRUE, recursive = FALSE)
lista2 <- lapply(files, extract)
lista3 <- bind_rows(lista2)

# Remove _norm from sample names
df_variants <- lista3 %>% mutate(cell = str_replace_all(cell, "_norm", ""))

tabla_AF <- df_variants %>%
  dplyr::select(MUT, cell, AF) %>%
  pivot_wider(names_from = "MUT", values_from = "AF")
tabla_AF[is.na(tabla_AF)] <- 0

tabla_DP <- df_variants %>%
  dplyr::select(MUT, cell, DP) %>%
  pivot_wider(names_from = "MUT", values_from = "DP")
tabla_DP[is.na(tabla_DP)] <- 0

write_csv(tabla_AF, file.path(out_dir, "tabla_AF_variants.csv"))
write_csv(tabla_DP, file.path(out_dir, "tabla_DP_variants.csv"))
