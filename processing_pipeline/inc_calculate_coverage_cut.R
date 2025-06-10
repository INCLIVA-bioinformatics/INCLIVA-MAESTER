#!/usr/bin/Rscript
# -*- coding: utf-8 -*-

# inc_calculate_coverage_cut.R

# Calculates the coverage cut-off point for declaring a cell to have good coverage. The comments are from
# UB281. Minimum good coverage is defined as the median coverage minus one SD.

# Author: UB INCLIVA (fmartinez@incliva.es)
# Date: 03/04/2024
#-------------------------------------------------------------------------------
library(argparse)
library(Cairo)
library(matrixStats)
library(tidyverse)

create_coverage_per_pos_table <- function(cov_path) {
  # Check if the file is empty
  if (file.info(cov_path)$size == 0) {
    stop("El archivo está vacío: ", cov_path)
  }

  sample <- sub(".*/([^/]+)_depth.tsv", "\\1", cov_path) # Extract the name of the sample
  header <- c("chr", "pos", paste0(sample, "_depth")) # Header is created since the tsv does not have one.
  cov_table <- tryCatch(
    {
      read.table(cov_path, header = FALSE, sep = "\t", check.names = FALSE, col.names = header)
    },
    error = function(e) {
      stop("Error reading file: ", e$message)
    }
  )

  # Ensure that the pos (position) column is numeric
  cov_table$pos <- as.numeric(cov_table$pos)

  # Ensure that the depth value is numeric
  cov_table[, 3] <- as.numeric(cov_table[, 3])

  return(cov_table)
}

# Define a function that merges dataframes by common columns (common_columns is HARDCODEADOWN)
fusionar_dataframes <- function(df1, df2) {
  columnas_comunes <- c("chr", "pos") # common columns in all df --> HARDCODEED
  merge(df1, df2, by = columnas_comunes, all = TRUE)
}

# Define the arguments
parser <- ArgumentParser()
parser$add_argument("--folder", help = "Path to the folder with the depth files", type = "character")
parser$add_argument("--output", help = "Path of the table file with the coverages to be created", type = "character")
parser$add_argument("--cutoff", help = "Path of the file to save the coverage breakpoint", type = "character")
args <- parser$parse_args()

if (is.null(args$folder)) {
  stop("You must provide the path to the folder using --folder")
}

if (is.null(args$output)) {
  stop("You must provide the path to the file to be created using --output")
}

if (is.null(args$cutoff)) {
  stop("You must provide the file path for the cut point using --cutoff")
}

coverages_path <- args$folder
ruta_archivo <- args$output
cutoff <- args$cutoff

print(paste("Path to folder:", coverages_path))
print(paste("Path of the file to create:", ruta_archivo))
print(paste("Path of the cut point file:", cutoff))

if (file.exists(ruta_archivo)) {
  print("The table with the coverages already exists. Reading the file...")

  df_c <- read.table(ruta_archivo, header = T, sep = "\t", check.names = F)

} else {
  print("The table with the coverages does NOT exist. We proceed to create it...")

  depth_files <- list.files(path = coverages_path, pattern = "*_depth.tsv$", full.names = T, recursive = F)
  lista2 <- lapply(depth_files, create_coverage_per_pos_table)

  dataframe_final <- Reduce(fusionar_dataframes, lista2) # It takes a long time (30-40 min, depending on the dataset).
  df_ordenado <- dataframe_final[order(dataframe_final$pos), ]

  rownames(df_ordenado) <- df_ordenado$pos
  df_c <- df_ordenado %>% select(-chr, -pos)

  write.table(df_c, file = ruta_archivo, sep = "\t", quote = FALSE, row.names = FALSE)
}

print(paste("There are", as.character(nrow(df_c)), "positions in the mitochondrial chromosome."))
print(paste("There are", as.character(ncol(df_c)), "cells in this sample."))

output_file <- file.path(coverages_path, "coverages_histogram.png")
CairoPNG(filename = output_file, width = 800, height = 600)
hist_to_plot <- colSums(df_c == 0) * 100 / nrow(df_c) # Positions that are NOT covered
hist(hist_to_plot, breaks = 1000, main = "Histogram of Uncovered Positions", xlab = "Percentage of Uncovered Positions", ylab = "Frequency")
dev.off()

print(paste("The histogram has been saved in:", output_file))

coverage_of_chrM <- data.frame(colSums(df_c > 10) * 100 / nrow(df_c)) # Draw, for each cell, the % of bases that are covered at more than 10x.
colnames(coverage_of_chrM) <- "chrM_coverage_(%)"

summary(coverage_of_chrM)
chrM_coverage_cut <- as.numeric(median(coverage_of_chrM$`chrM_coverage_(%)`) - sd(coverage_of_chrM$`chrM_coverage_(%)`))

writeLines(as.character(chrM_coverage_cut), cutoff)

print(paste("Coverage cutoff point stored in:", cutoff))

