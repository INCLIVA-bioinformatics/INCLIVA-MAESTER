#!/usr/bin/Rscript
# -*- coding: utf-8 -*-

# inc_assemble_fastq.R

# This script assembles the fastq files for high-quality cells from raw data.

# Author: UB INCLIVA (fmartinez@incliva.es)
# Date: 02/09/2024
#-------------------------------------------------------------------------------

# PACKAGES
# .libPaths("/home/fmartinez/R/x86_64-pc-linux-gnu-library/4.2") # user library path (comment out if using Singularity container)
library(argparse)
library(glue)
library(ShortRead)
library(Seurat)
library(MitoTrace)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(Rsamtools)
library(Matrix)
library(seqinr)

# PARSE ARGUMENTS
parser <- ArgumentParser(description = 'Script to analyze MAESTER in single cell data.')
parser$add_argument('--muestra', required = TRUE, help = 'Sample name')
parser$add_argument('--cell_barcodes', required = TRUE, help = 'Full path to CSV file containing high-quality cell barcodes')
parser$add_argument('--rawdata', required = TRUE, help = 'Full path to directory with raw data (R1 & R2)')
parser$add_argument('--output_dir', required = TRUE, help = 'Full path to output directory for assembled fastqs')
parser$add_argument('--CBlength', required = TRUE, help = 'Length of cell barcodes')
parser$add_argument('--UMIlength', required = TRUE, help = 'Length of UMIs')

# Parse arguments
args <- parser$parse_args()

Folder     <- args$rawdata
SampleName <- args$muestra
HQ_CBs     <- args$cell_barcodes
output_dir <- args$output_dir
CBlength   <- as.numeric(args$CBlength)
UMIlength  <- as.numeric(args$UMIlength)

# Print all variables
cat("Sample name:", SampleName, "\n")
cat("Raw data folder:", Folder, "\n")
cat("High-quality barcodes:", HQ_CBs, "\n")
cat("Output directory:", output_dir, "\n")
cat("Barcode length:", CBlength, "\n")
cat("UMI length:", UMIlength, "\n")

# FUNCTIONS
# Function to assemble FASTQ files
cutf <- function(x, f=1, d="/", ...) sapply(strsplit(x, d), function(i) i[f], ...)
assemble_fastq <- function(Folder, SampleName, CellBarcodes, CBlength, UMIlength, output_dir) {
    # Set options to handle large output and avoid converting strings to factors
    options(max.print = 500)
    options(stringsAsFactors = FALSE)
    options(scipen = 999)  # Disable scientific notation

    # Identify R1 FASTQ files (R2 files will be identified later by substitution)
    R1.ch <- list.files(Folder, pattern = paste0(SampleName, ".*_R1.fastq.gz$"), full.names = TRUE)

    # Display files to be processed
    message(Sys.time(), "\nLoading ", length(R1.ch) * 2, " FASTQ files:")
    message(cat(c(R1.ch, sub("_R1", "_R2", R1.ch)), sep = "\n"))

    # Stop if no FASTQ files were found
    if (length(R1.ch) == 0) stop("Did not find FASTQ files.")

    # Read and process cell barcodes
    cells.split <- unlist(strsplit(CellBarcodes, ","))
    cells.df <- do.call(rbind, lapply(cells.split, function(x) read.table(x)))
    cells.ch <- cutf(cells.df$V1, d = "-", f = 1)

    if (length(cells.ch) == 0) stop("No cells found.")

    message("Found ", length(cells.ch), " cells.\n")

    # Initialize list to store read counts
    report.ls <- list()  # List to store number of reads

    # Process each R1 FASTQ file
    for (f1 in R1.ch) {
        # Identify corresponding R2 file
        f2 <- sub("_R1", "_R2", f1)

        # Load files in chunks of 10 million reads
        message("file ", match(f1, R1.ch), "/", length(R1.ch), ": ", basename(f1), " ", appendLF = FALSE)
        strm1 <- FastqStreamer(f1, n = 1E7)  # 10 million reads by default
        strm2 <- FastqStreamer(f2, n = 1E7)

        # Process reads in chunks
        repeat {
            message("*", appendLF = FALSE)
            fq1 <- yield(strm1)
            fq2 <- yield(strm2)
            if (length(fq1) == 0 || length(fq2) == 0) break

            # Match reads to expected barcodes
            fq1.m <- is.element(as.vector(subseq(sread(fq1), 1, CBlength)), cells.ch)
            report.ls[[basename(f1)]] <- c(length(fq1.m), sum(fq1.m))

            # Filter unmatched reads
            fq1.f <- fq1[fq1.m]
            fq2.f <- fq2[fq1.m]

            # Extract barcode and UMI from Read1
            fq1.f.cell <- as.vector(subseq(sread(fq1.f), 1, CBlength))
            fq1.f.umi <- as.vector(subseq(sread(fq1.f), CBlength + 1, CBlength + UMIlength))

            # Update Read2 IDs with barcode and UMI
            fq2.f@id <- BStringSet(paste0(sub(" .:N:0:", "_", as.vector(fq2.f@id)), "_", fq1.f.cell, "_", fq1.f.umi))

            # Check if all IDs match
            if (!all(cutf(as.vector(fq1.f@id), d = " ", f = 1) == cutf(as.vector(fq2.f@id), d = "_", f = 1))) {
                stop("Read ID mismatch")
            }

            # Save filtered R2 reads to output file
            writeFastq(fq2.f, file = file.path(output_dir, paste0(SampleName, ".fastq.gz")), mode = "a")

            # Clean up memory
            invisible(gc())
        }

        close(strm1)
        close(strm2)
        message(" done")
        invisible(gc())
    }

    # Generate and save report
    report.mat <- do.call(rbind, report.ls)
    report.mat <- rbind(report.mat, colSums(report.mat))
    rownames(report.mat)[nrow(report.mat)] <- "total"
    report.mat <- cbind(report.mat, report.mat[, 2] / report.mat[, 1])
    colnames(report.mat) <- c("all", "filtered", "fraction")
    write.table(report.mat, file = file.path(output_dir, paste0(SampleName, ".stats.txt")), sep = "\t", quote = FALSE)

    invisible(gc())

    message("\nMaintained ", round(report.mat["total", "fraction"] * 100, 2), "% of reads.\n")
    sessionInfo()
    message("\nFinished!")
}

# SCRIPT
system(paste("mkdir -p", output_dir))

cat("Running 'assemble_fastq' function...", "\n\n")
assemble_fastq(Folder, SampleName, HQ_CBs, CBlength, UMIlength, output_dir)
