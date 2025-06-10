#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
inc_split_bam.py

Takes as input a BAM file with reads tagged with cell barcodes and a file listing the cell barcodes of interest.
It then creates a separate BAM file for each barcode of interest containing the corresponding reads.

Author: UB INCLIVA (fmartinez@incliva.es)
Date: 22/02/2024
"""
import argparse
import pysam

def create_individual_bam(bam_file, cell_barcodes):
    """
    Takes a BAM file with reads tagged with cell barcodes and a file listing barcodes of interest.
    Creates one BAM file per barcode containing its corresponding reads.

    INPUT
    - bam_file (str): Path to the BAM file containing reads tagged with cell barcodes.
    - cell_barcodes (set): Set of barcodes of interest (e.g., from a CSV file).
    """
    current_cell = None
    current_bam = None

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            cell_barcode = read.get_tag("CB") if read.has_tag("CB") else None

            # Check if we encountered a new cell barcode or an unmapped read
            if cell_barcode != current_cell or read.is_unmapped:
                if current_bam:
                    current_bam.close()
                current_cell = cell_barcode
                if current_cell in cell_barcodes:
                    current_bam = pysam.AlignmentFile(f"{current_cell}.bam", "wb", template=bam)
                else:
                    current_bam = None
            if current_bam:
                current_bam.write(read)

    if current_bam:
        current_bam.close()

def main(arguments):
    bam_file = arguments.bam_file  # e.g. /path/to/MAESTER_sorted.bam
    cell_barcodes_file = arguments.cell_barcodes_file  # e.g. /path/to/HQ_CBs.csv

    with open(cell_barcodes_file) as f:
        cell_barcodes = {line.strip() for line in f}

    create_individual_bam(bam_file, cell_barcodes)

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_bam', action='store', type=str, dest='bam_file', required=True,
                    help='Full path to the BAM file containing reads tagged with cell barcodes.')
parser.add_argument('-c', '--cell_barcodes', action='store', type=str, dest='cell_barcodes_file', required=True,
                    help='Full path to the CSV file listing single-cell barcodes of interest.')

arguments = parser.parse_args()

if __name__ == "__main__":
    main(arguments)
