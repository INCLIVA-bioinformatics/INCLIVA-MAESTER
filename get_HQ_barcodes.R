#!/usr/bin/Rscript
# -*- coding: utf-8 -*-

# inc_assemble_fastq.R

# This script processes single-cell RNA-seq datasets stored in .h5ad format,
# converting them to Seurat's .h5seurat format for analysis.
# It loads two datasets, performs cell counting by patient, biopsy site, and treatment stage,
# and exports cell barcodes to CSV files grouped by patient and experimental conditions.
# Parallelization is enabled to speed up data loading and processing.

# Author: UB INCLIVA (fmartinez@incliva.es, smzuniga@incliva.es)
# Date: 02/09/2024
#-------------------------------------------------------------------------------
library(hdf5r)
library(SeuratDisk)
library(future)

atlas <- LoadH5Seurat("/media/scratch0/20240212_UB281_MAESTER_piloto/nuevos_datasets/GSE215968_sc_AS_pre_vs_post.h5seurat")
# Calculate cell counts per patient, biopsy site (stroma or epithelium), and treatment stage (pre or post).
table(atlas$biospy_orig, atlas$Patient, atlas$Treatment_stage)

file_path <- "/media/scratch0/20240212_UB281_MAESTER_piloto/nuevos_datasets/GSE215968_sc_CD133.h5ad"
Convert(file_path, dest = "h5seurat", overwrite = T)
plan(multicore, workers = 18)
atlas <- LoadH5Seurat("/media/scratch0/20240212_UB281_MAESTER_piloto/nuevos_datasets/GSE215968_sc_AS_pre_vs_post.h5seurat")
cd133 <- LoadH5Seurat("/media/scratch0/20240212_UB281_MAESTER_piloto/nuevos_datasets/GSE215968_sc_CD133.h5seurat")

for(patient in unique(cd133$Patient)){
  one_patient_cells <- subset(cd133, subset = Patient == patient)
  one_patient_cells <- sub(".*CD133_","", colnames(one_patient_cells))
  write.table(one_patient_cells, paste0("/media/scratch0/20240212_UB281_MAESTER_piloto/nuevos_datasets/GSE215968_sc_CD133_",
                                        patient, "_cellBarcodes.csv"), quote = F, row.names = F,
                                        col.names = F)
}

for(patient in unique(atlas$Patient)){
  for (site in unique(atlas$biospy_orig)){
    for (treatment in unique(atlas$Treatment_stage)){
      one_patient_cells <- subset(atlas, subset = Patient == patient)
      one_patient_cells <- subset(one_patient_cells, subset = biospy_orig == site)
      if (nrow(one_patient_cells) != 0 ){
        one_patient_cells <- subset(one_patient_cells, subset = Treatment_stage == treatment)
        if(nrow(one_patient_cells) != 0){
          one_patient_cells <- sub(".*_","", colnames(one_patient_cells))
          write.table(one_patient_cells, paste0("/media/scratch0/20240212_UB281_MAESTER_piloto/nuevos_datasets/GSE215968_sc_",
                                                patient,"_", site,"_",treatment,"_cellBarcodes.csv"), quote = F, row.names = F, col.names = F)
        }
        else{
          cat (paste(patient, " ",site, " ", treatment, " is 0"))
        }
      }
      else{
        cat (paste(patient, " ",site, " is 0"))
      }
    }
  }
}
