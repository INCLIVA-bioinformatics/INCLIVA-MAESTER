#!/usr/bin/Rscript
# -*- coding: utf-8 -*-

# explore_variants_final.R

# Script to select informative variants from different samples and generate various plots.

# Author: UB INCLIVA (fmartinez@incliva.es)
# Date: 18/11/2024
#-------------------------------------------------------------------------------
# Packages
# .libPaths() # to see from where packages are loaded

library(cowplot)
library(ggplot2)
library(ggVennDiagram)
library(ggvenn)
library(gridExtra)
library(matrixStats)
library(RColorBrewer)
library(Seurat)
library(SeuratDisk)
library(reshape2)
library(tidyverse)
library(UpSetR)
library(vcfR)
library(viridis)

get_common_variants <- function(...) {
    # Function to get the common variants between several variant sets
    # Convert all arguments to a list
    variants_list <- list(...)
    # Ensure each list element is a vector by flattening nested lists
    variants_list <- lapply(variants_list, unlist)
    # Find the intersection of all variant sets
    common_variants <- Reduce(intersect, variants_list)

    return(common_variants)
}

get_first_exclusive_variants <- function(...) {
    # Convert all arguments to a list
    variants_list <- list(...)
    # Ensure each list element is a vector by flattening nested lists
    variants_list <- lapply(variants_list, unlist)
    # Find the variants exclusive to the first set
    common_variants <- Reduce(setdiff, variants_list)

    return(common_variants)
}

calcular_conteos_variantes_old <- function(data, conteo_cell_types, variantes) {
    # Function to calculate counts and percentages of variants
    for (variante in variantes) {
        # Create count table for the current variant
        conteo_variante <- as.data.frame(table(data[[variante]], useNA = "ifany"))

        # Define the search function
        buscar_conteo <- function(cell_type, conteo_df) {
            # Create pattern to search for (looking for entries starting with 1 indicating presence)
            patron <- paste0("^1_", cell_type, "$")
            # Extract the corresponding count if it exists
            freq <- conteo_df$Freq[grepl(patron, conteo_df$Var1)]
            # Return the count or NA if not found
            if (length(freq) == 0) {
                return(NA)
            } else {
                return(freq)
            }
        }

        # Count and add to the cell types counts dataframe
        conteo_cell_types[[paste0("count_", variante)]] <- sapply(conteo_cell_types$cell_type, buscar_conteo, conteo_df = conteo_variante)

        # Replace NA with 0
        conteo_cell_types[is.na(conteo_cell_types)] <- 0

        # Calculate percentage
        conteo_cell_types[[paste0("perc_", variante)]] <- (conteo_cell_types[[paste0("count_", variante)]] / sum(conteo_cell_types[[paste0("count_", variante)]])) * 100
    }

    return(conteo_cell_types)
}

calcular_conteos_variantes <- function(data, conteo_cell_types, variantes) {
    # Function to calculate counts, percentages and get barcodes of variants
    for (variante in variantes) {
        # Create count table for the current variant
        conteo_variante <- as.data.frame(table(data[[variante]], useNA = "ifany"))

        # Define the search function
        buscar_conteo <- function(cell_type, conteo_df) {
            # Create pattern to search for (looking for entries starting with 1 indicating presence)
            patron <- paste0("^1_", cell_type, "$")
            # Extract the corresponding count if it exists
            freq <- conteo_df$Freq[grepl(patron, conteo_df$Var1)]
            # Return the count or NA if not found
            if (length(freq) == 0) {
                return(NA)
            } else {
                return(freq)
            }
        }

        # Function to get barcodes of cells with the variant
        obtener_barcodes <- function(cell_type) {
            # Create pattern to find cells containing the variant in rownames
            patron <- paste0("^1_", cell_type, "$")
            barcodes <- rownames(data)[grepl(patron, data[[variante]])]
            # Return barcodes or NA if no cell found
            if (length(barcodes) == 0) {
                return(NA)
            } else {
                return(paste(barcodes, collapse = ", "))
            }
        }

        # Count and add to the cell types counts dataframe
        conteo_cell_types[[paste0("count_", variante)]] <- sapply(conteo_cell_types$cell_type, buscar_conteo, conteo_df = conteo_variante)

        # Replace NA with 0 in counts
        conteo_cell_types[is.na(conteo_cell_types)] <- 0

        # Calculate percentage
        conteo_cell_types[[paste0("perc_", variante)]] <- (conteo_cell_types[[paste0("count_", variante)]] / sum(conteo_cell_types[[paste0("count_", variante)]])) * 100

        # Add the HQcb column with barcodes of cells containing the variant
        conteo_cell_types[[paste0("HQcb_", variante)]] <- sapply(conteo_cell_types$cell_type, obtener_barcodes)
    }

    return(conteo_cell_types)
}

# Plots
plot_heatmap <- function(df, output_file = NULL) {
    # 1. Order columns by the sum of cells containing informative variants (descending)
    df_ordered <- df %>%
        select(names(sort(colSums(df != 0), decreasing = TRUE)))

    # 2. Order rows by variants in order
    df_ordered <- df_ordered[do.call(order, as.data.frame(df_ordered)), ]

    # 3. Remove rows where all variants are 0 (i.e., cells without informative variants)
    df_ordered <- df_ordered[rowSums(df_ordered == 0) != ncol(df_ordered), ]

    # 4. Get the number of remaining cells and variants
    n_celulas <- nrow(df_ordered)
    n_var <- ncol(df_ordered)

    # 5. Create a title for the heatmap
    titulo <- paste(n_celulas, "cells x", n_var, "variants")

    # 6. Generate the heatmap, properly handling whether an output file is provided or not
    if (is.null(output_file)) {
        # If no output file is provided, just display the heatmap
        heatmap <- pheatmap::pheatmap(t(df_ordered),
                                      show_colnames = FALSE,
                                      cluster_rows = FALSE,
                                      cluster_cols = FALSE,
                                      legend = TRUE,
                                      main = titulo,
                                      annotation_legend = FALSE,
                                      col = colorRampPalette(c("white", brewer.pal(5, "YlOrRd")))(255))
    } else {
        # If an output file is provided, save it
        heatmap <- pheatmap::pheatmap(t(df_ordered),
                                      show_colnames = FALSE,
                                      cluster_rows = FALSE,
                                      cluster_cols = FALSE,
                                      legend = TRUE,
                                      main = titulo,
                                      annotation_legend = FALSE,
                                      col = colorRampPalette(c("white", brewer.pal(5, "YlOrRd")))(255),
                                      filename = output_file)
    }

    # 7. Return the generated heatmap (though not necessary if printed immediately)
    return(heatmap)
}

plot_variant_ranked_vaf <- function(variant, af_df) {
  variant_row <- as.data.frame(af_df[, variant])
  rownames(variant_row) <- rownames(af_df)
  colnames(variant_row) <- "af"
  variant_row$celula <- rownames(variant_row) # silly, but without this it gives problems

  variant_row_ordenado <- variant_row[order(variant_row$af),]
  variant_row_ordenado$ranking <- seq_along(variant_row_ordenado$af)

  percentiles <- quantile(variant_row_ordenado$af, probs = c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99))

  ggplot(variant_row_ordenado, aes(x = ranking, y = af)) +
    geom_point() +
    labs(title = variant,
         y = "Variant Allele Frequency",
         x = "Rank sorted cells") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
}

plot_variant_ranked_vaf_q <- function(variant, af_df) {
  # Plot informative variants ordered by VAF showing percentiles
  variant_row <- as.data.frame(af_df[, variant])
  rownames(variant_row) <- rownames(af_df)
  colnames(variant_row) <- "af"
  variant_row$celula <- rownames(variant_row)

  variant_row_ordenado <- variant_row[order(variant_row$af),]
  variant_row_ordenado$ranking <- seq_along(variant_row_ordenado$af)

  percentiles <- quantile(variant_row_ordenado$af, probs = c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99))

  print(summary(variant_row_ordenado$af))  # Show summary of the data

  print(percentiles) # Show percentiles

  n <- nrow(variant_row_ordenado)
  q01_pos <- quantile(seq(1, n), 0.01)
  q10_pos <- quantile(seq(1, n), 0.10)
  q25_pos <- quantile(seq(1, n), 0.25)
  q50_pos <- quantile(seq(1, n), 0.50)
  q75_pos <- quantile(seq(1, n), 0.75)
  q90_pos <- quantile(seq(1, n), 0.90)
  q99_pos <- quantile(seq(1, n), 0.99)

  q01_y <- round(variant_row_ordenado$af[round(q01_pos)], 2)
  q10_y <- round(variant_row_ordenado$af[round(q10_pos)], 2)
  q25_y <- round(variant_row_ordenado$af[round(q25_pos)], 2)
  q50_y <- round(variant_row_ordenado$af[round(q50_pos)], 2)
  q75_y <- round(variant_row_ordenado$af[round(q75_pos)], 2)
  q90_y <- round(variant_row_ordenado$af[round(q90_pos)], 2)
  q99_y <- round(variant_row_ordenado$af[round(q99_pos)], 2)

  ggplot(variant_row_ordenado, aes(x = ranking, y = af)) +
    geom_point() +
    geom_vline(xintercept = c(q01_pos, q10_pos, q90_pos, q99_pos), linetype = "dashed", color = "firebrick3") +
    annotate("text", x = q01_pos, y = q01_y, label = paste("q01 =", q01_y), vjust = -3, color = "firebrick3") +
    annotate("text", x = q10_pos, y = q10_y, label = paste("q10 =", q10_y), vjust = -7, color = "firebrick3") +
    annotate("text", x = q90_pos, y = q90_y, label = paste("q90 =", q90_y), vjust = -1, color = "firebrick3") +
    annotate("text", x = q99_pos, y = q99_y, label = paste("q99 =", q99_y), vjust = 3, color = "firebrick3") +
    labs(title = variant,
         y = "Variant Allele Frequency",
         x = "Rank sorted cells") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
}

plot_scatterplot_informative_variants <- function(tabla_variantes, clone_sizes = seq(1, 200, by = 5), vaf_values = c(0.01, 0.05, 0.10, 0.50, 0.70, 0.90)) {
  # Function to generate a scatter plot to explore selection criteria for informative variants
  # Create a table with all combinations of clone_sizes and vaf_values
  combinations <- expand.grid(min_clone_size = clone_sizes, vaf = vaf_values)

  # Function to calculate the number of cells with informative variants given a minimum clone size and minimum VAF
  calc_informative_cells_vaf <- function(combo) {
    min_clone_size <- combo["min_clone_size"]
    minimum_vaf <- combo["vaf"]

    # Select informative variants for the given VAF
    info_var <- tabla_variantes[, colSums(tabla_variantes > minimum_vaf) >= min_clone_size, drop = FALSE]

    # Filter informative variants that do not exceed 90% of the cells
    info_var <- info_var[, colSums(info_var > 0) <= 0.1 * nrow(info_var), drop = FALSE]

    # Calculate how many cells have any informative variant
    cells_with_informative_variants <- rowSums(info_var > 0) > 0
    sum(cells_with_informative_variants)
  }

  # Apply calc_informative_cells_vaf function to each combination
  results_df <- as.data.frame(
    cbind(
      combinations,
      number_of_cells = apply(combinations, 1, calc_informative_cells_vaf)
    )
  )

  # Convert VAF to percentage for plot labels
  results_df$vaf_label <- paste0(results_df$vaf * 100, "%")

  # Create scatter plot
  scatterplot <- ggplot(results_df, aes(x = min_clone_size, y = number_of_cells, color = vaf_label)) +
    geom_line() +  # Use geom_line to connect points for each series
    geom_point() +  # Optional: add points for each value
    geom_hline(yintercept = nrow(tabla_variantes), color = "black", linetype = "dashed") +  # Black horizontal line
    labs(
      x = "Minimum Clone Size",
      y = "Number of Cells with Informative Variants",
      title = "Selection of informative variants",
      color = "VAF"
    ) +
    theme_minimal() +
    theme(
      text = element_text(size = 14),  # Increase general text size
      axis.title.x = element_text(size = 16),  # X axis title size
      axis.title.y = element_text(size = 16),  # Y axis title size
      plot.title = element_text(size = 18, hjust = 0.5),  # Plot title size and alignment
      legend.title = element_text(size = 14),  # Legend title size
      legend.text = element_text(size = 12)  # Legend text size
    )

  # Return results dataframe and plot
  list(results_df = results_df, scatterplot = scatterplot)
}

plot_variant_frequency <- function(data, column_name) {
    # Function to generate a bar plot of frequencies given a variant
    # Check if the column exists in the data.frame
    if (!column_name %in% colnames(data)) {
        stop("The specified column does not exist in the dataframe.")
    }

    # Extract the column
    col_data <- data[[column_name]]

    # Count frequencies of each value, including NAs
    frecuencias <- as.data.frame(table(col_data, useNA = "ifany"))

    # Rename columns for easier use in ggplot
    colnames(frecuencias) <- c("Value", "Frequency")

    # Create bar plot with ggplot
    ggplot(frecuencias, aes(x = Value, y = Frequency)) +
        geom_bar(stat = "identity", fill = "lightblue") +
        labs(title = paste("Frequencies of", column_name),
             x = column_name,
             y = "Frequency") +
        scale_y_continuous(limits = c(0, max(frecuencias$Frequency + 100, na.rm = TRUE)), expand = c(0, 0)) +  # Adjust Y axis limits
        theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate X axis labels
              plot.title = element_text(hjust = 0.5))  # Center title
}

# 115 - FINAL FIGURES ################################################################################################
## VARIANT TABLES 115 ###############################################################################################
tabla_variantes_path <- "/mnt/zonahpc/home/bioinformatica/servicios/20240724_UB375_MAESTER_Asherman/analysis/13_tablas_variantes"
plot_dir <- "/mnt/zonahpc/home/bioinformatica/servicios/20240724_UB375_MAESTER_Asherman/plots/variants"
patient <- "115"
samples <- c("CD133", "entire_E1", "entire_E2", "entire_S1", "entire_S2")

# Read the variant tables of all samples into a list
tabla_variantes_list <- lapply(samples, function(sample) {
    read.csv(file = paste0(tabla_variantes_path, "/", patient, "_", sample, "/tabla_AF_variantes.csv"), header = TRUE, check.names = FALSE, row.names = 1)
})

names(tabla_variantes_list) <- paste0(patient, "_", samples) # label the samples in the list

# Remove indels
filtered_list <- lapply(tabla_variantes_list, function(df) {
    # Create a vector with the names of columns to remove
    cols_to_remove <- colnames(df)[sapply(colnames(df), function(col) {
        # Remove numbers at the start
        new_col <- gsub("^[0-9]+", "", col)
        # Split by '>'
        parts <- strsplit(new_col, split = ">")[[1]]
        # Compare length
        nchar(parts[1]) != nchar(parts[2])
    })]

    # Remove the selected columns from the dataframe
    df_filtered <- df %>% select(-one_of(cols_to_remove))

    # Return the filtered dataframe
    return(df_filtered)
})

lapply(tabla_variantes_list, ncol)
lapply(filtered_list, ncol)

# Overwrite the list with filtered variant tables
tabla_variantes_list <- filtered_list

## VARIANTS OF 115_E2 ################################################################################################
# Informative variants: initially we will look for informative variants among those common to CD133 and entire_E2 BUT that are NOT present in entire_S1 (since we DON'T have entire_S1)
ids_variantes_list <- lapply(tabla_variantes_list, function(df) colnames(df)) # extract variant IDs

common_115_CD133_E2 <- get_common_variants(ids_variantes_list["115_CD133"], ids_variantes_list["115_entire_E2"]) # common between CD133 and entire_E2
common_115_CD133_E2_filt <- get_first_exclusive_variants(common_115_CD133_E2, ids_variantes_list["115_entire_E1"]) # remove those common with E1

# Checks
table(common_115_CD133_E2 %in% colnames(tabla_variantes_list[["115_CD133"]])) # 506
table(colnames(tabla_variantes_list[["115_CD133"]]) %in% common_115_CD133_E2)
table(common_115_CD133_E2 %in% colnames(tabla_variantes_list[["115_entire_E2"]])) # 506
table(colnames(tabla_variantes_list[["115_entire_E2"]]) %in% common_115_CD133_E2)

# PARAMETERS FOR INFORMATIVE VARIANTS
minimum_clone_size <- 5 # y, (e.g. 50)
minimum_vaf <- 0.10 # x, (e.g. 0.50)
vaf_2 <- 0 # w, (e.g., 0.01)
popul_percent <- 0.10 # z, (e.g. 0.75)

info_var_115_E2 <- tabla_variantes_list[["115_entire_E2"]][, colSums(tabla_variantes_list[["115_entire_E2"]] > minimum_vaf) >= minimum_clone_size]
info_var_115_E2 <- info_var_115_E2[, colSums(info_var_115_E2 > vaf_2) <= popul_percent * nrow(info_var_115_E2)]

ncol(info_var_115_E2) # 193 informative variants
colnames(info_var_115_E2)

nrow(info_var_115_E2)
sum(rowSums(info_var_115_E2 > 0) > 0) # 6715 tracked cells

table(colnames(info_var_115_E2) %in% common_115_CD133_E2_filt) # TRUE 19, 19 variants common between CD133 and E2 that are NOT in E1 are informative
info_var_115_E2_common <- info_var_115_E2[, colnames(info_var_115_E2) %in% common_115_CD133_E2_filt, drop = FALSE]
colnames(info_var_115_E2_common) # "152T>C"

nrow(info_var_115_E2_common)
sum(rowSums(info_var_115_E2_common > 0) > 0) # 400 tracked cells

## VARIANTS OF 115_S2 #################################################################################################
# Informative variants: initially we will look for informative variants among those common to CD133 and entire_S2 BUT that are NOT present in entire_S1
ids_variantes_list <- lapply(tabla_variantes_list, function(df) colnames(df)) # extract variant IDs

common_115_CD133_S2 <- get_common_variants(ids_variantes_list["115_CD133"], ids_variantes_list["115_entire_S2"]) # common between CD133 and entire_S2
common_115_CD133_S2_filt <- get_first_exclusive_variants(common_115_CD133_S2, ids_variantes_list["115_entire_S1"]) # remove those common with E1/S1

# Checks
table(common_115_CD133_S2 %in% colnames(tabla_variantes_list[["115_CD133"]])) # 429
table(colnames(tabla_variantes_list[["115_CD133"]]) %in% common_115_CD133_S2)
table(common_115_CD133_S2 %in% colnames(tabla_variantes_list[["115_entire_S2"]])) # 429
table(colnames(tabla_variantes_list[["115_entire_S2"]]) %in% common_115_CD133_S2)

# PARAMETERS FOR INFORMATIVE VARIANTS
minimum_clone_size <- 5 # y, (e.g. 50)
minimum_vaf <- 0.1 # x, (e.g. 0.50), maybe modify to 10%
vaf_2 <- 0 # w, (e.g., 0.01)
popul_percent <- 0.10 # z, (e.g. 0.75)

info_var_115_S2 <- tabla_variantes_list[["115_entire_S2"]][, colSums(tabla_variantes_list[["115_entire_S2"]] >= minimum_vaf) >= minimum_clone_size] # want it in at least these 50 cells
info_var_115_S2 <- info_var_115_S2[, colSums(info_var_115_S2 > vaf_2) <= popul_percent * nrow(info_var_115_S2)] # want it in less than 75%

ncol(info_var_115_S2) # 187 informative variants
colnames(info_var_115_S2)
sum(rowSums(info_var_115_S2 > 0) > 0) # 5690 tracked cells

table(colnames(info_var_115_S2) %in% common_115_CD133_S2_filt) # TRUE 1, 1 variant common between CD133 and S2 that is NOT in S1 is informative
info_var_115_S2_common <- info_var_115_S2[, colnames(info_var_115_S2) %in% common_115_CD133_S2_filt, drop = FALSE] # drop = FALSE so it does not become a vector
colnames(info_var_115_S2_common) # "11711G>A"
sum(rowSums(info_var_115_S2_common > 0) > 0) # 439 tracked cells

## SEURAT 115 E2 + S2 ##################################################################################################
datos_seurat_path <- "/mnt/zonahpc/home/bioinformatica/servicios/20240724_UB375_MAESTER_Asherman/GSE215968/AS_pre_post.H5Seurat"
seurat_GSE215968_raw <- LoadH5Seurat(datos_seurat_path) # load data into a Seurat object, takes a while

table(seurat_GSE215968_raw$cell_type_def) # AS_Epithelium, B_cells, CD8T_cells, Ciliated_Epithelium...
table(seurat_GSE215968_raw$Patient) # 01-05 01-07 01-08 01-09 01-10 01-12 01-13 01-14 01-15
table(seurat_GSE215968_raw$Cell_type) # Epithelium, Stroma, # rename this to biopsy_origin
table(seurat_GSE215968_raw$Treatment_stage) # A-Pre, B-Post

table(seurat_GSE215968_raw$Cell_type, seurat_GSE215968_raw$Patient, seurat_GSE215968_raw$Treatment_stage) # check table, is the number of cells expected? Remember that our 115 is actually 114 in the object.

# To add variants to the Seurat object
nrow(info_var_115_E2) # 12415
nrow(info_var_115_S2) # 12167
rownames(info_var_115_E2) <- paste0("ENTIRE-1-14-E2_", rownames(info_var_115_E2))
rownames(info_var_115_S2) <- paste0("ENTIRE-1-14-S2_", rownames(info_var_115_S2))

# Are there common variants?
length(intersect(colnames(info_var_115_E2), colnames(info_var_115_S2))) # 90, yes
# Common cells?
length(intersect(rownames(info_var_115_E2), rownames(info_var_115_S2))) # 0, no

# Merge variant tables E2 and S2
merged_samples_variants <- bind_rows(info_var_115_E2, info_var_115_S2)
ncol(merged_samples_variants) # 290
nrow(merged_samples_variants) # 24582

total_cells <- rownames(merged_samples_variants)
head(total_cells)
tail(total_cells)

# IMPORTANT: remember that 14 in the Seurat object is patient 115.
GSE215968_raw_115_post <- subset(seurat_GSE215968_raw, subset = Patient == "01-14" & Treatment_stage == "B-Post") # do not separate E2 and S2
GSE215968_raw_115_post <- subset(GSE215968_raw_115_post, cells = total_cells) # take post-115 cells
table(total_cells %in% colnames(GSE215968_raw_115_post)) # all TRUE

dimplot_115_post <- DimPlot(GSE215968_raw_115_post,
                    reduction = "umap",
                    group.by = "cell_type_def") +
                    labs(title = "N15",
                    color = "", # legend title
                    x = "",
                    y = "")

ggsave(paste0("/mnt/zonahpc/home/bioinformatica/servicios/20240724_UB375_MAESTER_Asherman/plots/umaps_", patient,"/umap_115_cell_type.png"), plot = dimplot_115_post, width = 11, height = 9, dpi = 300)

head(colnames(GSE215968_raw_115_post))
tail(colnames(GSE215968_raw_115_post))

merged_samples_variants[is.na(merged_samples_variants)] <- 0
head(merged_samples_variants)
tail(merged_samples_variants)

# UMAP
colnames(GSE215968_raw_115_post@meta.data)

trazado_115_E2 <- read.csv("/mnt/zonahpc/home/bioinformatica/servicios/20240724_UB375_MAESTER_Asherman/trazado_tipos_celulares/conteo_cell_types_115_E2_HQcb.csv", sep = ";")
trazado_115_S2 <- read.csv("/mnt/zonahpc/home/bioinformatica/servicios/20240724_UB375_MAESTER_Asherman/trazado_tipos_celulares/conteo_cell_types_115_S2_HQcb.csv", sep = ";")

trazado_115_E2[is.na(trazado_115_E2)] <- 0 # last column has NAs for some reason
trazado_115_S2[is.na(trazado_115_S2)] <- 0 # last column has NAs for some reason

head(trazado_115_E2)
head(trazado_115_S2)

trazado_115_E2_renamed <- trazado_115_E2 %>%
    select(starts_with("HQcb")) %>%
    rename_with(~ gsub("^HQcb_", "", .x)) %>%
    rename_with(~ paste0(.x, "_E2"))
trazado_115_S2_renamed <- trazado_115_S2 %>%
    select(starts_with("HQcb")) %>%
    rename_with(~ gsub("^HQcb_", "", .x)) %>%
    rename_with(~ paste0(.x, "_S2"))

trazado_115_E2_procesado <- trazado_115_E2_renamed %>%
    summarise(across(everything(), ~ paste(.x[.x != 0], collapse = ",")))
trazado_115_S2_procesado <- trazado_115_S2_renamed %>%
    summarise(across(everything(), ~ paste(.x[.x != 0], collapse = ",")))

# HARDCODED, REVIEW MANUALLY -> Remove variants that we saw do not trace CD133 cells
ncol(trazado_115_E2_procesado) # 9
trazado_115_E2_procesado <- trazado_115_E2_procesado %>%
    select(-c("2623A.G_E2", "13436C.G_E2"))
ncol(trazado_115_E2_procesado) # 7

ncol(trazado_115_S2_procesado) # 23
trazado_115_S2_procesado <- trazado_115_S2_procesado %>%
    select(-c("13436C.G_S2", "15731G.A_S2", "2833A.G_S2", "15734G.A_S2"))
ncol(trazado_115_S2_procesado) # 19

# Join both samples into a single dataframe with the HQcb
trazado_115_all <- cbind(trazado_115_E2_procesado, trazado_115_S2_procesado)
ncol(trazado_115_all) # 26

# Convert to a more usable list
barcode_lists_115 <- lapply(trazado_115_all, function(x) unlist(strsplit(x, ",")))
barcode_lists_115 <- lapply(barcode_lists_115, trimws) # Remove whitespace

cell_barcodes <- rownames(GSE215968_raw_115_post@meta.data)

# Create the highlight column
GSE215968_raw_115_post@meta.data$highlight <- sapply(cell_barcodes, function(cell) {
    # Find to which variants the cell belongs
    variantes <- names(barcode_lists_115)[sapply(barcode_lists_115, function(v) cell %in% v)]

    # If it belongs to at least one variant, concatenate names; else assign "0"
    if (length(variantes) > 0) {
        return(paste(variantes, collapse = ","))
    } else {
        return("0")
    }
})

# At this point we have the labels in the cells to plot with the "_" indicating from which cell they come
table(GSE215968_raw_115_post@meta.data$highlight)

GSE215968_raw_115_post@meta.data$highlight <- GSE215968_raw_115_post@meta.data$highlight %>%
    str_replace("_.*", "") %>%
    str_replace_all("\\.", ">")

table(GSE215968_raw_115_post@meta.data$highlight)

# Get unique levels of "highlight"
highlight_levels <- unique(GSE215968_raw_115_post@meta.data$highlight)

# Generate colors: "0" in grey, rest with extended palette
highlight_colors <- setNames(
    c("grey", viridis(length(highlight_levels) - 1)),
    c("0", highlight_levels[highlight_levels != "0"])
)

# Create the UMAP plot
umap_custom <- DimPlot(
    GSE215968_raw_115_post,
    reduction = "umap",
    group.by = "highlight",
    cols = highlight_colors,
    pt.size = 0.5
) +
    labs(
        title = "N15", # Change title
        color = "",    # Change legend title
        x = "",       # Remove X axis title
        y = ""        # Remove Y axis title
    ) +
    guides(
        color = guide_legend(override.aes = list(size = 3))  # Adjust legend point size
    ) +
    scale_color_manual(
        values = highlight_colors,
        breaks = names(highlight_colors)[names(highlight_colors) != "0"]  # Exclude "0" from legend
    ) +
    geom_point(
        data = GSE215968_raw_115_post@meta.data %>%
            mutate(UMAP_1 = Embeddings(GSE215968_raw_115_post, "umap")[, 1],
                   UMAP_2 = Embeddings(GSE215968_raw_115_post, "umap")[, 2]) %>%
            filter(highlight != "0"),  # Filter points to highlight
        aes(x = UMAP_1, y = UMAP_2, color = highlight), # Use highlight colors
        size = 1.25,                                    # Larger size for highlighted points
        inherit.aes = FALSE                             # Avoid conflicts with DimPlot mapping
    )

# Show custom plot
# print(umap_custom)

ggsave(paste0("/mnt/zonahpc/home/bioinformatica/servicios/20240724_UB375_MAESTER_Asherman/plots/umaps_", patient,"/umap_115_variants.png"), plot = umap_custom, width = 11, height = 9, dpi = 300)

# 109 - FINAL FIGURES ################################################################################################
## VARIANT TABLES 109 #############################################################################################
tabla_variantes_path <- "/mnt/zonahpc/home/bioinformatica/servicios/20240724_UB375_MAESTER_Asherman/analysis/13_tablas_variantes"
plot_dir <- "/mnt/zonahpc/home/bioinformatica/servicios/20240724_UB375_MAESTER_Asherman/plots/variants"
patient <- "109"
samples <- c("CD133", "entire_E1", "entire_E2", "entire_S1", "entire_S2")

# Read the variant table of all samples into a list
tabla_variantes_list <- lapply(samples, function(sample) {
    read.csv(file = paste0(tabla_variantes_path, "/", patient, "_", sample, "/tabla_AF_variantes.csv"), header = TRUE, check.names = FALSE, row.names = 1)
})

names(tabla_variantes_list) <- paste0(patient, "_", samples) # label the samples in the list

# Remove indels
filtered_list <- lapply(tabla_variantes_list, function(df) {
    # Create a vector with the names of the columns to remove
    cols_to_remove <- colnames(df)[sapply(colnames(df), function(col) {
        # Remove numbers from the start
        new_col <- gsub("^[0-9]+", "", col)
        # Split by '>'
        parts <- strsplit(new_col, split = ">")[[1]]
        # Compare length
        nchar(parts[1]) != nchar(parts[2])
    })]

    # Remove the selected columns from the dataframe
    df_filtered <- df %>% select(-one_of(cols_to_remove))

    # Return the filtered dataframe
    return(df_filtered)
})

lapply(tabla_variantes_list, ncol)
lapply(filtered_list, ncol)

# Overwrite the list with the variant tables
tabla_variantes_list <- filtered_list

#### VARIANTS OF 109_E2 ###############################################################################################
# Informative variants: initially we will look for informative variants among the variants common to CD133 and entire_E2 BUT that are NOT in entire_S1 (since we do NOT have entire_S1)
ids_variantes_list <- lapply(tabla_variantes_list, function(df) colnames(df)) # extract variant identifiers

common_109_CD133_E2 <- get_common_variants(ids_variantes_list["109_CD133"], ids_variantes_list["109_entire_E2"]) # common between CD133 and entire_E2
common_109_CD133_E2_filt <- get_first_exclusive_variants(common_109_CD133_E2, ids_variantes_list["109_entire_E1"]) # remove those common with E1

# PARAMETERS FOR INFORMATIVE VARIANTS
minimum_clone_size <- 5 # y, (e.g. 50)
minimum_vaf <- 0.10 # x, (e.g. 0.50)
vaf_2 <- 0 # w, (e.g., 0.01)
popul_percent <- 0.10 # z, (e.g. 0.75)

info_var_109_E2 <- tabla_variantes_list[["109_entire_E2"]][, colSums(tabla_variantes_list[["109_entire_E2"]] > minimum_vaf) >= minimum_clone_size]
info_var_109_E2 <- info_var_109_E2[, colSums(info_var_109_E2 > vaf_2) <= popul_percent * nrow(info_var_109_E2)]

ncol(info_var_109_E2) # 13 informative variants
colnames(info_var_109_E2)

nrow(info_var_109_E2)
sum(rowSums(info_var_109_E2 > 0) > 0) # 63 tracked cells

table(colnames(info_var_109_E2) %in% common_109_CD133_E2_filt) # TRUE 2, 2 variants common between CD133 and E2 that are NOT in E1 are informative
info_var_109_E2_common <- info_var_109_E2[, colnames(info_var_109_E2) %in% common_109_CD133_E2_filt, drop = FALSE]
colnames(info_var_109_E2_common) # "8276C>A"  "16278C>T"

nrow(info_var_109_E2_common)
sum(rowSums(info_var_109_E2_common > 0) > 0) # 14 tracked cells

## VARIANTS OF 109_S2 #################################################################################################
# Informative variants: initially we will look for informative variants among the variants common to CD133 and entire_S2 BUT that are NOT in entire_S1
ids_variantes_list <- lapply(tabla_variantes_list, function(df) colnames(df)) # extract variant identifiers

common_109_CD133_S2 <- get_common_variants(ids_variantes_list["109_CD133"], ids_variantes_list["109_entire_S2"]) # common between CD133 and entire_S2
common_109_CD133_S2_filt <- get_first_exclusive_variants(common_109_CD133_S2, ids_variantes_list["109_entire_S1"]) # remove those common with E1/S1

# Checks
table(common_109_CD133_S2 %in% colnames(tabla_variantes_list[["109_CD133"]])) # 429
table(colnames(tabla_variantes_list[["109_CD133"]]) %in% common_109_CD133_S2)
table(common_109_CD133_S2 %in% colnames(tabla_variantes_list[["109_entire_S2"]])) # 429
table(colnames(tabla_variantes_list[["109_entire_S2"]]) %in% common_109_CD133_S2)

# PARAMETERS FOR INFORMATIVE VARIANTS
minimum_clone_size <- 5 # y, (e.g. 50)
minimum_vaf <- 0.1 # x, (e.g. 0.50), maybe modify to 10%
vaf_2 <- 0 # w, (e.g., 0.01)
popul_percent <- 0.10 # z, (e.g. 0.75)

info_var_109_S2 <- tabla_variantes_list[["109_entire_S2"]][, colSums(tabla_variantes_list[["109_entire_S2"]] >= minimum_vaf) >= minimum_clone_size] # I want it to be in at least these 50 cells
info_var_109_S2 <- info_var_109_S2[, colSums(info_var_109_S2 > vaf_2) <= popul_percent * nrow(info_var_109_S2)] # I want it to be in less than 75%

ncol(info_var_109_S2) # informative variants
colnames(info_var_109_S2)
sum(rowSums(info_var_109_S2 > 0) > 0) # tracked cells

table(colnames(info_var_109_S2) %in% common_109_CD133_S2_filt) # TRUE 1, 1 variant common between CD133 and S2 that are NOT in S1 are informative
info_var_109_S2_common <- info_var_109_S2[, colnames(info_var_109_S2) %in% common_109_CD133_S2_filt, drop = FALSE] # drop = FALSE so it doesn't become a vector
colnames(info_var_109_S2_common) # ""
sum(rowSums(info_var_109_S2_common > 0) > 0) # 843 tracked cells

## SEURAT 109 S2 ##################################################################################################
datos_seurat_path <- "/mnt/zonahpc/home/bioinformatica/servicios/20240724_UB375_MAESTER_Asherman/GSE215968/AS_pre_post.H5Seurat"
seurat_GSE215968_raw <- LoadH5Seurat(datos_seurat_path) # load data into a Seurat object, it takes a little while

table(seurat_GSE215968_raw$cell_type_def) # AS_Epithelium, B_cells, CD8T_cells, Ciliated_Epithelium...
table(seurat_GSE215968_raw$Patient) # 01-05 01-07 01-08 01-09 01-10 01-12 01-13 01-14 01-15
table(seurat_GSE215968_raw$Cell_type) # Epithelium, Stroma, # change this name to biopsy_origin
table(seurat_GSE215968_raw$Treatment_stage) # A-Pre, B-Post

table(seurat_GSE215968_raw$Cell_type, seurat_GSE215968_raw$Patient, seurat_GSE215968_raw$Treatment_stage) # TODO - review table

# To add variants to the Seurat object
nrow(info_var_109_E2) # 143
nrow(info_var_109_S2) # 12167
rownames(info_var_109_E2) <- paste0("ENTIRE-1-09-E2_", rownames(info_var_109_E2))
rownames(info_var_109_S2) <- paste0("ENTIRE-1-09-S2_", rownames(info_var_109_S2))

# Are there common variants?
length(intersect(colnames(info_var_109_E2), colnames(info_var_109_S2))) # 90, yes
# Common cells?
length(intersect(rownames(info_var_109_E2), rownames(info_var_109_S2))) # 0, no

# Merge variant tables E2 and S2
merged_samples_variants <- bind_rows(info_var_109_E2, info_var_109_S2)
ncol(merged_samples_variants) # 59
nrow(merged_samples_variants) # 2474

total_cells <- rownames(merged_samples_variants)
head(total_cells)
tail(total_cells)

# IMPORTANT: remember that 14 in the Seurat object corresponds to patient 109.
GSE215968_raw_109_post <- subset(seurat_GSE215968_raw, subset = Patient == "01-09" & Treatment_stage == "B-Post" & Cell_type == "Stroma") # here we take only S2 because we don't have E2 variants
GSE215968_raw_109_post <- subset(GSE215968_raw_109_post, cells = total_cells) # take post-109 cells
table(total_cells %in% colnames(GSE215968_raw_109_post)) # all true

dimplot_109_post <- DimPlot(GSE215968_raw_109_post,
                    reduction = "umap",
                    group.by = "cell_type_def") +
                    labs(title = "N09",
                    color = "", # legend title
                    x = "",
                    y = "")

ggsave(paste0("/mnt/zonahpc/home/bioinformatica/servicios/20240724_UB375_MAESTER_Asherman/plots/umaps_", patient,"/umap_109_cell_type.png"), plot = dimplot_109_post, width = 11, height = 9, dpi = 300)

head(colnames(GSE215968_raw_109_post))
tail(colnames(GSE215968_raw_109_post))

merged_samples_variants[is.na(merged_samples_variants)] <- 0
head(merged_samples_variants)
tail(merged_samples_variants)

# UMAP
colnames(GSE215968_raw_109_post@meta.data)

# trazado_109_E2 <- read.csv("/mnt/zonahpc/home/bioinformatica/servicios/20240724_UB375_MAESTER_Asherman/trazado_tipos_celulares/conteo_cell_types_109_E2_HQcb.csv", sep = ";")
trazado_109_S2 <- read.csv("/mnt/zonahpc/home/bioinformatica/servicios/20240724_UB375_MAESTER_Asherman/trazado_tipos_celulares/conteo_cell_types_109_S2_HQcb.csv", sep = ";")

# trazado_109_E2[is.na(trazado_109_E2)] <- 0 # Last column has NAs for some reason
trazado_109_S2[is.na(trazado_109_S2)] <- 0 # Last column has NAs for some reason

# head(trazado_109_E2)
head(trazado_109_S2)

# trazado_109_E2_renamed <- trazado_109_E2 %>%
#     select(starts_with("HQcb")) %>%
#     rename_with(~ gsub("^HQcb_", "", .x)) %>%
#     rename_with(~ paste0(.x, "_E2"))
trazado_109_S2_renamed <- trazado_109_S2 %>%
    select(starts_with("HQcb")) %>%
    rename_with(~ gsub("^HQcb_", "", .x)) %>%
    rename_with(~ paste0(.x, "_S2"))

# trazado_109_E2_procesado <- trazado_109_E2_renamed %>%
#     summarise(across(everything(), ~ paste(.x[.x != 0], collapse = ",")))
trazado_109_S2_procesado <- trazado_109_S2_renamed %>%
    summarise(across(everything(), ~ paste(.x[.x != 0], collapse = ",")))

# HARDCODED, REVIEW MANUALLY -> Remove variants that we saw don't trace CD133 cells
# ncol(trazado_109_E2_procesado) # 9
# trazado_109_E2_procesado <- trazado_109_E2_procesado %>%
#     select(-c("2623A.G_E2", "13436C.G_E2"))
# ncol(trazado_109_E2_procesado) # 7

ncol(trazado_109_S2_procesado) # 9
trazado_109_S2_procesado <- trazado_109_S2_procesado %>%
    select(c("10326T.C_S2", "13218A.T_S2"))
ncol(trazado_109_S2_procesado) # 2

# Combine both samples into a single dataframe with the HQcb
# trazado_109_all <- cbind(trazado_109_E2_procesado, trazado_109_S2_procesado)
# ncol(trazado_109_all) # 26

trazado_109_all <- trazado_109_S2_procesado

# Convert to a more usable list
barcode_lists_109 <- lapply(trazado_109_all, function(x) unlist(strsplit(x, ",")))
barcode_lists_109 <- lapply(barcode_lists_109, trimws) # Remove white spaces

cell_barcodes <- rownames(GSE215968_raw_109_post@meta.data)

# Create the highlight column
GSE215968_raw_109_post@meta.data$highlight <- sapply(cell_barcodes, function(cell) {
    # Find which variants the cell belongs to
    variantes <- names(barcode_lists_109)[sapply(barcode_lists_109, function(v) cell %in% v)]

    # If it belongs to at least one variant, concatenate the names; otherwise, assign "0"
    if (length(variantes) > 0) {
        return(paste(variantes, collapse = ","))
    } else {
        return("0")
    }
})

# At this point we have the labels on the cells to be plotted with the "_" indicating which cell they come from
table(GSE215968_raw_109_post@meta.data$highlight)

GSE215968_raw_109_post@meta.data$highlight <- GSE215968_raw_109_post@meta.data$highlight %>%
    str_replace("_.*", "") %>%
    str_replace_all("\\.", ">")

table(GSE215968_raw_109_post@meta.data$highlight)

# Get unique levels of "highlight"
highlight_levels <- unique(GSE215968_raw_109_post@meta.data$highlight)

# Generate colors: "0" in grey, the rest with an extended palette
highlight_colors <- setNames(
    c("grey", viridis(length(highlight_levels) - 1)),
    c("0", highlight_levels[highlight_levels != "0"])
)

# Create the UMAP plot
umap_custom <- DimPlot(
    GSE215968_raw_109_post,
    reduction = "umap",
    group.by = "highlight",
    cols = highlight_colors,
    pt.size = 0.5
) +
    labs(
        title = "N09", # Change title
        color = "",    # Change legend title
        x = "",       # Remove x axis title
        y = ""        # Remove y axis title
    ) +
    guides(
        color = guide_legend(override.aes = list(size = 3))  # Adjust legend point size
    ) +
    scale_color_manual(
        values = highlight_colors,
        breaks = names(highlight_colors)[names(highlight_colors) != "0"]  # Exclude "0" from legend
    ) +
    geom_point(
        data = GSE215968_raw_109_post@meta.data %>%
            mutate(UMAP_1 = Embeddings(GSE215968_raw_109_post, "umap")[, 1],
                   UMAP_2 = Embeddings(GSE215968_raw_109_post, "umap")[, 2]) %>%
            filter(highlight != "0"),  # Filter points to highlight
        aes(x = UMAP_1, y = UMAP_2, color = highlight), # Respect "highlight" colors
        size = 1.25,                                    # Increased size for highlighted points
        inherit.aes = FALSE                             # Avoid conflicts with DimPlot mapping
    )

# Show the customized plot
# print(umap_custom)

ggsave(paste0("/mnt/zonahpc/home/bioinformatica/servicios/20240724_UB375_MAESTER_Asherman/plots/umaps_", patient,"/umap_109_variants.png"), plot = umap_custom, width = 11, height = 9, dpi = 300)

# 107 - FINAL FIGURES ################################################################################################
## VARIANT TABLES 107 #############################################################################################
tabla_variantes_path <- "/mnt/zonahpc/home/bioinformatica/servicios/20240724_UB375_MAESTER_Asherman/analysis/13_tablas_variantes"
plot_dir <- "/mnt/zonahpc/home/bioinformatica/servicios/20240724_UB375_MAESTER_Asherman/plots/variants"
patient <- "107"
samples <- c("CD133", "entire_E2", "entire_S1", "entire_S2")

# Read the variant table of all samples into a list
tabla_variantes_list <- lapply(samples, function(sample) {
    read.csv(file = paste0(tabla_variantes_path, "/", patient, "_", sample, "/tabla_AF_variantes.csv"), header = TRUE, check.names = FALSE, row.names = 1)
})

names(tabla_variantes_list) <- paste0(patient, "_", samples) # label the samples in the list

# Remove indels
filtered_list <- lapply(tabla_variantes_list, function(df) {
    # Create a vector with the names of the columns to remove
    cols_to_remove <- colnames(df)[sapply(colnames(df), function(col) {
        # Remove numbers from the beginning
        new_col <- gsub("^[0-9]+", "", col)
        # Split by '>'
        parts <- strsplit(new_col, split = ">")[[1]]
        # Compare length
        nchar(parts[1]) != nchar(parts[2])
    })]

    # Remove the selected columns from the dataframe
    df_filtered <- df %>% select(-one_of(cols_to_remove))

    # Return the filtered dataframe
    return(df_filtered)
})

lapply(tabla_variantes_list, ncol)
lapply(filtered_list, ncol)

# Overwrite the list with the filtered variant tables
tabla_variantes_list <- filtered_list

## VARIANTS OF 107_E2 #################################################################################################
# Informative variants: initially, we will look for informative variants among the variants common to CD133 and entire_E2 BUT that are NOT in entire_S1 (since we do NOT have entire_S1)
ids_variantes_list <- lapply(tabla_variantes_list, function(df) colnames(df)) # extract the variant identifiers

common_107_CD133_E2 <- get_common_variants(ids_variantes_list["107_CD133"], ids_variantes_list["107_entire_E2"]) # common between CD133 and entire_E2
common_107_CD133_E2_filt <- get_first_exclusive_variants(common_107_CD133_E2, ids_variantes_list["107_entire_S1"]) # remove those common with S1

# PARAMETERS FOR INFORMATIVE VARIANTS
minimum_clone_size <- 5 # y, (e.g. 50)
minimum_vaf <- 0.10 # x, (e.g. 0.50)
vaf_2 <- 0 # w, (e.g., 0.01)
popul_percent <- 0.10 # z, (e.g. 0.75)

info_var_107_E2 <- tabla_variantes_list[["107_entire_E2"]][, colSums(tabla_variantes_list[["107_entire_E2"]] > minimum_vaf) >= minimum_clone_size]
info_var_107_E2 <- info_var_107_E2[, colSums(info_var_107_E2 > vaf_2) <= popul_percent * nrow(info_var_107_E2)]

ncol(info_var_107_E2) # 72 informative variants
colnames(info_var_107_E2)

nrow(info_var_107_E2)
sum(rowSums(info_var_107_E2 > 0) > 0) # 2215 tracked cells

table(colnames(info_var_107_E2) %in% common_107_CD133_E2_filt) # TRUE 19, 19 variants common between CD133 and E2 that are NOT in S1 are informative
info_var_107_E2_common <- info_var_107_E2[, colnames(info_var_107_E2) %in% common_107_CD133_E2_filt, drop = FALSE]
colnames(info_var_107_E2_common) # "152T>C"

nrow(info_var_107_E2_common)
sum(rowSums(info_var_107_E2_common > 0) > 0) # 488 tracked cells

## VARIANTS OF 107_S2 #################################################################################################
# Informative variants: initially, we will look for informative variants among the variants common to CD133 and entire_S2 BUT that are NOT in entire_S1
ids_variantes_list <- lapply(tabla_variantes_list, function(df) colnames(df)) # extract the variant identifiers

common_107_CD133_S2 <- get_common_variants(ids_variantes_list["107_CD133"], ids_variantes_list["107_entire_S2"]) # common between CD133 and entire_S2
common_107_CD133_S2_filt <- get_first_exclusive_variants(common_107_CD133_S2, ids_variantes_list["107_entire_S1"]) # remove those common with E1/S1

# Checks
table(common_107_CD133_S2 %in% colnames(tabla_variantes_list[["107_CD133"]])) # 429
table(colnames(tabla_variantes_list[["107_CD133"]]) %in% common_107_CD133_S2)
table(common_107_CD133_S2 %in% colnames(tabla_variantes_list[["107_entire_S2"]])) # 429
table(colnames(tabla_variantes_list[["107_entire_S2"]]) %in% common_107_CD133_S2)

# PARAMETERS FOR INFORMATIVE VARIANTS
minimum_clone_size <- 5 # y, (e.g. 50)
minimum_vaf <- 0.1 # x, (e.g. 0.50), maybe change to 10%
vaf_2 <- 0 # w, (e.g., 0.01)
popul_percent <- 0.10 # z, (e.g. 0.75)

info_var_107_S2 <- tabla_variantes_list[["107_entire_S2"]][, colSums(tabla_variantes_list[["107_entire_S2"]] >= minimum_vaf) >= minimum_clone_size] # I want it to be in at least these 50 cells
info_var_107_S2 <- info_var_107_S2[, colSums(info_var_107_S2 > vaf_2) <= popul_percent * nrow(info_var_107_S2)] # I want it to be in less than 75%

ncol(info_var_107_S2) # informative variants
colnames(info_var_107_S2)
sum(rowSums(info_var_107_S2 > 0) > 0) # tracked cells

table(colnames(info_var_107_S2) %in% common_107_CD133_S2_filt) # TRUE 1, 1 variant common between CD133 and S2 that is NOT in S1 is informative
info_var_107_S2_common <- info_var_107_S2[, colnames(info_var_107_S2) %in% common_107_CD133_S2_filt, drop = FALSE] # drop = FALSE so it doesn't become a vector
colnames(info_var_107_S2_common) # ""
sum(rowSums(info_var_107_S2_common > 0) > 0) # 313 tracked cells

## SEURAT 107 E2 + S2 ##################################################################################################
datos_seurat_path <- "/mnt/zonahpc/home/bioinformatica/servicios/20240724_UB375_MAESTER_Asherman/GSE215968/AS_pre_post.H5Seurat"
seurat_GSE215968_raw <- LoadH5Seurat(datos_seurat_path) # load the data into a Seurat object, this takes a while

table(seurat_GSE215968_raw$cell_type_def) # AS_Epithelium, B_cells, CD8T_cells, Ciliated_Epithelium...
table(seurat_GSE215968_raw$Patient) # 01-05 01-07 01-08 01-09 01-10 01-12 01-13 01-14 01-15
table(seurat_GSE215968_raw$Cell_type) # Epithelium     Stroma, # change this name to biopsy_origin
table(seurat_GSE215968_raw$Treatment_stage) # A-Pre, B-Post

table(seurat_GSE215968_raw$Cell_type, seurat_GSE215968_raw$Patient, seurat_GSE215968_raw$Treatment_stage) # check table, is the number of cells as expected?

# To add the variants to the Seurat object
nrow(info_var_107_E2) # 3607
nrow(info_var_107_S2) # 2944
rownames(info_var_107_E2) <- paste0("ENTIRE-1-07-E2_", rownames(info_var_107_E2))
rownames(info_var_107_S2) <- paste0("ENTIRE-1-07-S2_", rownames(info_var_107_S2))

# Are there common variants?
length(intersect(colnames(info_var_107_E2), colnames(info_var_107_S2))) # 49, yes
# Common cells?
length(intersect(rownames(info_var_107_E2), rownames(info_var_107_S2))) # 0, no

# Merge E2 and S2 variant tables
merged_samples_variants <- bind_rows(info_var_107_E2, info_var_107_S2)
ncol(merged_samples_variants) # 80
nrow(merged_samples_variants) # 6551

total_cells <- rownames(merged_samples_variants)
head(total_cells)
tail(total_cells)

# IMPORTANT: remember that 14 in the Seurat object is patient 107.
GSE215968_raw_107_post <- subset(seurat_GSE215968_raw, subset = Patient == "01-07" & Treatment_stage == "B-Post") # we do not separate E2 from S2
GSE215968_raw_107_post <- subset(GSE215968_raw_107_post, cells = total_cells) # take post-107 cells
table(total_cells %in% colnames(GSE215968_raw_107_post)) # all true

dimplot_107_post <- DimPlot(GSE215968_raw_107_post,
                    reduction = "umap",
                    group.by = "cell_type_def") +
                    labs(title = "N07",
                    color = "", # Legend title
                    x = "",
                    y = "")

ggsave(paste0("/mnt/zonahpc/home/bioinformatica/servicios/20240724_UB375_MAESTER_Asherman/plots/umaps_", patient,"/umap_107_cell_type.png"), plot = dimplot_107_post, width = 11, height = 9, dpi = 300)

head(colnames(GSE215968_raw_107_post))
tail(colnames(GSE215968_raw_107_post))

merged_samples_variants[is.na(merged_samples_variants)] <- 0
head(merged_samples_variants)
tail(merged_samples_variants)

# UMAP
colnames(GSE215968_raw_107_post@meta.data)

trazado_107_E2 <- read.csv("/mnt/zonahpc/home/bioinformatica/servicios/20240724_UB375_MAESTER_Asherman/trazado_tipos_celulares/conteo_cell_types_107_E2_HQcb.csv", sep = ";")
trazado_107_S2 <- read.csv("/mnt/zonahpc/home/bioinformatica/servicios/20240724_UB375_MAESTER_Asherman/trazado_tipos_celulares/conteo_cell_types_107_S2_HQcb.csv", sep = ";")

trazado_107_E2[is.na(trazado_107_E2)] <- 0 # Last column has NAs for some reason
trazado_107_S2[is.na(trazado_107_S2)] <- 0 # Last column has NAs for some reason

head(trazado_107_E2)
head(trazado_107_S2)

trazado_107_E2_renamed <- trazado_107_E2 %>%
    select(starts_with("HQcb")) %>%
    rename_with(~ gsub("^HQcb_", "", .x)) %>%
    rename_with(~ paste0(.x, "_E2"))
trazado_107_S2_renamed <- trazado_107_S2 %>%
    select(starts_with("HQcb")) %>%
    rename_with(~ gsub("^HQcb_", "", .x)) %>%
    rename_with(~ paste0(.x, "_S2"))

trazado_107_E2_processed <- trazado_107_E2_renamed %>%
    summarise(across(everything(), ~ paste(.x[.x != 0], collapse = ",")))
trazado_107_S2_processed <- trazado_107_S2_renamed %>%
    summarise(across(everything(), ~ paste(.x[.x != 0], collapse = ",")))

# HARDCODED, CHECK MANUALLY: Remove those variants that we have seen do not trace CD133 cells
ncol(trazado_107_E2_processed) # 7
trazado_107_E2_processed <- trazado_107_E2_processed %>%
    select(-c("2623A.G_E2", "3085A.G_E2"))
ncol(trazado_107_E2_processed) # 5

ncol(trazado_107_S2_processed) # 3
trazado_107_S2_processed <- trazado_107_S2_processed %>%
    select(c("12416A.G_S2"))
ncol(trazado_107_S2_processed) # 1

# Combine both samples into a single dataframe with HQcb
trazado_107_all <- cbind(trazado_107_E2_processed, trazado_107_S2_processed)
ncol(trazado_107_all) # 6

# Convert to a more usable list
barcode_lists_107 <- lapply(trazado_107_all, function(x) unlist(strsplit(x, ",")))
barcode_lists_107 <- lapply(barcode_lists_107, trimws) # Remove whitespace

cell_barcodes <- rownames(GSE215968_raw_107_post@meta.data)

# Create the highlight column
GSE215968_raw_107_post@meta.data$highlight <- sapply(cell_barcodes, function(cell) {
    # Find which variants the cell belongs to
    variants <- names(barcode_lists_107)[sapply(barcode_lists_107, function(v) cell %in% v)]

    # If it belongs to at least one variant, concatenate the names; if not, assign "0"
    if (length(variants) > 0) {
        return(paste(variants, collapse = ","))
    } else {
        return("0")
    }
})

# At this point we already have the labels on the cells to plot, with "_" to indicate which cell they come from
table(GSE215968_raw_107_post@meta.data$highlight)

GSE215968_raw_107_post@meta.data$highlight <- GSE215968_raw_107_post@meta.data$highlight %>%
    str_replace("_.*", "") %>%
    str_replace_all("\\.", ">")

table(GSE215968_raw_107_post@meta.data$highlight)

# Get the unique levels of "highlight"
highlight_levels <- unique(GSE215968_raw_107_post@meta.data$highlight)

# Generate colors: "0" in grey, the rest with an extended palette
highlight_colors <- setNames(
    c("grey", viridis(length(highlight_levels) - 1)),
    c("0", highlight_levels[highlight_levels != "0"])
)

# Create the UMAP
umap_custom <- DimPlot(
    GSE215968_raw_107_post,
    reduction = "umap",
    group.by = "highlight",
    cols = highlight_colors,
    pt.size = 0.5
) +
    labs(
        title = "N07", # Change the title
        color = "",    # Change the legend title
        x = "",                # Remove X axis title
        y = ""                 # Remove Y axis title
    ) +
    guides(
        color = guide_legend(override.aes = list(size = 3))  # Adjust point size in the legend
    ) +
    scale_color_manual(
        values = highlight_colors,
        breaks = names(highlight_colors)[names(highlight_colors) != "0"]  # Exclude "0" from the legend
    ) +
    geom_point(
        data = GSE215968_raw_107_post@meta.data %>%
            mutate(UMAP_1 = Embeddings(GSE215968_raw_107_post, "umap")[, 1],
                   UMAP_2 = Embeddings(GSE215968_raw_107_post, "umap")[, 2]) %>%
            filter(highlight != "0"),  # Filter points to highlight
        aes(x = UMAP_1, y = UMAP_2, color = highlight), # Respect the "highlight" colors
        size = 1.25,                                    # Increased size for highlighted points
        inherit.aes = FALSE                             # Avoid conflicts with DimPlot mapping
    )

# Show the custom plot
# print(umap_custom)

ggsave(paste0("/mnt/zonahpc/home/bioinformatica/servicios/20240724_UB375_MAESTER_Asherman/plots/umaps_", patient,"/umap_107_variants.png"), plot = umap_custom, width = 11, height = 9, dpi = 300)

# SEURAT ALL #########################################################################################################
# Cell type
nrow(info_var_115_E2) # 12415
nrow(info_var_115_S2) # 12167
# rownames(info_var_115_E2) <- paste0("ENTIRE-1-14-E2_", rownames(info_var_115_E2))
# rownames(info_var_115_S2) <- paste0("ENTIRE-1-14-S2_", rownames(info_var_115_S2))

nrow(info_var_107_E2) # 3607
nrow(info_var_107_S2) # 2944
# rownames(info_var_107_E2) <- paste0("ENTIRE-1-07-E2_", rownames(info_var_107_E2))
# rownames(info_var_107_S2) <- paste0("ENTIRE-1-07-S2_", rownames(info_var_107_S2))

# nrow(info_var_109_E2) # 12415
nrow(info_var_109_S2) # 2331
# rownames(info_var_109_E2) <- paste0("ENTIRE-1-09-E2_", rownames(info_var_109_E2))
# rownames(info_var_109_S2) <- paste0("ENTIRE-1-09-S2_", rownames(info_var_109_S2))

merged_all_cells <- bind_rows(info_var_115_E2, info_var_115_S2, info_var_107_E2, info_var_107_S2, info_var_109_S2)
all_total_cells <- rownames(merged_all_cells)

length(all_total_cells) # 33464 cells

all(all_total_cells %in% colnames(seurat_GSE215968_raw))
GSE215968_selected_all_cells <- subset(seurat_GSE215968_raw, cells = all_total_cells)

table(GSE215968_selected_all_cells$Cell_type, GSE215968_selected_all_cells$Patient, GSE215968_selected_all_cells$Treatment_stage)

dimplot_all_post <- DimPlot(GSE215968_selected_all_cells,
                    reduction = "umap",
                    group.by = "cell_type_def") +
                    labs(
                        title = "",
                        color = "", # Legend title
                        x = "UMAP_1",
                        y = "UMAP_2"
                        )

dimplot_all_post_label <- LabelClusters(
  dimplot_all_post,
  id = "cell_type_def",
  repel = TRUE,
  size = 5,
  clusters = setdiff(unique(GSE215968_selected_all_cells$cell_type_def), "25") # Exclude the missing cluster
)

ggsave(paste0("/mnt/zonahpc/home/bioinformatica/servicios/20240724_UB375_MAESTER_Asherman/plots/umap_all_cell_type_label.png"), plot = dimplot_all_post_label, width = 11, height = 9, dpi = 300)

# ggsave(paste0("/mnt/zonahpc/home/bioinformatica/servicios/20240724_UB375_MAESTER_Asherman/plots/umap_all_cell_type_label.svg"), plot = dimplot_all_post_label, width = 11, height = 9, dpi = 300) # error due to missing png.h library

svg(filename = "/mnt/zonahpc/home/bioinformatica/servicios/20240724_UB375_MAESTER_Asherman/plots/umap_all_cell_type_label.svg",
    width = 11,
    height = 9)
print(dimplot_all_post_label)
dev.off()

# Variants
barcode_lists_107
barcode_lists_109
barcode_lists_115

barcode_lists_all <- c(barcode_lists_107, barcode_lists_109, barcode_lists_115)

cell_barcodes <- rownames(GSE215968_selected_all_cells@meta.data)

# We create the highlight column
GSE215968_selected_all_cells@meta.data$highlight <- sapply(cell_barcodes, function(cell) {
    # Search for which variants the cell belongs to
    variants <- names(barcode_lists_all)[sapply(barcode_lists_all, function(v) cell %in% v)]

    # If it belongs to at least one variant, concatenate the names; if not, assign "0"
    if (length(variants) > 0) {
        return(paste(variants, collapse = ","))
    } else {
        return("0")
    }
})

# At this point we already have the labels on the cells to be plotted, with "_" to indicate which cell they come from
table(GSE215968_selected_all_cells@meta.data$highlight)

table(GSE215968_selected_all_cells$cell_type_def)

GSE215968_selected_all_cells@meta.data$highlight <- GSE215968_selected_all_cells@meta.data$highlight %>%
    str_replace("_.*", "") %>%
    str_replace_all("\\.", ">")

table(GSE215968_selected_all_cells@meta.data$highlight)

# Get the unique levels of "highlight"
highlight_levels <- unique(GSE215968_selected_all_cells@meta.data$highlight)

# Generate colors: "0" in light grey, the rest with an extended palette
highlight_colors <- setNames(
    c("lightgrey", viridis(length(highlight_levels) - 1)),
    c("0", highlight_levels[highlight_levels != "0"])
)

# Create the UMAP
# group.by = "cell_type_def"
umap_custom <- DimPlot(
    GSE215968_selected_all_cells,
    reduction = "umap",
    group.by = "cell_type_def",
    cols = highlight_colors,
    pt.size = 0.5
) +
    labs(
        title = "", # Change the title
        color = "", # Change the legend title
        x = "UMAP_1",     # Remove X axis title
        y = "UMAP_2"      # Remove Y axis title
    ) +
    guides(
        color = guide_legend(override.aes = list(size = 3))  # Adjust point size in the legend
    ) +
    scale_color_manual(
        values = highlight_colors,
        breaks = names(highlight_colors)[names(highlight_colors) != "0"]  # Exclude "0" from the legend
    ) +
    geom_point(
        data = GSE215968_selected_all_cells@meta.data %>%
            mutate(UMAP_1 = Embeddings(GSE215968_selected_all_cells, "umap")[, 1],
                   UMAP_2 = Embeddings(GSE215968_selected_all_cells, "umap")[, 2]) %>%
            filter(highlight != "0"),  # Filter points to highlight
        aes(x = UMAP_1, y = UMAP_2, color = highlight), # Respect the "highlight" colors
        size = 1.25,                                    # Increased size for highlighted points
        inherit.aes = FALSE                             # Avoid conflicts with DimPlot mapping
    )

umap_variants_label <- LabelClusters(
  umap_custom,
  id = "cell_type_def",
  repel = TRUE,
  size = 5,
  clusters = setdiff(unique(GSE215968_selected_all_cells$cell_type_def), "25") # Exclude the missing cluster
)

# Show the custom plot
# print(umap_custom)

ggsave(paste0("/mnt/zonahpc/home/bioinformatica/servicios/20240724_UB375_MAESTER_Asherman/plots/umap_all_variants_label.png"), plot = umap_variants_label, width = 11, height = 9, dpi = 300)

svg(filename = "/mnt/zonahpc/home/bioinformatica/servicios/20240724_UB375_MAESTER_Asherman/plots/umap_all_variants_label.svg",
    width = 11,
    height = 9)
print(umap_variants_label)
dev.off()

## PLOTS ###############################################################################################################
# Made after what was the final submission
datos_seurat_path <- "/mnt/zonahpc/home/bioinformatica/servicios/20240724_UB375_MAESTER_Asherman/GSE215968/AS_pre_post.H5Seurat"
seurat_GSE215968_raw <- LoadH5Seurat(datos_seurat_path) # load the data into a Seurat object, this takes a while

table(seurat_GSE215968_raw$cell_type_def) # AS_Epithelium, B_cells, CD8T_cells, Ciliated_Epithelium...
table(seurat_GSE215968_raw$Patient) # 01-05 01-07 01-08 01-09 01-10 01-12 01-13 01-14 01-15
table(seurat_GSE215968_raw$Cell_type) # Epithelium     Stroma, # change this name to biopsy_origin
table(seurat_GSE215968_raw$Treatment_stage) # A-Pre, B-Post
table(seurat_GSE215968_raw$cell_type_def)

table(seurat_GSE215968_raw$Cell_type, seurat_GSE215968_raw$Patient, seurat_GSE215968_raw$Treatment_stage) # check table, is the number of cells as expected? Remember that our 115 is actually 114 in the object.

umap <- DimPlot(
  seurat_GSE215968_raw,
  reduction = "umap",
  group.by = "cell_type_def",
  split.by = "Treatment_stage",
  raster = FALSE,
) +
  labs(
    title = "cell_type_def",
    color = "", # Legend title
    x = "UMAP_1",
    y = "UMAP_2"
  ) +
  theme(
    strip.text = element_text(size = 16, face = "bold") # Size and style of facet titles
  )

# umap_label <- LabelClusters(
#   umap,
#   id = "cell_type_def",
#   repel = TRUE, # To avoid overlapping labels
#   size = 6,25
# )

umap_label <- LabelClusters(
  umap,
  id = "cell_type_def",
  repel = TRUE,
  size = 5,
  clusters = setdiff(unique(seurat_GSE215968_raw$cell_type_def), "25") # Exclude the missing cluster
)

umap_label <- umap_label +
  theme(
    plot.title = element_text(size = 18, face = "bold"), # Size and style of the main title
    axis.title.x = element_text(size = 14),             # Size of the X axis title
    axis.title.y = element_text(size = 14),             # Size of the Y axis title
    axis.text = element_text(size = 12)                 # Size of the axis text
  )

ggsave(paste0("/mnt/zonahpc/home/bioinformatica/servicios/20240724_UB375_MAESTER_Asherman/plots/umap_pre_vs_post_label.png"), plot = umap_label, width = 25, height = 15, dpi = 300)

svg(filename = "/mnt/zonahpc/home/bioinformatica/servicios/20240724_UB375_MAESTER_Asherman/plots/umap_pre_vs_post_label.svg",
    width = 25,
    height = 15)
print(umap_label)
dev.off()

## 03/03/2025: Check SOX9 gene expression ############################################################################
barcode_lists_todos # These are the cells we have been able to trace with maester
GSE215968_selected_all_cells # this Seurat object contains all cells for the selected patients and conditions

table(GSE215968_selected_all_cells$Cell_type, GSE215968_selected_all_cells$Patient, GSE215968_selected_all_cells$Treatment_stage)

# Merge all cells from the different groups into a single vector
todos_los_cells <- unique(unlist(barcode_lists_todos))
# From the Seurat object, keep the cells we have been able to trace:
subset_obj <- subset(GSE215968_selected_all_cells, cells = todos_los_cells)

# Untraced cells
not_trazadas_cells <- setdiff(colnames(GSE215968_selected_all_cells), todos_los_cells)
not_trazadas_obj <- subset(GSE215968_selected_all_cells, cells = not_trazadas_cells)

###### Search in the Seurat object
"SOX9" %in% rownames(subset_obj) # TRUE, the gene of interest is present
"SOX9" %in% rownames(subset_obj@assays$RNA@data) # TRUE
slotNames(subset_obj)
Assays(subset_obj)

sox9_umap_all <- FeaturePlot(GSE215968_selected_all_cells, features = "SOX9")
sox9_umap_trazadas <- FeaturePlot(subset_obj, features = "SOX9")
sox9_vln_all <- VlnPlot(GSE215968_selected_all_cells, features = "SOX9")
sox9_vln_trazadas <- VlnPlot(subset_obj, features = "SOX9")

ggsave(paste0("/mnt/zonahpc/home/bioinformatica/servicios/20240724_UB375_MAESTER_Asherman/plots/sox9_umap_all.png"), plot = sox9_umap_all, width = 25, height = 15, dpi = 300)
ggsave(paste0("/mnt/zonahpc/home/bioinformatica/servicios/20240724_UB375_MAESTER_Asherman/plots/sox9_umap_trazadas.png"), plot = sox9_umap_trazadas, width = 25, height = 15, dpi = 300)
ggsave(paste0("/mnt/zonahpc/home/bioinformatica/servicios/20240724_UB375_MAESTER_Asherman/plots/sox9_vln_all.png"), plot = sox9_vln_all, width = 25, height = 15, dpi = 300)
ggsave(paste0("/mnt/zonahpc/home/bioinformatica/servicios/20240724_UB375_MAESTER_Asherman/plots/sox9_vln_trazadas.png"), plot = sox9_vln_trazadas, width = 25, height = 15, dpi = 300)

# Checks to make sure everything is ok
actb_umap_all <- FeaturePlot(GSE215968_selected_all_cells, features = "ACTB")
actb_umap_trazadas <- FeaturePlot(subset_obj, features = "ACTB")
actb_vln_all <- VlnPlot(GSE215968_selected_all_cells, features = "ACTB")
actb_vln_trazadas <- VlnPlot(subset_obj, features = "ACTB")

ggsave(paste0("/mnt/zonahpc/home/bioinformatica/servicios/20240724_UB375_MAESTER_Asherman/plots/actb_umap_all.png"), plot = actb_umap_all, width = 25, height = 15, dpi = 300)
ggsave(paste0("/mnt/zonahpc/home/bioinformatica/servicios/20240724_UB375_MAESTER_Asherman/plots/actb_umap_trazadas.png"), plot = actb_umap_trazadas, width = 25, height = 15, dpi = 300)
ggsave(paste0("/mnt/zonahpc/home/bioinformatica/servicios/20240724_UB375_MAESTER_Asherman/plots/actb_vln_all.png"), plot = actb_vln_all, width = 25, height = 15, dpi = 300)
ggsave(paste0("/mnt/zonahpc/home/bioinformatica/servicios/20240724_UB375_MAESTER_Asherman/plots/actb_vln_trazadas.png"), plot = actb_vln_trazadas, width = 25, height = 15, dpi = 300)
#

subset_obj@assays$RNA@data %>% as.data.frame() %>% filter(rowname == "SOX9")

sox9_expression <- subset_obj@assays$RNA@data["SOX9", , drop = FALSE]
sox9_count <- subset_obj@assays$RNA@counts["SOX9", , drop = FALSE]

summary(as.data.frame(sox9_expression))
summary(as.data.frame(sox9_count))
