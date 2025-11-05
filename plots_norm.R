#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
input_path <- args[1]

library(Seurat)
library(scran)
library(Matrix)
library(ggplot2)

# Load normalized object
seurat_list <- readRDS(input_path)
sample_name <- names(seurat_list)[1]  # Just an example sample
obj <- seurat_list[[sample_name]]

# Extract raw counts and size factors 
raw_counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
scran_sf <- obj$size_factors

# Compute library size factor (total counts per cell)
lib_sf <- Matrix::colSums(raw_counts)
lib_sf <- lib_sf / mean(lib_sf)  # Normalize to mean 1

#  Scatter plot: Deconvolution vs Library size factors 
png("scran_vs_library_sf.png", width = 800, height = 600)
plot(lib_sf, scran_sf,
     xlab = "Library size factor",
     ylab = "Deconvolution size factor",
     log = "xy", pch = 16, col = "steelblue")
abline(a = 0, b = 1, col = "red", lwd = 2)
dev.off()

#  Histograms 
total_counts <- Matrix::colSums(raw_counts)
scran_counts <- Matrix::colSums(GetAssayData(obj, assay = "scran_norm", slot = "data"))

png("scran_histograms.png", width = 1000, height = 500)
par(mfrow = c(1,2))
hist(total_counts, breaks = 100, main = "Total counts", xlab = "total_counts", col = "lightblue")
hist(scran_counts, breaks = 100, main = "log1p with Scran estimated size factors", col = "lightgreen")
dev.off()

