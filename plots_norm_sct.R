#!/usr/bin/env Rscript

# Load required packages
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

# Get input file from command line
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript plots_norm_sct.R <path_to_seurat_list_sctransform_normalized.rds>")
}
input_path <- args[1]
stopifnot(file.exists(input_path))

# Load Seurat list
seurat_list <- readRDS(input_path)
message("Loaded Seurat list with ", length(seurat_list), " samples.")

# Initialize storage
all_metadata <- data.frame()

# Loop through all samples
for (sample_name in names(seurat_list)) {
  obj <- seurat_list[[sample_name]]
  
  # Extract data
  raw_total <- Matrix::colSums(GetAssayData(obj, assay = "RNA", slot = "counts"))
  norm_total <- Matrix::colSums(GetAssayData(obj, assay = "SCT", slot = "data"))  # Already log1p-normalized
  
  # Build table
  df <- data.frame(
    sample = sample_name,
    raw_total = raw_total,
    norm_total = norm_total,
    cell = names(raw_total)
  )
  
  all_metadata <- rbind(all_metadata, df)
}

# PLOT 1: Histograms before vs after normalization (log1p scale)
all_metadata$log1p_raw <- log1p(all_metadata$raw_total)
all_metadata$log1p_norm <- all_metadata$norm_total  # Already log1p from SCTransform

p1 <- ggplot(all_metadata, aes(x = log1p_raw)) +
  geom_histogram(fill = "steelblue", bins = 60) +
  labs(title = "Raw counts (log1p)", x = "log1p(Raw counts)", y = "Frequency") +
  theme_minimal()

p2 <- ggplot(all_metadata, aes(x = log1p_norm)) +
  geom_histogram(fill = "orchid", bins = 60) +
  labs(title = "SCTransform normalized counts (log1p)", x = "SCT normalized counts", y = "Frequency") +
  theme_minimal()

ggsave("sct_normalization_log1p_histograms.png", plot = p1 + p2, width = 10, height = 5)
message("Saved: sct_normalization_log1p_histograms.png")

# PLOT 2: Histograms in raw scale (no log1p)
real_raw <- all_metadata$raw_total
real_norm <- expm1(all_metadata$norm_total)

p_real_raw <- ggplot(data.frame(real_raw), aes(x = real_raw)) +
  geom_histogram(fill = "steelblue", bins = 60) +
  labs(title = "Total raw counts", x = "Raw counts", y = "Frequency") +
  scale_x_continuous(labels = comma_format()) +
  theme_minimal()

p_real_norm <- ggplot(data.frame(real_norm), aes(x = real_norm)) +
  geom_histogram(fill = "orchid", bins = 60) +
  labs(title = "SCTransform normalized counts", x = "Normalized counts", y = "Frequency") +
  scale_x_continuous(labels = comma_format()) +
  theme_minimal()

ggsave("sct_histograms_real_counts.png", plot = p_real_raw + p_real_norm, width = 10, height = 5)
message("Saved: sct_histograms_real_counts.png")
