#!/usr/bin/env Rscript

# Load required packages
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
})

# Get input file from command line
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript plots_norm.R <path_to_seurat_list_scran_normalized.rds>")
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
  
  # Skip if missing size_factors
  if (is.null(obj$size_factors)) next
  
  # Extract data
  raw_total <- Matrix::colSums(GetAssayData(obj, assay = "RNA", slot = "counts"))
  norm_total <- Matrix::colSums(GetAssayData(obj, assay = "scran_norm", slot = "data"))
  size_factors <- obj$size_factors
  
  # Build table
  df <- data.frame(
    sample = sample_name,
    raw_total = raw_total,
    norm_total = norm_total,
    size_factor = size_factors,
    cell = names(raw_total)
  )
  
  all_metadata <- rbind(all_metadata, df)
}

# ---- PLOT 1: Histograms before vs after normalization ----
all_metadata$log1p_raw <- log1p(all_metadata$raw_total)
all_metadata$log1p_norm <- log1p(all_metadata$norm_total)

p1 <- ggplot(all_metadata, aes(x = log1p_raw)) +
  geom_histogram(fill = "steelblue", bins = 60) +
  labs(title = "Total raw counts (log1p)", x = "log1p(Raw counts)", y = "Frequency") +
  theme_minimal()

p2 <- ggplot(all_metadata, aes(x = log1p_norm)) +
  geom_histogram(fill = "orchid", bins = 60) +
  labs(title = "Scran normalized counts (log1p)", x = "log1p(Scran counts)", y = "Frequency") +
  theme_minimal()

# Save histogram comparison
ggsave("scran_normalization_global_histograms.png", plot = p1 + p2, width = 10, height = 5)
message("Saved: scran_normalization_global_histograms.png")

# ---- PLOT 2: Library size factor vs scran size factor ----
lib_size_factor <- all_metadata$raw_total / mean(all_metadata$raw_total)

p3 <- ggplot(all_metadata, aes(x = lib_size_factor, y = size_factor)) +
  geom_point(alpha = 0.4, size = 1) +
  scale_x_log10() + scale_y_log10() +
  labs(
    title = "Library vs Scran size factor",
    x = "Library size factor (log10)",
    y = "Scran size factor (log10)"
  ) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  theme_minimal()

# Save scatter plot
ggsave("scran_size_factor_comparison.png", plot = p3, width = 6, height = 5)
message("Saved: scran_size_factor_comparison.png")
