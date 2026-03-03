#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(Matrix)
  library(future)
})

start_time <- Sys.time()

# 
option_list <- list(
  make_option(c("--seurat_rds"), type = "character", help = "Path to Seurat .rds", metavar = "FILE"),
  make_option(c("--out_dir"),    type = "character", default = "unassigned_only_final",
              help = "Output directory [default: %default]", metavar = "DIR"),
  make_option(c("--assay"),      type = "character", default = "RNA",
              help = "Assay to use if present (fallback to DefaultAssay if not) [default: %default]"),
  make_option(c("--layer_prefix"), type = "character", default = "data",
              help = "Layer family to join (e.g. 'data' joins data.1,data.2,...) [default: %default]"),
  make_option(c("--unassigned_mode"), type = "character", default = "auto",
              help = "How to define unassigned cells: auto | cell_type_na | cell_type_plot [default: %default]"),
  make_option(c("--unassigned_label"), type = "character", default = "Unassigned",
              help = "Label used in cell_type_plot for unassigned [default: %default]"),
  make_option(c("-w", "--workers"), type = "integer", default = 4,
              help = "Parallel workers (for future plan) [default: %default]"),
  make_option(c("--seed"), type = "integer", default = 42,
              help = "Random seed [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

stopifnot(!is.null(opt$seurat_rds))
stopifnot(file.exists(opt$seurat_rds))

set.seed(opt$seed)
plan(multicore, workers = opt$workers)

# ---- Output dir ----
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

cat("# Seurat rds: ", opt$seurat_rds, "\n")
cat("# Out dir:   ", normalizePath(opt$out_dir, winslash = "/", mustWork = FALSE), "\n")
cat("# Assay pref:", opt$assay, "\n")
cat("# Layer pref:", opt$layer_prefix, "\n")
cat("# Seed:      ", opt$seed, "\n")
cat("# Workers:   ", opt$workers, "\n\n")

#  Load object 
cat("# Loading Seurat object...\n")
obj <- readRDS(opt$seurat_rds)
cat("# Cells in Seurat: ", ncol(obj), "\n")

# 1) Select unassigned cells 
meta_cols <- colnames(obj@meta.data)

if (opt$unassigned_mode == "cell_type_na") {
  if (!("cell_type" %in% meta_cols)) stop("unassigned_mode=cell_type_na but 'cell_type' not found in meta.data")
  un_cells <- colnames(obj)[is.na(obj$cell_type)]
} else if (opt$unassigned_mode == "cell_type_plot") {
  if (!("cell_type_plot" %in% meta_cols)) stop("unassigned_mode=cell_type_plot but 'cell_type_plot' not found in meta.data")
  un_cells <- colnames(obj)[obj$cell_type_plot == opt$unassigned_label]
} else {
  # auto
  if ("cell_type" %in% meta_cols) {
    un_cells <- colnames(obj)[is.na(obj$cell_type)]
  } else if ("cell_type_plot" %in% meta_cols) {
    un_cells <- colnames(obj)[obj$cell_type_plot == opt$unassigned_label]
  } else {
    stop("Cannot find 'cell_type' or 'cell_type_plot' in meta.data")
  }
}

cat("# Unassigned cells: ", length(un_cells), "\n")
if (length(un_cells) == 0) {
  stop("No unassigned cells found. Nothing to export.")
}

obj_un <- subset(obj, cells = un_cells)

# 2) Use correct assay 
assay_use <- if (opt$assay %in% Assays(obj_un)) opt$assay else DefaultAssay(obj_un)
DefaultAssay(obj_un) <- assay_use

cat("# Using assay: ", assay_use, "\n")
cat("# Layers BEFORE join:\n")
print(Layers(obj_un[[assay_use]]))

# 3) Join layers (e.g., data.* -> data) 
cat("# Joining layers: ", opt$layer_prefix, ".* -> ", opt$layer_prefix, "\n", sep = "")
obj_un[[assay_use]] <- JoinLayers(obj_un[[assay_use]], layers = opt$layer_prefix)

cat("# Layers AFTER join:\n")
print(Layers(obj_un[[assay_use]]))

# 4) Extract joined layer and export 
cat("# Extracting layer: ", opt$layer_prefix, "\n", sep = "")
mat <- LayerData(obj_un[[assay_use]], layer = opt$layer_prefix)

cat("# Matrix dim (genes x cells): ", paste(dim(mat), collapse = " x "), "\n", sep = "")

mtx_path   <- file.path(opt$out_dir, "unassigned_logNorm.mtx")
genes_path <- file.path(opt$out_dir, "genes.txt")
cells_path <- file.path(opt$out_dir, "cells.txt")

Matrix::writeMM(mat, mtx_path)
writeLines(rownames(mat), genes_path)
writeLines(colnames(mat), cells_path)

cat("\n# Wrote:\n")
cat(" - ", normalizePath(mtx_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat(" - ", normalizePath(genes_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat(" - ", normalizePath(cells_path, winslash = "/", mustWork = FALSE), "\n", sep = "")

cat("# Files exist? ",
    file.exists(mtx_path), file.exists(genes_path), file.exists(cells_path), "\n")

end_time <- Sys.time()
cat("\n# DONE. Time taken: ", end_time - start_time, "\n")
