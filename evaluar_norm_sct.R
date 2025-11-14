#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Provide path to RDS file.")

input_path <- args[1]
x <- readRDS(input_path)

# Asegura que sea lista
if (inherits(x, "Seurat")) {
  seurat_list <- list(x)
  names(seurat_list) <- if ("sample_id" %in% colnames(x@meta.data)) unique(x$sample_id)[1] else "sample1"
} else if (is.list(x) && all(vapply(x, \(o) inherits(o, "Seurat"), logical(1)))) {
  seurat_list <- x
  if (is.null(names(seurat_list))) names(seurat_list) <- paste0("sample", seq_along(seurat_list))
} else {
  stop("Input must be a Seurat object or list of Seurat objects.")
}

# Evaluar objetos
results <- lapply(names(seurat_list), function(nm) {
  obj <- seurat_list[[nm]]
  sct_exists <- "SCT" %in% names(obj@assays)
  cells <- colnames(obj)
  has_dup_cells <- any(duplicated(cells))
  
  if (sct_exists) {
    assay <- obj[["SCT"]]
    n_cells_data <- ncol(GetAssayData(obj, assay = "SCT", layer = "data"))
    n_cells_metadata <- ncol(obj)
    list(
      sample = nm,
      has_SCT = TRUE,
      cells_in_SCT = n_cells_data,
      cells_meta = n_cells_metadata,
      n_genes = nrow(GetAssayData(obj, assay = "SCT", layer = "data")),
      duplicate_cells_within = has_dup_cells
    )
  } else {
    list(
      sample = nm,
      has_SCT = FALSE,
      cells_in_SCT = NA,
      cells_meta = length(cells),
      n_genes = NA,
      duplicate_cells_within = has_dup_cells
    )
  }
})

df <- bind_rows(results)

# Ver duplicados entre objetos
all_cells <- unlist(lapply(seurat_list, colnames))
duplicates_between <- any(duplicated(all_cells))

cat("# Evaluación completa de objetos:\n")
print(df, n = Inf)

cat("\n# ¿Hay duplicados entre objetos?:", duplicates_between, "\n")

