#!/usr/bin/env Rscript

set.seed(123)

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript normalize_scran.R <path_to_input.rds>\n",
       "Example: Rscript normalize_scran.R seurat_list_filtered.rds")
}
input_path <- args[1] # error if its not the argument 
stopifnot(file.exists(input_path))

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(SingleCellExperiment)
  library(scran)
})

message("# Loaded packages")

# Load input object
x <- readRDS(input_path)
message("# Loaded input: ", input_path)
message("# Class: ", paste(class(x), collapse = ", "))

# Normalize to list
if (inherits(x, "Seurat")) {
  seurat_list <- list(x)  # detet if its a seurat objet or a list of obj?
  nm <- if ("sample_id" %in% colnames(x@meta.data)) unique(x$sample_id)[1] else "sample1" # if it is just 1 objet -> convert to a list of a single obj 
  names(seurat_list) <- nm
} else if (is.list(x) && all(vapply(x, function(z) inherits(z, "Seurat"), logical(1)))) {
  seurat_list <- x
  if (is.null(names(seurat_list))) names(seurat_list) <- paste0("sample", seq_along(seurat_list))  # list of objets and name 
} else {
  stop("Input must be a Seurat object or list of Seurat objects.")
}

message("# Samples: ", length(seurat_list))
message("# Example sample names: ", paste(head(names(seurat_list)), collapse = ", "))

# Create output directory
stamp  <- format(Sys.time(), "%Y-%m-%d_%H-%M") # stamp
base   <- tools::file_path_sans_ext(basename(input_path)) # time 
outdir <- file.path("normalized_out", paste0(base, "-", stamp)) # save file in my stampp
dir.create(outdir, recursive = TRUE, showWarnings = FALSE) 
message("# Output directory: ", outdir)

# Helper to get counts safely (Seurat v5-compatible)

# funtion: Seurat v5 -> Obtein the count matriz thougt "layer "counts"" , 
get_counts_safe <- function(obj, assay = "RNA") {
  a <- obj[[assay]]
  if (!is.null(a@layers) && "counts" %in% names(a@layers)) {
    return(GetAssayData(obj, assay = assay, layer = "counts"))
  } else {  # if its not a seurat v5 -> thougt slot "counts"
    return(GetAssayData(obj, assay = assay, slot = "counts"))
  }
}

# Normalization function

normalize_with_scran <- function(obj, sample_name) {
  message("  >> Normalizing sample: ", sample_name)
  
  DefaultAssay(obj) <- "RNA" # ensured that it is RNA
  
  raw_counts <- tryCatch({  #o gets the count matrix 
    get_counts_safe(obj)
  }, error = function(e) {
    message("    ! Failed to extract counts: ", e$message)
    return(NULL)
  })
  
  if (is.null(raw_counts)) {
    message("    ! Skipping ", sample_name, ": no counts matrix.")
    return(NULL)
  }  
  
  # Convert to SingleCellExperiment
  sce <- tryCatch({
    SingleCellExperiment(assays = list(counts = raw_counts))
  }, error = function(e) {
    message("    ! Failed to convert to SCE: ", e$message)
    return(NULL)
  })
  if (is.null(sce)) return(NULL)
  
  # Compute size factors
  sce <- tryCatch({
    computeSumFactors(sce, min.mean = 0.1)
  }, error = function(e) {
    message("    ! Failed to compute size factors: ", e$message)
    return(NULL)
  })
  
  sf <- sizeFactors(sce)
  if (is.null(sf)) {
    message("    ! No size factors; skipping ", sample_name)
    return(NULL)
  }
  
  # Apply normalization
  norm_counts <- log1p(t(t(raw_counts) / sf))
  
  # Ensure rownames and colnames
  rownames(norm_counts) <- rownames(raw_counts)
  colnames(norm_counts) <- colnames(raw_counts)
  
  # Add normalized assay
  obj[["scran_norm"]] <- CreateAssayObject(data = norm_counts)
  DefaultAssay(obj) <- "scran_norm"
  obj$size_factors <- sf
  
  message("    ✓ Normalization complete for ", sample_name)
  return(obj)
}

# Run normalization
seurat_list_norm <- mapply(
  normalize_with_scran,
  seurat_list,
  names(seurat_list),
  SIMPLIFY = FALSE
)

# Drop failed samples
ok <- !vapply(seurat_list_norm, is.null, logical(1))
if (!all(ok)) {
  message("# Skipping ", sum(!ok), " failed samples.")
  seurat_list_norm <- seurat_list_norm[ok]
}

# Save output
output_file <- file.path(outdir, "seurat_list_scran_normalized.rds")
saveRDS(seurat_list_norm, file = output_file)
message("# DONE. Saved: ", output_file)

