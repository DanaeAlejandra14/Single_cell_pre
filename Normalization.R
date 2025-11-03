#Normalization

# Input: seurat_list_filtered : Just the qc 
# Output : seurat_list_normalized : Before the normalization

#!/usr/bin/env Rscript

set.seed(123)

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript normalize_scran.R <path_to_input.rds>\n",
       "Example: Rscript normalize_scran.R seurat_list_filtered.rds")
}
input_path <- args[1]
stopifnot(file.exists(input_path))

# Load required packages
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(SingleCellExperiment)
  library(scran)
  library(BiocParallel)
})

message("# Loaded packages")

# Load input Seurat object or list of Seurat objects
x <- readRDS(input_path)
message("# Loaded input: ", input_path)
message("# Class: ", paste(class(x), collapse = ", "))

# Normalize to a named list
if (inherits(x, "Seurat")) {
  seurat_list <- list(x)
  nm <- if ("sample_id" %in% colnames(x@meta.data)) unique(x$sample_id)[1] else "sample1"
  names(seurat_list) <- nm
} else if (is.list(x) && all(vapply(x, \(z) inherits(z, "Seurat"), logical(1)))) {
  seurat_list <- x
  if (is.null(names(seurat_list))) names(seurat_list) <- paste0("sample", seq_along(seurat_list))
} else {
  stop("Input must be a Seurat object or list of Seurat objects.")
}

message("# Samples: ", length(seurat_list))
message("# Example sample names: ", paste(head(names(seurat_list)), collapse = ", "))

# Create output directory
stamp  <- format(Sys.time(), "%Y-%m-%d_%H-%M")
base   <- tools::file_path_sans_ext(basename(input_path))
outdir <- file.path("normalized_out", paste0(base, "-", stamp))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
message("# Output directory: ", outdir)

# Function to normalize a Seurat object using scran
normalize_with_scran <- function(obj, sample_name) {
  message("  >> Normalizing sample: ", sample_name)
  
  # Convert to SingleCellExperiment
  sce <- tryCatch({
    as.SingleCellExperiment(obj, assay = "RNA")
  }, error = function(e) {
    message("    ! Error converting to SCE: ", e$message)
    return(NULL)
  })
  
  # Compute size factors using scran
  sce <- tryCatch({
    computeSumFactors(sce, min.mean = 0.1, BPPARAM = MulticoreParam())
  }, error = function(e) {
    message("    ! Error in computeSumFactors: ", e$message)
    return(NULL)
  })
  
  # Extract size factors
  size_factors <- sizeFactors(sce)
  obj$size_factors <- size_factors
  
  # Normalize counts and apply log1p transformation
  norm_counts <- log1p(t(t(counts(sce)) / size_factors))
  
  # Store normalized data in a new assay
  obj[["scran_norm"]] <- CreateAssayObject(data = norm_counts)
  DefaultAssay(obj) <- "scran_norm"
  
  return(obj)
}

# Apply normalization to each sample in the list
seurat_list_norm <- mapply(normalize_with_scran, seurat_list, names(seurat_list),
                           SIMPLIFY = FALSE)

# Save the output
saveRDS(seurat_list_norm, file = file.path(outdir, "seurat_list_scran_normalized.rds"))
message("# Saved: seurat_list_scran_normalized.rds")
