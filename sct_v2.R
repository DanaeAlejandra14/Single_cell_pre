#!/usr/bin/env Rscript

set.seed(123)

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript normalize_sctransform.R <path_to_input.rds>\n",
       "Example: Rscript normalize_sctransform.R seurat_list_filtered.rds")
}
input_path <- args[1]
stopifnot(file.exists(input_path))

# Load required packages
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
})

message("# Loaded packages")

# Load input object
x <- readRDS(input_path)
message("# Loaded input: ", input_path)
message("# Class: ", paste(class(x), collapse = ", "))

# Normalize to list
if (inherits(x, "Seurat")) {
  seurat_list <- list(x)
  nm <- if ("sample_id" %in% colnames(x@meta.data)) unique(x$sample_id)[1] else "sample1"
  names(seurat_list) <- nm
} else if (is.list(x) && all(vapply(x, function(z) inherits(z, "Seurat"), logical(1)))) {
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

# Normalization with SCTransform 

normalize_with_sct <- function(obj, sample_name) {
  message("  >> Normalizing sample: ", sample_name)
  
  DefaultAssay(obj) <- "RNA"
  
  tryCatch({
      obj <- SCTransform(
        obj,
        assay = "RNA",
        new.assay.name = "SCT",
        vars.to.regress = c("percent.mt", "nCount_RNA"),
        vst.flavor = "v2",
        return.only.var.genes = FALSE,
        verbose = FALSE
        )
    message("    ✓ Normalization complete for ", sample_name)
    return(obj)
  }, error = function(e) {
    message("    ! Failed SCTransform on ", sample_name, ": ", e$message)
    return(NULL)
  })
}

# Apply normalization
seurat_list_norm <- mapply(
  normalize_with_sct,
  seurat_list,
  names(seurat_list),
  SIMPLIFY = FALSE
)

# Drop failed
ok <- !vapply(seurat_list_norm, is.null, logical(1))
if (!all(ok)) {
  message("# Skipping ", sum(!ok), " failed samples.")
  seurat_list_norm <- seurat_list_norm[ok]
}

# Save output
output_file <- file.path(outdir, "seurat_list_sctransform_normalized.rds")
saveRDS(seurat_list_norm, file = output_file)
message("# DONE. Saved: ", output_file)