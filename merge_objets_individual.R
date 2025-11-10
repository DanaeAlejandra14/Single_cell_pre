#!/usr/bin/env Rscript

# Load required packages
suppressPackageStartupMessages({
  library(Seurat)
})

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript merge_objects_individual.R <seurat_list.rds> <biospecimen_metadata.csv>")
}

input_seurat <- args[1]
input_metadata <- args[2]

stopifnot(file.exists(input_seurat))
stopifnot(file.exists(input_metadata))

# Load Seurat list and metadata
seurat_list_filtered <- readRDS(input_seurat)
bio_meta <- read.csv(input_metadata)

cat("# Loaded Seurat objects:", length(seurat_list_filtered), "\n")
cat("# Loaded metadata rows:", nrow(bio_meta), "\n")

# Create correspondence table between sample and individual
id_map <- unique(bio_meta[, c("specimenID", "individualID")])

# Add individualID to each Seurat object
for (nm in names(seurat_list_filtered)) {
  obj <- seurat_list_filtered[[nm]]
  
  # Get the individual ID corresponding to this sample
  ind_id <- id_map$individualID[id_map$specimenID == nm]
  
  if (length(ind_id) == 1) {
    obj$individualID <- ind_id
    seurat_list_filtered[[nm]] <- obj
  } else {
    warning(paste0("No matching or multiple individualID found for: ", nm))
  }
}

# Group objects by individual ID
by_individual <- split(
  seurat_list_filtered,
  sapply(seurat_list_filtered, function(obj) unique(obj$individualID))
)

# Check how many were successfully mapped
mapped <- sapply(seurat_list_filtered, function(obj) !is.null(obj$individualID))
cat("# Successfully mapped objects:", sum(mapped), "out of", length(mapped), "\n")

# Show how many samples each individual has
individual_counts <- sapply(by_individual, length)
cat("# Number of samples per individual (first 5):\n")
print(head(individual_counts, 5))

# Merge all objects per individual
merged_by_individual <- lapply(by_individual, function(obj_list) {
  if (length(obj_list) == 1) {
    return(obj_list[[1]])  # No need to merge if only one object
  } else {
    Reduce(function(x, y) merge(x, y = y), obj_list)
  }
})

# Save merged objects to disk
output_path <- "/STORAGE/csbig/sc_ADers/qc_out/merged_by_individual.rds"
saveRDS(merged_by_individual, file = output_path)
cat("# Saved merged Seurat objects to:", output_path, "\n")

