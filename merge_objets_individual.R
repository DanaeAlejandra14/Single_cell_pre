# Create correspondence table between sample (specimenID) and individual (individualID)
id_map <- unique(bio_meta[, c("specimenID", "individualID")])  # remove duplicates if any

# Add the individualID field to each Seurat object in seurat_list_filtered
for (nm in names(seurat_list_filtered)) {
  obj <- seurat_list_filtered[[nm]]
  
  # Find the corresponding individualID
  ind_id <- id_map$individualID[id_map$specimenID == nm]
  
  if (length(ind_id) == 1) {
    obj$individualID <- ind_id
    seurat_list_filtered[[nm]] <- obj
  } else {
    warning(paste0("No matching or multiple individualID found for: ", nm))
  }
}

# Group the objects by individualID
by_individual <- split(
  seurat_list_filtered,
  sapply(seurat_list_filtered, function(obj) unique(obj$individualID))
)

# How many Seurat objects have an assigned individualID?
mapped <- sapply(seurat_list_filtered, function(obj) !is.null(obj$individualID))
table(mapped)  # We expect all to be TRUE

# See how many objects each individual has
sapply(by_individual, length)

# Show the first 3 individuals and the samples associated with each
lapply(by_individual[1:3], names)

# Merge objects per individual
merged_by_individual <- lapply(by_individual, function(obj_list) {
  if (length(obj_list) == 1) {
    return(obj_list[[1]])  # Already a single object
  } else {
    # Standard Seurat merge
    Reduce(function(x, y) merge(x, y = y), obj_list)
  }
})

# Save the merged objects to disk
saveRDS(merged_by_individual, file = "/datos/rosmap/single_cell/metadata/merged_by_individual.rds")

# Check for unmapped samples
ind_id <- id_map$individualID[id_map$specimenID == nm]

unmapped <- names(seurat_list_filtered)[!sapply(seurat_list_filtered, function(obj) !is.null(obj$individualID))]
length(unmapped)
print(unmapped)

# Total number of individuals expected from metadata
expected_individuals <- unique(bio_meta$individualID)
length(expected_individuals)  # Should be around 443 or 445

# Individuals that were actually mapped in your analysis
mapped_individuals <- names(by_individual)
length(mapped_individuals)  # Example: 438

# Which individuals are missing?
missing <- setdiff(expected_individuals, mapped_individuals)
length(missing)
print(missing)

# Read the original metadata file (if not already loaded)
bio_meta <- read.csv("/datos/rosmap/single_cell/metadata/ROSMAP_biospecimen_metadata.csv")

# Filter only snRNAseq samples and exclude flagged samples
bio_meta_filtered <- subset(bio_meta, assay == "scrnaSeq" & (is.na(exclude) | exclude == FALSE))

# Save the filtered metadata
write.csv(
  bio_meta_filtered,
  file = "/datos/rosmap/single_cell/metadata/ROSMAP_biospecimen_metadata_snRNAseq_filtered.csv",
  row.names = FALSE
)

# Get the specimenIDs used in your Seurat objects
specimen_ids_usados <- names(seurat_list_filtered)

# Filter the original metadata to only those specimenIDs
bio_meta_valid <- subset(bio_meta, specimenID %in% specimen_ids_usados)

# Save the valid metadata
write.csv(
  bio_meta_valid,
  file = "/datos/rosmap/single_cell/metadata/ROSMAP_biospecimen_metadata_final.csv",
  row.names = FALSE
)

# Check how many rows remain
cat("Number of filtered samples:", nrow(bio_meta_valid), "\n")
