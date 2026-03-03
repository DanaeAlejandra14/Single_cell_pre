# anotation celular + export Unassigned to h5ad

start_time <- Sys.time()

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(readr)
})

# =========================
# Inputs
# =========================
anno_path <- "/STORAGE/csbig/sc_ADers/metadata/Experiment2/cell-annotation.full-atlas.csv"
seurat_path <- "/STORAGE/csbig/sc_ADers/merge_integration_results_minimal/merged_by_individual_harmony.rds"

celular_annotation <- readr::read_csv(anno_path, show_col_types = TRUE)
seurat_list_integrated <- readRDS(seurat_path)

cat("# Loaded Seurat cells:", ncol(seurat_list_integrated), "\n")
cat("# Loaded annotation rows:", nrow(celular_annotation), "\n")

# =========================
# Output dir (NEW)
# =========================
out_dir <- "/STORAGE/csbig/sc_ADers/out_annotated_unassigned"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(out_dir)
cat("# Output dir:", out_dir, "\n")

# =========================
# Clean annotation
# =========================
celular_annotation_clean <- celular_annotation %>%
  dplyr::select(cell, cell.type, state) %>%
  dplyr::distinct()

annotation_map <- celular_annotation_clean %>%
  tibble::column_to_rownames("cell")

csv_keys <- rownames(annotation_map)

# =========================
# Build cell_key (batch_BARCODE)
# =========================
barcode_clean <- sub("_[0-9]+$", "", colnames(seurat_list_integrated))  # quita _1, _2...
batch <- seurat_list_integrated$libraryBatch                           # ej. 190403-B4-A
cell_key <- paste0(batch, "_", barcode_clean)

# =========================
# Annotate Seurat metadata
# =========================
seurat_list_integrated$cell_type  <- annotation_map[cell_key, "cell.type", drop = TRUE]
seurat_list_integrated$cell_state <- annotation_map[cell_key, "state", drop = TRUE]

seurat_list_integrated$cell_type_plot <- ifelse(
  is.na(seurat_list_integrated$cell_type),
  "Unassigned",
  seurat_list_integrated$cell_type
)

# =========================
# Summary + diagnostics
# =========================
total_cells <- ncol(seurat_list_integrated)
annotated_cells <- sum(!is.na(seurat_list_integrated$cell_type))
unannotated_cells <- sum(is.na(seurat_list_integrated$cell_type))

cat("\n# ===== Annotation summary =====\n")
cat("Total cells: ", total_cells, "\n")
cat("Annotated cells: ", annotated_cells, "\n")
cat("Unassigned cells (NA): ", unannotated_cells, "\n")
cat("Percent annotated: ", round(annotated_cells / total_cells * 100, 2), "%\n")
cat("Percent unassigned: ", round(unannotated_cells / total_cells * 100, 2), "%\n")


writeLines(
  c(
    paste0("Total cells: ", total_cells),
    paste0("Annotated cells: ", annotated_cells),
    paste0("Unassigned cells (NA): ", unannotated_cells),
    paste0("Percent annotated: ", round(annotated_cells / total_cells * 100, 2), "%"),
    paste0("Percent unassigned: ", round(unannotated_cells / total_cells * 100, 2), "%")
  ),
  con = "annotation_summary.txt"
)

cat("\n# Example NA keys:\n")
print(head(cell_key[is.na(seurat_list_integrated$cell_type)], 20))

cat("\n# NA per libraryBatch:\n")
print(sort(table(seurat_list_integrated$libraryBatch[is.na(seurat_list_integrated$cell_type)]), decreasing = TRUE))

# =========================
# Save annotated outputs (same as you had)
# =========================
write.csv(
  seurat_list_integrated@meta.data,
  file = "merged_harmony_integrated_cell_metadata_annotated.csv",
  row.names = TRUE
)

saveRDS(
  seurat_list_integrated,
  file = "merged_harmony_integrated_annotated.rds"
)

# =========================
# UMAPs
# =========================
pdf("umap_cell_type.pdf", width = 10, height = 8)
print(DimPlot(seurat_list_integrated, group.by = "cell_type", label = TRUE) + NoLegend())
dev.off()

pdf("umap_cell_type_with_unassigned.pdf", width = 10, height = 8)
print(DimPlot(seurat_list_integrated, group.by = "cell_type_plot", label = TRUE) + NoLegend())
dev.off()

pdf("umap_cell_state.pdf", width = 10, height = 8)
print(DimPlot(seurat_list_integrated, group.by = "cell_state", label = TRUE, repel = TRUE) + NoLegend())
dev.off()

# ============================================================
# NEW PART: Create Unassigned subset + export to h5ad
# ============================================================
cat("\n# ===== Creating Unassigned subset =====\n")

unassigned_cells <- colnames(seurat_list_integrated)[is.na(seurat_list_integrated$cell_type)]
cat("# Unassigned cells:", length(unassigned_cells), "\n")

seurat_unassigned <- subset(seurat_list_integrated, cells = unassigned_cells)

# trazabilidad para Python
seurat_unassigned$orig_cellname <- colnames(seurat_unassigned)
# guardar cell_key correspondiente (mismo orden por cellname)
cell_key_named <- setNames(cell_key, colnames(seurat_list_integrated))
seurat_unassigned$cell_key <- cell_key_named[colnames(seurat_unassigned)]

# guardar objetos/listas
saveRDS(seurat_unassigned, file = "unassigned_subset.rds")
writeLines(colnames(seurat_unassigned), con = "unassigned_cellnames.txt")
writeLines(seurat_unassigned$cell_key, con = "unassigned_cellkeys.txt")

readr::write_csv(
  tibble::tibble(
    cellname = colnames(seurat_unassigned),
    cell_key = seurat_unassigned$cell_key,
    libraryBatch = seurat_unassigned$libraryBatch
  ),
  "unassigned_mapping_for_python.csv"
)

# export h5ad usando SeuratDisk
cat("\n# ===== Exporting to h5ad =====\n")

if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
  cat("# SeuratDisk not found. Installing...\n")
  install.packages("SeuratDisk")
}
suppressPackageStartupMessages(library(SeuratDisk))

# asegúrate de exportar desde RNA
if (!("RNA" %in% Assays(seurat_unassigned))) {
  stop("No RNA assay found in seurat_unassigned. Available assays: ",
       paste(Assays(seurat_unassigned), collapse = ", "))
}
DefaultAssay(seurat_unassigned) <- "RNA"

# Nota: si tu objeto NO tiene slot data, puedes descomentar NormalizeData:
# seurat_unassigned <- NormalizeData(seurat_unassigned)

SaveH5Seurat(seurat_unassigned, filename = "unassigned.h5Seurat", overwrite = TRUE)
Convert("unassigned.h5Seurat", dest = "h5ad", overwrite = TRUE)

cat("# Wrote: unassigned_subset.rds, unassigned.h5Seurat, unassigned.h5ad\n")

end_time <- Sys.time()
cat("\n# DONE. Time taken: ", end_time - start_time, "\n")
# ENDD
