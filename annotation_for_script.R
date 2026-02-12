#anotation celular 

start_time <- Sys.time()

#inputs

celular_annotation <- "/STORAGE/csbig/sc_ADers/metadata/Experiment2/cell-annotation.full-atlas.csv"


celular_annotation <- readr::read_csv(celular_annotation, show_col_types = TRUE)
head(celular_annotation)

seurat_list_integrated <- readRDS("/STORAGE/csbig/sc_ADers/merge_integration_results_minimal/merged_by_individual_harmony.rds"
)
head(colnames(seurat_list_integrated))

#output
out_dir <- "/STORAGE/csbig/sc_ADers/out_annotated"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(out_dir)

# Clean annotation - Just keeping the cell type and state
celular_annotation_clean <- celular_annotation %>%
  dplyr::select(cell, cell.type, state) %>%
  dplyr::distinct()

# tablita
annotation_map <- celular_annotation_clean %>%
  tibble::column_to_rownames("cell")


# key to annotate the cells  

# En Seurat las celdas son ACTACGA...-1_1 y en la anno son : 190403-B4-A_ACTACGA...-1
# Entonces: quitamos _1 del colname y le pegamos el batch al inicio

barcode_clean <- sub("_[0-9]+$", "", colnames(seurat_list_integrated))  # quita _1, _2, etc.
batch <- seurat_list_integrated$libraryBatch                           # batch por célula (ej. 190403-B4-A)
cell_key <- paste0(batch, "_", barcode_clean)                           # el barcode final para la celula 


# Pegar anotación al metadata del Seurat 

seurat_list_integrated$cell_type  <- annotation_map[cell_key, "cell.type", drop = TRUE]
seurat_list_integrated$cell_state <- annotation_map[cell_key, "state", drop = TRUE]

# Crear una columna SOLO para graficar, donde NA se vuelve "Unassigned"
# (Esto no cambia tu cell_type original; solo ayuda a visualizar las no anotadas)
seurat_list_integrated$cell_type_plot <- ifelse(
  is.na(seurat_list_integrated$cell_type),
  "Unassigned",
  seurat_list_integrated$cell_type
)


# Revisar qué porcentaje se anotó y diagnosticar NA
total_cells <- ncol(seurat_list_integrated)
mean(!is.na(seurat_list_integrated$cell_type))
table(is.na(seurat_list_integrated$cell_type))

# Ejemplos de celulas que NO hicieron match (NA)
head(cell_key[is.na(seurat_list_integrated$cell_type)])

# En qué batches se concentran las NA
table(seurat_list_integrated$libraryBatch[is.na(seurat_list_integrated$cell_type)])

# Exportar tabla anotada (metadata completo)
write.csv(
  seurat_list_integrated@meta.data,
  file = "merged_harmony_integrated_cell_metadata_annotated.csv",
  row.names = TRUE
)

# Guardar Seurat ya anotado

saveRDS(
  seurat_list_integrated,
  file = "merged_harmony_integrated_annotated.rds"
)

#UMAPs anotados 

# UMAP con cell_type 
pdf("umap_cell_type.pdf", width = 10, height = 8)
print(DimPlot(seurat_list_integrated, group.by = "cell_type", label = TRUE) + NoLegend())
dev.off()

# UMAP con NA 
pdf("umap_cell_type_with_unassigned.pdf", width = 10, height = 8)
print(DimPlot(seurat_list_integrated, group.by = "cell_type_plot", label = TRUE) + NoLegend())
dev.off()

# UMAP por cell_state
pdf("umap_cell_state.pdf", width = 10, height = 8)
print(DimPlot(seurat_list_integrated, group.by = "cell_state", label = TRUE,   repel = TRUE) + NoLegend())
dev.off()


end_time <- Sys.time()
print(end_time - start_time)

#ENDD
