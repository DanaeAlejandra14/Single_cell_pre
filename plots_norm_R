#!/usr/bin/env Rscript

# ==== Evaluar normalización scran de objetos Seurat ====

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
})

# Leer argumento
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Uso: Rscript evaluar_normalizacion.R /ruta/al/seurat_list_scran_normalized.rds")
}
input_path <- args[1]
if (!file.exists(input_path)) stop("No se encuentra el archivo: ", input_path)

# Leer objeto
seurat_list <- readRDS(input_path)

# Elegir una muestra para evaluar
sample_name <- names(seurat_list)[1]
obj <- seurat_list[[sample_name]]
message("Evaluando muestra: ", sample_name)

# Calcular counts
obj$nCount_RNA   <- Matrix::colSums(GetAssayData(obj, assay = "RNA", layer = "counts"))
obj$nCount_scran <- Matrix::colSums(GetAssayData(obj, assay = "scran_norm", slot = "data"))

# Agregar al metadata
obj <- AddMetaData(obj, metadata = obj$nCount_RNA, col.name = "raw_counts")
obj <- AddMetaData(obj, metadata = obj$nCount_scran, col.name = "scran_norm_counts")

# Gráfico 1: size factors vs raw counts
p1 <- ggplot(obj@meta.data, aes(x = raw_counts, y = size_factors)) +
  geom_point(alpha = 0.3) +
  scale_x_log10() + scale_y_log10() +
  labs(title = "Size factors vs. Raw counts", x = "Raw counts", y = "Size factors") +
  theme_minimal()

# Gráfico 2: densidad de raw vs scran_norm
p2 <- ggplot(obj@meta.data, aes(x = log1p(raw_counts))) +
  geom_density(fill = "blue", alpha = 0.4) +
  geom_density(aes(x = log1p(scran_norm_counts)), fill = "red", alpha = 0.4) +
  labs(title = "Distribución Raw vs Scran (log1p)", x = "log1p(counts)", y = "Densidad") +
  theme_minimal()

# Guardar gráfico combinado
output_file <- "diagnostic_plot.png"
ggsave(output_file, p1 + p2, width = 10, height = 5, dpi = 300)
message("Gráfico guardado en: ", output_file)
