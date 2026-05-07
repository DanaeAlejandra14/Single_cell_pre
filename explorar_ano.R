#!/usr/bin/env Rscript

# ============================================================

#   Rscript explorar_objeto.R /ruta/a/tu/objeto.rds
#
# Genera: reporte_exploratorio.html en el directorio actual
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(scales)
})

# ── Leer argumento de línea de comandos ──────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  ruta <- "merged_harmony_integrated_annotated_plus_celltypist.rds"
} else {
  ruta <- args[1]
}

if (!file.exists(ruta)) stop(paste("ERROR: No se encontró el archivo:", ruta))

cat("============================================================\n")
cat("INICIANDO EXPLORACIÓN DEL OBJETO\n")
cat("============================================================\n")
cat("Archivo:", ruta, "\n")
cat("Inicio: ", format(Sys.time()), "\n\n")

# ── 1. CARGAR OBJETO ─────────────────────────────────────────
cat("[1/7] Cargando objeto...\n")
obj  <- readRDS(ruta)
meta <- obj@meta.data
cat(sprintf("      OK %s células | %s genes | %d columnas metadata\n",
            format(ncol(obj), big.mark=","),
            format(nrow(obj), big.mark=","),
            ncol(meta)))

# ── 2. HELPERS HTML ──────────────────────────────────────────
cat("[2/7] Preparando helpers...\n")

tabla_html <- function(df, caption="") {
  header <- paste0("<th>", colnames(df), "</th>", collapse="")
  filas  <- apply(df, 1, function(r) {
    paste0("<tr>", paste0("<td>", r, "</td>", collapse=""), "</tr>")
  })
  paste0(
    ifelse(caption!="", paste0("<h4>",caption,"</h4>"),""),
    "<div style='overflow-x:auto'>",
    "<table class='tabla'><thead><tr>",header,"</tr></thead><tbody>",
    paste(filas, collapse="\n"),
    "</tbody></table></div>"
  )
}

badges <- function(cols, clase="badge") {
  if (length(cols)==0) return("<em>Ninguna detectada</em>")
  paste0('<span class="',clase,'">',cols,'</span>', collapse=" ")
}

# ── 3. INFO GENERAL ──────────────────────────────────────────
cat("[3/7] Calculando info general...\n")

info_df <- data.frame(
  Parametro = c("Células totales","Genes totales",
                "Assays","Reducciones","Columnas metadata"),
  Valor = c(
    format(ncol(obj), big.mark=","),
    format(nrow(obj), big.mark=","),
    paste(Assays(obj),     collapse=", "),
    paste(Reductions(obj), collapse=", "),
    as.character(ncol(meta))
  )
)

# ── 4. RESUMEN METADATA ──────────────────────────────────────
cat("[4/7] Analizando columnas de metadata...\n")

resumen_meta <- data.frame(
  Columna  = colnames(meta),
  Tipo     = sapply(meta, function(x) class(x)[1]),
  N_unicos = sapply(meta, function(x) length(unique(x))),
  N_NA     = sapply(meta, function(x) sum(is.na(x))),
  Pct_NA   = paste0(sapply(meta, function(x) round(100*mean(is.na(x)),1)),"%")
)

# ── 5. DETECCIÓN DE COLUMNAS ─────────────────────────────────
cat("[5/7] Detectando columnas de anotación...\n")

cols_lower <- tolower(colnames(meta))

detectar <- function(keywords) {
  colnames(meta)[sapply(cols_lower, function(x)
    any(sapply(keywords, function(k) grepl(k, x, fixed=TRUE))))]
}

canon_cols <- detectar(c("cell_type","celltype","subtype","cluster",
                          "annotation","label","beyond","synapse",
                          "subpopulation","class","seurat_clusters"))

ct_cols    <- detectar(c("celltypist","predicted","majority","over_clust",
                          "conf_score","probability","majority_voting"))

clin_cols  <- detectar(c("age","sex","gender","diagnosis","braak","cerad",
                          "amyloid","tau","cogn","decline","batch","donor",
                          "participant","sample","pmi","rin","apoe",
                          "pathol","disease","condition","group","projid","study"))

# ── 6. CONTEOS POR ANOTACIÓN ─────────────────────────────────
cat("[6/7] Calculando conteos por anotación...\n")

bloque_anotacion <- function(col) {
  v       <- meta[[col]]
  n_total <- length(v)
  n_sin   <- sum(is.na(v)) + sum(v=="", na.rm=TRUE)
  n_con   <- n_total - n_sin

  resumen <- data.frame(
    Metrica = c("Total células","Con anotación","Sin anotación"),
    N       = c(format(n_total,big.mark=","),
                format(n_con,  big.mark=","),
                format(n_sin,  big.mark=",")),
    Pct     = c("100%",
                paste0(round(100*n_con/n_total,1),"%"),
                paste0(round(100*n_sin/n_total,1),"%"))
  )

  top <- head(as.data.frame(sort(table(v, useNA="ifany"), decreasing=TRUE)), 25)
  colnames(top) <- c("Tipo_celular","N_celulas")
  top$Porcentaje <- paste0(round(100*top$N_celulas/n_total,2),"%")

  list(resumen=resumen, top=top)
}

bloques_canon <- lapply(setNames(canon_cols, canon_cols), bloque_anotacion)
bloques_ct    <- lapply(setNames(ct_cols,    ct_cols),    bloque_anotacion)

# Cruce
html_cruce <- "<p><em>No se detectaron ambos tipos de anotación para el cruce.</em></p>"
if (length(canon_cols)>0 && length(ct_cols)>0) {
  col_c <- canon_cols[1]; col_t <- ct_cols[1]
  meta$tiene_canon <- !is.na(meta[[col_c]]) & meta[[col_c]] != ""
  meta$tiene_ct    <- !is.na(meta[[col_t]]) & meta[[col_t]] != ""
  cruce <- meta %>%
    mutate(Fuente = case_when(
       tiene_canon &  tiene_ct ~ "Ambas anotaciones",
       tiene_canon & !tiene_ct ~ "Solo canonica (BEYOND)",
      !tiene_canon &  tiene_ct ~ "Solo CellTypist",
       TRUE                    ~ "Sin anotacion"
    )) %>%
    count(Fuente) %>%
    mutate(Porcentaje = paste0(round(100*n/sum(n),1),"%")) %>%
    rename(N_celulas = n)
  html_cruce <- paste0('<div class="caja">', tabla_html(cruce, "Origen de la anotación por célula"), '</div>')
}

# Metadata clínica
clin_df <- if (length(clin_cols)>0) {
  data.frame(
    Columna = clin_cols,
    Tipo    = sapply(clin_cols, function(c) class(meta[[c]])[1]),
    N_NA    = sapply(clin_cols, function(c) sum(is.na(meta[[c]]))),
    Pct_NA  = paste0(sapply(clin_cols, function(c) round(100*mean(is.na(meta[[c]])),1)),"%")
  )
} else { data.frame(Mensaje="No se detectaron variables clínicas") }

# ── 7. CONSTRUIR Y ESCRIBIR HTML ─────────────────────────────
cat("[7/7] Escribiendo reporte HTML...\n")

css <- "
<style>
body{font-family:'Segoe UI',Arial,sans-serif;margin:40px;color:#2c3e50;background:#f8f9fa}
h1{color:#2c3e50;border-bottom:3px solid #3498db;padding-bottom:10px}
h2{color:#2980b9;border-left:4px solid #3498db;padding-left:10px;margin-top:40px}
h3{color:#34495e}
.tabla{border-collapse:collapse;width:100%;font-size:13px;margin-bottom:20px}
.tabla th{background:#3498db;color:white;padding:8px 12px;text-align:left}
.tabla td{padding:6px 12px;border-bottom:1px solid #dee2e6}
.tabla tr:nth-child(even){background:#ecf0f1}
.tabla tr:hover{background:#d5e8f7}
.caja{background:white;border-radius:8px;padding:20px;box-shadow:0 2px 6px rgba(0,0,0,0.08);margin-bottom:25px}
.badge{display:inline-block;background:#3498db;color:white;border-radius:4px;padding:2px 8px;font-size:12px;margin:2px}
.badge-green{background:#27ae60}
.badge-orange{background:#e67e22}
.grid2{display:grid;grid-template-columns:1fr 1fr;gap:20px}
@media(max-width:768px){.grid2{grid-template-columns:1fr}}
</style>"

# Bloques de anotación
render_bloques <- function(bloques) {
  html <- ""
  for (col in names(bloques)) {
    b <- bloques[[col]]
    html <- paste0(html,
      '<div class="caja"><h3>Columna: <code>',col,'</code></h3>',
      '<div class="grid2">',
      tabla_html(b$resumen, "Resumen"),
      tabla_html(b$top,     "Top 25 tipos celulares"),
      '</div></div>')
  }
  if (html=="") html <- "<p><em>No se detectaron columnas de este tipo.</em></p>"
  html
}

html_final <- paste0(
'<!DOCTYPE html><html lang="es"><head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Reporte Exploratorio Seurat</title>', css, '
</head><body>

<h1>Reporte Exploratorio del Objeto Seurat</h1>
<p><strong>Archivo:</strong> <code>',ruta,'</code><br>
<strong>Generado:</strong> ',format(Sys.time()),'</p>

<h2>1. Información General</h2>
<div class="caja">',tabla_html(info_df),'</div>

<h2>2. Columnas de Metadata</h2>
<div class="caja">',tabla_html(resumen_meta,"Todas las columnas en @meta.data"),'</div>

<h2>3. Columnas Detectadas por Tipo</h2>
<div class="caja">
<p><strong>Anotacion canonica (BEYOND/Synapse):</strong><br>',badges(canon_cols),'</p>
<p><strong>CellTypist:</strong><br>',badges(ct_cols,"badge badge-green"),'</p>
<p><strong>Variables clinicas:</strong><br>',badges(clin_cols,"badge badge-orange"),'</p>
</div>

<h2>4. Anotación Canónica (BEYOND / Synapse)</h2>',
render_bloques(bloques_canon),

'<h2>5. Anotación CellTypist</h2>',
render_bloques(bloques_ct),

'<h2>6. Cruce de Anotaciones</h2>',
html_cruce,

'<h2>7. Variables Clínicas Detectadas</h2>
<div class="caja">',tabla_html(clin_df),'</div>

<hr><p style="color:#999;font-size:12px;">
Generado con Rscript | ',format(Sys.time()),'
</p></body></html>'
)

salida <- "reporte_exploratorio.html"
writeLines(html_final, salida)

cat("\n============================================================\n")
cat("REPORTE GENERADO EXITOSAMENTE\n")
cat("============================================================\n")
cat(sprintf("  Archivo de salida:  %s\n", salida))
cat(sprintf("  Células:            %s\n", format(ncol(obj), big.mark=",")))
cat(sprintf("  Genes:              %s\n", format(nrow(obj), big.mark=",")))
cat(sprintf("  Cols canon.:        %s\n", paste(canon_cols, collapse=", ")))
cat(sprintf("  Cols CellTypist:    %s\n", paste(ct_cols,    collapse=", ")))
cat(sprintf("  Cols clínicas:      %s\n", paste(clin_cols,  collapse=", ")))
cat(sprintf("  Fin:                %s\n", format(Sys.time())))
cat("============================================================\n")
