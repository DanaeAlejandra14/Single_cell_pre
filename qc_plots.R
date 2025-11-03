#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage:\n  Rscript qc_plots_final.R <RUN_DIR> [min_feat] [min_count] [mt_max] [down_per]\n",
       "Ex:\n  Rscript qc_plots_final.R /path/qc_out/XXXX 700 1500 13 10000\n")
}
run_dir   <- normalizePath(args[1], mustWork = TRUE)
min_feat  <- ifelse(length(args) >= 2, as.integer(args[2]), 700L)
min_count <- ifelse(length(args) >= 3, as.integer(args[3]), 1500L)
mt_max    <- ifelse(length(args) >= 4, as.numeric(args[4]), 13)
down_per  <- ifelse(length(args) >= 5, as.integer(args[5]), 10000L)

cat("# Reading from:\n ", run_dir, "\n")

# Input files
f_pre <- file.path(run_dir, "qc_cells_pre.csv")
f_dbl <- file.path(run_dir, "doublet_summary.csv")
stopifnot(file.exists(f_pre))

# Output dir
plots_root <- file.path(dirname(run_dir), "qc_plots_final")
out_dir    <- file.path(plots_root, basename(run_dir))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat("# Saving to:\n ", out_dir, "\n")

# Load PRE 
pre <- readr::read_csv(f_pre, show_col_types = FALSE)

# Normalize mt column name if needed
if ("percent.mt" %in% names(pre) && !("percent_mt" %in% names(pre))) {
  pre <- dplyr::rename(pre, percent_mt = "percent.mt")
}

# Basic checks
need_cols <- c("sample_id","cell_barcode","nFeature_RNA","nCount_RNA","percent_mt")
miss <- setdiff(need_cols, names(pre))
if (length(miss)) stop("Missing columns in qc_cells_pre.csv: ", paste(miss, collapse = ", "))

# Flags 
# Keep after counts & features thresholds (no mt yet)
pre <- pre %>%
  mutate(keep_cf = if_else(nFeature_RNA >= min_feat & nCount_RNA >= min_count, "keep", "remove"),
         keep_mt = if_else(percent_mt <= mt_max, "keep", "remove")) %>%
  mutate(keep_all = if_else(keep_cf == "keep" & keep_mt == "keep", "keep", "remove"))

# Robust downsample per sample 
downsample_per_sample <- function(df, by = "sample_id", cap = 10000L) {
  cap <- as.integer(cap)
  if (is.na(cap) || cap <= 0L) return(df)
  parts <- split(df, df[[by]], drop = TRUE)
  sampled <- lapply(parts, function(dd) {
    m <- min(cap, nrow(dd))
    if (m >= nrow(dd)) return(dd)
    dd[sample.int(nrow(dd), size = m), , drop = FALSE]
  })
  bind_rows(sampled)
}

set.seed(123)
pre_ds <- downsample_per_sample(pre, by = "sample_id", cap = down_per)

# 1) Global histograms with threshold lines 
mk_hist <- function(df, x, cut = NULL, title = x, bins = 80, logx = FALSE) {
  p <- ggplot(df, aes(x = .data[[x]])) +
    geom_histogram(bins = bins, fill = "grey50") +
    labs(title = title, x = x, y = "cells") +
    theme_bw(11)
  if (!is.null(cut)) p <- p + geom_vline(xintercept = cut, color = "red", linetype = 2, linewidth = 0.7)
  if (isTRUE(logx))  p <- p + scale_x_continuous(trans = "log10", labels = label_number(accuracy = 1))
  p
}

p_hist_counts  <- mk_hist(pre, "nCount_RNA",   cut = min_count, title = "Global nCount_RNA (red = min)",     bins = 80, logx = TRUE)
p_hist_feat    <- mk_hist(pre, "nFeature_RNA", cut = min_feat,  title = "Global nFeature_RNA (red = min)",   bins = 80, logx = FALSE)
p_hist_mt      <- mk_hist(pre, "percent_mt",   cut = mt_max,    title = "Global percent_mt (red = max)",     bins = 80, logx = FALSE)

ggsave(file.path(out_dir, "01_hist_nCount_global.png"),   p_hist_counts, width = 9, height = 6, dpi = 140)
ggsave(file.path(out_dir, "02_hist_nFeature_global.png"), p_hist_feat,   width = 9, height = 6, dpi = 140)
ggsave(file.path(out_dir, "03_hist_percent_mt_global.png"), p_hist_mt,   width = 9, height = 6, dpi = 140)

# 2) Scatter counts vs features (log scales), colored by keep/remove
p_sc_cf <- ggplot(pre_ds, aes(x = nCount_RNA, y = nFeature_RNA, color = keep_cf)) +
  geom_point(alpha = 0.4, size = 0.4) +
  scale_x_continuous(trans = "log10", labels = label_number(accuracy = 1)) +
  scale_y_continuous(trans = "log10", labels = label_number(accuracy = 1)) +
  scale_color_manual(values = c(keep = "#1b9e77", remove = "#d95f02")) +
  labs(title = "Counts vs Features (log–log) — keep/remove by counts+features",
       x = "nCount_RNA (log10)", y = "nFeature_RNA (log10)", color = "") +
  theme_bw(11)
ggsave(file.path(out_dir, "04_scatter_counts_vs_features_keepCF.png"), p_sc_cf, width = 9, height = 7, dpi = 140)

# 3) Violins per sample with threshold lines 
vln_theme <- theme_bw(11) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_vln_nc <- ggplot(pre_ds, aes(x = sample_id, y = nCount_RNA)) +
  geom_violin(fill = "#8cc5f2", alpha = 0.8, scale = "width") +
  geom_hline(yintercept = min_count, linetype = 2, color = "red") +
  coord_flip() + labs(title = "nCount_RNA by sample", x = "sample", y = "nCount_RNA") + vln_theme

p_vln_nf <- ggplot(pre_ds, aes(x = sample_id, y = nFeature_RNA)) +
  geom_violin(fill = "#f28e8c", alpha = 0.8, scale = "width") +
  geom_hline(yintercept = min_feat, linetype = 2, color = "red") +
  coord_flip() + labs(title = "nFeature_RNA by sample", x = "sample", y = "nFeature_RNA") + vln_theme

p_vln_mt <- ggplot(pre_ds, aes(x = sample_id, y = percent_mt)) +
  geom_violin(fill = "#b8e3a1", alpha = 0.8, scale = "width") +
  geom_hline(yintercept = mt_max, linetype = 2, color = "red") +
  coord_flip() + labs(title = "percent_mt by sample", x = "sample", y = "percent_mt") + vln_theme

ggsave(file.path(out_dir, "05_violins_by_sample.png"),
       (p_vln_nc / p_vln_nf / p_vln_mt), width = 14, height = 18, dpi = 130)

# 4) Doublets (optional) 
if (file.exists(f_dbl)) {
  dbl <- readr::read_csv(f_dbl, show_col_types = FALSE)
  # Expected columns: sample, singlets, doublets, frac_doublet
  if (all(c("sample","singlets","doublets","frac_doublet") %in% names(dbl))) {
    dbl$frac_doublet <- as.numeric(dbl$frac_doublet)
    p_dbl <- dbl %>%
      arrange(desc(frac_doublet)) %>%
      slice_head(n = 30) %>%
      mutate(sample = factor(sample, levels = rev(sample))) %>%
      ggplot(aes(x = sample, y = frac_doublet)) +
      geom_col(fill = "#7f7fff") +
      coord_flip() +
      scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
      labs(title = "Doublet fraction per sample (top 30)", x = "sample", y = "doublet fraction") +
      theme_bw(11)
    ggsave(file.path(out_dir, "06_doublets_per_sample.png"), p_dbl, width = 10, height = 10, dpi = 140)
  } else {
    message("# doublet_summary.csv present but with unexpected columns; skipping plot.")
  }
} else {
  message("# No doublet_summary.csv; skipping doublets plot.")
}

# 5) Scatter counts vs features colored by MT keep/remove 
p_sc_mt <- ggplot(pre_ds, aes(x = nCount_RNA, y = nFeature_RNA, color = keep_mt)) +
  geom_point(alpha = 0.4, size = 0.4) +
  scale_x_continuous(trans = "log10", labels = label_number(accuracy = 1)) +
  scale_y_continuous(trans = "log10", labels = label_number(accuracy = 1)) +
  scale_color_manual(values = c(keep = "#1b9e77", remove = "#d95f02")) +
  labs(title = "Counts vs Features (log–log) — keep/remove by mitochondrial %",
       x = "nCount_RNA (log10)", y = "nFeature_RNA (log10)", color = "") +
  theme_bw(11)
ggsave(file.path(out_dir, "07_scatter_counts_vs_features_keepMT.png"), p_sc_mt, width = 9, height = 7, dpi = 140)

#  6) Mini summary (pre, post counts+features, post doublets, post mt) 
total_pre <- nrow(pre)
after_cf  <- sum(pre$keep_cf == "keep")

if (file.exists(f_dbl)) {
  dbl <- readr::read_csv(f_dbl, show_col_types = FALSE)
  if (all(c("singlets","doublets") %in% names(dbl))) {
    after_doublets <- sum(dbl$singlets, na.rm = TRUE)
  } else after_doublets <- NA_integer_
} else {
  after_doublets <- NA_integer_
}

# Apply mt on top of counts+features (approximate final)
after_mt_on_cf <- sum(pre$keep_cf == "keep" & pre$percent_mt <= mt_max)

mini <- tibble(
  metric = c("cells_pre", "cells_after_counts_features", "cells_after_doublets", "cells_after_counts_features_and_mt"),
  value  = c(total_pre, after_cf, after_doublets, after_mt_on_cf)
)
readr::write_csv(mini, file.path(out_dir, "00_qc_mini_summary.csv"))
cat("# Mini summary written: ", file.path(out_dir, "00_qc_mini_summary.csv"), "\n")

cat("# Done. Plots in: ", out_dir, "\n")
