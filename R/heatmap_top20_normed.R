suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
})
source("R/utils_normed.R")

out_dir <- "outputs_normed"
dir.create(out_dir, showWarnings = FALSE)

obj <- load_normed_expr("data/normed_cpms_filtered_annot.csv")
expr_log <- obj$expr_log
grp      <- obj$grp
samples  <- colnames(expr_log)

# --- Top-20 by NA-aware variance (no renaming) ---
N <- 20L
row_var <- apply(expr_log, 1, function(x) var(x[is.finite(x)], na.rm = TRUE))
row_var[!is.finite(row_var)] <- -Inf
ord <- order(row_var, decreasing = TRUE)
sel_n <- min(N, sum(is.finite(row_var) & row_var > -Inf))
top_idx <- ord[seq_len(sel_n)]
mat_top <- expr_log[top_idx, , drop = FALSE]

# --- Long + group-wise median imputation (fallback = gene median) ---
df <- as.data.frame(mat_top) |>
  tibble::rownames_to_column("Gene") |>
  tidyr::pivot_longer(-Gene, names_to = "Sample", values_to = "logCPM") |>
  dplyr::mutate(Group = grp[match(Sample, samples)])

df_imp <- df |>
  dplyr::group_by(Gene) |>
  dplyr::mutate(gene_med = median(logCPM[is.finite(logCPM)], na.rm = TRUE)) |>
  dplyr::group_by(Gene, Group, .add = TRUE) |>
  dplyr::mutate(
    grp_med = if (any(is.finite(logCPM))) median(logCPM[is.finite(logCPM)], na.rm = TRUE) else NA_real_,
    logCPM_imp = dplyr::case_when(
      is.finite(logCPM) ~ logCPM,
      is.finite(grp_med) ~ grp_med,
      TRUE ~ gene_med
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::select(Gene, Sample, logCPM = logCPM_imp)

# collapse duplicate (Gene, Sample) if any
df_imp <- df_imp |>
  dplyr::group_by(Gene, Sample) |>
  dplyr::summarise(logCPM = mean(logCPM, na.rm = TRUE), .groups = "drop")

# back to matrix, then row-wise Z
mat_imp <- df_imp |>
  tidyr::pivot_wider(names_from = Sample, values_from = logCPM) |>
  tibble::column_to_rownames("Gene") |>
  as.matrix()

row_z <- function(m){
  t(apply(m, 1, function(x){
    mu <- mean(x, na.rm = TRUE); sdv <- sd(x, na.rm = TRUE)
    if (!is.finite(sdv) || sdv == 0) return(rep(0, length(x)))
    (x - mu) / sdv
  }))
}
mat_z <- row_z(mat_imp)

# annotation + palette
ann <- data.frame(Group = grp); rownames(ann) <- colnames(mat_z)
pal <- colorRampPalette(c("#0d0887","#6a00a8","#b12a90","#e16462","#fca636","#f0f921"))

# outputs
readr::write_csv(tibble(Gene = rownames(mat_imp)),
                 file.path(out_dir, "Top20_variable_genes.csv"))

png(file.path(out_dir, "Heatmap_Top20_rowZ.png"), width = 1100, height = 800)
pheatmap(
  mat_z,
  color = pal(101),
  annotation_col = ann,
  show_rownames = TRUE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = sprintf("Top %d most variable genes (row Z of log2-CPM, imputed)", sel_n),
  fontsize = 11,
  border_color = NA
)
dev.off()

message("Saved: ",
        file.path(out_dir, "Heatmap_Top20_rowZ.png"), "; ",
        file.path(out_dir, "Top20_variable_genes.csv"))
