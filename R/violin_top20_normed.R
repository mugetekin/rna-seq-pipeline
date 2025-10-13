suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})
source("R/utils_normed.R")

out_dir <- "outputs_normed"
dir.create(out_dir, showWarnings = FALSE)

# load matrix + groups + reuse Top20 selection logic
obj <- load_normed_expr("data/normed_cpms_filtered_annot.csv")
expr_log <- obj$expr_log
grp      <- obj$grp
samples  <- colnames(expr_log)

row_var <- apply(expr_log, 1, function(x) var(x[is.finite(x)], na.rm = TRUE))
row_var[!is.finite(row_var)] <- -Inf
top_idx <- order(row_var, decreasing = TRUE)[seq_len(min(20L, sum(is.finite(row_var)&row_var>-Inf)))]
mat20   <- expr_log[top_idx, , drop = FALSE]

# long form + group-wise median imputation for plotting only
df <- as.data.frame(mat20) |>
  tibble::rownames_to_column("Gene") |>
  tidyr::pivot_longer(-Gene, names_to = "Sample", values_to = "logCPM") |>
  dplyr::mutate(Group = grp[match(Sample, samples)])

df_imp <- df |>
  dplyr::group_by(Gene) |>
  dplyr::mutate(gene_med = median(logCPM[is.finite(logCPM)], na.rm = TRUE)) |>
  dplyr::group_by(Gene, Group, .add = TRUE) |>
  dplyr::mutate(
    grp_med = if (any(is.finite(logCPM))) median(logCPM[is.finite(logCPM)], na.rm = TRUE) else NA_real_,
    logCPM_plot = dplyr::case_when(
      is.finite(logCPM) ~ logCPM,
      is.finite(grp_med) ~ grp_med,
      TRUE ~ gene_med
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::select(Gene, Sample, Group, logCPM = logCPM_plot) |>
  dplyr::mutate(Gene = factor(Gene, levels = rownames(mat20)))

set.seed(42)
p <- ggplot(df_imp, aes(x = Group, y = logCPM, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = "black") +
  geom_jitter(width = 0.10, size = 1.8, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2.2, fill = "white") +
  facet_wrap(~ Gene, ncol = 5, scales = "free_y") +
  labs(title = "Top 20 variable genes (normalized log2-CPM)",
       x = "", y = "log2(CPM + 1)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

ggsave(file.path(out_dir, "Violin_Top20_logCPM.png"), p, width = 16, height = 12, dpi = 300)
message("Saved: ", file.path(out_dir, "Violin_Top20_logCPM.png"))
