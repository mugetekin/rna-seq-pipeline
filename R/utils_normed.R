suppressPackageStartupMessages({
  library(tidyverse); library(ggplot2)             # wrangling + plots
})

# Simple helper: cap extreme values (useful before plotting log fold-changes)
cap_values <- function(x, lower = -6, upper = 6) {
  pmax(lower, pmin(upper, x))                      # clamp between bounds
}

# Build a volcano-friendly tibble from a limma top table
volcano_df <- function(tt, lfc_col = "logFC", p_col = "adj.P.Val", gene_col = "gene",
                       lfc_thresh = 1, fdr_thresh = 0.05) {
  tt %>%
    transmute(
      gene = .data[[gene_col]],                    # gene name/label
      logFC = .data[[lfc_col]],                    # effect size (log2 fold change)
      negLog10FDR = -log10(pmax(.data[[p_col]], 1e-300)),  # display-friendly scale
      sig = abs(.data[[lfc_col]]) >= lfc_thresh & .data[[p_col]] <= fdr_thresh  # highlight
    )
}

# Minimal volcano plot factory (returns a ggplot object)
make_volcano <- function(df, title = NULL) {
  ggplot(df, aes(x = logFC, y = negLog10FDR)) +
    geom_point(aes(alpha = !sig), size = 1.2) +    # auto deemphasize non-significant points
    geom_vline(xintercept = c(-1, 1), linetype = 2) +
    geom_hline(yintercept = -log10(0.05), linetype = 2) +
    labs(title = title, x = "log2 fold-change", y = "-log10 FDR") +
    guides(alpha = "none") +
    theme_minimal(base_size = 12)
}

# Safe writer for figures
save_png <- function(plot, path, width = 1800, height = 1400, res = 150) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE) # ensure folder
  png(path, width = width, height = height, res = res)
  print(plot)
  dev.off()
}
