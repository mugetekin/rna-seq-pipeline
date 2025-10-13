suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
})

source("R/utils_normed.R")
out_dir <- "outputs_normed"
dir.create(out_dir, showWarnings = FALSE)

# load matrix + groups
obj <- load_normed_expr("data/normed_cpms_filtered_annot.csv")
expr_log <- obj$expr_log
grp      <- obj$grp

# choose top N by variance (handle small matrices gracefully)
N <- 50L
n_avail <- nrow(expr_log)
sel_n <- min(N, n_avail)
ord <- order(apply(expr_log, 1, var), decreasing = TRUE)
top_idx <- ord[seq_len(sel_n)]
mat_top <- expr_log[top_idx, , drop = FALSE]

# column annotation
ann <- data.frame(Group = grp)
rownames(ann) <- colnames(expr_log)

# optional palette (nice contrast but not too aggressive)
pal <- colorRampPalette(c("#0d0887","#6a00a8","#b12a90","#e16462","#fca636","#f0f921"))

# save the list of selected genes
readr::write_csv(
  tibble(Gene = rownames(mat_top)),
  file.path(out_dir, "Top50_variable_genes.csv")
)

# draw heatmap
png(file.path(out_dir, "Heatmap_Top50_logCPM.png"), width = 1200, height = 950)
pheatmap(
  mat_top,
  color = pal(101),
  annotation_col = ann,
  show_rownames = FALSE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = sprintf("Top %d most variable genes (log2-CPM)", sel_n),
  fontsize = 11,
  border_color = NA
)
dev.off()

message("Saved: ",
        file.path(out_dir, "Heatmap_Top50_logCPM.png"),
        " and Top50 list CSV.")
