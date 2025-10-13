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

# diagnostics: show groups and sample counts
message("Groups present and counts:"); print(table(grp))
message("Samples seen: ", paste(colnames(expr_log), collapse = ", "))

# choose top N by variance (NA-aware)
N <- 50L
row_var <- apply(expr_log, 1, function(x) var(x[is.finite(x)], na.rm = TRUE))
row_var[!is.finite(row_var)] <- -Inf
ord <- order(row_var, decreasing = TRUE)
sel_n <- min(N, sum(is.finite(row_var) & row_var > -Inf))
top_idx <- ord[seq_len(sel_n)]
mat_top <- expr_log[top_idx, , drop = FALSE]

# impute NAs by row-median for visualization
mat_top <- t(apply(mat_top, 1, function(x){
  if (all(!is.finite(x))) return(rep(0, length(x)))
  xr <- x
  xr[!is.finite(xr)] <- median(xr[is.finite(xr)], na.rm = TRUE)
  xr
}))
rownames(mat_top) <- rownames(expr_log)[top_idx]

# column annotation
ann <- data.frame(Group = grp)
rownames(ann) <- colnames(expr_log)

# palette
pal <- colorRampPalette(c("#0d0887","#6a00a8","#b12a90","#e16462","#fca636","#f0f921"))

# save Top list
readr::write_csv(tibble(Gene = rownames(mat_top)),
                 file.path(out_dir, "Top50_variable_genes.csv"))

# plot
png(file.path(out_dir, "Heatmap_Top50_logCPM.png"), width = 1200, height = 950)
pheatmap(
  mat_top,
  color = pal(101),
  annotation_col = ann,
  show_rownames = FALSE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = sprintf("Top %d most variable genes (log2-CPM, NA-imputed row median)", sel_n),
  fontsize = 11,
  border_color = NA
)
dev.off()

message("Saved: ",
        file.path(out_dir, "Heatmap_Top50_logCPM.png"),
        " and Top50_variable_genes.csv")
