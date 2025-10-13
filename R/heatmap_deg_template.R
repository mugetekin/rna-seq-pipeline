suppressPackageStartupMessages({ library(tidyverse); library(pheatmap) })
source("R/utils_normed.R")
out_dir <- "outputs_normed"; dir.create(out_dir, showWarnings = FALSE)

# === EDIT ME LATER ===
de_csv      <- file.path(out_dir, "DE_results_example.csv")  # put your real DESeq2/edgeR table here
contrast_id <- "Hi_vs_PBS"                                   # which contrast to plot
top_each    <- 25                                            # top up + top down
# =====================

obj <- load_normed_expr("data/normed_cpms_filtered_annot.csv")
X   <- obj$expr_log; grp <- obj$grp

stopifnot(file.exists(de_csv))
de <- readr::read_csv(de_csv, show_col_types = FALSE) |>
  dplyr::filter(Contrast == contrast_id) |>
  dplyr::mutate(direction = dplyr::case_when(
    log2FC > 0 ~ "Up", log2FC < 0 ~ "Down", TRUE ~ "Neutral"
  )) |>
  dplyr::arrange(padj)

up   <- de |> dplyr::filter(direction=="Up")   |> dplyr::slice_head(n = top_each)
down <- de |> dplyr::filter(direction=="Down") |> dplyr::slice_head(n = top_each)
sel  <- dplyr::bind_rows(up, down)
genes <- unique(sel$Gene)
genes <- genes[genes %in% rownames(X)]
stopifnot(length(genes) > 1)

mat <- X[genes, , drop=FALSE]

# row Z for visualization
row_z <- function(m) t(apply(m,1,function(x){mu<-mean(x); sdv<-sd(x); if(!is.finite(sdv)||sdv==0) rep(0,length(x)) else (x-mu)/sdv}))
mat_z <- row_z(mat)

ann <- data.frame(Group = grp); rownames(ann) <- colnames(X)
png(file.path(out_dir, sprintf("Heatmap_DE_%s.png", contrast_id)), width = 1200, height = 950)
pheatmap(mat_z, annotation_col = ann, show_rownames = TRUE, cluster_rows = TRUE,
         cluster_cols = TRUE, main = sprintf("DE heatmap: %s (row Z)", contrast_id),
         border_color = NA)
dev.off()
message("Saved: ", file.path(out_dir, sprintf("Heatmap_DE_%s.png", contrast_id)))
