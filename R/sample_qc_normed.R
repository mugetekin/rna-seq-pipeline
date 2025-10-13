suppressPackageStartupMessages({ library(tidyverse); library(pheatmap) })
source("R/utils_normed.R")
out_dir <- "outputs_normed"; dir.create(out_dir, showWarnings = FALSE)

obj <- load_normed_expr("data/normed_cpms_filtered_annot.csv")
X   <- obj$expr_log; grp <- obj$grp

# pairwise-complete correlation across genes
cors <- cor(X, use = "pairwise.complete.obs", method = "pearson")
ann  <- data.frame(Group = grp); rownames(ann) <- colnames(X)

png(file.path(out_dir, "QC_SampleCorrelation.png"), width = 1100, height = 950)
pheatmap(cors, annotation_col = ann, annotation_row = ann,
         main = "Sample–sample correlation (Pearson)",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         border_color = NA)
dev.off()
message("Saved: ", file.path(out_dir, "QC_SampleCorrelation.png"))
