suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(factoextra)
})

# --- paths ---
in_csv <- "data/normed_cpms_filtered_annot.csv"
out_dir <- "outputs_normed"
dir.create(out_dir, showWarnings = FALSE)

# --- load normalized matrix ---
dat <- data.table::fread(in_csv, data.table = FALSE)

# --- build expression matrix + rownames (prefer SYMBOL, fallback to id) ---
gene_cols <- c("id","SYMBOL","gene_type")
expr <- dat[ , !(names(dat) %in% gene_cols), drop = FALSE]
genes <- ifelse(is.na(dat$SYMBOL) | dat$SYMBOL=="", dat$id, dat$SYMBOL)
rownames(expr) <- make.unique(genes)

# --- numeric, log2 transform, clean ---
expr_num <- as.matrix(data.frame(lapply(as.data.frame(expr), as.numeric),
                                 check.names = FALSE, row.names = rownames(expr)))
expr_log <- log2(expr_num + 1)
keep_rows <- apply(is.finite(expr_log), 1, all)
expr_log  <- expr_log[keep_rows, , drop = FALSE]
keep_cols <- apply(is.finite(expr_log), 2, all)
expr_log  <- expr_log[, keep_cols, drop = FALSE]
nzv <- apply(expr_log, 1, sd) > 0
expr_log <- expr_log[nzv, , drop = FALSE]

# --- groups from sample IDs ---
sample_ids <- colnames(expr_log)
grp <- dplyr::case_when(
  startsWith(sample_ids,"PBS") ~ "PBS",
  startsWith(sample_ids,"Lo")  ~ "Lo",
  startsWith(sample_ids,"Med") ~ "Med",
  startsWith(sample_ids,"Hi")  ~ "Hi",
  TRUE ~ "Other"
)
grp <- factor(grp, levels = c("PBS","Lo","Med","Hi","Other"))

# --- PCA ---
pca <- prcomp(t(expr_log), scale. = TRUE)

# --- plot -> outputs_normed/PCA_normed_CPM.png ---
png(file.path(out_dir, "PCA_normed_CPM.png"), width = 900, height = 700)
print(
  fviz_pca_ind(
    pca,
    geom.ind = "point",
    col.ind  = grp,
    pointsize = 3,
    pointshape = 16,
    addEllipses = TRUE,
    ellipse.level = 0.95,
    legend.title = list(color = "Group", shape = "Group")
  ) +
    guides(linetype = "none") +
    theme_minimal(base_size = 12)
)
dev.off()
message("Saved: ", file.path(out_dir, "PCA_normed_CPM.png"))
