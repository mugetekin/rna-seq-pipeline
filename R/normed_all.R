# R/normed_all.R — consolidated normalized-data visualizations
# Run: Rscript R/normed_all.R
# This file inlines several *_normed.R scripts so you can run everything at once.
# Each section is wrapped in local({ ... }) to avoid object collisions.

suppressPackageStartupMessages({ library(tidyverse) })
message("=== Running consolidated normalized-data visualizations ===")


# --------------------[ BEGIN: R/sample_qc_normed.R ]--------------------
local({
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

})
# ---------------------[ END: R/sample_qc_normed.R ]---------------------


# --------------------[ BEGIN: R/pca_normed.R ]--------------------
local({
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

})
# ---------------------[ END: R/pca_normed.R ]---------------------


# --------------------[ BEGIN: R/violin_top20_normed.R ]--------------------
local({
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

})
# ---------------------[ END: R/violin_top20_normed.R ]---------------------


# --------------------[ BEGIN: R/heatmap_top20_normed.R ]--------------------
local({
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

})
# ---------------------[ END: R/heatmap_top20_normed.R ]---------------------


# --------------------[ BEGIN: R/heatmap_top50_normed.R ]--------------------
local({
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

message("Groups present and counts:"); print(table(grp))
message("Samples: ", paste(samples, collapse = ", "))

## --- Top-N by NA-aware variance (no name stripping) ---
N <- 50L
row_var <- apply(expr_log, 1, function(x) var(x[is.finite(x)], na.rm = TRUE))
row_var[!is.finite(row_var)] <- -Inf
ord <- order(row_var, decreasing = TRUE)
sel_n <- min(N, sum(is.finite(row_var) & row_var > -Inf))
top_idx <- ord[seq_len(sel_n)]
mat_top <- expr_log[top_idx, , drop = FALSE]

## --- Long form + diagnostics of NAs (no renaming of Gene!) ---
df <- as.data.frame(mat_top) |>
  tibble::rownames_to_column("Gene") |>
  tidyr::pivot_longer(-Gene, names_to = "Sample", values_to = "logCPM") |>
  dplyr::mutate(Group = grp[match(Sample, samples)])

na_diag <- df |>
  dplyr::group_by(Gene) |>
  dplyr::summarise(
    nNA_total   = sum(!is.finite(logCPM)),
    nNA_by_group= paste(tapply(!is.finite(logCPM), Group, sum), collapse = ";"),
    .groups = "drop"
  )

## --- Impute within Gene x Group (fallback to gene median) ---
df_imputed <- df |>
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
  dplyr::select(Gene, Sample, Group, logCPM = logCPM_imp)

## --- Collapse any duplicate (Gene, Sample) pairs by mean BEFORE widening ---
df_imputed <- df_imputed |>
  dplyr::group_by(Gene, Sample) |>
  dplyr::summarise(logCPM = mean(logCPM, na.rm = TRUE), .groups = "drop")

## --- Back to matrix safely (no stripping); enforce unique rownames if needed ---
mat_imp <- df_imputed |>
  tidyr::pivot_wider(names_from = Sample, values_from = logCPM) |>
  tibble::column_to_rownames("Gene") |>
  as.data.frame()

if (any(duplicated(rownames(mat_imp)))) {
  dups <- unique(rownames(mat_imp)[duplicated(rownames(mat_imp))])
  message("WARNING: duplicate gene names found after widening; making unique: ",
          paste(dups, collapse = ", "))
  rownames(mat_imp) <- make.unique(rownames(mat_imp))
}
mat_imp <- as.matrix(mat_imp)

## --- Column annotation
ann <- data.frame(Group = grp); rownames(ann) <- colnames(mat_imp)

## --- Palette
pal <- colorRampPalette(c("#0d0887","#6a00a8","#b12a90","#e16462","#fca636","#f0f921"))

## --- Outputs: Top list, diagnostics, heatmap
readr::write_csv(
  tibble(Gene = rownames(mat_top)),  # the selected rows (pre-imputation)
  file.path(out_dir, "Top50_variable_genes.csv")
)
readr::write_csv(
  na_diag |> dplyr::arrange(desc(nNA_total), Gene),
  file.path(out_dir, "Top50_imputation_diagnostics.csv")
)

png(file.path(out_dir, "Heatmap_Top50_logCPM.png"), width = 1200, height = 950)
pheatmap(
  mat_imp,
  color = pal(101),
  annotation_col = ann,
  show_rownames = FALSE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = sprintf("Top %d most variable genes (log2-CPM, group-wise median imputation)", sel_n),
  fontsize = 11,
  border_color = NA
)
dev.off()

message("Saved: ",
        file.path(out_dir, "Heatmap_Top50_logCPM.png"), "; ",
        file.path(out_dir, "Top50_variable_genes.csv"), "; ",
        file.path(out_dir, "Top50_imputation_diagnostics.csv"))

})
# ---------------------[ END: R/heatmap_top50_normed.R ]---------------------


# --------------------[ BEGIN: R/trend_genes_normed.R ]--------------------
local({
suppressPackageStartupMessages({ library(tidyverse); library(ggplot2) })
source("R/utils_normed.R")
out_dir <- "outputs_normed"; dir.create(out_dir, showWarnings = FALSE)

# --- REQUIRE Top20 list ---
top20_path <- file.path(out_dir, "Top20_variable_genes.csv")
if (!file.exists(top20_path)) {
  stop("Top20 list not found at: ", top20_path,
       "\nGenerate it first (e.g., run R/heatmap_top20_normed.R).")
}

# load data + groups
obj <- load_normed_expr("data/normed_cpms_filtered_annot.csv")
X   <- obj$expr_log
grp <- factor(obj$grp, levels = c("PBS","Lo","Med","Hi"))
samples <- colnames(X)

# read Top20 and preserve order
top20 <- readr::read_csv(top20_path, show_col_types = FALSE)$Gene
top20 <- unique(top20)  # safety
# subset to available genes (warn if any missing)
missing <- setdiff(top20, rownames(X))
if (length(missing)) message("WARNING: missing genes in matrix (will skip): ",
                             paste(missing, collapse = ", "))
genes <- top20[top20 %in% rownames(X)]
stopifnot(length(genes) > 0)

# long-form; keep Top20 order in facets
df <- as.data.frame(X[genes, , drop = FALSE]) |>
  tibble::rownames_to_column("Gene") |>
  tidyr::pivot_longer(-Gene, names_to = "Sample", values_to = "logCPM") |>
  dplyr::mutate(
    Group = grp[match(Sample, samples)],
    Gene  = factor(Gene, levels = genes)  # facet order = Top20 order
  )

# mean ± SEM
sem <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x)))
sumdf <- df |>
  dplyr::group_by(Gene, Group) |>
  dplyr::summarise(
    mean = mean(logCPM, na.rm = TRUE),
    sem  = sem(logCPM),
    .groups = "drop"
  )

# plot
p <- ggplot(sumdf, aes(Group, mean, group = Gene)) +
  geom_line(aes(color = Gene)) +
  geom_point(aes(color = Gene)) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.15) +
  facet_wrap(~ Gene, ncol = 4, scales = "free_y") +
  labs(title = "Top 20 genes — Mean ± SEM across dose groups",
       y = "log2(CPM+1)", x = "") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

ggsave(file.path(out_dir, "Trend_Genes_MeanSEM_Top20.png"), p,
       width = 14, height = 9, dpi = 300)
message("Saved: ", file.path(out_dir, "Trend_Genes_MeanSEM_Top20.png"))

})
# ---------------------[ END: R/trend_genes_normed.R ]---------------------


# --------------------[ BEGIN: R/trend_dose_normed.R ]--------------------
local({
suppressPackageStartupMessages({ library(tidyverse); library(ggplot2) })
source("R/utils_normed.R")
out_dir <- "outputs_normed"; dir.create(out_dir, showWarnings = FALSE)

obj <- load_normed_expr("data/normed_cpms_filtered_annot.csv")
X   <- obj$expr_log; grp <- factor(obj$grp, levels=c("PBS","Lo","Med","Hi"))
samples <- colnames(X)

# pick top K variable genes
K <- 12L
top_idx <- order(apply(X,1,var), decreasing = TRUE)[1:min(K, nrow(X))]
mat     <- X[top_idx, , drop=FALSE]

df <- as.data.frame(mat) |>
  tibble::rownames_to_column("Gene") |>
  tidyr::pivot_longer(-Gene, names_to="Sample", values_to="logCPM") |>
  dplyr::mutate(Group = grp[match(Sample, samples)])

sem <- function(x) sd(x)/sqrt(sum(is.finite(x)))
sumdf <- df |>
  dplyr::group_by(Gene, Group) |>
  dplyr::summarise(mean = mean(logCPM, na.rm=TRUE), sem = sem(logCPM), .groups="drop")

p <- ggplot(sumdf, aes(Group, mean, group = 1)) +
  geom_line() + geom_point() +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.15) +
  facet_wrap(~ Gene, ncol = 4, scales = "free_y") +
  labs(title = "Top variable genes — mean ± SEM across dose groups",
       y = "log2(CPM+1)", x = "") +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "Trend_TopVar_MeanSEM.png"), p, width = 14, height = 9, dpi = 300)
message("Saved: ", file.path(out_dir, "Trend_TopVar_MeanSEM.png"))

})
# ---------------------[ END: R/trend_dose_normed.R ]---------------------


# --------------------[ BEGIN: R/de_summary_normed.R ]--------------------
local({
suppressPackageStartupMessages({ library(tidyverse); library(ggplot2) })
source("R/utils_normed.R")
out_dir <- "outputs_normed"; dir.create(out_dir, showWarnings = FALSE)

obj <- load_normed_expr("data/normed_cpms_filtered_annot.csv")
X   <- obj$expr_log; grp <- obj$grp; samples <- colnames(X)
stopifnot(all(c("PBS","Hi") %in% levels(grp)))

pbs <- which(grp=="PBS"); hi <- which(grp=="Hi")

# per-gene Welch t-test with NA-safety
safe_t <- function(x){
  a <- x[pbs]; b <- x[hi]
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (length(a) < 2 || length(b) < 2) return(NA_real_)
  t.test(a, b, alternative="two.sided", var.equal=FALSE)$p.value
}

# compute stats
mu_pbs <- apply(X[, pbs, drop=FALSE], 1, function(v) mean(v[is.finite(v)], na.rm=TRUE))
mu_hi  <- apply(X[, hi,  drop=FALSE],  1, function(v) mean(v[is.finite(v)], na.rm=TRUE))
log2FC <- mu_hi - mu_pbs
pvals  <- apply(X, 1, safe_t)
padj   <- p.adjust(pvals, method = "BH")

res <- tibble(Gene = rownames(X), log2FC = log2FC, pval = pvals, padj = padj)

# how many skipped
nskip <- sum(!is.finite(res$pval))
message("Skipped genes (insufficient data for t-test): ", nskip)

# Volcano (drop NA p-values)
res_plot <- res |> dplyr::filter(is.finite(padj))
p1 <- ggplot(res_plot, aes(x = log2FC, y = -log10(padj))) +
  geom_point(alpha = 0.6, size = 1.4) +
  geom_vline(xintercept = c(-1, 1), linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  labs(title = "Mock Volcano (Hi vs PBS)", x = "log2FC (Hi - PBS)", y = "-log10(adj p)") +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "DE_Volcano_mock_Hi_vs_PBS.png"), p1, width = 8, height = 6, dpi = 300)

# Bar of up/down counts (ignore NA padj)
up   <- sum(res$padj < 0.05 & res$log2FC > 0, na.rm = TRUE)
down <- sum(res$padj < 0.05 & res$log2FC < 0, na.rm = TRUE)
bar  <- tibble(Direction=c("Up","Down"), Count=c(up,down))
p2 <- ggplot(bar, aes(Direction, Count, fill = Direction)) +
  geom_col(width = 0.6) +
  theme_minimal(base_size = 12) + theme(legend.position = "none") +
  labs(title = "Mock DEG counts (Hi vs PBS)")
ggsave(file.path(out_dir, "DE_BarCounts_mock_Hi_vs_PBS.png"), p2, width = 6, height = 5, dpi = 300)

# write table (keep NAs so downstream knows which were skipped)
readr::write_csv(res, file.path(out_dir, "DE_Mock_Hi_vs_PBS_results.csv"))
message("Saved: DE_Mock_Hi_vs_PBS_results.csv + volcano + bar counts")

})
# ---------------------[ END: R/de_summary_normed.R ]---------------------

message("All normalized plots finished. Check outputs_normed/")
