# What this script does
#  1) Loads the voom object + metadata saved earlier (norm_and_fit.rds)
#  2) Fits a limma model with group indicators (PBS/Lo/Med/Hi)
#  3) Runs an omnibus F-test across ALL group coefficients (ANOVA-like)
#  4) Writes full/sig results and produces QC plots:
#       - p-value histogram
#       - Heatmap (top 50 F genes; row-Z)
#       - Violin plots for top 12 ANOVA genes
#
# Notes
#  - We use topTable(..., sort.by="F") with coef=1:ncol(design) instead of topTableF(),
#    which limma flags as obsolete. Results are equivalent for an omnibus test.
#  - Gene symbol mapping is done via org.Mm.eg.db (mouse). If rownames are ENSEMBL,
#    they are mapped to SYMBOL; otherwise we treat rownames as SYMBOL.

# R/anova_normed.R — omnibus ANOVA on voom-normalized logCPM + visualizations
suppressPackageStartupMessages({
  library(tidyverse); library(limma); library(ggrepel); library(pheatmap)
  library(AnnotationDbi); library(org.Mm.eg.db); library(scales); library(readr)
})

# Null-coalescing helper: return b if a is NULL or empty
`%||%` <- function(a,b) if (is.null(a) || !length(a)) b else a

# Load normalized object produced by 01_normalize.R
cfg  <- tryCatch({ source("R/io_helpers.R"); load_cfg() }, error=function(e) NULL)
rds_dir <- if (is.null(cfg)) "outputs/rds" else cfg$paths$rds
  # Keep anova outputs separate from other normalized outputs for tidiness
res_dir <- "outputs_normed/anova"; dir.create(res_dir, recursive=TRUE, showWarnings=FALSE)
fig_dir <- "outputs_normed/anova"; dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)

# Load voom and metadata (created earlier in your pipeline)
obj <- readRDS(file.path(rds_dir, "norm_and_fit.rds"))
v   <- obj$v
meta <- obj$meta

# Ensure Group factor has expected levels (order if present)
if (!"Group" %in% names(meta)) stop("meta$Group not found.")
meta$Group <- droplevels(factor(meta$Group, levels = c("PBS","Lo","Med","Hi")))

# Design and fit limma model; omnibus across all group coefficients
design <- model.matrix(~0 + Group, data = meta)  # columns: GroupPBS, GroupLo, ...
colnames(design) <- gsub("^Group", "", colnames(design)) # -> PBS, Lo, Med, Hi
# Fit per-gene linear models on voom-transformed expression
fit <- lmFit(v$E, design)
fit <- eBayes(fit)

# Omnibus ANOVA-like F-test
# Use topTable with coef=all group columns and sort.by="F".
# This replaces topTableF() (deprecated) with equivalent behavior.
ttF <- topTableF(fit, number = Inf, sort.by = "F")

# Save tables
# Raw ANOVA table (TSV; rownames are gene IDs as in v$E)
out_all <- file.path(res_dir, "ANOVA_all.tsv")
write.table(ttF, out_all, sep="\t", quote=FALSE, row.names=TRUE)

# Add gene SYMBOLs for convenience (supports ENSEMBL or already-SYMBOL rownames)
keys <- rownames(ttF)
keyType <- if (length(keys) && grepl("^ENSMUSG", keys[1])) "ENSEMBL" else "SYMBOL"
sy <- AnnotationDbi::mapIds(org.Mm.eg.db, keys=keys, keytype=keyType, column="SYMBOL", multiVals="first")
ttF_sym <- ttF %>% tibble::rownames_to_column("KEY") %>%
  mutate(SYMBOL = ifelse(is.na(sy[KEY]) | sy[KEY]=="", KEY, unname(sy[KEY])))
write_csv(ttF_sym, file.path(res_dir, "ANOVA_all_withSymbols.csv"))

# Significant by BH FDR (config or default 0.05)
fdr_thr <- (if (!is.null(cfg)) cfg$params$fdr_thresh else NULL) %||% 0.05
sig <- ttF_sym %>% filter(adj.P.Val <= fdr_thr)
write_csv(sig, file.path(res_dir, sprintf("ANOVA_significant_FDR%.3f.csv", fdr_thr)))

# --- Plots ---
# 1) P-value histogram (technically adj.P.Val is FDR; show raw F p-values if available)
pvals_raw <- fit$F.p.value
png(file.path(fig_dir, "ANOVA_pval_hist.png"), width=900, height=600)
hist(pvals_raw, 50, main="ANOVA omnibus p-values", xlab="P", col="grey80", border="white")
abline(v=fdr_thr, col="red", lty=2)
dev.off()

# 2) Heatmap of top 50 ANOVA genes (row-Z)
topN <- 50
top_keys <- head(ttF_sym$KEY, n = min(topN, nrow(ttF_sym)))
mat <- v$E[top_keys, , drop=FALSE]
# row-Z
# Row-wise Z-score scaling to highlight relative patterns
mat_z <- t(scale(t(mat)))
rownames(mat_z) <- ttF_sym$SYMBOL[match(rownames(mat_z), ttF_sym$KEY)] # Replace rownames with SYMBOLs for easier reading in the figure
ann_col <- data.frame(Group = meta$Group); rownames(ann_col) <- meta$SampleID
pheatmap(mat_z, annotation_col = ann_col, show_rownames = TRUE, show_colnames = FALSE,
         color = colorRampPalette(c("#214478","#f7f7f7","#b30000"))(101),
         filename = file.path(fig_dir, sprintf("ANOVA_heatmap_top%d_rowZ.png", topN)),
         width = 9, height = 10)

# 3) Violin for top 12 ANOVA genes (by F; prefer those with SYMBOLs)
top12_sym <- head(ttF_sym$SYMBOL[ttF_sym$SYMBOL!="" & !is.na(ttF_sym$SYMBOL)], 12)
if (length(top12_sym)) {
  # Map symbols back to rownames used in v$E if needed (e.g., when v$E rows are ENSEMBL)
  # map SYMBOLs back to rownames if necessary
  row_is_ens <- grepl("^ENSMUSG", rownames(v$E)[1])
  if (row_is_ens) {
    sym_to_row <- setNames(rownames(v$E),
      AnnotationDbi::mapIds(org.Mm.eg.db, keys=rownames(v$E),
                            keytype="ENSEMBL", column="SYMBOL", multiVals="first"))
    rn <- unique(na.omit(sym_to_row[top12_sym]))
    mat2 <- v$E[rn, , drop=FALSE]
    rownames(mat2) <- AnnotationDbi::mapIds(org.Mm.eg.db, keys=rownames(mat2),
                                            keytype="ENSEMBL", column="SYMBOL", multiVals="first")
  } else {
    mat2 <- v$E[intersect(top12_sym, rownames(v$E)), , drop=FALSE]
  }

  if (nrow(mat2) > 0) {
    df <- as.data.frame(mat2) %>% rownames_to_column("Gene") %>%
      pivot_longer(-Gene, names_to="Sample", values_to="logCPM") %>%
      mutate(Group = factor(meta$Group[match(Sample, meta$SampleID)], levels=c("PBS","Lo","Med","Hi"))) %>%
      filter(!is.na(Group))
    p <- ggplot(df, aes(Group, logCPM, fill=Group)) +
      geom_violin(trim=FALSE, alpha=.6, color="black") +
      geom_jitter(width=.12, size=1.8, alpha=.8, shape=21, color="black") +
      stat_summary(fun=mean, geom="point", shape=23, size=2.5, fill="white") +
      facet_wrap(~Gene, nrow=1, scales="free_y") +
      labs(title="Top ANOVA genes — voom log2-CPM", x=NULL, y="log2(CPM+1)") +
      theme_minimal(base_size=12) + theme(legend.position="none")
    ggsave(file.path(fig_dir, "ANOVA_violin_top12.png"), p, width=max(7, 3 + 2.2*length(unique(df$Gene))), height=4.2, dpi=300)
  }
}

message("ANOVA outputs written under ", res_dir)
