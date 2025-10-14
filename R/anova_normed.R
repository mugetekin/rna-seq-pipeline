# R/anova_normed.R — omnibus ANOVA on voom-normalized logCPM + visualizations
suppressPackageStartupMessages({
  library(tidyverse); library(limma); library(ggrepel); library(pheatmap)
  library(AnnotationDbi); library(org.Mm.eg.db); library(scales); library(readr)
})

`%||%` <- function(a,b) if (is.null(a) || !length(a)) b else a

# Load normalized object produced by 01_normalize.R
cfg  <- tryCatch({ source("R/io_helpers.R"); load_cfg() }, error=function(e) NULL)
rds_dir <- if (is.null(cfg)) "outputs/rds" else cfg$paths$rds
res_dir <- "outputs_normed/anova"; dir.create(res_dir, recursive=TRUE, showWarnings=FALSE)
fig_dir <- "outputs_normed/anova"; dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)

obj <- readRDS(file.path(rds_dir, "norm_and_fit.rds"))
v   <- obj$v
meta <- obj$meta

# Ensure Group factor has expected levels (order if present)
if (!"Group" %in% names(meta)) stop("meta$Group not found.")
meta$Group <- droplevels(factor(meta$Group, levels = c("PBS","Lo","Med","Hi")))

# Design and fit limma model; omnibus across all group coefficients
design <- model.matrix(~0 + Group, data = meta)  # PBS/Lo/Med/Hi columns
colnames(design) <- gsub("^Group", "", colnames(design))
fit <- lmFit(v$E, design)
fit <- eBayes(fit)

# Omnibus F-test table (sorted by F)
ttF <- topTableF(fit, number = Inf, sort.by = "F")

# Save tables
out_all <- file.path(res_dir, "ANOVA_all.tsv")
write.table(ttF, out_all, sep="\t", quote=FALSE, row.names=TRUE)

# Add SYMBOLs for convenience
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
mat_z <- t(scale(t(mat)))
rownames(mat_z) <- ttF_sym$SYMBOL[match(rownames(mat_z), ttF_sym$KEY)]
ann_col <- data.frame(Group = meta$Group); rownames(ann_col) <- meta$SampleID
pheatmap(mat_z, annotation_col = ann_col, show_rownames = TRUE, show_colnames = FALSE,
         color = colorRampPalette(c("#214478","#f7f7f7","#b30000"))(101),
         filename = file.path(fig_dir, sprintf("ANOVA_heatmap_top%d_rowZ.png", topN)),
         width = 9, height = 10)

# 3) Violin for top 12 ANOVA genes
top12_sym <- head(ttF_sym$SYMBOL[ttF_sym$SYMBOL!="" & !is.na(ttF_sym$SYMBOL)], 12)
if (length(top12_sym)) {
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
