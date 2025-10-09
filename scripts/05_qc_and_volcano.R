# QC + volcano + MA + PCA + hist
source("R/io_helpers.R")
suppressPackageStartupMessages({
  library(edgeR); library(limma); library(tidyverse)
  library(ggrepel); library(factoextra); library(ragg)
  library(AnnotationDbi); library(org.Mm.eg.db)
})

cfg <- load_cfg()
obj <- readRDS(file.path(cfg$paths$rds,"norm_and_fit.rds"))
annot <- obj$annot; meta <- obj$meta; dge <- obj$dge; v <- obj$v; fit2 <- obj$fit2

dir.create(cfg$paths$figures, showWarnings=FALSE, recursive=TRUE)

png(file.path(cfg$paths$figures,"library_sizes.png"), width=900, height=600)
barplot(colSums(dge$counts), las=2, main="Library sizes", ylab="Total counts")
dev.off()

png(file.path(cfg$paths$figures,"MDS_groups.png"), width=700, height=700)
plotMDS(dge, labels = meta$Group); title("MDS (TMM-normalized)")
dev.off()

pdf(file.path(cfg$paths$figures,"voom_meanvar.pdf"), 6, 5)
voom(dge, model.matrix(~0+Group, data=meta), plot=TRUE)
dev.off()

tt_lo  <- topTable(fit2, coef="Lo_vs_PBS",  number=Inf, sort.by="P")
tt_med <- topTable(fit2, coef="Med_vs_PBS", number=Inf, sort.by="P")
tt_hi  <- topTable(fit2, coef="Hi_vs_PBS",  number=Inf, sort.by="P")

plot_volcano <- function(tt, title_txt, lfc_thresh=1, fdr=0.05){
  df <- as.data.frame(tt) %>%
    mutate(sig = adj.P.Val < fdr & abs(logFC) > lfc_thresh,
           color = case_when(sig & logFC>0 ~ "Up",
                             sig & logFC<0 ~ "Down",
                             TRUE ~ "NS"))
  ggplot(df, aes(logFC, -log10(P.Value), color=color)) +
    geom_point(alpha=.6, size=1) +
    scale_color_manual(values=c(Up="red",Down="blue",NS="grey70")) +
    geom_vline(xintercept=c(-lfc_thresh, lfc_thresh), linetype="dashed") +
    geom_hline(yintercept=-log10(fdr), linetype="dashed") +
    labs(title=paste0("Volcano: ", title_txt), x="log2 Fold Change", y="-log10 P-value") +
    theme_minimal()
}

plot_volcano_labeled <- function(tt, title_txt, lfc_thresh=1, fdr=0.05, label_n=12){
  keyType <- if (grepl("^ENSMUSG", rownames(tt)[1])) "ENSEMBL" else "SYMBOL"
  keys <- rownames(tt)
  symbol_map <- AnnotationDbi::mapIds(org.Mm.eg.db, keys=keys, keytype=keyType, column="SYMBOL", multiVals="first")
  df <- as.data.frame(tt) %>%
    tibble::rownames_to_column("KEY") %>%
    mutate(gene = ifelse(is.na(symbol_map[KEY]), KEY, unname(symbol_map[KEY])),
           sig = adj.P.Val < fdr & abs(logFC) > lfc_thresh,
           color = case_when(sig & logFC>0 ~ "Up",
                             sig & logFC<0 ~ "Down",
                             TRUE ~ "NS"),
           neglog10p = -log10(P.Value))
  lab_up   <- df %>% filter(color=="Up")   %>% slice_max(logFC, n=label_n, with_ties=FALSE)
  lab_down <- df %>% filter(color=="Down") %>% slice_min(logFC, n=label_n, with_ties=FALSE)
  labs_df  <- bind_rows(lab_up, lab_down)

  ggplot(df, aes(logFC, neglog10p, color=color)) +
    geom_point(alpha=.6, size=1) +
    geom_vline(xintercept=c(-lfc_thresh, lfc_thresh), linetype="dashed") +
    geom_hline(yintercept=-log10(fdr), linetype="dashed") +
    scale_color_manual(values=c(Up="red",Down="blue",NS="grey70")) +
    ggrepel::geom_text_repel(data=labs_df, aes(label=gene),
                             size=3, max.overlaps=Inf, box.padding=.3, point.padding=.2,
                             segment.alpha=.5, min.segment.length=0) +
    labs(title=paste0("Volcano (labeled): ", title_txt),
         x="log2 Fold Change", y="-log10 P-value") +
    theme_minimal()
}

ggsave(file.path(cfg$paths$figures,"Volcano_Lo_vs_PBS.png"),  plot_volcano(tt_lo, "Lo vs PBS",  cfg$params$lfc_thresh, cfg$params$fdr_thresh), width=6, height=5, dpi=300)
ggsave(file.path(cfg$paths$figures,"Volcano_Med_vs_PBS.png"), plot_volcano(tt_med,"Med vs PBS", cfg$params$lfc_thresh, cfg$params$fdr_thresh), width=6, height=5, dpi=300)
ggsave(file.path(cfg$paths$figures,"Volcano_Hi_vs_PBS.png"),  plot_volcano(tt_hi, "Hi vs PBS",  cfg$params$lfc_thresh, cfg$params$fdr_thresh), width=6, height=5, dpi=300)

ggsave(file.path(cfg$paths$figures,"Volcano_Lo_vs_PBS_labeled.png"),
       plot_volcano_labeled(tt_lo, "Lo vs PBS", cfg$params$lfc_thresh, cfg$params$fdr_thresh), width=6, height=5, dpi=300)

png(file.path(cfg$paths$figures,"pval_histograms.png"), width=1200, height=400)
par(mfrow=c(1,3))
hist(tt_lo$P.Value,  50, main="Lo vs PBS",  col="lightblue",  xlab="P-value")
hist(tt_med$P.Value, 50, main="Med vs PBS", col="lightgreen", xlab="P-value")
hist(tt_hi$P.Value,  50, main="Hi vs PBS",  col="lightcoral", xlab="P-value")
dev.off()

png(file.path(cfg$paths$figures,"MA_Lo_Med_Hi.png"), width=1200, height=400)
par(mfrow=c(1,3))
limma::plotMA(fit2, coef="Lo_vs_PBS",  main="MA: Lo vs PBS");  abline(h=c(-1,1), col="red", lty=2)
limma::plotMA(fit2, coef="Med_vs_PBS", main="MA: Med vs PBS"); abline(h=c(-1,1), col="red", lty=2)
limma::plotMA(fit2, coef="Hi_vs_PBS",  main="MA: Hi vs PBS");  abline(h=c(-1,1), col="red", lty=2)
dev.off()

pca <- prcomp(t(v$E), scale.=TRUE)
p <- factoextra::fviz_pca_ind(
       pca,
       geom.ind = "point",
       habillage = meta$Group,   # << doğru kullanım
       addEllipses = TRUE,
       ellipse.level = 0.95,
       legend.title = "Group"
     ) + theme_minimal()
ggsave(file.path(cfg$paths$figures,"PCA_voomE.png"), p, width=6.5, height=5.5, dpi=300)
