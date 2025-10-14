# 05_qc_and_volcano.R — QC, volcano/MA/PCA/histograms (dynamic coef names)
source("R/io_helpers.R")
suppressPackageStartupMessages({
  library(edgeR); library(limma); library(tidyverse)
  library(ggrepel); library(factoextra); library(ragg)
  library(AnnotationDbi); library(org.Mm.eg.db)
})

`%||%` <- function(a, b) if (is.null(a) || !length(a)) b else a

cfg <- load_cfg()
obj <- readRDS(file.path(cfg$paths$rds,"norm_and_fit.rds"))
annot <- obj$annot; meta <- obj$meta; dge <- obj$dge; v <- obj$v; fit2 <- obj$fit2

dir.create(cfg$paths$figures, showWarnings = FALSE, recursive = TRUE)

# --- coefficient name helpers: map "Lo - PBS" -> "Lo_vs_PBS" ---
norm_name   <- function(x) gsub(" - ", "_vs_", x, fixed = TRUE)
pretty_name <- function(x) gsub("_vs_", " vs ", x, fixed = TRUE)
avail_raw   <- colnames(fit2$coefficients)
avail_norm  <- norm_name(avail_raw)

get_tt <- function(cn_norm) {
  idx <- match(cn_norm, avail_norm)
  if (is.na(idx)) return(NULL)
  topTable(fit2, coef = avail_raw[idx], number = Inf, sort.by = "P")
}

# --- 1) Library sizes ---
png(file.path(cfg$paths$figures,"library_sizes.png"), width=900, height=600)
barplot(colSums(dge$counts), las=2, main="Library sizes", ylab="Total counts")
dev.off()

# --- 2) MDS (TMM-normalized) ---
png(file.path(cfg$paths$figures,"MDS_groups.png"), width=700, height=700)
plotMDS(dge, labels = meta$Group); title("MDS (TMM-normalized)")
dev.off()

# --- 3) voom mean–variance ---
pdf(file.path(cfg$paths$figures,"voom_meanvar.pdf"), 6, 5)
voom(dge, model.matrix(~0+Group, data=meta), plot=TRUE)
dev.off()

# --- thresholds ---
lfc_thr <- cfg$params$lfc_thresh %||% 1
fdr_thr <- cfg$params$fdr_thresh %||% 0.05

# --- helpers ---
plot_volcano <- function(tt, title_txt, lfc_thresh=lfc_thr, fdr=fdr_thr){
  req <- c("logFC","P.Value","adj.P.Val"); if (is.null(tt) || !all(req %in% names(tt))) return(NULL)
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

plot_volcano_labeled <- function(tt, title_txt, lfc_thresh=lfc_thr, fdr=fdr_thr, label_n=12){
  if (is.null(tt) || !nrow(tt)) return(NULL)
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

# --- 4) Decide which contrasts to plot ---
want <- cfg$params$contrasts %||% avail_norm   # if config lists some, respect it
to_plot <- intersect(want, avail_norm)
if (!length(to_plot)) {
  message("No contrasts to plot; available: ", paste(avail_norm, collapse=", "))
}

# --- 5) Per-contrast volcanoes / labeled volcanoes / MA / p-value hists ---
for (cn in to_plot) {
  tt <- get_tt(cn)
  title_txt <- pretty_name(cn)

  # Volcano
  pv <- plot_volcano(tt, title_txt)
  if (!is.null(pv))
    ggsave(file.path(cfg$paths$figures, paste0("Volcano_", cn, ".png")), pv, width=6, height=5, dpi=300)

  # Labeled Volcano
  pvl <- plot_volcano_labeled(tt, title_txt)
  if (!is.null(pvl))
    ggsave(file.path(cfg$paths$figures, paste0("Volcano_", cn, "_labeled.png")), pvl, width=6, height=5, dpi=300)

  # MA plot
  png(file.path(cfg$paths$figures, paste0("MA_", cn, ".png")), width=700, height=550)
  idx <- match(cn, avail_norm)
  if (!is.na(idx)) { limma::plotMA(fit2, coef=avail_raw[idx], main=paste0("MA: ", title_txt)); abline(h=c(-1,1), col="red", lty=2) } else plot.new()
  dev.off()

  # P-value histogram
  png(file.path(cfg$paths$figures, paste0("PvalHist_", cn, ".png")), width=700, height=500)
  if (!is.null(tt) && "P.Value" %in% names(tt)) hist(tt$P.Value, 50, main=title_txt, col="grey80", xlab="P-value")
  dev.off()
}

# --- 6) Classic Lo/Med/Hi panels (kept for continuity if present) ---
tt_lo  <- get_tt("Lo_vs_PBS")
tt_med <- get_tt("Med_vs_PBS")
tt_hi  <- get_tt("Hi_vs_PBS")

if (!is.null(tt_lo))
  ggsave(file.path(cfg$paths$figures,"Volcano_Lo_vs_PBS.png"),
         plot_volcano(tt_lo, "Lo vs PBS"), width=6, height=5, dpi=300)
if (!is.null(tt_med))
  ggsave(file.path(cfg$paths$figures,"Volcano_Med_vs_PBS.png"),
         plot_volcano(tt_med,"Med vs PBS"), width=6, height=5, dpi=300)
if (!is.null(tt_hi))
  ggsave(file.path(cfg$paths$figures,"Volcano_Hi_vs_PBS.png"),
         plot_volcano(tt_hi, "Hi vs PBS"), width=6, height=5, dpi=300)

pv <- plot_volcano_labeled(tt_lo, "Lo vs PBS")
if (!is.null(pv))
  ggsave(file.path(cfg$paths$figures,"Volcano_Lo_vs_PBS_labeled.png"), pv, width=6, height=5, dpi=300)

png(file.path(cfg$paths$figures,"pval_histograms.png"), width=1200, height=400)
par(mfrow=c(1,3))
if (!is.null(tt_lo))  hist(tt_lo$P.Value,  50, main="Lo vs PBS",  col="lightblue",  xlab="P-value") else plot.new()
if (!is.null(tt_med)) hist(tt_med$P.Value, 50, main="Med vs PBS", col="lightgreen", xlab="P-value") else plot.new()
if (!is.null(tt_hi))  hist(tt_hi$P.Value,  50, main="Hi vs PBS",  col="lightcoral", xlab="P-value") else plot.new()
dev.off()

png(file.path(cfg$paths$figures,"MA_Lo_Med_Hi.png"), width=1200, height=400)
par(mfrow=c(1,3))
if (!is.null(tt_lo))  { limma::plotMA(fit2, coef=avail_raw[match("Lo_vs_PBS", avail_norm)],  main="MA: Lo vs PBS");  abline(h=c(-1,1), col="red", lty=2) } else plot.new()
if (!is.null(tt_med)) { limma::plotMA(fit2, coef=avail_raw[match("Med_vs_PBS",avail_norm)], main="MA: Med vs PBS"); abline(h=c(-1,1), col="red", lty=2) } else plot.new()
if (!is.null(tt_hi))  { limma::plotMA(fit2, coef=avail_raw[match("Hi_vs_PBS", avail_norm)],  main="MA: Hi vs PBS");  abline(h=c(-1,1), col="red", lty=2) } else plot.new()
dev.off()

# --- 7) PCA (voom E) ---
pca <- prcomp(t(v$E), scale.=TRUE)
p <- factoextra::fviz_pca_ind(
       pca,
       geom.ind = "point",
       habillage = meta$Group,
       addEllipses = TRUE,
       ellipse.level = 0.95,
       legend.title = "Group"
     ) + theme_minimal()
ggsave(file.path(cfg$paths$figures,"PCA_voomE.png"), p, width=6.5, height=5.5, dpi=300)
