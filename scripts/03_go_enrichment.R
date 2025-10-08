source("R/io_helpers.R")
suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(ragg)
  library(tibble)
})

cfg  <- load_cfg()
obj  <- readRDS(file.path(cfg$paths$rds, "de_tables.rds"))

# ---- Build ranked vectors from limma topTable (SYMBOL or ENSEMBL rownames) ----
make_rank_from_tt <- function(tt, fc_col = "logFC") {
  stopifnot(fc_col %in% names(tt))
  keyType <- if (grepl("^ENSMUSG", rownames(tt)[1])) "ENSEMBL" else "SYMBOL"
  ids <- AnnotationDbi::select(org.Mm.eg.db, keys = rownames(tt),
                               keytype = keyType, columns = "ENTREZID")
  ids <- ids[!is.na(ids$ENTREZID), c(keyType, "ENTREZID")]
  colnames(ids)[1] <- "KEY"
  df <- tibble(KEY = rownames(tt), FC = tt[[fc_col]])
  m  <- merge(df, ids, by = "KEY", all.x = TRUE)
  agg <- stats::aggregate(FC ~ ENTREZID, m, function(x) mean(x, na.rm = TRUE))
  v <- setNames(agg$FC, agg$ENTREZID)
  v <- v[is.finite(v)]
  sort(v, decreasing = TRUE)
}

r_lo  <- make_rank_from_tt(obj$tt_lo)
r_med <- make_rank_from_tt(obj$tt_med)
r_hi  <- make_rank_from_tt(obj$tt_hi)

# ---- gseGO (BP) ----
if (requireNamespace("BiocParallel", quietly = TRUE)) {
  BiocParallel::register(BiocParallel::SerialParam())
}
set.seed(42)
g_lo  <- gseGO(geneList = r_lo,  OrgDb = org.Mm.eg.db, ont = cfg$params$go_ont, eps = 0, verbose = FALSE)
g_med <- gseGO(geneList = r_med, OrgDb = org.Mm.eg.db, ont = cfg$params$go_ont, eps = 0, verbose = FALSE)
g_hi  <- gseGO(geneList = r_hi,  OrgDb = org.Mm.eg.db, ont = cfg$params$go_ont, eps = 0, verbose = FALSE)

# ---- write tables ----
dir.create(file.path(cfg$paths$results, "go"), showWarnings = FALSE, recursive = TRUE)
readr::write_csv(as.data.frame(g_lo),  file.path(cfg$paths$results, "go", "gseGO_Lo.csv"))
readr::write_csv(as.data.frame(g_med), file.path(cfg$paths$results, "go", "gseGO_Med.csv"))
readr::write_csv(as.data.frame(g_hi),  file.path(cfg$paths$results, "go", "gseGO_Hi.csv"))

# ---- dotplots (split by .sign) ----
dir.create(file.path(cfg$paths$figures, "go"), showWarnings = FALSE, recursive = TRUE)
save_plot <- function(p, name, w=8, h=5){
  ggplot2::ggsave(file.path(cfg$paths$figures, "go", paste0(name, ".png")),
                  plot = p, width = w, height = h, dpi = 300, device = ragg::agg_png)
}

p_lo  <- enrichplot::dotplot(g_lo,  showCategory = cfg$params$go_n_show, split=".sign") + ggtitle("Lo vs PBS — GO:BP")
p_med <- enrichplot::dotplot(g_med, showCategory = cfg$params$go_n_show, split=".sign") + ggtitle("Med vs PBS — GO:BP")
p_hi  <- enrichplot::dotplot(g_hi,  showCategory = cfg$params$go_n_show, split=".sign") + ggtitle("Hi vs PBS — GO:BP")

save_plot(p_lo,  "dotplot_gseGO_Lo_BP")
save_plot(p_med, "dotplot_gseGO_Med_BP")
save_plot(p_hi,  "dotplot_gseGO_Hi_BP")

# ---- Save R objects for downstream (GOI shortlist step) ----
saveRDS(list(g_lo=g_lo, g_med=g_med, g_hi=g_hi, r_lo=r_lo, r_med=r_med, r_hi=r_hi),
        file.path(cfg$paths$rds, "go_objs.rds"))
