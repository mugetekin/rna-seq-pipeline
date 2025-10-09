# 03_go_enrichment.R — Gene Ontology GSEA (voom+limma ranked genes)

source("R/io_helpers.R")
suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(ragg)
  library(tibble)
  library(dplyr)
  library(readr)
})

cfg  <- load_cfg()
obj  <- readRDS(file.path(cfg$paths$rds, "de_tables.rds"))

# --- 1) Prepare ranked gene vectors (logFC) from limma results ---
make_rank_from_tt <- function(tt, fc_col = "logFC") {
  stopifnot(fc_col %in% names(tt))
  keyType <- if (grepl("^ENSMUSG", rownames(tt)[1])) "ENSEMBL" else "SYMBOL"
  ids <- AnnotationDbi::select(org.Mm.eg.db, keys = rownames(tt),
                               keytype = keyType, columns = "ENTREZID")
  ids <- ids[!is.na(ids$ENTREZID), c(keyType, "ENTREZID")]
  colnames(ids)[1] <- "KEY"
  df <- tibble(KEY = rownames(tt), FC = tt[[fc_col]])
  m  <- merge(df, ids, by = "KEY", all.x = TRUE)
  agg <- stats::aggregate(FC ~ ENTREZID, m, mean, na.rm = TRUE)
  v <- setNames(agg$FC, agg$ENTREZID)
  v <- v[is.finite(v)]
  sort(v, decreasing = TRUE)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

tt_lo  <- obj$tt_lo  %||% obj$Lo_vs_PBS
tt_med <- obj$tt_med %||% obj$Med_vs_PBS
tt_hi  <- obj$tt_hi  %||% obj$Hi_vs_PBS

r_lo  <- make_rank_from_tt(tt_lo)
r_med <- make_rank_from_tt(tt_med)
r_hi  <- make_rank_from_tt(tt_hi)

# --- 2) Safe wrapper for GSEA ---
set.seed(42)
gse_safe <- function(rnk) {
  if (is.null(rnk) || length(rnk) < 50) return(NULL)
  gseGO(
    geneList = rnk,
    OrgDb = org.Mm.eg.db,
    ont = cfg$params$go_ont,
    minGSSize = 10, maxGSSize = 500,
    pAdjustMethod = "BH",
    eps = 0,
    verbose = FALSE,
    seed = TRUE
  )
}

g_lo  <- gse_safe(r_lo)
g_med <- gse_safe(r_med)
g_hi  <- gse_safe(r_hi)

# --- 3) Export raw results ---
dir.create(file.path(cfg$paths$results, "go"), showWarnings = FALSE, recursive = TRUE)
if (!is.null(g_lo))  readr::write_csv(as.data.frame(g_lo),  file.path(cfg$paths$results, "go", "gseGO_Lo.csv"))
if (!is.null(g_med)) readr::write_csv(as.data.frame(g_med), file.path(cfg$paths$results, "go", "gseGO_Med.csv"))
if (!is.null(g_hi))  readr::write_csv(as.data.frame(g_hi),  file.path(cfg$paths$results, "go", "gseGO_Hi.csv"))

# --- 4) Filter for visualization ---
flt <- function(g) {
  if (is.null(g)) return(NULL)
  df <- as.data.frame(g)
  if (!nrow(df)) return(NULL)
  df <- df %>%
    mutate(size_ok = between(as.numeric(setSize), 10, 500)) %>%
    filter(size_ok, p.adjust < 0.05)
  if (!nrow(df)) return(NULL)
  g@result <- df
  g
}

g_lo_p  <- flt(g_lo)
g_med_p <- flt(g_med)
g_hi_p  <- flt(g_hi)

# --- 5) Dotplots ---
dir.create(file.path(cfg$paths$figures, "go"), showWarnings = FALSE, recursive = TRUE)
save_plot <- function(p, name, w=8, h=5){
  ggplot2::ggsave(file.path(cfg$paths$figures, "go", paste0(name, ".png")),
                  plot = p, width = w, height = h, dpi = 300, device = ragg::agg_png)
}

dplot <- function(g, title_txt) {
  if (is.null(g)) return(NULL)
  suppressWarnings(
    enrichplot::dotplot(g, showCategory = cfg$params$go_n_show, split = ".sign") +
      ggtitle(title_txt)
  )
}

if (!is.null(g_lo_p))  save_plot(dplot(g_lo_p,  "Lo vs PBS — GO:BP"),  "dotplot_gseGO_Lo_BP")
if (!is.null(g_med_p)) save_plot(dplot(g_med_p, "Med vs PBS — GO:BP"), "dotplot_gseGO_Med_BP")
if (!is.null(g_hi_p))  save_plot(dplot(g_hi_p,  "Hi vs PBS — GO:BP"),  "dotplot_gseGO_Hi_BP")

# --- 6) Save objects for downstream steps ---
saveRDS(
  list(g_lo = g_lo, g_med = g_med, g_hi = g_hi,
       r_lo = r_lo, r_med = r_med, r_hi = r_hi),
  file.path(cfg$paths$rds, "go_objs.rds")
)
