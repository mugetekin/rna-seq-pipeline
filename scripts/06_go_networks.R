# GO emap + cnet plots
source("R/io_helpers.R")
suppressPackageStartupMessages({
  library(clusterProfiler); library(enrichplot); library(org.Mm.eg.db)
  library(ggplot2); library(ragg)
})

cfg   <- load_cfg()
g_obj <- readRDS(file.path(cfg$paths$rds, "go_objs.rds"))
g_lo <- g_obj$g_lo; g_med <- g_obj$g_med; g_hi <- g_obj$g_hi
r_lo <- g_obj$r_lo; r_med <- g_obj$r_med; r_hi <- g_obj$r_hi

dir.create(file.path(cfg$paths$figures,"go"), showWarnings=FALSE, recursive=TRUE)

make_emap <- function(g, label, show_n=25, w=9, h=7){
  if (is.null(g) || nrow(as.data.frame(g))==0) return(invisible(NULL))
  set.seed(42)
  g_sim <- enrichplot::pairwise_termsim(g)
  p <- enrichplot::emapplot(g_sim, showCategory=show_n, layout.params=list(layout="fr")) +
       ggtitle(paste0(label," — GO:BP (emap)"))
  ggsave(file.path(cfg$paths$figures,"go", paste0("emap_", label, "_BP.png")),
         plot=p, width=w, height=h, dpi=300, device=ragg::agg_png)
  invisible(p)
}

make_cnet <- function(g, label, fc_vec, show_n=15, w=10, h=7){
  if (is.null(g) || nrow(as.data.frame(g))==0) return(invisible(NULL))
  set.seed(42)
  p <- enrichplot::cnetplot(
         g,
         showCategory = show_n,
         circular     = FALSE,
         node_label   = "category",
         layout.params= list(layout="fr"),
         color.params = list(foldChange = fc_vec, edge = TRUE)  # yeni API
       ) + ggtitle(paste0(label," — gene–GO network"))
  ggsave(file.path(cfg$paths$figures,"go", paste0("cnet_", label, "_BP.png")),
         plot=p, width=w, height=h, dpi=300, device=ragg::agg_png)
  invisible(p)
}

make_emap(g_lo,  "Lo_vs_PBS")
make_emap(g_med, "Med_vs_PBS")
make_emap(g_hi,  "Hi_vs_PBS")

make_cnet(g_lo,  "Lo_vs_PBS",  r_lo)
make_cnet(g_med, "Med_vs_PBS", r_med)
make_cnet(g_hi,  "Hi_vs_PBS",  r_hi)
