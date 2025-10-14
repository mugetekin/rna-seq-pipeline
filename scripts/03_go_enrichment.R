# 03_go_enrichment.R — GO GSEA + ORA with robust ID handling (mouse)
source("R/io_helpers.R")
suppressPackageStartupMessages({
  library(tidyverse)
  library(AnnotationDbi); library(org.Mm.eg.db)
  library(clusterProfiler); library(DOSE)
  library(readr)
})

`%||%` <- function(a, b) if (is.null(a) || !length(a)) b else a

cfg <- load_cfg()
rds_dir <- cfg$paths$rds
res_dir <- cfg$paths$results
dir.create(rds_dir, showWarnings = FALSE, recursive = TRUE)

# --- inputs ---
obj    <- readRDS(file.path(rds_dir, "norm_and_fit.rds"))
de_obj <- readRDS(file.path(rds_dir, "de_tables.rds"))
tt_lo  <- de_obj$tt_lo
tt_med <- de_obj$tt_med
tt_hi  <- de_obj$tt_hi

# --- helpers ---
strip_ver <- function(x) sub("\\.\\d+$", "", x)  # ENSMUSG00000000001.1 -> ENSMUSG00000000001

safe_map <- function(keys, keytype){
  # Return a vector of same length as keys; NAs where unmapped or on error.
  out <- tryCatch(
    AnnotationDbi::mapIds(org.Mm.eg.db, keys = keys, keytype = keytype,
                          column = "ENTREZID", multiVals = "first"),
    error = function(e) {
      # when *none* are valid for this keytype, mapIds errors; return all NA
      rep(NA_character_, length(keys))
    }
  )
  # mapIds may return a named vector shorter than input; align to input order
  if (length(out) != length(keys) || is.null(names(out))) {
    nm <- names(out)
    aligned <- rep(NA_character_, length(keys))
    names(aligned) <- keys
    if (!is.null(nm)) aligned[nm] <- out
    unname(aligned)
  } else {
    unname(out)
  }
}

map_to_entrez_any <- function(keys){
  if (!length(keys)) return(character(0))
  # Try multiple keytypes safely; prefer SYMBOL/ENSEMBL/ALIAS/ENTREZID
  tries <- c("SYMBOL", "ENSEMBL", "ALIAS", "ENTREZID")
  candidates <- list(
    SYMBOL   = keys,
    ENSEMBL  = strip_ver(keys),
    ALIAS    = keys,
    ENTREZID = keys
  )
  results <- lapply(tries, function(kt) safe_map(candidates[[kt]], kt))
  names(results) <- tries
  success_counts <- sapply(results, function(v) sum(!is.na(v) & nzchar(v)))
  best <- names(which.max(success_counts))
  message(sprintf("ID mapping summary — best=%s  counts: %s",
                  best, paste(sprintf("%s=%d", names(success_counts), success_counts), collapse=", ")))
  results[[best]]
}

make_rank <- function(tt){
  if (is.null(tt) || !nrow(tt)) return(NULL)
  s  <- as.data.frame(tt)
  rn <- rownames(s)
  eg <- map_to_entrez_any(rn)
  metric <- if ("t" %in% colnames(s)) s$t else sign(s$logFC) * -log10(s$P.Value %||% s$adj.P.Val)
  names(metric) <- eg
  metric <- metric[!is.na(names(metric)) & names(metric) != ""]
  if (!length(metric)) return(NULL)
  tapply(metric, INDEX = names(metric), FUN = mean) |> sort(decreasing = TRUE)
}

deg_vector <- function(tt, lfc = cfg$params$lfc_thresh %||% 1, fdr = cfg$params$fdr_thresh %||% 0.05){
  if (is.null(tt) || !nrow(tt)) return(character(0))
  s  <- as.data.frame(tt)
  rn <- rownames(s)
  sig <- s$adj.P.Val <= fdr & abs(s$logFC) >= lfc
  eg  <- map_to_entrez_any(rn[sig])
  unique(eg[!is.na(eg) & nzchar(eg)])
}

# Universe for ORA: all genes in voom object that can map to ENTREZ
univ_entrez <- {
  rn <- rownames(obj$v$E)
  eg <- map_to_entrez_any(rn)
  unique(eg[!is.na(eg) & nzchar(eg)])
}

run_gsea_go <- function(rank_vec, ont = cfg$params$go_ont %||% "BP"){
  if (is.null(rank_vec) || !length(rank_vec)) return(NULL)
  suppressMessages(
    tryCatch(
      clusterProfiler::gseGO(
        geneList     = rank_vec,
        OrgDb        = org.Mm.eg.db,
        ont          = ont,
        keyType      = "ENTREZID",
        minGSSize    = 10,
        maxGSSize    = 2000,
        pvalueCutoff = 1,
        verbose      = FALSE
      ),
      error = function(e) { message("gseGO error: ", e$message); NULL }
    )
  )
}

run_ora_go <- function(gene_ids, universe_ids = NULL, ont = cfg$params$go_ont %||% "BP"){
  if (!length(gene_ids)) return(NULL)
  suppressMessages(
    tryCatch(
      clusterProfiler::enrichGO(
        gene          = gene_ids,
        universe      = universe_ids,
        OrgDb         = org.Mm.eg.db,
        keyType       = "ENTREZID",
        ont           = ont,
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.2,
        readable      = TRUE
      ),
      error = function(e) { message("enrichGO error: ", e$message); NULL }
    )
  )
}

# --- build ranked lists and DEG sets ---
rank_lo  <- make_rank(tt_lo)
rank_med <- make_rank(tt_med)
rank_hi  <- make_rank(tt_hi)

deg_lo  <- deg_vector(tt_lo)
deg_med <- deg_vector(tt_med)
deg_hi  <- deg_vector(tt_hi)

# --- run GSEA/ORA ---
g_lo  <- run_gsea_go(rank_lo)
g_med <- run_gsea_go(rank_med)
g_hi  <- run_gsea_go(rank_hi)

r_lo  <- run_ora_go(deg_lo,  universe_ids = univ_entrez)
r_med <- run_ora_go(deg_med, universe_ids = univ_entrez)
r_hi  <- run_ora_go(deg_hi,  universe_ids = univ_entrez)

# --- save ---
go_obj <- list(g_lo=g_lo, g_med=g_med, g_hi=g_hi, r_lo=r_lo, r_med=r_med, r_hi=r_hi)
saveRDS(go_obj, file.path(rds_dir, "go_objs.rds"))

# optional: export top tables
outdir <- file.path(res_dir, "go")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

save_if <- function(obj, nm){
  if (is.null(obj)) return(invisible())
  df <- as.data.frame(obj)
  if (!is.null(df) && nrow(df)) readr::write_tsv(df, file.path(outdir, paste0(nm, ".tsv")))
}

save_if(g_lo,  "GO_GSEA_lo")
save_if(g_med, "GO_GSEA_med")
save_if(g_hi,  "GO_GSEA_hi")
save_if(r_lo,  "GO_ORA_lo")
save_if(r_med, "GO_ORA_med")
save_if(r_hi,  "GO_ORA_hi")

message("GO objects saved to ", file.path(rds_dir, "go_objs.rds"))
message("GO result tables (if any) in ", outdir)
