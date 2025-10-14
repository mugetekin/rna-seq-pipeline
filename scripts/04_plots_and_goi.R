# GOI shortlist + generic plots (resilient even if GO results empty)
source("R/io_helpers.R")
suppressPackageStartupMessages({
  library(tidyverse); library(ggrepel); library(scales)
  library(AnnotationDbi); library(org.Mm.eg.db)
  library(emmeans); library(DOSE); library(readr)
})

`%||%` <- function(a, b) if (is.null(a) || !length(a)) b else a

cfg      <- load_cfg()
norm_obj <- readRDS(file.path(cfg$paths$rds, "norm_and_fit.rds"))
de_obj   <- readRDS(file.path(cfg$paths$rds, "de_tables.rds"))
go_path  <- file.path(cfg$paths$rds, "go_objs.rds")
go_obj   <- if (file.exists(go_path)) readRDS(go_path) else list()

annot <- norm_obj$annot; meta <- norm_obj$meta; v <- norm_obj$v
tt_lo <- de_obj$tt_lo;   tt_med <- de_obj$tt_med; tt_hi <- de_obj$tt_hi
g_lo  <- go_obj$g_lo  %||% NULL; g_med <- go_obj$g_med %||% NULL; g_hi <- go_obj$g_hi %||% NULL

dir.create(file.path(cfg$paths$results, "goi"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(cfg$paths$figures, "goi"), showWarnings = FALSE, recursive = TRUE)

# --- Helper functions ---
safe_symbol_map <- function(keys){
  keys <- as.character(keys)
  if (!length(keys)) return(character(0))
  tryCatch(
    AnnotationDbi::mapIds(org.Mm.eg.db, keys = keys, keytype = "SYMBOL",
                          column = "SYMBOL", multiVals = "first"),
    error = function(e) setNames(keys, keys)
  )
}
find_col <- function(nms, candidates) {
  cand <- candidates[candidates %in% nms]
  if (length(cand)) cand[1] else NA_character_
}

pick_de <- function(tt, lbl){
  if (is.null(tt) || !nrow(tt)) return(tibble(KEY=character(), !!paste0("logFC_", lbl):=numeric(), !!paste0("FDR_", lbl):=numeric()))
  df <- as.data.frame(tt)
  lfc_col <- find_col(names(df), c("logFC","log2FoldChange"))
  fdr_col <- find_col(names(df), c("adj.P.Val","padj","FDR","qvalue"))
  rn <- if (!is.null(rownames(df)) && all(nzchar(rownames(df)))) rownames(df) else as.character(df$Gene %||% seq_len(nrow(df)))
  tibble(KEY = rn, !!paste0("logFC_", lbl):=df[[lfc_col]], !!paste0("FDR_", lbl):=df[[fdr_col]])
}

# --- merge DEs ---
de_tbl <- list(lo=tt_lo, med=tt_med, hi=tt_hi) |>
  purrr::imap(~pick_de(.x, .y)) |>
  purrr::reduce(full_join, by="KEY")

sym_map <- safe_symbol_map(de_tbl$KEY)
de_tbl <- de_tbl |>
  mutate(SYMBOL = unname(sym_map[KEY]),
         pass_lo  = !is.na(FDR_lo)  & FDR_lo  <= (cfg$params$fdr_thresh %||% 0.05) & abs(logFC_lo)  >= (cfg$params$lfc_thresh %||% 1),
         pass_med = !is.na(FDR_med) & FDR_med <= (cfg$params$fdr_thresh %||% 0.05) & abs(logFC_med) >= (cfg$params$lfc_thresh %||% 1),
         pass_hi  = !is.na(FDR_hi)  & FDR_hi  <= (cfg$params$fdr_thresh %||% 0.05) & abs(logFC_hi)  >= (cfg$params$lfc_thresh %||% 1),
         pass_any = pass_lo | pass_med | pass_hi,
         mean_abs_logFC = rowMeans(cbind(abs(logFC_lo), abs(logFC_med), abs(logFC_hi)), na.rm = TRUE),
         dir_consistent = {
           s <- c(logFC_lo, logFC_med, logFC_hi); non_na <- s[!is.na(s)]
           length(non_na) >= 2 && all(sign(non_na) == sign(non_na[1]))
         })

# --- extract GO focus ---
extract_symbols_from_gsea <- function(g, n_terms = 12){
  df <- as.data.frame(g)
  if (is.null(df) || !nrow(df)) return(character(0))
  rank_col <- if ("qvalues" %in% names(df)) "qvalues" else if ("p.adjust" %in% names(df)) "p.adjust" else names(df)[1]
  df <- df |>
    arrange(.data[[rank_col]], desc(abs(.data$NES))) |>
    slice_head(n = min(n_terms, nrow(df)))
  ent <- unique(unlist(strsplit(paste(df$core_enrichment, collapse="/"), "/")))
  if (!length(ent)) return(character(0))
  sy <- AnnotationDbi::mapIds(org.Mm.eg.db, keys=ent, keytype="ENTREZID", column="SYMBOL", multiVals="first")
  unique(unname(sy[!is.na(sy)]))
}
go_focus <- unique(c(
  extract_symbols_from_gsea(g_lo, 12),
  extract_symbols_from_gsea(g_med,12),
  extract_symbols_from_gsea(g_hi, 12)
))

# --- Fallback: if GO empty, fill with top DE genes ---
if (length(go_focus) == 0) {
  message("No GO terms found; falling back to top DE genes for GOI focus.")
  go_focus <- de_tbl |>
    arrange(desc(mean_abs_logFC)) |>
    slice_head(n=300) |>
    pull(SYMBOL) |>
    unique()
}

# --- rank and shortlist ---
rng <- range(de_tbl$mean_abs_logFC, na.rm = TRUE); if (diff(rng)==0) rng <- c(0,1)
de_ranked <- de_tbl |>
  mutate(in_GO_focus = SYMBOL %in% go_focus,
         score = 2*in_GO_focus + 1*dir_consistent +
                 scales::rescale(mean_abs_logFC, to=c(0,1), from=rng)) |>
  arrange(desc(score), desc(mean_abs_logFC))

# always produce populated shortlists
short_strict <- de_ranked |> filter(pass_any | in_GO_focus) |> slice_head(n = 150)
short_broad  <- de_ranked |> filter(in_GO_focus | dir_consistent | pass_any) |> slice_head(n = 300)

write_csv(de_ranked,    file.path(cfg$paths$results, "goi", "GOI_ranked_all.csv"))
write_csv(short_strict, file.path(cfg$paths$results, "goi", "GOI_shortlist_strict.csv"))
write_csv(short_broad,  file.path(cfg$paths$results, "goi", "GOI_shortlist_broad.csv"))

message(sprintf("Generated GOI shortlists: strict=%d genes, broad=%d genes",
                nrow(short_strict), nrow(short_broad)))
