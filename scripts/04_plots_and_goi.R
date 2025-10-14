# GOI shortlist + generic plots
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
go_obj   <- readRDS(file.path(cfg$paths$rds, "go_objs.rds"))

annot <- norm_obj$annot; meta <- norm_obj$meta; v <- norm_obj$v
tt_lo <- de_obj$tt_lo;   tt_med <- de_obj$tt_med; tt_hi <- de_obj$tt_hi
g_lo  <- go_obj$g_lo;    g_med  <- go_obj$g_med;  g_hi  <- go_obj$g_hi
r_lo  <- go_obj$r_lo;    r_med  <- go_obj$r_med;  r_hi  <- go_obj$r_hi

dir.create(file.path(cfg$paths$results, "goi"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(cfg$paths$figures, "goi"), showWarnings = FALSE, recursive = TRUE)

# -------- helpers --------
find_col <- function(nms, candidates) {
  cand <- candidates[candidates %in% nms]
  if (length(cand)) cand[1] else NA_character_
}

# Safe SYMBOL mapper: if keys look numeric or empty, just return them
safe_symbol_map <- function(keys){
  keys <- as.character(keys)
  if (!length(keys)) return(character(0))
  # If most keys are purely numeric, skip mapping
  frac_numeric <- mean(grepl("^\\d+$", keys))
  if (is.na(frac_numeric)) frac_numeric <- 0
  if (frac_numeric > 0.5) {
    names(keys) <- keys
    return(keys)
  }
  # try SYMBOL first, then ENSEMBL, but do not error
  try_map <- function(keytype){
    tryCatch(
      AnnotationDbi::mapIds(org.Mm.eg.db, keys = keys, keytype = keytype,
                            column = "SYMBOL", multiVals = "first"),
      error = function(e) setNames(rep(NA_character_, length(keys)), keys)
    )
  }
  m1 <- try_map("SYMBOL")
  # if nothing mapped, try ENSEMBL
  if (all(is.na(m1))) {
    keys2 <- sub("\\.\\d+$","", keys)
    m2 <- try_map("ENSEMBL")
    if (!is.null(names(m2))) names(m2) <- keys
    sy <- unname(ifelse(is.na(m2), keys, m2))
  } else {
    sy <- unname(ifelse(is.na(m1) | m1=="", keys, m1))
  }
  names(sy) <- keys
  sy
}

pick_de <- function(tt, lbl){
  if (is.null(tt) || !nrow(tt)) {
    return(tibble(KEY=character(),
                  !!paste0("logFC_", lbl) := numeric(),
                  !!paste0("FDR_",   lbl) := numeric()))
  }
  df <- as.data.frame(tt)

  lfc_col <- find_col(names(df), c("logFC","log2FoldChange","LFC","coef"))
  fdr_col <- find_col(names(df), c("adj.P.Val","adj.P.Val.","FDR","qvalue","padj","q.val","qvalueBH"))
  if (is.na(lfc_col) || is.na(fdr_col)) {
    stop(sprintf("Could not find LFC/FDR columns for '%s'. Available columns: %s",
                 lbl, paste(names(df), collapse=", ")))
  }

  # Prefer rownames if they exist & are not plain indices; else Gene column if present
  rn <- rownames(df)
  use_rn <- !is.null(rn) && length(rn) == nrow(df) &&
            mean(grepl("^\\d+$", rn)) < 0.5 && all(nzchar(rn))
  if (use_rn) {
    key <- rn
  } else if ("Gene" %in% names(df)) {
    key <- as.character(df$Gene)
  } else {
    key <- as.character(seq_len(nrow(df)))
  }

  tibble(KEY = key) |>
    mutate(
      !!paste0("logFC_", lbl) := as.numeric(df[[lfc_col]]),
      !!paste0("FDR_",   lbl) := as.numeric(df[[fdr_col]])
    )
}

row_z <- function(m){
  t(apply(m, 1, function(x){
    mu <- mean(x, na.rm = TRUE); sdv <- sd(x, na.rm = TRUE)
    if (!is.finite(sdv) || sdv == 0) return(rep(0, length(x)))
    (x - mu) / sdv
  }))
}

# -------- unified DE table + flags --------
de_tbl <- list(lo=tt_lo, med=tt_med, hi=tt_hi) |>
  purrr::imap(~pick_de(.x, .y)) |>
  purrr::reduce(full_join, by="KEY")

sym_map <- safe_symbol_map(de_tbl$KEY)
de_tbl  <- de_tbl |>
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

# -------- GSEA leading-edge → SYMBOL list --------
extract_symbols_from_gsea <- function(g, n_terms = 12){
  df <- as.data.frame(g)
  if (is.null(df) || !nrow(df)) return(character(0))
  rank_col <- if ("qvalues" %in% names(df)) "qvalues" else if ("p.adjust" %in% names(df)) "p.adjust" else names(df)[1]
  df <- df |>
    dplyr::arrange(.data[[rank_col]], dplyr::desc(abs(.data$NES))) |>
    dplyr::slice_head(n = min(n_terms, nrow(df)))
  ent <- unique(unlist(strsplit(paste(df$core_enrichment, collapse = "/"), "/")))
  ent <- ent[ent != "" & !is.na(ent)]
  if (!length(ent)) return(character(0))
  sy <- AnnotationDbi::mapIds(org.Mm.eg.db, keys = ent, keytype = "ENTREZID",
                              column = "SYMBOL", multiVals = "first")
  unique(unname(sy[!is.na(sy)]))
}
go_focus <- unique(c(
  extract_symbols_from_gsea(g_lo, 12),
  extract_symbols_from_gsea(g_med,12),
  extract_symbols_from_gsea(g_hi, 12)
))

# -------- score & shortlists --------
rng <- range(de_tbl$mean_abs_logFC, na.rm = TRUE); if (diff(rng)==0) rng <- c(0,1)
de_ranked <- de_tbl |>
  mutate(in_GO_focus = SYMBOL %in% go_focus,
         score = 2*in_GO_focus + 1*dir_consistent +
                 scales::rescale(mean_abs_logFC, to=c(0,1), from=rng)) |>
  arrange(desc(score), desc(mean_abs_logFC))

short_strict <- de_ranked |> filter(pass_any & in_GO_focus) |> slice_head(n = 150)
short_broad  <- de_ranked |> filter(in_GO_focus | dir_consistent) |> slice_head(n = 300)

readr::write_csv(de_ranked,   file.path(cfg$paths$results, "goi", "GOI_ranked_all.csv"))
readr::write_csv(short_strict,file.path(cfg$paths$results, "goi", "GOI_shortlist_strict.csv"))
readr::write_csv(short_broad, file.path(cfg$paths$results, "goi", "GOI_shortlist_broad.csv"))

# -------- plotting helpers --------
plot_violin_goi <- function(goi, fname){
  sy_map <- safe_symbol_map(rownames(v$E))
  keep <- tolower(unname(sy_map)) %in% tolower(goi)
  stopifnot(any(keep))
  mat <- v$E[keep,, drop=FALSE]; rownames(mat) <- unname(sy_map[rownames(mat)])
  grp <- meta$Group[match(colnames(mat), meta$SampleID)]
  grp <- factor(grp, levels=c("PBS","Lo","Med","Hi"))
  df <- as.data.frame(mat) |>
    rownames_to_column("Gene") |>
    pivot_longer(-Gene, names_to="Sample", values_to="logCPM") |>
    mutate(Group = grp[match(Sample, colnames(mat))]) |>
    filter(!is.na(Group))
  p <- ggplot(df, aes(Group, logCPM, fill=Group)) +
    geom_violin(trim=FALSE, alpha=.6, color="black") +
    geom_jitter(width=.12, size=2, alpha=.9, shape=21, color="black") +
    stat_summary(fun=mean, geom="point", shape=23, size=3, fill="white") +
    facet_wrap(~Gene, scales="free_y", nrow=1) +
    labs(title="GOI — voom log2-CPM (violin)", x=NULL, y="log2(CPM+1)") +
    theme_minimal(base_size=12) + theme(legend.position = "none")
  ggsave(fname, p, width = max(7, 3 + 2.5*length(unique(df$Gene))), height=4, dpi=300)
}

plot_heatmap_goi <- function(goi, fname){
  sy_map <- safe_symbol_map(rownames(v$E))
  keep <- tolower(unname(sy_map)) %in% tolower(goi)
  stopifnot(any(keep))
  mat <- v$E[keep,, drop=FALSE]; rownames(mat) <- unname(sy_map[rownames(mat)])
  mat_z <- row_z(mat)
  sid <- colnames(mat); grp <- case_when(
    startsWith(sid,"PBS")~"PBS", startsWith(sid,"Lo")~"Lo",
    startsWith(sid,"Med")~"Med", startsWith(sid,"Hi")~"Hi", TRUE~"Other")
  grp <- factor(grp, levels=c("PBS","Lo","Med","Hi","Other"))
  order_within_group <- function(M, idx){
    if (length(idx) <= 2) return(idx)
    Msub <- M[, idx, drop=FALSE]
    Msub_imp <- apply(Msub, 2, function(col){ if (all(!is.finite(col))) return(rep(0,length(col)))
      col[!is.finite(col)] <- median(col[is.finite(col)], na.rm=TRUE); col })
    idx[ hclust(dist(t(Msub_imp)), method="average")$order ]
  }
  ord <- unlist(lapply(levels(grp), function(g){
    idx <- which(grp==g); if (length(idx)) order_within_group(mat_z, idx) else integer(0)
  }), use.names=FALSE)
  mat_z <- mat_z[, ord, drop=FALSE]; sid <- sid[ord]; grp <- droplevels(grp[ord])

  df <- as.data.frame(mat_z) |>
    rownames_to_column("Gene") |>
    pivot_longer(-Gene, names_to="Sample", values_to="z") |>
    mutate(Sample=factor(Sample, levels=sid), Gene=factor(Gene, levels=rev(rownames(mat_z))))
  pal <- scales::gradient_n_pal(c("#214478","#f7f7f7","#b30000"))
  p <- ggplot(df, aes(Sample, Gene, fill=z)) +
    geom_tile() +
    scale_fill_gradientn(colors = pal(seq(0,1,length.out=101)), limits=c(-2.5,2.5), oob=squish, name="row Z") +
    labs(title="GOI — row Z heatmap (voom log2-CPM)", x=NULL, y=NULL) +
    theme_minimal(base_size=12) +
    theme(panel.grid=element_blank(), axis.text.x=element_text(angle=90, vjust=.5, hjust=1))
  ggsave(fname, p, width = max(8, 1.2*ncol(mat_z)), height = max(3.5, 0.35*nrow(mat_z)+1), dpi=200)
}

plot_barjit_dunnett_goi <- function(goi, fname){
  sy_map <- safe_symbol_map(rownames(v$E))
  keep <- tolower(unname(sy_map)) %in% tolower(goi)
  stopifnot(any(keep))
  mat <- v$E[keep,, drop=FALSE]; rownames(mat) <- unname(sy_map[rownames(mat)])
  df <- as.data.frame(mat) |>
    rownames_to_column("Gene") |>
    pivot_longer(-Gene, names_to="Sample", values_to="logCPM") |>
    mutate(Group = meta$Group[match(Sample, meta$SampleID)]) |>
    filter(!is.na(Group)) |>
    mutate(Group=factor(Group, levels=c("PBS","Lo","Med","Hi")))

  sem <- function(x){ x <- x[is.finite(x)]; if (!length(x)) NA_real_ else sd(x)/sqrt(length(x)) }
  get_stars <- function(p) ifelse(is.na(p),"", ifelse(p<0.001,"***", ifelse(p<0.01,"**", ifelse(p<0.05,"*",""))))

  dunnett_df <- lapply(split(df, df$Gene), function(sub){
    if (dplyr::n_distinct(sub$Group) < 2) return(NULL)
    fit <- aov(logCPM ~ Group, data=sub)
    emm <- emmeans(fit, ~ Group)
    dtab <- as.data.frame(summary(contrast(emm, method="dunnett", ref="PBS"), infer=c(TRUE,TRUE)))
    dtab$Gene  <- unique(sub$Gene)
    dtab$Group <- gsub(" - PBS","", dtab$contrast, fixed=TRUE)
    transmute(dtab, Gene, Group=factor(Group, levels=c("PBS","Lo","Med","Hi")), pval=p.value, stars=get_stars(p.value))
  }) |> bind_rows()

  pbs_stub <- expand.grid(Gene=unique(df$Gene), Group="PBS", stringsAsFactors=FALSE) |>
    mutate(pval=NA_real_, stars="")
  dunnett_df <- bind_rows(dunnett_df, pbs_stub)

  y_pos <- df |>
    group_by(Gene, Group) |>
    summarise(y = mean(logCPM, na.rm=TRUE) + sem(logCPM) + 0.12, .groups="drop")
  labels_df <- left_join(y_pos, dunnett_df, by=c("Gene","Group"))

  p <- ggplot(df, aes(Group, logCPM, fill=Group)) +
    stat_summary(fun=mean, geom="col", width=.62, alpha=.8, color="black") +
    stat_summary(fun.data=function(.x){ m<-mean(.x,na.rm=TRUE); s<-sem(.x); data.frame(y=m, ymin=m-s, ymax=m+s) },
                 geom="errorbar", width=.22, linewidth=.7) +
    geom_point(position=position_jitter(width=.10, height=0, seed=42),
               size=2.3, shape=21, color="black", alpha=.9) +
    geom_text(data=labels_df, aes(x=Group, y=y, label=stars), inherit.aes=FALSE, vjust=0, size=4) +
    facet_wrap(~ Gene, nrow=1, scales="free_y") +
    labs(title="GOI — mean±SEM (Dunnett vs PBS)", y="log2(CPM+1)", x=NULL) +
    theme_minimal(base_size=12) + theme(legend.position="none", panel.grid.major.x=element_blank())
  ggsave(fname, p, width = max(7, 3 + 2.5*length(unique(df$Gene))), height=4, dpi=300)
}

# -------- run --------
go_focus <- go_focus  # keep in scope

rng <- range(de_tbl$mean_abs_logFC, na.rm = TRUE); if (diff(rng)==0) rng <- c(0,1)
de_ranked <- de_tbl |>
  mutate(in_GO_focus = SYMBOL %in% go_focus,
         score = 2*in_GO_focus + 1*dir_consistent +
                 scales::rescale(mean_abs_logFC, to=c(0,1), from=rng)) |>
  arrange(desc(score), desc(mean_abs_logFC))

short_strict <- de_ranked |> filter(pass_any & in_GO_focus) |> slice_head(n = 150)
short_broad  <- de_ranked |> filter(in_GO_focus | dir_consistent) |> slice_head(n = 300)

readr::write_csv(de_ranked,   file.path(cfg$paths$results, "goi", "GOI_ranked_all.csv"))
readr::write_csv(short_strict,file.path(cfg$paths$results, "goi", "GOI_shortlist_strict.csv"))
readr::write_csv(short_broad, file.path(cfg$paths$results, "goi", "GOI_shortlist_broad.csv"))

auto_goi <- short_strict |> filter(!is.na(SYMBOL) & SYMBOL != "") |> slice_head(n=12) |> pull(SYMBOL) |> unique()
if (length(auto_goi) == 0) {
  auto_goi <- de_ranked |> filter(!is.na(SYMBOL) & SYMBOL != "") |> slice_head(n=12) |> pull(SYMBOL) |> unique()
}

plot_violin_goi(auto_goi,   file.path(cfg$paths$figures, "goi", "GOI_violin_auto12.png"))
plot_heatmap_goi(auto_goi,  file.path(cfg$paths$figures, "goi", "GOI_heatmap_auto12.png"))
plot_barjit_dunnett_goi(auto_goi, file.path(cfg$paths$figures, "goi", "GOI_barjit_dunnett_auto12.png"))
message("GOI CSVs saved under outputs/results/goi/ and figures under outputs/figures/goi/")
