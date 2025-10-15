#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse); library(data.table); library(yaml)
  library(edgeR)
  library(msigdbr); library(fgsea)
  library(ComplexHeatmap); library(circlize)
  library(UpSetR); library(VennDiagram)
  library(RColorBrewer); library(grid)
})

`%||%` <- function(a, b) if (is.null(a)) b else a
safe_dir <- function(...) { d <- file.path(...); dir.create(d, showWarnings = FALSE, recursive = TRUE); d }
say <- function(...) cat(paste0(..., "\n"))
load_cfg <- function(path = "config.yaml") yaml::read_yaml(path)

# ---------- plotting helpers (embedded) ----------
plot_enrich_heatmap <- function(NES, FDR, out_png, title = "Pathway NES across contrasts") {
  if (nrow(NES) == 0 || ncol(NES) == 0) { say("Heatmap skipped: empty NES matrix"); return(invisible(FALSE)) }
  pal <- colorRamp2(c(-2, 0, 2), c("#3b4cc0", "#f7f7f7", "#b40426"))
  sig <- ifelse(FDR < 0.05, "*", "")
  hm <- Heatmap(
    NES, name = "NES", col = pal, cluster_rows = TRUE, cluster_columns = TRUE,
    column_title = title,
    cell_fun = function(j, i, x, y, w, h, fill) { grid.text(sig[i, j], x, y) }
  )
  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
  png(out_png, width = 1600, height = 1400, res = 180); draw(hm); dev.off()
  say("Saved heatmap: ", out_png); invisible(TRUE)
}

plot_enrich_upset <- function(sig_list, out_png) {
  if (length(sig_list) < 2) { say("UpSet skipped: <2 contrasts"); return(invisible(FALSE)) }
  if (all(lengths(sig_list) == 0)) { say("UpSet skipped: no significant pathways in any contrast"); return(invisible(FALSE)) }
  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
  png(out_png, width = 1400, height = 900, res = 150)
  suppressWarnings(
    UpSetR::upset(UpSetR::fromList(sig_list),
                  nsets = length(sig_list), nintersects = 40, order.by = "freq")
  )
  dev.off()
  say("Saved UpSet:   ", out_png); invisible(TRUE)
}

plot_enrich_venn <- function(sig_list, out_png) {
  if (length(sig_list) > 5 || length(sig_list) < 2) { say("Venn skipped: need 2–5 contrasts"); return(invisible(FALSE)) }
  if (all(lengths(sig_list) == 0)) { say("Venn skipped: no significant pathways"); return(invisible(FALSE)) }
  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
  png(out_png, width = 1200, height = 1000, res = 160)
  g <- try(VennDiagram::venn.diagram(sig_list, filename = NULL,
                                     fill = brewer.pal(max(3, length(sig_list)), "Set2")),
           silent = TRUE)
  if(!inherits(g,"try-error")) grid.draw(g)
  dev.off()
  say("Saved Venn:    ", out_png); invisible(TRUE)
}

plot_enrich_violin <- function(cpm_long, genes, out_png, title = "Leading-edge mean CPM") {
  if (length(genes) == 0) { say("Violin skipped: empty gene set"); return(invisible(FALSE)) }
  df <- cpm_long %>% filter(gene %in% genes) %>% group_by(sample, group) %>%
    summarize(meanCPM = mean(value, na.rm = TRUE), .groups = "drop")
  if (nrow(df) == 0) { say("Violin skipped: no CPM rows for gene set"); return(invisible(FALSE)) }
  p <- ggplot(df, aes(group, meanCPM)) +
    geom_violin(trim = FALSE) + geom_jitter(width = 0.15, alpha = 0.6, size = 1.1) +
    labs(title = title, x = NULL, y = "Mean CPM of leading-edge genes") + theme_bw(base_size = 12)
  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
  ggsave(out_png, p, width = 8, height = 5, dpi = 200)
  say("Saved violin:  ", out_png); invisible(TRUE)
}

# ---------- I/O helpers ----------
read_counts <- function(path) {
  dt <- fread(path, data.table = FALSE)
  nms <- names(dt); lname <- tolower(nms)
  if (!"gene" %in% nms) {
    cand_symbol <- which(lname %in% c("symbol","gene_symbol","genesymbol","hgnc_symbol","gene.name","genename"))
    cand_id     <- which(lname %in% c("gene","geneid","gene_id","ensembl","ensembl_id","id"))
    if (length(cand_symbol))      dt$gene <- dt[[ nms[cand_symbol[1]] ]]
    else if (length(cand_id))     dt$gene <- dt[[ nms[cand_id[1]] ]]
  }
  if (!"gene" %in% names(dt)) stop("Could not find a gene column (symbol/gene/id-like) in: ", path)
  dt
}

infer_group <- function(nms) {
  case_when(
    str_detect(nms, regex("^PBS[_-]", ignore_case=TRUE)) ~ "PBS",
    str_detect(nms, regex("^Lo[_-]",  ignore_case=TRUE)) ~ "Lo",
    str_detect(nms, regex("^Med[_-]", ignore_case=TRUE)) ~ "Med",
    str_detect(nms, regex("^Hi[_-]",  ignore_case=TRUE)) ~ "Hi",
    TRUE ~ "UNK"
  )
}

# Mouse-native sets: Hallmark = MH; GO:BP = M5 + subcollection "GO:BP"
get_genesets <- function(species="Mus musculus", kind=c("H","GO")) {
  kind <- match.arg(kind)
  if (kind=="H") {
    gs <- msigdbr(species=species, db_species="MM", collection="MH") %>% dplyr::select(gs_name, gene_symbol)
    label <- "Hallmark"
  } else {
    gs <- msigdbr(species=species, db_species="MM", collection="M5", subcollection="GO:BP") %>% dplyr::select(gs_name, gene_symbol)
    label <- "GOBP"
  }
  list(pathways = split(gs$gene_symbol, gs$gs_name), label=label)
}

read_de <- function(dir="outputs/results", contrast) {
  f <- file.path(dir, sprintf("DE_%s.tsv", contrast))
  if (!file.exists(f)) stop("Missing DE results: ", f)
  df <- fread(f, data.table=FALSE, sep="\t")
  lower <- tolower(names(df))
  find1 <- function(cands) { ix <- which(lower %in% tolower(cands)); if (length(ix)) names(df)[ix[1]] else NULL }

  gcol <- find1(c("gene","Gene","symbol","SYMBOL","gene_symbol","genesymbol","hgnc_symbol","gene.name","genename","id"))
  if (!is.null(gcol) && !"gene" %in% names(df)) df$gene <- df[[gcol]]

  pcol   <- find1(c("P.Value","p.value","pvalue","pval","pr(>|t|)","p.value...wald","pval_wald"))
  adjcol <- find1(c("adj.P.Val","adj.p.val","fdr","padj","qval","adjp","adj.p.value"))
  tcol   <- find1(c("t","T","stat","statistic","wald.stat","score"))

  if (!"gene" %in% names(df)) stop("No gene column in ", f, " (looked for Gene/gene/SYMBOL/symbol/id).")
  if (!"logFC" %in% names(df)) {
    lcol <- find1(c("logFC","logfc","log2fc","log2foldchange","log2_fold_change"))
    if (!is.null(lcol)) names(df)[match(lcol, names(df))] <- "logFC"
  }
  if (!is.null(pcol)   && !"P.Value"   %in% names(df)) names(df)[match(pcol,   names(df))] <- "P.Value"
  if (!is.null(adjcol) && !"adj.P.Val" %in% names(df)) names(df)[match(adjcol, names(df))] <- "adj.P.Val"
  if (!is.null(tcol)   && !"t"         %in% names(df)) names(df)[match(tcol,   names(df))] <- "t"

  if ("P.Value" %in% names(df)) {
    df$rank_metric <- with(df, sign(logFC) * -log10(pmax(P.Value, 1e-300)))
  } else if ("adj.P.Val" %in% names(df)) {
    df$rank_metric <- with(df, sign(logFC) * -log10(pmax(adj.P.Val, 1e-300)))
  } else if ("t" %in% names(df)) {
    df$rank_metric <- with(df, sign(logFC) * abs(t))
  } else {
    warning("No P.Value/adj.P.Val/statistic found in ", f, " — ranking by |logFC|.")
    df$rank_metric <- with(df, sign(logFC) * abs(logFC))
  }
  df %>% filter(!is.na(gene) & gene != "")
}

make_ranks <- function(de_df) {
  de_df %>% group_by(gene) %>% slice_max(order_by=abs(rank_metric), n=1, with_ties=FALSE) %>%
    ungroup() %>% arrange(desc(rank_metric)) %>% { setNames(.$rank_metric, .$gene) }
}

# ----------------------------- MAIN ---------------------------------------
cfg <- load_cfg()
counts_path <- cfg$paths$counts %||% "data/raw_annotated_combined.counts"
de_dir      <- cfg$paths$results %||% "outputs/results"
fig_root    <- cfg$paths$figures %||% "outputs/figures"
fig_dir     <- safe_dir(fig_root, "enrich")
out_tab_dir <- safe_dir("outputs","enrichment")

contrasts   <- cfg$params$contrasts
fdr_thr     <- cfg$params$fdr_thresh %||% 0.05
n_show      <- cfg$params$go_n_show %||% 30
species     <- cfg$params$species %||% "Mus musculus"
gs_kind     <- cfg$params$gene_set_kind %||% "H"

say("Loading raw counts from: ", counts_path)
raw <- read_counts(counts_path)

sample_cols <- setdiff(
  names(raw),
  c("gene","Gene","GeneID","ENTREZID","SYMBOL","symbol","id","ID","Ensembl","ensembl","EnsemblID","ensembl_id")
)
counts_mat <- as.matrix(raw[, sample_cols, drop=FALSE]); mode(counts_mat) <- "numeric"
rownames(counts_mat) <- raw$gene

say("Computing CPMs for violin plots...")
cpm_mat <- edgeR::cpm(counts_mat, log=FALSE)
cpm_long <- as.data.frame(cpm_mat) %>% rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to="sample", values_to="value") %>%
  mutate(group = infer_group(sample))

say("Loading gene sets (mouse-native, ", ifelse(gs_kind=='H','MH','M5/GO:BP'), ") ...")
gs <- get_genesets(species=species, kind=gs_kind)
pathways <- gs$pathways; gs_label <- gs$label

say("Running fgsea for contrasts: ", paste(contrasts, collapse=", "))
all_res <- list(); sig_sets <- list()
for (ct in contrasts) {
  say("  - ", ct)
  de <- read_de(dir=de_dir, contrast=ct)
  ranks <- make_ranks(de)

  # Prefer fgseaMultilevel with eps=0 to better estimate tiny p-values; fallback to fgsea(nperm=10000)
  fg <- try({
    fgseaMultilevel(pathways=pathways, stats=ranks, minSize=10, maxSize=500, eps=0)
  }, silent=TRUE)
  if (inherits(fg, "try-error")) {
    fg <- fgsea(pathways=pathways, stats=ranks, minSize=10, maxSize=500, nperm=10000)
  }

  fg <- fg %>% arrange(padj, desc(NES)) %>% mutate(contrast=ct)
  fwrite(fg, file.path(out_tab_dir, sprintf("fgsea_%s_%s.csv", gs_label, ct)))
  saveRDS(fg, file.path(out_tab_dir, sprintf("fgsea_%s_%s.rds", gs_label, ct)))
  all_res[[ct]] <- fg
  sig_sets[[ct]] <- fg %>% filter(padj < fdr_thr) %>% pull(pathway)
}

say("Building summary and visualizations...")
comb <- bind_rows(all_res); comb$padj[is.na(comb$padj)] <- 1
minpadj <- comb %>% group_by(pathway) %>% summarize(min_padj=min(padj, na.rm=TRUE), .groups="drop")
top_paths <- minpadj %>% arrange(min_padj) %>% slice_head(n=max(n_show,30)) %>% pull(pathway)

NES <- matrix(0, nrow=length(top_paths), ncol=length(contrasts), dimnames=list(top_paths, contrasts))
FDR <- matrix(1, nrow=length(top_paths), ncol=length(contrasts), dimnames=list(top_paths, contrasts))
for (ct in contrasts) {
  fg <- all_res[[ct]]; idx <- match(top_paths, fg$pathway); keep <- !is.na(idx)
  NES[keep, ct] <- fg$NES[idx[keep]]; FDR[keep, ct] <- fg$padj[idx[keep]]
}

fwrite(comb, file.path(out_tab_dir, "summary_enrichment_long.csv"))
say("Saved table:   ", file.path(out_tab_dir, "summary_enrichment_long.csv"))

# ---- FIGURES ----
plot_enrich_heatmap(NES, FDR, file.path(fig_dir, sprintf("Heatmap_NES_%s.png", gs_label)),
                    title = paste0(gs_label, " pathways (NES)"))
plot_enrich_upset(sig_sets, file.path(fig_dir, sprintf("UpSet_%s.png", gs_label)))
plot_enrich_venn(sig_sets,  file.path(fig_dir, sprintf("Venn_%s.png",  gs_label)))

# Leading-edge violins
shared_counts <- sort(table(unlist(sig_sets)), decreasing=TRUE)
if (length(shared_counts)) {
  shared_term <- names(shared_counts)[1]
  owner <- names(Filter(function(x) shared_term %in% x, sig_sets))[1]
  le <- all_res[[owner]] %>% filter(pathway==shared_term) %>% pull(leadingEdge) %>% .[[1]]
  if (length(le)>0) {
    le_csv <- file.path(out_tab_dir, sprintf("LeadingEdge_%s.csv", make.names(shared_term)))
    fwrite(data.frame(pathway=shared_term, gene=le), le_csv)
    say("Saved leading-edge genes: ", le_csv)
    plot_enrich_violin(cpm_long, le,
                       file.path(fig_dir, sprintf("Violin_%s_shared.png", make.names(shared_term))),
                       title = sprintf("Leading-edge mean CPM: %s (shared)", shared_term))
  }
}
for (ct in contrasts) {
  uniq <- setdiff(sig_sets[[ct]], unique(unlist(sig_sets[names(sig_sets)!=ct])))
  pick <- if(length(uniq)>0) uniq[1] else (all_res[[ct]] %>% arrange(padj) %>% slice(1) %>% pull(pathway))
  fe <- all_res[[ct]] %>% filter(pathway==pick) %>% pull(leadingEdge) %>% .[[1]]
  if (length(fe)>0) {
    le_csv <- file.path(out_tab_dir, sprintf("LeadingEdge_%s_%s.csv", make.names(pick), ct))
    fwrite(data.frame(pathway=pick, gene=fe), le_csv)
    say("Saved leading-edge genes: ", le_csv)
    plot_enrich_violin(cpm_long, fe,
                       file.path(fig_dir, sprintf("Violin_%s_%s.png", make.names(pick), ct)),
                       title = sprintf("Leading-edge mean CPM: %s (%s)", pick, ct))
  }
}

say("Tables dir: ", out_tab_dir)
say("Figures dir:", fig_dir)
