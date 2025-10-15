#!/usr/bin/env Rscript

## Multi-contrast pathway enrichment (RAW dataset; TSV DE files in outputs/results/)
## - Mouse-native MSigDB (db_species="MM"): Hallmark = MH, GO:BP = M5/GO:BP
## - Reads DE_<Contrast>.tsv from outputs/results/ (tolerant Gene/gene/SYMBOL/etc.)
## - Outputs heatmap (NES), UpSet/Venn, and violins for leading-edge genes

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
load_cfg <- function(path = "config.yaml") yaml::read_yaml(path)

# ---------- I/O helpers ----------
read_counts <- function(path) {
  dt <- fread(path, data.table = FALSE)
  nms <- names(dt); lname <- tolower(nms)
  if (!"gene" %in% nms) {
    cand_symbol <- which(lname %in% c("symbol","gene_symbol","genesymbol","hgnc_symbol","gene.name","genename"))
    cand_id     <- which(lname %in% c("gene","geneid","gene_id","ensembl","ensembl_id","id"))
    if (length(cand_symbol)) {
      dt$gene <- dt[[ nms[cand_symbol[1]] ]]
    } else if (length(cand_id)) {
      dt$gene <- dt[[ nms[cand_id[1]] ]]
    }
  }
  if (!"gene" %in% names(dt)) stop("Could not find a gene column (symbol/gene/id-like) in: ", path)
  dt
}

infer_group <- function(nms) {
  case_when(
    str_detect(nms, regex("\\bPBS\\b", ignore_case=TRUE)) ~ "PBS",
    str_detect(nms, regex("\\bLo\\b",  ignore_case=TRUE)) ~ "Lo",
    str_detect(nms, regex("\\bMed\\b", ignore_case=TRUE)) ~ "Med",
    str_detect(nms, regex("\\bHi\\b",  ignore_case=TRUE)) ~ "Hi",
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

  pcol <- find1(c("P.Value","p.value","pvalue","pval","pr(>|t|)","p.value...wald","pval_wald"))
  if (!is.null(pcol) && !"P.Value" %in% names(df)) names(df)[match(pcol, names(df))] <- "P.Value"

  fcol <- find1(c("adj.P.Val","adj.p.val","fdr","padj","qval","adjp","adj.p.value"))
  if (!is.null(fcol) && !"adj.P.Val" %in% names(df)) names(df)[match(fcol, names(df))] <- "adj.P.Val"

  if (!"gene" %in% names(df)) stop("No gene column in ", f, " (looked for Gene/gene/SYMBOL/symbol/id).")
  if (!"logFC" %in% names(df)) {
    lcol <- find1(c("logFC","logfc","log2fc","log2foldchange","log2_fold_change"))
    if (!is.null(lcol)) names(df)[match(lcol, names(df))] <- "logFC"
  }
  if (!"P.Value" %in% names(df)) stop("No P.Value-like column in ", f, ".")

  df$rank_metric <- with(df, sign(logFC) * -log10(pmax(P.Value, 1e-300)))
  df %>% filter(!is.na(gene) & gene != "")
}

make_ranks <- function(de_df) {
  de_df %>%
    group_by(gene) %>%
    slice_max(order_by=abs(rank_metric), n=1, with_ties=FALSE) %>%
    ungroup() %>%
    arrange(desc(rank_metric)) %>%
    { setNames(.$rank_metric, .$gene) }
}

plot_heatmap <- function(mat, fdr_mat, out_png) {
  pal <- colorRamp2(c(-2,0,2), c("#3b4cc0","#f7f7f7","#b40426"))
  sig <- ifelse(fdr_mat < 0.05, "*", "")
  hm <- Heatmap(mat, name="NES", col=pal, cluster_rows=TRUE, cluster_columns=TRUE,
                cell_fun=function(j,i,x,y,w,h,fill){ grid.text(sig[i,j], x, y) })
  png(out_png, width=1600, height=1400, res=180); draw(hm); dev.off()
}

plot_upset <- function(sig_list, out_png) {
  if (length(sig_list) < 2) return(invisible())
  png(out_png, width=1400, height=900, res=150)
  UpSetR::upset(UpSetR::fromList(sig_list), nsets=length(sig_list), nintersects=40, order.by="freq")
  dev.off()
}

plot_venn <- function(sig_list, out_png) {
  if (length(sig_list) > 5 || length(sig_list) < 2) return(invisible())
  png(out_png, width=1200, height=1000, res=160)
  g <- try(VennDiagram::venn.diagram(sig_list, filename=NULL, fill=brewer.pal(max(3,length(sig_list)),"Set2")), silent=TRUE)
  if(!inherits(g,"try-error")) grid.draw(g)
  dev.off()
}

violin_le <- function(cpm_long, le_genes, title, out_png) {
  df <- cpm_long %>% filter(gene %in% le_genes) %>%
    group_by(sample, group) %>% summarize(meanCPM = mean(value, na.rm=TRUE), .groups="drop")
  p <- ggplot(df, aes(group, meanCPM)) +
    geom_violin(trim=FALSE) + geom_jitter(width=0.15, alpha=0.6, size=1.1) +
    labs(title=title, x=NULL, y="Mean CPM of leading-edge genes") + theme_bw(base_size=12)
  ggsave(out_png, p, width=8, height=5, dpi=200)
}

# ---------- MAIN ----------
cfg <- load_cfg()
counts_path <- cfg$paths$counts %||% "data/raw_annotated_combined.counts"
de_dir      <- cfg$paths$results %||% "outputs/results"
contrasts   <- cfg$params$contrasts
fdr_thr     <- cfg$params$fdr_thresh %||% 0.05
n_show      <- cfg$params$go_n_show %||% 30
species     <- cfg$params$species %||% "Mus musculus"
gs_kind     <- cfg$params$gene_set_kind %||% "H"

outdir <- safe_dir("outputs","enrichment")

message("Loading raw counts from: ", counts_path)
raw <- read_counts(counts_path)

sample_cols <- setdiff(
  names(raw),
  c("gene","Gene","GeneID","ENTREZID","SYMBOL","symbol","id","ID","Ensembl","ensembl","EnsemblID","ensembl_id")
)
counts_mat <- as.matrix(raw[, sample_cols, drop=FALSE]); mode(counts_mat) <- "numeric"
rownames(counts_mat) <- raw$gene

message("Computing CPMs for violin plots...")
cpm_mat <- edgeR::cpm(counts_mat, log=FALSE)
cpm_long <- as.data.frame(cpm_mat) %>% rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to="sample", values_to="value") %>%
  mutate(group = infer_group(sample))

message("Loading gene sets (mouse-native, ", ifelse(gs_kind=='H','MH','M5/GO:BP'), ") ...")
gs <- get_genesets(species=species, kind=gs_kind)
pathways <- gs$pathways; gs_label <- gs$label

message("Running fgsea for contrasts: ", paste(contrasts, collapse=", "))
all_res <- list(); sig_sets <- list()
for (ct in contrasts) {
  message("  - ", ct)
  de <- read_de(dir=de_dir, contrast=ct)
  ranks <- make_ranks(de)
  # >>> Robust to fgsea versions: use fgsea() with nperm <<<
  fg <- fgsea(pathways=pathways, stats=ranks, minSize=10, maxSize=500, nperm=10000) %>%
        arrange(padj, desc(NES)) %>% mutate(contrast=ct)
  fwrite(fg, file.path(outdir, sprintf("fgsea_%s_%s.csv", gs_label, ct)))
  saveRDS(fg, file.path(outdir, sprintf("fgsea_%s_%s.rds", gs_label, ct)))
  all_res[[ct]] <- fg
  sig_sets[[ct]] <- fg %>% filter(padj < fdr_thr) %>% pull(pathway)
}

message("Building summary and visualizations...")
comb <- bind_rows(all_res); comb$padj[is.na(comb$padj)] <- 1
minpadj <- comb %>% group_by(pathway) %>% summarize(min_padj=min(padj, na.rm=TRUE), .groups="drop")
top_paths <- minpadj %>% arrange(min_padj) %>% slice_head(n=max(n_show,30)) %>% pull(pathway)

NES <- matrix(0, nrow=length(top_paths), ncol=length(contrasts), dimnames=list(top_paths, contrasts))
FDR <- matrix(1, nrow=length(top_paths), ncol=length(contrasts), dimnames=list(top_paths, contrasts))
for (ct in contrasts) {
  fg <- all_res[[ct]]; idx <- match(top_paths, fg$pathway); keep <- !is.na(idx)
  NES[keep, ct] <- fg$NES[idx[keep]]; FDR[keep, ct] <- fg$padj[idx[keep]]
}

fwrite(comb, file.path(outdir, "summary_enrichment_long.csv"))
plot_heatmap(NES, FDR, file.path(outdir, sprintf("Heatmap_NES_%s.png", gs_label)))
plot_upset(sig_sets, file.path(outdir, sprintf("UpSet_%s.png", gs_label)))
plot_venn(sig_sets,  file.path(outdir, sprintf("Venn_%s.png", gs_label)))

# Leading-edge violins
shared_counts <- sort(table(unlist(sig_sets)), decreasing=TRUE)
if (length(shared_counts)) {
  shared_term <- names(shared_counts)[1]
  owner <- names(Filter(function(x) shared_term %in% x, sig_sets))[1]
  le <- all_res[[owner]] %>% filter(pathway==shared_term) %>% pull(leadingEdge) %>% .[[1]]
  if (length(le)>0) {
    fwrite(data.frame(pathway=shared_term, gene=le),
           file.path(outdir, sprintf("LeadingEdge_%s.csv", make.names(shared_term))))
    violin_le(cpm_long, le, sprintf("Leading-edge mean CPM: %s (shared)", shared_term),
              file.path(outdir, sprintf("Violin_%s_shared.png", make.names(shared_term))))
  }
}
for (ct in contrasts) {
  uniq <- setdiff(sig_sets[[ct]], unique(unlist(sig_sets[names(sig_sets)!=ct])))
  pick <- if(length(uniq)>0) uniq[1] else (all_res[[ct]] %>% arrange(padj) %>% slice(1) %>% pull(pathway))
  fe <- all_res[[ct]] %>% filter(pathway==pick) %>% pull(leadingEdge) %>% .[[1]]
  if (length(fe)>0) {
    fwrite(data.frame(pathway=pick, gene=fe),
           file.path(outdir, sprintf("LeadingEdge_%s_%s.csv", make.names(pick), ct)))
    violin_le(cpm_long, fe, sprintf("Leading-edge mean CPM: %s (%s)", pick, ct),
              file.path(outdir, sprintf("Violin_%s_%s.png", make.names(pick), ct)))
  }
}

message("RAW-based enrichment complete. Outputs in: ", outdir)
