# 07_goi_zoom.R — Deep dive for selected GOI

source("R/io_helpers.R")
suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
  library(edgeR)
  library(ggrepel)
  library(ragg)
  library(AnnotationDbi); library(org.Mm.eg.db)
  library(grid)
  # optional
  suppressWarnings({
    if (requireNamespace("ggVennDiagram", quietly = TRUE)) library(ggVennDiagram)
    if (requireNamespace("VennDiagram",  quietly = TRUE)) library(VennDiagram)
    if (requireNamespace("GGally",        quietly = TRUE)) library(GGally)
  })
  # needed to coerce enrichResult objects saved in go_objs.rds
  if (requireNamespace("DOSE", quietly = TRUE)) library(DOSE)
})

`%||%` <- function(a,b) if (!is.null(a)) a else b

cfg   <- load_cfg()
obj1  <- readRDS(file.path(cfg$paths$rds, "norm_and_fit.rds"))
obj2  <- readRDS(file.path(cfg$paths$rds, "de_tables.rds"))
goes  <- file.path(cfg$paths$rds, "go_objs.rds")
obj3  <- if (file.exists(goes)) readRDS(goes) else NULL

annot <- obj1$annot; meta <- obj1$meta; dge <- obj1$dge; v <- obj1$v; fit2 <- obj1$fit2

# ---------- parameters ----------
goi     <- c("Bub1","Cenpi","Esco2","Rpl36-ps12")
fdr_thr <- cfg$params$fdr_thresh %||% 0.05
lfc_thr <- cfg$params$lfc_thresh %||% 1
outF    <- file.path(cfg$paths$figures, "goi_zoom")
outR    <- file.path(cfg$paths$results, "goi_zoom")
dir.create(outF, showWarnings = FALSE, recursive = TRUE)
dir.create(outR, showWarnings = FALSE, recursive = TRUE)

# ---------- 1) Robust GOI matching in voom E ----------
rnames <- rownames(v$E)
match_sym <- match(goi, rnames)
if (any(is.na(match_sym))) {
  for (i in which(is.na(match_sym))) {
    cand <- which(tolower(rnames) == tolower(goi[i]))
    match_sym[i] <- if (length(cand)) cand[1] else NA_integer_
  }
}
present <- !is.na(match_sym)
stopifnot(any(present))
goi_found <- goi[present]

mat <- v$E[goi_found, , drop = FALSE]
grp <- factor(meta$Group[match(colnames(mat), meta$SampleID)],
              levels = c("PBS","Lo","Med","Hi"))

df_long <- as.data.frame(mat) |>
  tibble::rownames_to_column("Gene") |>
  tidyr::pivot_longer(-Gene, names_to = "Sample", values_to = "logCPM") |>
  dplyr::mutate(Group = grp[match(Sample, colnames(mat))]) |>
  dplyr::filter(!is.na(Group)) |>
  dplyr::mutate(Gene = factor(Gene, levels = goi_found),
                Group = factor(Group, levels = c("PBS","Lo","Med","Hi")))

# ---------- 2) limma topTables for Lo/Med/Hi (dynamic coef names) ----------
norm_name  <- function(x) gsub(" - ", "_vs_", x, fixed = TRUE)
avail_raw  <- colnames(fit2$coefficients)
avail_norm <- norm_name(avail_raw)
get_tt <- function(name_norm) {
  idx <- match(name_norm, avail_norm)
  if (is.na(idx)) return(NULL)
  topTable(fit2, coef = avail_raw[idx], number = Inf, sort.by = "P", confint = TRUE)
}
tt_list <- list(
  Lo  = get_tt("Lo_vs_PBS"),
  Med = get_tt("Med_vs_PBS"),
  Hi  = get_tt("Hi_vs_PBS")
)

# ---------- 3) Lollipop: GOI logFC (±CI) ----------
fc_tab <- bind_rows(lapply(names(tt_list), function(k){
  tt <- tt_list[[k]]; if (is.null(tt)) return(NULL)
  sub <- tt[rownames(tt) %in% goi_found, , drop = FALSE]
  if (!nrow(sub)) return(NULL)
  tibble(Contrast = k, Gene = rownames(sub),
         logFC = sub$logFC,
         CI.L  = sub$CI.L %||% NA_real_,
         CI.R  = sub$CI.R %||% NA_real_)
}), .id = NULL)

if (nrow(fc_tab)) {
  fc_tab$Contrast <- factor(fc_tab$Contrast, levels = c("Lo","Med","Hi"))
  p_lolli <- ggplot(fc_tab, aes(x = Gene, y = logFC, color = Contrast)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_pointrange(aes(ymin = CI.L, ymax = CI.R),
                    position = position_dodge(width = 0.5)) +
    coord_flip() +
    labs(title = "GOI logFC by contrast (limma)",
         y = "log2 fold-change (vs PBS)", x = "") +
    theme_minimal()
  ggsave(file.path(outF, "GOI_lollipop_logFC.png"), p_lolli, width = 7.5, height = 4.5, dpi = 300)
  readr::write_csv(fc_tab, file.path(outR, "GOI_logFC_table.csv"))
}

# ---------- 3a) Significant UP sets & lists (Lo/Med/Hi) ----------
sig_up <- lapply(names(tt_list), function(k){
  tt <- tt_list[[k]]
  if (is.null(tt)) return(character(0))
  rownames(tt)[which(tt$adj.P.Val < fdr_thr & tt$logFC > lfc_thr)]
})
names(sig_up) <- names(tt_list)
for (k in names(sig_up)) {
  readr::write_lines(sig_up[[k]], file.path(outR, paste0("UP_", k, "_sig_genes.txt")))
}

# helper: draw Venn safely (2–3 sets)
draw_venn <- function(sets_named, title_txt, outfile){
  if (length(sets_named) < 2) return(invisible(NULL))
  # ensure named list via setNames (avoid data.table's :=)
  sets_named <- setNames(lapply(sets_named, unique), names(sets_named))
  if (requireNamespace("ggVennDiagram", quietly = TRUE)) {
    p <- ggVennDiagram::ggVennDiagram(sets_named) + ggplot2::labs(title = title_txt)
    ggplot2::ggsave(outfile, p, width = 6, height = 5, dpi = 300)
  } else if (requireNamespace("VennDiagram", quietly = TRUE)) {
    png(outfile, width = 1200, height = 900, res = 150)
    grid::grid.newpage()
    VennDiagram::grid.draw(
      VennDiagram::venn.diagram(x = sets_named, fill = c("#66c2a5","#8da0cb","#fc8d62")[seq_along(sets_named)],
                                alpha = .5, cex = 1.4, cat.cex = 1.3, lwd = 2, filename = NULL)
    )
    dev.off()
  }
}

# draw Lo/Med/Hi UP venn if at least two non-empty
non_empty <- names(sig_up)[vapply(sig_up, length, 1L) > 0]
if (length(non_empty) >= 2) {
  draw_venn(sig_up[non_empty], "Significant UP genes (FDR<0.05 & |LFC|>1)",
            file.path(outF, "Venn_UP_Lo_Med_Hi.png"))
}

# ---------- 3b) UpSet-like bar plot (no extra deps) ----------
L <- unique(sig_up$Lo  %||% character(0))
M <- unique(sig_up$Med %||% character(0))
H <- unique(sig_up$Hi  %||% character(0))
lab_counts <- tibble::tibble(
  pattern = c("Lo","Med","Hi","Lo∩Med","Lo∩Hi","Med∩Hi","Lo∩Med∩Hi"),
  n = c(
    length(setdiff(L, union(M, H))),
    length(setdiff(M, union(L, H))),
    length(setdiff(H, union(L, M))),
    length(intersect(L, setdiff(M, H))),
    length(intersect(L, setdiff(H, M))),
    length(intersect(M, setdiff(H, L))),
    length(Reduce(intersect, list(L, M, H)))
  )
)
if (sum(lab_counts$n > 0) == 0) lab_counts$n[lab_counts$pattern=="Lo"] <- length(L)
p_upset <- ggplot(lab_counts %>% dplyr::filter(n > 0), aes(x = reorder(pattern, -n), y = n)) +
  geom_col(fill = "steelblue") + geom_text(aes(label = n), vjust = -0.2, size = 3.5) +
  labs(title = "Significant UP genes — overlaps (FDR<0.05 & |LFC|>1)",
       x = "Pattern", y = "Count") + theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(outF, "UP_overlaps_bar.png"), p_upset, width = 7, height = 4.2, dpi = 300)

# ---------- 4) Correlation heatmap (gene x gene) ----------
cors <- cor(t(mat), use = "pairwise.complete.obs")
cors_df <- as.data.frame(cors) |>
  tibble::rownames_to_column("Gene1") |>
  tidyr::pivot_longer(-Gene1, names_to = "Gene2", values_to = "r")
p_cor <- ggplot(cors_df, aes(Gene1, Gene2, fill = r)) +
  geom_tile() +
  scale_fill_gradient2(limits = c(-1,1)) +
  coord_equal() +
  labs(title = "GOI correlation (voom logCPM across samples)", x = "", y = "", fill = "r") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(outF, "GOI_correlation_heatmap.png"), p_cor, width = 5.5, height = 5, dpi = 300)

# ---------- 5) Pairwise scatter (if GGally available) ----------
if (requireNamespace("GGally", quietly = TRUE)) {
  df_pairs <- df_long |>
    dplyr::select(Gene, Sample, logCPM) |>
    tidyr::pivot_wider(names_from = Gene, values_from = logCPM) |>
    dplyr::mutate(Group = grp[match(Sample, colnames(mat))]) |>
    dplyr::filter(!is.na(Group))
  p_pairs <- GGally::ggpairs(df_pairs, columns = colnames(df_pairs)[colnames(df_pairs) %in% goi_found],
                             aes(color = Group, alpha = 0.85))
  ggsave(file.path(outF, "GOI_pairwise_scatter.png"), p_pairs, width = 8.5, height = 8, dpi = 250)
}

# ---------- 6) Leading-edge membership (if gse objects exist) ----------
make_le_table <- function(g_obj, gene_symbols, label){
  if (is.null(g_obj)) return(NULL)
  df <- as.data.frame(g_obj); if (!nrow(df)) return(NULL)
  ids <- AnnotationDbi::select(org.Mm.eg.db, keys = gene_symbols, keytype = "SYMBOL", columns = "ENTREZID")
  ids <- ids[!is.na(ids$ENTREZID),]
  if (!nrow(ids)) return(NULL)
  out <- list()
  for (i in seq_len(nrow(ids))) {
    sym <- ids$SYMBOL[i]; ent <- ids$ENTREZID[i]
    hit <- df %>% dplyr::filter(grepl(paste0("\\b", ent, "\\b"), core_enrichment), p.adjust < 0.05)
    if (nrow(hit)) {
      hit$Gene <- sym
      out[[sym]] <- hit[, c("Gene","ID","Description","NES","p.adjust","core_enrichment")]
    }
  }
  if (!length(out)) return(NULL)
  res <- dplyr::bind_rows(out); res$Contrast <- label; res
}
le_tab <- NULL
if (!is.null(obj3)) {
  le_lo  <- make_le_table(obj3$g_lo,  goi_found, "Lo_vs_PBS")
  le_med <- make_le_table(obj3$g_med, goi_found, "Med_vs_PBS")
  le_hi  <- make_le_table(obj3$g_hi,  goi_found, "Hi_vs_PBS")
  le_tab <- dplyr::bind_rows(le_lo, le_med, le_hi)
  if (!is.null(le_tab) && nrow(le_tab)) {
    readr::write_csv(le_tab, file.path(outR, "GOI_leading_edge_terms.csv"))
  }
}

# ---------- 7) Focused bar + jitter (mean±SEM) ----------
sem <- function(x){ x <- x[is.finite(x)]; if (!length(x)) NA_real_ else sd(x)/sqrt(length(x)) }
p_barjit <- ggplot(df_long, aes(x = Group, y = logCPM, fill = Group)) +
  stat_summary(fun = mean, geom = "col", width = 0.62, alpha = 0.75, color = "black") +
  stat_summary(fun.data = function(.x){
    m <- mean(.x, na.rm=TRUE); s <- sem(.x)
    data.frame(y=m, ymin=m-s, ymax=m+s)
  }, geom="errorbar", width=0.22, linewidth=0.7) +
  geom_point(position = position_jitter(width = 0.10, height = 0, seed=42),
             size = 2.2, shape = 21, color = "black", alpha = 0.9) +
  facet_wrap(~ Gene, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c(PBS="#9e9e9e", Lo="#2ca02c", Med="#1f77b4", Hi="#d62728")) +
  labs(title = "GOI (voom log2-CPM) — mean±SEM with per-sample jitter",
       y = "log2(CPM+1)", x = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"),
        panel.grid.major.x = element_blank())
ggsave(file.path(outF, "GOI_bar_jitter_focus.png"), p_barjit, width = 12, height = 3.8, dpi = 300)

# ---------- 8) LO-only: GOI vs Lo-UP/DOWN venn + status table ----------
tt_lo_only <- tt_list$Lo
if (!is.null(tt_lo_only)) {
  lo_up   <- rownames(tt_lo_only)[which(tt_lo_only$adj.P.Val < fdr_thr & tt_lo_only$logFC >  lfc_thr)]
  lo_down <- rownames(tt_lo_only)[which(tt_lo_only$adj.P.Val < fdr_thr & tt_lo_only$logFC < -lfc_thr)]

  stat_lo <- tt_lo_only[rownames(tt_lo_only) %in% goi_found, , drop = FALSE]
  status_tbl <- tibble::tibble(
    Gene      = goi_found,
    logFC     = stat_lo[goi_found, "logFC"] %||% NA_real_,
    adj.P.Val = stat_lo[goi_found, "adj.P.Val"] %||% NA_real_,
    sig_up    = goi_found %in% lo_up,
    sig_down  = goi_found %in% lo_down
  )
  readr::write_csv(status_tbl, file.path(outR, "GOI_status_Lo.csv"))

  draw_2set_venn <- function(setA, setB, nameA, nameB, outfile){
    nl <- setNames(list(unique(setA), unique(setB)), c(nameA, nameB))
    if (requireNamespace("ggVennDiagram", quietly = TRUE)) {
      p <- ggVennDiagram::ggVennDiagram(nl) + ggplot2::labs(title = paste0(nameA, " ∩ ", nameB))
      ggplot2::ggsave(outfile, p, width = 5.5, height = 4.5, dpi = 300)
    } else if (requireNamespace("VennDiagram", quietly = TRUE)) {
      png(outfile, width = 1100, height = 900, res = 150)
      grid::grid.newpage()
      VennDiagram::grid.draw(
        VennDiagram::venn.diagram(x = nl, fill = c("#66c2a5","#fc8d62"),
                                  alpha = .5, cex = 1.5, cat.cex = 1.5,
                                  lwd = 2, filename = NULL)
      )
      dev.off()
    }
  }

  draw_2set_venn(goi_found, lo_up,  "GOI", "Lo-UP (FDR<0.05 & |LFC|>1)",
                  file.path(outF, "Venn_GOI_vs_Lo_UP.png"))
  draw_2set_venn(goi_found, lo_down,"GOI", "Lo-DOWN (FDR<0.05 & |LFC|>1)",
                  file.path(outF, "Venn_GOI_vs_Lo_DOWN.png"))
}

message("GOI deep-dive outputs -> ", outF, " ; tables -> ", outR)
