suppressPackageStartupMessages({
  library(edgeR); library(limma); library(tidyverse)
})

# ---- Normalization: CPM filter + TMM ----
# Matches your monolithic script:
#   - keep genes with CPM >= cpm_min in at least cpm_min_samples samples
#   - TMM normalization
do_normalize <- function(counts, meta, cpm_min = 1, cpm_min_samples = 3) {
  stopifnot(is.data.frame(meta), "Group" %in% names(meta))
  dge <- DGEList(counts = as.matrix(counts), samples = as.data.frame(meta))
  keep <- rowSums(cpm(dge) >= cpm_min) >= cpm_min_samples
  dge  <- dge[keep,, keep.lib.sizes = FALSE]
  dge  <- calcNormFactors(dge, method = "TMM")
  dge
}

# ---- limma-voom fit + contrasts (Lo/Med/Hi vs PBS) ----
# Design: ~ 0 + Group, columns ordered by factor levels of meta$Group
fit_limma <- function(dge, meta) {
  design <- model.matrix(~ 0 + Group, data = meta)
  colnames(design) <- levels(meta$Group)

  v <- voom(dge, design, plot = FALSE)

  contr <- makeContrasts(
    Lo_vs_PBS  = Lo  - PBS,
    Med_vs_PBS = Med - PBS,
    Hi_vs_PBS  = Hi  - PBS,
    levels = design
  )

  fit  <- lmFit(v, design)
  fit2 <- eBayes(contrasts.fit(fit, contr))
  list(v = v, fit2 = fit2)
}

# ---- Snapshot writers (CSV of CPM/voom, plus kept raw counts) ----
write_post_norm <- function(dge, v, annot_keep, cfg) {
  rk <- rownames(dge$counts)

  name_key <- make.unique(ifelse(is.na(annot_keep$SYMBOL) | annot_keep$SYMBOL=="",
                                 annot_keep$id, annot_keep$SYMBOL))
  idx <- match(rk, name_key)

  keep_annot <- tibble(id = annot_keep$id[idx], SYMBOL = annot_keep$SYMBOL[idx])

  # (A) TMM-CPM
  cpm_mat <- edgeR::cpm(dge, normalized.lib.sizes = TRUE)
  readr::write_csv(bind_cols(keep_annot, as.data.frame(cpm_mat, check.names=FALSE)),
                   file.path(cfg$paths$results, "CPM_filtered_annot_POST.csv"))

  # (B) voom log-CPM
  voomE <- v$E[match(rk, rownames(v$E)),, drop = FALSE]
  readr::write_csv(bind_cols(keep_annot, as.data.frame(voomE, check.names=FALSE)),
                   file.path(cfg$paths$results, "voom_logCPM_filtered_annot_POST.csv"))

  # (C) Kept raw counts
  readr::write_csv(bind_cols(keep_annot, as.data.frame(dge$counts, check.names=FALSE)),
                   file.path(cfg$paths$results, "raw_counts_filtered_annot_POST.csv"))
}
