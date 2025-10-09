suppressPackageStartupMessages({
  library(edgeR); library(limma); library(tidyverse)
})

# ---- edgeR önerisi: filterByExpr + TMM ----
do_normalize <- function(counts, meta) {
  meta <- meta %>% mutate(Group = droplevels(Group)) # kullanılmayan seviyeleri düşür
  dge <- DGEList(counts = as.matrix(counts), samples = as.data.frame(meta))
  keep <- filterByExpr(dge, group = meta$Group)
  dge  <- dge[keep,, keep.lib.sizes = FALSE]
  dge  <- calcNormFactors(dge, method = "TMM")
  dge
}

# ---- limma-voom + contrasts (seviyeler daraltılmış) ----
fit_limma <- function(dge, meta) {
  meta <- meta %>% mutate(Group = droplevels(Group))
  design <- model.matrix(~ 0 + Group, data = meta)
  colnames(design) <- levels(meta$Group)

  v <- voom(dge, design, plot = FALSE)

  # mevcut olmayan seviyeler otomatik düşer; contrasts olanları varsa kur
  has <- colnames(design)
  mk <- function(a,b) if (all(c(a,b) %in% has)) paste0(a," - ",b) else NA_character_

  contrs <- c(
    Lo_vs_PBS  = mk("Lo","PBS"),
    Med_vs_PBS = mk("Med","PBS"),
    Hi_vs_PBS  = mk("Hi","PBS")
  )
  contrs <- contrs[!is.na(contrs)]
  if (!length(contrs)) stop("No valid contrasts with available groups: ", paste(has, collapse=", "))

  cm <- makeContrasts(contrasts = unname(contrs), levels = design)
  names(cm) <- names(contrs)

  fit  <- lmFit(v, design)
  fit2 <- eBayes(contrasts.fit(fit, cm))
  list(v = v, fit2 = fit2)
}

# ---- Snapshot writers ----
write_post_norm <- function(dge, v, annot_keep, cfg) {
  rk <- rownames(dge$counts)
  name_key <- make.unique(ifelse(is.na(annot_keep$SYMBOL) | annot_keep$SYMBOL=="",
                                 annot_keep$id, annot_keep$SYMBOL))
  idx <- match(rk, name_key)
  keep_annot <- tibble(id = annot_keep$id[idx], SYMBOL = annot_keep$SYMBOL[idx])

  cpm_mat <- edgeR::cpm(dge, normalized.lib.sizes = TRUE)
  readr::write_csv(bind_cols(keep_annot, as.data.frame(cpm_mat, check.names=FALSE)),
                   file.path(cfg$paths$results, "CPM_filtered_annot_POST.csv"))

  voomE <- v$E[match(rk, rownames(v$E)),, drop = FALSE]
  readr::write_csv(bind_cols(keep_annot, as.data.frame(voomE, check.names=FALSE)),
                   file.path(cfg$paths$results, "voom_logCPM_filtered_annot_POST.csv"))

  readr::write_csv(bind_cols(keep_annot, as.data.frame(dge$counts, check.names=FALSE)),
                   file.path(cfg$paths$results, "raw_counts_filtered_annot_POST.csv"))
}
