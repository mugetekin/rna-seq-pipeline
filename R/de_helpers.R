suppressPackageStartupMessages({
  library(edgeR); library(limma); library(tidyverse)
})

# edgeR best practice: filterByExpr + TMM normalization
do_normalize <- function(counts, meta) {
  # Drop unused factor levels
  meta <- meta %>% mutate(Group = droplevels(Group))
  dge <- DGEList(counts = as.matrix(counts), samples = as.data.frame(meta))
  
  # Filter low-expressed genes and normalize library sizes
  keep <- filterByExpr(dge, group = meta$Group)
  dge  <- dge[keep,, keep.lib.sizes = FALSE]
  dge  <- calcNormFactors(dge, method = "TMM")
  dge
}

# limma-voom model fitting + contrasts (handles reduced factor levels)
fit_limma <- function(dge, meta) {
  meta <- meta %>% mutate(Group = droplevels(Group))
  design <- model.matrix(~ 0 + Group, data = meta)
  colnames(design) <- levels(meta$Group)

  v <- voom(dge, design, plot = FALSE)

  # Automatically create contrasts among available groups
  has <- colnames(design)                 # e.g., PBS, Lo, Med, Hi present in this dataset
  mk  <- function(a,b) if (all(c(a,b) %in% has)) paste0(a," - ",b) else NA_character_
  
  # 1) classic vs-PBS contrasts
  contrs_vs_ref <- c(
    Lo_vs_PBS  = mk("Lo","PBS"),
    Med_vs_PBS = mk("Med","PBS"),
    Hi_vs_PBS  = mk("Hi","PBS")
  )
  
  # 2) all pairwise among non-reference levels (Lo, Med, Hi)
  non_ref <- setdiff(has, "PBS")
  if (length(non_ref) >= 2) {
    pw <- combn(non_ref, 2, simplify = FALSE)
    contrs_pw <- setNames(
      vapply(pw, function(p) mk(p[1], p[2]), character(1)),
      vapply(pw, function(p) paste0(p[1], "_vs_", p[2]), character(1))
    )
  } else {
    contrs_pw <- c()
  }
  
  contrs <- c(contrs_vs_ref, contrs_pw)
  contrs <- contrs[!is.na(contrs)]
  
  if (!length(contrs))
    stop("No valid contrasts with available groups: ", paste(has, collapse = ", "))
  
  cm <- makeContrasts(contrasts = unname(contrs), levels = design)
  colnames(cm) <- names(contrs)           # keep friendly names like Lo_vs_Med
  
  fit  <- lmFit(v, design)
  fit2 <- eBayes(contrasts.fit(fit, cm), trend = TRUE)
  list(v = v, fit2 = fit2)
  
}

# Write post-normalization snapshots
write_post_norm <- function(dge, v, annot_keep, cfg) {
  rk <- rownames(dge$counts)
  name_key <- make.unique(ifelse(is.na(annot_keep$SYMBOL) | annot_keep$SYMBOL=="",
                                 annot_keep$id, annot_keep$SYMBOL))
  idx <- match(rk, name_key)
  keep_annot <- tibble(id = annot_keep$id[idx], SYMBOL = annot_keep$SYMBOL[idx])

  # Export normalized CPMs
  cpm_mat <- edgeR::cpm(dge, normalized.lib.sizes = TRUE)
  readr::write_csv(bind_cols(keep_annot, as.data.frame(cpm_mat, check.names=FALSE)),
                   file.path(cfg$paths$results, "CPM_filtered_annot_POST.csv"))

  # Export voom-transformed logCPMs
  voomE <- v$E[match(rk, rownames(v$E)),, drop = FALSE]
  readr::write_csv(bind_cols(keep_annot, as.data.frame(voomE, check.names=FALSE)),
                   file.path(cfg$paths$results, "voom_logCPM_filtered_annot_POST.csv"))

  # Export filtered raw counts
  readr::write_csv(bind_cols(keep_annot, as.data.frame(dge$counts, check.names=FALSE)),
                   file.path(cfg$paths$results, "raw_counts_filtered_annot_POST.csv"))
}
