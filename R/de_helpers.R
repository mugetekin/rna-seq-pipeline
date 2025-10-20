#' This script implements the main preprocessing and modeling pipeline for bulk RNA-seq:
#'   1. edgeR-based filtering and TMM normalization
#'   2. limma-voom linear modeling and contrast generation
#'   3. Export of normalized matrices for inspection and downstream analyses
#'
#' Key methods:
#' - edgeR::filterByExpr : removes low-expressed genes safely based on design
#' - edgeR::calcNormFactors (TMM): normalizes for library size differences
#' - limma::voom : transforms counts to logCPM for linear modeling
#' - limma::eBayes : moderates standard errors for small-sample inference

suppressPackageStartupMessages({
  library(edgeR); library(limma); library(tidyverse)
})



#' Performs filtering and TMM normalization using edgeR.
#' Input:
#'   - counts: numeric matrix of raw read counts (genes × samples)
#'   - meta: sample metadata (must include a 'Group' column)
#' Output:
#'   - edgeR DGEList object, filtered and normalized
#'
#' Steps:
#'   1. Convert to DGEList (edgeR’s main container for count data)
#'   2. Filter genes with low expression (filterByExpr)
#'   3. Recalculate library sizes and normalize (TMM)
do_normalize <- function(counts, meta) {
  # Drop unused factor levels (important for edgeR design construction)
  meta <- meta %>% mutate(Group = droplevels(Group))
  dge <- DGEList(counts = as.matrix(counts), samples = as.data.frame(meta))
  
  # Filter low-expressed genes and normalize library sizes
  keep <- filterByExpr(dge, group = meta$Group)
  dge  <- dge[keep,, keep.lib.sizes = FALSE]
  # Compute TMM normalization factors to scale library sizes
  dge  <- calcNormFactors(dge, method = "TMM")
  dge
}

#' Fits a linear model using limma-voom, automatically building contrasts.
#' Input:
#'   - dge: normalized DGEList from edgeR
#'   - meta: sample metadata with factor Group
#' Output:
#'   - list(v, fit2): voom-transformed expression object + limma fit with eBayes moderation
#'
#' Workflow:
#'   1. Create a design matrix without intercept (~0 + Group)
#'   2. Transform counts with voom (mean-variance modeling)
#'   3. Automatically generate contrasts:
#'        a) each treatment vs PBS (reference)
#'        b) all pairwise among non-reference groups (e.g., Lo vs Med)
#'   4. Fit linear model and apply empirical Bayes shrinkage
fit_limma <- function(dge, meta) {
  meta <- meta %>% mutate(Group = droplevels(Group))
  design <- model.matrix(~ 0 + Group, data = meta)
  colnames(design) <- levels(meta$Group)  # e.g., PBS, Lo, Med, Hi

  # Transform counts to logCPM + precision weights for linear modeling
  v <- voom(dge, design, plot = FALSE)

  # Automatically create contrasts among available groups
  has <- colnames(design)                 # e.g., PBS, Lo, Med, Hi present in this dataset
  mk  <- function(a,b) if (all(c(a,b) %in% has)) paste0(a," - ",b) else NA_character_
  
  # classic vs-PBS contrasts
  contrs_vs_ref <- c(
    Lo_vs_PBS  = mk("Lo","PBS"),
    Med_vs_PBS = mk("Med","PBS"),
    Hi_vs_PBS  = mk("Hi","PBS")
  )
  
  # all pairwise among non-reference levels (Lo, Med, Hi)
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

  # Combine valid contrasts and remove NAs
  contrs <- c(contrs_vs_ref, contrs_pw)
  contrs <- contrs[!is.na(contrs)]

  # Sanity check: stop if no valid contrasts exist (e.g., single-group dataset)
  if (!length(contrs))
    stop("No valid contrasts with available groups: ", paste(has, collapse = ", "))

  # Build contrast matrix for limma
  cm <- makeContrasts(contrasts = unname(contrs), levels = design)
  colnames(cm) <- names(contrs)           # keep friendly names like Lo_vs_Med

  # Fit model and apply empirical Bayes moderation
  fit  <- lmFit(v, design)
  fit2 <- eBayes(contrasts.fit(fit, cm), trend = TRUE)
  # Return both voom object and moderated fit
  list(v = v, fit2 = fit2)
  
}

             
#' Saves normalized data for transparency and downstream QC.
#' Input:
#'   - dge: edgeR DGEList (filtered + normalized)
#'   - v: voom object with logCPM
#'   - annot_keep: annotation table (id, SYMBOL columns)
#'   - cfg: configuration with output paths
#'
#' Output files (CSV under cfg$paths$results):
#'   - CPM_filtered_annot_POST.csv      : normalized CPMs
#'   - voom_logCPM_filtered_annot_POST.csv : voom-transformed logCPMs
#'   - raw_counts_filtered_annot_POST.csv  : filtered raw counts
write_post_norm <- function(dge, v, annot_keep, cfg) {
  rk <- rownames(dge$counts)
  # Safely reconstruct gene names from annotation
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
