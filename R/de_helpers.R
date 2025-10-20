# Quietly load required R packages so the console isn't spammed
suppressPackageStartupMessages({
  library(edgeR); library(limma); library(tidyverse)  # RNA-seq DE pipeline + data wrangling
})

# edgeR best practice: filter lowly expressed genes, then TMM-normalize library sizes
do_normalize <- function(counts, meta) {
  meta <- meta %>% mutate(Group = droplevels(Group))               # drop unused groups (clean factors)
  dge  <- DGEList(counts = as.matrix(counts), samples = as.data.frame(meta))  # edgeR container
  keep <- filterByExpr(dge, group = meta$Group)                    # keep genes expressed enough to be informative
  dge  <- dge[keep, , keep.lib.sizes = FALSE]                      # subset to kept genes; reset lib sizes to recompute
  dge  <- calcNormFactors(dge, method = "TMM")                     # TMM normalization (makes libraries comparable)
  dge                                                                # return normalized DGEList
}

# Build the design matrix (which groups we compare) for linear modeling
make_design <- function(meta, ref_level = NULL) {
  g <- droplevels(meta$Group)                                      # ensure only present groups are kept
  if (!is.null(ref_level) && ref_level %in% levels(g)) {
    g <- relevel(g, ref = ref_level)                               # set reference (e.g., "PBS") so contrasts are intuitive
  }
  design <- model.matrix(~ 0 + g)                                  # no intercept; get one column per group
  colnames(design) <- levels(g)                                    # name columns as group names
  design
}

# voom-transform (logCPM + precision weights) and fit linear model
voom_fit <- function(dge, design) {
  v   <- voom(dge, design, plot = FALSE)                           # convert counts to logCPM with mean-variance weights
  fit <- lmFit(v, design)                                          # fit linear model per gene
  list(v = v, fit = fit)
}

# Create contrast matrix from vector of "A_vs_B" strings like "Lo_vs_PBS"
make_contrast_matrix <- function(contrasts, design) {
  # parse "A_vs_B" into "A - B", ensuring the groups exist in design
  to_expr <- function(x) {
    parts <- strsplit(x, "_vs_")[[1]]
    stopifnot(length(parts) == 2)
    paste0(parts[1], " - ", parts[2])
  }
  contrast_exprs <- sapply(contrasts, to_expr, USE.NAMES = TRUE)
  makeContrasts(contrasts = contrast_exprs, levels = design)       # limma contrast matrix
}

# Apply contrasts, borrow strength with eBayes, return fitted object
fit_contrasts_ebayes <- function(fit, contrast_matrix) {
  fit2 <- contrasts.fit(fit, contrast_matrix)                       # apply requested comparisons
  fit2 <- eBayes(fit2)                                              # empirical Bayes: stabilizes variance estimates
  fit2
}

# Get ranked tables per contrast with useful statistics (logFC, FDR, etc.)
top_tables <- function(fit2, lfc_thresh = 1, fdr_thresh = 0.05, n = Inf) {
  # Return a named list of data frames, one per contrast
  res <- lapply(colnames(fit2), function(ct) {
    tt <- topTable(fit2, coef = ct, number = n, sort.by = "P")      # ranked by p-value
    tt <- tt %>%
      rownames_to_column("gene") %>%
      mutate(sig = if_else(adj.P.Val <= fdr_thresh & abs(logFC) >= lfc_thresh, TRUE, FALSE))  # simple flag used in plots
    tt
  })
  names(res) <- colnames(fit2)
  res
}

# Attach gene annotation (SYMBOL/ID) to top tables for human-readable outputs
add_annotation_to_toptables <- function(tt_list, annot) {
  # Expect annot to have columns: id, SYMBOL; rownames(tt) are gene names from normalization step
  lapply(tt_list, function(tt) {
    # try to match by SYMBOL first, then by id as backup
    left_join(tt, annot %>% select(id, SYMBOL), by = c("gene" = "SYMBOL")) %>%
      mutate(GeneLabel = if_else(!is.na(gene), gene, id)) %>%       # prefer the readable gene symbol
      relocate(GeneLabel, .before = gene)
  })
}

# Write normalized CPMs and voom logCPMs for downstream inspection and plotting
write_post_norm <- function(dge, v, annot_keep, cfg) {
  rk <- rownames(dge$counts)                                       # final kept gene names (post-filter)
  # Build a stable "name key": prefer SYMBOL; fall back to ID; ensure uniqueness
  name_key <- make.unique(ifelse(is.na(annot_keep$SYMBOL) | annot_keep$SYMBOL == "",
                                 annot_keep$id, annot_keep$SYMBOL))
  idx <- match(rk, name_key)                                       # align annotation to kept rows
  keep_annot <- tibble(id = annot_keep$id[idx], SYMBOL = annot_keep$SYMBOL[idx])

  # Export normalized CPMs (not log) — easy to read/report
  cpm_mat <- edgeR::cpm(dge, normalized.lib.sizes = TRUE)          # CPM using TMM-normalized library sizes
  readr::write_csv(bind_cols(keep_annot, as.data.frame(cpm_mat, check.names = FALSE)),
                   file.path(cfg$paths$results, "CPM_filtered_annot_POST.csv"))  # CSV under outputs/results

  # Export voom-transformed logCPMs (used by limma models/plots)
  voomE <- v$E[match(rk, rownames(v$E)), , drop = FALSE]           # align voom matrix to kept genes
  readr::write_csv(bind_cols(keep_annot, as.data.frame(voomE, check.names = FALSE)),
                   file.path(cfg$paths$results, "voom_logCPM_filtered_annot_POST.csv"))
}

# Convenience wrapper: complete DE workflow from counts+meta to top tables
run_de_pipeline <- function(counts, meta, contrasts, ref_level = "PBS",
                            lfc_thresh = 1, fdr_thresh = 0.05) {
  dge          <- do_normalize(counts, meta)                        # filter + TMM
  design       <- make_design(dge$samples, ref_level = ref_level)   # design matrix (one column per group)
  vf           <- voom_fit(dge, design)                             # voom + lmFit
  cm           <- make_contrast_matrix(contrasts, design)           # construct A-B contrasts
  fit2         <- fit_contrasts_ebayes(vf$fit, cm)                  # apply contrasts + eBayes
  tt_list      <- top_tables(fit2, lfc_thresh = lfc_thresh, fdr_thresh = fdr_thresh)  # ranked gene lists
  list(dge = dge, design = design, v = vf$v, fit = vf$fit, fit2 = fit2, tables = tt_list)
}

# Optional: write DE result tables to disk (one CSV per contrast)
write_de_tables <- function(tt_list, out_dir) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)       # ensure destination exists
  for (nm in names(tt_list)) {
    fn <- file.path(out_dir, paste0("DE_", nm, "_results.csv"))     # e.g., DE_Lo_vs_PBS_results.csv
    readr::write_csv(tt_list[[nm]], fn)
  }
}
