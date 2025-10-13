suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
})
source("R/utils_normed.R")

out_dir <- "outputs_normed"
dir.create(out_dir, showWarnings = FALSE)

# --- load normalized log2-CPM and groups (PBS/Lo/Med/Hi) ---
obj <- load_normed_expr("data/normed_cpms_filtered_annot.csv")
X    <- obj$expr_log             # log2(CPM+1)
grp  <- factor(obj$grp, levels = c("PBS","Lo","Med","Hi"))
samples <- colnames(X)

# --- collapse duplicate gene names by mean (DE should not have duplicate rows) ---
# keep original rownames as "Gene", average by Gene across samples
dd <- as.data.frame(X) |>
  tibble::rownames_to_column("Gene") |>
  dplyr::group_by(Gene) |>
  dplyr::summarise(across(all_of(samples), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

# back to matrix with unique genes
mat <- dd |>
  tibble::column_to_rownames("Gene") |>
  as.matrix()

# --- filter: require >= 2 finite values per group (limma can handle NAs, but this is sane QC) ---
ok_by_grp <- sapply(levels(grp), function(g){
  idx <- which(grp == g)
  apply(mat, 1, function(x) sum(is.finite(x[idx])) >= 2)
})
keep <- apply(ok_by_grp, 1, all)
mat  <- mat[keep, , drop = FALSE]

# --- design matrix and contrasts ---
design <- model.matrix(~ 0 + grp)        # no intercept
colnames(design) <- levels(grp)

fit  <- lmFit(mat, design)
ct   <- makeContrasts(
  Lo_vs_PBS  = Lo  - PBS,
  Med_vs_PBS = Med - PBS,
  Hi_vs_PBS  = Hi  - PBS,
  levels = design
)
fit2 <- contrasts.fit(fit, ct)
fit2 <- eBayes(fit2)

# helper to extract and tidy a contrast
extract_contrast <- function(fit2, coef_name) {
  tt <- topTable(fit2, coef = coef_name, number = Inf, sort.by = "P")
  tt <- tt |>
    tibble::rownames_to_column("Gene") |>
    dplyr::transmute(
      Gene,
      log2FC  = logFC,
      AveExpr = AveExpr,
      t       = t,
      P.Value = P.Value,
      padj    = adj.P.Val,
      Contrast = coef_name
    )
  tt
}

res <- bind_rows(
  extract_contrast(fit2, "Lo_vs_PBS"),
  extract_contrast(fit2, "Med_vs_PBS"),
  extract_contrast(fit2, "Hi_vs_PBS")
)

# write combined results
out_csv <- file.path(out_dir, "DE_results_example.csv")
readr::write_csv(res, out_csv)

message("Saved: ", out_csv)
# optional: quick summary to console
print(res |>
  mutate(sig = padj < 0.05) |>
  group_by(Contrast, Direction = if_else(log2FC>0,"Up","Down")) |>
  summarise(n = sum(sig), .groups="drop"))
