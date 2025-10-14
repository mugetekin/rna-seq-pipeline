# scripts/02_de_limma.R
source("R/io_helpers.R")
suppressPackageStartupMessages({ library(limma); library(tidyverse); library(readr) })

`%||%` <- function(a, b) if (is.null(a) || !length(a)) b else a

cfg <- load_cfg()

# read voom+fit from step 01
obj  <- readRDS(file.path(cfg$paths$rds, "norm_and_fit.rds"))
fit2 <- obj$fit2

# Helper: normalize names like "Lo - PBS" -> "Lo_vs_PBS"
norm_name <- function(x) gsub(" - ", "_vs_", x, fixed = TRUE)

# Available coefficients in the fitted model
avail_raw  <- colnames(fit2$coefficients)              # e.g., "Lo - PBS", "Lo - Med", ...
avail_norm <- norm_name(avail_raw)                      # e.g., "Lo_vs_PBS", "Lo_vs_Med", ...
names(avail_raw) <- avail_norm                          # map normalized -> raw

# Requested contrasts (from config) or default (ALL available)
want <- cfg$params$contrasts
if (is.null(want) || !length(want)) {
  message("No contrasts specified in config; exporting ALL available: ",
          paste(avail_norm, collapse = ", "))
  want <- avail_norm
}

# Keep only those that actually exist
use_norm <- intersect(want, avail_norm)
if (!length(use_norm)) {
  stop("No requested contrasts found in fit2.\nRequested: ",
       paste(want, collapse = ", "),
       "\nAvailable: ", paste(avail_norm, collapse = ", "))
}

# thresholds from config (with safe defaults)
lfc_thr <- cfg$params$lfc_thresh %||% 1
fdr_thr <- cfg$params$fdr_thresh %||% 0.05
adj_mth <- "BH"

# Prepare outputs
dir.create(cfg$paths$results, showWarnings = FALSE, recursive = TRUE)

# Compute & write each DE table
de_list <- setNames(vector("list", length(use_norm)), use_norm)
for (cn_norm in use_norm) {
  coef_raw <- avail_raw[[cn_norm]]   # real coef name in fit2 (e.g., "Lo - Med")
  tt <- topTable(
    fit2,
    coef          = coef_raw,
    number        = Inf,
    sort.by       = "P",
    lfc           = lfc_thr,
    p.value       = fdr_thr,
    adjust.method = adj_mth
  )
  # add Gene column from rownames for tidy TSV
  if (!tibble::has_name(tt, "Gene")) {
    tt <- tt %>% tibble::rownames_to_column("Gene")
  }
  out_tsv <- file.path(cfg$paths$results, paste0("DE_", cn_norm, ".tsv"))
  readr::write_tsv(tt, out_tsv)
  de_list[[cn_norm]] <- tt
  message("Wrote: ", out_tsv, "  (coef: ", coef_raw, ")")
}

# Back-compat for 03/04: expose tt_lo/tt_med/tt_hi if present
tt_lo  <- if ("Lo_vs_PBS"  %in% names(de_list)) de_list[["Lo_vs_PBS"]]  else NULL
tt_med <- if ("Med_vs_PBS" %in% names(de_list)) de_list[["Med_vs_PBS"]] else NULL
tt_hi  <- if ("Hi_vs_PBS"  %in% names(de_list)) de_list[["Hi_vs_PBS"]]  else NULL

# Save RDS with both styles
saveRDS(
  c(list(tt_lo = tt_lo, tt_med = tt_med, tt_hi = tt_hi), de_list),
  file.path(cfg$paths$rds, "de_tables.rds")
)

message("Done. Exported contrasts: ", paste(names(de_list), collapse = ", "))
