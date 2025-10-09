source("R/io_helpers.R")
suppressPackageStartupMessages({ library(limma); library(tidyverse) })

cfg <- load_cfg()

# read voom+fit from step 01
obj  <- readRDS(file.path(cfg$paths$rds, "norm_and_fit.rds"))
fit2 <- obj$fit2

# Helper: normalize names like "Lo - PBS" -> "Lo_vs_PBS"
norm_name <- function(x) gsub(" - ", "_vs_", x, fixed = TRUE)

avail_raw  <- colnames(fit2$coefficients)
avail_norm <- norm_name(avail_raw)

# Requested contrasts: from config or default
want <- cfg$params$contrasts
if (is.null(want) || !length(want)) want <- c("Lo_vs_PBS","Med_vs_PBS","Hi_vs_PBS")

# Match requested against available (normalized)
use_norm <- intersect(want, avail_norm)
if (!length(use_norm)) {
  stop("No requested contrasts found in fit2. Available (normalized): ",
       paste(avail_norm, collapse=", "))
}

# Map normalized -> raw index
idx <- match(use_norm, avail_norm)

# Prepare outputs
dir.create(cfg$paths$results, showWarnings = FALSE, recursive = TRUE)

# Compute & write each DE table
de_list <- setNames(vector("list", length(use_norm)), use_norm)
for (i in seq_along(use_norm)) {
  coef_raw <- avail_raw[idx[i]]   # real coef name in fit2
  cn_norm  <- use_norm[i]         # normalized name for files/keys
  tt <- topTable(fit2, coef = coef_raw, number = Inf, sort.by = "P")
  out_tsv <- file.path(cfg$paths$results, paste0("DE_", cn_norm, ".tsv"))
  write.table(tt, out_tsv, sep = "\t", quote = FALSE)
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
