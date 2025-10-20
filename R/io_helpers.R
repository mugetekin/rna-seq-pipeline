#' - "counts": raw read counts per gene from RNA-seq (integer table).
#' - "CPM": counts per million; a way to compare samples after scaling for library size.
#' - "TMM": trimmed mean of M-values; a robust way to normalize RNA-seq libraries.
#' - "voom": transforms counts to logCPM with precision weights for linear modeling.
#' -----------

suppressPackageStartupMessages({
  library(data.table); library(tidyverse); library(yaml)
})

#' Load configuration file (YAML format)
#' Reads paths and parameters (e.g., contrasts, thresholds) from config.yaml.
load_cfg <- function(path = "config.yaml") yaml::read_yaml(path)


#' This function reads a tab-delimited table containing gene-level counts and annotations.
#' It automatically detects ID, gene symbol, and biotype columns and extracts sample columns.
#' Returns a list with:
#'   - annot: annotation dataframe (gene IDs, symbols, types)
#'   - counts: numeric matrix of counts (genes × samples)
read_counts_annot <- function(counts_path) {
  # Read the entire count table as a data.frame
  dt <- data.table::fread(counts_path, data.table = FALSE)
  nms <- names(dt)

  # Helper: find column names case-insensitively
  pick <- function(cands) {
    ix <- which(tolower(nms) %in% tolower(cands))
    if (length(ix)) nms[ix[1]] else NA_character_
  }

  # Detect typical annotation columns
  id_col   <- pick(c("id","gene_id","gene","ensembl_id","Gene"))
  sym_col  <- pick(c("SYMBOL","symbol","gene_name"))
  type_col <- pick(c("gene_type","biotype","gene_biotype"))
  stopifnot(!is.na(id_col) | !is.na(sym_col))

  # Extract annotation columns
  annot <- tibble(
    id        = if (!is.na(id_col))   dt[[id_col]]   else NA_character_,
    SYMBOL    = if (!is.na(sym_col))  dt[[sym_col]]  else NA_character_,
    gene_type = if (!is.na(type_col)) dt[[type_col]] else NA_character_
  )

  # Define rownames for counts:
  # Prefer gene symbols; fall back to IDs if symbol missing.
  genes <- annot$SYMBOL
  genes[is.na(genes) | genes == ""] <- annot$id
  genes <- make.unique(genes)

  # Identify all sample columns (everything not part of annotation)
  sample_cols <- setdiff(nms, unique(na.omit(c(id_col, sym_col, type_col))))
  counts <- dt[, sample_cols, drop = FALSE]

  # Ensure numeric matrix for counts (coerce factors/characters safely)
  counts[] <- lapply(counts, function(x) as.numeric(as.character(x)))
  rownames(counts) <- genes
  # Sanity check: ensure all count values are numeric and finite
  stopifnot(all(is.finite(as.matrix(counts))))

  list(annot = annot, counts = counts)
}


#' Generates a metadata table from column names of the count matrix.
#' Automatically assigns experimental groups based on name prefixes:
#'   - "PBS", "Lo", "Med", "Hi" → mapped to factors
#' Any other prefix defaults to "Other".
#' Returns a tibble with columns:
#'   - SampleID: sanitized sample name
#'   - Group: factor of experimental group
make_meta_from_colnames <- function(cn) {
  raw <- tibble(
    SampleID = make.names(cn),
    Group = dplyr::case_when(
      startsWith(cn, "PBS") ~ "PBS",
      startsWith(cn, "Lo")  ~ "Lo",
      startsWith(cn, "Med") ~ "Med",
      startsWith(cn, "Hi")  ~ "Hi",
      TRUE ~ "Other"
    )
  )

  # Drop unused levels; keep only levels that actually appear
  present <- unique(raw$Group)
  wanted  <- intersect(c("PBS","Lo","Med","Hi","Other"), present)
  raw$Group <- factor(raw$Group, levels = wanted)
  raw
}


#' Creates all necessary directories defined in the config.yaml file.
#' This ensures that subsequent steps (saving plots, RDS files, results)
#' do not fail due to missing folders.
                     
dir_prep <- function(cfg) {
  dirs <- c(cfg$paths$outputs, cfg$paths$rds, cfg$paths$results, cfg$paths$figures)
  for (d in dirs) dir.create(d, showWarnings = FALSE, recursive = TRUE)
}
