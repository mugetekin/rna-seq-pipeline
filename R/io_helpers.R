#' - "counts": raw read counts per gene from RNA-seq (integer table).
#' - "CPM": counts per million; a way to compare samples after scaling for library size.
#' - "TMM": trimmed mean of M-values; a robust way to normalize RNA-seq libraries.
#' - "voom": transforms counts to logCPM with precision weights for linear modeling.
#' -----------

# Quietly load required R packages so the console isn't spammed.
suppressPackageStartupMessages({
# Load an external package that provides functions we use below.
  library(data.table); library(tidyverse); library(yaml)  # basic data wrangling + YAML parsing
})

# Load analysis configuration (paths, parameters) from config.yaml
load_cfg <- function(path = "config.yaml") yaml::read_yaml(path)  # returns a list like cfg$paths, cfg$params

# Read annotated raw counts file (genes × samples) and split into annotation + counts matrix
read_counts_annot <- function(counts_path) {
  dt <- data.table::fread(counts_path, data.table = FALSE)  # fast import of large tab-delimited table
  nms <- names(dt)  # store column names for convenience

  # helper to find which column name matches known candidates (case-insensitive)
  pick <- function(cands) {
    ix <- which(tolower(nms) %in% tolower(cands))
    if (length(ix)) nms[ix[1]] else NA_character_
  }
  
  # try to locate annotation columns (these vary between datasets)
  id_col   <- pick(c("id","gene_id","gene","ensembl_id","Gene"))
  sym_col  <- pick(c("SYMBOL","symbol","gene_name"))
  type_col <- pick(c("gene_type","biotype","gene_biotype"))
  stopifnot(!is.na(id_col) | !is.na(sym_col))

  # make a clean annotation table
  annot <- tibble(
    id        = if (!is.na(id_col))   dt[[id_col]]   else NA_character_,
    SYMBOL    = if (!is.na(sym_col))  dt[[sym_col]]  else NA_character_,
    gene_type = if (!is.na(type_col)) dt[[type_col]] else NA_character_
  )

  # pick SYMBOLs as gene names, fall back to IDs when missing; ensure uniqueness
  genes <- annot$SYMBOL
  genes[is.na(genes) | genes == ""] <- annot$id
  genes <- make.unique(genes)

  # everything else are the count columns (samples)
  sample_cols <- setdiff(nms, unique(na.omit(c(id_col, sym_col, type_col))))
  counts <- dt[, sample_cols, drop = FALSE]

  # Ensure numeric matrix for counts (coerce factors/characters safely)
  counts[] <- lapply(counts, function(x) as.numeric(as.character(x)))
  rownames(counts) <- genes
  stopifnot(all(is.finite(as.matrix(counts))))  # sanity check: no NA/Inf

  list(annot = annot, counts = counts)  # return both pieces
}

# Create sample metadata (sample ID + group) from column names
make_meta_from_colnames <- function(cn) {
  raw <- tibble(
    SampleID = make.names(cn),  # valid R names for sample IDs
    Group = dplyr::case_when(   # infer group from prefix convention
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

# Make sure all output directories defined in config exist (creates them if missing)
dir_prep <- function(cfg) {
  dirs <- c(cfg$paths$outputs, cfg$paths$rds, cfg$paths$results, cfg$paths$figures)
  for (d in dirs) dir.create(d, showWarnings = FALSE, recursive = TRUE)
}
