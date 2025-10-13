suppressPackageStartupMessages({
  library(data.table); library(tidyverse); library(yaml)
})

load_cfg <- function(path = "config.yaml") yaml::read_yaml(path)

read_counts_annot <- function(counts_path) {
  dt <- data.table::fread(counts_path, data.table = FALSE)
  nms <- names(dt)

  pick <- function(cands) {
    ix <- which(tolower(nms) %in% tolower(cands))
    if (length(ix)) nms[ix[1]] else NA_character_
  }

  id_col   <- pick(c("id","gene_id","gene","ensembl_id","Gene"))
  sym_col  <- pick(c("SYMBOL","symbol","gene_name"))
  type_col <- pick(c("gene_type","biotype","gene_biotype"))
  stopifnot(!is.na(id_col) | !is.na(sym_col))

  annot <- tibble(
    id        = if (!is.na(id_col))   dt[[id_col]]   else NA_character_,
    SYMBOL    = if (!is.na(sym_col))  dt[[sym_col]]  else NA_character_,
    gene_type = if (!is.na(type_col)) dt[[type_col]] else NA_character_
  )

  genes <- annot$SYMBOL
  genes[is.na(genes) | genes == ""] <- annot$id
  genes <- make.unique(genes)

  sample_cols <- setdiff(nms, unique(na.omit(c(id_col, sym_col, type_col))))
  counts <- dt[, sample_cols, drop = FALSE]

  # Ensure numeric matrix for counts (coerce factors/characters safely)
  counts[] <- lapply(counts, function(x) as.numeric(as.character(x)))
  rownames(counts) <- genes
  stopifnot(all(is.finite(as.matrix(counts))))

  list(annot = annot, counts = counts)
}

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

dir_prep <- function(cfg) {
  dirs <- c(cfg$paths$outputs, cfg$paths$rds, cfg$paths$results, cfg$paths$figures)
  for (d in dirs) dir.create(d, showWarnings = FALSE, recursive = TRUE)
}
