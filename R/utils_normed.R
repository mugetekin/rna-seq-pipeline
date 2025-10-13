suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
})

# read normalized CPM file, build log2 matrix and groups
load_normed_expr <- function(csv = "data/normed_cpms_filtered_annot.csv") {
  stopifnot(file.exists(csv))
  dat <- data.table::fread(csv, data.table = FALSE)

  # expression matrix with SYMBOL fallback
  gene_cols <- c("id","SYMBOL","gene_type")
  expr <- dat[, !(names(dat) %in% gene_cols), drop = FALSE]
  genes <- ifelse(is.na(dat$SYMBOL) | dat$SYMBOL=="", dat$id, dat$SYMBOL)
  rownames(expr) <- make.unique(genes)

  # numeric + log2(CPM+1)
  expr_num <- as.matrix(data.frame(lapply(as.data.frame(expr), as.numeric),
                                   check.names = FALSE, row.names = rownames(expr)))
  expr_log <- log2(expr_num + 1)

  # ---- KEEP columns even if they have some NAs ----
  # keep rows (genes) that have at least 2 finite values (so var is defined)
  keep_rows <- apply(expr_log, 1, function(x) sum(is.finite(x)) >= 2)
  expr_log  <- expr_log[keep_rows, , drop = FALSE]

  # keep columns that have at least 1 finite value
  keep_cols <- apply(expr_log, 2, function(x) any(is.finite(x)))
  expr_log  <- expr_log[, keep_cols, drop = FALSE]

  # groups from sample id prefix
  sample_ids <- colnames(expr_log)
  grp <- dplyr::case_when(
    startsWith(sample_ids,"PBS") ~ "PBS",
    startsWith(sample_ids,"Lo")  ~ "Lo",
    startsWith(sample_ids,"Med") ~ "Med",
    startsWith(sample_ids,"Hi")  ~ "Hi",
    TRUE ~ "Other"
  )
  grp <- factor(grp, levels = c("PBS","Lo","Med","Hi","Other"))

  list(expr_log = expr_log, grp = grp)
}
