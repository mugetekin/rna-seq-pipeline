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

  # numeric + log2(CPM+1) + basic cleaning
  expr_num <- as.matrix(data.frame(lapply(as.data.frame(expr), as.numeric),
                                   check.names = FALSE, row.names = rownames(expr)))
  expr_log <- log2(expr_num + 1)
  keep_rows <- apply(is.finite(expr_log), 1, all)
  keep_cols <- apply(is.finite(expr_log), 2, all)
  expr_log  <- expr_log[keep_rows, keep_cols, drop = FALSE]
  nzv <- apply(expr_log, 1, sd) > 0
  expr_log <- expr_log[nzv, , drop = FALSE]

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
