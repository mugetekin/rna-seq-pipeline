#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(yaml)
  library(edgeR)
  library(limma)
})

`%||%` <- function(a,b) if (is.null(a)) b else a
safe_dir <- function(...) { d <- file.path(...); dir.create(d, showWarnings=FALSE, recursive=TRUE); d }
load_cfg <- function(path="config.yaml") yaml::read_yaml(path)

read_counts <- function(path) {
  dt <- fread(path, data.table=FALSE)
  nms <- names(dt); lname <- tolower(nms)
  if (!"gene" %in% nms) {
    sym <- which(lname %in% c("symbol","gene","gene_symbol","genesymbol","hgnc_symbol","gene.name","genename"))
    ens <- which(lname %in% c("ensembl","ensembl_id","id","geneid","gene_id"))
    if (length(sym)) dt$gene <- dt[[ nms[sym[1]] ]]
    else if (length(ens)) dt$gene <- dt[[ nms[ens[1]] ]]
  }
  stopifnot("gene" %in% names(dt))
  dt
}

# Infer group from the prefix before the first underscore (PBS_*, Lo_*, Med_*, Hi_*)
infer_group <- function(vec) {
  pref <- toupper(sub("_.*$", "", vec))
  dplyr::recode(pref,
                "PBS" = "PBS",
                "LO"  = "Lo",
                "MED" = "Med",
                "HI"  = "Hi",
                .default = NA_character_)
}

make_contrasts <- function(levels, wanted) {
  vapply(wanted, function(w){
    parts <- strsplit(w, "_vs_")[[1]]
    stopifnot(length(parts)==2)
    sprintf("%s - %s", parts[1], parts[2])
  }, character(1))
}

cfg <- load_cfg()
counts_path <- cfg$paths$counts %||% "data/raw_annotated_combined.counts"
outdir <- safe_dir(cfg$paths$results %||% "outputs/results")
treat_levels_cfg <- cfg$factors$treatment_levels %||% c("PBS","Lo","Med","Hi")
contrasts_cfg <- cfg$params$contrasts
cpm_min <- cfg$params$cpm_min %||% 1
cpm_min_samples <- cfg$params$cpm_min_samples %||% 3

message("Reading counts: ", counts_path)
raw <- read_counts(counts_path)
sample_cols <- setdiff(names(raw),
  c("gene","Gene","SYMBOL","symbol","GeneID","ENTREZID","id","ID","Ensembl","ensembl","EnsemblID","ensembl_id"))
stopifnot(length(sample_cols) >= 2)

# Build group factor from sample name prefixes
sample_groups <- infer_group(sample_cols)

# Print mapping preview
mapping_df <- tibble(sample = sample_cols,
                     prefix = sub("_.*$", "", sample_cols),
                     group  = sample_groups)
print(mapping_df, n = nrow(mapping_df))
tab <- table(sample_groups, useNA="ifany"); print(tab)

# Ensure groups in levels and order as config
group <- factor(sample_groups, levels = treat_levels_cfg)
if (any(is.na(group))) {
  bad <- mapping_df %>% filter(is.na(group)) %>% pull(sample)
  stop("Some samples could not be assigned to any of ",
       paste(treat_levels_cfg, collapse=", "),
       ". Offending samples: ", paste(bad, collapse=", "))
}

# DGE setup
y <- DGEList(counts = as.matrix(raw[, sample_cols, drop=FALSE]), genes = raw["gene"])
keep <- filterByExpr(y, group = group, min.count = cpm_min, min.total.count = 10)
keep <- keep | rowSums(cpm(y) >= cpm_min) >= cpm_min_samples
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y, method="TMM")

# Design and voom
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
v <- voom(y, design, plot=FALSE)
fit <- lmFit(v, design)

# Contrasts
contr_formulas <- make_contrasts(levels = levels(group), wanted = contrasts_cfg)
contr_mat <- makeContrasts(contrasts = contr_formulas, levels = design)
fit2 <- contrasts.fit(fit, contr_mat)
fit2 <- eBayes(fit2)

# Write all requested contrasts (full tables)
for (i in seq_along(contrasts_cfg)) {
  cname <- contrasts_cfg[i]
  tt <- topTable(fit2, coef = i, n = Inf, sort.by="P")
  if (!"gene" %in% names(tt) && "ID" %in% names(tt)) tt$gene <- tt$ID
  if (!"gene" %in% names(tt)) tt$gene <- y$genes$gene[match(rownames(tt), rownames(y$counts))]
  cols <- intersect(c("gene","logFC","AveExpr","t","P.Value","adj.P.Val","B"), names(tt))
  tt <- tt[, unique(c(cols, setdiff(names(tt), cols)))]
  out <- file.path(outdir, sprintf("DE_%s.tsv", cname))
  fwrite(tt, out, sep="\t", quote=FALSE)
  message("Wrote: ", out, " (", nrow(tt), " rows)")
}

message("DE rebuild complete. Results in: ", outdir)
