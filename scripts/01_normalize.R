# Load helper modules (I/O + DE utilities)
source("R/io_helpers.R"); source("R/de_helpers.R")
suppressPackageStartupMessages({ library(edgeR); library(limma) })

# Read config and ensure outputs/rds/results/figures exist
cfg  <- load_cfg(); dir_prep(cfg)

# Returns:
#   io$annot  : tibble(id, SYMBOL, gene_type)
#   io$counts : matrix-like counts (genes × samples), rownames are unique gene labels
io   <- read_counts_annot(cfg$paths$counts)
# Build metadata from column names (Group inferred from prefixes PBS/Lo/Med/Hi)
meta <- make_meta_from_colnames(colnames(io$counts))
# Sanity check: metadata SampleID matches sanitized column names
stopifnot(all(meta$SampleID == make.names(colnames(io$counts))))

# do_normalize now only takes (counts, meta); filterByExpr is handled inside
dge  <- do_normalize(io$counts, meta)

# fit_limma(dge, meta): voom -> design (~0+Group) -> contrasts -> eBayes
fit  <- fit_limma(dge, meta)
v    <- fit$v; fit2 <- fit$fit2

# Bundle key objects for downstream scripts (ANOVA, enrichment, etc.)
saveRDS(list(annot = io$annot, meta = meta, dge = dge, v = v, fit2 = fit2),
        file.path(cfg$paths$rds, "norm_and_fit.rds"))

# Export CPM, voom logCPM, and filtered counts with annotations
write_post_norm(dge, v, io$annot, cfg)
