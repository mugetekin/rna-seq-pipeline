source("R/io_helpers.R"); source("R/de_helpers.R")
suppressPackageStartupMessages({ library(edgeR); library(limma) })

cfg  <- load_cfg(); dir_prep(cfg)

io   <- read_counts_annot(cfg$paths$counts)
meta <- make_meta_from_colnames(colnames(io$counts))
stopifnot(all(meta$SampleID == make.names(colnames(io$counts))))

# do_normalize artık sadece (counts, meta) alıyor; filterByExpr içeride
dge  <- do_normalize(io$counts, meta)

fit  <- fit_limma(dge, meta)
v    <- fit$v; fit2 <- fit$fit2

saveRDS(list(annot = io$annot, meta = meta, dge = dge, v = v, fit2 = fit2),
        file.path(cfg$paths$rds, "norm_and_fit.rds"))

write_post_norm(dge, v, io$annot, cfg)
