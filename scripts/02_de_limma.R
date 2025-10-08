source("R/io_helpers.R")
suppressPackageStartupMessages({ library(limma); library(tidyverse) })

cfg <- load_cfg()

# read voom+fit from step 01
obj <- readRDS(file.path(cfg$paths$rds, "norm_and_fit.rds"))
fit2 <- obj$fit2

# DE tables
tt_lo  <- topTable(fit2, coef = "Lo_vs_PBS",  number = Inf, sort.by = "P")
tt_med <- topTable(fit2, coef = "Med_vs_PBS", number = Inf, sort.by = "P")
tt_hi  <- topTable(fit2, coef = "Hi_vs_PBS",  number = Inf, sort.by = "P")

# write results
dir.create(cfg$paths$results, showWarnings = FALSE, recursive = TRUE)
write.table(tt_lo,  file.path(cfg$paths$results, "DE_Lo_vs_PBS.tsv"),  sep="\t", quote=FALSE)
write.table(tt_med, file.path(cfg$paths$results, "DE_Med_vs_PBS.tsv"), sep="\t", quote=FALSE)
write.table(tt_hi,  file.path(cfg$paths$results, "DE_Hi_vs_PBS.tsv"),  sep="\t", quote=FALSE)

# save RDS for downstream steps
saveRDS(list(tt_lo=tt_lo, tt_med=tt_med, tt_hi=tt_hi),
        file.path(cfg$paths$rds, "de_tables.rds"))
