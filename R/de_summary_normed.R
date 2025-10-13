suppressPackageStartupMessages({ library(tidyverse); library(ggplot2) })
source("R/utils_normed.R")
out_dir <- "outputs_normed"; dir.create(out_dir, showWarnings = FALSE)

obj <- load_normed_expr("data/normed_cpms_filtered_annot.csv")
X   <- obj$expr_log; grp <- obj$grp; samples <- colnames(X)
stopifnot(all(c("PBS","Hi") %in% levels(grp)))

pbs <- which(grp=="PBS"); hi <- which(grp=="Hi")

# per-gene Welch t-test with NA-safety
safe_t <- function(x){
  a <- x[pbs]; b <- x[hi]
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (length(a) < 2 || length(b) < 2) return(NA_real_)
  t.test(a, b, alternative="two.sided", var.equal=FALSE)$p.value
}

# compute stats
mu_pbs <- apply(X[, pbs, drop=FALSE], 1, function(v) mean(v[is.finite(v)], na.rm=TRUE))
mu_hi  <- apply(X[, hi,  drop=FALSE],  1, function(v) mean(v[is.finite(v)], na.rm=TRUE))
log2FC <- mu_hi - mu_pbs
pvals  <- apply(X, 1, safe_t)
padj   <- p.adjust(pvals, method = "BH")

res <- tibble(Gene = rownames(X), log2FC = log2FC, pval = pvals, padj = padj)

# how many skipped
nskip <- sum(!is.finite(res$pval))
message("Skipped genes (insufficient data for t-test): ", nskip)

# Volcano (drop NA p-values)
res_plot <- res |> dplyr::filter(is.finite(padj))
p1 <- ggplot(res_plot, aes(x = log2FC, y = -log10(padj))) +
  geom_point(alpha = 0.6, size = 1.4) +
  geom_vline(xintercept = c(-1, 1), linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  labs(title = "Mock Volcano (Hi vs PBS)", x = "log2FC (Hi - PBS)", y = "-log10(adj p)") +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "DE_Volcano_mock_Hi_vs_PBS.png"), p1, width = 8, height = 6, dpi = 300)

# Bar of up/down counts (ignore NA padj)
up   <- sum(res$padj < 0.05 & res$log2FC > 0, na.rm = TRUE)
down <- sum(res$padj < 0.05 & res$log2FC < 0, na.rm = TRUE)
bar  <- tibble(Direction=c("Up","Down"), Count=c(up,down))
p2 <- ggplot(bar, aes(Direction, Count, fill = Direction)) +
  geom_col(width = 0.6) +
  theme_minimal(base_size = 12) + theme(legend.position = "none") +
  labs(title = "Mock DEG counts (Hi vs PBS)")
ggsave(file.path(out_dir, "DE_BarCounts_mock_Hi_vs_PBS.png"), p2, width = 6, height = 5, dpi = 300)

# write table (keep NAs so downstream knows which were skipped)
readr::write_csv(res, file.path(out_dir, "DE_Mock_Hi_vs_PBS_results.csv"))
message("Saved: DE_Mock_Hi_vs_PBS_results.csv + volcano + bar counts")
