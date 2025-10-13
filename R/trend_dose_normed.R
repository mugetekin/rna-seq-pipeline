suppressPackageStartupMessages({ library(tidyverse); library(ggplot2) })
source("R/utils_normed.R")
out_dir <- "outputs_normed"; dir.create(out_dir, showWarnings = FALSE)

obj <- load_normed_expr("data/normed_cpms_filtered_annot.csv")
X   <- obj$expr_log; grp <- factor(obj$grp, levels=c("PBS","Lo","Med","Hi"))
samples <- colnames(X)

# pick top K variable genes
K <- 12L
top_idx <- order(apply(X,1,var), decreasing = TRUE)[1:min(K, nrow(X))]
mat     <- X[top_idx, , drop=FALSE]

df <- as.data.frame(mat) |>
  tibble::rownames_to_column("Gene") |>
  tidyr::pivot_longer(-Gene, names_to="Sample", values_to="logCPM") |>
  dplyr::mutate(Group = grp[match(Sample, samples)])

sem <- function(x) sd(x)/sqrt(sum(is.finite(x)))
sumdf <- df |>
  dplyr::group_by(Gene, Group) |>
  dplyr::summarise(mean = mean(logCPM, na.rm=TRUE), sem = sem(logCPM), .groups="drop")

p <- ggplot(sumdf, aes(Group, mean, group = 1)) +
  geom_line() + geom_point() +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.15) +
  facet_wrap(~ Gene, ncol = 4, scales = "free_y") +
  labs(title = "Top variable genes — mean ± SEM across dose groups",
       y = "log2(CPM+1)", x = "") +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "Trend_TopVar_MeanSEM.png"), p, width = 14, height = 9, dpi = 300)
message("Saved: ", file.path(out_dir, "Trend_TopVar_MeanSEM.png"))
