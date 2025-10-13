suppressPackageStartupMessages({ library(tidyverse); library(ggplot2) })
source("R/utils_normed.R")
out_dir <- "outputs_normed"; dir.create(out_dir, showWarnings = FALSE)

# --- REQUIRE Top20 list ---
top20_path <- file.path(out_dir, "Top20_variable_genes.csv")
if (!file.exists(top20_path)) {
  stop("Top20 list not found at: ", top20_path,
       "\nGenerate it first (e.g., run R/heatmap_top20_normed.R).")
}

# load data + groups
obj <- load_normed_expr("data/normed_cpms_filtered_annot.csv")
X   <- obj$expr_log
grp <- factor(obj$grp, levels = c("PBS","Lo","Med","Hi"))
samples <- colnames(X)

# read Top20 and preserve order
top20 <- readr::read_csv(top20_path, show_col_types = FALSE)$Gene
top20 <- unique(top20)  # safety
# subset to available genes (warn if any missing)
missing <- setdiff(top20, rownames(X))
if (length(missing)) message("WARNING: missing genes in matrix (will skip): ",
                             paste(missing, collapse = ", "))
genes <- top20[top20 %in% rownames(X)]
stopifnot(length(genes) > 0)

# long-form; keep Top20 order in facets
df <- as.data.frame(X[genes, , drop = FALSE]) |>
  tibble::rownames_to_column("Gene") |>
  tidyr::pivot_longer(-Gene, names_to = "Sample", values_to = "logCPM") |>
  dplyr::mutate(
    Group = grp[match(Sample, samples)],
    Gene  = factor(Gene, levels = genes)  # facet order = Top20 order
  )

# mean ± SEM
sem <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x)))
sumdf <- df |>
  dplyr::group_by(Gene, Group) |>
  dplyr::summarise(
    mean = mean(logCPM, na.rm = TRUE),
    sem  = sem(logCPM),
    .groups = "drop"
  )

# plot
p <- ggplot(sumdf, aes(Group, mean, group = Gene)) +
  geom_line(aes(color = Gene)) +
  geom_point(aes(color = Gene)) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.15) +
  facet_wrap(~ Gene, ncol = 4, scales = "free_y") +
  labs(title = "Top 20 genes — Mean ± SEM across dose groups",
       y = "log2(CPM+1)", x = "") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

ggsave(file.path(out_dir, "Trend_Genes_MeanSEM_Top20.png"), p,
       width = 14, height = 9, dpi = 300)
message("Saved: ", file.path(out_dir, "Trend_Genes_MeanSEM_Top20.png"))
