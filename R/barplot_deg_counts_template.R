suppressPackageStartupMessages({ library(tidyverse); library(ggplot2) })
out_dir <- "outputs_normed"; dir.create(out_dir, showWarnings = FALSE)

# === EDIT ME LATER ===
de_csv   <- file.path(out_dir, "DE_results_example.csv")
alpha    <- 0.05
# =====================

stopifnot(file.exists(de_csv))
de <- readr::read_csv(de_csv, show_col_types = FALSE) |>
  dplyr::mutate(Direction = dplyr::case_when(
    padj < alpha & log2FC > 0 ~ "Up",
    padj < alpha & log2FC < 0 ~ "Down",
    TRUE ~ "NS"
  )) |>
  dplyr::filter(Direction != "NS")

counts <- de |>
  dplyr::count(Contrast, Direction) |>
  tidyr::complete(Contrast, Direction, fill = list(n = 0))

p <- ggplot(counts, aes(x = Contrast, y = n, fill = Direction)) +
  geom_col(position = "stack", width = 0.65) +
  coord_flip() +
  labs(title = "DEG counts per contrast (padj < 0.05)", x = "", y = "Count") +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "DE_BarCounts_byContrast.png"), p, width = 8, height = 5, dpi = 300)
message("Saved: ", file.path(out_dir, "DE_BarCounts_byContrast.png"))
