# ANOVA (voom-normalized) — Figures and Explanations
This folder contains visualizations produced by `R/anova_normed.R`, which fits a limma model on voom-normalized logCPM and performs an **omnibus F-test across PBS/Lo/Med/Hi**.

## ANOVA_pval_hist.png
![ANOVA P-value Histogram](ANOVA_pval_hist.png)
- **What it shows:** Histogram of **raw** omnibus F-test p-values across all genes.
- **Reference line:** The red dashed line marks the FDR threshold from `config.yaml` (**0.05**) for orientation (note: coloring in results uses **adjusted p-values**, but this plot shows raw p-values).
- **Interpretation:** A **left-skew** (excess near 0) indicates many genes vary across groups; a near-uniform shape suggests few global differences.

## ANOVA_heatmap_top50_rowZ.png
![ANOVA Heatmap Top50](ANOVA_heatmap_top50_rowZ.png)
- **What it shows:** Heatmap of the **top 50 genes ranked by the ANOVA F-statistic**.
- **Scaling:** Rows are **Z-scored** (row-wise) so patterns reflect **relative up/down** per gene rather than absolute expression.
- **Annotations:** Columns annotated by **Group** (PBS/Lo/Med/Hi). Gene labels are **SYMBOLs** (ENSMUSG IDs mapped to symbols).
- **Interpretation:** Stripe-like blocks across columns indicate group-specific expression programs; tight blocks within a group suggest consistent replicates.

## ANOVA_violin_top12.png
![ANOVA Violin Plots Top12](ANOVA_violin_top12.png)
- **What it shows:** **Violin plots** for the **top 12 ANOVA genes** (highest F-statistics), showing **per-sample log2(CPM+1)** across groups.
- **Overlays:** Mean per group (white diamond) and **jittered** sample points (black outline) to show variability.
- **Selection logic:** Genes are ranked by **F-statistic** from the omnibus test; the top 12 by **F** are displayed.
- **Interpretation:** Shifts in violin location/width across groups indicate changes in mean expression and dispersion; large separation between PBS and treated groups suggests treatment-associated differences.

