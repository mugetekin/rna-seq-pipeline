Stepwise RNA-seq pipeline split into scripts:
- `scripts/01_normalize.R`: filtering (CPM), TMM, voom; snapshot exports
- `scripts/02_de_limma.R`: contrasts + limma DE tables
- `scripts/03_go_enrichment.R`: GSEA (gseGO) + figures
- `scripts/04_plots_and_goi.R`: GOI shortlist + TGFBR & generic plots

Edit `config.yaml` to point to your raw counts and tweak params.
