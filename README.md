# RNA-seq Analysis Pipeline (Mouse Pituitary: PBS / Lo / Med / Hi Design)

A modular **R-based RNA-seq analysis pipeline** designed for stepwise or end-to-end processing of **pituitary transcriptomics data** from a **PBS/Lo/Med/Hi LPS treatment experiment** (Garcia *et al.*, *Front. Endocrinol.*, 2024).

This pipeline performs normalization, differential expression, enrichment analysis, QC, and visualization, including **gene-of-interest (GOI)** summaries and **GSEA-driven insights**.

---

## Project Overview

```
R/                     # Reusable functions (I/O, normalization, modeling)
scripts/               # Analysis steps (01–07)
├── 01_normalize.R       # Filtering, TMM normalization, voom
├── 02_de_limma.R        # Differential expression (Lo/Med/Hi vs PBS)
├── 03_go_enrichment.R   # GSEA (GO:BP)
├── 04_plots_and_goi.R   # GOI ranking & plots
├── 05_qc_and_volcano.R  # QC & volcano plots
├── 06_go_networks.R     # GO term networks
└── 07_goi_zoom.R        # Detailed GOI & Venn analyses
run_all.R               # Master script (runs all steps sequentially)
Makefile                # For automation (`make all`)
config.yaml             # User parameters and file paths
data/                   # Raw input data (counts + annotations)
outputs/                # Results, figures, RDS intermediates
```

---

## Requirements

- **R 4.3+** (Windows/macOS/Linux)
- **Git LFS** (large files)
- Core R packages:
  ```
  data.table, tidyverse, edgeR, limma, yaml,
  clusterProfiler, enrichplot, DOSE, org.Mm.eg.db,
  GO.db, ggrepel, factoextra, ragg, ggVennDiagram
  ```

> The scripts auto-install missing packages when possible.

**Windows note:**  
If `Rscript` is not in your PATH:
```bash
"/c/Program Files/R/R-4.3.3/bin/x64/Rscript.exe" --vanilla run_all.R
```

---

## Input Data

**Raw counts:** `data/raw_annotated_combined.counts`  
Tab-delimited file with:
- 1+ annotation columns (e.g., `id`, `SYMBOL`, `gene_type`)
- Followed by raw count columns per sample

**Sample naming convention:**
| Prefix | Group | Notes |
|--------|--------|-------|
| PBS*   | PBS (reference) | control group |
| Lo*    | Low dose LPS    | treatment 1 |
| Med*   | Medium dose LPS | treatment 2 |
| Hi*    | High dose LPS   | treatment 3 |

Samples outside these prefixes are labeled “Other” and excluded from contrasts.

**Normalized CPMs:** `data/normed_cpms_filtered_annot.csv` Used for visualization, QC checks, and validation of main output trends.

---

## Configuration

Edit `config.yaml` to adjust:
```yaml
paths:
  counts:  data/raw_annotated_combined.counts
  outputs: outputs

params:
  ref_group: PBS
  contrasts: [Lo_vs_PBS, Med_vs_PBS, Hi_vs_PBS]
  cpm_min: 1
  cpm_min_samples: 3
  lfc_thresh: 1
  fdr_thresh: 0.05
  go_ont: BP
  go_n_show: 10
```

---

## Quick Start

### Option 1 — One command
```bash
Rscript --vanilla run_all.R
```

### Option 2 — With Make
```bash
make all
```
Or run individual stages:
```bash
make normalize de go goi qc networks goi_zoom
```

All results are stored in `outputs/`, with intermediates in `outputs/rds/`.

---

## Main Outputs

| Output | Description |
|--------|--------------|
| `outputs/results/DE_*_vs_PBS.tsv` | DE tables (limma-voom) |
| `outputs/results/go/gseGO_*.csv` | GO:BP enrichment results |
| `outputs/results/goi/GOI_ranked_all.csv` | Ranked GOI summary |
| `outputs/figures/` | QC, volcanoes, dotplots, GO networks, GOI visuals |
| `outputs/rds/` | Serialized intermediate R objects |

---

## Pipeline Summary

| Step | Script | Description |
|------|---------|-------------|
| **01** | `normalize.R` | CPM filter, TMM, voom normalization |
| **02** | `de_limma.R` | DE modeling (PBS as reference) |
| **03** | `go_enrichment.R` | GSEA via `clusterProfiler` |
| **04** | `plots_and_goi.R` | GOI scoring + violin, heatmap, bar+SEM |
| **05** | `qc_and_volcano.R` | Library QC, PCA/MDS, volcano plots |
| **06** | `go_networks.R` | GO term & gene-term networks |
| **07** | `goi_zoom.R` | Focused per-GOI plots & Venn diagrams |

---

## Reproducibility Notes

- **Normalization:** TMM + `filterByExpr` (edgeR best practices)  
- **Modeling:** `limma-voom` with empirical Bayes  
- **Multiple testing:** Benjamini–Hochberg FDR  
- **GSEA:** `clusterProfiler::gseGO` (GO:BP ontology)  
- **GOI logic:** Combines DE significance and enrichment context  

---

## Troubleshooting

| Issue | Explanation |
|--------|-------------|
| “Coefficients not estimable: Other” | Samples not matching PBS/Lo/Med/Hi are excluded automatically. |
| Contrast name mismatch | Check contrast names via: `readRDS("outputs/rds/norm_and_fit.rds") |> \`colnames(obj$fit2$coefficients)` |
| Line-ending warnings (CRLF) | Safe to ignore on Windows. |
| Memory usage | Intermediate `.rds` files reduce memory load between steps. |

---

## Citation

If using this pipeline or its outputs, please cite:

**Primary dataset:**  
Garcia *et al.* (2024). *Lipopolysaccharide-induced chronic inflammation increases female serum gonadotropins and shifts the pituitary transcriptomic landscape.*  
*Front. Endocrinol.* 14:1279878. [DOI: 10.3389/fendo.2023.1279878](https://doi.org/10.3389/fendo.2023.1279878)

**Software & methods:**  
- **edgeR:** Robinson *et al.*, *Bioinformatics* (2010)  
- **limma/voom:** Law *et al.*, *Genome Biol* (2014); Ritchie *et al.*, *NAR* (2015)  
- **clusterProfiler/enrichplot/DOSE:** Yu *et al.*  
- **org.Mm.eg.db / GO.db:** Bioconductor annotation packages  
