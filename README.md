# RNA-seq Analysis Pipeline (mouse pituitary, PBS/Lo/Med/Hi design)

This repository provides a stepwise, script-oriented RNA-seq pipeline for differential expression, gene set enrichment, quality control, and “gene-of-interest” (GOI) visualization. The code is written in R and organized so each stage can be run independently or end-to-end.

The pipeline assumes a PBS/Lo/Med/Hi treatment design (PBS as reference) and a single raw counts file with gene annotations. All outputs are written to versioned subfolders and large artifacts are tracked via Git LFS.

Contents

├── R/

│   ├── io_helpers.R         # I/O, config, annotation parsing, output helpers

│   └── de_helpers.R         # normalization, design/contrasts, model fitting

├── scripts/

│   ├── 01_normalize.R       # CPM filter, TMM, voom; post-normalization exports

│   ├── 02_de_limma.R        # limma-voom DE; contrast tables (Lo/Med/Hi vs PBS)

│   ├── 03_go_enrichment.R   # GSEA over GO (BP); dotplots and CSVs

│   ├── 04_plots_and_goi.R   # GOI ranking/shortlists; violin/heatmap/bar+SEM

│   ├── 05_qc_and_volcano.R  # QC (library size, MDS/PCA, MA) and volcano plots

│   ├── 06_go_networks.R     # GO term networks (emapplot, cnetplot)

│   └── 07_goi_zoom.R        # Focused analysis for selected GOIs (incl. Venn)

├── run_all.R                # Orchestrates scripts 01–07 with error trapping

├── Makefile                 # `make all` (01→07), or run steps individually

├── config.yaml              # Paths and analysis parameters

├── data/

│   └── raw_annotated_combined.counts   # raw counts+annotation

└── outputs/

│   ├── rds/                 # intermediate R objects

│   ├── results/             # TSV/CSV analysis tables

│  └── figures/             # PNG/PDF figures


## Requirements

R 4.3+ on Windows/macOS/Linux.

**R packages:** data.table, tidyverse, edgeR, limma, yaml, clusterProfiler, enrichplot, DOSE, org.Mm.eg.db, GO.db, ggrepel, factoextra, ragg, ggVennDiagram (the scripts check and install when feasible).

Git LFS for large binary outputs (already configured via .gitattributes).

**Windows users:** If Rscript is not on PATH, invoke it explicitly, e.g.
"/c/Program Files/R/R-4.3.3/bin/x64/Rscript.exe" --vanilla ...

## Input specification

**Counts file:** data/raw_annotated_combined.counts
A tab-delimited matrix with at least one annotation column (e.g. id, SYMBOL, gene_type) and then sample columns (counts).

Row identifiers are built as SYMBOL (fallback to id) with make.unique.

**Sample naming:** Group is inferred from sample name prefix:

PBS* → PBS (reference) 

Lo* → Lo

Med* → Med

Hi* → Hi

Samples not matching these prefixes are labeled “Other” and safely excluded from contrasts.

## Configuration

**Edit config.yaml:**

**paths:**

  counts:  "data/raw_annotated_combined.counts"

  outputs: "outputs"

  rds:     "outputs/rds"

  results: "outputs/results"

  figures: "outputs/figures"


**params:**

  cpm_min: 1

  cpm_min_samples: 3

  ref_group: "PBS"

  contrasts: ["Lo_vs_PBS","Med_vs_PBS","Hi_vs_PBS"]

  lfc_thresh: 1

  fdr_thresh: 0.05

  go_ont: "BP"

  go_n_show: 10


## Quick start

**Option A — one command (recommended):**

Rscript --vanilla run_all.R

## on Windows Git Bash, if Rscript is not on PATH:
"/c/Program Files/R/R-4.3.3/bin/x64/Rscript.exe" --vanilla run_all.R


**Option B — with Make:**

make all      # runs 01→07

## or individual targets:
make normalize de go goi qc networks goi_zoom


All outputs are written under outputs/ and intermediate objects under outputs/rds/.

What each step does

**01_normalize.R:** Filters low counts (CPM ≥ cpm_min in ≥ cpm_min_samples; plus filterByExpr). TMM normalization and voom. Exports: filtered raw counts, CPM (TMM), voom logCPM.

**02_de_limma.R:** Builds design (PBS/Lo/Med/Hi), computes contrasts vs PBS. Fits lmFit + eBayes; writes full DE tables (.tsv) for each contrast.

**03_go_enrichment.R:** Ranks genes by logFC, maps IDs to Entrez, performs GSEA (GO:BP). Saves dotplots and CSVs per contrast; robust to small gene lists.

**04_plots_and_goi.R:** Data-driven GOI ranking combining DE and GSEA signals. Produces: GOI violin plots, heatmap (row-Z), bar+jitter with mean±SEM and Dunnett stars vs PBS. Writes GOI_ranked_all.csv and two shortlists (broad/strict).

**05_qc_and_volcano.R:** QC: library sizes, MDS, PCA, MA plots; Volcano plots per contrast (with optional labels). Dynamic handling of available contrast names.

**06_go_networks.R:** Term–term similarity networks (emapplot) and gene–term networks (cnetplot). Uses fold-change colors where available; saves PNG/PDF.

**07_goi_zoom.R:** Focused analysis for specified GOIs (default: Bub1, Cenpi, Esco2, Rpl36-ps12). Plots violin / bar+jitter for these genes; Venn diagrams for leading-edge genes (global) and for Lo vs PBS; Exports per-gene summary tables.

## Outputs (selected)

outputs/results/DE_<Lo|Med|Hi>_vs_PBS.tsv — full limma DE results.

outputs/results/go/gseGO_<Lo|Med|Hi>.csv — GSEA results (GO:BP).

outputs/results/goi/GOI_ranked_all.csv — global GOI ranking.

outputs/figures/ — QC (library size, MDS/PCA/MA), volcanoes, GO dotplots, GO networks, GOI panels, and GOI zoom (incl. Venn).

outputs/rds/ — intermediate objects (norm_and_fit.rds, de_tables.rds, go_objs.rds) for reuse.

All images and large artifacts are tracked by Git LFS.

## Reproducibility notes

**Normalization:** CPM filtering + TMM; filterByExpr is applied to match limma/edgeR best practice.

**Modeling:** limma-voom with empirical Bayes moderation; PBS as reference.

**Multiple testing:** Benjamini–Hochberg FDR unless stated.

**GSEA:** GO:BP with clusterProfiler::gseGO; robust ID mapping via org.Mm.eg.db.

**GOI logic:** Integrates DE magnitude/significance and enrichment context; shortlists are parameterized.

## Troubleshooting

“Coefficients not estimable: Other”

Informational: samples not matching PBS/Lo/Med/Hi are labeled “Other” and excluded from contrasts. Safe to ignore.

Contrast name mismatch

Scripts detect available contrast names and will error with a clear message if a requested contrast is missing. Check colnames(obj$fit2$coefficients) via:

obj <- readRDS("outputs/rds/norm_and_fit.rds"); colnames(obj$fit2$coefficients)


### Windows line endings warnings (CRLF)

Benign; the repo is configured to normalize line endings. You can suppress the warning or set git config core.autocrlf input if preferred.

### Memory usage

The pipeline exports intermediate RDS/CSV snapshots to avoid re-computation and to keep per-session memory bounded.

### Extending the pipeline

Change GOI list: Edit the goi vector in scripts/07_goi_zoom.R.

Add contrasts: Update params.contrasts in config.yaml and re-run 02_de_limma.R.

Different ontology: Set params.go_ont to MF or CC and re-run 03_go_enrichment.R.

New design: If your study design deviates from PBS/Lo/Med/Hi, adapt make_meta_from_colnames() in R/io_helpers.R and the model matrix in R/de_helpers.R.

### How to cite

Please cite the software/libraries that enable this workflow:

edgeR — Robinson MD, McCarthy DJ, Smyth GK. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics (2010).

limma/voom — Law CW, Chen Y, Shi W, Smyth GK. voom: Precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biol (2014); Ritchie ME et al., Nucleic Acids Res (2015).

clusterProfiler / enrichplot / DOSE — Yu G et al. and related packages for enrichment analysis and visualization.

org.Mm.eg.db / GO.db — Bioconductor annotation packages for Mus musculus and Gene Ontology.

If you analyze the published LPS pituitary dataset, also cite the corresponding data source and paper, and include the Dryad/SRA accession as appropriate.
