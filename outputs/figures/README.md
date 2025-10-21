# Core Figures — Detailed, Plot-Specific Explanations

This README documents each plot in **`outputs/figures/`** with *plot-specific* notes.  
Where relevant, thresholds are pulled from `config.yaml` and exact rules are taken from `scripts/05_qc_and_volcano.R`.

**Key thresholds (from `config.yaml`):**
- `lfc_thresh = 1`
- `fdr_thresh = 0.05`

**Important nuance used by the script** (`scripts/05_qc_and_volcano.R`):
- **Coloring & significance** use **adjusted p-values (`adj.P.Val`)** **and** LFC: a point is called **Up** (red) if `adj.P.Val < 0.05` and `logFC > +1`; **Down** (blue) if `adj.P.Val < 0.05` and `logFC < -1`; otherwise **NS** (grey).
- The **y-axis** is **`-log10(P.Value)` (raw p-value)**. The dashed horizontal line is drawn at `-log10(fdr_thresh) ≈ 1.3` for reference, **even though coloring uses adjusted p-values**.  
  → This explains why, for some contrasts, many points sit above the dashed line but remain **grey**: their **raw p < 0.05**, but **adjusted p ≥ 0.05**.

---

## Volcano plots (by contrast)

### Volcano: Lo vs PBS &nbsp;(`Volcano_Lo_vs_PBS.png` and `_labeled`)
![Volcano Lo vs PBS Labeled](Volcano_Lo_vs_PBS_labeled.png)
- **What you see:** A dense red cloud on the right (positive log2FC), plus some blue on the left.  
- **Why:** Many genes satisfy both **|log2FC| ≥ 1** and **adj.P.Val < 0.05**, hence colored **Up/Down**.  
- **Labeled version:** The script labels the **top 12** DE genes among significant ones by **extreme logFC** (ties broken by effect size), split between up and down.  
- **Interpretation:** Lo vs PBS shows **broad upregulation** (numerically more Up than Down), suggesting a strong activation relative to PBS.

### Volcano: Med vs PBS &nbsp;(`Volcano_Med_vs_PBS.png` and `_labeled`)
- **What you see:** Fewer colored points than Lo vs PBS, but still a visible tail of Up genes.  
- **Why:** A moderate number pass **adj.P.Val < 0.05** *and* **|log2FC| ≥ 1**; others remain grey (raw p may be <0.05 but FDR is not).  
- **Labeled version:** Top 12 significant genes, prioritizing larger |log2FC|.

### Volcano: Hi vs PBS &nbsp;(`Volcano_Hi_vs_PBS.png` and `_labeled`)
![Volcano Hi vs PBS Labeled](Volcano_Hi_vs_PBS_labeled.png)
- **What you see:** **Most points are grey**, with only a **handful labeled** in the `_labeled` plot.  
- **Real reason (why it differs from Lo):**  
  - Coloring uses **adjusted** p-values; here, **few genes meet `adj.P.Val < 0.05` + |log2FC| ≥ 1**.  
  - You may still see many points above the dashed horizontal line (raw **p < 0.05**), but they remain grey because **FDR ≥ 0.05**.  
  - The labeled genes are the **top 12** that **do** pass the stricter adjusted threshold.
- **Interpretation:** Under Hi vs PBS, **effect sizes and/or variance** yield fewer FDR-significant genes; biology or sample variability can both contribute.

### Volcano: Lo vs Med, Lo vs Hi, Med vs Hi &nbsp;(`Volcano_*_vs_*`)
- **What you see:** Asymmetry varies by pair; typically fewer colored points than Lo vs PBS.  
- **Why:** Between-treatment comparisons often have smaller differential signal than treatment vs control, hence fewer points passing **adj.P.Val + logFC** thresholds.  
- **Labeled versions:** Again, **top 12** significant hits by extreme logFC are annotated.

**Legend note:** The legend “Down / NS / Up” corresponds to the **coloring rule using `adj.P.Val`**; the y-axis and dashed line relate to **raw p** only.

---

## MA plots (`MA_*`)
- **Axes:** M = log2 fold-change; A = average log-expression.  
- **Lo vs PBS:** Broad positive M tail (more **upregulated** genes at higher expression).  
- **Med vs PBS:** Intermediate pattern; fewer extremes than Lo vs PBS.  
- **Hi vs PBS:** Tail is thinner; fewer points exceed the |log2FC| cutoff at high confidence.  
- **Between treatments (Lo vs Med / Lo vs Hi / Med vs Hi):** Generally fewer strong shifts than treatment vs PBS.

---

## P-value histograms (`PvalHist_*`, `pval_histograms.png`)
![P-Value Histograms](pval_histograms.png)
- **Lo vs PBS:** Left-skewed (excess of small p-values), consistent with many DE genes.  
- **Med vs PBS:** Left skew present but milder.  
- **Hi vs PBS:** Closer to uniform, matching the scarcity of FDR-significant hits.  
- **Pairwise treatments:** Often closer to uniform than treatment vs PBS.

---

## PCA & MDS (`PCA_voomE.png`, `MDS_groups.png`)
- **What you see:** Samples cluster by group; 95% ellipses (PCA) capture within-group spread.  
- **Interpretation:** Tighter grouping implies consistent replicates; separation between groups suggests treatment effects.

---

## Library sizes (`library_sizes.png`)
- Bars per sample showing total reads; look for large imbalances pre-normalization. TMM in the pipeline accounts for this but extreme outliers are flagged for QC.

---

### How labeled points are chosen (applies to all `_labeled` volcanoes)
1. Compute **significance** using **adjusted p-values (`adj.P.Val`)** + **|log2FC| ≥ 1**.  
2. Split into **Up** and **Down**.  
3. Pick the **top 12** by **most extreme logFC** (separately for Up and Down), and annotate their gene symbols.

> Source: `scripts/05_qc_and_volcano.R` — functions `plot_volcano` and `plot_volcano_labeled` (color rule: `adj.P.Val < fdr & abs(logFC) > lfc`; y-axis uses `-log10(P.Value)`).

