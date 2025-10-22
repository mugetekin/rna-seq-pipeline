# GOI Zoom Figures — RNA-seq Analysis
This folder provides detailed visualizations ('zoomed-in') for selected **Genes of Interest (GOI)** and their overlaps across contrasts. Selected genes are Bub1, Cenpi, Esco2, Rpl36-ps12.

---
### Goi Bar Jitter Focus

![Goi Bar Jitter Focus](GOI_bar_jitter_focus.png)
**Bar plot with jittered points** Mean ± SEM log₂(CPM + 1) expression is shown for Bub1, Cenpi, Esco2, and Rpl36-ps12 with individual sample points (jittered). These genes, associated with cell-cycle progression and ribosomal function, display higher expression in LPS-treated groups compared to PBS, consistent with the activation of mitotic and translational pathways observed in the GO:BP enrichment and broader GOI analyses.

---
### Goi Correlation Heatmap

![Goi Correlation Heatmap](GOI_correlation_heatmap.png)
**Correlation heatmap** Pairwise Pearson correlations show strong positive co-expression among Bub1, Cenpi, and Esco2, consistent with their shared roles in mitotic regulation. In contrast, Rpl36-ps12 exhibits weak or inverse correlation with these genes, suggesting its expression follows a distinct translational pattern. These relationships parallel the coordinated and opposing trends observed in the GOI bar and violin plots.

---
### Goi Lollipop Logfc

![Goi Lollipop Logfc](GOI_lollipop_logFC.png)
**Lollipop plot** Estimated expression changes (±95% CI) for Bub1, Cenpi, Esco2, and Rpl36-ps12 are shown relative to PBS. Bub1, Cenpi, and Esco2 show consistent up-regulation across Lo–Hi contrasts, reflecting activation of mitotic and cell-cycle pathways under LPS treatment. In contrast, Rpl36-ps12 displays a weaker or inverse trend, likely due to its pseudogene origin and partial decoupling from active ribosomal transcription, resulting in less coordinated regulation with the core cell-cycle genes.

---
### Up Overlaps Bar

![Up Overlaps Bar](UP_overlaps_bar.png)
**Bar plot with jittered points** Bars represent the number of genes with FDR < 0.05 and |log₂FC| > 1 that are significantly upregulated in each condition. A large set of genes (n = 396) is unique to the Lo group, while only a small subset (n = 11) remains consistently upregulated across Lo and Hi. This indicates a strong, early transcriptional activation at low-dose LPS that narrows to a more selective response under high-dose stimulation.

---
### Venn Goi Vs Lo Down

![Venn Goi Vs Lo Down](Venn_GOI_vs_Lo_DOWN.png)
**Venn diagram** showing overlaps between gene sets (e.g., DEGs or GOI subsets) across indicated contrasts. Shared regions represent genes common to multiple groups.

---
### Venn Goi Vs Lo Up

![Venn Goi Vs Lo Up](Venn_GOI_vs_Lo_UP.png)
**Venn diagram** showing overlaps between gene sets (e.g., DEGs or GOI subsets) across indicated contrasts. Shared regions represent genes common to multiple groups.

---
### Venn Up Lo Med Hi

![Venn Up Lo Med Hi](Venn_UP_Lo_Med_Hi.png)
**Venn diagram** showing overlaps between gene sets (e.g., DEGs or GOI subsets) across indicated contrasts. Shared regions represent genes common to multiple groups.

---
