# 08 — Cell-type annotation (first pass, broad categories)

**Question.** Which broad cell type is each cluster, in each arm?

**Method.** By hand, the way the authors did it. Each cluster was read off its evidence
row and feature plots and assigned one of the authors' seven categories. Cluster
numbers do not carry across arms.

**Result.** Annotated on the rebuilt PCA arm (SCT residuals, 12 clusters) and the current
`scgpt_harmony` labels; the earlier annotation is superseded. The two arms agree:

| identity | PCA-15 Harmony | scGPT Harmony | difference |
| --- | --- | --- | --- |
| T Cells | 54,803 | 54,684 | 119 |
| Monocytes | 7,387 | 7,470 | 83 |
| DC | 4,794 | 4,849 | 55 |
| NK Cells | 2,369 | 2,305 | 64 |
| B Cells | 217 | 208 | 9 |
| Undetermined | 821 | 875 | 54 |

PCA cluster 3 (macrophage-like: APOE, C1Q genes, LYVE1) is filed under Monocytes. PCA
cluster 8 and scGPT cluster 7 are the same cells, dominated by XIST and long-transcript
genes with no canonical lineage marker, and are Undetermined. PCA's one extra cluster
subdivides the T cells (cycling cells). `Plasma Cells` went unused in both arms.

**This is the first of two passes.** Figure 1D's finer labels (CD4/CD8/Treg and the
monocyte subtypes) come from reclustering the T-cell and monocyte subsets. That work is
[09](../09_T_CELL_SUBCLUSTERS/) and [10](../10_MONOCYTE_SUBCLUSTERS/), and neither has
been done yet.

**Code.** [`02_CODE/15_annotate_clusters.R`](../../02_CODE/15_annotate_clusters.R).

**Canonical folder.** [`Age_trajectory/results/15_annotation/`](../../Age_trajectory/results/15_annotation/)

| linked file | what it is |
| --- | --- |
| `15_celltypes_pca15_harmony_cluster.png` | Broad cell types over the PCA-15 layout. |
| `15_celltypes_scgpt_harmony_cluster.png` | The same over the scGPT layout. |
| `annotated_cluster_labels.csv.gz` | Cell-type identity per cell, both arms. |
| `annotation_rationale.csv` | The evidence behind each hand assignment. |

Not linked: `annotation_audit.json`.
