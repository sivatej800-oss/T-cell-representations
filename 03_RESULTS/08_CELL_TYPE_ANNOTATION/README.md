# 08 — Cell-type annotation (first pass, broad categories)

**Question.** Which broad cell type is each cluster, in each arm?

**Method.** By hand, the way the authors did it. Each cluster was read off its evidence
row and feature plots and assigned one of the authors' seven categories. Cluster
numbers do not carry across arms.

**Result.** The two arms agree at this resolution:

| identity | PCA-30 Harmony | scGPT Harmony | difference |
| --- | --- | --- | --- |
| T Cells | 55,485 | 55,554 | 69 |
| Monocytes | 7,392 | 7,471 | 79 |
| DC | 4,931 | 4,848 | 83 |
| NK Cells | 2,362 | 2,311 | 51 |
| B Cells | 221 | 207 | 14 |

PCA's two extra clusters subdivide the T cells (a Treg cluster and a cycling cluster);
they are not a lineage that scGPT missed. `Undetermined` and `Plasma Cells` went unused
in both arms.

**This is the first of two passes.** Figure 1D's finer labels (CD4/CD8/Treg and the
monocyte subtypes) come from reclustering the T-cell and monocyte subsets. That work is
[09](../09_T_CELL_SUBCLUSTERS/) and [10](../10_MONOCYTE_SUBCLUSTERS/), and neither has
been done yet.

**Code.** [`02_CODE/15_annotate_clusters.R`](../../02_CODE/15_annotate_clusters.R).

**Canonical folder.** [`Age_trajectory/results/15_annotation/`](../../Age_trajectory/results/15_annotation/)

| linked file | what it is |
| --- | --- |
| `15_celltypes_pca30_harmony_cluster.png` | Broad cell types over the PCA-30 layout. |
| `15_celltypes_scgpt_harmony_cluster.png` | The same over the scGPT layout. |
| `annotated_cluster_labels.csv.gz` | Cell-type identity per cell, both arms. |
| `annotation_rationale.csv` | The evidence behind each hand assignment. |

Not linked: `annotation_audit.json`.
