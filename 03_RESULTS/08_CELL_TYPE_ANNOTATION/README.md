# 08 — Cell-type annotation (first pass, broad categories)

**Question.** Which broad cell type is each cluster, in each arm?

**Method.** By hand, the way the authors did it. Each cluster was read off its evidence
row and feature plots and assigned one of the authors' seven categories. Cluster
numbers do not carry across arms.

**Status: pending.** The earlier assignments were made on a PCA arm built from the SCT
`data` layer (13 clusters) and are superseded. The PCA arm has been rebuilt from SCT
residuals (12 clusters), and `scgpt_harmony` shifted slightly, so Stage 15 will be rerun
once the proposed assignments are approved. The earlier outputs are in
`Age_trajectory/results_superseded_2026-11-22/`. The links below stay dangling until then.

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
