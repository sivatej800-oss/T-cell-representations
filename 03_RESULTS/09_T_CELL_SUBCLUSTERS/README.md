# 09 — T-cell subclusters

> **PLANNED / NOT YET COMPLETED.** No code has been written and no results exist.

**Goal.** Subset the T-cell clusters from [08_CELL_TYPE_ANNOTATION](../08_CELL_TYPE_ANNOTATION/)
and recluster them, following the authors: PCA, neighbours, clustering, and UMAP within
the subset. This gives the Figure 1D labels CD4+, CD8+, CD4+/CD8+ and regulatory T
cells, for each arm.

**Watch for.** CD4 is not a T-cell marker in the full dataset. It is expressed in
68–88% of myeloid-cluster cells, so it only becomes usable after subsetting.
