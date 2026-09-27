# 02 — scGPT baseline (the frozen-embedding arm)

**Question.** What does the same downstream workflow produce when the cell
representation is a frozen scGPT embedding instead of PCA?

**Result.** 13 clusters, with no batch correction of any kind. It splits finer at the
small end: five clusters under 1,000 cells, down to 103.

**Gene coverage.** 15,960 of 19,689 genes (81.1%) are in scGPT's vocabulary. Those
genes carry **98.86% of total expression**, and none of 65 canonical immune markers
is missing. The dropped genes are mostly clone-based loci, lncRNAs and renamed
symbols, concentrated in histones and tRNA synthetases.

**Code.** [`02_CODE/04_scgpt_figure1d.py`](../../02_CODE/04_scgpt_figure1d.py) → `04_scgpt_figure1d.R`;
coverage audit by [`02_CODE/05_audit_scgpt_gene_coverage.py`](../../02_CODE/05_audit_scgpt_gene_coverage.py).

**Canonical folder.** [`Age_trajectory/results/04_scgpt_embedding/`](../../Age_trajectory/results/04_scgpt_embedding/)

| linked file | what it is |
| --- | --- |
| `04_scgpt_raw_umap_clusters.png` | The 13-cluster UMAP, same size and palette as the PCA baseline. |
| `GSE200164_scgpt_umap_clusters.csv.gz` | Cluster label and UMAP coordinates per cell. |
| `scgpt_gene_coverage_audit.json` | Which genes matched the vocabulary, and the expression behind the gap. |
| `scgpt_environment_audit.json` | Package versions and seeds for this run. |

Not linked: `scgpt_clustered_umap.rds` (final Seurat object) and `r_input/` (the
embedding handed to R).
