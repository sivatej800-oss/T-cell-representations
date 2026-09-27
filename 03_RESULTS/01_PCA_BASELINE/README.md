# 01 — PCA baseline (the paper reproduction)

**Question.** Does the clustering workflow behind Figure 1D of Piehl et al. reproduce
from the public GSE200164 matrix?

**Result.** All 70,391 cells in 12 unsupervised clusters (the paper reports 13).
Clusters run from 18,533 cells down to 787. Requesting the full ~2,816-gene
integration pool instead of 2,000 is what keeps the rare populations. The earlier
2,000-feature run bottomed out at 671 cells.

**Method.** SCTransform per sample (v1, 1,000 variable genes), anchor integration,
PCs 1–15, Louvain at resolution 0.3.

**Code.** [`02_CODE/03_reproduce_figure1d.py`](../../02_CODE/03_reproduce_figure1d.py)
→ `03_reproduce_figure1d.R` (after Stages 01–02).

**Canonical folder.** [`Age_trajectory/results/03_seurat_reproduction/`](../../Age_trajectory/results/03_seurat_reproduction/)

| linked file | what it is |
| --- | --- |
| `03_seurat_umap_clusters.png` | The 12-cluster UMAP. |
| `GSE200164_umap_clusters.csv.gz` | Cluster label and UMAP coordinates per cell. |
| `GSE200164_sample_metadata_mapping.csv` | GEO accessions reconciled to the authors' A1–H8 sample names. |
| `figure1d_environment_audit.json` | Package versions and seeds for this run. |

Not linked, in the canonical folder: `integrated_clustered_umap.rds` (final Seurat
object), `GSE200164_after_PCA.rds` and `GSE200164_after_qc_audit.h5ad` (resume
checkpoints), `sct_samples/` (59 per-sample SCTransform objects), `seurat_input/`
(the Python→R hand-off), and `GSE200164_umap_clusters_BEFORE_2000features.csv.gz`
(the superseded 2,000-feature run, kept as the counter-example).
