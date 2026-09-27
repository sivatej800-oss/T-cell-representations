# 04 — Silhouette at native widths

**Question.** Which partition is more compact and better separated, with each arm at
its native width (PCA-30 against scGPT-512)?

**Result.** Each partition is scored in both spaces, and comparisons only run down a
column. Comparing across a row is invalid because distances in 30 and 512 dimensions
are not the same quantity.

| partition | scored in PCA-30 | scored in scGPT-512 |
| --- | --- | --- |
| PCA labels (13) | 0.006 | 0.103 |
| scGPT labels (11) | 0.008 | **0.207** |

In PCA space the two partitions are indistinguishable. In scGPT space they separate,
and every partition scores higher there. On its own, this pass cannot tell whether
the embedding organizes the cells better or whether silhouette just behaves
differently at 512 dimensions. [06_DIMENSION_CONTROL](../06_DIMENSION_CONTROL/) is
what settles that.

**Code.** [`02_CODE/07_silhouette_comparison.py`](../../02_CODE/07_silhouette_comparison.py) → `.R`.

**Canonical folder.** [`Age_trajectory/results/07_silhouette/`](../../Age_trajectory/results/07_silhouette/)

| linked file | what it is |
| --- | --- |
| `07_silhouette_per_cluster.png` | Per-cluster silhouette, both partitions in both spaces. |
| `silhouette_summary.csv` | Mean silhouette and negative fraction, per partition per space. |
| `silhouette_per_cluster.csv` | The same by cluster, with cluster sizes. |

Not linked: `silhouette_environment_audit.json`.
