# 04 — Silhouette at native widths

**Question.** Which partition is more compact and better separated, with each arm at
its native width (PCA-15 against scGPT-512)?

**Result.** Each partition is scored in both spaces, and comparisons only run down a
column. Comparing across a row is invalid because distances in 15 and 512 dimensions
are not the same quantity.

| partition | scored in PCA-15 | scored in scGPT-512 |
| --- | --- | --- |
| PCA labels (12) | **0.176** | 0.100 |
| scGPT labels (11) | 0.151 | **0.220** |

Each partition scores highest in its own space, so neither description wins in both.
That says each favours the clusters it produced, not that either is better defined.

**Code.** [`02_CODE/07_silhouette_comparison.py`](../../02_CODE/07_silhouette_comparison.py) → `.R`.

**Canonical folder.** [`Age_trajectory/results/07_silhouette/`](../../Age_trajectory/results/07_silhouette/)

| linked file | what it is |
| --- | --- |
| `07_silhouette_per_cluster.png` | Per-cluster silhouette, both partitions in both spaces. |
| `silhouette_summary.csv` | Mean silhouette and negative fraction, per partition per space. |
| `silhouette_per_cluster.csv` | The same by cluster, with cluster sizes. |

Not linked: `silhouette_environment_audit.json`.
