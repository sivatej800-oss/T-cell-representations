# 06 — Dimension control (PCA-512 against scGPT-512)

**Question.** Is scGPT's silhouette advantage just a consequence of having 512
dimensions instead of 30?

**Result: no.** Given the embedding's exact width, PCA gets worse rather than better:

| partition | scored in PCA-512 | scored in scGPT-512 |
| --- | --- | --- |
| PCA-512 labels (9 clusters) | **−0.033** | 0.160 |
| scGPT labels (11 clusters) | −0.048 | **0.219** |

PCA-512 is also a worse partition in its own right. It has 9 clusters with nothing
below 378 cells, so the rare populations are lost. Five of its nine clusters score
negative. All eleven scGPT clusters score positive, and the smallest ones score
highest. PCs 31–512 hold 59.7% of the variance, spread across many near-noise
directions, and that dilutes the signal that separates rare populations. Because
PCA-512 has the coarser partition, the cluster-count bias favours PCA, and PCA still
loses.

**Note.** The PCA-512 arm needs exact SVD (`approx = FALSE`). With the default method
the decomposition fails to converge, and nothing downstream reports it.

**Code.** [`02_CODE/10_pca512_arm.py`](../../02_CODE/10_pca512_arm.py) → `.R`, then
[`02_CODE/11_silhouette_dimension_matched.py`](../../02_CODE/11_silhouette_dimension_matched.py) → `.R`.

**Canonical folders.** [`Age_trajectory/results/06_harmony_benchmark/`](../../Age_trajectory/results/06_harmony_benchmark/) (the `pca512` files) and
[`Age_trajectory/results/11_silhouette_matched/`](../../Age_trajectory/results/11_silhouette_matched/)

| linked file | what it is |
| --- | --- |
| `10_pca512_harmony_umap.png` | PCA-512 UMAP, 9 clusters. |
| `11_silhouette_matched.png` | The matched 2×2. |
| `11_silhouette_matched_per_cluster.png` | Per-cluster width for both partitions in their own spaces. |
| `benchmark_pca512_summary.csv` | PCA-512 cluster counts and the variance split across components. |
| `silhouette_matched_summary.csv` | Mean silhouette and negative fraction, per partition per space. |
| `silhouette_matched_per_cluster.csv` | The same by cluster, with sizes. |

Not linked: `benchmark_pca512_arm.rds`, `benchmark_pca512_assignments.csv.gz`,
`pca512_environment_audit.json`, `silhouette_matched_environment_audit.json`.
