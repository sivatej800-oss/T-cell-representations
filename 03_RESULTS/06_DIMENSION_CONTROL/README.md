# 06 — Dimension control (PCA-512 against scGPT-512)

**Question.** Is the silhouette comparison between PCA and scGPT just a consequence of
the arms having different numbers of dimensions (15 against 512)?

**Result.** At the embedding's exact width, each partition still scores best in its
own space:

| partition | scored in PCA-512 | scored in scGPT-512 |
| --- | --- | --- |
| PCA-512 labels (8 clusters) | **0.053** | 0.136 |
| scGPT labels (11 clusters) | −0.022 | **0.221** |
| *PCA-15 labels (12, carried)* | *0.029* | *0.100* |

The pattern matches native width, so dimensionality does not explain it. PCA-512 is a
different partition from PCA-15: 8 clusters with the smallest at 76 cells, against 12
and 217. Seven of its eight clusters score positive, and all eleven scGPT clusters do.
PCs 16–512 hold 64.7% of the variance, spread across many near-noise directions. PCA-512
has the coarser partition, so the cluster-count bias favours PCA, not scGPT.

**Note.** The PCA-512 arm needs exact SVD (`approx = FALSE`). With the default method
the decomposition fails to converge, and nothing downstream reports it.

**Code.** [`02_CODE/10_pca512_arm.py`](../../02_CODE/10_pca512_arm.py) → `.R`, then
[`02_CODE/11_silhouette_dimension_matched.py`](../../02_CODE/11_silhouette_dimension_matched.py) → `.R`.

**Canonical folders.** [`Age_trajectory/results/06_harmony_benchmark/`](../../Age_trajectory/results/06_harmony_benchmark/) (the `pca512` files) and
[`Age_trajectory/results/11_silhouette_matched/`](../../Age_trajectory/results/11_silhouette_matched/)

| linked file | what it is |
| --- | --- |
| `10_pca512_harmony_umap.png` | PCA-512 UMAP, 8 clusters. |
| `11_silhouette_matched.png` | The matched 2×2. |
| `11_silhouette_matched_per_cluster.png` | Per-cluster width for both partitions in their own spaces. |
| `benchmark_pca512_summary.csv` | PCA-512 cluster counts and the variance split across components. |
| `silhouette_matched_summary.csv` | Mean silhouette and negative fraction, per partition per space. |
| `silhouette_matched_per_cluster.csv` | The same by cluster, with sizes. |

Not linked: `benchmark_pca512_arm.rds`, `benchmark_pca512_assignments.csv.gz`,
`pca512_environment_audit.json`, `silhouette_matched_environment_audit.json`.
