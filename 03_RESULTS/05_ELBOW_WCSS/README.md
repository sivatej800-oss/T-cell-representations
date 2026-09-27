# 05 — Elbow / WCSS sweep *(supporting analysis)*

**Question.** Resolution 0.3 was inherited from the authors, not tuned here. Does it
sit at a sensible place on each representation's compactness curve?

**Result.** For PCA-Harmony (30 dims) it lands on the bend, at 12–13 clusters. For
scGPT-Harmony (512 dims) it lands *before* the bend, at 11 clusters, while the curve
is still descending. scGPT keeps improving out to about 16–18 clusters.

**Read the shapes, not the heights.** The two curves live in spaces of different
dimensionality. The PCA curve sitting lower is not evidence that PCA clusters better.
This analysis supports the verdict; it does not carry it.

**Code.** [`02_CODE/08_wcss_resolution_sweep.py`](../../02_CODE/08_wcss_resolution_sweep.py) → `.R`.

**Canonical folder.** [`Age_trajectory/results/08_wcss_sweep/`](../../Age_trajectory/results/08_wcss_sweep/)

| linked file | what it is |
| --- | --- |
| `08_wcss_compactness.png` | Normalized WCSS against cluster count — the elbow. |
| `08_wcss_dispersion_explained.png` | The same sweep as dispersion explained. |
| `wcss_resolution_sweep.csv` | Resolution, cluster count, normalized WCSS, per representation. |

Not linked: `wcss_sweep_audit.json`.
