# 05 — Elbow / WCSS sweep *(supporting analysis)*

**Question.** Resolution 0.3 was inherited from the authors, not tuned here. Does it
sit at a sensible place on each representation's compactness curve?

**Result.** Resolution 0.3 gives 12 clusters for PCA-Harmony (15 dims) and 11 for
scGPT-Harmony (512 dims). Both curves keep descending past that point: PCA from 0.336
at 12 clusters to 0.278 at 20, scGPT from 0.443 at 11 to 0.382 at 21. By 13 clusters
PCA has covered 69% of its total drop across the sweep and scGPT 62%. No bend was
computed; the figure shows the shape.

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
