# 03 — Batch correction (Harmony on both arms)

**Question.** The two baselines are not comparable: PCA was anchor-integrated, scGPT
was not corrected at all. With Harmony applied equally to both, how much patient
identity does each still carry?

**PCA input.** SCT residuals built the way the authors' `PrepSCTIntegration` builds
them, the closest match to their normalization that works without anchors. The arm
without Harmony is therefore *normalized per patient, no Harmony*, not raw: residuals
are centred within each patient. `percent.mt` is reduced, not removed.

**Result.**

| arm | dims | clusters | smallest | donor mixing | donor NMI |
| --- | --- | --- | --- | --- | --- |
| `pca_raw` (per patient, no Harmony) | 15 | 10 | 206 | 0.883 | 0.049 |
| `pca_harmony` | 15 | 12 | 217 | **0.969** | **0.038** |
| `pca_harmony_30dims` | 30 | 12 | 226 | 0.967 | 0.038 |
| `scgpt_raw` | 512 | 12 | 103 | 0.806 | 0.072 |
| `scgpt_harmony` | 512 | 11 | 105 | 0.866 | 0.051 |

After Harmony PCA has the lower donor NMI and the better local mixing. Before Harmony
the arms are not a like-for-like batch-robustness comparison, since the PCA input is
already partly donor-adjusted, so no claim is made about it. The two corrected arms
agree moderately: ARI 0.543 and NMI 0.641 between `pca_harmony` and `scgpt_harmony`.
`scgpt_harmony` shifted slightly from the earlier run (975 cells, ARI 0.958); `scgpt_raw`
is identical.

**Caveat.** Harmony removes all between-patient variation, and with it every age and
diagnosis signal. These coordinates are for grouping cells by type only.

**Code.** [`02_CODE/06_benchmark_representations.py`](../../02_CODE/06_benchmark_representations.py) → `.R`;
the uncorrected donor-coloured layout comes from [`02_CODE/09_pca_raw_umap.R`](../../02_CODE/09_pca_raw_umap.R).

**Canonical folder.** [`Age_trajectory/results/06_harmony_benchmark/`](../../Age_trajectory/results/06_harmony_benchmark/)

| linked file | what it is |
| --- | --- |
| `06_pca_raw_umap.png` | PCA without Harmony, 10 clusters. |
| `06_pca_raw_umap_by_donor.png` | The same layout coloured by patient — donor NMI 0.049. |
| `06_pca_harmony_umap.png` | PCA after Harmony, 12 clusters. |
| `06_scgpt_harmony_umap.png` | scGPT after Harmony, 11 clusters. |
| `benchmark_arm_summary.csv` | Cluster counts, sizes, donor mixing, per arm. |
| `benchmark_arm_agreement.csv` | ARI and NMI between every pair of arms. |
| `benchmark_donor_leakage.csv` | How much of each clustering patient identity explains. |
| `benchmark_cluster_assignments.csv.gz` | All five arms' labels plus the corrected UMAP layouts, per cell. |

Not linked: `benchmark_representations.rds` (all five reductions),
`benchmark_pca_uncorrected.rds`, `benchmark_pca_raw_umap_coordinates.csv.gz`,
`r_input/`, `benchmark_environment_audit.json`. The PCA-512 files in this folder
belong to [06_DIMENSION_CONTROL](../06_DIMENSION_CONTROL/).
