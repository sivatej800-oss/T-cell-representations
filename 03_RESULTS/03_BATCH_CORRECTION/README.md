# 03 — Batch correction (Harmony on both arms)

**Question.** The two baselines are not comparable: PCA was anchor-integrated, scGPT
was not corrected at all. With Harmony applied equally to both, how much patient
identity does each still carry?

**Result.**

| arm | dims | clusters | smallest | donor mixing | donor NMI |
| --- | --- | --- | --- | --- | --- |
| `pca_raw` | 15 | 11 | 490 | 0.801 | **0.167** |
| `pca_harmony` | 15 | 12 | 109 | **0.959** | 0.051 |
| `pca_harmony_30dims` | 30 | 13 | 139 | 0.950 | 0.051 |
| `scgpt_raw` | 512 | 12 | 103 | 0.806 | **0.072** |
| `scgpt_harmony` | 512 | 11 | 110 | 0.867 | 0.051 |

Before correction, scGPT was the more batch-robust arm (patient NMI 0.072 against
0.167). Harmony brings both to 0.051, and correcting PCA brought back its rare
populations. The two corrected arms agree only moderately: ARI 0.584 and NMI 0.694
between `pca_harmony` and `scgpt_harmony`.

**Caveat.** Harmony removes all between-patient variation, and with it every age and
diagnosis signal. These coordinates are for grouping cells by type only.

**Code.** [`02_CODE/06_benchmark_representations.py`](../../02_CODE/06_benchmark_representations.py) → `.R`;
the uncorrected donor-coloured layout comes from [`02_CODE/09_pca_raw_umap.R`](../../02_CODE/09_pca_raw_umap.R).

**Canonical folder.** [`Age_trajectory/results/06_harmony_benchmark/`](../../Age_trajectory/results/06_harmony_benchmark/)

| linked file | what it is |
| --- | --- |
| `06_pca_raw_umap.png` | Uncorrected PCA control, 11 clusters. |
| `06_pca_raw_umap_by_donor.png` | The same layout coloured by patient — what NMI 0.167 looks like. |
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
