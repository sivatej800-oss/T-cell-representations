# 07 — Cell-type markers

**Question.** Geometry cannot tell whether a cluster is a cell type. What does the
RNA say each cluster is?

**Design.** The labels go in, and the RNA provides the evidence:
- **Freeze the labels (Stage 12).** The two partitions under test are
  `pca30_harmony_cluster` (13 clusters) and `scgpt_harmony_cluster` (11). They are
  frozen and checksummed, and are never re-clustered after this point.
- **Find markers (Stage 13).** One shared RNA object is normalized once.
  `FindAllMarkers` then runs twice on it with identical parameters (the authors'
  thresholds, Wilcoxon via presto), so the only thing that differs between the two
  runs is which cells carry which label.
- **Plot the evidence (Stage 14).** The authors' 14-gene canonical panel is drawn over
  each arm's own UMAP, with the percentages behind every colour.

**Why RNA and not SCT.** The per-sample SCT assays share only 9,013 of 19,689 genes,
and the genes that drop out are the rare-population markers.

**Code.** [`02_CODE/12_freeze_cluster_labels.py`](../../02_CODE/12_freeze_cluster_labels.py),
[`02_CODE/13_marker_discovery.py`](../../02_CODE/13_marker_discovery.py) → `.R`,
[`02_CODE/14_marker_evidence.py`](../../02_CODE/14_marker_evidence.py) → `.R`.

**Canonical folders.** [`12_frozen_labels/`](../../Age_trajectory/results/12_frozen_labels/),
[`13_marker_discovery/`](../../Age_trajectory/results/13_marker_discovery/),
[`14_marker_evidence/`](../../Age_trajectory/results/14_marker_evidence/)

| linked file | what it is |
| --- | --- |
| `frozen_cluster_labels.csv.gz` | Both frozen partitions, per cell. |
| `markers_pca30_harmony_cluster.csv.gz` | FindAllMarkers output, PCA-30 partition. |
| `markers_scgpt_harmony_cluster.csv.gz` | The same, scGPT partition. |
| `cluster_evidence_table.csv` | One row per cluster: size, donors, top markers, panel genes present. This is what annotation is read from. |
| `top_markers_*.csv` | Top markers per cluster, each arm. |
| `14_featureplots_*.png` | The 14-gene panel over each arm's UMAP. |
| `14_clusters_*.png` | Each arm's cluster numbers on its own layout. |
| `14_feature_plots/` | One plot per marker per arm. |

Not linked: `frozen_labels_audit.json` (holds the checksum), `canonical_panel_by_cluster.csv`,
`canonical_pct_expressing_*.csv`, and the marker and environment audits.
