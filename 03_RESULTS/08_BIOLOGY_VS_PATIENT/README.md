# 08 — Biology vs. patient (SingleR reference labels)

**Question.** Does each arm's clustering follow cell type or the patient a cell came
from? The arms are the authors' published pipeline (Stage 03) and scGPT used out of
the box (`scgpt_raw`).

**Method.** SingleR labelled every cell against purified immune profiles
(`celldex::MonacoImmuneData`, main labels). The labels come from expression alone and
belong to neither arm. 530 cells (0.75%) stayed unassigned and are left out. Each
arm's clusters are then scored against SingleR (biology NMI) and against patient
(patient NMI), and each hand annotation from [07](../07_CELL_TYPE_ANNOTATION/) is
checked cell by cell against SingleR.

The authors' pipeline includes batch correction (anchor integration) and scGPT raw
has none. That is intended, so a lower patient NMI for the authors' arm is expected.

**Result.**
- **Biology dominates both arms.** Biology NMI is 0.50–0.60; patient NMI is 0.07 or
  less.
- **scGPT follows cell type more closely** (biology NMI 0.587 against 0.501).
- **The authors' pipeline carries less patient identity** (0.039 against 0.072).
- **Agreement with SingleR is about equal** (95.8% against 95.3%; 99.7% against
  98.9% with monocytes and DCs merged). scGPT leads on NK cells and DCs, and is the
  only arm with a B-cell cluster.
- **Both arms disagree with SingleR on the same boundary:** about 2,500 cells on the
  monocyte/DC line.

**Code.** `02_CODE/16_singler_reference_labels.py`, `02_CODE/17_cluster_structure_vs_singler.py` → `.R`.

**Canonical folders.** [`Age_trajectory/results/16_singler_labels/`](../../Age_trajectory/results/16_singler_labels/), [`Age_trajectory/results/17_biology_vs_donor/`](../../Age_trajectory/results/17_biology_vs_donor/)

| file | what it is |
| --- | --- |
| `17_patient_vs_biology_scatter.png` | Patient NMI against biology NMI, one point per arm. |
| `17_agreement_heatmap_authors_pipeline.png` | Hand annotation against SingleR, authors' pipeline. |
| `17_agreement_heatmap_scgpt_raw.png` | Hand annotation against SingleR, scGPT raw. |
| `17_umap_singler_vs_annotation_*.png` | Each arm's UMAP coloured by SingleR, by the annotation, and by kind of disagreement. |
| `disagreement_kinds_*.csv` | Every kind of disagreement (annotation → SingleR) with its cell count, per arm. |
| `nmi_by_arm.csv` | Patient and biology NMI per arm. |
| `annotation_vs_singler_agreement.csv` | Overall agreement, monocyte/DC apart and merged. |
| `agreement_per_type_*.csv` | Agreement per broad cell type. |
| `confusion_*.csv`, `confusion_merged_myeloid_*.csv` | Hand annotation against SingleR, cell counts, per arm. |
| `monocyte_dc_marker_groups.csv` | Marker percentages for the disputed monocyte/DC cells. |
| `label_counts_broad.csv` | SingleR's broad label counts. |

The earlier four-arm run with the Harmony arms is kept in
`Age_trajectory/results/17_biology_vs_donor_harmony_record/`.
