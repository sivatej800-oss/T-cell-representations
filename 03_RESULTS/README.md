# Results

Every result this project has produced, grouped by the scientific question it answers.

**This folder is a navigation layer.** Each file here is a relative link into the
pipeline workspace, [`Age_trajectory/results/`](../Age_trajectory/results/), where the
canonical files live in folders named after the stage that wrote them. Nothing here is
a copy, so nothing here can drift out of date. Large checkpoints (`.rds`,
`sct_samples/`, `seurat_input/`, `r_input/`) are deliberately not linked; they stay
in the workspace for resuming the pipeline.

The links only resolve on a machine where the pipeline has been run. On a fresh clone
they are dangling until it has.

| folder | question | stages | canonical folder(s) |
| --- | --- | --- | --- |
| [01_PCA_BASELINE](01_PCA_BASELINE/) | Does the Figure 1D clustering workflow reproduce from the public data? | 01–03 | `03_seurat_reproduction` |
| [02_SCGPT_BASELINE](02_SCGPT_BASELINE/) | What does the same workflow give on a frozen scGPT embedding? | 04–05 | `04_scgpt_embedding` |
| [03_BATCH_CORRECTION](03_BATCH_CORRECTION/) | With Harmony applied to both, which leaks less patient identity? | 06, 09 | `06_harmony_benchmark` |
| [04_SILHOUETTE](04_SILHOUETTE/) | Which partition is more compact and separated, at native widths? | 07 | `07_silhouette` |
| [05_ELBOW_WCSS](05_ELBOW_WCSS/) | Where does resolution 0.3 sit on each compactness curve? *(supporting)* | 08 | `08_wcss_sweep` |
| [06_DIMENSION_CONTROL](06_DIMENSION_CONTROL/) | Is scGPT's advantage just its 512 dimensions? | 10–11 | `06_harmony_benchmark`, `11_silhouette_matched` |
| [07_CELL_TYPE_MARKERS](07_CELL_TYPE_MARKERS/) | What does the RNA say each frozen cluster is? | 12–14 | `12_frozen_labels`, `13_marker_discovery`, `14_marker_evidence` |
| [08_CELL_TYPE_ANNOTATION](08_CELL_TYPE_ANNOTATION/) | Which broad cell type is each cluster, in each arm? | 15 | `15_annotation` |
| [09_T_CELL_SUBCLUSTERS](09_T_CELL_SUBCLUSTERS/) | **PLANNED / NOT YET COMPLETED** — CD4, CD8, Treg by reclustering | — | — |
| [10_MONOCYTE_SUBCLUSTERS](10_MONOCYTE_SUBCLUSTERS/) | **PLANNED / NOT YET COMPLETED** — monocyte subtypes by reclustering | — | — |
| [FINAL_FIGURES](FINAL_FIGURES/) | The figures chosen for the paper | — | — |

The full scientific narrative behind every number is in
[`01_PROJECT_OVERVIEW/ANALYSIS_WALKTHROUGH.md`](../01_PROJECT_OVERVIEW/ANALYSIS_WALKTHROUGH.md).

## Reproducing any of it

Each stage reads the previous stage's saved output, so a stage can be re-run on its
own. The Python driver of each pair sets up paths and calls the R script. Run from the
repository root.

```bash
python 02_CODE/03_reproduce_figure1d.py         # 01_PCA_BASELINE
python 02_CODE/04_scgpt_figure1d.py             # 02_SCGPT_BASELINE
python 02_CODE/05_audit_scgpt_gene_coverage.py
python 02_CODE/06_benchmark_representations.py  # 03_BATCH_CORRECTION
Rscript 02_CODE/09_pca_raw_umap.R \
  Age_trajectory/results/06_harmony_benchmark \
  Age_trajectory/results/figures 42
python 02_CODE/07_silhouette_comparison.py      # 04_SILHOUETTE
python 02_CODE/08_wcss_resolution_sweep.py      # 05_ELBOW_WCSS
python 02_CODE/10_pca512_arm.py                 # 06_DIMENSION_CONTROL
python 02_CODE/11_silhouette_dimension_matched.py
python 02_CODE/12_freeze_cluster_labels.py      # 07_CELL_TYPE_MARKERS
python 02_CODE/13_marker_discovery.py
python 02_CODE/14_marker_evidence.py
Rscript 02_CODE/15_annotate_clusters.R \
  Age_trajectory/results/12_frozen_labels \
  Age_trajectory/results/14_marker_evidence \
  Age_trajectory/results/06_harmony_benchmark \
  Age_trajectory/results/15_annotation \
  Age_trajectory/results/figures                # 08_CELL_TYPE_ANNOTATION
```
