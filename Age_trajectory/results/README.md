# Results

Every figure and table this project has produced, grouped by the stage that wrote it.

Each stage directory holds that stage's data — checkpoints, saved objects, audits.
**All figures are collected in [`figures/`](figures/)** rather than sitting beside
that data, so the report can link them from one place.

The headline figures are embedded in the [main README](../../README.md); this page
is the full inventory.

## Figures

| figure | stage | what it shows |
| --- | --- | --- |
| [`03_seurat_umap_clusters.png`](figures/03_seurat_umap_clusters.png) | 03 | The paper reproduction — 70,391 cells, 12 clusters, anchor-integrated PCA. |
| [`04_scgpt_raw_umap_clusters.png`](figures/04_scgpt_raw_umap_clusters.png) | 04 | The same cells from a frozen scGPT embedding, 13 clusters, no batch correction. |
| [`06_pca_raw_umap.png`](figures/06_pca_raw_umap.png) | 06 | The uncorrected PCA control, 11 clusters. Read with the donor panel below. |
| [`06_pca_raw_umap_by_donor.png`](figures/06_pca_raw_umap_by_donor.png) | 06 | The same layout coloured by patient — what "uncorrected" costs. |
| [`06_pca_harmony_umap.png`](figures/06_pca_harmony_umap.png) | 06 | PCA after Harmony, 12 clusters. |
| [`06_scgpt_harmony_umap.png`](figures/06_scgpt_harmony_umap.png) | 06 | scGPT after Harmony, 11 clusters. |
| [`07_silhouette_per_cluster.png`](figures/07_silhouette_per_cluster.png) | 07 | Per-cluster silhouette width, both partitions in both spaces. |
| [`10_pca512_harmony_umap.png`](figures/10_pca512_harmony_umap.png) | 10 | PCA carried to 512 dimensions: 9 clusters, rare populations lost. |
| [`11_silhouette_matched.png`](figures/11_silhouette_matched.png) | 11 | The matched 2×2 — 512 vs 512, dimensionality controlled. |
| [`11_silhouette_matched_per_cluster.png`](figures/11_silhouette_matched_per_cluster.png) | 11 | Per-cluster width for both matched partitions in their own spaces. |
| [`08_wcss_compactness.png`](figures/08_wcss_compactness.png) | 08 | Normalized WCSS against cluster count — the elbow. |
| [`08_wcss_dispersion_explained.png`](figures/08_wcss_dispersion_explained.png) | 08 | The same sweep as dispersion explained. |

## Stages

### [`03_seurat_reproduction/`](03_seurat_reproduction/) — the paper reproduction

SCTransform per sample, anchor integration, PCA, Louvain at 0.3. Written by
`scripts/03_reproduce_figure1d.R`.

| file | what it is |
| --- | --- |
| `GSE200164_umap_clusters.csv.gz` | Cluster labels and UMAP coordinates per cell. |
| `GSE200164_umap_clusters_BEFORE_2000features.csv.gz` | The superseded 2,000-feature run, kept as the counter-example. |
| `integrated_clustered_umap.rds` | The final Seurat object. |
| `GSE200164_after_PCA.rds`, `GSE200164_after_qc_audit.h5ad` | Stage checkpoints, for resuming without re-running. |
| `sct_samples/` | Per-sample SCTransform output, 59 files. |
| `seurat_input/` | The Matrix Market hand-off from Python to R. |
| `figure1d_environment_audit.json` | Package versions and seeds for this run. |

### [`04_scgpt_embedding/`](04_scgpt_embedding/) — the frozen-embedding arm

The same downstream code on a 512-dimension scGPT embedding instead of PCA.
Written by `scripts/04_scgpt_figure1d.R` and audited by `scripts/05_audit_scgpt_gene_coverage.py`.

| file | what it is |
| --- | --- |
| `GSE200164_scgpt_umap_clusters.csv.gz` | Cluster labels and UMAP coordinates per cell. |
| `scgpt_clustered_umap.rds` | The final Seurat object. |
| `scgpt_gene_coverage_audit.json` | Which genes scGPT's vocabulary matched, and the expression behind the gap. |
| `r_input/` | The embedding and metadata handed to R. |
| `scgpt_environment_audit.json` | Package versions and seeds for this run. |

### [`06_harmony_benchmark/`](06_harmony_benchmark/) — the five-arm comparison

Harmony applied to both representations equally, so the two can be compared without
one being batch-corrected and the other not. Written by `scripts/06_benchmark_representations.R`.

| file | what it is |
| --- | --- |
| `benchmark_arm_summary.csv` | Cluster counts, smallest and largest cluster, donor mixing, per arm. |
| `benchmark_arm_agreement.csv` | ARI and NMI between every pair of arms. |
| `benchmark_donor_leakage.csv` | How much of each clustering is explained by patient identity. |
| `benchmark_cluster_assignments.csv.gz` | All five arms' labels plus the two corrected UMAP layouts, per cell. |
| `benchmark_pca_raw_umap_coordinates.csv.gz` | The `pca_raw` layout, added by Stage 09. |
| `benchmark_pca512_assignments.csv.gz` | Stage 10's PCA-512 labels and UMAP layout, per cell. |
| `benchmark_pca512_summary.csv` | PCA-512 cluster counts and the variance split across components. |
| `benchmark_pca512_arm.rds` | The Stage 10 Seurat object, carrying `pca512` and `pca512_harmony`. |
| `benchmark_representations.rds` | The Seurat object carrying all five reductions. |
| `benchmark_pca_uncorrected.rds` | The merged, unintegrated PCA checkpoint. |
| `r_input/` | The embedding and metadata handed to R. |
| `benchmark_environment_audit.json` | Package versions and seeds for this run. |

### [`07_silhouette/`](07_silhouette/) — per-cell separation

Each partition scored in both spaces, since distances in 30 principal components and
512 learned dimensions are not the same quantity. Written by `scripts/07_silhouette_comparison.R`.

| file | what it is |
| --- | --- |
| `silhouette_summary.csv` | Mean silhouette and negative fraction, per partition per space. |
| `silhouette_per_cluster.csv` | The same broken out by cluster, with cluster sizes. |
| `silhouette_environment_audit.json` | Package versions, seeds, and draw size for this run. |

### [`11_silhouette_matched/`](11_silhouette_matched/) — dimensionality-matched silhouette

PCA-512 against scGPT-512, so the comparison varies representation only. Written by
`scripts/11_silhouette_dimension_matched.R`; the PCA-512 arm it scores comes from
`scripts/10_pca512_arm.R` and lands in `06_harmony_benchmark/`.

| file | what it is |
| --- | --- |
| `silhouette_matched_summary.csv` | Mean silhouette and negative fraction, per partition per space. |
| `silhouette_matched_per_cluster.csv` | The same broken out by cluster, with cluster sizes. |
| `silhouette_matched_environment_audit.json` | Package versions, seeds, and draw size for this run. |

### [`08_wcss_sweep/`](08_wcss_sweep/) — the elbow

Louvain resolution swept 0.1 → 1.0, plotted against the cluster count each setting
happened to produce. Written by `scripts/08_wcss_resolution_sweep.R`.

| file | what it is |
| --- | --- |
| `wcss_resolution_sweep.csv` | Resolution, resulting cluster count, and normalized WCSS per representation. |
| `wcss_sweep_audit.json` | Package versions and seeds for this run. |

### [`12_frozen_labels/`](12_frozen_labels/) — the partitions under test

The two arms Stage 06 produced, copied out and frozen so the biological comparison
judges a fixed target. Written by `scripts/12_freeze_cluster_labels.py`.

| file | what it is |
| --- | --- |
| `frozen_cluster_labels.csv.gz` | `pca30_harmony_cluster` and `scgpt_harmony_cluster` per cell. |
| `frozen_labels_audit.json` | Cluster sizes and the SHA-256 later stages check the labels against. |

### [`13_marker_discovery/`](13_marker_discovery/) — what the RNA says the clusters are

Both partitions marker-tested against one shared expression object, so the only thing
differing between the two runs is which cells carry which label. Written by
`scripts/13_marker_discovery.R`.

| file | what it is |
| --- | --- |
| `markers_pca30_harmony_cluster.csv.gz` | FindAllMarkers output for the PCA-30 partition. |
| `markers_scgpt_harmony_cluster.csv.gz` | The same for the scGPT partition, identical parameters. |
| `canonical_panel_by_cluster.csv` | The paper's 13 marker genes, mean expression and percent expressing, per cluster. |
| `marker_discovery_audit.json` | Assay, thresholds, Wilcoxon implementation, and the expression fingerprint both runs shared. |

## Reproducing any of it

Each stage reads the previous stage's saved output, so a stage can be re-run on its
own. The Python driver of each pair sets up paths and calls the R script.

```bash
python scripts/03_reproduce_figure1d.py    # the paper reproduction
python scripts/04_scgpt_figure1d.py        # the embedding arm
python scripts/05_audit_scgpt_gene_coverage.py
python scripts/06_benchmark_representations.py
python scripts/07_silhouette_comparison.py
python scripts/08_wcss_resolution_sweep.py
Rscript scripts/09_pca_raw_umap.R \
  Age_trajectory/results/06_harmony_benchmark \
  Age_trajectory/results/figures 42
python scripts/10_pca512_arm.py              # the PCA-512 control
python scripts/11_silhouette_dimension_matched.py
python scripts/12_freeze_cluster_labels.py   # freeze the two partitions
python scripts/13_marker_discovery.py        # markers, both arms, one expression object
```
