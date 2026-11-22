# Status

*Last updated 2026-11-22.*

## Research question

Which way of representing a cell, **PCA on SCTransformed expression** or a **frozen
scGPT embedding**, organizes CSF immune cells into cleaner, donor-independent cell
populations? The data are GSE200164 (Piehl et al. 2022, *Cell*): 70,391 cells from 59
people. The baseline is a reproduction of the clustering workflow behind the paper's
Figure 1D.

## Completed analyses

| # | analysis | results |
| --- | --- | --- |
| 1 | Figure 1D clustering reproduced (PCA arm) | [01_PCA_BASELINE](../03_RESULTS/01_PCA_BASELINE/) |
| 2 | Same workflow on a frozen scGPT embedding, plus a gene-coverage audit | [02_SCGPT_BASELINE](../03_RESULTS/02_SCGPT_BASELINE/) |
| 3 | Harmony on both arms, measuring donor mixing and leakage | [03_BATCH_CORRECTION](../03_RESULTS/03_BATCH_CORRECTION/) |
| 4 | Silhouette at native widths | [04_SILHOUETTE](../03_RESULTS/04_SILHOUETTE/) |
| 5 | WCSS elbow sweep *(supporting)* | [05_ELBOW_WCSS](../03_RESULTS/05_ELBOW_WCSS/) |
| 6 | PCA-512 control and silhouette at matched width | [06_DIMENSION_CONTROL](../03_RESULTS/06_DIMENSION_CONTROL/) |
| 7 | Frozen labels, marker discovery, marker evidence | [07_CELL_TYPE_MARKERS](../03_RESULTS/07_CELL_TYPE_MARKERS/) |
| 8 | Broad cell-type annotation, both arms, by hand | [08_CELL_TYPE_ANNOTATION](../03_RESULTS/08_CELL_TYPE_ANNOTATION/) |

## Headline results

The PCA arm was rebuilt on 2026-11-22 from SCT residuals, the quantity the authors'
integration consumes, so it follows their normalization as closely as a comparison
without anchors allows. The earlier arm used the SCT `data` layer and was dominated by
depth and gene count; every PCA number below is from the rebuilt arm. Earlier results
are in `Age_trajectory/results_superseded_2026-11-22/`.

- **Reproduction:** 12 clusters (the paper reports 13). Requesting the full
  integration-feature pool is what recovers the rare populations.
- **Donor leakage:** the PCA arm *normalized per patient, no Harmony* has donor NMI
  0.049 against scGPT's 0.072. That is not a like-for-like comparison: residuals are
  centred within each patient, so the PCA arm is partly donor-adjusted. After Harmony,
  donor NMI is 0.038 for PCA and 0.051 for scGPT, and PCA mixes donors better locally
  (0.969 against 0.866).
- **Separation:** each partition scores best in its own space. At native width, PCA
  labels score 0.176 in PCA space against scGPT labels' 0.151, and scGPT labels score
  0.220 in scGPT space against PCA labels' 0.100. At matched 512 dimensions it is the
  same: PCA-512 labels 0.053 against −0.022 in PCA-512 space, scGPT labels 0.221
  against 0.136 in scGPT space. Neither description wins in both spaces.
- **Biology:** both arms resolve to the same broad cell types (T, monocyte, DC, NK,
  B), within 119 cells for every identity. 821 PCA cells (875 in scGPT) are
  Undetermined: a cluster dominated by XIST and long-transcript genes. PCA cluster 3
  is macrophage-like and filed under Monocytes.

## Current interpretation

With PCA built the way the authors normalized, the earlier claim that scGPT separates
clusters better and is more batch-robust does not stand. The two descriptions give
comparable geometry in their own spaces, and PCA is better at donor mixing after
Harmony. Whether they differ in biology, at broad level or in T-cell and monocyte
subtypes, is untested.

## Next analyses

1. **T-cell subclustering:** CD4, CD8, CD4/CD8, Treg, for each arm. See
   [09_T_CELL_SUBCLUSTERS](../03_RESULTS/09_T_CELL_SUBCLUSTERS/).
2. **Monocyte subclustering:** classical, intermediate, non-classical, for each arm.
   See [10_MONOCYTE_SUBCLUSTERS](../03_RESULTS/10_MONOCYTE_SUBCLUSTERS/).
3. Write the refined labels back onto the full UMAP to produce Figure 1D itself.
4. Pick the figures for [FINAL_FIGURES](../03_RESULTS/FINAL_FIGURES/) and draft the
   [paper](../04_PAPER/). A paper outline and abstract already exist in
   [`04_PAPER/`](../04_PAPER/).

## Project scope

**In scope:** the Figure 1D reproduction, and the PCA-vs-scGPT comparison on
cell-state structure (donor mixing, separation, the dimensionality control, and marker
and cell-type agreement).

**Out of scope:** patient-level aging and cognitive-impairment analyses, meaning
age trajectories and deviation of cognitively impaired patients from healthy
expectation. This was decided on 2026-09-27. Also out of scope: re-running the
authors' upstream preprocessing (Cell Ranger, SoupX, DoubletFinder). The public
matrix is already its output.
