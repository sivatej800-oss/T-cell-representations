# Status

*Last updated 2026-09-27.*

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

- **Reproduction:** 12 clusters (the paper reports 13). Requesting the full
  integration-feature pool is what recovers the rare populations.
- **Batch robustness:** before correction, scGPT leaked less patient identity than PCA
  (donor NMI 0.072 against 0.167). Harmony brings both to 0.051.
- **Separation:** at matched 512 dimensions, scGPT's partition scores 0.219 in its own
  space and PCA-512's scores −0.033 in its own. All 11 scGPT clusters score positive,
  against 4 of 9 for PCA-512. Adding dimensions to PCA loses the rare populations.
- **Biology:** both arms resolve to the same five broad cell types (T, monocyte, DC,
  NK, B), with counts within 83 cells of each other. PCA's two extra clusters are
  subdivisions of T cells.

## Current interpretation

At the broad cell-type level, PCA and scGPT agree on the biology. They differ in
geometry: scGPT gives a better-separated, more batch-robust space, and that advantage
is not an artifact of its dimensionality. Whether the arms differ in finer biology,
meaning T-cell and monocyte subtypes, is still untested.

## Next analyses

1. **T-cell subclustering:** CD4, CD8, CD4/CD8, Treg, for each arm. See
   [09_T_CELL_SUBCLUSTERS](../03_RESULTS/09_T_CELL_SUBCLUSTERS/).
2. **Monocyte subclustering:** classical, intermediate, non-classical, for each arm.
   See [10_MONOCYTE_SUBCLUSTERS](../03_RESULTS/10_MONOCYTE_SUBCLUSTERS/).
3. Write the refined labels back onto the full UMAP to produce Figure 1D itself.
4. Pick the figures for [FINAL_FIGURES](../03_RESULTS/FINAL_FIGURES/) and draft the
   [paper](../04_PAPER/).

## Project scope

**In scope:** the Figure 1D reproduction, and the PCA-vs-scGPT comparison on
cell-state structure (donor mixing, separation, the dimensionality control, and marker
and cell-type agreement).

**Out of scope:** patient-level aging and cognitive-impairment analyses, meaning
age trajectories and deviation of cognitively impaired patients from healthy
expectation. This was decided on 2026-09-27. Also out of scope: re-running the
authors' upstream preprocessing (Cell Ranger, SoupX, DoubletFinder). The public
matrix is already its output.
