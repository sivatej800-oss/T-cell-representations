# GSE200164 Figure 1D Reproduction

This repository reproduces the clustering workflow used for Figure 1D of Piehl et al. using the public GSE200164 single-cell RNA-seq dataset.

The dataset contains **70,391 CSF immune cells from 59 individuals**. The current reproduction generates an integrated UMAP with **11 unsupervised clusters**. The published figure is a later, annotated result in which the authors used marker genes and separate reclustering of T cells and monocytes to assign final cell-type labels.

![Current unsupervised UMAP](Age_trajectory/results/figure1d_python/figure1d_unsupervised_clusters.png)

## What the authors did

### 1. Preprocess the sequencing data

The authors began with raw sequencing reads and:

1. aligned reads to GRCh38 with Cell Ranger 6.0.0;
2. removed empty droplets;
3. corrected ambient RNA contamination with SoupX;
4. created one Seurat object per sample;
5. removed doublets with DoubletFinder;
6. removed cells with fewer than 200 detected genes or more than 10% mitochondrial reads.

The public expression matrix used here is already the output of these preprocessing steps.

### 2. Normalize each sample independently

The authors split the merged Seurat object by sample and ran `SCTransform` separately on each sample.

Important settings:

- 1,000 variable genes per sample;
- regression of `nCount_RNA`, `nFeature_RNA`, and `percent.mt`;
- up to 10,000 cells used during model fitting.

### 3. Integrate the 59 samples

The authors selected shared integration features, prepared the SCT assays, found integration anchors, and combined the samples into one integrated assay.

Samples A4 and A6 were used as references to reduce memory use during anchor finding.

### 4. Generate the first cell representation with PCA

PCA was run on the integrated assay. The first **15 principal components** were used for both clustering and UMAP.

### 5. Build the cell-neighbor graph and cluster the cells

The authors constructed a shared nearest-neighbor graph from PCs 1–15 and clustered it with the Louvain algorithm at:

```text
resolution = 0.3
```

Their software environment produced 13 first-pass numerical Seurat clusters.

### 6. Generate the UMAP

UMAP was calculated from PCs 1–15 using:

```text
n.neighbors = 30
```

UMAP is used to visualize the cells. The clustering itself is determined from the neighbor graph, not from the two-dimensional UMAP coordinates.

### 7. Identify broad immune-cell populations

The authors calculated cluster marker genes with `FindAllMarkers` and inspected canonical markers including:

```text
CD3D, CD4, CD8A, CD8B, FOXP3
NKG7, GNLY
FCER1A
CD68, CD14, FCGR3A
CD19, JCHAIN
```

The initial numerical clusters were grouped into broad T-cell, monocyte, NK-cell, dendritic-cell, B-cell, plasma-cell, and undetermined compartments.

### 8. Recluster T cells and monocytes

The authors then isolated the broad T-cell and monocyte compartments and repeated PCA, neighbor finding, clustering, and UMAP within each subset.

This produced more specific labels, including:

- CD4+ T cells;
- CD8+ T cells;
- CD4+/CD8+ T cells;
- regulatory T cells;
- classical, intermediate, and non-classical monocyte populations.

These refined labels were written back onto the original full-dataset UMAP to create Figure 1D.

## What this repository currently does

Our reproduction starts from the processed GEO expression matrix rather than the raw FASTQ files. Because that matrix is already the output of the authors' step 1, this pipeline picks up at their step 2 and carries through their step 6.

The work is split across three scripts, run in order. Each writes its result to disk and the next reads that file, so a stage can be re-run without repeating the ones before it.

```text
  Stage 01  ──▶  Stage 02  ──▶  Stage 03  ──▶  Figure 1D
  acquire        prepare        analyze        11 clusters
  (Python)       (Python)       (R, via Python)
```

### Stage 01 — acquire the data and prove it is intact

`scripts/01_download_audit_build_anndata.py`

1. downloads the GEO expression matrix and the series metadata file, verifying each archive's gzip integrity and recording a SHA-256 digest so a later run can prove it started from identical source data;
2. converts the very large dense CSV into a sparse matrix, storing only nonzero measurements, written in both row and column orientations;
3. re-reads the original download and verifies the conversion changed no value — everything downstream depends on this holding;
4. binds the measurements, gene names, and cell barcodes into a single AnnData object;
5. parses the patient table from the series matrix and copies each person's age, sex, and diagnosis onto every one of their cells.

The SoupX-adjusted values are fractional, not integer counts. They are preserved exactly at every step; rounding them would silently undo the authors' ambient-RNA correction.

### Stage 02 — measure cell quality and prepare the hand-off

`scripts/02_prepare_figure1d.py`

1. opens Stage 01's result and confirms its dimensions are unchanged;
2. computes the three quality measures the authors regressed out — `nCount_RNA`, `nFeature_RNA`, and `percent.mt`;
3. reconciles the two naming systems: GEO's accession-based sample names against the authors' `A1`–`H8` labels, and writes a 59-row summary of the study;
4. saves its own checkpoint so a re-run need not repeat the above;
5. exports the matrix in Matrix Market form for R — skipped when Stage 03 already has a checkpoint to resume from, since the export is roughly 4 GB;
6. confirms R and its packages are present now, rather than discovering they are not once Stage 03 is hours into the analysis.

**No cells are removed.** The thresholds are recomputed and reported, but nothing is filtered: `percent.mt` maxes at 9.9979 and `nFeature_RNA` bottoms at 202, the fingerprints of the authors' filters having already run. Filtering again would diverge from the published dataset.

### Stage 03 — the clustering workflow

`scripts/03_reproduce_figure1d.R`, started automatically by `scripts/03_reproduce_figure1d.py`

This is the direct counterpart of the authors' steps 2 through 6:

1. rebuilds the Seurat object with `min.cells = 0, min.features = 0` — deliberately no thresholds, since the authors' already ran;
2. splits the cells into 59 per-sample objects, asserting none are lost;
3. runs `SCTransform` per sample with 1,000 variable genes, `ncells = 10000`, and regression of the three quality measures — matching the authors, with `vst.flavor = "v1"` set explicitly (see below);
4. selects 2,000 shared integration features and prepares the SCT assays;
5. finds integration anchors using samples **A4 and A6** as references, located by name rather than by position;
6. integrates the samples into one assay, asserting the cell count is unchanged;
7. runs PCA on the integrated assay, 50 components;
8. builds the shared nearest-neighbor graph from PCs 1–15;
9. clusters with Louvain at `resolution = 0.3`;
10. runs UMAP from PCs 1–15 with `n.neighbors = 30`;
11. writes the full Seurat object, a per-cell CSV of coordinates and cluster assignments, and the figure.

Because steps 3 through 6 take hours, each is checkpointed. `--restart-from {none,clustering,pca,integration,all}` controls how far back a re-run goes: re-tuning the clustering costs minutes, while `all` redoes 59 `SCTransform` runs and costs hours.

### Result and the unexplained cluster-count difference

The result contains all 70,391 cells but currently produces **11 numerical clusters** where the authors reported 13.

**Why this happens is currently unknown.** What can be stated is only the following:

- every published parameter is reproduced as written — 1,000 variable genes, the same three regressed covariates, 2,000 integration features, references A4 and A6, PCs 1–15, resolution 0.3, `n.neighbors = 30`;
- the pipeline starts from the same deposited matrix the authors' own Figure 1D script starts from;
- all 70,391 cells are present, so the difference is not caused by dropped or filtered cells.

A plausible explanation is that the Seurat implementations of `SCTransform` and Louvain clustering have changed since October 2021, moving where cluster boundaries fall. **That is a hypothesis, not a finding — it has not been tested here.** Doing so would mean re-running the pipeline pinned to the authors' original package versions and comparing, which has not been attempted. Other explanations have not been excluded: an undocumented parameter, a difference in the object the authors actually started from, or an unrecorded step between their preprocessing and their Figure 1D script.

Two decisions follow from not knowing:

- **`vst.flavor = "v1"` is set deliberately.** The authors wrote their script before SCTransform v2 existed, so their unspecified call ran v1; modern Seurat defaults to v2. Removing this line would look like a closer match to their code while running a different algorithm. This narrows one known source of divergence — it does not resolve the count.
- **The resolution is not tuned to force 13 clusters.** Raising it until the number matched would report the published count obtained by a different method, and would bury the discrepancy rather than record it.

**We are proceeding on this basis.** 11 clusters is the honest current output, the downstream comparison work does not depend on the count matching, and the difference is documented here rather than tuned away. If the cause is identified later, this section should be updated to say so.

The authors' steps 7 and 8 — marker-gene identification and reclustering of T cells and monocytes — are not yet reproduced. Those are the subject of the labeling work described below.


## Future scGPT comparison

The next comparison will replace the PCA-based cell representation with a frozen scGPT embedding while keeping the downstream workflow as similar as possible.

Planned steps:

1. load the same AnnData object and retain the same cells and metadata;
2. generate one scGPT embedding vector for each cell;
3. store the vectors in `adata.obsm["X_scgpt"]`;
4. build a nearest-neighbor graph directly from the scGPT vectors;
5. cluster the graph using the same resolution target;
6. generate a UMAP using the same neighbor setting;
7. compare the PCA and scGPT results on the same cells.

The comparison will evaluate:

- cluster stability;
- separation of broad immune populations;
- mixing of cells from different donors;
- agreement with marker-gene expression;
- recovery of rare or transitional populations.

SCTransform integration and PCA define the original representation. In the scGPT branch, those coordinates are replaced by the scGPT vectors; the graph clustering, UMAP, and biological interpretation steps remain parallel.

## Future cluster labeling from gene expression

Cluster labels should be assigned from the original normalized RNA expression values, not from UMAP coordinates or the embedding dimensions alone.

For both the PCA and scGPT clusters, the planned annotation workflow is:

1. run differential-expression analysis for every cluster against all remaining cells;
2. rank the positively enriched genes for each cluster;
3. inspect marker heatmaps and expression plots;
4. assign broad labels from combinations of canonical markers;
5. validate that marker expression is consistent across multiple donors;
6. isolate broad T-cell and monocyte populations and recluster them for finer labels;
7. keep uncertain populations labeled as ambiguous rather than forcing a cell type.

This will allow the PCA and scGPT clusterings to be compared using the same biological evidence.

## References

- Piehl N, et al. *Cerebrospinal fluid immune dysregulation during healthy brain aging and cognitive impairment.* Cell. 2022. DOI: 10.1016/j.cell.2022.11.019.
- Dataset: GEO accession **GSE200164**.
- Author code: `gatelabNW/csf_aging`.
