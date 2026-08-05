# GSE200164 Figure 1D Reproduction

This repository reproduces the clustering workflow used for Figure 1D of Piehl et al. from the public GSE200164 single-cell RNA-seq dataset: **70,391 CSF immune cells from 59 individuals**.

The reproduction generates an integrated UMAP of unsupervised clusters. The published figure is a later, annotated result, in which the authors used marker genes and separate reclustering of T cells and monocytes to assign final cell-type labels.

The image below is from the earlier 11-cluster run, before the integration-feature correction described in [Result and the cluster-count difference](#result-and-the-cluster-count-difference); it will be replaced when the corrected pipeline is re-run.

![Current unsupervised UMAP](Age_trajectory/results/figure1d_python/figure1d_unsupervised_clusters.png)

## What the authors did

```mermaid
flowchart TD
    A["<b>1. Preprocess</b><br/>Cell Ranger 6.0.0 → GRCh38<br/>empty droplets removed<br/>SoupX ambient RNA correction<br/>DoubletFinder<br/>filter: ≥200 genes, &lt;10% mito"]
    B["<b>2. Normalize per sample</b><br/>SCTransform on each of 59 samples<br/>1,000 variable genes<br/>regress nCount_RNA, nFeature_RNA, percent.mt<br/>ncells = 10,000"]
    C["<b>3. Integrate</b><br/>shared features → anchors → one assay<br/>references: samples A4, A6"]
    D["<b>4. PCA</b><br/>PCs 1–15 retained"]
    E["<b>5. Cluster</b><br/>shared nearest-neighbor graph<br/>Louvain, resolution = 0.3<br/>→ 13 clusters"]
    F["<b>6. UMAP</b><br/>PCs 1–15, n.neighbors = 30"]
    G["<b>7. Identify populations</b><br/>FindAllMarkers + canonical markers<br/>→ T, mono, NK, DC, B, plasma"]
    H["<b>8. Recluster T and mono</b><br/>PCA, neighbors, clustering, UMAP<br/>within each subset"]
    I(["<b>Figure 1D</b><br/>refined labels written back<br/>onto the full-dataset UMAP"])

    A --> B --> C --> D
    D --> E
    D --> F
    E --> G --> H --> I
    F --> I

    style A fill:#e8e8e8,stroke:#666
    style I fill:#d4e8d4,stroke:#4a7a4a
```

The public expression matrix used here is already the output of step 1.

**PCA feeds clustering and UMAP separately.** UMAP visualizes the cells; the clusters come from the neighbor graph, not from the two-dimensional UMAP coordinates.

Markers used in step 7:

```text
CD3D, CD4, CD8A, CD8B, FOXP3    T cells and regulatory T cells
NKG7, GNLY                      NK cells
FCER1A                          dendritic cells
CD68, CD14, FCGR3A              monocytes
CD19, JCHAIN                    B cells and plasma cells
```

Step 8 produced the finer labels: CD4+, CD8+, CD4+/CD8+, and regulatory T cells; classical, intermediate, and non-classical monocytes.

## What this repository currently does

Our reproduction starts from the processed GEO matrix rather than the raw FASTQ files, so it picks up at the authors' step 2 and carries through their step 6.

```mermaid
flowchart LR
    subgraph authors ["Authors' workflow"]
        direction TB
        S1["1. Preprocess"]
        S2["2. Normalize"]
        S3["3. Integrate"]
        S4["4. PCA"]
        S5["5. Cluster"]
        S6["6. UMAP"]
        S7["7. Markers"]
        S8["8. Recluster"]
        S1 --> S2 --> S3 --> S4 --> S5 --> S6 --> S7 --> S8
    end

    subgraph ours ["This repository"]
        direction TB
        P1["deposited GEO matrix<br/><i>already the output of step 1</i>"]
        P2["Stages 01–03"]
        P3["not yet reproduced<br/><i>the labeling work below</i>"]
        P1 --> P2 --> P3
    end

    S1 -.-> P1
    S2 -.-> P2
    S6 -.-> P2
    S7 -.-> P3

    style S1 fill:#e8e8e8,stroke:#666
    style S7 fill:#f5e6d3,stroke:#a08050
    style S8 fill:#f5e6d3,stroke:#a08050
    style P1 fill:#e8e8e8,stroke:#666
    style P3 fill:#f5e6d3,stroke:#a08050
```

The work is split across three scripts, run in order. Each writes its result to disk and the next reads that file, so a stage can be re-run without repeating the ones before it.

```mermaid
flowchart LR
    GEO[("GEO<br/>GSE200164")]
    S1["<b>Stage 01</b> · Python<br/>acquire and audit"]
    S2["<b>Stage 02</b> · Python<br/>measure quality"]
    S3["<b>Stage 03</b> · R<br/>clustering workflow"]
    FIG(["<b>Figure 1D</b><br/>clustered UMAP"])

    C1[/"AnnData<br/>70,391 × genes"/]
    C2[/"Matrix Market<br/>~4 GB, for R"/]
    C3[/"per-stage checkpoints<br/>SCT · integration · PCA · clustering"/]

    GEO --> S1 --> C1 --> S2 --> C2 --> S3 --> FIG
    C3 -.resume.-> S3
    S3 -.-> C3

    style FIG fill:#d4e8d4,stroke:#4a7a4a
```

### Stage 01 — acquire the data and prove it is intact

`scripts/01_download_audit_build_anndata.py`

1. downloads the GEO expression matrix and series metadata, verifying each archive's gzip integrity and recording a SHA-256 digest, so a later run can prove it started from identical source data;
2. converts the very large dense CSV into a sparse matrix, storing only nonzero measurements, in both row and column orientations;
3. re-reads the original download and verifies the conversion changed no value — everything downstream depends on this holding;
4. binds measurements, gene names, and cell barcodes into a single AnnData object;
5. parses the patient table from the series matrix and copies each person's age, sex, and diagnosis onto every one of their cells.

The SoupX-adjusted values are fractional, not integer counts. They are preserved exactly at every step; rounding them would silently undo the authors' ambient-RNA correction.

### Stage 02 — measure cell quality and prepare the hand-off

`scripts/02_prepare_figure1d.py`

1. opens Stage 01's result and confirms its dimensions are unchanged;
2. computes the three quality measures the authors regressed out — `nCount_RNA`, `nFeature_RNA`, `percent.mt`;
3. reconciles GEO's accession-based sample names against the authors' `A1`–`H8` labels, and writes a 59-row summary of the study;
4. saves its own checkpoint so a re-run need not repeat the above;
5. exports the matrix in Matrix Market form for R — skipped when Stage 03 already has a checkpoint to resume from, since the export is roughly 4 GB;
6. confirms R and its packages are present now, rather than discovering they are not once Stage 03 is hours into the analysis.

**No cells are removed.** Thresholds are recomputed and reported, but nothing is filtered: `percent.mt` maxes at 9.9979 and `nFeature_RNA` bottoms at 202 — the fingerprints of the authors' filters having already run. Filtering again would diverge from the published dataset.

### Stage 03 — the clustering workflow

`scripts/03_reproduce_figure1d.R`, started automatically by `scripts/03_reproduce_figure1d.py`

The direct counterpart of the authors' steps 2 through 6:

```mermaid
flowchart TD
    R1["rebuild Seurat object<br/><code>min.cells = 0, min.features = 0</code><br/><i>no thresholds — the authors' already ran</i>"]
    R2["split into 59 per-sample objects<br/><i>assert none lost</i>"]
    R3["<b>SCTransform</b> per sample<br/>1,000 variable genes · ncells = 10,000<br/>regress the three quality measures<br/><code>vst.flavor = &quot;v1&quot;</code> set explicitly"]
    R4["select integration features<br/>request 3,000 → the full pool of ~2,816<br/>prepare SCT assays"]
    R5["<b>find anchors</b><br/>references A4 and A6, located by name"]
    R6["<b>integrate</b> into one assay<br/><i>assert cell count unchanged</i>"]
    R7["<b>PCA</b> · 30 components"]
    R8["shared nearest-neighbor graph<br/>PCs 1–15"]
    R9["<b>Louvain</b> · resolution = 0.3"]
    R10["<b>UMAP</b> · PCs 1–15, n.neighbors = 30"]
    R11(["write Seurat object,<br/>per-cell CSV, figure"])

    R1 --> R2 --> R3 --> R4 --> R5 --> R6 --> R7 --> R8 --> R9 --> R10 --> R11

    style R11 fill:#d4e8d4,stroke:#4a7a4a
```

The expensive steps are checkpointed, so `--restart-from` controls how far back a re-run goes:

| `--restart-from` | redoes | keeps | cost |
| --- | --- | --- | --- |
| `none` *(default)* | nothing | the most advanced checkpoint | — |
| `clustering` | neighbors, Louvain, UMAP | the PCA checkpoint | minutes |
| `pca` | PCA onward | the integrated object | ~30 min |
| `integration` | anchors, `IntegrateData` | the 59 SCT checkpoints | hours |
| `all` | everything, including per-sample SCTransform | nothing | many hours |

Each level discards its own checkpoint and every later one, so asking for a cheap re-run never triggers an expensive step.

### Result and the cluster-count difference

The result contains all 70,391 cells. The authors reported 13 clusters.

**An earlier version of this pipeline produced 11, and most of that gap now has an identified cause.** It requested 2,000 integration features, matching the authors' code as written. But `SelectIntegrationFeatures` keeps only genes measured in every one of the 59 samples, and the smallest sample holds 151 cells, so on Seurat 5.5.1 the candidate pool is roughly 2,816 genes rather than the tens of thousands the request implies. Asking for 2,000 discarded about 800 of them.

Which 800 is what matters. Candidates rank by how many samples each gene varies in, so the cut keeps genes that vary everywhere and drops genes that vary in only a few — and a plasma-cell or B-cell marker can only look variable in the samples that contained those cells. The markers of rare populations were removed before the merge, leaving those cells nothing to separate on:

| | clusters | smallest cluster |
| --- | --- | --- |
| authors | 13 | B cells, Undetermined |
| pipeline requesting 2,000 features | 11 | 671 cells |
| notebook run on the full pool | 12 | 120 cells |

Stage 03 now requests 3,000, which returns the pool entire. **The corrected pipeline has not yet been re-run end to end** — the 12-cluster row above comes from a notebook run over the identical matrix (`edited_good_reproduction.md`), not from this pipeline.

Two further differences were corrected at the same time. PCA computed 50 components where that reference run computed 30; the decomposition is iterative rather than exact, so the count perturbs the 15 components actually used. And `DimPlot` painted all 70,391 cells onto a 512 × 512 grid before stretching it across a 2,700 × 2,100 image, washing the figure out for reasons that have nothing to do with the analysis.

**Matching 13 exactly is still not expected.** The remaining gap is still attributed to Seurat implementation drift since October 2021 — a hypothesis that has not been tested, since testing it means re-running pinned to the authors' original package versions. Other explanations remain open: an undocumented parameter, or an unrecorded step between their preprocessing and their Figure 1D script.

Two decisions stand:

- **`vst.flavor = "v1"` is set deliberately.** The authors wrote their script before SCTransform v2 existed, so their unspecified call ran v1; modern Seurat defaults to v2. Removing this line would look like a closer match to their code while running a different algorithm. This narrows one known source of divergence — it does not resolve the count.
- **The resolution is not tuned to force 13 clusters.** Raising it until the number matched would report the published count obtained by a different method, burying the discrepancy rather than recording it.

Cluster count is treated as a diagnostic rather than a target. A run that returns no small clusters is evidence that rare cell types were lost during feature selection, and that is where to look.

## Future scGPT comparison

The next comparison replaces the PCA-based cell representation with a frozen scGPT embedding, keeping everything downstream as similar as possible.

```mermaid
flowchart TD
    A[("same AnnData object<br/>same cells, same metadata")]

    B["SCTransform integration"]
    C["<b>PCA</b><br/>PCs 1–15"]

    D["frozen <b>scGPT</b> embedding<br/>one vector per cell<br/><code>adata.obsm[&quot;X_scgpt&quot;]</code>"]

    E["nearest-neighbor graph"]
    F["nearest-neighbor graph"]
    G["clustering · same resolution target"]
    H["clustering · same resolution target"]
    I["UMAP · same n.neighbors"]
    J["UMAP · same n.neighbors"]

    K(["<b>compare on the same cells</b><br/>cluster stability · separation of broad populations<br/>donor mixing · marker-gene agreement<br/>recovery of rare and transitional populations"])

    A --> B --> C --> E --> G --> I --> K
    A --> D --> F --> H --> J --> K

    style C fill:#dbe8f5,stroke:#4a6a8a
    style D fill:#f0dbe8,stroke:#8a4a6a
    style K fill:#d4e8d4,stroke:#4a7a4a
```

Only the representation differs. The graph clustering, UMAP, and biological interpretation steps stay parallel.

## Future cluster labeling from gene expression

Cluster labels should be assigned from the original normalized RNA expression values — not from UMAP coordinates or embedding dimensions alone. The same workflow applies to both the PCA and scGPT clusters, so the two can be compared on the same biological evidence.

```mermaid
flowchart TD
    A["differential expression<br/>each cluster vs. all remaining cells"]
    B["rank positively enriched genes per cluster"]
    C["inspect marker heatmaps and expression plots"]
    D["assign broad labels<br/>from combinations of canonical markers"]
    E{"marker expression<br/>consistent across<br/>multiple donors?"}
    F["isolate broad T-cell and monocyte populations<br/>recluster for finer labels"]
    G(["labeled clusters"])
    H["keep labeled <b>ambiguous</b><br/><i>rather than forcing a cell type</i>"]

    A --> B --> C --> D --> E
    E -->|yes| F --> G
    E -->|no| H --> G

    style G fill:#d4e8d4,stroke:#4a7a4a
    style H fill:#f5e6d3,stroke:#a08050
```

## References

- Piehl N, et al. *Cerebrospinal fluid immune dysregulation during healthy brain aging and cognitive impairment.* Cell. 2022. DOI: 10.1016/j.cell.2022.11.019.
- Dataset: GEO accession **GSE200164**.
- Author code: `gatelabNW/csf_aging`.
