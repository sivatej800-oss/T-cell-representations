# GSE200164 Figure 1D reproduction

This repository contains a resumable three-stage pipeline that downloads the public
GSE200164 single-cell RNA-sequencing dataset, builds a canonical AnnData object, and
reproduces the unsupervised clustering and UMAP portion of Figure 1D from Piehl et al.
(2022).

If you are new to the project, start here. The files you run are the three numbered
Python scripts in `scripts/`. The notebooks, authors' code, and scGPT snapshot are
historical/reference material under `archive/`; they are not part of the working
pipeline.

## What the pipeline produces

The maintained endpoint is an unannotated Seurat UMAP of 70,391 cells from 59 CSF
samples. A successful run creates:

- `Age_trajectory/results/figure1d_python/integrated_clustered_umap.rds`: the complete
  integrated Seurat object with PCA, graph, clusters, and UMAP;
- `Age_trajectory/results/figure1d_python/GSE200164_umap_clusters.csv.gz`: a compact
  cell-level table of barcodes, UMAP coordinates, sample IDs, and clusters;
- `Age_trajectory/results/figure1d_python/figure1d_unsupervised_clusters.png`: the
  cluster UMAP that can be inspected without opening the RDS.

The maintained endpoint does **not** include the paper's manual cell-type labels,
T-cell/monocyte subclustering, TCR cleanup, or supplementary plots. Those analyses
depend partly on private metadata and unavailable helper code. The supplied author
workflow is retained at `archive/references/figure1d_author_workflow.R` to document
that boundary.

## Pipeline at a glance

Run these entry points in order from the repository root:

```bash
python3 scripts/01_download_audit_build_anndata.py
python3 scripts/02_prepare_figure1d.py
python3 scripts/03_reproduce_figure1d.py
```

```text
GEO downloads
    |
    v
01_download_audit_build_anndata.py
    |  canonical sparse AnnData + conversion audits
    v
02_prepare_figure1d.py
    |  QC fields + 59-sample mapping + optional Seurat export
    v
03_reproduce_figure1d.py
    |  R/Seurat integration or checkpoint resume
    v
RDS + cluster CSV + UMAP PNG
```

The scripts are safe to rerun. Each stage validates and reuses completed work. Do not
add `--force` unless you intentionally want to replace checkpoints and repeat expensive
steps.

## Requirements

### Software

- Python 3.10 or newer;
- the Python packages in `requirements.txt`;
- R and `Rscript` on `PATH`;
- R packages `Seurat`, `Matrix`, `glmGamPoi`, and `ggplot2`;
- internet access for a from-scratch Stage 01 download.

Create the Python environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

Install the R dependencies from an R session if they are not already available:

```r
install.packages(c("Seurat", "Matrix", "ggplot2"))
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("glmGamPoi")
```

Confirm both runtimes before starting:

```bash
python -c "import anndata, numpy, pandas, psutil, scipy; print('Python dependencies OK')"
Rscript -e 'for (p in c("Seurat", "Matrix", "glmGamPoi", "ggplot2")) stopifnot(requireNamespace(p, quietly=TRUE)); cat("R dependencies OK\n")'
```

### Storage and compute

This is a large single-cell workflow. The existing complete workspace occupies about
47 GB, including several 4-6 GB Seurat checkpoints. For a new run, use a filesystem
with substantially more free space than the final workspace and a machine with ample
RAM. The full SCTransform integration is the expensive part and can take hours;
runtime depends strongly on hardware and installed Seurat versions. Keep matrices
sparse—converting the full matrix to a dense array is impractical.

## How to run

### Option A: reuse existing local checkpoints

If generated files have been supplied under `Age_trajectory/` (as they are in the
current working copy), install the dependencies and run all three commands normally:

```bash
source .venv/bin/activate
python scripts/01_download_audit_build_anndata.py
python scripts/02_prepare_figure1d.py
python scripts/03_reproduce_figure1d.py
```

Expected behavior:

1. Stage 01 validates and reuses the canonical H5AD instead of downloading and
   rebuilding it.
2. Stage 02 recalculates the lightweight QC/sample audit and avoids a large Matrix
   Market export when a downstream R checkpoint is already available.
3. Stage 03 returns the completed result, or resumes from the post-PCA/integrated
   checkpoint instead of repeating SCTransform integration.

The final command should print the path to `integrated_clustered_umap.rds`. Repeating
Stage 03 with an intact final result should be a quick validation/no-op.

### Option B: run from scratch in a separate workspace

Use `--workspace` when you want to leave the existing results untouched:

```bash
python scripts/01_download_audit_build_anndata.py --workspace /path/to/new_workspace
python scripts/02_prepare_figure1d.py --workspace /path/to/new_workspace
python scripts/03_reproduce_figure1d.py --workspace /path/to/new_workspace
```

Use the same workspace path for all three commands. The scripts create its directory
structure automatically. Stage 01 downloads approximately 1.1 GB of compressed source
data, verifies the gzip stream, converts the matrix without discarding fractional
values, audits selected entries, and constructs the H5AD. Stages 02 and 03 then create
the Seurat inputs and checkpoints.

### Command-line options

Every stage accepts:

- `--workspace PATH`: choose the generated-data workspace; defaults to
  `Age_trajectory/` in this checkout;
- `--force`: replace that stage's generated results and recompute expensive work.

Stage 01 also accepts `--chunk-size`; Stages 02 and 03 accept `--rscript` for an
explicit Rscript path; Stage 03 accepts `--seed` (default `1234`). Use, for example:

```bash
python scripts/03_reproduce_figure1d.py --help
```

### Why there are entry-point files and implementation files

There is only one pipeline. The files have two different jobs:

- `scripts/01_*.py`, `scripts/02_*.py`, and `scripts/03_*.py` are the public commands.
  **These are the files a user runs.** They parse command-line options and call the
  appropriate implementation.
- `scripts/figure1d_pipeline/` is an internal Python package. **Do not run its Python
  files directly.** It keeps shared paths, validation, and reusable functions out of
  the command-line wrappers so logic is not duplicated.
- `run_seurat_figure1d.R` is an internal backend invoked automatically by Stage 03.
  A normal user does not call it directly.

The exact call mapping is:

| Command to run | Internal code it calls |
| --- | --- |
| `python scripts/01_download_audit_build_anndata.py` | `figure1d_pipeline/data_preparation.py` |
| `python scripts/02_prepare_figure1d.py` | `figure1d_pipeline/figure1d_preparation.py` and `r_environment.py` |
| `python scripts/03_reproduce_figure1d.py` | `figure1d_pipeline/seurat_runner.py`, which invokes `run_seurat_figure1d.R` |

`config.py` and `common.py` are shared by multiple stages; this is what prevents paths,
downloads, checksums, and file-validation logic from being repeated.

## What each stage does

| Stage | Reads | Work performed | Main output | Resume behavior |
| --- | --- | --- | --- | --- |
| `01` | GEO counts CSV and series matrix | Download validation, chunked sparse conversion, value audit, AnnData creation, clinical metadata integration | `data/processed/GSE200164_canonical_input.h5ad` | Validates and reuses an existing canonical H5AD |
| `02` | Canonical H5AD | Cell QC, count-like layer preservation, barcode-to-sample recovery, 59-sample clinical mapping, optional sparse Seurat export | `results/figure1d_python/GSE200164_after_qc_audit.h5ad` | Skips multi-GB Seurat export if an integrated/PCA checkpoint exists |
| `03` | Seurat export or reusable R checkpoint | SCTransform integration, PCA, neighbors, Louvain clustering, UMAP, plot/table export | `results/figure1d_python/integrated_clustered_umap.rds` | Prefers final result, then post-PCA, then integrated checkpoint; otherwise rebuilds |

All paths in the table are relative to the chosen workspace.

The R backend follows the authors' core settings: per-sample SCTransform with 1,000
variable features and the published QC regressors, A4/A6 integration references,
2,000 integration features, PCs 1-15 for neighbors and UMAP, clustering resolution
0.3, and 30 UMAP neighbors.

## Repository organization

```text
.
├── README.md                              # Start here: purpose, setup, run, outputs
├── requirements.txt                       # Active Python pipeline dependencies
├── scripts/                               # Maintained pipeline code
│   ├── 01_download_audit_build_anndata.py # RUN THIS: Stage 01 command
│   ├── 02_prepare_figure1d.py              # RUN THIS: Stage 02 command
│   ├── 03_reproduce_figure1d.py            # RUN THIS: Stage 03 command
│   ├── check_repository.py                 # Repository/notebook/file-size validator
│   ├── README.md                           # Developer-oriented module details
│   └── figure1d_pipeline/                  # INTERNAL: imported by the 3 commands
│       ├── config.py                       # Shared URLs and workspace paths
│       ├── common.py                       # Shared I/O and validation helpers
│       ├── data_preparation.py             # Called by the Stage 01 command
│       ├── figure1d_preparation.py         # Called by the Stage 02 command
│       ├── r_environment.py                # Called by Stages 02 and 03
│       ├── seurat_runner.py                # Called by the Stage 03 command
│       └── run_seurat_figure1d.R           # Invoked automatically by seurat_runner.py
├── Age_trajectory/                         # Generated/resumable workspace; not source
│   ├── README.md                           # Generated directory and retention guide
│   ├── data/                               # Downloads, sparse intermediates, metadata, H5AD
│   ├── logs/                               # Checksums and reproducibility summaries
│   ├── tables/                             # Conversion and QC audit tables
│   └── results/figure1d_python/             # Seurat inputs, checkpoints, and final outputs
├── archive/                                # Not executed by the maintained pipeline
│   ├── README.md                           # Archive inventory and status
│   ├── notebooks/                          # Superseded numbered notebooks/experiments
│   ├── references/                         # Authors' supplied non-runnable R workflow
│   ├── vendor/scGPT-main/                  # Unused upstream scGPT snapshot
│   └── legacy_data_directory/              # Obsolete root data-folder documentation
├── .github/                                # Active issue/PR templates and CI validation
├── .vscode/                                # Optional shared editor settings
├── CONTRIBUTING.md                         # Rules for changing active and archived code
└── SECURITY.md                             # Security and sensitive-data policy
```

The distinction is deliberate:

- run code only from `scripts/`;
- treat `Age_trajectory/` as generated, resumable state;
- treat `archive/` as read-only provenance, not a source of pipeline commands.

## Generated workspace contents

The default `Age_trajectory/` workspace is excluded from Git except for its README.
Its important paths are:

| Path | Created by | Meaning |
| --- | --- | --- |
| `data/raw/GSE200164_counts.csv.gz` | Stage 01 | Original GEO expression archive |
| `data/raw/GSE200164_series_matrix.txt.gz` | Stage 01 | GEO sample/clinical metadata header |
| `data/interim/*.npz` | Stage 01 | Sparse matrices in both required orientations |
| `data/interim/GSE200164_genes.csv` | Stage 01 | Matrix-aligned gene identifiers/annotations |
| `data/interim/GSE200164_barcodes.txt` | Stage 01 | Matrix-aligned cell barcodes |
| `data/metadata/*` | Stage 01 | Sample- and cell-level metadata exports |
| `data/processed/GSE200164_canonical_input.h5ad` | Stage 01 | Supported input to Stage 02 |
| `tables/*` and `logs/*` | Stages 01-02 | Conversion, QC, checksum, and run audits |
| `results/figure1d_python/seurat_input/` | Stage 02 | Raw sparse R boundary, only needed without R checkpoints |
| `results/figure1d_python/sct_samples/` | Stage 03 | Per-sample SCTransform resume points created by a fresh maintained run |
| `results/figure1d_python/GSE200164_after_PCA.rds` | Stage 03 | Preferred downstream resume point |
| `results/figure1d_python/integrated_clustered_umap.rds` | Stage 03 | Canonical final Seurat object |
| `results/figure1d_python/GSE200164_umap_clusters.csv.gz` | Stage 03 | Portable final cell/UMAP/cluster table |
| `results/figure1d_python/figure1d_unsupervised_clusters.png` | Stage 03 | Final visual output |

Temporary files ending in `.tmp` indicate interrupted atomic writes. Remove one only
after confirming that its completed sibling exists. See `Age_trajectory/README.md` for
the complete retention guide. Historical notebook-generated plots, tables, and unused
R checkpoints are preserved separately under
`Age_trajectory/archive/legacy_notebook_outputs/`.

## Expected scientific result

The current environment produced 11 unsupervised clusters and a 70,391-row cluster
CSV. The authors' historical environment reported 13 clusters. Differences in modern
Seurat/UMAP implementations can change the exact graph and cluster count even when the
published parameters are preserved, so exact cluster-number identity is not asserted.
The output is intentionally described as a pre-annotation reconstruction.

## Troubleshooting

- **`Rscript was not found on PATH`**: install R or pass
  `--rscript /absolute/path/to/Rscript` to Stages 02 and 03.
- **Missing R packages**: install the four packages listed under Requirements. The
  pipeline also writes `figure1d_environment_audit.json` to the results directory.
- **Missing raw Seurat input in Stage 03**: run Stage 02 using the same workspace.
- **Unexpected full recomputation**: check that checkpoint filenames remain inside
  `results/figure1d_python/` and that you did not pass `--force`.
- **Disk-space or memory failure**: move to a larger workspace with `--workspace`; do
  not densify the matrices or delete valid checkpoints needed for resume.
- **An interrupted `.tmp` RDS**: keep it until you know whether a valid completed RDS
  exists; Stage 03 uses completed checkpoints only.

## Validation

Before committing source or documentation changes, run:

```bash
python3 scripts/check_repository.py
git diff --check
```

The repository check validates notebook JSON, rejects accidentally tracked generated
artifacts, and enforces GitHub's practical per-file size boundary.

## Citation

If this repository contributes to published work, cite the underlying dataset and
methods paper:

> Piehl N, et al. Cerebrospinal fluid immune dysregulation during healthy brain aging
> and cognitive impairment. *Cell*. 2022.
> https://doi.org/10.1016/j.cell.2022.11.019

No license has been granted yet. Until a license is added, reuse requires permission
from the repository owner.
