# Pipeline scripts

The numbered Python files in this directory are the supported end-to-end workflow. Run them in numeric order from the repository root.

```bash
python3 scripts/01_download_audit_build_anndata.py
python3 scripts/02_prepare_figure1d.py
python3 scripts/03_reproduce_figure1d.py
```

Each entry point accepts `--workspace PATH`; the default is the repository's `Age_trajectory/` directory. Run any stage with `--help` for its complete options. `--force` intentionally recomputes or replaces that stage's generated outputs and should normally be unnecessary.

## Package layout

The numbered scripts are the public commands and are the only pipeline files users
should execute directly. They are intentionally thin: each parses its command-line
arguments and delegates to reusable code in `figure1d_pipeline/`.

`figure1d_pipeline/` is an internal package, not a second pipeline and not another set
of commands:

- `config.py` defines URLs and every workspace-relative path.
- `common.py` provides shared downloads, checksums, validation, and atomic writes.
- `data_preparation.py` implements Stage 01.
- `figure1d_preparation.py` implements Stage 02.
- `r_environment.py` audits the R packages used by Stages 02 and 03.
- `seurat_runner.py` is the Python-to-R boundary for Stage 03.
- `run_seurat_figure1d.R` contains the Seurat computation and checkpoint-resume logic.

Stage 03 invokes the R file automatically through `seurat_runner.py`; users do not
need to run the R file themselves.

## Checkpoint behavior

- Stage 01 reopens and validates an existing canonical H5AD instead of repeating the download and sparse conversion.
- Stage 02 recalculates QC and sample mapping, but skips the multi-gigabyte Matrix Market export when an integrated or post-PCA Seurat checkpoint is already available.
- Stage 03 prefers the completed result, then the post-PCA checkpoint, then the integrated checkpoint. It only rebuilds from Matrix Market input when no usable Seurat checkpoint exists.
- Per-sample SCTransform objects are checkpointed when a full rebuild is required.
- The final RDS is written to a temporary sibling file and renamed only after a successful save.

## Main outputs

Stage 03 writes these files under `Age_trajectory/results/figure1d_python/`:

- `integrated_clustered_umap.rds`: complete integrated Seurat object.
- `GSE200164_umap_clusters.csv.gz`: cell barcodes, UMAP coordinates, and cluster IDs.
- `figure1d_unsupervised_clusters.png`: unsupervised cluster UMAP.

The plot is a pre-annotation reconstruction. Historical manual cell-type annotation and T-cell/monocyte subclustering from the authors' analysis are intentionally outside the required pipeline.

The authors' supplied workflow is archived separately at
`archive/references/figure1d_author_workflow.R`; reference material does not live
among executable scripts.
