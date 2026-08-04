# Sequential pipeline files

The pipeline is intentionally flat. Run the three Python files in order:

```bash
python3 scripts/01_download_audit_build_anndata.py
python3 scripts/02_prepare_figure1d.py
python3 scripts/03_reproduce_figure1d.py
```

Each Python file contains its own arguments, paths, validation helpers, and stage logic.
The Python files do not import one another, and there is no internal pipeline package.

The stages connect only through saved files:

```text
01_download_audit_build_anndata.py
    |
    | writes
    v
Age_trajectory/data/processed/GSE200164_canonical_input.h5ad
    |
    | read by
    v
02_prepare_figure1d.py
    |
    | writes
    v
Age_trajectory/results/figure1d_python/GSE200164_after_qc_audit.h5ad
Age_trajectory/results/figure1d_python/seurat_input/
    |
    | read by
    v
03_reproduce_figure1d.py
    |
    | starts the adjacent 03_reproduce_figure1d.R because Seurat is an R library
    v
integrated_clustered_umap.rds + cluster CSV + UMAP PNG
```

## File 01: data preparation

`01_download_audit_build_anndata.py` contains the complete Python logic for:

1. creating the workspace directories;
2. reusing a valid canonical H5AD when one already exists;
3. downloading and checking the two GEO source files;
4. streaming the large expression CSV into sparse matrices;
5. auditing selected source and sparse-matrix values;
6. creating AnnData;
7. attaching sample metadata; and
8. saving `GSE200164_canonical_input.h5ad`.

The `main()` function is at the bottom of the file. It reads arguments and calls
`run_data_preparation()`, which lists the steps above in execution order.

## File 02: Figure 1D preparation

`02_prepare_figure1d.py` contains the complete Python logic for:

1. opening the canonical H5AD written by File 01;
2. confirming the expected 70,391 × 19,689 sparse matrix;
3. calculating `nCount_RNA`, `nFeature_RNA`, and `percent.mt`;
4. mapping GEO accessions to the paper's A1-H8 sample IDs;
5. saving its QC H5AD checkpoint;
6. exporting sparse Seurat files when no later checkpoint exists; and
7. checking R and its required packages.

File 02 never runs File 01. If File 01's H5AD does not exist in the selected workspace,
File 02 stops with a missing-file error.

## File 03: Seurat reproduction

`03_reproduce_figure1d.py` contains all Python checkpoint and environment logic. It:

1. looks for a completed final RDS;
2. otherwise looks for a post-PCA or integrated RDS checkpoint;
3. otherwise requires the sparse files written by File 02;
4. verifies R and its packages;
5. starts `03_reproduce_figure1d.R`; and
6. verifies the final RDS exists before reporting success.

`03_reproduce_figure1d.R` is the only companion implementation file. It is necessary
because the published method uses Seurat, which is an R package. Users still run only
the Python file; File 03 supplies the R input directory, output directory, seed, and
force setting automatically.

## Shared workspace rule

Use the same workspace for all three files. With no option, all use
`Age_trajectory/`. For another location:

```bash
python3 scripts/01_download_audit_build_anndata.py --workspace /path/to/workspace
python3 scripts/02_prepare_figure1d.py --workspace /path/to/workspace
python3 scripts/03_reproduce_figure1d.py --workspace /path/to/workspace
```

The files repeat their small path and validation helpers deliberately. This keeps each
stage readable and runnable by itself without hidden cross-file Python imports.

## Resume behavior

- File 01 reuses a valid canonical H5AD.
- File 02 reuses its H5AD and skips an unnecessary multi-gigabyte export.
- File 03 uses the most advanced completed R checkpoint.

Do not use `--force` unless you intend to overwrite checkpoints and repeat expensive
work.
