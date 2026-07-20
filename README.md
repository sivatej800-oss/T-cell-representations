# Summer Project 2026

This repo expects large research data files to be uploaded directly into local ignored folders:

- `data/ad_datasets/`
- `data/download_family/`

The large data folders are intentionally ignored by git.

## Folder layout

- `code/` contains source code, scripts, notebooks, tests, and the vendored `scGPT-main` project.
- `data/` contains local research datasets and downloaded source files.

## Expected files

Put these files in `data/ad_datasets/`:

- `GSE200164_metadata.csv`
- `GSE200164_counts.csv`

Put these files in `data/download_family/`:

- `GSE200164_series_matrix.txt`
- `GSE200164_family.soft`

`GSE200164_lognorm_expression.csv` is intentionally excluded from the local dataset upload.
