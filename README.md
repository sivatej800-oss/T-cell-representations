# Summer Project 2026

This repository contains a reproducible recreation of Figure 1D from Piehl et
al., *Cell* 2022, using the deposited GSE200164 CSF single-cell expression data.

## Figure 1D recreation

The corrected workflow is implemented in:

- `code/test_code/prepare_gse200164_sparse.py` — losslessly converts the wide
  SoupX-adjusted GEO table to a reusable sparse matrix without rounding its
  fractional values.
- `code/test_code/generate_fig1d_umap.R` — applies the authors' per-sample SCT
  integration, PCA, Louvain clustering, UMAP, hierarchical T-cell/monocyte
  reclustering, marker validation, and plotting workflow.
- `FIGURE1D_REPRODUCTION.md` — records the audit findings, methods mapping,
  expected outputs, and limitations of the deposited processed matrix.

Run the complete 70,391-cell analysis from the repository root:

```bash
Rscript code/test_code/generate_fig1d_umap.R
```

The first run creates a sparse cache and analysis checkpoints under
`results/figure1d/`. Later runs reuse compatible checkpoints. To test the
workflow quickly while retaining all 59 samples, use a separate output folder:

```bash
Rscript code/test_code/generate_fig1d_umap.R \
  --max-cells-per-sample=50 \
  --sparse-dir=results/figure1d/sparse_cache \
  --out-dir=results/figure1d_smoke
```

Use `--rebuild-sparse` only when the source matrix changes. Use
`--rebuild-analysis` when analysis parameters change and existing checkpoints
must be replaced.

## Local data layout

The scripts recognize either of these local layouts:

- `data/AD datasets/GSE200164_counts.csv` and
  `data/AD datasets/GSE200164_metadata.csv`
- `data/ad_datasets/GSE200164_counts.csv` and
  `data/ad_datasets/GSE200164_metadata.csv`

Large data and generated results are intentionally ignored by Git.

The authors' uploaded reference repository is available under
`code/csf_aging-main/`; the recreation does not modify it.
