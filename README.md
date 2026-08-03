# T Cell Representations

Research workflows for constructing and evaluating single-cell T-cell representations from cerebrospinal-fluid (CSF) transcriptomic data. The current work focuses on the public [GSE200164](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200164) dataset and the Figure 1D analysis from Piehl et al. (2022).

## Project status

This is an active research repository. The notebooks capture data auditing, AnnData construction, Seurat-based Figure 1D reproduction, and early representation-learning experiments. Generated matrices and model checkpoints are intentionally not versioned.

## Repository layout

```text
.
├── notebooks/
│   ├── data_preparation/   # Download, audit, and AnnData construction
│   ├── figure1d/           # Figure 1D preparation and reproduction
│   ├── experiments/        # Superseded or diagnostic experiments
│   └── utilities/          # Data-transfer and integrity checks
├── Age_trajectory/         # Local generated data, results, and QC artifacts
├── code/scGPT-main/        # Vendored scGPT research code used by experiments
├── scripts/                # Repository integrity checks
├── data/                   # Local-only data area; see data/README.md
└── requirements.txt        # Core Python notebook dependencies
```

## Notebook workflow

Run the notebooks in this order when rebuilding the analysis:

1. `notebooks/data_preparation/01_download_audit_build_anndata.ipynb`
2. `notebooks/figure1d/02_prepare_figure1d.ipynb`
3. `notebooks/figure1d/03_reproduce_figure1d.ipynb`
4. `notebooks/utilities/validate_canonical_data_transfer.ipynb` when moving a canonical AnnData file between systems

The notebook under `notebooks/experiments/` is retained as a documented failed attempt and is not part of the canonical workflow.

## Setup

Python 3.10 or newer is recommended.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python -m ipykernel install --user --name t-cell-representations
```

The Figure 1D reproduction also requires R with `Seurat`, `Matrix`, `glmGamPoi`, `future`, `dplyr`, `ggplot2`, `patchwork`, and `readr`.

Several notebooks currently define an absolute project path near the top. Update that single configuration cell for your checkout before running the workflow. Keep the full expression matrix sparse; converting the complete matrix to a dense array requires impractical memory.

## Data

The analysis uses GSE200164, comprising 70,391 cells and 19,689 genes across 59 samples in the processed matrix. Raw downloads, sparse matrices, AnnData objects, Seurat objects, and checkpoints are excluded from Git. See [data/README.md](data/README.md) for the expected local data policy.

## Validation

Run the repository check before committing:

```bash
python scripts/check_repository.py
```

It verifies notebook JSON, checks for accidentally tracked generated artifacts, and enforces GitHub's practical per-file size boundary.

## Contributing and security

See [CONTRIBUTING.md](CONTRIBUTING.md) for the research workflow and [SECURITY.md](SECURITY.md) for responsible disclosure and data-handling guidance.

## Citation

If this repository contributes to published work, cite the underlying dataset and methods paper:

> Piehl N, et al. Cerebrospinal fluid immune dysregulation during healthy brain aging and cognitive impairment. *Cell*. 2022. https://doi.org/10.1016/j.cell.2022.11.019

No license has been granted yet. Until a license is added, reuse requires permission from the repository owner.
