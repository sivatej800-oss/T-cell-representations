# Archived notebooks

This directory preserves notebook implementations that have been superseded by the
modular command-line pipeline under `../../scripts/`.

The notebooks are retained for provenance, saved outputs, and debugging history. They are not the supported way to rebuild the analysis because they contain duplicated cells, exploratory state, and machine-specific absolute paths.

The corresponding maintained entry points are:

| Archived notebook | Maintained entry point |
| --- | --- |
| `legacy_pipeline/data_preparation/01_download_audit_build_anndata.ipynb` | `../../scripts/01_download_audit_build_anndata.py` |
| `legacy_pipeline/figure1d/02_prepare_figure1d.ipynb` | `../../scripts/02_prepare_figure1d.py` |
| `legacy_pipeline/figure1d/03_reproduce_figure1d.ipynb` | `../../scripts/03_reproduce_figure1d.py` |

`experiments/figure1d_failed_attempt.ipynb` documents the earlier integration/checkpoint failure and its recovery. It has no maintained pipeline counterpart and should not be run as a canonical stage.

Do not delete the archived notebooks when changing the maintained pipeline; they document how earlier results were produced.
