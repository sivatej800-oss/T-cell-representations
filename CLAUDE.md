# Notes for Claude

A bioinformatics reproduction: rebuild Figure 1D of Piehl et al. 2022 (*Cell*) from
the public GSE200164 data. The scientific result is the deliverable; the code exists
to produce it. Read `README.md` first — it explains the biology and the analysis.

## Working here

Judge changes by whether they make the analysis more correct, more reproducible, or
easier to follow. This is not a deployed system and has no users but its author, so
speculative error handling, extra options, and abstraction layers are not wanted.

**Checks that protect scientific correctness are the exception and should stay.** The
pipeline verifies dataset dimensions at every boundary, confirms the sparse conversion
preserved values exactly, and compares one known patient's metadata against
hand-verified values. A mistake in single-cell data rarely crashes anything — it
produces a plausible-looking, wrong figure. That is the failure mode worth guarding.

Comments explain the data and the science, never Python or R syntax. The `.py` files
are organized as notebook-style blocks with a purpose header above each step.

## Facts that are easy to get wrong

**The GEO matrix is already filtered.** It is the *output* of the authors'
preprocessing, not the input. `percent.mt` maxes at 9.9979 and `nFeature_RNA` bottoms
at 202 — the fingerprints of their `percent.mt < 10` and `min.features = 200` filters
having already run. Doublets were removed before publication too. Do not add QC
filtering; Stage 02 measures and reports, deliberately removing nothing.

**Expression values are fractional on purpose.** The data is SoupX-corrected for
ambient RNA, which leaves non-integer values. Never round or cast to int — much
single-cell tooling assumes whole counts, and that assumption is wrong here. Stage 01
verifies this survived conversion.

**`vst.flavor = "v1"` in the R script is load-bearing.** The authors wrote their code
in October 2021, before v2 existed, so their unspecified call ran v1. Modern Seurat
defaults to v2. Removing this line would look like a closer match to their script but
would actually break the reproduction.

**11 clusters is the expected result, not a bug.** The paper reports 13. Every
published parameter is reproduced exactly (verified against
`figure1d_author_workflow.R` and `author_preprocessing/`), and the pipeline starts
from the same object their Figure 1D script starts from. The difference is four years
of Seurat implementation drift, mostly in SCTransform. Do not "fix" it by tuning
resolution.

**Stage 03 parameters that look like magic numbers are Seurat defaults** — 2,000
integration features, `dims = 1:30`, 50 PCs, `algorithm = 1`. The authors left them
unspecified; they are written out here for readability. The genuinely chosen values
are 1,000 variable features, resolution 0.3, and `dims = 1:15`.

**The three stages repeat their path and validation helpers on purpose.** They import
nothing from each other, so each file can be read and run on its own. This trades some
duplication for the ability to understand one stage without tracing across files.
Factoring the shared parts into a module would be the conventional move; it is not
wanted here.

## Environment

- Python: `.venv/` in this directory. Note `~/.venvs/jupyter` also exists and is the
  notebook kernel — not the pipeline environment.
- R: system install, packages under `~/R`. Seurat 5.5.1.
- Disk is tight (~38 GB free of 99 GB). The workspace is already ~47 GB, and a
  from-scratch run needs 60–80 GB free. Do not delete anything under
  `Age_trajectory/` without asking.

## Running it

```bash
python scripts/01_download_audit_build_anndata.py
python scripts/02_prepare_figure1d.py
python scripts/03_reproduce_figure1d.py
```

All three must share a workspace (default `Age_trajectory/`). Every stage resumes
from saved work, so re-running after an interruption is cheap.

To redo part of Stage 03, use `--restart-from {none,clustering,pca,integration,all}`
rather than `--force`. Re-tuning clustering costs minutes; `all` redoes 59
SCTransform runs and costs hours.

## Layout

- `scripts/` — the pipeline; the only code that runs. Stage 03's real work is in the
  `.R` file, started automatically by the `.py` file.
- `Age_trajectory/` — generated output, resumable, git-ignored. Not source.
- `figure1d_author_workflow.R` — the authors' Figure 1D code, reference only.
- `author_preprocessing/` — the authors' upstream steps, reference only. Confirms
  nothing analytical is missing from this pipeline.
- `archive/` — superseded notebooks and an unused scGPT snapshot. Never run.

## Before committing

```bash
python3 scripts/check_repository.py
```

Guards against generated data entering git.
