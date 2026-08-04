# GSE200164 Figure 1D reproduction

This project reproduces one figure from a published study of the immune cells found
in cerebrospinal fluid — the liquid surrounding the brain and spinal cord — during
healthy ageing and cognitive decline.

The study measured, one cell at a time, which genes were active in 70,391 immune
cells drawn from 59 people. Some were cognitively healthy, others had mild cognitive
impairment or Alzheimer's disease. Figure 1D of that paper is the overview: every one
of those cells drawn as a single point on a two-dimensional map, positioned so that
cells behaving similarly sit near each other, and coloured by the group they were
sorted into. That map is what this pipeline rebuilds from the public data.

## What the result looks like

![Unsupervised clusters](Age_trajectory/results/figure1d_python/figure1d_unsupervised_clusters.png)

Every point is one cell. The map separates into distinct territories, and those are
the different kinds of immune cell present — the large landmasses on the right are
the most abundant populations, while small isolated islands such as clusters 8, 9,
and 10 are rare cell types making up a few hundred cells each.

This run found **11 clusters** across the 70,391 cells. They are very unevenly sized,
which is expected of immune populations: the largest holds 22,660 cells and the
smallest 671.

**The clusters here are numbered, not named.** The published figure labels each
territory with the cell type it represents — T cells, monocytes, and so on. Assigning
those names is a separate interpretive step that this pipeline stops short of, for
reasons described under [Where this reproduction stops](#where-this-reproduction-stops).

## How faithful is this to the paper?

The paper reported 13 clusters; this reproduction finds 11.

That difference is not a mistake in the settings. Every published parameter is
reproduced exactly — the same normalization method, the same 2,000 genes used to
align samples, the same two reference patients, the same 15 dimensions, the same
clustering resolution. What differs is the software: the clustering and layout
algorithms have been revised across Seurat versions since the paper was written, and
those revisions change where the boundaries fall. Two territories the authors' version
separated may be joined here.

So this should be read as a faithful reconstruction of the paper's *method*, not a
claim to have regenerated its exact output. The overall structure — the number of
major populations, their relative sizes, their arrangement — is what carries over.
The exact cluster count does not, and is not asserted.

## What the pipeline does

The work runs in three stages, each a single file you run in order. Each stage saves
its result to disk, and the next stage reads that file. They never call one another.

```text
    The two public files published with the study
                        │
   ┌────────────────────▼────────────────────┐
   │ Stage 01  download, check, and assemble │
   │                                          │
   │ Fetches the data and proves it arrived   │
   │ intact. Converts a very large table into │
   │ a compact form, then verifies nothing    │
   │ changed in the conversion. Attaches each │
   │ patient's age, sex, and diagnosis to     │
   │ every one of their cells.                │
   └────────────────────┬────────────────────┘
                        │  70,391 cells, each knowing whose it is
   ┌────────────────────▼────────────────────┐
   │ Stage 02  measure cell quality           │
   │                                          │
   │ Calculates three measurements per cell   │
   │ describing how well it was captured, and │
   │ matches the two systems of patient names │
   │ used by the repository and the paper.    │
   │ Writes the data out in a form R reads.   │
   └────────────────────┬────────────────────┘
                        │  the same cells, now described
   ┌────────────────────▼────────────────────┐
   │ Stage 03  the analysis itself (in R)     │
   │                                          │
   │ Normalizes each patient's cells, merges  │
   │ all 59 samples so they can be compared,  │
   │ groups the cells into clusters, and lays │
   │ them out on the two-dimensional map.     │
   └────────────────────┬────────────────────┘
                        │
        the map above, a cluster table, and the full result
```

Stage 03 is written in R rather than Python, because the published method uses
Seurat, an R library. You never run the R file yourself — the third Python file
starts it with the right settings.

### Why the analysis works the way it does

The central problem is that the 59 samples cannot simply be pooled. Each was
collected and processed separately, and those differences in handling leave a mark on
the measurements. Pooled naively, cells would group by which patient they came from
rather than by what type of cell they are, which is the opposite of what the figure
is meant to show.

The analysis therefore proceeds in stages that each remove a source of unwanted
variation:

1. **Normalize within each patient.** Raw measurements are not comparable between
   cells, because a cell that yielded more material shows higher numbers for every
   gene regardless of its biology. This step corrects for that, and simultaneously
   removes three known technical influences — how much RNA was captured, how many
   genes were detected, and how much of the signal came from mitochondria, which is
   a signature of a dying cell.

2. **Merge the samples.** Pairs of cells in different patients that closely resemble
   each other are taken to be the same kind of cell. Those pairings reveal how each
   sample is systematically shifted relative to the others, which is then corrected.
   Two patients act as the fixed reference everything else is aligned to.

3. **Condense to the main patterns.** Each cell is described by thousands of gene
   measurements, most of which rise and fall together. These are reduced to a smaller
   number of combined dimensions, ordered strongest first, on the reasoning that the
   leading ones carry the biology while the rest carry noise.

4. **Group the cells.** Each cell is linked to the cells it most resembles, and that
   network is divided into clusters where connections are densest.

5. **Lay them out.** Finally the cells are arranged on a flat map for drawing, keeping
   similar cells near each other.

A note on reading the map: it is a way of *viewing* structure, not a measurement.
Distances across it should not be treated as quantities, and the clustering in step 4
does not depend on it.

## Running it

### What you need

- Python 3.10 or newer, and the packages in `requirements.txt`
- R, with `Seurat`, `Matrix`, `glmGamPoi`, and `ggplot2`
- Internet access, if starting from scratch

Set up Python:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

Set up R, from an R session:

```r
install.packages(c("Seurat", "Matrix", "ggplot2"))
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("glmGamPoi")
```

Confirm both are ready:

```bash
python -c "import anndata, numpy, pandas, psutil, scipy; print('Python dependencies OK')"
Rscript -e 'for (p in c("Seurat", "Matrix", "glmGamPoi", "ggplot2")) stopifnot(requireNamespace(p, quietly=TRUE)); cat("R dependencies OK\n")'
```

### The three commands

```bash
source .venv/bin/activate
python scripts/01_download_audit_build_anndata.py
python scripts/02_prepare_figure1d.py
python scripts/03_reproduce_figure1d.py
```

All three must use the same workspace. With no options they all use
`Age_trajectory/` inside this checkout. To work somewhere else, pass the same
`--workspace PATH` to each of the three.

Every command is safe to run again. Each checks for finished work before starting,
so re-running after an interruption resumes rather than repeating. With results
already present, all three finish in seconds.

### What it costs to run from scratch

The figures below are from an actual complete run on a machine with 63 GB of RAM.
Yours will differ, but the proportions are informative — nearly all the time is spent
in the two steps that normalize and merge the samples.

| | |
| --- | --- |
| Download | about 1.1 GB compressed |
| Finished workspace | about 47 GB |
| Free disk needed | **60–80 GB**, more than the final size, because intermediate files coexist |
| Stage 01 | roughly an hour, dominated by converting and verifying the large table |
| Stage 02 | minutes, plus roughly 15 minutes if it writes the 4 GB hand-off file |
| Stage 03 | **several hours** — normalizing and merging the 59 samples is the bulk of it |

Some individual saved files are 4–4.5 GB. Keep the data sparse; converting the full
matrix to a dense one is not feasible on ordinary hardware.

### Redoing part of the analysis

Because Stage 03 takes hours, it saves progress along the way and can restart from
any point rather than from the beginning:

```bash
python scripts/03_reproduce_figure1d.py --restart-from clustering
```

| Level | Redoes | Roughly |
| --- | --- | --- |
| `none` | nothing; resumes from the furthest saved point | seconds |
| `clustering` | the grouping and the map | minutes |
| `pca` | condensing to main patterns, onward | ~30 minutes |
| `integration` | merging the samples, onward | hours |
| `all` | everything, including per-patient normalization | many hours |

This matters most when adjusting the settings that shape the figure — the clustering
resolution or the number of dimensions. Those live in `scripts/03_reproduce_figure1d.R`,
and changing them only requires `--restart-from clustering`, not a full re-run.

Other options: `--workspace PATH` on all three stages, `--chunk-size` on Stage 01,
`--rscript` on Stages 02 and 03 when R is not on your `PATH`, and `--seed` on Stage 03.
`--force` still works everywhere and means "redo this stage's work".

## What you get

Everything lands in `Age_trajectory/results/figure1d_python/`:

| File | What it is |
| --- | --- |
| `figure1d_unsupervised_clusters.png` | The figure shown above |
| `GSE200164_umap_clusters.csv.gz` | One row per cell: its map position, cluster, and patient. A few megabytes, opens anywhere — **start here if you want to work with the result** |
| `integrated_clustered_umap.rds` | The complete analysis, about 4.5 GB. What any further analysis would build on |
| `GSE200164_sample_metadata_mapping.csv` | The 59 patients: both their names, age, sex, diagnosis, and cell count |
| `figure1d_environment_audit.json` | Which software versions produced this result |

Files ending in `.tmp` are interrupted writes. Delete one only after confirming its
finished counterpart exists.

## Where this reproduction stops

This pipeline produces the map and the clusters. It does not go on to name the
clusters by cell type, subdivide the T cells and monocytes, clean up the T-cell
receptor data, or generate the supplementary panels.

That boundary is not arbitrary. Those analyses depend partly on metadata that was
never made public and on helper code that was not distributed with the paper. Rather
than guess at them, the pipeline stops at the last point that can be reproduced
honestly from public data.

The authors' original workflow is kept at `figure1d_author_workflow.R` in the project
root. It is not runnable here and is not part of the pipeline — it is reference
material, showing what the published analysis did and therefore exactly where this
reproduction parts company with it.

## How the project is organized

```text
.
├── README.md                               # This file
├── requirements.txt                        # Python packages the pipeline needs
├── figure1d_author_workflow.R              # The authors' original code, for reference only
├── scripts/                                # The pipeline — the only code you run
│   ├── 01_download_audit_build_anndata.py  # Stage 01, complete and standalone
│   ├── 02_prepare_figure1d.py              # Stage 02, complete and standalone
│   ├── 03_reproduce_figure1d.py            # Stage 03, starts the R analysis
│   ├── 03_reproduce_figure1d.R             # The analysis itself
│   ├── check_repository.py                 # Guards against committing large data
│   └── README.md                           # Step-by-step detail on each stage
├── Age_trajectory/                         # Everything the pipeline generates
│   ├── data/                               # Downloads and intermediate files
│   ├── logs/                               # Checksums and run summaries
│   ├── tables/                             # Verification results
│   └── results/figure1d_python/            # The outputs listed above
├── archive/                                # Kept for reference, never run
│   ├── notebooks/                          # Earlier notebook versions of this work
│   ├── vendor/scGPT-main/                  # An unused library snapshot
│   └── legacy_data_directory/              # Notes on a retired folder layout
├── CONTRIBUTING.md
└── SECURITY.md
```

Three rules keep this straight: run only what is in `scripts/`, treat
`Age_trajectory/` as regenerable output rather than source, and read `archive/` as
history rather than instructions.

`Age_trajectory/` is excluded from version control apart from its README, because
single-cell data is far too large for Git.

## Notes for anyone extending this

A few decisions are worth knowing before changing anything:

**The measurements are deliberately not whole numbers.** The published data was
corrected for ambient RNA — free-floating genetic material that contaminates droplets
— and that correction leaves fractional values behind. Rounding them to integers, as
much single-cell software expects, would destroy the correction. The pipeline
preserves them at every step, and Stage 01 verifies this after converting the data.

**No cells are filtered out.** Stage 02 calculates the paper's quality thresholds and
reports how many cells fall outside them, but removes nothing. In practice none do,
because the authors published an already-filtered dataset. Filtering again would
silently diverge from what was published.

**Nothing is trusted without checking.** The dataset's dimensions are verified at
every boundary, the conversion in Stage 01 is checked against the original download,
and one patient's details are compared against values confirmed by hand. This is
deliberate: a mistake in single-cell data does not usually crash anything, it just
produces a plausible-looking and wrong figure.

## Validation

Before committing changes:

```bash
python3 scripts/check_repository.py
git diff --check
```

This confirms no generated data has crept into version control and no file exceeds
what GitHub accepts.

## Citation

If this work contributes to something published, cite the original study:

> Piehl N, et al. Cerebrospinal fluid immune dysregulation during healthy brain aging
> and cognitive impairment. *Cell*. 2022.
> https://doi.org/10.1016/j.cell.2022.11.019

The data is public under accession
[GSE200164](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200164).

No license has been granted for this repository yet. Until one is added, reuse
requires permission from the owner.
