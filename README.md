# START HERE

PCA vs. a frozen scGPT embedding on 70,391 CSF immune cells (GSE200164, Piehl et al.
2022), starting from a reproduction of the paper's Figure 1D clustering.

Go through the folders in this order:

| step | where | what you'll find |
| --- | --- | --- |
| 1 | [`01_PROJECT_OVERVIEW/STATUS.md`](01_PROJECT_OVERVIEW/STATUS.md) | Research question, what's done, headline results, what's next, scope. **Read this first.** |
| | [`01_PROJECT_OVERVIEW/ANALYSIS_WALKTHROUGH.md`](01_PROJECT_OVERVIEW/ANALYSIS_WALKTHROUGH.md) | The full scientific write-up of every stage: methods, diagrams, numbers, caveats. |
| 2 | [`02_CODE/`](02_CODE/) | The pipeline, numbered by stage. Run from the repository root. Each `.py` drives its matching `.R`. |
| 3 | [`03_RESULTS/`](03_RESULTS/) | Results by scientific question (`01_PCA_BASELINE` … `08_CELL_TYPE_ANNOTATION`), each with a README. |
| 4 | [`04_PAPER/`](04_PAPER/) | Outline and abstract. |
| 5 | [`05_REFERENCES/`](05_REFERENCES/) | The original paper, its published figure, and the authors' code. |

Other folders, which you normally don't need to open:

| where | what |
| --- | --- |
| `Age_trajectory/` | The pipeline workspace: data, checkpoints, and the canonical result files that `03_RESULTS/` links to. Generated, git-ignored, about 64 GB. Don't delete anything here without checking. |
| `ARCHIVE/`, `archive/` | Superseded work. `archive/` has to stay where it is because Stage 05 reads the scGPT gene vocabulary from `archive/vendor/`. |

## Running

```bash
python 02_CODE/01_download_audit_build_anndata.py
python 02_CODE/02_prepare_figure1d.py
python 02_CODE/03_reproduce_figure1d.py
# … later stages: see 03_RESULTS/README.md
```

Every stage resumes from saved work. Before committing, run
`python 02_CODE/check_repository.py`.
