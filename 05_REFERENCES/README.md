# References

Reference material only. None of it is run, and no pipeline script reads from here.

| where | what |
| --- | --- |
| `paper/Original_paper.pdf` | Piehl N, et al. *Cerebrospinal fluid immune dysregulation during healthy brain aging and cognitive impairment.* Cell. 2022. DOI: 10.1016/j.cell.2022.11.019 |
| `paper/Author's figure.png` | The published Figure 1D. |
| `author_code/figure1d_author_workflow.R` | The authors' integration, clustering, manual annotation, and reclustering code. This is the one to read. |
| `author_code/00_…` – `06_…` | The authors' upstream preprocessing: Cell Ranger, SoupX, TCR, DoubletFinder, pseudobulk (from `gatelabNW/csf_aging`). |

`figure1d_author_workflow.R` is not an executable entry point. It has placeholder
paths, depends on unpublished helpers and metadata fields, and includes downstream
analyses outside this project. Its core parameters are already built into
`02_CODE/03_reproduce_figure1d.R`. The authors' upstream preprocessing is already
applied to the public GEO matrix.

The scGPT source code, which is also reference material, stays at
`archive/vendor/scGPT-main/` because Stage 05 reads its gene vocabulary from that path.
