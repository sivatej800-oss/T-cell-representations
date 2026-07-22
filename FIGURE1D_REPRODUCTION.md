# Figure 1D reproduction audit

## Biological and computational target

Figure 1D in Piehl et al., *Cerebrospinal fluid immune dysregulation during
healthy brain aging and cognitive impairment* (Cell 185, 2022; DOI:
10.1016/j.cell.2022.11.019) is an all-cell UMAP colored by annotated CSF immune
cell type. The paper's deposited source is available in
`code/csf_aging-main/code/1_figures/Fig1D-F_Fig4H_SuppFig1C_1G_celltype_annotation.R`.

The released code and STAR Methods specify this analysis:

1. Remove ambient RNA with SoupX 1.5.2.
2. Remove doublets per sample with DoubletFinder 2.0.3, using 12 PCs, `pN =
   0.25`, an estimated `pK`, and an expected doublet rate of 1%.
3. Retain cells with at least 200 detected features, features detected in at
   least three cells, and cells with less than 10% mitochondrial reads.
4. Run `SCTransform` separately for each sample with 1,000 variable features,
   at most 10,000 sampled cells, and regression of `nCount_RNA`,
   `nFeature_RNA`, and `percent.mt`.
5. Integrate samples with SCT anchors, using healthy sample A4 and impaired
   sample A6 as reference samples to reduce memory use.
6. Run PCA on the integrated assay; use PCs 1–15 for a shared-nearest-neighbor
   graph, Louvain clustering at resolution 0.3, and UMAP with 30 neighbors.
7. Inspect canonical markers to assign broad T, monocyte, NK, dendritic, B,
   plasma, and undetermined clusters.
8. Recluster broad T-cell and monocyte compartments with the same PCA/graph/
   clustering/UMAP settings, then transfer the subtypes back to the global UMAP.

## Problems found in the previous local workflow

1. **Lossy count conversion.** `GSE200164_counts.csv` contains SoupX-adjusted
   fractional values (97%+ of sampled nonzero entries are non-integers). The
   previous Python converter used `dtype=np.int64`, truncating expression before
   Seurat analysis. The corrected cache is `float64` and includes a provenance
   manifest stating whether fractional values were observed.
2. **Wrong mitochondrial cutoff.** The local script retained cells at or below
   20% mitochondrial reads. The authors used a strict `<10%` threshold.
3. **Wrong normalization/integration workflow.** The local script applied one
   global `LogNormalize`/variable-feature/scale sequence and did not integrate
   samples. Figure 1D was built from per-sample SCTransform followed by SCT
   anchor integration. The corrected workflow implements that sequence.
4. **Missing hierarchical clustering.** The local script clustered all cells
   once. The authors separately reclustered broad T cells and monocytes before
   assigning the Figure 1D subtypes. The corrected workflow reproduces both
   compartment-specific analyses.
5. **Unreliable label scoring.** The previous annotation selected the largest
   unscaled mean from marker lists containing broadly shared genes such as
   `NKG7`, `LYZ`, `LST1`, and `CST3`. Those scores are not directly comparable
   across lineages and can systematically favor lists with abundant shared
   genes. The corrected workflow uses the authors' reproduced cluster mapping
   when the expected cluster solution is present and independently calculates
   standardized cluster-level canonical-marker scores as validation. It falls
   back to those scores only if software-version changes alter the cluster IDs.
6. **A requested cell type could never be produced.** `CD4+/CD8+ T Cells` was
   included in the legend but absent from the earlier marker-set classifier.
   Hierarchical T-cell reclustering now produces this subtype through the
   authors' intended mapping.
7. **Input path mismatch.** The previous script expected `data/ad_datasets`, but
   the uploaded files are in `data/AD datasets`. Both layouts are now resolved
   explicitly.

## Problems found in the authors' public Figure 1D script

These do not invalidate the paper's biological workflow, but they prevent the
released script from running reliably in a clean R session:

1. `color_vector` is used for Figure 1D before it is defined later in the file.
   The recreation defines the named palette before plotting.
2. The reclustering loop uses `cell_type` while the loop variable is named
   `celltype`. The recreation uses a tested helper function instead.
3. T-cell and monocyte cluster factors are assigned into character columns.
   Base R inserts the factors' one-based integer codes rather than labels such
   as `"0"`, shifting every cluster by +1. The authors noticed the shift in code
   comments and mapped the shifted values. The recreation transfers
   `as.character(seurat_clusters)` and converts their intended mappings back to
   the actual zero-based cluster IDs.
4. Cell-type composition divides by `sum(celltype_composition_pct)` while that
   object is still being created. The recreation divides counts by their real
   total.
5. The public label contains the typo `Monoctyes`; outputs use `Monocytes` while
   preserving the intended biological category and original color.

## What can and cannot be repeated from the deposited matrix

The GEO matrix contains 19,689 genes and 70,391 cells, exactly matching the sum
of the per-sample `n_cells` metadata. Its fractional values are consistent with
the paper's per-sample SoupX subtraction. The matrix was deposited after the
authors' empty-droplet, SoupX, DoubletFinder, and QC processing, so those steps
cannot be rerun from this combined table: removed droplets/doublets are no
longer present and raw Cell Ranger directories are required for SoupX. The
recreation therefore starts at the scientifically correct processed-count
entry point, verifies the published QC limits, and exactly/functionally repeats
the Figure 1D integration and annotation pipeline.

## Expected outputs

The full run writes:

- `results/figure1d/figure1d_recreated.pdf`
- `results/figure1d/figure1d_recreated.png`
- `results/figure1d/gse200164_figure1d_seurat.rds`
- `results/figure1d/figure1d_umap_coordinates.csv.gz`
- `results/figure1d/celltype_composition.csv`
- broad, T-cell, and monocyte marker-score and mapping-validation CSV files
- per-sample QC, exact run parameters, package versions, and diagnostic plots

Every plotted point is produced from the integrated expression data. No UMAP
coordinate, cluster, cell label, or visual element is placed manually.
