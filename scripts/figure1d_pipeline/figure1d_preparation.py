"""Implementation of Stage 02: prepare canonical AnnData for Seurat.

This stage validates Stage 01, calculates the three QC fields expected by the authors'
Seurat workflow, recovers the internal sample IDs, and exports a sparse R interchange
format only when Stage 03 cannot resume from an existing checkpoint.
"""

from __future__ import annotations

import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
from scipy import sparse

from .common import atomic_h5ad_write, require_files
from .config import PipelinePaths
from .data_preparation import EXPECTED_CELLS, EXPECTED_GENES, EXPECTED_SAMPLES
from .r_environment import audit_r_environment


def load_and_audit_canonical(paths: PipelinePaths) -> ad.AnnData:
    """Load Stage 01 and stop early if its shape, sparsity, or names are unexpected."""

    require_files([paths.canonical_h5ad])
    adata = ad.read_h5ad(paths.canonical_h5ad)
    if adata.shape != (EXPECTED_CELLS, EXPECTED_GENES):
        raise ValueError(
            f"Canonical AnnData has shape {adata.shape}; "
            f"expected {(EXPECTED_CELLS, EXPECTED_GENES)}"
        )
    if not sparse.issparse(adata.X):
        raise TypeError("Canonical AnnData.X must remain sparse")
    # CSR is efficient for the per-cell row calculations performed below.
    if not sparse.isspmatrix_csr(adata.X):
        adata.X = adata.X.tocsr()
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ValueError("Canonical cell and gene names must be unique")

    # Sparse matrices keep only nonzero values in ``data``. Inspecting that array is
    # much cheaper than expanding the complete cells-by-genes matrix.
    values = adata.X.data
    negative_values = int(np.count_nonzero(values < 0))
    integer_fraction = float(
        np.isclose(values, np.round(values), rtol=0, atol=1e-8).mean()
    )
    if negative_values:
        raise ValueError(f"Expression matrix contains {negative_values:,} negative values")
    print(f"Loaded canonical AnnData: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    print(f"Sparse nonzeros: {adata.X.nnz:,}; integer-like fraction: {integer_fraction:.6f}")
    if integer_fraction < 0.99:
        print(
            "WARNING: the matrix is predominantly fractional. The pipeline preserves "
            "the dataset's SoupX-adjusted count-like interpretation."
        )
    return adata


def add_figure1d_qc(adata: ad.AnnData, *, sample_key: str = "gsm") -> None:
    """Preserve count-like data and calculate Seurat-compatible cell QC fields.

    ``nCount_RNA`` is total expression per cell, ``nFeature_RNA`` is the number of
    detected genes, and ``percent.mt`` is the mitochondrial-expression percentage.
    """

    if sample_key not in adata.obs:
        raise KeyError(f"Sample key {sample_key!r} is missing from AnnData.obs")
    if adata.obs[sample_key].isna().any():
        raise ValueError(f"Sample key {sample_key!r} contains missing values")
    if adata.obs[sample_key].nunique() != EXPECTED_SAMPLES:
        raise ValueError(
            f"Expected {EXPECTED_SAMPLES} samples in {sample_key!r}, "
            f"found {adata.obs[sample_key].nunique()}"
        )
    adata.obs[sample_key] = adata.obs[sample_key].astype("category")

    # AnnData layers are alternate versions of the same matrix. Keeping a named counts
    # layer prevents later transformations from overwriting the Stage 03 input.
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
    counts = adata.layers["counts"]
    if not sparse.issparse(counts) or counts.shape != adata.shape:
        raise TypeError("AnnData counts layer must be sparse and match AnnData.X")
    if not sparse.isspmatrix_csr(counts):
        counts = counts.tocsr()
        adata.layers["counts"] = counts

    # axis=1 reduces across genes, producing one QC value for each cell.
    adata.obs["nCount_RNA"] = np.asarray(counts.sum(axis=1)).ravel()
    adata.obs["nFeature_RNA"] = np.diff(counts.indptr)
    mitochondrial = adata.var_names.astype(str).str.startswith("MT-")
    if int(mitochondrial.sum()) == 0:
        raise ValueError("No MT- genes were found; percent.mt cannot be calculated")
    mitochondrial_counts = np.asarray(counts[:, mitochondrial].sum(axis=1)).ravel()
    totals = adata.obs["nCount_RNA"].to_numpy()
    # np.divide's ``where`` avoids division-by-zero warnings for an empty cell.
    adata.obs["percent.mt"] = np.divide(
        mitochondrial_counts,
        totals,
        out=np.zeros(adata.n_obs, dtype=np.float64),
        where=totals > 0,
    ) * 100

    below_200 = int((adata.obs["nFeature_RNA"] < 200).sum())
    high_mitochondrial = int((adata.obs["percent.mt"] >= 10).sum())
    print(adata.obs[["nCount_RNA", "nFeature_RNA", "percent.mt"]].describe().T)
    print(f"Cells below 200 features: {below_200:,}")
    print(f"Cells at or above 10% mitochondrial expression: {high_mitochondrial:,}")
    if below_200 or high_mitochondrial:
        print("WARNING: some cells do not meet the paper's measurable QC thresholds")

    # Record the scientific parameters beside the data so the checkpoint is
    # self-describing when opened outside this repository.
    adata.uns["figure1d_reproduction"] = {
        "dataset": "GSE200164",
        "sample_key": sample_key,
        "expected_samples": EXPECTED_SAMPLES,
        "expected_cells": EXPECTED_CELLS,
        "expected_genes": EXPECTED_GENES,
        "input_layer": "counts",
        "input_interpretation": "Processed SoupX-adjusted quality-controlled count-like matrix",
        "authors_initial_pcs": 15,
        "authors_umap_neighbors": 30,
        "authors_cluster_resolution": 0.3,
        "authors_variable_features_per_sample": 1000,
    }


def recover_sample_mapping(adata: ad.AnnData, paths: PipelinePaths) -> pd.DataFrame:
    """Recover paper-style sample IDs and verify a one-to-one mapping to GEO IDs.

    The authors identify samples as A1, A2, ..., H8, while GEO uses GSM accessions.
    Stage 03 needs the internal IDs because A4 and A6 are integration references.
    """

    if "sample_prefix" in adata.obs:
        prefixes = adata.obs["sample_prefix"].astype(str)
    else:
        prefixes = adata.obs_names.astype(str).str.split("_").str[0]
    valid = prefixes.str.match(r"^[A-Z]\d+$")
    if not bool(valid.all()):
        examples = prefixes[~valid].unique()[:10].tolist()
        raise ValueError(f"Invalid internal sample prefixes: {examples}")
    adata.obs["internal_sample_id"] = pd.Categorical(prefixes)

    # Every cell in one GSM sample must agree on its sample-level metadata.
    consistency_columns = [
        column
        for column in ("internal_sample_id", "title", "diagnosis", "age", "sex")
        if column in adata.obs
    ]
    for column in consistency_columns:
        counts = adata.obs.groupby("gsm", observed=True)[column].nunique(dropna=False)
        if (counts > 1).any():
            raise ValueError(f"Metadata column {column!r} varies within a GSM sample")

    table_columns = [
        column
        for column in (
            "gsm",
            "internal_sample_id",
            "patient_id",
            "title",
            "diagnosis",
            "cognitive_group",
            "age",
            "sex",
        )
        if column in adata.obs
    ]
    # Collapse the cell-level AnnData metadata to one representative row per sample,
    # then add the observed number of cells as an additional alignment check.
    sample_table = (
        adata.obs[table_columns]
        .reset_index(names="cell_barcode")
        .groupby("gsm", observed=True, as_index=False)
        .first()
    )
    cell_counts = adata.obs["gsm"].value_counts().rename("observed_cells")
    sample_table = sample_table.merge(
        cell_counts, left_on="gsm", right_index=True, validate="one_to_one"
    )
    sample_table["internal_sample_id"] = sample_table["internal_sample_id"].astype(str)
    if (
        sample_table["gsm"].nunique() != EXPECTED_SAMPLES
        or sample_table["internal_sample_id"].nunique() != EXPECTED_SAMPLES
        or not sample_table["internal_sample_id"].str.match(re.compile(r"^[A-Z]\d+$")).all()
    ):
        raise ValueError("GSM-to-internal-sample mapping is not one-to-one")
    # These references come directly from the supplied author workflow.
    for reference in ("A4", "A6"):
        if int((sample_table["internal_sample_id"] == reference).sum()) != 1:
            raise ValueError(f"Integration reference {reference} is not uniquely mapped")

    sample_table = sample_table.sort_values("internal_sample_id").reset_index(drop=True)
    output = paths.figure1d_dir / "GSE200164_sample_metadata_mapping.csv"
    sample_table.to_csv(output, index=False)
    print(f"Saved sample mapping: {output}")
    return sample_table


def export_seurat_input(
    adata: ad.AnnData,
    paths: PipelinePaths,
    *,
    force: bool = False,
) -> dict[str, Path]:
    """Export the four aligned files used to transfer AnnData into Seurat.

    Matrix Market stores the sparse numeric matrix. Separate text files preserve the
    exact gene and cell order, and metadata.csv carries cell annotations into R.
    """

    outputs = {
        "matrix": paths.seurat_input_dir / "counts.mtx",
        "genes": paths.seurat_input_dir / "genes.tsv",
        "barcodes": paths.seurat_input_dir / "barcodes.tsv",
        "metadata": paths.seurat_input_dir / "metadata.csv",
    }
    if all(path.exists() for path in outputs.values()) and not force:
        print(f"Using existing Seurat export: {paths.seurat_input_dir}")
        return outputs
    if "counts" not in adata.layers:
        raise KeyError("AnnData is missing the counts layer")
    counts = adata.layers["counts"]
    if not sparse.issparse(counts):
        raise TypeError("A sparse AnnData counts layer is required for Seurat export")
    # AnnData is cells-by-genes; Seurat expects genes-by-cells. COO is the sparse format
    # accepted by SciPy's Matrix Market writer.
    counts_for_seurat = counts.transpose().tocoo()
    if counts_for_seurat.shape != (adata.n_vars, adata.n_obs):
        raise ValueError("Transposed Seurat matrix has unexpected dimensions")

    scipy.io.mmwrite(outputs["matrix"], counts_for_seurat, field="real")
    pd.Series(adata.var_names.astype(str)).to_csv(
        outputs["genes"], sep="\t", index=False, header=False
    )
    pd.Series(adata.obs_names.astype(str)).to_csv(
        outputs["barcodes"], sep="\t", index=False, header=False
    )
    # Include barcodes as an explicit column so R can align metadata without relying
    # on pandas' index serialization.
    metadata = adata.obs.copy()
    metadata.insert(0, "cell_barcode", adata.obs_names.astype(str))
    metadata.to_csv(outputs["metadata"], index=False)
    for label, path in outputs.items():
        print(f"{label:9s}: {path} ({path.stat().st_size / 1024**3:.3f} GB)")
    return outputs


def prepare_figure1d(
    paths: PipelinePaths,
    *,
    force: bool = False,
    rscript: str | None = None,
) -> None:
    """Run Stage 02 and avoid exporting raw R inputs when checkpoints can resume."""

    paths.create_directories()
    adata = load_and_audit_canonical(paths)
    add_figure1d_qc(adata)
    recover_sample_mapping(adata, paths)
    atomic_h5ad_write(adata, paths.figure1d_checkpoint_h5ad, force=force)
    # Any of these objects lets Stage 03 start after the raw Matrix Market boundary.
    # Skipping that export saves substantial disk I/O on a resumed run.
    reusable_r_checkpoints = (
        paths.figure1d_dir / "GSE200164_after_PCA.rds",
        paths.figure1d_dir / "GSE200164_seurat_integrated.rds",
        paths.figure1d_dir / "integrated_clustered_umap.rds",
    )
    if any(path.exists() for path in reusable_r_checkpoints) and not force:
        print("Reusable Seurat checkpoint found; large Matrix Market export skipped.")
    else:
        export_seurat_input(adata, paths, force=force)
    audit_r_environment(paths, rscript=rscript)
    print(f"Figure 1D input checkpoint: {paths.figure1d_checkpoint_h5ad}")
