#!/usr/bin/env python3
"""Stage 02: read Stage 01, calculate QC, and prepare Seurat input.

This stage validates Stage 01, calculates the three QC fields expected by the authors'
Seurat workflow, recovers the internal sample IDs, and exports a sparse R interchange
format only when Stage 03 cannot resume from an existing checkpoint. Everything the
Python stage needs is defined in this file.
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import anndata as ad
import numpy as np
import pandas as pd
import psutil
import scipy.io
from scipy import sparse

# --- Dataset and software expectations --------------------------------------
# The same fixed dimensions Stage 01 checked, repeated here so a truncated or
# mismatched hand-off is caught the moment this stage opens the file.
#
# The R packages listed are what the published analysis needs: Seurat performs the
# integration and clustering, Matrix carries the sparse data into it, glmGamPoi
# makes the per-sample normalization tractable, and ggplot2 draws the final figure.
EXPECTED_CELLS = 70_391
EXPECTED_GENES = 19_689
EXPECTED_SAMPLES = 59
REQUIRED_R_PACKAGES = ("Seurat", "Matrix", "glmGamPoi", "ggplot2")


# --- Workspace paths ---------------------------------------------------------
# Stage 02 reads exactly one file from Stage 01 and writes everything it produces
# into the results folder, where Stage 03 will look for it.
@dataclass(frozen=True)
class PipelinePaths:
    """All Stage 02 input and output paths below one workspace directory."""

    workspace: Path

    @classmethod
    def from_workspace(cls, workspace: str | Path) -> "PipelinePaths":
        return cls(Path(workspace).expanduser().resolve())

    # The finished object Stage 01 produced: 70,391 cells with patient details.
    @property
    def canonical_h5ad(self) -> Path:
        return self.workspace / "data" / "processed" / "GSE200164_canonical_input.h5ad"

    @property
    def figure1d_dir(self) -> Path:
        return self.workspace / "results" / "figure1d_python"

    # Where the data gets handed across to R. Only needed when Stage 03 is starting
    # from nothing; a resumed run skips writing these entirely.
    @property
    def seurat_input_dir(self) -> Path:
        return self.figure1d_dir / "seurat_input"

    # Stage 02's own saved result: the same cells, now carrying quality measurements.
    @property
    def figure1d_checkpoint_h5ad(self) -> Path:
        return self.figure1d_dir / "GSE200164_after_qc_audit.h5ad"

    def create_directories(self) -> None:
        """Create this stage's folders, leaving any existing results in place."""

        for directory in (self.workspace, self.figure1d_dir, self.seurat_input_dir):
            directory.mkdir(parents=True, exist_ok=True)


def default_workspace() -> Path:
    """Use the repository's Age_trajectory directory unless --workspace is given."""

    return Path(__file__).resolve().parents[1] / "Age_trajectory"


# --- General file helpers ----------------------------------------------------
def require_files(paths: Iterable[Path]) -> None:
    """Stop with a readable list when an expected input file is missing."""

    missing = [path for path in paths if not path.is_file()]
    if missing:
        formatted = "\n".join(f"  - {path}" for path in missing)
        raise FileNotFoundError(f"Required files were not found:\n{formatted}")


def atomic_h5ad_write(adata: Any, destination: Path, *, force: bool = False) -> None:
    """Save a result so a crash mid-write cannot destroy a good existing file."""

    if destination.exists() and not force:
        print(f"Using existing AnnData file: {destination}")
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    # Write under a temporary name and rename only on success, so an interrupted
    # save leaves the previous good file intact rather than a half-written one.
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    if temporary.exists():
        temporary.unlink()
    adata.write_h5ad(
        temporary,
        compression="gzip",
        convert_strings_to_categoricals=False,
    )
    temporary.replace(destination)


# --- R and computer-resource preflight --------------------------------------
# The heavy analysis in Stage 03 happens in R, because the published method uses
# Seurat. This check runs R briefly to confirm it is installed with the packages
# that analysis needs, and records how much memory and disk the machine has.
#
# The point is to fail in seconds rather than hours in. Without this, a missing
# package would surface only after Stage 03 had already spent a long time
# normalizing all 59 samples, and that work would be wasted.
def audit_r_environment(
    paths: PipelinePaths,
    *,
    rscript: str | None = None,
) -> dict[str, object]:
    """Confirm R and its required packages are available before any long work."""

    executable = rscript or shutil.which("Rscript")
    if executable is None:
        raise RuntimeError("Rscript was not found on PATH")
    # Ask R to report, for each required package, whether it is installed and at
    # what version. Versions are worth recording: Seurat's results shift between
    # releases, which is part of why this reproduction differs from the paper's.
    package_vector = ", ".join(json.dumps(package) for package in REQUIRED_R_PACKAGES)
    code = f"""
packages <- c({package_vector})
for (pkg in packages) {{
  installed <- requireNamespace(pkg, quietly = TRUE)
  version <- if (installed) as.character(packageVersion(pkg)) else NA_character_
  cat(pkg, installed, version, sep = "\\t")
  cat("\\n")
}}
"""
    result = subprocess.run(
        [executable, "-e", code], check=True, capture_output=True, text=True
    )
    packages: dict[str, dict[str, object]] = {}
    for line in result.stdout.splitlines():
        name, installed, version = line.split("\t")
        packages[name] = {
            "installed": installed == "TRUE",
            "version": None if version == "NA" else version,
        }
    missing = [name for name in REQUIRED_R_PACKAGES if not packages[name]["installed"]]
    # Note the machine's memory and free disk alongside the package versions. The
    # integration in Stage 03 is demanding on both, so if a run later fails this
    # record shows what the machine had available at the time.
    memory = psutil.virtual_memory()
    disk = shutil.disk_usage(paths.figure1d_dir)
    audit = {
        "python_version": sys.version,
        "rscript_path": executable,
        "required_r_packages": packages,
        "missing_r_packages": missing,
        "total_ram_gb": memory.total / 1024**3,
        "available_ram_gb": memory.available / 1024**3,
        "disk_free_gb": disk.free / 1024**3,
        "ready": not missing,
    }
    # Save it as plain text next to the results, so the environment that produced
    # a given figure can be inspected later without re-running anything.
    output = paths.figure1d_dir / "figure1d_environment_audit.json"
    output.write_text(json.dumps(audit, indent=2, default=str) + "\n", encoding="utf-8")
    print(f"R environment audit: missing packages = {missing or 'none'}")
    return audit


# --- Step 1: load and validate Stage 01 --------------------------------------
# Open the object Stage 01 built and confirm it is what this stage expects: the
# right number of cells and genes, no cell or gene appearing twice, and no negative
# measurements, which would mean the data had been corrupted or transformed.
#
# This step also reports how many of the values are whole numbers. It is expected
# to be a low fraction, and that is not a problem — this dataset was corrected for
# ambient RNA contamination, which leaves fractional measurements behind. The
# printed figure is there so the character of the data is visible rather than
# assumed, since much single-cell software takes whole-number counts for granted.
def load_and_audit_canonical(paths: PipelinePaths) -> ad.AnnData:
    """Open Stage 01's result and confirm it is intact before building on it."""

    require_files([paths.canonical_h5ad])
    adata = ad.read_h5ad(paths.canonical_h5ad)
    if adata.shape != (EXPECTED_CELLS, EXPECTED_GENES):
        raise ValueError(
            f"Canonical AnnData has shape {adata.shape}; "
            f"expected {(EXPECTED_CELLS, EXPECTED_GENES)}"
        )
    if not sparse.issparse(adata.X):
        raise TypeError("Canonical AnnData.X must remain sparse")
    # Arrange the data so that reading one cell's full set of measurements is fast,
    # which is what the per-cell quality calculations in the next step need.
    if not sparse.isspmatrix_csr(adata.X):
        adata.X = adata.X.tocsr()
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ValueError("Canonical cell and gene names must be unique")

    # Examine only the recorded measurements rather than the full grid, since the
    # vast majority of gene-cell combinations were never observed at all.
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


# --- Step 2: calculate the quality-control columns used by the authors -------
# Not every cell in a single-cell experiment is a healthy intact cell. Some droplets
# captured debris, some captured a dying cell, some captured almost nothing. The
# standard way to tell them apart is three measurements per cell:
#
#   nCount_RNA    how much RNA was measured in total. Very low means little was
#                 captured; very high can mean two cells were caught together.
#   nFeature_RNA  how many distinct genes were detected at all. A healthy cell
#                 expresses thousands; a few hundred suggests an empty droplet.
#   percent.mt    what share of the cell's RNA came from mitochondrial genes. When
#                 a cell is dying its ordinary RNA leaks out while mitochondrial
#                 RNA remains, so a high share is a signature of a damaged cell.
#
# The paper treated cells with fewer than 200 detected genes or 10% or more
# mitochondrial content as failing. This step calculates all three and reports how
# many cells fall outside those bounds, but deliberately removes nothing. The
# published dataset was already filtered by its authors, so any cells outside the
# thresholds are a fact about the source data worth surfacing rather than quietly
# discarding. Stage 03 receives every cell and the decision stays visible.
def add_figure1d_qc(adata: ad.AnnData, *, sample_key: str = "gsm") -> None:
    """Measure each cell's RNA total, gene count, and mitochondrial share."""

    # Confirm every cell traces back to one of the 59 patients before measuring
    # anything, since the per-sample normalization in Stage 03 depends on it.
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

    # The first two measurements come straight from each cell's row: the sum of
    # everything measured, and how many genes were detected at all.
    counts = adata.X
    adata.obs["nCount_RNA"] = np.asarray(counts.sum(axis=1)).ravel()
    adata.obs["nFeature_RNA"] = np.diff(counts.indptr)
    # Mitochondrial genes are identifiable by name — in human data they all begin
    # with "MT-". Their share of each cell's total is the third measurement.
    mitochondrial = adata.var_names.astype(str).str.startswith("MT-")
    if int(mitochondrial.sum()) == 0:
        raise ValueError("No MT- genes were found; percent.mt cannot be calculated")
    mitochondrial_counts = np.asarray(counts[:, mitochondrial].sum(axis=1)).ravel()
    totals = adata.obs["nCount_RNA"].to_numpy()
    # A cell with nothing measured would give an undefined share; those are recorded
    # as zero rather than left as an error value.
    adata.obs["percent.mt"] = np.divide(
        mitochondrial_counts,
        totals,
        out=np.zeros(adata.n_obs, dtype=np.float64),
        where=totals > 0,
    ) * 100

    # Report how many cells fall outside the paper's thresholds, without removing
    # them, as described above.
    below_200 = int((adata.obs["nFeature_RNA"] < 200).sum())
    high_mitochondrial = int((adata.obs["percent.mt"] >= 10).sum())
    print(adata.obs[["nCount_RNA", "nFeature_RNA", "percent.mt"]].describe().T)
    print(f"Cells below 200 features: {below_200:,}")
    print(f"Cells at or above 10% mitochondrial expression: {high_mitochondrial:,}")
    if below_200 or high_mitochondrial:
        print("WARNING: some cells do not meet the paper's measurable QC thresholds")

    # Store the paper's analysis settings inside the file itself. These are the
    # numbers Stage 03 reproduces — how many dimensions to summarize the data down
    # to, how finely to divide the cells into clusters, and so on. Keeping them with
    # the data means anyone opening this file can see the settings the figure was
    # made with, without needing this code alongside it.
    adata.uns["figure1d_reproduction"] = {
        "dataset": "GSE200164",
        "sample_key": sample_key,
        "expected_samples": EXPECTED_SAMPLES,
        "expected_cells": EXPECTED_CELLS,
        "expected_genes": EXPECTED_GENES,
        "input_interpretation": "Processed SoupX-adjusted quality-controlled count-like matrix",
        "authors_initial_pcs": 15,
        "authors_umap_neighbors": 30,
        "authors_cluster_resolution": 0.3,
        "authors_variable_features_per_sample": 1000,
    }


# --- Step 3: recover the authors' A1-H8 sample labels ------------------------
# The same patient goes by two different names. The public repository calls them
# things like "GSM6016346", while the paper and its code call them "A1" through
# "H8". Both appear in the data, and this step establishes that they refer to the
# same 59 people, one for one.
#
# This matters for a specific reason. When Stage 03 merges the 59 samples together,
# it needs two of them to act as the reference the others are aligned against, and
# the authors chose A4 and A6. Those names only exist in the paper's numbering, so
# without this mapping the merge could not be reproduced as published.
#
# The step also confirms that every cell belonging to one patient agrees about that
# patient's age, sex, and diagnosis. Disagreement would mean the details had been
# attached wrongly somewhere upstream.
def recover_sample_mapping(adata: ad.AnnData, paths: PipelinePaths) -> pd.DataFrame:
    """Match the repository's patient names to the paper's, and check they agree."""

    # The paper's identifier is the prefix already carried on every cell's label.
    if "sample_prefix" in adata.obs:
        prefixes = adata.obs["sample_prefix"].astype(str)
    else:
        prefixes = adata.obs_names.astype(str).str.split("_").str[0]
    # These identifiers always read as a letter followed by digits, so anything else
    # means the labels are not what this code expects.
    valid = prefixes.str.match(r"^[A-Z]\d+$")
    if not bool(valid.all()):
        examples = prefixes[~valid].unique()[:10].tolist()
        raise ValueError(f"Invalid internal sample prefixes: {examples}")
    adata.obs["internal_sample_id"] = pd.Categorical(prefixes)

    # Check that all of one patient's cells report the same details for that patient.
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
    # Condense the 70,391 cell records down to one row per patient, and record how
    # many cells each contributed — a useful sanity figure, since a patient with far
    # fewer cells than the rest will carry less weight in the merged result.
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
    # Confirm the two naming systems line up exactly: 59 patients under each, with
    # no name shared between two patients and none left unmatched.
    if (
        sample_table["gsm"].nunique() != EXPECTED_SAMPLES
        or sample_table["internal_sample_id"].nunique() != EXPECTED_SAMPLES
        or not sample_table["internal_sample_id"].str.match(re.compile(r"^[A-Z]\d+$")).all()
    ):
        raise ValueError("GSM-to-internal-sample mapping is not one-to-one")
    # Confirm the two reference patients the authors merged against are present and
    # unambiguous, since Stage 03 cannot reproduce the published merge without them.
    for reference in ("A4", "A6"):
        if int((sample_table["internal_sample_id"] == reference).sum()) != 1:
            raise ValueError(f"Integration reference {reference} is not uniquely mapped")

    # Save the mapping as a small spreadsheet-readable table. It is the clearest
    # single view of the study: 59 rows showing each patient under both names, with
    # their age, sex, diagnosis, and how many cells they contributed.
    sample_table = sample_table.sort_values("internal_sample_id").reset_index(drop=True)
    output = paths.figure1d_dir / "GSE200164_sample_metadata_mapping.csv"
    sample_table.to_csv(output, index=False)
    print(f"Saved sample mapping: {output}")
    return sample_table


# --- Step 4/5: export aligned files for R/Seurat when needed -----------------
# This is the hand-off from Python to R. The analysis continues in Seurat, which
# cannot read the format used up to now, so the data is written out in a form R
# understands: the measurements in one file, and the gene names, cell names, and
# per-cell details in three companion files beside it.
#
# The matrix is flipped on its way out. Everything so far has treated cells as rows;
# Seurat expects genes as rows. The three companion files record the exact ordering
# so R can put the names back on the right rows and columns.
#
# This export is around 4 GB and takes a while, which is why the controller below
# skips it whenever Stage 03 already has a saved checkpoint to resume from.
def export_seurat_input(
    adata: ad.AnnData,
    paths: PipelinePaths,
    *,
    force: bool = False,
) -> dict[str, Path]:
    """Write the measurements and labels out in the form R/Seurat can read."""

    # The four files are only meaningful together, since their ordering is what
    # links names back to the matrix. A partial export is treated as no export.
    outputs = {
        "matrix": paths.seurat_input_dir / "counts.mtx",
        "genes": paths.seurat_input_dir / "genes.tsv",
        "barcodes": paths.seurat_input_dir / "barcodes.tsv",
        "metadata": paths.seurat_input_dir / "metadata.csv",
    }
    if all(path.exists() for path in outputs.values()) and not force:
        print(f"Using existing Seurat export: {paths.seurat_input_dir}")
        return outputs
    counts = adata.X
    if not sparse.issparse(counts):
        raise TypeError("A sparse expression matrix is required for Seurat export")
    # Flip the matrix so genes are rows, as described above.
    counts_for_seurat = counts.transpose().tocoo()
    if counts_for_seurat.shape != (adata.n_vars, adata.n_obs):
        raise ValueError("Transposed Seurat matrix has unexpected dimensions")

    # Write the measurements, then the gene names, cell names, and per-cell details.
    # The measurements are written as decimal numbers rather than whole ones, which
    # is essential here: rounding would destroy the ambient-RNA correction that the
    # whole pipeline has been careful to preserve.
    scipy.io.mmwrite(outputs["matrix"], counts_for_seurat, field="real")
    pd.Series(adata.var_names.astype(str)).to_csv(
        outputs["genes"], sep="\t", index=False, header=False
    )
    pd.Series(adata.obs_names.astype(str)).to_csv(
        outputs["barcodes"], sep="\t", index=False, header=False
    )
    # Repeat each cell's name as an ordinary column so R can match the details to
    # the matrix by name rather than trusting that row order survived the transfer.
    metadata = adata.obs.copy()
    metadata.insert(0, "cell_barcode", adata.obs_names.astype(str))
    metadata.to_csv(outputs["metadata"], index=False)
    for label, path in outputs.items():
        print(f"{label:9s}: {path} ({path.stat().st_size / 1024**3:.3f} GB)")
    return outputs


# --- Stage 02 controller -----------------------------------------------------
# The whole stage in order. Read this function to see what Stage 02 does; read the
# steps above for how each part works.
#
# Going in: Stage 01's object, 70,391 cells knowing which patient they came from.
# Coming out: the same cells with three quality measurements each, both patient
# naming systems matched up, and the data written out in a form R can read.
#
# Nothing here changes a single measurement. This stage describes and repackages
# the data; the analysis itself is Stage 03's job.
def prepare_figure1d(
    paths: PipelinePaths,
    *,
    force: bool = False,
    rscript: str | None = None,
) -> None:
    """Measure cell quality, match the patient names, and hand off to R."""

    paths.create_directories()

    # Step 1: open Stage 01's result and confirm it arrived intact.
    adata = load_and_audit_canonical(paths)

    # Step 2: measure each cell's RNA total, gene count, and mitochondrial share.
    add_figure1d_qc(adata)

    # Step 3: match the repository's patient names to the paper's, and save the
    # 59-row summary table of the study.
    recover_sample_mapping(adata, paths)

    # Step 4: save this stage's own result, so a rerun need not repeat the above.
    atomic_h5ad_write(adata, paths.figure1d_checkpoint_h5ad, force=force)

    # Step 5: write the data out for R, but only if it will actually be needed.
    # Once Stage 03 has produced any of its saved checkpoints it resumes from those
    # instead, and this 4 GB export would be written only to sit unread.
    reusable_r_checkpoints = (
        paths.figure1d_dir / "GSE200164_after_PCA.rds",
        paths.figure1d_dir / "GSE200164_seurat_integrated.rds",
        paths.figure1d_dir / "integrated_clustered_umap.rds",
    )
    if any(path.exists() for path in reusable_r_checkpoints) and not force:
        print("Reusable Seurat checkpoint found; large Matrix Market export skipped.")
    else:
        export_seurat_input(adata, paths, force=force)

    # Step 6: confirm R is ready now rather than discovering it is not once Stage 03
    # is already hours into the analysis.
    audit_r_environment(paths, rscript=rscript)
    print(f"Figure 1D input checkpoint: {paths.figure1d_checkpoint_h5ad}")


# --- Command-line interface -------------------------------------------------
# What can be typed after the filename when running this stage: where the generated
# files live, which R installation to use, and whether to redo finished work. The
# scientific settings are not options; they live in the steps above, in context.
def parse_args() -> argparse.Namespace:
    """Read Stage 02 command-line options."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--workspace",
        default=default_workspace(),
        help="Generated-data workspace (default: checkout-local Age_trajectory)",
    )
    parser.add_argument("--rscript", help="Optional explicit path to Rscript")
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite the Figure 1D checkpoint and Seurat export",
    )
    return parser.parse_args()


def main() -> None:
    """Run the complete second stage from this one file."""

    # Example: python scripts/02_prepare_figure1d.py
    args = parse_args()
    paths = PipelinePaths.from_workspace(args.workspace)
    prepare_figure1d(paths, force=args.force, rscript=args.rscript)


if __name__ == "__main__":
    main()
