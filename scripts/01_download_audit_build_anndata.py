#!/usr/bin/env python3
"""Stage 01: download, audit, and build the canonical GSE200164 AnnData file.

The source expression file looks like a dense CSV, but loading it all at once would
need excessive memory. This module reads a small group of genes at a time, builds a
sparse matrix, verifies that conversion against the source, and attaches GEO sample
metadata before saving the canonical ``.h5ad`` file. It is self-contained: run this
file directly and do not import another project file.
"""

from __future__ import annotations

import argparse
import csv
import gc
import gzip
import hashlib
import json
import random
import re
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

# --- Dataset configuration ---------------------------------------------------
# The two public files this stage is built around: one holds the expression
# measurements for every cell, the other describes the 59 patients they came from.
GEO_ACCESSION = "GSE200164"
COUNTS_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE200nnn/"
    "GSE200164/suppl/GSE200164_counts.csv.gz"
)
SERIES_MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE200nnn/"
    "GSE200164/matrix/GSE200164_series_matrix.txt.gz"
)


# --- Workspace paths ---------------------------------------------------------
# Every file this stage reads or writes lives under one workspace directory, so a
# run can be pointed at a scratch location without touching finished results.
# Naming them all here means a filename is defined in exactly one place.
@dataclass(frozen=True)
class PipelinePaths:
    """All Stage 01 input and output paths below one workspace directory."""

    workspace: Path

    @classmethod
    def from_workspace(cls, workspace: str | Path) -> "PipelinePaths":
        return cls(Path(workspace).expanduser().resolve())

    # Folders: raw downloads, intermediate matrices, the finished H5AD, metadata
    # exports, and the audit trail of logs and tables.
    @property
    def raw_dir(self) -> Path:
        return self.workspace / "data" / "raw"

    @property
    def interim_dir(self) -> Path:
        return self.workspace / "data" / "interim"

    @property
    def processed_dir(self) -> Path:
        return self.workspace / "data" / "processed"

    @property
    def metadata_dir(self) -> Path:
        return self.workspace / "data" / "metadata"

    @property
    def logs_dir(self) -> Path:
        return self.workspace / "logs"

    @property
    def tables_dir(self) -> Path:
        return self.workspace / "tables"

    @property
    def figure1d_dir(self) -> Path:
        return self.workspace / "results" / "03_seurat_reproduction"

    @property
    def seurat_input_dir(self) -> Path:
        return self.figure1d_dir / "seurat_input"

    # Individual files, in the order the stage produces them: the two downloads,
    # the sparse matrices and their row/column labels, then the H5AD outputs.
    @property
    def counts_gzip(self) -> Path:
        return self.raw_dir / "GSE200164_counts.csv.gz"

    @property
    def series_matrix_gzip(self) -> Path:
        return self.raw_dir / "GSE200164_series_matrix.txt.gz"

    @property
    def genes_by_cells_npz(self) -> Path:
        return self.interim_dir / "GSE200164_expression_genes_by_cells_csr_float32.npz"

    @property
    def cells_by_genes_npz(self) -> Path:
        return self.interim_dir / "GSE200164_expression_cells_by_genes_csr_float32.npz"

    @property
    def genes_csv(self) -> Path:
        return self.interim_dir / "GSE200164_genes.csv"

    @property
    def barcodes_txt(self) -> Path:
        return self.interim_dir / "GSE200164_barcodes.txt"

    # The measurements alone, before patient details are attached. Kept as a resume
    # point so an interrupted run need not repeat the conversion and audit.
    @property
    def unannotated_h5ad(self) -> Path:
        return self.processed_dir / "GSE200164_expression_unannotated.h5ad"

    # Stage 01's finished output, and the only file Stage 02 reads.
    @property
    def canonical_h5ad(self) -> Path:
        return self.processed_dir / "GSE200164_canonical_input.h5ad"

    def create_directories(self) -> None:
        """Create this stage's folders, leaving any existing data in place."""

        for directory in (
            self.workspace,
            self.raw_dir,
            self.interim_dir,
            self.processed_dir,
            self.metadata_dir,
            self.logs_dir,
            self.tables_dir,
            self.figure1d_dir,
            self.seurat_input_dir,
        ):
            directory.mkdir(parents=True, exist_ok=True)


def default_workspace() -> Path:
    """Use the repository's Age_trajectory directory unless --workspace is given."""

    return Path(__file__).resolve().parents[1] / "Age_trajectory"


# --- General file helpers ----------------------------------------------------
# Small utilities shared by the steps below. They exist so the scientific steps can
# read as a sequence of decisions rather than a mix of analysis and file plumbing.
#
# The recurring theme is safety on a long-running job: a download or save that is
# interrupted must never leave behind a file that later looks complete and valid.
def require_files(paths: Iterable[Path]) -> None:
    """Stop with a readable list when an expected input file is missing."""

    missing = [path for path in paths if not path.is_file()]
    if missing:
        formatted = "\n".join(f"  - {path}" for path in missing)
        raise FileNotFoundError(f"Required files were not found:\n{formatted}")


def write_json(path: Path, payload: dict[str, Any]) -> None:
    """Record what the pipeline observed, as plain text anyone can open."""

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    """Fingerprint the download so a future run can prove it got the same data."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def validate_gzip(path: Path) -> None:
    """Catch a truncated or corrupted download before hours are spent reading it."""

    # Reading the stream to its end forces gzip to verify its own checksum.
    with gzip.open(path, "rb") as handle:
        while handle.read(8 * 1024 * 1024):
            pass


def download_with_resume(url: str, destination: Path, *, force: bool = False) -> None:
    """Fetch a GEO file, resuming rather than restarting if the connection drops."""

    if destination.exists() and not force:
        print(f"Using existing download: {destination}")
        return
    curl = shutil.which("curl")
    if curl is None:
        raise RuntimeError("curl is required for GEO downloads but was not found on PATH")
    destination.parent.mkdir(parents=True, exist_ok=True)
    # Partial bytes accumulate under a ``.part`` name and are renamed only on
    # success, so an interrupted 1.1 GB download can never masquerade as complete.
    partial = destination.with_suffix(destination.suffix + ".part")
    subprocess.run(
        [
            curl,
            "--location",
            "--fail",
            "--retry",
            "5",
            "--retry-delay",
            "5",
            "--continue-at",
            "-",
            "--output",
            str(partial),
            url,
        ],
        check=True,
    )
    partial.replace(destination)


def atomic_h5ad_write(adata: Any, destination: Path, *, force: bool = False) -> None:
    """Save an H5AD so a crash mid-write cannot destroy a good existing file."""

    if destination.exists() and not force:
        print(f"Using existing AnnData file: {destination}")
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    # Same guarantee as the download: write under a temporary name, then rename.
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    if temporary.exists():
        temporary.unlink()
    adata.write_h5ad(
        temporary,
        compression="gzip",
        convert_strings_to_categoricals=False,
    )
    temporary.replace(destination)


# --- Known dataset layout ----------------------------------------------------
# The published GEO file has a fixed shape: six descriptive columns per gene,
# followed by one column per cell. Hard-coding the expected dimensions turns a
# partial download, a revised GEO release, or a wrong input file into an immediate
# error instead of a subtly wrong figure discovered hours later.
ANNOTATION_COLUMNS = ["gene_symbol", "chromosome", "start", "end", "width", "strand"]

EXPECTED_CELLS = 70_391
EXPECTED_GENES = 19_689
EXPECTED_SAMPLES = 59


# --- Step 1: download the two GEO inputs -------------------------------------
# Fetch both public files and prove they arrived intact, recording a checksum so a
# future run can demonstrate it started from identical source data.
def download_source_files(paths: PipelinePaths, *, force: bool = False) -> None:
    """Download, integrity-check, and record provenance for both GEO source files."""

    download_with_resume(COUNTS_URL, paths.counts_gzip, force=force)
    validate_gzip(paths.counts_gzip)
    digest = sha256_file(paths.counts_gzip)
    write_json(
        paths.logs_dir / "GSE200164_counts_download_metadata.json",
        {
            "geo_accession": GEO_ACCESSION,
            "source_url": COUNTS_URL,
            "local_file": str(paths.counts_gzip),
            "compressed_size_bytes": paths.counts_gzip.stat().st_size,
            "sha256": digest,
            "gzip_integrity_passed": True,
        },
    )
    print(f"Validated expression archive ({paths.counts_gzip.stat().st_size / 1024**3:.3f} GB)")
    print(f"SHA-256: {digest}")

    # The series matrix is small but carries the clinical variables the study is
    # about: each sample's age, diagnosis, and sex.
    download_with_resume(SERIES_MATRIX_URL, paths.series_matrix_gzip, force=force)
    validate_gzip(paths.series_matrix_gzip)


def _read_source_header(path: Path) -> list[str]:
    """Read the first CSV row, which names the gene columns and all 70,391 cells."""

    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        return next(csv.reader(handle))


# --- Step 2: convert the large dense CSV into sparse matrices ----------------
# GEO ships the expression data as a dense CSV where most entries are zero, which
# is far too large to hold in memory. This step rewrites it as a sparse matrix,
# storing only the nonzero measurements.
#
# Two orientations are saved: GEO lists genes as rows, AnnData expects cells as
# rows. Keeping both makes the transposition auditable and means the expensive
# conversion never has to run twice.
def convert_expression_to_sparse(
    paths: PipelinePaths,
    *,
    chunk_size: int = 100,
    force: bool = False,
) -> None:
    """Convert the dense GEO CSV stream into two sparse matrix orientations."""

    # Treat the step as done only if all four aligned outputs survived the last run.
    outputs = (
        paths.genes_by_cells_npz,
        paths.cells_by_genes_npz,
        paths.genes_csv,
        paths.barcodes_txt,
    )
    if all(path.exists() for path in outputs) and not force:
        print("Sparse conversion outputs already exist; conversion skipped.")
        return

    require_files([paths.counts_gzip])
    # Learn the dataset's true shape before committing to a long conversion: the
    # header names the cells, and a fast scan counts the genes.
    header = _read_source_header(paths.counts_gzip)
    cell_barcodes = header[len(ANNOTATION_COLUMNS) :]
    with gzip.open(paths.counts_gzip, "rt", encoding="utf-8") as handle:
        n_genes = sum(1 for _ in handle) - 1

    if len(cell_barcodes) != EXPECTED_CELLS or n_genes != EXPECTED_GENES:
        raise ValueError(
            "Unexpected GEO dimensions: "
            f"{len(cell_barcodes):,} cells and {n_genes:,} genes"
        )
    # Save the barcodes in their original order. Everything downstream assumes row i
    # of the matrix is this cell, so this ordering is the contract for the whole run.
    paths.barcodes_txt.write_text(
        "".join(f"{barcode}\n" for barcode in cell_barcodes), encoding="utf-8"
    )

    # Read the file a few hundred genes at a time rather than all at once. Smaller
    # chunks lower peak memory at the cost of more passes.
    chunks = pd.read_csv(
        paths.counts_gzip,
        compression="gzip",
        header=0,
        names=ANNOTATION_COLUMNS + cell_barcodes,
        chunksize=chunk_size,
        engine="c",
        low_memory=False,
    )
    sparse_blocks: list[sparse.csr_matrix] = []
    metadata_blocks: list[pd.DataFrame] = []
    processed = 0
    fractional_nonzero = 0
    minimum = np.inf
    maximum = -np.inf
    started = time.monotonic()

    # Convert one chunk at a time: split each into gene descriptions and expression
    # values, reject anything negative or non-numeric, and keep only the nonzero
    # measurements. Running tallies of the value range and how many values are
    # fractional become the conversion audit, at no extra cost in passes.
    #
    # Values stay fractional on purpose. This dataset is SoupX-corrected for ambient
    # RNA, so its "counts" are not whole numbers, and rounding would corrupt them.
    for chunk_number, chunk in enumerate(chunks, start=1):
        metadata = chunk.iloc[:, : len(ANNOTATION_COLUMNS)].copy()
        metadata.columns = ANNOTATION_COLUMNS
        values = chunk.iloc[:, len(ANNOTATION_COLUMNS) :].to_numpy(
            dtype=np.float32, copy=True
        )
        if not np.isfinite(values).all() or np.any(values < 0):
            raise ValueError(f"Invalid expression values in chunk {chunk_number}")
        minimum = min(minimum, float(values.min()))
        maximum = max(maximum, float(values.max()))
        nonzero = values[values != 0]
        fractional_nonzero += int(np.count_nonzero(nonzero != np.floor(nonzero)))

        block = sparse.csr_matrix(values, dtype=np.float32)
        block.eliminate_zeros()
        sparse_blocks.append(block)
        metadata_blocks.append(metadata)
        processed += block.shape[0]
        elapsed = (time.monotonic() - started) / 60
        print(
            f"\rChunk {chunk_number:4d}: {processed:,}/{n_genes:,} genes "
            f"({elapsed:.1f} min)",
            end="",
            flush=True,
        )
        # Release each chunk immediately; this loop runs for thousands of iterations
        # and peak memory is the limiting factor on an ordinary machine.
        del values, nonzero, chunk
        gc.collect()
    print()

    if processed != n_genes:
        raise RuntimeError(f"Processed {processed:,} genes; expected {n_genes:,}")
    # Reassemble the chunks in source order, then flip the matrix so cells are rows,
    # which is the orientation AnnData and every later stage expect.
    genes_by_cells = sparse.vstack(sparse_blocks, format="csr", dtype=np.float32)
    cells_by_genes = genes_by_cells.transpose().tocsr().astype(np.float32)
    cells_by_genes.sort_indices()
    cells_by_genes.eliminate_zeros()
    if cells_by_genes.shape != (EXPECTED_CELLS, EXPECTED_GENES):
        raise RuntimeError(f"Unexpected sparse matrix shape: {cells_by_genes.shape}")

    # Save both orientations and the labels that go with them, so this conversion
    # never needs to be repeated and its row/column ordering stays checkable.
    sparse.save_npz(paths.genes_by_cells_npz, genes_by_cells, compressed=True)
    sparse.save_npz(paths.cells_by_genes_npz, cells_by_genes, compressed=True)
    pd.concat(metadata_blocks, ignore_index=True).to_csv(paths.genes_csv, index=False)

    # Summarize the finished matrix once, so anyone judging whether the conversion
    # looks sane can read a small file instead of reloading gigabytes.
    cell_sums = np.asarray(cells_by_genes.sum(axis=1)).ravel()
    features_per_cell = np.diff(cells_by_genes.indptr)
    write_json(
        paths.logs_dir / "GSE200164_float32_sparse_conversion.json",
        {
            "matrix_shape_cells_by_genes": list(cells_by_genes.shape),
            "dtype": str(cells_by_genes.dtype),
            "nonzero_entries": int(cells_by_genes.nnz),
            "matrix_density": float(cells_by_genes.nnz / np.prod(cells_by_genes.shape)),
            "minimum_value": minimum,
            "maximum_value": maximum,
            "fractional_nonzero_values": fractional_nonzero,
            "mean_expression_sum_per_cell": float(cell_sums.mean()),
            "median_expression_sum_per_cell": float(np.median(cell_sums)),
            "mean_detected_genes_per_cell": float(features_per_cell.mean()),
            "median_detected_genes_per_cell": float(np.median(features_per_cell)),
        },
    )
    print(f"Saved sparse matrix: {cells_by_genes.shape}, {cells_by_genes.nnz:,} nonzeros")


# --- Step 3: prove the sparse conversion preserved values --------------------
# The previous step rewrote the entire dataset into a new format and transposed it.
# Every result in this project rests on that rewrite being faithful, so this step
# goes back to the original downloaded file and re-checks a sample of it by hand.
#
# Three things are checked, from coarsest to finest:
#
#   1. Structure. Does the saved matrix still have 70,391 cells and 19,689 genes,
#      are the values still fractional rather than rounded, and do the saved cell
#      labels appear in exactly the order the original file listed them? A silent
#      reordering here would attach every cell's measurements to the wrong patient.
#
#   2. Whole genes. For 14 chosen genes, every one of the 70,391 cell measurements
#      is compared against the original file. This catches an error affecting one
#      gene's entire row, which spot-checking individual values could easily miss.
#
#   3. Individual measurements. For those genes crossed with 22 chosen cells, each
#      single value is compared. The chosen genes and cells always include the very
#      first and very last of each, since off-by-one errors surface at the edges.
#
# The selection is random but seeded, so the audit examines the same entries on
# every run and its published results stay reproducible.
#
# Both comparison tables are written to disk before any pass/fail decision, so the
# evidence survives even when the audit fails.
def audit_sparse_conversion(paths: PipelinePaths) -> None:
    """Re-read the original download and confirm the converted matrix matches it."""

    require_files(
        [
            paths.counts_gzip,
            paths.cells_by_genes_npz,
            paths.genes_csv,
            paths.barcodes_txt,
        ]
    )
    # Check 1: the converted matrix has the expected size, its values were never
    # rounded to integers, and its cell ordering still matches the original file.
    expression = sparse.load_npz(paths.cells_by_genes_npz).tocsr()
    if expression.dtype != np.float32:
        raise TypeError(f"Expected float32 sparse data, found {expression.dtype}")
    genes = pd.read_csv(paths.genes_csv)["gene_symbol"].astype(str).tolist()
    barcodes = paths.barcodes_txt.read_text(encoding="utf-8").splitlines()
    if expression.shape != (len(barcodes), len(genes)):
        raise ValueError("Sparse matrix dimensions do not match genes and barcodes")

    source_header = _read_source_header(paths.counts_gzip)
    if source_header[len(ANNOTATION_COLUMNS) :] != barcodes:
        raise ValueError("Saved barcode order does not match the GEO CSV")

    # Choose which genes and cells to examine: the first and last of each, plus a
    # random spread. Seeding keeps the same selection across runs, and using a
    # private generator means this choice cannot disturb randomness elsewhere.
    generator = random.Random(42)
    gene_indices = sorted({0, len(genes) - 1, *generator.sample(range(len(genes)), 12)})
    cell_indices = sorted(
        {0, len(barcodes) - 1, *generator.sample(range(len(barcodes)), 20)}
    )
    # Pull just those genes back out of the original download, reading only as far
    # into the file as the last one needed rather than re-reading all 1.1 GB.
    selected_rows: dict[int, list[str]] = {}
    with gzip.open(paths.counts_gzip, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        next(reader)
        for gene_index, row in enumerate(reader):
            if gene_index in gene_indices:
                selected_rows[gene_index] = row
            if gene_index >= gene_indices[-1]:
                break

    # The genes-as-rows copy saved during conversion is used here. It holds exactly
    # the same numbers as the cells-as-rows copy, but reading one gene's full set of
    # measurements out of it is near-instant, whereas pulling that same gene out of
    # the cells-as-rows copy would mean touching all 70,391 cell records.
    by_gene = sparse.load_npz(paths.genes_by_cells_npz).tocsr()
    if by_gene.shape != (len(genes), len(barcodes)):
        raise ValueError("Genes-as-rows matrix does not match genes and barcodes")

    # Now run checks 2 and 3. For each chosen gene: compare its complete set of
    # measurements against the original file, then compare individual values at the
    # chosen cells. Both records become audit tables saved alongside the results.
    pair_records: list[dict[str, object]] = []
    vector_records: list[dict[str, object]] = []
    for gene_index in gene_indices:
        row = selected_rows[gene_index]
        if row[0] != genes[gene_index]:
            raise ValueError(f"Gene-order mismatch at index {gene_index}")
        # Comparisons are made at the same precision the data is stored in, so a
        # difference here means a real conversion error rather than rounding noise.
        source_vector = np.asarray(row[len(ANNOTATION_COLUMNS) :], dtype=np.float32)
        sparse_vector = by_gene[gene_index].toarray().ravel().astype(np.float32)
        vector_matches = np.isclose(source_vector, sparse_vector, rtol=1e-6, atol=1e-7)
        vector_records.append(
            {
                "gene_index": gene_index,
                "gene_symbol": genes[gene_index],
                "all_cells_match": bool(vector_matches.all()),
                "mismatched_positions": int((~vector_matches).sum()),
                "maximum_absolute_difference": float(
                    np.abs(source_vector.astype(np.float64) - sparse_vector).max()
                ),
                "original_nonzero_cells": int(np.count_nonzero(source_vector)),
                "sparse_nonzero_cells": int(np.count_nonzero(sparse_vector)),
                "original_total_expression": float(source_vector.sum(dtype=np.float64)),
                "sparse_total_expression": float(sparse_vector.sum(dtype=np.float64)),
            }
        )
        # Read individual values from the cells-as-rows matrix, the copy every later
        # stage actually consumes, so this confirms the file that will really be used.
        for cell_index in cell_indices:
            source_value = source_vector[cell_index]
            sparse_value = np.float32(expression[cell_index, gene_index])
            pair_records.append(
                {
                    "gene_index": gene_index,
                    "gene_symbol": genes[gene_index],
                    "cell_index": cell_index,
                    "cell_barcode": barcodes[cell_index],
                    "original_csv_float32": float(source_value),
                    "sparse_matrix_float32": float(sparse_value),
                    "absolute_difference": float(abs(np.float64(source_value) - sparse_value)),
                    "match": bool(np.isclose(source_value, sparse_value, rtol=1e-6, atol=1e-7)),
                }
            )

    # Save the evidence first, then judge it. If the audit fails, the tables showing
    # exactly which values disagreed are already on disk to be examined.
    pair_table = pd.DataFrame(pair_records)
    vector_table = pd.DataFrame(vector_records)
    pair_table.to_csv(paths.tables_dir / "GSE200164_sparse_conversion_audit.csv", index=False)
    vector_table.to_csv(
        paths.tables_dir / "GSE200164_sparse_full_vector_audit.csv", index=False
    )
    if not pair_table["match"].all() or not vector_table["all_cells_match"].all():
        raise AssertionError("Sparse conversion audit failed")
    print(f"Sparse audit passed: {len(pair_table):,} pairs and {len(vector_table)} vectors")


# --- Step 4: place expression and labels into an AnnData object --------------
# So far the measurements, the gene descriptions, and the cell labels are three
# separate files that only line up by convention. This step binds them into a
# single object where that alignment is guaranteed and travels with the data.
#
# The result is a table of 70,391 cells by 19,689 genes carrying, per cell, which
# sample it came from, and per gene, its name and genomic location. Clinical fields
# such as age and diagnosis are created here but left as "unknown"; the next step
# fills them from the patient metadata file.
#
# One detail is worth knowing: each cell's label in the original data looks like
# "A1_AAACCTGGTTGGTGGA-1", which packs two facts together — "A1" identifies the
# patient sample, and the rest identifies the individual cell within that sample.
# Splitting them apart is what makes it possible to attach patient information to
# every cell later, and to group cells by sample during integration in Stage 03.
def build_unannotated_anndata(paths: PipelinePaths, *, force: bool = False) -> None:
    """Combine the matrix, gene descriptions, and cell labels into one object."""

    if paths.unannotated_h5ad.exists() and not force:
        print(f"Using existing AnnData file: {paths.unannotated_h5ad}")
        return
    require_files([paths.cells_by_genes_npz, paths.genes_csv, paths.barcodes_txt])
    # Reload the three separate pieces and confirm they still describe each other:
    # one label per matrix row, one gene description per matrix column.
    expression = sparse.load_npz(paths.cells_by_genes_npz).tocsr().astype(np.float32)
    genes = pd.read_csv(paths.genes_csv, low_memory=False)
    barcodes = paths.barcodes_txt.read_text(encoding="utf-8").splitlines()
    if expression.shape != (len(barcodes), len(genes)):
        raise ValueError("AnnData inputs have inconsistent dimensions")

    # Tidy the gene descriptions: names and chromosomes become plain text, and the
    # genomic coordinates become numbers, so later filtering behaves predictably.
    for column in ("gene_symbol", "chromosome", "strand"):
        genes[column] = genes[column].fillna("unknown").astype(str)
    for column in ("start", "end", "width"):
        genes[column] = pd.to_numeric(genes[column], errors="coerce")

    # Separate each cell label into its sample identifier and its within-sample
    # identifier, as described above.
    split_barcodes = pd.Series(barcodes).str.split("_", n=1, expand=True)
    if split_barcodes.shape[1] != 2:
        raise ValueError("Cell barcodes do not contain sample prefixes")
    # Build the per-cell description. The clinical fields are deliberately declared
    # here as empty placeholders so the object's structure is complete and obvious
    # from the start, even though the real values only arrive in the next step.
    obs = pd.DataFrame(index=pd.Index(barcodes, name="cell_barcode"))
    obs["sample_prefix"] = split_barcodes.iloc[:, 0].to_numpy(dtype=str)
    obs["tenx_barcode"] = split_barcodes.iloc[:, 1].to_numpy(dtype=str)
    obs["patient_id"] = "unknown"
    obs["diagnosis"] = "unknown"
    obs["sex"] = "unknown"
    obs["cell_type"] = "unknown"
    obs["age"] = np.full(len(obs), np.nan, dtype=np.float32)

    # Join the measurements, the gene descriptions, and the cell descriptions into
    # the finished object. A handful of gene symbols appear twice in the source; the
    # repeats are given distinguishing names so each column can be referred to
    # unambiguously, while the original symbol stays recorded alongside.
    var = genes.copy()
    var.index = pd.Index(var["gene_symbol"].astype(str), name="gene_id")
    adata = ad.AnnData(X=expression, obs=obs, var=var)
    adata.var_names_make_unique(join="-duplicate-")
    # Record facts about the dataset as a whole, so anyone opening this file later
    # can tell where it came from and what has and has not been done to it. Stating
    # that no normalization has been applied matters: the values look unusual
    # (fractional rather than whole counts) and this confirms that is the source
    # data's own character, not something this pipeline introduced.
    adata.uns["dataset"] = {
        "geo_accession": GEO_ACCESSION,
        "source_filename": paths.counts_gzip.name,
        "matrix_description": "Audited GEO expression matrix stored as sparse float32.",
        "expression_status": "Fractional count-like values; untreated raw UMI counts are not assumed.",
    }
    adata.uns["processing"] = {
        "matrix_orientation": "cells_by_genes",
        "expression_modified": False,
        "normalization_applied_by_this_pipeline": False,
        "log_transform_applied_by_this_pipeline": False,
        "cell_metadata_attached": False,
    }
    # Save the object, immediately reopen it to confirm it was written correctly,
    # and note its size and contents in the run log.
    atomic_h5ad_write(adata, paths.unannotated_h5ad, force=True)
    validate_h5ad(paths.unannotated_h5ad, expected_shape=adata.shape)
    write_json(
        paths.logs_dir / "GSE200164_anndata_construction_summary.json",
        {
            "output_file": str(paths.unannotated_h5ad),
            "n_cells": adata.n_obs,
            "n_genes": adata.n_vars,
            "n_nonzero": int(adata.X.nnz),
            "dtype": str(adata.X.dtype),
            "sparse_format": adata.X.getformat(),
            "save_validation_passed": True,
        },
    )


# --- Step 5a: parse sample metadata from GEO's series matrix -----------------
# The patient information arrives in an awkward shape. Instead of one row per
# patient, the file is turned on its side: each line describes a single attribute
# across all 59 patients at once, so ages sit on one line, diagnoses on another.
# Worse, several different attributes share the same line label, distinguished only
# by text inside each value like "age: 67" or "diagnosis: Healthy".
#
# This step turns that sideways layout into an ordinary table of 59 rows, one per
# patient, with columns for sample identifier, age, sex, and diagnosis.
#
# Two cleanups happen along the way. Diagnosis labels are spelled inconsistently
# in the source ("AD", "Alzheimer's disease"), so they are mapped onto one spelling
# each — capitalization only, never inventing a diagnosis that was not recorded.
# And patients are grouped into healthy versus cognitively impaired, which is the
# comparison the study is built around.
#
# Because a shift in this file would silently mismatch patients to their data, the
# parse ends by re-checking one known patient against values verified by hand.
def parse_geo_sample_metadata(series_matrix: Path) -> pd.DataFrame:
    """Turn GEO's sideways metadata layout into one clean row per patient."""

    # Collect every attribute line, keyed by its label. Labels that appear more than
    # once are kept together, since that is how age, diagnosis, and sex arrive.
    metadata_rows: dict[str, list[list[str]]] = {}
    with gzip.open(series_matrix, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("!Sample_"):
                continue
            fields = line.rstrip("\n").split("\t")
            name = fields[0].removeprefix("!Sample_")
            metadata_rows.setdefault(name, []).append(
                [value.strip().strip('"') for value in fields[1:]]
            )

    def single_row(name: str) -> list[str]:
        # Attributes like the accession number must appear exactly once. Finding
        # none or several means the file's format has changed from what this code
        # was written against, and guessing would risk mismatching patients.
        rows = metadata_rows.get(name)
        if rows is None or len(rows) != 1:
            raise ValueError(f"Expected exactly one GEO Sample_{name} row")
        return rows[0]

    # Pull the attributes that identify each patient, and confirm the file really
    # describes 59 of them with every attribute line the same length. A short line
    # would mean the columns no longer align and values belong to the wrong patient.
    accessions = single_row("geo_accession")
    titles = single_row("title")
    descriptions = single_row("description")
    n_samples = len(accessions)
    characteristic_rows = metadata_rows.get("characteristics_ch1", [])
    if n_samples != EXPECTED_SAMPLES or any(
        len(row) != n_samples for row in characteristic_rows
    ):
        raise ValueError("Unexpected or misaligned GEO sample metadata")

    # Unpack the shared attribute lines. Each entry reads like "age: 67", so
    # splitting on the colon recovers which attribute it is and what its value is,
    # gathering everything known about each patient into one place.
    characteristics: list[dict[str, str]] = [{} for _ in range(n_samples)]
    for row in characteristic_rows:
        for index, text in enumerate(row):
            key, separator, value = text.partition(":")
            characteristics[index][key.strip().lower()] = value.strip() if separator else ""

    # Build the finished record for each patient, one per eventual table row.
    records: list[dict[str, object]] = []
    for index in range(n_samples):
        title = titles[index]
        # Find this patient's short identifier (A1, B3, and so on) — the same code
        # that prefixes every one of their cells. It is normally stated outright;
        # when it is missing, it can be read off the end of the sample's title.
        sample_prefix = descriptions[index].strip()
        if not sample_prefix or sample_prefix.lower() == "unknown":
            match = re.search(r"_([A-Z]\d+)$", title)
            if match is None:
                raise ValueError(f"Cannot recover sample prefix from {title!r}")
            sample_prefix = match.group(1)
        values = characteristics[index]
        diagnosis_original = values.get("diagnosis", "unknown").strip()
        # Settle on one spelling per diagnosis. This only unifies capitalization and
        # wording that already mean the same thing; the untouched original is kept
        # alongside so nothing is lost. A diagnosis never seen before passes through
        # unchanged rather than being forced into one of these categories.
        diagnosis = {
            "healthy": "Healthy",
            "mci": "MCI",
            "ad": "AD",
            "alzheimer's disease": "AD",
            "alzheimers disease": "AD",
            "mci/ad": "MCI/AD",
        }.get(diagnosis_original.lower(), diagnosis_original)
        # Sex is recorded under either of two labels depending on the sample, and in
        # either full or abbreviated form; both are accepted and written one way.
        sex_raw = values.get("sex", values.get("gender", "unknown")).lower()
        sex = {"female": "Female", "f": "Female", "male": "Male", "m": "Male"}.get(
            sex_raw, sex_raw
        )
        records.append(
            {
                "sample_prefix": sample_prefix,
                "patient_id": sample_prefix,
                "gsm": accessions[index],
                "title": title,
                "description": descriptions[index],
                "sex": sex,
                "age": pd.to_numeric(values.get("age", ""), errors="coerce"),
                "diagnosis": diagnosis,
                "diagnosis_original": diagnosis_original,
                "cognitive_group": (
                    "Healthy" if diagnosis == "Healthy" else "Cognitively impaired"
                ),
            }
        )
    # Check the finished table before any of it reaches the cells. Two patients
    # sharing an identifier, or an age that failed to read as a number, would spread
    # silently across tens of thousands of cells and be very hard to trace back.
    metadata = pd.DataFrame(records)
    if not metadata["sample_prefix"].is_unique or not metadata["gsm"].is_unique:
        raise ValueError("GEO sample identifiers are not unique")
    if metadata["age"].isna().any():
        raise ValueError("One or more sample ages could not be parsed")
    # Finally, compare one patient against values confirmed by hand from the GEO
    # record. If the sideways layout above were ever misread — a shift of one
    # position, an attribute landing under the wrong label — the resulting table
    # would still look perfectly well-formed. This single known-good patient is what
    # would catch that, because their details would come back wrong.
    a1 = metadata.loc[metadata["sample_prefix"].eq("A1")].iloc[0]
    if (a1["gsm"], a1["sex"], a1["diagnosis"], float(a1["age"])) != (
        "GSM6016346",
        "Female",
        "Healthy",
        67.0,
    ):
        raise AssertionError("Known A1 GEO metadata did not match expected values")
    return metadata.sort_values("sample_prefix").reset_index(drop=True)


# --- Step 5b: attach sample metadata to cells and save Stage 01 ---------------
# The last step of Stage 01 brings the two halves together: the measurements for
# 70,391 cells, and the table describing the 59 patients they came from.
#
# The link between them is the sample identifier carried in every cell's label. A
# patient's age, sex, and diagnosis get copied onto each of their cells, so that a
# single cell record answers both "what was this cell expressing" and "whose cell
# was it". That is what makes it possible, in later stages, to colour a plot by
# diagnosis or compare healthy against impaired patients.
#
# The copying is deliberately strict: if any cell carries a sample identifier that
# does not appear in the patient table, or any attribute fails to attach to every
# cell, the run stops. A cell silently left without a diagnosis would quietly drop
# out of exactly the comparison this project exists to make.
#
# The finished object is Stage 01's single output and the only file Stage 02 reads.
def attach_metadata_and_write_canonical(paths: PipelinePaths, *, force: bool = False) -> None:
    """Copy each patient's details onto their cells and save Stage 01's output."""

    if paths.canonical_h5ad.exists() and not force:
        validate_h5ad(paths.canonical_h5ad, expected_shape=(EXPECTED_CELLS, EXPECTED_GENES))
        print(f"Using existing canonical AnnData: {paths.canonical_h5ad}")
        return
    require_files([paths.unannotated_h5ad, paths.series_matrix_gzip])
    # Read the patient table, and also save it on its own. It is small and readable
    # in any spreadsheet, which makes it the easiest way to see the study's makeup.
    sample_metadata = parse_geo_sample_metadata(paths.series_matrix_gzip)
    sample_metadata.to_csv(
        paths.metadata_dir / "GSE200164_sample_metadata.csv", index=False
    )

    # Confirm every cell's sample identifier is one of the 59 known patients before
    # copying anything, so a mismatch is reported up front rather than part-way.
    adata = ad.read_h5ad(paths.unannotated_h5ad)
    prefixes = adata.obs["sample_prefix"].astype(str)
    lookup = sample_metadata.set_index("sample_prefix")
    unmatched = sorted(set(prefixes) - set(lookup.index))
    if unmatched:
        raise ValueError(f"Cell prefixes missing from GEO metadata: {unmatched}")
    # Copy each patient attribute down onto every one of that patient's cells.
    for column in (
        "patient_id",
        "gsm",
        "title",
        "sex",
        "age",
        "diagnosis",
        "diagnosis_original",
        "cognitive_group",
    ):
        mapped = prefixes.map(lookup[column])
        if mapped.isna().any():
            raise ValueError(f"Missing cell metadata after mapping {column!r}")
        adata.obs[column] = mapped.to_numpy()
    # Age stays a number so it can be compared and plotted; everything else becomes
    # text, with anything genuinely absent labelled "unknown" rather than left blank.
    adata.obs["age"] = pd.to_numeric(adata.obs["age"], errors="raise").astype(np.float32)
    for column in adata.obs.columns:
        if column != "age":
            adata.obs[column] = adata.obs[column].fillna("unknown").astype(str)

    # Update the object's own description now that it carries patient information,
    # and spell out what each clinical field means, so the file explains itself to
    # someone opening it without this code in front of them.
    adata.uns["dataset"]["expression_status"] = (
        "SoupX-adjusted count-like values; fractional values are expected after ambient-RNA correction."
    )
    adata.uns["processing"]["cell_metadata_attached"] = True
    adata.uns["processing"]["metadata_source"] = "GSE200164 GEO series matrix"
    adata.uns["clinical_variable_definitions"] = {
        "patient_id": "Study sample identifier, equivalent to sample prefix.",
        "age": "Age in years at CSF collection.",
        "sex": "Sex reported in GEO metadata.",
        "diagnosis": "Diagnosis label reported in GEO metadata.",
        "cognitive_group": "Healthy versus cognitively impaired.",
    }
    # Export the per-cell information as a plain table too, so it can be read by
    # tools that cannot open the main file format.
    adata.obs.to_csv(
        paths.metadata_dir / "GSE200164_cell_metadata.csv.gz", compression="gzip"
    )

    # Save Stage 01's finished object and immediately reopen it to confirm the save
    # succeeded, then log how the study broke down between the two patient groups.
    atomic_h5ad_write(adata, paths.canonical_h5ad, force=True)
    validate_h5ad(paths.canonical_h5ad, expected_shape=adata.shape)
    write_json(
        paths.logs_dir / "GSE200164_metadata_integration_summary.json",
        {
            "canonical_h5ad": str(paths.canonical_h5ad),
            "n_cells": adata.n_obs,
            "n_genes": adata.n_vars,
            "n_samples": int(sample_metadata.shape[0]),
            "n_healthy_samples": int((sample_metadata["cognitive_group"] == "Healthy").sum()),
            "n_cognitively_impaired_samples": int(
                (sample_metadata["cognitive_group"] == "Cognitively impaired").sum()
            ),
            "save_validation_passed": True,
        },
    )
    print(f"Canonical AnnData saved to {paths.canonical_h5ad}")


def validate_h5ad(path: Path, *, expected_shape: tuple[int, int]) -> None:
    """Validate a saved H5AD without loading its full expression matrix into RAM."""

    # Inspect the file's structure while leaving the bulk of the data on disk, so
    # checking a multi-gigabyte result costs almost nothing.
    check = ad.read_h5ad(path, backed="r")
    try:
        if check.shape != expected_shape:
            raise AssertionError(f"{path} has shape {check.shape}; expected {expected_shape}")
        if not check.obs_names.is_unique or not check.var_names.is_unique:
            raise AssertionError(f"{path} contains duplicate observation or variable names")
    finally:
        check.file.close()


# --- Stage 01 controller -----------------------------------------------------
# The whole stage in order. Read this function to see what Stage 01 does; read the
# five steps above for how each part works.
#
# Going in: two public files, roughly 1.1 GB compressed.
# Coming out: one object holding 70,391 cells by 19,689 genes, where every cell
# knows which patient it came from and that patient's age, sex, and diagnosis.
def run_data_preparation(
    paths: PipelinePaths,
    *,
    chunk_size: int = 100,
    force: bool = False,
) -> None:
    """Run Stage 01 end to end, reusing a valid finished file when one exists."""

    paths.create_directories()

    # If the finished object is already present and passes inspection, everything
    # below has been done before. Stage 01 takes hours from scratch, so re-running
    # the command after an interruption should cost seconds, not a repeat run.
    if paths.canonical_h5ad.exists() and not force:
        validate_h5ad(
            paths.canonical_h5ad,
            expected_shape=(EXPECTED_CELLS, EXPECTED_GENES),
        )
        print(
            "Stage 01 complete: reusing validated canonical AnnData at "
            f"{paths.canonical_h5ad}"
        )
        return

    # Step 1: fetch the expression measurements and the patient metadata, and prove
    # both arrived undamaged.
    download_source_files(paths, force=force)

    # Step 2: rewrite the huge dense expression file as a compact sparse matrix,
    # keeping only the nonzero measurements, in both row orientations.
    convert_expression_to_sparse(paths, chunk_size=chunk_size, force=force)

    # Step 3: go back to the original download and verify the rewrite changed none
    # of the values. Everything downstream depends on this holding.
    audit_sparse_conversion(paths)

    # Step 4: bind the measurements, gene descriptions, and cell labels into one
    # object, with clinical fields present but not yet filled.
    build_unannotated_anndata(paths, force=force)

    # Step 5: read the patient table, copy each patient's details onto their cells,
    # and save the result — the file Stage 02 will pick up.
    attach_metadata_and_write_canonical(paths, force=force)


# --- Command-line interface -------------------------------------------------
# What can be typed after the filename when running this stage. Only two things are
# adjustable: where the generated files go, and whether to redo finished work. The
# scientific settings are not options — they live in the steps above, where they
# can be read in context.
def parse_args() -> argparse.Namespace:
    """Read Stage 01 command-line options."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--workspace",
        default=default_workspace(),
        help="Generated-data workspace (default: checkout-local Age_trajectory)",
    )
    parser.add_argument("--chunk-size", type=int, default=100, help="Gene rows per CSV chunk")
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-download and overwrite existing generated outputs",
    )
    return parser.parse_args()


def main() -> None:
    """Run the complete first stage from this one file."""

    # Example: python scripts/01_download_audit_build_anndata.py --chunk-size 100
    args = parse_args()
    paths = PipelinePaths.from_workspace(args.workspace)
    run_data_preparation(paths, chunk_size=args.chunk_size, force=args.force)


if __name__ == "__main__":
    main()
