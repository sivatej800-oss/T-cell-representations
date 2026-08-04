"""Implementation of Stage 01: build the canonical GSE200164 AnnData object.

The source expression file looks like a dense CSV, but loading it all at once would
need excessive memory. This module reads a small group of genes at a time, builds a
sparse matrix, verifies that conversion against the source, and attaches GEO sample
metadata before saving the canonical ``.h5ad`` file.
"""

from __future__ import annotations

import csv
import gc
import gzip
import random
import re
import shutil
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from .common import (
    atomic_h5ad_write,
    download_with_resume,
    require_files,
    sha256_file,
    validate_gzip,
    write_json,
)
from .config import COUNTS_URL, GEO_ACCESSION, SERIES_MATRIX_URL, PipelinePaths


# The first six CSV columns describe each gene. Every remaining column is one cell.
ANNOTATION_COLUMNS = ["gene_symbol", "chromosome", "start", "end", "width", "strand"]

# Known public-dataset dimensions provide early protection against partial downloads,
# unexpected GEO revisions, or an accidentally selected input file.
EXPECTED_CELLS = 70_391
EXPECTED_GENES = 19_689
EXPECTED_SAMPLES = 59


def download_source_files(paths: PipelinePaths, *, force: bool = False) -> None:
    """Download, integrity-check, and record provenance for both GEO source files."""

    # The large counts file is validated end to end before any conversion begins.
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

    # The much smaller series matrix contains sample ages, diagnoses, and sexes.
    download_with_resume(SERIES_MATRIX_URL, paths.series_matrix_gzip, force=force)
    validate_gzip(paths.series_matrix_gzip)


def _read_source_header(path: Path) -> list[str]:
    """Read only the first CSV row, which contains annotations and cell barcodes."""

    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        return next(csv.reader(handle))


def convert_expression_to_sparse(
    paths: PipelinePaths,
    *,
    chunk_size: int = 100,
    force: bool = False,
) -> None:
    """Convert the dense GEO CSV stream into two sparse matrix orientations.

    The GEO file is genes-by-cells, while AnnData expects cells-by-genes. Both forms
    are saved so the conversion is auditable and later steps never repeat this work.
    """

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
    # Split the header into six gene-annotation labels and 70,391 cell barcodes.
    header = _read_source_header(paths.counts_gzip)
    cell_barcodes = header[len(ANNOTATION_COLUMNS) :]
    # This quick first pass counts rows without materializing the expression values.
    with gzip.open(paths.counts_gzip, "rt", encoding="utf-8") as handle:
        n_genes = sum(1 for _ in handle) - 1

    if len(cell_barcodes) != EXPECTED_CELLS or n_genes != EXPECTED_GENES:
        raise ValueError(
            "Unexpected GEO dimensions: "
            f"{len(cell_barcodes):,} cells and {n_genes:,} genes"
        )
    paths.barcodes_txt.write_text(
        "".join(f"{barcode}\n" for barcode in cell_barcodes), encoding="utf-8"
    )

    # pandas now streams only ``chunk_size`` gene rows at a time. Smaller chunks use
    # less peak RAM but require more loop iterations.
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

    for chunk_number, chunk in enumerate(chunks, start=1):
        # Keep gene annotations as a table and expression values as float32. Float32
        # is intentional: SoupX-adjusted values are fractional and must not be rounded.
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

        # CSR stores only nonzero values and their positions, avoiding a huge dense
        # matrix in memory and on disk.
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
        # Explicit cleanup helps keep peak memory stable while streaming this file.
        del values, nonzero, chunk
        gc.collect()
    print()

    if processed != n_genes:
        raise RuntimeError(f"Processed {processed:,} genes; expected {n_genes:,}")
    # Stack gene chunks in source order, then transpose once into AnnData's expected
    # cells-by-genes orientation.
    genes_by_cells = sparse.vstack(sparse_blocks, format="csr", dtype=np.float32)
    cells_by_genes = genes_by_cells.transpose().tocsr().astype(np.float32)
    cells_by_genes.sort_indices()
    cells_by_genes.eliminate_zeros()
    if cells_by_genes.shape != (EXPECTED_CELLS, EXPECTED_GENES):
        raise RuntimeError(f"Unexpected sparse matrix shape: {cells_by_genes.shape}")

    sparse.save_npz(paths.genes_by_cells_npz, genes_by_cells, compressed=True)
    sparse.save_npz(paths.cells_by_genes_npz, cells_by_genes, compressed=True)
    pd.concat(metadata_blocks, ignore_index=True).to_csv(paths.genes_csv, index=False)

    # Store summary statistics as a small audit record rather than recomputing them
    # whenever someone wants to inspect the conversion.
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


def audit_sparse_conversion(paths: PipelinePaths) -> None:
    """Verify selected source values survived sparse conversion unchanged.

    First/last indices cover boundaries. A fixed random seed selects additional genes
    and cells reproducibly. Full cell vectors are also compared for selected genes.
    """

    require_files(
        [
            paths.counts_gzip,
            paths.cells_by_genes_npz,
            paths.genes_csv,
            paths.barcodes_txt,
        ]
    )
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

    # A local seeded generator makes the audit repeatable without changing NumPy's or
    # Python's global random state.
    generator = random.Random(42)
    gene_indices = sorted({0, len(genes) - 1, *generator.sample(range(len(genes)), 12)})
    cell_indices = sorted(
        {0, len(barcodes) - 1, *generator.sample(range(len(barcodes)), 20)}
    )
    # Stream the source again, retaining only the small set of gene rows being audited.
    selected_rows: dict[int, list[str]] = {}
    with gzip.open(paths.counts_gzip, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        next(reader)
        for gene_index, row in enumerate(reader):
            if gene_index in gene_indices:
                selected_rows[gene_index] = row
            if gene_index >= gene_indices[-1]:
                break

    pair_records: list[dict[str, object]] = []
    vector_records: list[dict[str, object]] = []
    for gene_index in gene_indices:
        row = selected_rows[gene_index]
        if row[0] != genes[gene_index]:
            raise ValueError(f"Gene-order mismatch at index {gene_index}")
        # Compare the entire selected gene across all cells before checking individual
        # gene/cell pairs. Both comparisons use the same float32 precision as storage.
        source_vector = np.asarray(row[len(ANNOTATION_COLUMNS) :], dtype=np.float32)
        sparse_vector = expression[:, gene_index].toarray().ravel().astype(np.float32)
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

    pair_table = pd.DataFrame(pair_records)
    vector_table = pd.DataFrame(vector_records)
    pair_table.to_csv(paths.tables_dir / "GSE200164_sparse_conversion_audit.csv", index=False)
    vector_table.to_csv(
        paths.tables_dir / "GSE200164_sparse_full_vector_audit.csv", index=False
    )
    if not pair_table["match"].all() or not vector_table["all_cells_match"].all():
        raise AssertionError("Sparse conversion audit failed")
    print(f"Sparse audit passed: {len(pair_table):,} pairs and {len(vector_table)} vectors")


def build_unannotated_anndata(paths: PipelinePaths, *, force: bool = False) -> None:
    """Build AnnData before clinical sample metadata is attached.

    AnnData stores cells in ``obs`` (rows), genes in ``var`` (columns), and the sparse
    expression matrix in ``X``. Placeholder clinical fields make the schema explicit
    until the GEO series matrix is joined in the next step.
    """

    if paths.unannotated_h5ad.exists() and not force:
        print(f"Using existing AnnData file: {paths.unannotated_h5ad}")
        return
    require_files([paths.cells_by_genes_npz, paths.genes_csv, paths.barcodes_txt])
    expression = sparse.load_npz(paths.cells_by_genes_npz).tocsr().astype(np.float32)
    genes = pd.read_csv(paths.genes_csv, low_memory=False)
    barcodes = paths.barcodes_txt.read_text(encoding="utf-8").splitlines()
    if expression.shape != (len(barcodes), len(genes)):
        raise ValueError("AnnData inputs have inconsistent dimensions")

    for column in ("gene_symbol", "chromosome", "strand"):
        genes[column] = genes[column].fillna("unknown").astype(str)
    for column in ("start", "end", "width"):
        genes[column] = pd.to_numeric(genes[column], errors="coerce")

    # A source barcode such as ``A1_AAAC...`` contains both the study sample prefix
    # and the within-sample 10x barcode.
    split_barcodes = pd.Series(barcodes).str.split("_", n=1, expand=True)
    if split_barcodes.shape[1] != 2:
        raise ValueError("Cell barcodes do not contain sample prefixes")
    obs = pd.DataFrame(index=pd.Index(barcodes, name="cell_barcode"))
    obs["sample_prefix"] = split_barcodes.iloc[:, 0].to_numpy(dtype=str)
    obs["tenx_barcode"] = split_barcodes.iloc[:, 1].to_numpy(dtype=str)
    obs["patient_id"] = "unknown"
    obs["diagnosis"] = "unknown"
    obs["sex"] = "unknown"
    obs["cell_type"] = "unknown"
    obs["age"] = np.full(len(obs), np.nan, dtype=np.float32)

    # AnnData requires unique variable names. Original symbols remain in ``var`` even
    # when a duplicate suffix is added to the index.
    var = genes.copy()
    var.index = pd.Index(var["gene_symbol"].astype(str), name="gene_id")
    adata = ad.AnnData(X=expression, obs=obs, var=var)
    adata.var_names_make_unique(join="-duplicate-")
    # ``uns`` holds dataset-level provenance that does not belong to one cell or gene.
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


def parse_geo_sample_metadata(series_matrix: Path) -> pd.DataFrame:
    """Turn GEO's repeated header rows into one clean record per biological sample."""

    # GEO encodes each field as a tab-separated !Sample_* row. Some field names, such
    # as characteristics_ch1, occur multiple times for age, diagnosis, and sex.
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
        rows = metadata_rows.get(name)
        if rows is None or len(rows) != 1:
            raise ValueError(f"Expected exactly one GEO Sample_{name} row")
        return rows[0]

    accessions = single_row("geo_accession")
    titles = single_row("title")
    descriptions = single_row("description")
    n_samples = len(accessions)
    characteristic_rows = metadata_rows.get("characteristics_ch1", [])
    if n_samples != EXPECTED_SAMPLES or any(
        len(row) != n_samples for row in characteristic_rows
    ):
        raise ValueError("Unexpected or misaligned GEO sample metadata")

    # Reorganize repeated ``key: value`` characteristics from row-oriented GEO text
    # into one dictionary per sample.
    characteristics: list[dict[str, str]] = [{} for _ in range(n_samples)]
    for row in characteristic_rows:
        for index, text in enumerate(row):
            key, separator, value = text.partition(":")
            characteristics[index][key.strip().lower()] = value.strip() if separator else ""

    records: list[dict[str, object]] = []
    for index in range(n_samples):
        title = titles[index]
        sample_prefix = descriptions[index].strip()
        if not sample_prefix or sample_prefix.lower() == "unknown":
            match = re.search(r"_([A-Z]\d+)$", title)
            if match is None:
                raise ValueError(f"Cannot recover sample prefix from {title!r}")
            sample_prefix = match.group(1)
        values = characteristics[index]
        diagnosis_original = values.get("diagnosis", "unknown").strip()
        # Normalize spelling/capitalization only; do not infer a diagnosis not present
        # in the public metadata.
        diagnosis = {
            "healthy": "Healthy",
            "mci": "MCI",
            "ad": "AD",
            "alzheimer's disease": "AD",
            "alzheimers disease": "AD",
            "mci/ad": "MCI/AD",
        }.get(diagnosis_original.lower(), diagnosis_original)
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
    metadata = pd.DataFrame(records)
    if not metadata["sample_prefix"].is_unique or not metadata["gsm"].is_unique:
        raise ValueError("GEO sample identifiers are not unique")
    if metadata["age"].isna().any():
        raise ValueError("One or more sample ages could not be parsed")
    # A known sample acts as a sentinel for column shifts or parsing mistakes.
    a1 = metadata.loc[metadata["sample_prefix"].eq("A1")].iloc[0]
    if (a1["gsm"], a1["sex"], a1["diagnosis"], float(a1["age"])) != (
        "GSM6016346",
        "Female",
        "Healthy",
        67.0,
    ):
        raise AssertionError("Known A1 GEO metadata did not match expected values")
    return metadata.sort_values("sample_prefix").reset_index(drop=True)


def attach_metadata_and_write_canonical(paths: PipelinePaths, *, force: bool = False) -> None:
    """Map sample metadata onto every cell and save the canonical Stage 01 result."""

    if paths.canonical_h5ad.exists() and not force:
        validate_h5ad(paths.canonical_h5ad, expected_shape=(EXPECTED_CELLS, EXPECTED_GENES))
        print(f"Using existing canonical AnnData: {paths.canonical_h5ad}")
        return
    require_files([paths.unannotated_h5ad, paths.series_matrix_gzip])
    sample_metadata = parse_geo_sample_metadata(paths.series_matrix_gzip)
    sample_metadata.to_csv(
        paths.metadata_dir / "GSE200164_sample_metadata.csv", index=False
    )

    adata = ad.read_h5ad(paths.unannotated_h5ad)
    # Map through the short A1-H8-style prefix already recovered from each barcode.
    # This repeats sample-level values for every cell belonging to that sample.
    prefixes = adata.obs["sample_prefix"].astype(str)
    lookup = sample_metadata.set_index("sample_prefix")
    unmatched = sorted(set(prefixes) - set(lookup.index))
    if unmatched:
        raise ValueError(f"Cell prefixes missing from GEO metadata: {unmatched}")
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
    adata.obs["age"] = pd.to_numeric(adata.obs["age"], errors="raise").astype(np.float32)
    for column in adata.obs.columns:
        if column != "age":
            adata.obs[column] = adata.obs[column].fillna("unknown").astype(str)

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
    adata.obs.to_csv(
        paths.metadata_dir / "GSE200164_cell_metadata.csv.gz", compression="gzip"
    )

    # The canonical filename is the supported Stage 02 input. The annotated filename
    # is retained as a compatibility copy for earlier project work.
    atomic_h5ad_write(adata, paths.canonical_h5ad, force=True)
    validate_h5ad(paths.canonical_h5ad, expected_shape=adata.shape)
    if paths.annotated_h5ad.exists():
        if force:
            paths.annotated_h5ad.unlink()
    if not paths.annotated_h5ad.exists():
        shutil.copy2(paths.canonical_h5ad, paths.annotated_h5ad)
    write_json(
        paths.logs_dir / "GSE200164_metadata_integration_summary.json",
        {
            "canonical_h5ad": str(paths.canonical_h5ad),
            "annotated_h5ad": str(paths.annotated_h5ad),
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

    # backed="r" reads metadata from disk while leaving the large matrix on disk.
    check = ad.read_h5ad(path, backed="r")
    try:
        if check.shape != expected_shape:
            raise AssertionError(f"{path} has shape {check.shape}; expected {expected_shape}")
        if not check.obs_names.is_unique or not check.var_names.is_unique:
            raise AssertionError(f"{path} contains duplicate observation or variable names")
    finally:
        check.file.close()


def run_data_preparation(
    paths: PipelinePaths,
    *,
    chunk_size: int = 100,
    force: bool = False,
) -> None:
    """Run Stage 01 in dependency order, reusing a valid final file when possible."""

    paths.create_directories()
    # The canonical H5AD is the terminal Stage 01 artifact. If it validates, every
    # expensive upstream operation can be skipped safely.
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
    # From-scratch order: source -> sparse conversion -> audit -> AnnData -> metadata.
    download_source_files(paths, force=force)
    convert_expression_to_sparse(paths, chunk_size=chunk_size, force=force)
    audit_sparse_conversion(paths)
    build_unannotated_anndata(paths, force=force)
    attach_metadata_and_write_canonical(paths, force=force)
