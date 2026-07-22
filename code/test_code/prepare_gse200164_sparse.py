#!/usr/bin/env python3
"""Convert the published GSE200164 table to a lossless sparse matrix cache.

GEO calls ``GSE200164_counts.csv`` a raw-count matrix, but its values are the
authors' SoupX-adjusted counts. SoupX subtraction legitimately produces
fractional counts. The converter therefore stores floating-point values without
rounding them to integers. The output is a SciPy CSR ``.npz`` matrix, which is
substantially faster and smaller than a multi-gigabyte text Matrix Market file.
"""

# Allow modern type annotations while retaining compatibility with Python 3.9.
from __future__ import annotations

# argparse defines the command-line interface used by the R analysis script.
import argparse
# json records provenance and matrix statistics for cache validation.
import json
# os provides atomic file replacement after every output is complete.
import os
# Path keeps all source and output paths explicit and platform-independent.
from pathlib import Path
# tempfile creates an isolated staging directory beside the final cache.
import tempfile

# NumPy performs numeric validation and detects fractional SoupX values.
import numpy as np
# pandas streams the 4+ GB wide CSV without loading the full dense table.
import pandas as pd
# SciPy constructs and serializes a standard sparse gene-by-cell matrix.
from scipy import sparse


# These are the six gene-annotation columns that precede the cell barcodes.
GENE_METADATA_COLUMNS = ["Unnamed: 0", "chromosome", "start", "end", "width", "strand"]


def parse_args() -> argparse.Namespace:
    """Read and validate command-line settings for the conversion stage."""
    # Create the parser shown by ``--help``.
    parser = argparse.ArgumentParser(
        description="Create a lossless SciPy sparse cache from GSE200164_counts.csv."
    )
    # Point to the authors' SoupX-adjusted, doublet/QC-filtered count matrix.
    parser.add_argument("--counts", required=True, type=Path)
    # Point to the sample-level metadata assembled from GEO.
    parser.add_argument("--metadata", required=True, type=Path)
    # Choose the directory that will hold the reusable sparse cache.
    parser.add_argument("--out-dir", required=True, type=Path)
    # Keep only a modest number of dense gene rows in memory at one time.
    parser.add_argument("--chunksize", default=32, type=int)
    # Recreate a valid cache when the user explicitly requests it.
    parser.add_argument("--force", action="store_true")
    # Return the parsed namespace to the main conversion function.
    return parser.parse_args()


def make_unique(values: list[str]) -> list[str]:
    """Make gene symbols unique using the same suffix convention as R."""
    # Track how many times each original symbol has appeared.
    seen: dict[str, int] = {}
    # Accumulate unique symbols in exactly the source row order.
    unique: list[str] = []
    # Process one source symbol at a time so the first occurrence is unchanged.
    for value in values:
        # Convert missing-looking or numeric labels to a stable string.
        value = str(value)
        # Read the number of earlier occurrences of this symbol.
        count = seen.get(value, 0)
        # Update the occurrence count for subsequent duplicate rows.
        seen[value] = count + 1
        # Add ``.2``, ``.3``, and so on only to later duplicate symbols.
        unique.append(value if count == 0 else f"{value}.{count + 1}")
    # Return names suitable for a Seurat feature-name vector.
    return unique


def source_signature(path: Path) -> dict[str, object]:
    """Describe an input file well enough to reject a stale sparse cache."""
    # Resolve the path so equivalent relative spellings compare identically.
    resolved = path.resolve()
    # Read size and nanosecond modification time without scanning the large file.
    stat = resolved.stat()
    # Return only stable, JSON-serializable provenance fields.
    return {
        "path": str(resolved),
        "size_bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }


def cache_is_current(
    manifest_path: Path,
    expected_outputs: list[Path],
    counts_path: Path,
    metadata_path: Path,
) -> bool:
    """Return True only when every output matches the current source files."""
    # A missing output means the previous conversion was absent or interrupted.
    if not all(path.exists() for path in expected_outputs):
        return False
    # Treat an unreadable manifest as invalid instead of risking silent reuse.
    try:
        manifest = json.loads(manifest_path.read_text())
    except (OSError, json.JSONDecodeError):
        return False
    # Compare both large counts and small metadata inputs to their saved stamps.
    return (
        manifest.get("counts_source") == source_signature(counts_path)
        and manifest.get("metadata_source") == source_signature(metadata_path)
        and manifest.get("matrix_format") == "scipy_csr_npz_float64"
    )


def build_cell_metadata(
    barcodes: list[str], sample_metadata_path: Path
) -> pd.DataFrame:
    """Expand one GEO metadata row per sample to one row per cell barcode."""
    # Read the compact sample-level table created from the GEO family record.
    sample_metadata = pd.read_csv(sample_metadata_path)
    # Fail early if the cell-to-sample join key is unavailable.
    if "sample_id" not in sample_metadata.columns:
        raise ValueError(f"Expected a sample_id column in {sample_metadata_path}")
    # Require unique sample IDs to prevent accidental one-to-many row expansion.
    if sample_metadata["sample_id"].duplicated().any():
        raise ValueError("sample_id values must be unique in the metadata table")

    # Start a one-row-per-cell table in count-matrix column order.
    cell_metadata = pd.DataFrame({"cell": barcodes})
    # The prefix before the first underscore is the authors' sample identifier.
    cell_metadata["sample_id"] = cell_metadata["cell"].str.split("_", n=1).str[0]
    # Join all diagnosis, sex, age, GSM, and expected cell-count fields.
    cell_metadata = cell_metadata.merge(
        sample_metadata, on="sample_id", how="left", validate="many_to_one", sort=False
    )
    # Restore matrix order explicitly because downstream row names must align.
    cell_metadata = cell_metadata.set_index("cell").loc[barcodes].reset_index()
    # Match Seurat's original identity name used in the authors' integration.
    cell_metadata["orig.ident"] = cell_metadata["sample_id"]

    # Detect cell prefixes that did not match any sample-level GEO record.
    missing = sorted(
        cell_metadata.loc[cell_metadata["gsm"].isna(), "sample_id"].unique()
    )
    # Stop rather than allow cells with incomplete covariates into integration.
    if missing:
        raise ValueError(
            "Cell barcodes had sample prefixes missing from metadata: "
            + ", ".join(missing[:20])
        )

    # Confirm that the published per-sample cell totals match matrix barcodes.
    observed = cell_metadata.groupby("sample_id", sort=False).size()
    expected = sample_metadata.set_index("sample_id")["n_cells"]
    # Find any sample whose expected and observed counts are not identical.
    mismatched = [
        sample
        for sample in expected.index
        if int(observed.get(sample, 0)) != int(expected.loc[sample])
    ]
    # Reject a mismatched matrix/metadata pair before expensive analysis begins.
    if mismatched:
        details = [
            f"{sample}: observed={int(observed.get(sample, 0))}, "
            f"expected={int(expected.loc[sample])}"
            for sample in mismatched[:20]
        ]
        raise ValueError("Per-sample cell-count mismatch: " + "; ".join(details))

    # Return the validated one-row-per-cell metadata table.
    return cell_metadata


def main() -> None:
    """Stream, validate, sparsify, and atomically cache the GEO matrix."""
    # Parse the settings supplied by the R pipeline or command line.
    args = parse_args()
    # Require a positive chunk size because pandas cannot stream zero rows.
    if args.chunksize < 1:
        raise ValueError("--chunksize must be at least 1")
    # Resolve inputs before recording their provenance in the manifest.
    args.counts = args.counts.resolve()
    args.metadata = args.metadata.resolve()
    # Ensure the destination exists before creating its staging directory.
    args.out_dir.mkdir(parents=True, exist_ok=True)

    # Name every reusable cache component consumed by the R workflow.
    matrix_path = args.out_dir / "GSE200164_counts_csr.npz"
    genes_path = args.out_dir / "GSE200164_genes.txt"
    barcodes_path = args.out_dir / "GSE200164_barcodes.txt"
    cell_metadata_path = args.out_dir / "GSE200164_cell_metadata.csv"
    manifest_path = args.out_dir / "GSE200164_sparse_manifest.json"
    done_path = args.out_dir / ".sparse_cache_complete"
    # Require the sentinel and all data files before considering reuse.
    expected_outputs = [
        matrix_path,
        genes_path,
        barcodes_path,
        cell_metadata_path,
        manifest_path,
        done_path,
    ]

    # Reuse a cache only when its inputs and lossless format still match.
    if not args.force and cache_is_current(
        manifest_path, expected_outputs, args.counts, args.metadata
    ):
        print(f"Sparse cache is current: {args.out_dir}", flush=True)
        return

    # Read only the wide header to discover cell barcodes inexpensively.
    header = pd.read_csv(args.counts, nrows=0)
    # Preserve the exact source-column order throughout the conversion.
    columns = list(header.columns)
    # Confirm the expected GEO gene-annotation schema before selecting counts.
    missing_gene_meta = [
        column for column in GENE_METADATA_COLUMNS if column not in columns
    ]
    # A schema mismatch could otherwise make annotations look like cell counts.
    if missing_gene_meta:
        raise ValueError(f"Missing count metadata columns: {missing_gene_meta}")
    # Every non-annotation column represents one cell barcode.
    count_columns = [
        column for column in columns if column not in GENE_METADATA_COLUMNS
    ]
    # Duplicate barcodes would make cell metadata and matrix columns ambiguous.
    if len(count_columns) != len(set(count_columns)):
        raise ValueError("Count-matrix cell barcodes are not unique")

    # Validate and expand sample metadata before scanning the 4+ GB matrix body.
    cell_metadata = build_cell_metadata(count_columns, args.metadata)
    # Accumulate sparse row blocks rather than the impossible dense full matrix.
    sparse_blocks: list[sparse.csr_matrix] = []
    # Accumulate feature names in the same order as sparse matrix rows.
    genes: list[str] = []
    # Track whether values prove this is SoupX-adjusted rather than integer UMI data.
    fractional_values_seen = False
    # Track extrema for an auditable cache manifest.
    minimum_nonzero = np.inf
    maximum_value = 0.0
    # Count processed rows for progress reports and final shape validation.
    n_genes = 0
    # Count stored values without rescanning sparse blocks later.
    nnz = 0

    # Stream manageable gene blocks from the uncompressed GEO CSV.
    for chunk in pd.read_csv(args.counts, chunksize=args.chunksize):
        # Append gene symbols before discarding the annotation columns.
        genes.extend(chunk["Unnamed: 0"].astype(str).tolist())
        # Convert only cell columns to float64 so SoupX decimals remain unchanged.
        dense_block = chunk[count_columns].to_numpy(dtype=np.float64, copy=False)
        # Reject malformed, infinite, or missing numeric values immediately.
        if not np.isfinite(dense_block).all():
            raise ValueError(f"Non-finite count encountered near gene row {n_genes + 1}")
        # SoupX-adjusted expression should never contain negative abundance.
        if np.any(dense_block < 0):
            raise ValueError(f"Negative count encountered near gene row {n_genes + 1}")
        # Record evidence that integer conversion would corrupt this dataset.
        nonzero_values = dense_block[dense_block != 0]
        # Update summary statistics only when this block contains expressed genes.
        if nonzero_values.size:
            fractional_values_seen |= bool(
                np.any(~np.isclose(nonzero_values, np.rint(nonzero_values)))
            )
            minimum_nonzero = min(minimum_nonzero, float(nonzero_values.min()))
            maximum_value = max(maximum_value, float(nonzero_values.max()))
        # Convert this gene block to CSR without changing its numeric values.
        sparse_block = sparse.csr_matrix(dense_block)
        # Remove any explicit zeros defensively before serialization.
        sparse_block.eliminate_zeros()
        # Save the block for one efficient final sparse concatenation.
        sparse_blocks.append(sparse_block)
        # Advance matrix and progress counters using the block's actual shape.
        n_genes += sparse_block.shape[0]
        nnz += sparse_block.nnz
        # Report progress frequently enough to distinguish work from a stalled run.
        print(
            f"processed {n_genes:,} genes; nonzero values so far: {nnz:,}",
            flush=True,
        )

    # Fail clearly if a truncated or empty source produced no feature rows.
    if not sparse_blocks:
        raise ValueError("The count matrix contained no gene rows")
    # Stack row blocks into the final gene-by-cell CSR matrix.
    counts = sparse.vstack(sparse_blocks, format="csr")
    # Release the list's references before compressed serialization peaks in memory.
    sparse_blocks.clear()
    # Verify the final shape against independently counted genes and barcodes.
    if counts.shape != (n_genes, len(count_columns)):
        raise RuntimeError(
            f"Sparse matrix shape {counts.shape} did not match "
            f"({n_genes}, {len(count_columns)})"
        )
    # Make duplicate symbols unique only after preserving every matrix row.
    unique_genes = make_unique(genes)

    # Build a provenance record before staging any final output names.
    manifest = {
        "matrix_format": "scipy_csr_npz_float64",
        "matrix_orientation": "genes_by_cells",
        "n_genes": n_genes,
        "n_cells": len(count_columns),
        "nnz": int(counts.nnz),
        "fractional_values_seen": fractional_values_seen,
        "minimum_nonzero": None if not np.isfinite(minimum_nonzero) else minimum_nonzero,
        "maximum_value": maximum_value,
        "counts_source": source_signature(args.counts),
        "metadata_source": source_signature(args.metadata),
    }

    # Stage every output beside its destination so final renames are atomic.
    with tempfile.TemporaryDirectory(dir=args.out_dir, prefix="sparse-build-") as tmp:
        # Convert the temporary directory string into convenient Path objects.
        staging_dir = Path(tmp)
        # Define staging paths with the same base names as final cache outputs.
        staged_matrix = staging_dir / matrix_path.name
        staged_genes = staging_dir / genes_path.name
        staged_barcodes = staging_dir / barcodes_path.name
        staged_metadata = staging_dir / cell_metadata_path.name
        staged_manifest = staging_dir / manifest_path.name
        staged_done = staging_dir / done_path.name

        # Save a compressed standard SciPy sparse archive without integer coercion.
        sparse.save_npz(staged_matrix, counts, compressed=True)
        # Write feature names one per line in sparse row order.
        staged_genes.write_text("\n".join(unique_genes) + "\n")
        # Write cell barcodes one per line in sparse column order.
        staged_barcodes.write_text("\n".join(count_columns) + "\n")
        # Write aligned cell metadata without an extraneous pandas index.
        cell_metadata.to_csv(staged_metadata, index=False)
        # Pretty-print provenance so it is easy to audit outside Python.
        staged_manifest.write_text(json.dumps(manifest, indent=2) + "\n")
        # Create the completion sentinel only after all staged writes succeed.
        staged_done.write_text("ok\n")

        # Move data files into place before publishing the completion sentinel.
        for staged, final in [
            (staged_matrix, matrix_path),
            (staged_genes, genes_path),
            (staged_barcodes, barcodes_path),
            (staged_metadata, cell_metadata_path),
            (staged_manifest, manifest_path),
        ]:
            os.replace(staged, final)
        # Publish the sentinel last so interrupted builds are never reused.
        os.replace(staged_done, done_path)

    # Print a concise completion record for both humans and calling scripts.
    print(
        f"Wrote {counts.shape[0]:,} genes x {counts.shape[1]:,} cells "
        f"with {counts.nnz:,} nonzero values to {matrix_path}",
        flush=True,
    )


# Run conversion only when the file is executed, not when imported by tests.
if __name__ == "__main__":
    main()
