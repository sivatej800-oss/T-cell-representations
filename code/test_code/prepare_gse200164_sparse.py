#!/usr/bin/env python3
"""Convert the wide GSE200164 count CSV into sparse Matrix Market files.

The GEO count matrix is gene-by-cell in a very wide CSV. Reading it directly
inside R is expensive, so this script streams row chunks and writes a reusable
sparse cache for Seurat.
"""

# Let Python understand modern type hints such as list[str].
from __future__ import annotations

# argparse reads command-line options passed from the R script.
import argparse
# csv writes space-separated sparse-matrix coordinate rows.
import csv
# Path represents file and folder locations safely.
from pathlib import Path
# tempfile creates a temporary coordinate file during conversion.
import tempfile

# numpy efficiently finds nonzero numeric values.
import numpy as np
# pandas reads the CSV in manageable blocks of rows.
import pandas as pd


# The source CSV has these gene-description columns before its many cell columns.
GENE_METADATA_COLUMNS = ["Unnamed: 0", "chromosome", "start", "end", "width", "strand"]


def parse_args() -> argparse.Namespace:
    """Read settings provided on the command line."""
    parser = argparse.ArgumentParser()
    # The original wide gene-by-cell count table.
    parser.add_argument("--counts", required=True, type=Path)
    # Sample-level metadata that will be attached to each cell.
    parser.add_argument("--metadata", required=True, type=Path)
    # Folder where the reusable sparse cache will be created.
    parser.add_argument("--out-dir", required=True, type=Path)
    # Number of gene rows to keep in memory at one time.
    parser.add_argument("--chunksize", default=128, type=int)
    # Rebuild outputs even when an earlier completed cache is available.
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def make_unique(values: list[str]) -> list[str]:
    # Matrix row names must be unique; preserve the first occurrence and suffix
    # later duplicates using an R-like `.2`, `.3`, ... convention.
    seen: dict[str, int] = {}
    unique: list[str] = []
    for value in values:
        value = str(value)
        count = seen.get(value, 0)
        seen[value] = count + 1
        unique.append(value if count == 0 else f"{value}.{count + 1}")
    return unique


def write_lines(path: Path, values: list[str]) -> None:
    """Write names to a text file, one name per line."""
    path.write_text("\n".join(values) + "\n")


def write_cell_metadata(barcodes: list[str], sample_metadata_path: Path, output_path: Path) -> None:
    """Make a metadata row for each cell barcode in the count matrix."""
    # Read the original sample-level metadata CSV.
    sample_metadata = pd.read_csv(sample_metadata_path)
    if "sample_id" not in sample_metadata.columns:
        raise ValueError(f"Expected a sample_id column in {sample_metadata_path}")

    # Start a new table with one row per cell barcode.
    cell_metadata = pd.DataFrame({"cell": barcodes})
    # GSE200164 barcodes are prefixed with their sample ID, which provides the
    # join key needed to attach GEO-level sample metadata to every cell.
    cell_metadata["sample_id"] = cell_metadata["cell"].str.split("_", n=1).str[0]
    # Add matching sample information to each cell. A left join keeps every cell
    # barcode, even temporarily if its sample information is missing.
    cell_metadata = cell_metadata.merge(sample_metadata, on="sample_id", how="left")

    # Check for missing joins before writing output, rather than allowing a
    # downstream analysis to proceed with incomplete metadata.
    missing = sorted(cell_metadata.loc[cell_metadata["gsm"].isna(), "sample_id"].unique())
    if missing:
        raise ValueError(
            "Some cell barcodes had sample prefixes missing from metadata: "
            + ", ".join(missing[:20])
        )

    # Write the completed per-cell metadata table without an extra row-number column.
    cell_metadata.to_csv(output_path, index=False)


def main() -> None:
    """Convert the wide input CSV into files that R can load efficiently."""
    # Read command-line settings and make the cache folder if needed.
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    # These paths are the four data files in the sparse cache.
    matrix_path = args.out_dir / "GSE200164_counts.mtx"
    genes_path = args.out_dir / "GSE200164_genes.txt"
    barcodes_path = args.out_dir / "GSE200164_barcodes.txt"
    cell_metadata_path = args.out_dir / "GSE200164_cell_metadata.csv"
    # This small sentinel file means all cache outputs finished writing.
    done_path = args.out_dir / ".sparse_cache_complete"

    # The completion sentinel prevents a partially written cache from being
    # treated as usable after an interrupted conversion.
    expected_outputs = [matrix_path, genes_path, barcodes_path, cell_metadata_path, done_path]
    if not args.force and all(path.exists() for path in expected_outputs):
        print(f"Sparse cache already exists in {args.out_dir}")
        return

    # Read only the column names first. This avoids loading the very large CSV.
    header = pd.read_csv(args.counts, nrows=0)
    columns = list(header.columns)
    # Confirm that the input has the expected gene-description columns.
    missing_gene_meta = [column for column in GENE_METADATA_COLUMNS if column not in columns]
    if missing_gene_meta:
        raise ValueError(f"Missing expected count metadata columns: {missing_gene_meta}")

    # Every remaining column is a cell barcode with one count per gene.
    count_columns = [column for column in columns if column not in GENE_METADATA_COLUMNS]
    # Keep track of matrix dimensions and the number of stored nonzero entries.
    n_cells = len(count_columns)
    n_genes = 0
    nnz = 0
    genes: list[str] = []

    # Save column names and create matching metadata before streaming counts.
    write_lines(barcodes_path, count_columns)
    write_cell_metadata(count_columns, args.metadata, cell_metadata_path)

    # Matrix Market requires dimensions and total nonzero entries in its header.
    # Stream coordinates to a temporary file first, then prepend that header
    # once the complete matrix shape and nnz are known.
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", newline="", dir=args.out_dir, delete=False
    ) as coords_file:
        coords_path = Path(coords_file.name)
        writer = csv.writer(coords_file, delimiter=" ")

        # Process a small block of genes, then release it before the next block.
        for chunk in pd.read_csv(args.counts, chunksize=args.chunksize):
            # Save the gene names in exactly the same order as matrix rows.
            chunk_genes = chunk["Unnamed: 0"].astype(str).tolist()
            genes.extend(chunk_genes)

            # Select just count columns and turn this block into a numeric array.
            counts = chunk[count_columns].to_numpy(dtype=np.int64, copy=False)
            for row_offset in range(counts.shape[0]):
                # Sparse matrices record nonzero values only; zeros are implied.
                row_values = counts[row_offset]
                nonzero_cols = np.flatnonzero(row_values)
                # Matrix Market uses 1-based row/column positions, unlike Python.
                row_index = n_genes + row_offset + 1
                for col_index in nonzero_cols:
                    writer.writerow((row_index, col_index + 1, int(row_values[col_index])))
                nnz += int(nonzero_cols.size)

            # Update progress after this full block of gene rows was handled.
            n_genes += counts.shape[0]
            print(f"processed {n_genes:,} genes; nonzero entries so far: {nnz:,}", flush=True)

    # Make duplicate gene names unique before R uses them as row names.
    genes = make_unique(genes)
    write_lines(genes_path, genes)

    # Write the Matrix Market header followed by every saved coordinate row.
    # The header needs the final dimensions and nnz, which are known only now.
    with matrix_path.open("w", encoding="utf-8") as matrix_file:
        matrix_file.write("%%MatrixMarket matrix coordinate integer general\n")
        matrix_file.write("% GSE200164_counts.csv converted by prepare_gse200164_sparse.py\n")
        matrix_file.write(f"{n_genes} {n_cells} {nnz}\n")
        # Copy the temporary coordinate records into the final matrix file.
        with coords_path.open("r", encoding="utf-8") as coords_file:
            for line in coords_file:
                matrix_file.write(line)

    # Remove the temporary file now that the final matrix has been written.
    coords_path.unlink()
    # Create the completion marker last, so interrupted conversions are not reused.
    done_path.write_text("ok\n")
    print(f"Wrote sparse cache to {args.out_dir}")


if __name__ == "__main__":
    # This allows the functions above to be imported without starting conversion.
    main()
