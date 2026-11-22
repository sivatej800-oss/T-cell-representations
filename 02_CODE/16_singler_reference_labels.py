#!/usr/bin/env python3
"""Stage 16: an independent per-cell cell-type label from SingleR and Monaco.

Stages 13-15 named clusters by reading marker genes. This stage labels every cell on its
own, by correlating its expression against purified immune-cell profiles (celldex's
MonacoImmuneData), with no clusters involved. Because the label for a cell does not
depend on the PCA partition or the scGPT partition, it can later grade both without
belonging to either.

The RNA is the same unchanged seurat_input/ that Stage 13 reads, log-normalized the same
way. Cells SingleR cannot place confidently stay "unassigned" in every label column.

This file checks inputs and hands them to R. It performs no analysis.
"""

from __future__ import annotations

import argparse
import json
import resource
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


REQUIRED_R_PACKAGES = ("Seurat", "Matrix", "jsonlite", "SingleR", "celldex", "BiocParallel")

EXPECTED_CELLS = 70391
EXPECTED_GENES = 19689


@dataclass(frozen=True)
class PipelinePaths:
    """All Stage 16 input and output paths below one workspace directory."""

    workspace: Path

    @classmethod
    def from_workspace(cls, workspace: str | Path) -> "PipelinePaths":
        return cls(Path(workspace).expanduser().resolve())

    @property
    def seurat_input_dir(self) -> Path:
        return self.workspace / "results" / "03_seurat_reproduction" / "seurat_input"

    @property
    def stage13_audit(self) -> Path:
        return self.workspace / "results" / "13_marker_discovery" / "marker_discovery_audit.json"

    @property
    def output_dir(self) -> Path:
        return self.workspace / "results" / "16_singler_labels"

    def create_directories(self) -> None:
        self.output_dir.mkdir(parents=True, exist_ok=True)


def default_workspace() -> Path:
    """Use the repository's Age_trajectory directory unless --workspace is given."""

    return Path(__file__).resolve().parents[1] / "Age_trajectory"


def require_files(paths: Iterable[Path]) -> None:
    """Stop before doing any work if an input this stage depends on is absent."""

    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required inputs:\n  " + "\n  ".join(missing))


# --- Step 1: confirm the expression hand-off is the expected size -------------------
# Labels attach to cells by barcode, so the barcode and gene lists must be the ones the
# whole pipeline shares. A different matrix would give plausible labels for other cells.
def verify_expression(paths: PipelinePaths) -> None:
    """Check the RNA hand-off holds the expected genes and unique cells."""

    barcodes_file = paths.seurat_input_dir / "barcodes.tsv"
    genes_file = paths.seurat_input_dir / "genes.tsv"
    require_files([barcodes_file, genes_file, paths.seurat_input_dir / "counts.mtx"])

    barcodes = barcodes_file.read_text(encoding="utf-8").split()
    genes = genes_file.read_text(encoding="utf-8").split()
    if len(barcodes) != EXPECTED_CELLS or len(set(barcodes)) != EXPECTED_CELLS:
        raise ValueError(f"Expected {EXPECTED_CELLS} unique cells, found {len(barcodes)}")
    if len(genes) != EXPECTED_GENES or len(set(genes)) != EXPECTED_GENES:
        raise ValueError(f"Expected {EXPECTED_GENES} unique genes, found {len(genes)}")
    print(f"Expression: {len(genes)} genes x {len(barcodes)} cells.")


def audit_r_environment(rscript: str | None = None) -> dict[str, str]:
    """Confirm the R packages exist and record what versions ran."""

    executable = rscript or shutil.which("Rscript")
    if executable is None:
        raise RuntimeError("Rscript was not found on PATH")

    program = (
        "versions <- vapply(c("
        + ", ".join(f'"{package}"' for package in REQUIRED_R_PACKAGES)
        + "), function(p) tryCatch(as.character(packageVersion(p)), "
        'error = function(e) "MISSING"), character(1));'
        "cat(paste(names(versions), versions, sep = '=', collapse = '\\n'))"
    )
    result = subprocess.run(
        [executable, "-e", program], check=True, capture_output=True, text=True
    )
    versions = dict(
        line.split("=", 1) for line in result.stdout.strip().splitlines() if "=" in line
    )
    missing = [name for name, version in versions.items() if version == "MISSING"]
    if missing:
        raise RuntimeError("Missing required R packages: " + ", ".join(missing))
    print("R packages: " + ", ".join(f"{k} {v}" for k, v in versions.items()))
    return versions


def bundled_r_script() -> Path:
    """The R script that sits beside this one and does the actual work."""

    script = Path(__file__).resolve().with_suffix(".R")
    if not script.exists():
        raise FileNotFoundError(f"Companion R script not found at {script}")
    return script


# --- Step 2: run R and record peak memory -------------------------------------------
# The R script reports its own peak; the largest process R started (a SingleR worker or
# R itself) is added here from the operating system, since forked workers are not
# visible from inside the parent.
def run_singler(paths: PipelinePaths, rscript: str | None = None) -> dict[str, object]:
    """Start the R stage, confirm its outputs, and add the OS-measured peak memory."""

    executable = rscript or shutil.which("Rscript")
    if executable is None:
        raise RuntimeError("Rscript was not found on PATH")

    command = [
        executable,
        "--vanilla",
        str(bundled_r_script()),
        str(paths.seurat_input_dir),
        str(paths.output_dir),
        str(paths.stage13_audit),
    ]
    print("\nStarting R:\n  " + " ".join(command) + "\n")
    subprocess.run(command, check=True)

    expected = [
        "singler_cell_labels.csv.gz",
        "label_counts_broad.csv",
        "label_counts_fine.csv",
        "mapping_table.csv",
        "singler_audit.json",
    ]
    require_files(paths.output_dir / name for name in expected)

    audit_path = paths.output_dir / "singler_audit.json"
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    audit["peak_memory_gb_largest_process_os"] = round(
        resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024**2, 3
    )
    audit_path.write_text(json.dumps(audit, indent=2), encoding="utf-8")
    return audit


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--workspace",
        default=str(default_workspace()),
        help="Directory holding data/ and results/ (default: Age_trajectory)",
    )
    parser.add_argument("--rscript", help="Optional explicit path to Rscript")
    parser.add_argument(
        "--verify-only",
        action="store_true",
        help="Run the input checks, then stop before the analysis",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    paths = PipelinePaths.from_workspace(args.workspace)

    verify_expression(paths)
    if args.verify_only:
        print("\nVerification complete; stopping before the analysis as requested.")
        return

    paths.create_directories()
    audit_r_environment(rscript=args.rscript)
    audit = run_singler(paths, rscript=args.rscript)

    print(
        f"\nStage 16 complete: {audit['n_cells']} cells, {audit['n_unassigned']} unassigned "
        f"({audit['pct_unassigned']:.2f}%).\n  Outputs: {paths.output_dir}"
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as error:  # noqa: BLE001
        print(f"Stage 16 failed: {error}", file=sys.stderr)
        raise SystemExit(1) from error
