#!/usr/bin/env python3
"""Stage 17: does cluster structure follow biology or the patient?

For the two arms of the main comparison -- the authors' published pipeline (Stage 03) and
scGPT used out of the box (Stage 06's scgpt_raw) -- this stage measures how strongly each
clustering tracks the patient a cell came from and how strongly it tracks cell type, the
latter against SingleR labels (Stage 16) that belong to neither arm. It also compares each
arm's hand annotation (Stage 15) with SingleR cell by cell, and looks at marker genes for
the Monocyte/DC cells the two disagree about.

The earlier four-arm version of this stage (PCA and scGPT, each with and without Harmony)
is in git history; its outputs are kept in results/17_biology_vs_donor_harmony_record/.

This file checks inputs and hands them to R. It performs no biological analysis.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd


REQUIRED_R_PACKAGES = ("Seurat", "Matrix", "jsonlite", "SingleR", "celldex", "BiocParallel")

EXPECTED_CELLS = 70391
EXPECTED_GENES = 19689
TOTAL_CORES = 16

# The two arms, named as the R script names them.
ARMS = ("authors_pipeline", "scgpt_raw")

@dataclass(frozen=True)
class PipelinePaths:
    """All Stage 17 input and output paths below one workspace directory."""

    workspace: Path

    @classmethod
    def from_workspace(cls, workspace: str | Path) -> "PipelinePaths":
        return cls(Path(workspace).expanduser().resolve())

    @property
    def seurat_input_dir(self) -> Path:
        return self.workspace / "results" / "03_seurat_reproduction" / "seurat_input"

    @property
    def clusters(self) -> Path:
        return self.workspace / "results" / "06_harmony_benchmark" / "benchmark_cluster_assignments.csv.gz"

    @property
    def donor_leakage(self) -> Path:
        return self.workspace / "results" / "06_harmony_benchmark" / "benchmark_donor_leakage.csv"

    @property
    def singler(self) -> Path:
        return self.workspace / "results" / "16_singler_labels" / "singler_cell_labels.csv.gz"

    @property
    def annotation(self) -> Path:
        return self.workspace / "results" / "15_annotation" / "annotated_cluster_labels.csv.gz"

    @property
    def figure_dir(self) -> Path:
        return self.workspace / "results" / "figures"

    @property
    def output_dir(self) -> Path:
        return self.workspace / "results" / "17_biology_vs_donor"

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


# --- Step 1: every input must describe the same cells -------------------------------
# Labels are joined by barcode downstream, so every source must hold exactly the matrix's
# 70,391 unique barcodes. A source that silently held a different set would attach one
# cell's label to another cell's expression.
def verify_same_cells(paths: PipelinePaths) -> dict[str, pd.DataFrame]:
    """Load every label source and confirm they all cover the matrix's cells."""

    require_files(
        [
            paths.seurat_input_dir / "barcodes.tsv",
            paths.seurat_input_dir / "genes.tsv",
            paths.seurat_input_dir / "counts.mtx",
            paths.clusters,
            paths.donor_leakage,
            paths.singler,
            paths.annotation,
        ]
    )
    barcodes = (paths.seurat_input_dir / "barcodes.tsv").read_text(encoding="utf-8").split()
    genes = (paths.seurat_input_dir / "genes.tsv").read_text(encoding="utf-8").split()
    if len(barcodes) != EXPECTED_CELLS or len(set(barcodes)) != EXPECTED_CELLS:
        raise ValueError(f"Expected {EXPECTED_CELLS} unique cells, found {len(barcodes)}")
    if len(genes) != EXPECTED_GENES or len(set(genes)) != EXPECTED_GENES:
        raise ValueError(f"Expected {EXPECTED_GENES} unique genes, found {len(genes)}")

    sources = {
        "stage06_clusters": pd.read_csv(paths.clusters, dtype=str),
        "singler": pd.read_csv(paths.singler, dtype=str),
        "annotation": pd.read_csv(paths.annotation, dtype=str),
    }
    for name, frame in sources.items():
        if len(frame) != EXPECTED_CELLS or frame["cell_barcode"].duplicated().any():
            raise ValueError(f"{name}: expected {EXPECTED_CELLS} unique barcodes")
        if set(frame["cell_barcode"]) != set(barcodes):
            raise ValueError(f"{name}: barcodes differ from seurat_input/barcodes.tsv")
    print(f"All inputs cover the same {EXPECTED_CELLS} unique barcodes.")
    return {name: frame.set_index("cell_barcode") for name, frame in sources.items()}


# --- Step 2: choose how many cores to take ---------------------------------------------
# Other work on this machine comes in bursts, so the busier of the 1- and 5-minute load
# averages is used: the 1-minute figure alone can catch a lull. Free cores are the
# machine's 16 minus the whole cores busy, with one left spare.
def choose_workers() -> tuple[int, tuple[float, float]]:
    """Return (workers, (load_1min, load_5min)) using the load at launch."""

    load_1, load_5, _ = os.getloadavg()
    free = TOTAL_CORES - math.floor(max(load_1, load_5))
    return max(1, free - 1), (load_1, load_5)


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


def run_r(paths: PipelinePaths, workers: int, rscript: str | None = None) -> dict[str, object]:
    """Start the R stage and confirm it produced every table."""

    executable = rscript or shutil.which("Rscript")
    if executable is None:
        raise RuntimeError("Rscript was not found on PATH")
    command = [
        executable,
        "--vanilla",
        str(bundled_r_script()),
        str(paths.seurat_input_dir),
        str(paths.clusters),
        str(paths.donor_leakage),
        str(paths.singler),
        str(paths.annotation),
        str(paths.output_dir),
        str(workers),
        str(paths.figure_dir),
    ]
    print("\nStarting R:\n  " + " ".join(command) + "\n")
    subprocess.run(command, check=True)

    expected = [
        "nmi_by_arm.csv",
        "annotation_vs_singler_agreement.csv",
        *[
            f"{table}_{arm}.csv"
            for arm in ARMS
            for table in ("agreement_per_type", "confusion", "confusion_merged_myeloid")
        ],
        "myeloid_merge_table.csv",
        "monocyte_dc_marker_groups.csv",
        "confidence_by_broad_label.csv",
        "broad_level_confidence.csv.gz",
        "annotation_table_used.csv",
        "stage17_audit_r.json",
    ]
    require_files(paths.output_dir / name for name in expected)
    require_files(
        paths.figure_dir / name
        for name in (
            "17_patient_vs_biology_scatter.png",
            *[f"17_agreement_heatmap_{arm}.png" for arm in ARMS],
            *[f"17_umap_singler_vs_annotation_{arm}.png" for arm in ARMS],
        )
    )
    return json.loads((paths.output_dir / "stage17_audit_r.json").read_text(encoding="utf-8"))


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

    sources = verify_same_cells(paths)
    paths.create_directories()
    if args.verify_only:
        print("\nVerification complete; stopping before the analysis as requested.")
        return

    versions = audit_r_environment(rscript=args.rscript)
    workers, load = choose_workers()
    print(f"Load average {load[0]:.2f} (1 min), {load[1]:.2f} (5 min) -> {workers} SingleR workers.")
    r_audit = run_r(paths, workers, rscript=args.rscript)

    audit = {
        "stage": "17_cluster_structure_vs_singler",
        "inputs": {
            "workspace": str(paths.workspace),
            "clusters": str(paths.clusters),
            "donor_leakage": str(paths.donor_leakage),
            "singler_labels": str(paths.singler),
            "seurat_input": str(paths.seurat_input_dir),
            "stage15_annotation": str(paths.annotation),
        },
        "same_cells_in_every_input": EXPECTED_CELLS,
        "workers": workers,
        "load_average_at_launch": {"1min": load[0], "5min": load[1]},
        "r_packages": versions,
        "r_stage": r_audit,
    }
    audit_path = paths.output_dir / "stage17_audit.json"
    audit_path.write_text(json.dumps(audit, indent=2), encoding="utf-8")
    (paths.output_dir / "stage17_audit_r.json").unlink()

    print(f"\nStage 17 complete.\n  Outputs: {paths.output_dir}\n  Audit: {audit_path}")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:  # noqa: BLE001
        print(f"Stage 17 failed: {error}", file=sys.stderr)
        raise SystemExit(1) from error
