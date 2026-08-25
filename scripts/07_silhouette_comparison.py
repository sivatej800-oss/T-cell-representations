#!/usr/bin/env python3
"""Stage 07: measure how well-defined each arm's Louvain clusters are.

Stage 06 showed that the PCA and scGPT descriptions find related but not identical
structure, and left open which one draws better-separated groups. This stage answers
that with silhouette widths, computed on the clusters Stage 06 already produced.
Nothing is re-clustered and no k is chosen: the Louvain partitions are the input.

The comparison is not a matter of scoring each arm in its own space and reading off
the larger number. Silhouette is built from distances, and 30 principal components
and 512 learned dimensions are not comparable — distances concentrate as dimension
grows, which lowers silhouette for reasons unrelated to how good the grouping is. So
each partition is scored in both spaces. Within one space the geometry is fixed and
the two partitions can be compared directly; a partition that scores higher in both
spaces is better defined in a way that does not depend on where it is measured.

This file checks inputs and hands them to R. It performs no analysis.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd


# --- Required R software -----------------------------------------------------
# cluster provides silhouette; Seurat is needed only to read the saved object.
REQUIRED_R_PACKAGES = ("Seurat", "cluster", "ggplot2")

# The two arms compared, named as Stage 06 named them.
PCA_ARM = "pca_harmony_30dims"
SCGPT_ARM = "scgpt_harmony"

EXPECTED_CELLS = 70391
DEFAULT_SEED = 42

# Silhouette needs every pairwise distance within the scored set, so the full 70,391
# cells are out of reach and a subsample is scored instead. Stage 07 draws several
# times; the spread across draws is reported so that a difference between the arms
# can be checked against it rather than read as exact.
DEFAULT_SAMPLE_SIZE = 12000


@dataclass(frozen=True)
class PipelinePaths:
    """All Stage 07 input and output paths below one workspace directory."""

    workspace: Path

    @classmethod
    def from_workspace(cls, workspace: str | Path) -> "PipelinePaths":
        return cls(Path(workspace).expanduser().resolve())

    @property
    def benchmark_dir(self) -> Path:
        return self.workspace / "results" / "06_harmony_benchmark"

    @property
    def silhouette_dir(self) -> Path:
        return self.workspace / "results" / "07_silhouette"

    def create_directories(self) -> None:
        self.silhouette_dir.mkdir(parents=True, exist_ok=True)


def default_workspace() -> Path:
    """Use the repository's Age_trajectory directory unless --workspace is given."""

    return Path(__file__).resolve().parents[1] / "Age_trajectory"


def require_files(paths: Iterable[Path]) -> None:
    """Stop before doing any work if an input this stage depends on is absent."""

    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required inputs:\n  " + "\n  ".join(missing))


# --- Step 1: confirm Stage 06's output is present and complete ---------------
# Both files are needed: the saved object holds the coordinates distances are
# measured in, and the assignments table holds the clusters being scored. They must
# describe the same cells in the same order, which the R script re-checks before
# using them.
def verify_benchmark_output(paths: PipelinePaths) -> pd.DataFrame:
    """Check that Stage 06 finished and produced both arms this stage compares."""

    object_file = paths.benchmark_dir / "benchmark_representations.rds"
    assignments_file = paths.benchmark_dir / "benchmark_cluster_assignments.csv.gz"
    require_files([object_file, assignments_file])

    assignments = pd.read_csv(assignments_file)
    if len(assignments) != EXPECTED_CELLS:
        raise ValueError(
            f"Expected {EXPECTED_CELLS} cells in Stage 06's assignments, "
            f"found {len(assignments)}"
        )
    for arm in (PCA_ARM, SCGPT_ARM):
        if arm not in assignments.columns:
            raise ValueError(f"Stage 06 output does not contain the arm {arm}")

    for arm in (PCA_ARM, SCGPT_ARM):
        sizes = assignments[arm].value_counts()
        print(
            f"{arm}: {len(sizes)} clusters, "
            f"{sizes.min()}-{sizes.max()} cells per cluster"
        )
    return assignments


# --- Step 2: report what the sample can and cannot resolve -------------------
# The smallest cluster decides whether the comparison can say anything about rare
# populations. Drawing proportionally, a cluster holding 110 of 70,391 cells
# contributes about 19 cells to a 12,000-cell draw — enough for a mean, but with
# enough spread that its per-cluster number should be read as approximate.
def report_sampling_resolution(assignments: pd.DataFrame, sample_size: int) -> dict[str, object]:
    """State how many cells the smallest clusters contribute to each draw."""

    report: dict[str, object] = {}
    print(f"\nAt {sample_size} cells per draw, drawn proportionally:")
    for arm in (PCA_ARM, SCGPT_ARM):
        sizes = assignments[arm].value_counts()
        smallest = int(sizes.min())
        drawn = max(2, round(sample_size * smallest / len(assignments)))
        print(
            f"  {arm}: smallest cluster has {smallest} cells, "
            f"about {drawn} per draw"
        )
        if drawn < 10:
            print(
                "    (fewer than 10 drawn; this cluster's own score is indicative "
                "only, though it still contributes to the overall mean)"
            )
        report[arm] = {"smallest_cluster": smallest, "drawn_per_repeat": int(drawn)}
    return report


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


def run_silhouette(
    paths: PipelinePaths, seed: int, sample_size: int, rscript: str | None = None
) -> Path:
    """Start the R stage and confirm it produced the result it claims to."""

    executable = rscript or shutil.which("Rscript")
    if executable is None:
        raise RuntimeError("Rscript was not found on PATH")

    command = [
        executable,
        "--vanilla",
        str(bundled_r_script()),
        str(paths.benchmark_dir),
        str(paths.silhouette_dir),
        str(seed),
        str(sample_size),
    ]
    print("\nStarting R:\n  " + " ".join(command) + "\n")
    subprocess.run(command, check=True)

    summary = paths.silhouette_dir / "silhouette_summary.csv"
    if not summary.exists():
        raise RuntimeError("R finished without writing the silhouette summary")
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--workspace",
        default=str(default_workspace()),
        help="Directory holding data/ and results/ (default: Age_trajectory)",
    )
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED, help="Random seed")
    parser.add_argument(
        "--sample-size",
        type=int,
        default=DEFAULT_SAMPLE_SIZE,
        help=f"Cells per silhouette draw (default: {DEFAULT_SAMPLE_SIZE})",
    )
    parser.add_argument("--rscript", help="Optional explicit path to Rscript")
    parser.add_argument(
        "--verify-only",
        action="store_true",
        help="Run the input checks and the sampling report, then stop",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    paths = PipelinePaths.from_workspace(args.workspace)
    paths.create_directories()

    assignments = verify_benchmark_output(paths)
    resolution = report_sampling_resolution(assignments, args.sample_size)

    if args.verify_only:
        print("\nVerification complete; stopping before the analysis as requested.")
        return

    versions = audit_r_environment(rscript=args.rscript)
    summary = run_silhouette(
        paths, seed=args.seed, sample_size=args.sample_size, rscript=args.rscript
    )

    audit = {
        "stage": "07_silhouette_comparison",
        "seed": args.seed,
        "sample_size": args.sample_size,
        "arms_compared": [PCA_ARM, SCGPT_ARM],
        "r_packages": versions,
        "sampling_resolution": resolution,
    }
    audit_path = paths.silhouette_dir / "silhouette_environment_audit.json"
    audit_path.write_text(json.dumps(audit, indent=2), encoding="utf-8")

    print(f"\nStage 07 complete.\n  Summary: {summary}\n  Audit: {audit_path}")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:  # noqa: BLE001
        print(f"Stage 07 failed: {error}", file=sys.stderr)
        raise SystemExit(1) from error
