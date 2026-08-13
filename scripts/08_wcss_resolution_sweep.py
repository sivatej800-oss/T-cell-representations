#!/usr/bin/env python3
"""Stage 08: an elbow-style compactness curve for each representation.

The elbow method plots how much of a dataset's spread a grouping accounts for
against the number of groups, and looks for the point where adding groups stops
helping. It is normally run with k-means, which takes the number of groups as an
argument. Louvain does not: it takes a resolution, and the number of clusters is an
outcome. So this stage sweeps resolution from 0.1 to 1.0, records whatever cluster
count each setting produced, and plots against that count.

Compactness is the within-cluster sum of squared distances to each cluster centroid,
divided by the total sum of squares about the representation's overall centroid. The
complement of that ratio is the fraction of total dispersion the grouping explains.

The two curves are computed in different spaces — 30 principal components against
512 learned dimensions — so their heights are not directly comparable. What is
comparable is how steeply each falls as clusters are added, which is a property of
how structure is organized within a space rather than of the space's size.

Nothing is recomputed: the coordinates are Stage 06's Harmony-corrected output and
the clustering is the same Seurat procedure at the same seed. Stage 06's files are
read only. This file checks inputs and hands them to R; it performs no analysis.
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


REQUIRED_R_PACKAGES = ("Seurat", "ggplot2")

# The two arms swept, and the Stage 06 columns whose cluster counts the sweep must
# reproduce at resolution 0.3.
REPRESENTATIONS = {
    "pca_harmony_30dims": "pca_harmony_30dims",
    "scgpt_harmony": "scgpt_harmony",
}

EXPECTED_CELLS = 70391
DEFAULT_SEED = 42
STAGE06_RESOLUTION = 0.3


@dataclass(frozen=True)
class PipelinePaths:
    """All Stage 08 input and output paths below one workspace directory."""

    workspace: Path

    @classmethod
    def from_workspace(cls, workspace: str | Path) -> "PipelinePaths":
        return cls(Path(workspace).expanduser().resolve())

    @property
    def benchmark_dir(self) -> Path:
        return self.workspace / "results" / "benchmark"

    @property
    def sweep_dir(self) -> Path:
        return self.workspace / "results" / "wcss_sweep"

    def create_directories(self) -> None:
        self.sweep_dir.mkdir(parents=True, exist_ok=True)


def default_workspace() -> Path:
    """Use the repository's Age_trajectory directory unless --workspace is given."""

    return Path(__file__).resolve().parents[1] / "Age_trajectory"


def require_files(paths: Iterable[Path]) -> None:
    """Stop before doing any work if an input this stage depends on is absent."""

    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required inputs:\n  " + "\n  ".join(missing))


# --- Step 1: confirm Stage 06's output is present ----------------------------
# The saved object holds the corrected coordinates the sweep clusters on, and the
# assignments table holds the cluster counts the sweep must reproduce at resolution
# 0.3. Reproducing them is what proves this sweep clusters the way Stage 06 did.
def verify_benchmark_output(paths: PipelinePaths) -> dict[str, int]:
    """Check Stage 06 finished, and record the counts the sweep has to match."""

    object_file = paths.benchmark_dir / "benchmark_representations.rds"
    assignments_file = paths.benchmark_dir / "benchmark_cluster_assignments.csv.gz"
    require_files([object_file, assignments_file])

    assignments = pd.read_csv(assignments_file)
    if len(assignments) != EXPECTED_CELLS:
        raise ValueError(
            f"Expected {EXPECTED_CELLS} cells in Stage 06's assignments, "
            f"found {len(assignments)}"
        )

    expected: dict[str, int] = {}
    for name, column in REPRESENTATIONS.items():
        if column not in assignments.columns:
            raise ValueError(f"Stage 06 output does not contain the arm {column}")
        expected[name] = int(assignments[column].nunique())
        print(f"{name}: Stage 06 recorded {expected[name]} clusters at resolution {STAGE06_RESOLUTION}")
    return expected


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


def run_sweep(paths: PipelinePaths, seed: int, rscript: str | None = None) -> Path:
    """Start the R stage and confirm it produced the result it claims to."""

    executable = rscript or shutil.which("Rscript")
    if executable is None:
        raise RuntimeError("Rscript was not found on PATH")

    command = [
        executable,
        "--vanilla",
        str(bundled_r_script()),
        str(paths.benchmark_dir),
        str(paths.sweep_dir),
        str(seed),
    ]
    print("\nStarting R:\n  " + " ".join(command) + "\n")
    subprocess.run(command, check=True)

    table = paths.sweep_dir / "wcss_resolution_sweep.csv"
    if not table.exists():
        raise RuntimeError("R finished without writing the sweep table")
    return table


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--workspace",
        default=str(default_workspace()),
        help="Directory holding data/ and results/ (default: Age_trajectory)",
    )
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED, help="Random seed")
    parser.add_argument("--rscript", help="Optional explicit path to Rscript")
    parser.add_argument(
        "--verify-only",
        action="store_true",
        help="Run the input checks, then stop",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    paths = PipelinePaths.from_workspace(args.workspace)
    paths.create_directories()

    expected = verify_benchmark_output(paths)

    if args.verify_only:
        print("\nVerification complete; stopping before the analysis as requested.")
        return

    versions = audit_r_environment(rscript=args.rscript)
    table = run_sweep(paths, seed=args.seed, rscript=args.rscript)

    # Report where Stage 06's setting falls on each curve, which is the point that
    # connects this sweep to every other result in the project.
    sweep = pd.read_csv(table)
    at_stage06 = sweep[sweep["resolution"] == STAGE06_RESOLUTION]
    print("\nAt Stage 06's resolution:")
    for _, row in at_stage06.iterrows():
        print(
            f"  {row['representation']}: {int(row['n_clusters'])} clusters, "
            f"normalized WCSS {row['normalized_wcss']:.4f}, "
            f"dispersion explained {1 - row['normalized_wcss']:.4f}"
        )

    audit = {
        "stage": "08_wcss_resolution_sweep",
        "seed": args.seed,
        "resolutions": [round(0.1 * i, 1) for i in range(1, 11)],
        "representations": list(REPRESENTATIONS),
        "stage06_cluster_counts": expected,
        "r_packages": versions,
    }
    audit_path = paths.sweep_dir / "wcss_sweep_audit.json"
    audit_path.write_text(json.dumps(audit, indent=2), encoding="utf-8")

    print(f"\nStage 08 complete.\n  Table: {table}\n  Audit: {audit_path}")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:  # noqa: BLE001
        print(f"Stage 08 failed: {error}", file=sys.stderr)
        raise SystemExit(1) from error
