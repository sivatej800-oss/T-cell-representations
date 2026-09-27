#!/usr/bin/env python3
"""Stage 11: silhouette with dimensionality held fixed.

Stage 07 found that every partition scores higher in scGPT coordinates than in PCA
coordinates, including the PCA partition itself. That has two readings — the
embedding organizes these cells more cleanly, or silhouette simply behaves
differently in 512 dimensions than in 30 — and Stage 07 cannot separate them,
because its two spaces differ in dimension as well as in kind.

This stage removes the confound by scoring PCA-512 against scGPT-512, where the
only remaining difference is how the coordinates were produced. The PCA-512 labels
come from Stage 10, which clustered that arm fresh rather than reusing the
30-dimension labels.

The PCA-30 labels are carried along as a third partition, scored in both
512-dimensional spaces. They are not part of the matched comparison; they show what
changes when the scoring space grows and the labels stay put.

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

REQUIRED_R_PACKAGES = ("Seurat", "cluster", "ggplot2")

PCA30_ARM = "pca_harmony_30dims"
PCA512_ARM = "pca512_harmony"
SCGPT_ARM = "scgpt_harmony"

EXPECTED_CELLS = 70391
EXPECTED_DIMS = 512
DEFAULT_SEED = 42
DEFAULT_SAMPLE_SIZE = 25000


@dataclass(frozen=True)
class PipelinePaths:
    """All Stage 11 input and output paths below one workspace directory."""

    workspace: Path

    @classmethod
    def from_workspace(cls, workspace: str | Path) -> "PipelinePaths":
        return cls(Path(workspace).expanduser().resolve())

    @property
    def benchmark_dir(self) -> Path:
        return self.workspace / "results" / "06_harmony_benchmark"

    @property
    def matched_dir(self) -> Path:
        return self.workspace / "results" / "11_silhouette_matched"

    def create_directories(self) -> None:
        self.matched_dir.mkdir(parents=True, exist_ok=True)


def default_workspace() -> Path:
    """Use the repository's Age_trajectory directory unless --workspace is given."""

    return Path(__file__).resolve().parents[1] / "Age_trajectory"


def bundled_r_script() -> Path:
    script = Path(__file__).resolve().with_suffix(".R")
    if not script.exists():
        raise FileNotFoundError(f"Companion R script is missing: {script}")
    return script


def require_files(paths: Iterable[Path]) -> None:
    """Stop before doing any work if an input this stage depends on is absent."""

    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required inputs:\n  " + "\n  ".join(missing))


# --- Step 1: confirm both stages left what this comparison needs -------------
# Three files describing the same 70,391 cells, written by two different runs. The R
# script checks their cell order matches; this checks they exist and hold the arms
# they are supposed to before an hour of scoring starts.
def verify_inputs(paths: PipelinePaths) -> pd.DataFrame:
    stage06_assignments = paths.benchmark_dir / "benchmark_cluster_assignments.csv.gz"
    stage10_assignments = paths.benchmark_dir / "benchmark_pca512_assignments.csv.gz"
    require_files(
        [
            paths.benchmark_dir / "benchmark_representations.rds",
            paths.benchmark_dir / "benchmark_pca512_arm.rds",
            stage06_assignments,
            stage10_assignments,
        ]
    )

    stage06 = pd.read_csv(stage06_assignments)
    stage10 = pd.read_csv(stage10_assignments)

    for frame, name, arms in (
        (stage06, "Stage 06", (PCA30_ARM, SCGPT_ARM)),
        (stage10, "Stage 10", (PCA512_ARM,)),
    ):
        if len(frame) != EXPECTED_CELLS:
            raise ValueError(
                f"{name} holds {len(frame)} cells, expected {EXPECTED_CELLS}"
            )
        for arm in arms:
            if arm not in frame.columns:
                raise ValueError(f"{name} output does not contain the arm {arm}")

    if not stage06["cell_barcode"].equals(stage10["cell_barcode"]):
        raise ValueError(
            "Stage 06 and Stage 10 do not describe the same cells in the same order"
        )

    merged = stage06[["cell_barcode", PCA30_ARM, SCGPT_ARM]].copy()
    merged[PCA512_ARM] = stage10[PCA512_ARM]

    print("Partitions to be scored:")
    for arm in (PCA512_ARM, SCGPT_ARM, PCA30_ARM):
        sizes = merged[arm].value_counts()
        print(
            f"  {arm}: {len(sizes)} clusters, "
            f"smallest {sizes.min()}, largest {sizes.max()}"
        )
    return merged


# --- Step 2: report what each draw can actually resolve ----------------------
# A cluster drawn down to one or two cells has no meaningful within-cluster
# distance. Printing this before the run makes an unusable sample size visible up
# front rather than as a strange score afterwards.
def report_sampling_resolution(
    assignments: pd.DataFrame, sample_size: int
) -> dict[str, dict[str, int]]:
    resolution = {}
    print(f"\nAt {sample_size} cells per draw:")
    for arm in (PCA512_ARM, SCGPT_ARM, PCA30_ARM):
        sizes = assignments[arm].value_counts()
        smallest = int(sizes.min())
        drawn = max(2, round(sample_size * smallest / len(assignments)))
        resolution[arm] = {"smallest_cluster": smallest, "drawn_per_repeat": drawn}
        print(f"  {arm}: smallest cluster {smallest} -> about {drawn} cells drawn")
    return resolution


def audit_r_environment(rscript: str | None = None) -> dict[str, str]:
    executable = rscript or shutil.which("Rscript")
    if executable is None:
        raise RuntimeError("Rscript was not found on PATH")
    names = ", ".join(f'"{package}"' for package in REQUIRED_R_PACKAGES)
    code = (
        f"packages <- c({names}); "
        "cat(paste(packages, sapply(packages, function(p) "
        'as.character(packageVersion(p))), collapse="\\n"))'
    )
    result = subprocess.run(
        [executable, "--vanilla", "-e", code],
        check=True,
        capture_output=True,
        text=True,
    )
    versions = {}
    for line in result.stdout.strip().splitlines():
        name, _, version = line.partition(" ")
        versions[name] = version
    return versions


def run_matched_silhouette(
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
        str(paths.matched_dir),
        str(seed),
        str(sample_size),
    ]
    print("\nStarting R:\n  " + " ".join(command) + "\n")
    subprocess.run(command, check=True)

    summary = paths.matched_dir / "silhouette_matched_summary.csv"
    if not summary.exists():
        raise RuntimeError("R finished without writing the matched summary")
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

    assignments = verify_inputs(paths)
    resolution = report_sampling_resolution(assignments, args.sample_size)

    if args.verify_only:
        print("\nVerification complete; stopping before the analysis as requested.")
        return

    versions = audit_r_environment(rscript=args.rscript)
    summary = run_matched_silhouette(
        paths, seed=args.seed, sample_size=args.sample_size, rscript=args.rscript
    )

    audit = {
        "stage": "11_silhouette_dimension_matched",
        "seed": args.seed,
        "sample_size": args.sample_size,
        "matched_dims": EXPECTED_DIMS,
        "arms_compared": [PCA512_ARM, SCGPT_ARM],
        "arms_carried": [PCA30_ARM],
        "r_packages": versions,
        "sampling_resolution": resolution,
    }
    audit_path = paths.matched_dir / "silhouette_matched_environment_audit.json"
    audit_path.write_text(json.dumps(audit, indent=2), encoding="utf-8")

    print(f"\nStage 11 complete.\n  Summary: {summary}\n  Audit: {audit_path}")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:  # noqa: BLE001
        print(f"\nStage 11 failed: {error}", file=sys.stderr)
        raise SystemExit(1) from error
