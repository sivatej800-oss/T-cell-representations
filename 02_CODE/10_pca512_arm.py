"""Build a PCA arm at the embedding's dimensionality, to control the comparison.

Stage 07 scored a 13-cluster PCA partition in 30 dimensions against an 11-cluster
scGPT partition in 512, and found every partition scores higher in scGPT
coordinates than in PCA coordinates. Two explanations fit that: the embedding
organizes the cells more cleanly, or silhouette simply behaves differently in 512
dimensions. Distances concentrate as dimension grows, and silhouette is built from
distances, so dimensionality moves the score by itself.

This stage builds the control that separates them — PCA carried to 512 components,
Harmony-corrected, and clustered fresh at Stage 06's settings. The clustering is
new rather than the 30-dimension labels re-scored, because the comparison needs a
PCA-512 partition, not the PCA-30 partition viewed from further away.

The heavy work is in 10_pca512_arm.R; this file checks the inputs, records the
environment, and reports what the arm produced.
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

DEFAULT_SEED = 42
DEFAULT_DIMS = 512
EXPECTED_CELLS = 70391


@dataclass(frozen=True)
class PipelinePaths:
    workspace: Path

    @classmethod
    def from_workspace(cls, workspace: str | Path) -> "PipelinePaths":
        return cls(Path(workspace).expanduser().resolve())

    @property
    def benchmark_dir(self) -> Path:
        return self.workspace / "results" / "06_harmony_benchmark"

    @property
    def figures_dir(self) -> Path:
        return self.workspace / "results" / "figures"

    def create_directories(self) -> None:
        self.figures_dir.mkdir(parents=True, exist_ok=True)


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


# --- Step 1: confirm Stage 06 left the checkpoint this arm rebuilds from -----
# The uncorrected checkpoint carries the merged, scaled expression with no donor
# correction applied. Rebuilding from it is what keeps this arm corrected exactly
# once, by the same Harmony the scGPT arm gets.
def verify_inputs(paths: PipelinePaths) -> None:
    checkpoint = paths.benchmark_dir / "benchmark_pca_uncorrected.rds"
    require_files([checkpoint])
    size_gb = checkpoint.stat().st_size / 1024**3
    print(f"Stage 06 uncorrected checkpoint: {checkpoint} ({size_gb:.1f} GB)")


def audit_r_environment(rscript: str | None = None) -> dict[str, str]:
    """Record the package versions this run used, as the other stages do."""

    executable = rscript or shutil.which("Rscript")
    if executable is None:
        raise RuntimeError("Rscript was not found on PATH")
    code = (
        'cat(paste(c("Seurat", "harmony", "ggplot2"), '
        'sapply(c("Seurat", "harmony", "ggplot2"), '
        'function(p) as.character(packageVersion(p))), collapse="\\n"))'
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


def run_arm(
    paths: PipelinePaths, seed: int, n_dims: int, rscript: str | None = None
) -> Path:
    """Start the R stage and confirm it produced the labels it claims to."""

    executable = rscript or shutil.which("Rscript")
    if executable is None:
        raise RuntimeError("Rscript was not found on PATH")

    command = [
        executable,
        "--vanilla",
        str(bundled_r_script()),
        str(paths.benchmark_dir),
        str(paths.figures_dir),
        str(seed),
        str(n_dims),
        str(EXPECTED_CELLS),
    ]
    print("\nStarting R:\n  " + " ".join(command) + "\n")
    subprocess.run(command, check=True)

    assignments = paths.benchmark_dir / "benchmark_pca512_assignments.csv.gz"
    if not assignments.exists():
        raise RuntimeError("R finished without writing the PCA-512 assignments")
    return assignments


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--workspace",
        default=str(default_workspace()),
        help="Directory holding data/ and results/ (default: Age_trajectory)",
    )
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED, help="Random seed")
    parser.add_argument(
        "--dims",
        type=int,
        default=DEFAULT_DIMS,
        help=f"Principal components to compute (default: {DEFAULT_DIMS})",
    )
    parser.add_argument("--rscript", help="Optional explicit path to Rscript")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    paths = PipelinePaths.from_workspace(args.workspace)
    paths.create_directories()

    verify_inputs(paths)
    versions = audit_r_environment(rscript=args.rscript)
    assignments = run_arm(
        paths, seed=args.seed, n_dims=args.dims, rscript=args.rscript
    )

    summary_path = paths.benchmark_dir / "benchmark_pca512_summary.csv"
    summary = summary_path.read_text(encoding="utf-8") if summary_path.exists() else ""

    audit = {
        "stage": "10_pca512_arm",
        "seed": args.seed,
        "n_dims": args.dims,
        "resolution": 0.3,
        "r_packages": versions,
    }
    audit_path = paths.benchmark_dir / "pca512_environment_audit.json"
    audit_path.write_text(json.dumps(audit, indent=2), encoding="utf-8")

    print(f"\nStage 10 complete.\n  Labels: {assignments}\n  Audit: {audit_path}")
    if summary:
        print("\n" + summary)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:  # noqa: BLE001
        print(f"\nStage 10 failed: {error}", file=sys.stderr)
        raise SystemExit(1) from error
