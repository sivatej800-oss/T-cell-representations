#!/usr/bin/env python3
"""Stage 03: run or resume the R/Seurat Figure 1D analysis.

The actual analysis happens in R, because the published method uses Seurat. This
file's job is everything around that: work out whether there is already saved
progress worth resuming from, confirm R is installed and capable, start the R
script with the right settings, and check afterwards that it really produced what
it promised. It performs no analysis of its own.
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

import psutil


# --- Required R software -----------------------------------------------------
# Seurat performs the merge and clustering, Matrix carries the sparse data into it,
# glmGamPoi makes the per-sample normalization fast enough to be practical, and
# ggplot2 draws the final figure.
REQUIRED_R_PACKAGES = ("Seurat", "Matrix", "glmGamPoi", "ggplot2")

# How much of the R workflow to recompute, cheapest first. Each level discards its
# own checkpoint and every later one, so asking to redo clustering never triggers
# the multi-hour SCTransform step. The R script accepts these same names.
RESTART_LEVELS = ("none", "clustering", "pca", "integration", "all")


# --- Workspace paths ---------------------------------------------------------
# Both Python and R need to agree on where the data is, so paths are settled here
# once and handed to R as complete, unambiguous locations.
@dataclass(frozen=True)
class PipelinePaths:
    """All Stage 03 input and output paths below one workspace directory."""

    workspace: Path

    @classmethod
    def from_workspace(cls, workspace: str | Path) -> "PipelinePaths":
        return cls(Path(workspace).expanduser().resolve())

    # Where the analysis writes everything: its saved progress, the final object,
    # the cluster table, and the figure.
    @property
    def figure1d_dir(self) -> Path:
        return self.workspace / "results" / "figure1d_python"

    # Where Stage 02 left the data for R to pick up, if it needed to write it.
    @property
    def seurat_input_dir(self) -> Path:
        return self.figure1d_dir / "seurat_input"

    def create_directories(self) -> None:
        """Create this stage's folders, leaving saved progress in place."""

        for directory in (self.workspace, self.figure1d_dir, self.seurat_input_dir):
            directory.mkdir(parents=True, exist_ok=True)


def default_workspace() -> Path:
    """Use the repository's Age_trajectory directory unless --workspace is given."""

    return Path(__file__).resolve().parents[1] / "Age_trajectory"


# --- Small validation helper ------------------------------------------------
def require_files(paths: Iterable[Path]) -> None:
    """Stop with a readable list when an expected input file is missing."""

    missing = [path for path in paths if not path.is_file()]
    if missing:
        formatted = "\n".join(f"  - {path}" for path in missing)
        raise FileNotFoundError(f"Required files were not found:\n{formatted}")


# --- R and computer-resource preflight --------------------------------------
# The same check Stage 02 performs, repeated here because this stage can be run on
# its own or on a different machine, and because the analysis it is about to start
# is by far the most expensive part of the pipeline.
#
# Recording the package versions matters for interpreting the result: Seurat's
# clustering and UMAP change between releases, which is why this reproduction finds
# a different number of clusters than the paper reported.
def audit_r_environment(
    paths: PipelinePaths,
    *,
    rscript: str | None = None,
) -> dict[str, object]:
    """Confirm R and its required packages are available before any long work."""

    executable = rscript or shutil.which("Rscript")
    if executable is None:
        raise RuntimeError("Rscript was not found on PATH")
    # Ask R which of the required packages it has, and at what versions.
    package_vector = ", ".join(json.dumps(package) for package in REQUIRED_R_PACKAGES)
    code = f"""
packages <- c({package_vector})
for (pkg in packages) {{
  installed <- requireNamespace(pkg, quietly = TRUE)
  version <- if (installed) as.character(packageVersion(pkg)) else NA_character_
  cat(pkg, installed, version, sep = "\\t")
  cat("\\n")
}}
"""
    result = subprocess.run(
        [executable, "-e", code], check=True, capture_output=True, text=True
    )
    packages: dict[str, dict[str, object]] = {}
    for line in result.stdout.splitlines():
        name, installed, version = line.split("\t")
        packages[name] = {
            "installed": installed == "TRUE",
            "version": None if version == "NA" else version,
        }
    missing = [name for name in REQUIRED_R_PACKAGES if not packages[name]["installed"]]
    # Note the machine's memory and free disk too. Merging 59 samples is demanding on
    # both, so if a run later fails this records what was available at the time.
    memory = psutil.virtual_memory()
    disk = shutil.disk_usage(paths.figure1d_dir)
    audit = {
        "python_version": sys.version,
        "rscript_path": executable,
        "required_r_packages": packages,
        "missing_r_packages": missing,
        "total_ram_gb": memory.total / 1024**3,
        "available_ram_gb": memory.available / 1024**3,
        "disk_free_gb": disk.free / 1024**3,
        "ready": not missing,
    }
    # Save it beside the results, so the environment that produced a given figure
    # can be inspected later without re-running anything.
    output = paths.figure1d_dir / "figure1d_environment_audit.json"
    output.write_text(json.dumps(audit, indent=2, default=str) + "\n", encoding="utf-8")
    print(f"R environment audit: missing packages = {missing or 'none'}")
    return audit


# --- Locate the only companion file -----------------------------------------
def bundled_r_script() -> Path:
    """Find the R script that does the actual analysis, stored in this folder."""

    return Path(__file__).with_name("03_reproduce_figure1d.R")


# --- Stage 03 controller -----------------------------------------------------
# The whole stage in order. Everything here is preparation and checking; the
# analysis itself happens in the R script this function starts.
#
# Going in: either the data Stage 02 wrote out for R, or saved progress from an
# earlier attempt at this stage.
# Coming out: the finished analysis — every cell placed on a two-dimensional map
# and assigned to a cluster — plus a cluster table and the figure.
#
# The resume logic exists because this stage can take many hours. R decides which
# saved progress to use; this function only works out whether the original data
# files are needed at all, so it can complain early if they are missing.
def run_seurat_figure1d(
    paths: PipelinePaths,
    *,
    rscript: str | None = None,
    seed: int = 1234,
    restart_from: str = "none",
) -> Path:
    """Check everything is in place, run the R analysis, and confirm it succeeded."""

    paths.create_directories()

    # Step 1: the four files Stage 02 writes when R is starting from nothing. They
    # hold the measurements plus the gene names, cell names, and per-cell details.
    raw_inputs = [
        paths.seurat_input_dir / "counts.mtx",
        paths.seurat_input_dir / "genes.tsv",
        paths.seurat_input_dir / "barcodes.tsv",
        paths.seurat_input_dir / "metadata.csv",
    ]
    # Step 2: if any saved progress exists, R will resume from it and never look at
    # those four files, so their absence is fine. With no saved progress, they are
    # the only way in, and a missing one should be reported now.
    reusable_checkpoints = [
        paths.figure1d_dir / "integrated_clustered_umap.rds",
        paths.figure1d_dir / "GSE200164_after_PCA.rds",
        paths.figure1d_dir / "GSE200164_seurat_integrated.rds",
    ]
    if not any(path.exists() for path in reusable_checkpoints):
        require_files(raw_inputs)
    else:
        print("Existing final/integrated/PCA checkpoint found; raw Seurat export is not required.")
    # Step 3: find R and confirm it can do the job before committing to a long run.
    executable = rscript or shutil.which("Rscript")
    if executable is None:
        raise RuntimeError("Rscript was not found on PATH")
    audit = audit_r_environment(paths, rscript=executable)
    if audit["missing_r_packages"]:
        raise RuntimeError(
            "Cannot run Seurat; missing R packages: "
            + ", ".join(audit["missing_r_packages"])
        )

    # Step 4: locate the R script that performs the analysis.
    script = bundled_r_script()
    require_files([script])
    # Step 5: start it, telling it where to read, where to write, which random seed
    # to use so the result is repeatable, and how much work to redo. Each value is
    # passed separately rather than as one typed-out command line, so a folder name
    # containing a space cannot break the call.
    command = [
        executable,
        str(script),
        str(paths.seurat_input_dir),
        str(paths.figure1d_dir),
        str(seed),
        restart_from,
    ]
    print("Running Seurat Figure 1D workflow...")
    subprocess.run(command, check=True)

    # Step 6: R finishing without complaint is not proof it produced anything. Insist
    # on seeing the finished file before reporting success.
    result = paths.figure1d_dir / "integrated_clustered_umap.rds"
    require_files([result])
    return result


# --- Command-line interface -------------------------------------------------
# What can be typed after the filename when running this stage.
def parse_args() -> argparse.Namespace:
    """Read Stage 03 command-line options."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--workspace",
        default=default_workspace(),
        help="Generated-data workspace (default: checkout-local Age_trajectory)",
    )
    parser.add_argument("--rscript", help="Optional explicit path to Rscript")
    # Several steps ahead involve randomness — where the clustering starts, how the
    # map is initialized. Fixing the seed makes the same run give the same figure.
    parser.add_argument("--seed", type=int, default=1234, help="Random seed")
    # Restarting names the earliest step to recompute, so re-tuning clustering
    # does not throw away hours of SCTransform and integration work.
    parser.add_argument(
        "--restart-from",
        choices=RESTART_LEVELS,
        default="none",
        help=(
            "Earliest step to recompute: none (resume, default); "
            "clustering (minutes); pca (~30 min); integration (hours); "
            "all (many hours, redoes per-sample SCTransform)"
        ),
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Deprecated alias for --restart-from all",
    )
    return parser.parse_args()


def main() -> None:
    """Run the complete third stage from this one Python file."""

    # Example: python scripts/03_reproduce_figure1d.py --restart-from clustering
    args = parse_args()
    # --force predates --restart-from and meant "recompute everything".
    restart_from = "all" if args.force else args.restart_from
    if args.force and args.restart_from != "none":
        raise SystemExit("Use either --force or --restart-from, not both")
    paths = PipelinePaths.from_workspace(args.workspace)
    result = run_seurat_figure1d(
        paths,
        rscript=args.rscript,
        seed=args.seed,
        restart_from=restart_from,
    )
    print(f"Integrated Seurat result: {result}")


if __name__ == "__main__":
    main()
