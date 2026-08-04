"""Implementation of Stage 03's Python-to-R handoff.

Python decides whether inputs/checkpoints are sufficient and audits the R environment.
The bundled R script then performs or resumes the actual Seurat computation.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from .common import require_files
from .config import PipelinePaths
from .r_environment import audit_r_environment


def bundled_r_script() -> Path:
    """Return the R backend stored beside this Python module."""

    return Path(__file__).with_name("run_seurat_figure1d.R")


def run_seurat_figure1d(
    paths: PipelinePaths,
    *,
    rscript: str | None = None,
    seed: int = 1234,
    force: bool = False,
) -> Path:
    """Validate prerequisites, invoke R, and return the completed Seurat RDS path."""

    paths.create_directories()
    # A from-scratch run crosses from Python to R through these four aligned files.
    raw_inputs = [
        paths.seurat_input_dir / "counts.mtx",
        paths.seurat_input_dir / "genes.tsv",
        paths.seurat_input_dir / "barcodes.tsv",
        paths.seurat_input_dir / "metadata.csv",
    ]
    # Checkpoints are ordered here for readability, but R makes the final choice:
    # completed result first, then post-PCA, then integrated, then raw inputs.
    reusable_checkpoints = [
        paths.figure1d_dir / "integrated_clustered_umap.rds",
        paths.figure1d_dir / "GSE200164_after_PCA.rds",
        paths.figure1d_dir / "GSE200164_seurat_integrated.rds",
    ]
    if not any(path.exists() for path in reusable_checkpoints):
        require_files(raw_inputs)
    else:
        print("Existing final/integrated/PCA checkpoint found; raw Seurat export is not required.")
    # An explicit --rscript value is useful when multiple R installations exist.
    executable = rscript or shutil.which("Rscript")
    if executable is None:
        raise RuntimeError("Rscript was not found on PATH")
    audit = audit_r_environment(paths, rscript=executable)
    if audit["missing_r_packages"]:
        raise RuntimeError(
            "Cannot run Seurat; missing R packages: "
            + ", ".join(audit["missing_r_packages"])
        )

    script = bundled_r_script()
    require_files([script])
    # Use an argument list rather than a shell command. This preserves paths containing
    # spaces and avoids shell-quoting or command-injection problems.
    command = [
        executable,
        str(script),
        str(paths.seurat_input_dir),
        str(paths.figure1d_dir),
        str(seed),
        str(force).lower(),
    ]
    print("Running Seurat Figure 1D workflow...")
    subprocess.run(command, check=True)
    # Do not report success merely because R exited; require the promised final file.
    result = paths.figure1d_dir / "integrated_clustered_umap.rds"
    require_files([result])
    return result
