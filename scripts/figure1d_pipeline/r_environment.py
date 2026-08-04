"""Check that R, required packages, memory, and disk are ready for Seurat.

Stages 02 and 03 both use this module so environment failures are reported before a
long-running integration begins.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys

import psutil

from .common import write_json
from .config import PipelinePaths


# Keep this list synchronized with imports in run_seurat_figure1d.R.
REQUIRED_R_PACKAGES = ("Seurat", "Matrix", "glmGamPoi", "ggplot2")


def audit_r_environment(
    paths: PipelinePaths,
    *,
    rscript: str | None = None,
) -> dict[str, object]:
    """Audit the runtime and save the result as reproducibility metadata."""

    executable = rscript or shutil.which("Rscript")
    if executable is None:
        raise RuntimeError("Rscript was not found on PATH")
    # Build a small R expression that reports package availability and version. JSON
    # quoting safely turns each Python package name into an R string literal.
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
    # check=True converts an R execution failure into a Python exception immediately.
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
    # Memory and free-disk values do not enforce a guessed minimum; they document the
    # environment and help diagnose failures of this unusually large analysis.
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
    write_json(paths.figure1d_dir / "figure1d_environment_audit.json", audit)
    print(f"R environment audit: missing packages = {missing or 'none'}")
    return audit
