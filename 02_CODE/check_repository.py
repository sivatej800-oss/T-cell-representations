#!/usr/bin/env python3
"""Validate repository hygiene without loading scientific project dependencies.

This lightweight check is suitable for both local use and GitHub Actions. It prevents
large generated data from entering Git and confirms archived notebooks remain valid
JSON documents.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


# GitHub rejects files at 100 MiB, so 95 MiB leaves a small safety margin.
ROOT = Path(__file__).resolve().parents[1]
MAX_FILE_BYTES = 95 * 1024 * 1024
GENERATED_SUFFIXES = {
    ".h5ad",
    ".hdf5",
    ".loom",
    ".mtx",
    ".npy",
    ".npz",
    ".rda",
    ".rdata",
    ".rds",
}


def tracked_files() -> list[Path]:
    """Return staged/tracked files plus visible untracked files for pre-commit checks."""

    result = subprocess.run(
        ["git", "ls-files", "--cached", "--others", "--exclude-standard"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return [ROOT / line for line in result.stdout.splitlines() if line]


def main() -> int:
    """Collect all hygiene errors, print them together, and return a shell status."""

    errors: list[str] = []
    files = tracked_files()

    for path in files:
        if not path.is_file():
            continue
        relative = path.relative_to(ROOT)
        size = path.stat().st_size
        if size > MAX_FILE_BYTES:
            errors.append(f"oversized file ({size / 1024**2:.1f} MiB): {relative}")
        if path.suffix.lower() in GENERATED_SUFFIXES:
            errors.append(f"generated data artifact: {relative}")

    # A notebook can be valid JSON but still be incomplete, so check its two required
    # top-level fields after parsing.
    notebooks = sorted(
        path for path in files if path.suffix == ".ipynb" and path.is_file()
    )
    for notebook in notebooks:
        try:
            payload = json.loads(notebook.read_text(encoding="utf-8"))
        except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
            errors.append(f"invalid notebook {notebook.relative_to(ROOT)}: {exc}")
            continue
        if "cells" not in payload or "nbformat" not in payload:
            errors.append(f"incomplete notebook structure: {notebook.relative_to(ROOT)}")

    if errors:
        print("Repository check failed:")
        for error in errors:
            print(f"  - {error}")
        return 1

    print(f"Repository check passed ({len(files)} files, {len(notebooks)} notebooks).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
