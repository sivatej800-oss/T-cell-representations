"""Small file-handling helpers shared by all three pipeline stages.

Keeping these operations here ensures every stage checks files, writes JSON, and saves
large artifacts in the same safe way.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import shutil
import subprocess
from pathlib import Path
from typing import Any, Iterable


def require_files(paths: Iterable[Path]) -> None:
    """Raise one readable error listing every required file that is missing."""

    missing = [path for path in paths if not path.is_file()]
    if missing:
        formatted = "\n".join(f"  - {path}" for path in missing)
        raise FileNotFoundError(f"Required files were not found:\n{formatted}")


def write_json(path: Path, payload: dict[str, Any]) -> None:
    """Write a human-readable JSON audit record, creating its parent folder first."""

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    """Calculate a file checksum without reading a large file into memory at once."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def validate_gzip(path: Path) -> None:
    """Read the complete gzip stream so CRC and truncation errors are detected."""

    with gzip.open(path, "rb") as handle:
        while handle.read(8 * 1024 * 1024):
            pass


def download_with_resume(url: str, destination: Path, *, force: bool = False) -> None:
    """Download with curl, retaining a partial file if transfer is interrupted.

    curl's ``--continue-at -`` option resumes the ``.part`` file on the next run.
    The partial file is renamed to the requested destination only after curl succeeds.
    """

    if destination.exists() and not force:
        print(f"Using existing download: {destination}")
        return

    curl = shutil.which("curl")
    if curl is None:
        raise RuntimeError("curl is required for GEO downloads but was not found on PATH")

    destination.parent.mkdir(parents=True, exist_ok=True)
    # Never expose an incomplete download under the final filename.
    partial = destination.with_suffix(destination.suffix + ".part")
    command = [
        curl,
        "--location",
        "--fail",
        "--retry",
        "5",
        "--retry-delay",
        "5",
        "--continue-at",
        "-",
        "--output",
        str(partial),
        url,
    ]
    subprocess.run(command, check=True)
    partial.replace(destination)


def atomic_h5ad_write(adata: Any, destination: Path, *, force: bool = False) -> None:
    """Save AnnData without leaving a corrupt final file after an interruption.

    Data are first written to a neighboring ``.tmp`` file. A successful write is then
    renamed into place, which is an atomic operation on the same filesystem.
    """

    if destination.exists() and not force:
        print(f"Using existing AnnData file: {destination}")
        return

    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    # A leftover temporary file means an earlier save did not finish; it is never used
    # as a checkpoint and can be replaced safely here.
    if temporary.exists():
        temporary.unlink()
    adata.write_h5ad(
        temporary,
        compression="gzip",
        convert_strings_to_categoricals=False,
    )
    temporary.replace(destination)
