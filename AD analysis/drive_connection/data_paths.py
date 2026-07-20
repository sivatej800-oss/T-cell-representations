from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Iterable


CONNECTION_ROOT = Path(__file__).resolve().parent
PROJECT_ROOT = CONNECTION_ROOT.parent
DEFAULT_DATA_DIR = PROJECT_ROOT / "AD datasets"
ENV_DATA_DIR = "AD_DATA_DIR"
LOCAL_MANIFEST_PATH = CONNECTION_ROOT / "drive_data_manifest.local.json"
EXAMPLE_MANIFEST_PATH = CONNECTION_ROOT / "drive_data_manifest.example.json"
FALLBACK_DATASET_FILES = (
    "GSE200164_metadata.csv",
    "GSE200164_counts.csv",
    "GSE200164_lognorm_expression.csv",
)


def get_data_dir() -> Path:
    configured = os.environ.get(ENV_DATA_DIR)
    if configured:
        return Path(configured).expanduser()
    return DEFAULT_DATA_DIR


def dataset_file(file_name: str) -> Path:
    return get_data_dir() / file_name


def manifest_path() -> Path:
    if LOCAL_MANIFEST_PATH.exists():
        return LOCAL_MANIFEST_PATH
    return EXAMPLE_MANIFEST_PATH


def load_drive_manifest() -> dict:
    path = manifest_path()
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _folder_manifest(title: str, manifest: dict | None = None) -> dict:
    manifest = manifest or load_drive_manifest()
    for folder in manifest.get("folders", []):
        if folder.get("title") == title:
            return folder
    return {}


def expected_dataset_files(manifest: dict | None = None) -> list[str]:
    folder = _folder_manifest("AD datasets", manifest)
    files = [file_info["title"] for file_info in folder.get("files", [])]
    return files or list(FALLBACK_DATASET_FILES)


def missing_dataset_files(file_names: Iterable[str] | None = None) -> list[Path]:
    names = list(file_names or expected_dataset_files())
    return [dataset_file(name) for name in names if not dataset_file(name).exists()]


def require_dataset_files(file_names: Iterable[str] | None = None) -> None:
    missing = missing_dataset_files(file_names)
    if not missing:
        return

    manifest = load_drive_manifest()
    drive_folder = _folder_manifest("AD datasets", manifest)
    drive_url = drive_folder.get("url") or manifest.get("root_folder", {}).get("url")
    missing_names = "\n".join(f"- {path}" for path in missing)
    help_lines = [
        "Required dataset files are not available locally.",
        f"Expected data directory: {get_data_dir()}",
        "Missing:",
        missing_names,
        f"Set {ENV_DATA_DIR} to your synced Google Drive 'AD datasets' folder, or sync/link it to the expected directory.",
    ]
    if drive_url:
        help_lines.append(f"Drive folder: {drive_url}")
    raise FileNotFoundError("\n".join(help_lines))
