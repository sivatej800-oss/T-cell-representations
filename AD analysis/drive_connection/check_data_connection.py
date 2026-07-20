from __future__ import annotations

from data_paths import (
    DEFAULT_DATA_DIR,
    ENV_DATA_DIR,
    expected_dataset_files,
    get_data_dir,
    load_drive_manifest,
    missing_dataset_files,
)


def _human_size(size: int | str | None) -> str:
    if size in (None, ""):
        return "unknown size"
    size_float = float(size)
    units = ("B", "KB", "MB", "GB", "TB")
    unit_index = 0
    while size_float >= 1024 and unit_index < len(units) - 1:
        size_float /= 1024
        unit_index += 1
    if unit_index == 0:
        return f"{int(size_float)} {units[unit_index]}"
    return f"{size_float:.2f} {units[unit_index]}"


def _dataset_folder(manifest: dict) -> dict:
    for folder in manifest.get("folders", []):
        if folder.get("title") == "AD datasets":
            return folder
    return {}


def main() -> int:
    manifest = load_drive_manifest()
    root = manifest.get("root_folder", {})
    dataset_folder = _dataset_folder(manifest)
    data_dir = get_data_dir()

    print(f"Drive root: {root.get('title', 'not configured')}")
    if root.get("url"):
        print(f"Drive URL: {root['url']}")
    print(f"Local data directory: {data_dir}")
    print()

    known_files = {
        file_info.get("title"): file_info
        for file_info in dataset_folder.get("files", [])
    }
    missing = set(missing_dataset_files())
    for file_name in expected_dataset_files(manifest):
        path = data_dir / file_name
        status = "OK" if path.exists() else "MISSING"
        size = _human_size(known_files.get(file_name, {}).get("size"))
        print(f"{status:7} {file_name} ({size})")

    if missing:
        print()
        print("Data is connected in the manifest, but the files are not synced locally yet.")
        print(f"Sync the Drive folder named '{root.get('title', 'Summer 2026_ AD_Research')}' with Google Drive for desktop.")
        print(f"Then either place/link its 'AD datasets' folder at '{DEFAULT_DATA_DIR}'")
        print(f"or set {ENV_DATA_DIR} to the synced folder path before running analysis.")
        return 1

    print()
    print("Data connection is ready.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
