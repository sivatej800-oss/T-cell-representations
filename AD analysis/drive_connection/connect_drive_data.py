from __future__ import annotations

import os
from pathlib import Path

from data_paths import PROJECT_ROOT, load_drive_manifest


def _google_drive_project_dirs(root_title: str) -> list[Path]:
    cloud_storage = Path.home() / "Library" / "CloudStorage"
    if not cloud_storage.exists():
        return []
    candidates = []
    for drive_root in sorted(cloud_storage.glob("GoogleDrive-*")):
        project_dir = drive_root / "My Drive" / root_title
        if project_dir.exists():
            candidates.append(project_dir)
    return candidates


def _link_folder(source: Path, destination: Path) -> bool:
    if destination.is_symlink():
        current = destination.resolve(strict=False)
        if current == source:
            print(f"OK      {destination} -> {source}")
            return True
        destination.unlink()
    elif destination.exists():
        print(f"SKIP    {destination} already exists and is not a symlink")
        return False

    os.symlink(source, destination, target_is_directory=True)
    print(f"LINKED  {destination} -> {source}")
    return True


def main() -> int:
    manifest = load_drive_manifest()
    root_title = manifest.get("root_folder", {}).get("title", "Summer 2026_ AD_Research")
    project_dirs = _google_drive_project_dirs(root_title)
    if not project_dirs:
        print(f"Could not find a local Google Drive folder named '{root_title}'.")
        print("Open Google Drive for desktop and make that Drive folder available locally, then rerun this script.")
        return 1

    drive_project_dir = project_dirs[0]
    all_linked = True
    for folder in manifest.get("folders", []):
        folder_title = folder.get("title")
        local_path = folder.get("local_path")
        if not folder_title or not local_path:
            continue
        source = drive_project_dir / folder_title
        destination = PROJECT_ROOT / local_path
        if not source.exists():
            print(f"MISSING {source}")
            all_linked = False
            continue
        all_linked = _link_folder(source, destination) and all_linked

    return 0 if all_linked else 1


if __name__ == "__main__":
    raise SystemExit(main())
