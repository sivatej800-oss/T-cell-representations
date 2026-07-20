from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd


ANALYSIS_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ANALYSIS_ROOT))

from drive_connection.data_paths import dataset_file, require_dataset_files  # noqa: E402


METADATA_FILE = "GSE200164_metadata.csv"
RAW_COUNTS_FILE = "GSE200164_counts.csv"


def main() -> None:
    require_dataset_files([METADATA_FILE, RAW_COUNTS_FILE])

    metadata_preview = pd.read_csv(dataset_file(METADATA_FILE), nrows=5)
    raw_counts_preview = pd.read_csv(dataset_file(RAW_COUNTS_FILE), nrows=5)

    print("Metadata first 5 rows:")
    print(metadata_preview)
    print()
    print("Raw counts first 5 rows:")
    print(raw_counts_preview)


if __name__ == "__main__":
    main()
