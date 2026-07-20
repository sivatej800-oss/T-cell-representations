from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = PROJECT_ROOT / "data" / "ad_datasets"


METADATA_FILE = "GSE200164_metadata.csv"
RAW_COUNTS_FILE = "GSE200164_counts.csv"


def dataset_file(file_name: str) -> Path:
    return DATA_DIR / file_name


def require_dataset_files(file_names: list[str]) -> None:
    missing = [dataset_file(file_name) for file_name in file_names if not dataset_file(file_name).exists()]
    if not missing:
        return

    missing_list = "\n".join(f"- {path}" for path in missing)
    raise FileNotFoundError(
        "\n".join(
            [
                "Required dataset files are not available locally.",
                f"Expected data directory: {DATA_DIR}",
                "Missing:",
                missing_list,
                "Upload the missing files into data/ad_datasets/.",
            ]
        )
    )


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
