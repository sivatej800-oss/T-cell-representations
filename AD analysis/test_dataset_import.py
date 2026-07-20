from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent
DATA_DIR = PROJECT_ROOT / "AD datasets"

metadata_path = DATA_DIR / "GSE200164_metadata.csv"
counts_path = DATA_DIR / "GSE200164_counts.csv"
lognorm_path = DATA_DIR / "GSE200164_lognorm_expression.csv"

metadata = pd.read_csv(metadata_path)
counts_preview = pd.read_csv(counts_path, nrows=5)
lognorm_preview = pd.read_csv(lognorm_path, nrows=5)

print("metadata shape:", metadata.shape)
print(metadata.head())
print("counts preview shape:", counts_preview.shape)
print("log-normalized expression preview shape:", lognorm_preview.shape)
