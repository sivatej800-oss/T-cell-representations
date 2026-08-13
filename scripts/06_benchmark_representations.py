#!/usr/bin/env python3
"""Stage 06: compare the PCA and scGPT descriptions with the batch confound removed.

Stages 03 and 04 are not directly comparable. Stage 03 corrects for patient before
reducing dimensions, using integration anchors. Stage 04 does not correct at all,
because scGPT expects count-scale input and pre-correcting would feed the model
values unlike anything it was trained on. Comparing them therefore compares a
corrected representation against an uncorrected one, and a difference in clustering
could be patient structure surviving in the second rather than anything about the
descriptions.

This stage removes that confound by correcting after the description exists instead
of before, with Harmony, on both arms equally. It produces a 2x2 — each description
corrected and uncorrected — so the effect of the correction can be read separately
from the effect of the description, plus a fifth arm that repeats the corrected PCA
on 30 dimensions rather than 15 to show whether that cut matters.

Neither Stage 03 nor Stage 04 is modified or re-run. Stage 03's clustering remains
the paper reproduction; this stage's PCA arm is a different analysis that happens to
start from the same per-sample normalizations, and its cluster count should not be
compared back to the paper's.

This file verifies inputs and hands them to R. It performs no analysis of its own.
The verification matters because the embedding is 70,391 anonymous rows of 512
numbers: if its rows were ordered differently from the cells, every step downstream
would still run and would still produce a plausible figure, and it would be wrong.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


# --- Required R software -----------------------------------------------------
# Seurat for the clustering, harmony for the correction. glmGamPoi is not needed:
# this stage reads the per-sample normalizations Stage 03 already produced rather
# than recomputing them.
REQUIRED_R_PACKAGES = ("Seurat", "Matrix", "harmony", "ggplot2")

EXPECTED_EMBEDDING_DIMS = 512
EXPECTED_CELLS = 70391
EXPECTED_SAMPLES = 59

# The seed Stages 03 and 04 use. Sharing it means a difference between arms cannot
# be a difference in where the random numbers happened to start.
DEFAULT_SEED = 42


# --- Workspace paths ---------------------------------------------------------
# Stages 03 and 04 each repeat these helpers rather than importing them, so that any
# stage can be read and run on its own. This stage follows that decision.
@dataclass(frozen=True)
class PipelinePaths:
    """All Stage 06 input and output paths below one workspace directory."""

    workspace: Path

    @classmethod
    def from_workspace(cls, workspace: str | Path) -> "PipelinePaths":
        return cls(Path(workspace).expanduser().resolve())

    # The unpacked scGPT bundle: the embedding, its manifest, and its run metadata.
    @property
    def scgpt_dir(self) -> Path:
        return self.workspace / "data" / "scgpt"

    # Stage 03's directory. Read only: the per-sample normalizations are reused as
    # this stage's starting point, and its cell list is used to check correspondence.
    @property
    def pca_arm_dir(self) -> Path:
        return self.workspace / "results" / "figure1d_python"

    @property
    def sct_dir(self) -> Path:
        return self.pca_arm_dir / "sct_samples"

    # This stage's results, kept apart from both earlier arms so nothing collides.
    @property
    def benchmark_dir(self) -> Path:
        return self.workspace / "results" / "benchmark"

    @property
    def r_input_dir(self) -> Path:
        return self.benchmark_dir / "r_input"

    def create_directories(self) -> None:
        """Create this stage's folders, leaving any previous result in place."""

        for directory in (self.benchmark_dir, self.r_input_dir):
            directory.mkdir(parents=True, exist_ok=True)


def default_workspace() -> Path:
    """Use the repository's Age_trajectory directory unless --workspace is given."""

    return Path(__file__).resolve().parents[1] / "Age_trajectory"


def require_files(paths: Iterable[Path]) -> None:
    """Stop before doing any work if an input this stage depends on is absent."""

    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required inputs:\n  " + "\n  ".join(missing))


# --- Step 1: confirm Stage 03's normalizations are present -------------------
# The PCA arm starts from the 59 per-sample SCTransform checkpoints rather than
# recomputing them, which is what keeps this stage to about an hour instead of the
# many hours Stage 03 takes. All 59 must be there: a missing sample would silently
# shrink the dataset and change every count this stage reports.
def verify_sct_checkpoints(paths: PipelinePaths) -> list[str]:
    """Check that every per-sample normalization Stage 03 wrote is available."""

    if not paths.sct_dir.is_dir():
        raise FileNotFoundError(
            f"Stage 03's per-sample normalizations were not found at {paths.sct_dir}. "
            "Run Stage 03 at least as far as SCTransform before this stage."
        )
    checkpoints = sorted(paths.sct_dir.glob("*_sct.rds"))
    if len(checkpoints) != EXPECTED_SAMPLES:
        raise ValueError(
            f"Expected {EXPECTED_SAMPLES} per-sample normalizations, "
            f"found {len(checkpoints)} in {paths.sct_dir}"
        )
    empty = [str(path) for path in checkpoints if path.stat().st_size == 0]
    if empty:
        raise ValueError("Empty normalization checkpoints:\n  " + "\n  ".join(empty))
    sample_names = [path.name.removesuffix("_sct.rds") for path in checkpoints]
    print(f"Per-sample normalizations verified: {len(sample_names)} samples")
    return sample_names


# --- Step 2: prove the embedding describes these cells, in this order --------
def verify_cell_correspondence(
    paths: PipelinePaths,
) -> tuple[np.ndarray, pd.DataFrame, dict[str, object]]:
    """Load the embedding and prove its rows are this pipeline's cells, in order."""

    metadata_file = paths.scgpt_dir / "scgpt_run_metadata.json"
    embedding_file = paths.scgpt_dir / "X_scgpt.npy"
    manifest_file = paths.scgpt_dir / "scgpt_cell_manifest.csv.gz"
    require_files([metadata_file, embedding_file, manifest_file])
    run_metadata = json.loads(metadata_file.read_text(encoding="utf-8"))

    # A non-finite value would propagate silently through the neighbor graph, and a
    # zero-length row would have no direction at all.
    embedding = np.load(embedding_file)
    if list(embedding.shape) != run_metadata["output_shape"]:
        raise ValueError("Embedding shape disagrees with the recorded run metadata")
    if embedding.shape != (EXPECTED_CELLS, EXPECTED_EMBEDDING_DIMS):
        raise ValueError(f"Expected a {EXPECTED_CELLS} x {EXPECTED_EMBEDDING_DIMS} embedding")
    if not np.isfinite(embedding).all():
        raise ValueError("Embedding contains NaN or infinite values")
    if np.linalg.norm(embedding, axis=1).min() <= 0:
        raise ValueError("At least one cell has an all-zero embedding")

    manifest = pd.read_csv(manifest_file)
    if len(manifest) != embedding.shape[0]:
        raise ValueError("Manifest row count does not match the embedding")
    if not (manifest["embedding_row"].to_numpy() == np.arange(len(manifest))).all():
        raise ValueError("Manifest is not in embedding row order")
    if not manifest["cell_id"].is_unique:
        raise ValueError("Manifest contains duplicate cell identifiers")

    # The recorded hash of the embedded cell order. Reproducing it confirms row i
    # belongs to cell i, rather than to some other arrangement of the same cells.
    cell_ids = manifest["cell_id"].to_numpy()
    order_hash = hashlib.sha256("\n".join(cell_ids).encode("utf-8")).hexdigest()
    if order_hash != run_metadata["cell_order_sha256"]:
        raise ValueError(
            "Recorded cell-order hash does not match the manifest; the embedding "
            "rows may not correspond to these cells in this order"
        )

    # The cells Stage 03 ran on, which are also the cells behind the per-sample
    # normalizations this stage's PCA arm rebuilds from.
    pca_barcodes_file = paths.pca_arm_dir / "seurat_input" / "barcodes.tsv"
    require_files([pca_barcodes_file])
    pca_barcodes = pca_barcodes_file.read_text(encoding="utf-8").split()
    if set(pca_barcodes) != set(cell_ids):
        raise ValueError("Embedded cells are not the same cells Stage 03 used")
    if manifest["internal_sample_id"].nunique() != EXPECTED_SAMPLES:
        raise ValueError(f"Expected {EXPECTED_SAMPLES} internal sample IDs")

    print(
        f"Cell correspondence verified: {embedding.shape[0]} cells, "
        f"{embedding.shape[1]} dimensions, order hash reproduced"
    )
    return embedding, manifest, run_metadata


# --- Step 3: report the structure the correction will remove -----------------
# Each sample here is one patient, and age and diagnosis are recorded per patient,
# so they are constant across all of that patient's cells. Correcting for patient
# therefore removes all age and all diagnosis variation along with the batch
# differences — necessarily, because those variables have no within-patient
# variation left to survive it.
#
# That is correct for grouping cells by type, which is all this stage does. It is
# printed rather than assumed so that nobody later reuses these coordinates for the
# age trajectory question, where the signal has been removed by construction.
def report_donor_confounding(manifest: pd.DataFrame) -> dict[str, object]:
    """State plainly which variables patient identity determines, and warn."""

    per_donor = manifest.groupby("internal_sample_id")
    nested = [
        column
        for column in ("age", "diagnosis", "sex", "cognitive_group")
        if column in manifest.columns and (per_donor[column].nunique() == 1).all()
    ]
    sizes = manifest["internal_sample_id"].value_counts()
    proportions = (sizes / len(manifest)).to_numpy()
    mixing_ceiling = float(1.0 - np.sum(proportions**2))

    print(
        f"\nDonor structure: {len(sizes)} patients, "
        f"{sizes.min()}-{sizes.max()} cells each"
    )
    print(f"Constant within patient, so removed by the correction: {', '.join(nested)}")
    print(
        f"Cross-patient neighbor fraction under perfect mixing: {mixing_ceiling:.4f} "
        "(the ceiling the mixing score is measured against)"
    )
    print(
        "WARNING: the Harmony arms are valid for cell-type grouping only. Age and "
        "diagnosis analysis must use uncorrected expression within clusters."
    )
    return {
        "n_donors": int(len(sizes)),
        "smallest_donor": int(sizes.min()),
        "largest_donor": int(sizes.max()),
        "variables_removed_by_correction": nested,
        "donor_mixing_ceiling": mixing_ceiling,
    }


# --- Step 4: hand the embedding to R -----------------------------------------
# R has no reader for NumPy's format, so the embedding is rewritten as a flat block
# of 32-bit floats in row order. The values are copied unchanged; this converts how
# the numbers are packaged, not what they are.
def write_r_inputs(
    paths: PipelinePaths, embedding: np.ndarray, manifest: pd.DataFrame
) -> None:
    """Write the embedding and per-cell details in the form the R script reads."""

    embedding_path = paths.r_input_dir / "embedding_float32.bin"
    embedding.astype("<f4", copy=False).tofile(embedding_path)
    expected_bytes = embedding.size * 4
    if embedding_path.stat().st_size != expected_bytes:
        raise ValueError("Embedding binary was not written at its expected size")

    # The same column name Stages 03 and 04 use for a cell, so the output tables of
    # all three stages line up without translation.
    metadata = manifest.rename(columns={"cell_id": "cell_barcode"})
    metadata.to_csv(paths.r_input_dir / "metadata.csv", index=False)
    print(f"Wrote R inputs to {paths.r_input_dir}")


# --- R preflight -------------------------------------------------------------
# Both arms must run the same Seurat, or a difference between them could come from
# the software rather than from the representation. The versions are recorded for
# the same reason Stages 03 and 04 record them.
def audit_r_environment(paths: PipelinePaths, rscript: str | None = None) -> dict[str, str]:
    """Confirm the R packages exist and record what versions ran."""

    executable = rscript or shutil.which("Rscript")
    if executable is None:
        raise RuntimeError("Rscript was not found on PATH")

    program = (
        "versions <- vapply(c("
        + ", ".join(f'"{package}"' for package in REQUIRED_R_PACKAGES)
        + "), function(p) tryCatch(as.character(packageVersion(p)), "
        'error = function(e) "MISSING"), character(1));'
        "cat(paste(names(versions), versions, sep = '=', collapse = '\\n'))"
    )
    result = subprocess.run(
        [executable, "-e", program], check=True, capture_output=True, text=True
    )
    versions = dict(
        line.split("=", 1) for line in result.stdout.strip().splitlines() if "=" in line
    )
    missing = [name for name, version in versions.items() if version == "MISSING"]
    if missing:
        raise RuntimeError(
            "Missing required R packages: "
            + ", ".join(missing)
            + '\nInstall with: install.packages(c("harmony"))'
        )
    print("R packages: " + ", ".join(f"{k} {v}" for k, v in versions.items()))
    return versions


def bundled_r_script() -> Path:
    """The R script that sits beside this one and does the actual work."""

    script = Path(__file__).resolve().with_suffix(".R")
    if not script.exists():
        raise FileNotFoundError(f"Companion R script not found at {script}")
    return script


def run_benchmark(paths: PipelinePaths, seed: int, rscript: str | None = None) -> Path:
    """Start the R stage and confirm it produced the result it claims to."""

    executable = rscript or shutil.which("Rscript")
    if executable is None:
        raise RuntimeError("Rscript was not found on PATH")

    command = [
        executable,
        "--vanilla",
        str(bundled_r_script()),
        str(paths.r_input_dir),
        str(paths.sct_dir),
        str(paths.benchmark_dir),
        str(seed),
        str(EXPECTED_CELLS),
        str(EXPECTED_EMBEDDING_DIMS),
    ]
    print("\nStarting R:\n  " + " ".join(command) + "\n")
    subprocess.run(command, check=True)

    summary = paths.benchmark_dir / "benchmark_arm_summary.csv"
    if not summary.exists():
        raise RuntimeError("R finished without writing the arm summary")
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--workspace",
        default=str(default_workspace()),
        help="Directory holding data/ and results/ (default: Age_trajectory)",
    )
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED, help="Random seed")
    parser.add_argument("--rscript", help="Optional explicit path to Rscript")
    parser.add_argument(
        "--verify-only",
        action="store_true",
        help="Run the input checks and the confounding report, then stop",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    paths = PipelinePaths.from_workspace(args.workspace)
    paths.create_directories()

    verify_sct_checkpoints(paths)
    embedding, manifest, run_metadata = verify_cell_correspondence(paths)
    confounding = report_donor_confounding(manifest)

    if args.verify_only:
        print("\nVerification complete; stopping before the analysis as requested.")
        return

    versions = audit_r_environment(paths, rscript=args.rscript)
    write_r_inputs(paths, embedding, manifest)
    summary = run_benchmark(paths, seed=args.seed, rscript=args.rscript)

    audit = {
        "stage": "06_benchmark_representations",
        "seed": args.seed,
        "r_packages": versions,
        "scgpt_model": run_metadata.get("model"),
        "scgpt_cell_order_sha256": run_metadata.get("cell_order_sha256"),
        "scgpt_gene_vocabulary_coverage": run_metadata.get("gene_vocabulary_coverage"),
        "donor_structure": confounding,
    }
    audit_path = paths.benchmark_dir / "benchmark_environment_audit.json"
    audit_path.write_text(json.dumps(audit, indent=2), encoding="utf-8")

    print(f"\nStage 06 complete.\n  Summary: {summary}\n  Audit: {audit_path}")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:  # noqa: BLE001
        print(f"Stage 06 failed: {error}", file=sys.stderr)
        raise SystemExit(1) from error
