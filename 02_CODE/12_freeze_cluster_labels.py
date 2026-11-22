#!/usr/bin/env python3
"""Stage 12: freeze the cluster partitions the biological comparison will judge.

Stages 06 through 11 asked which representation draws better-separated groups, and
answered in the representations' own terms — silhouette, WCSS, donor mixing. None of
those say what the groups *are*. A partition can be geometrically clean and still cut
across cell types, and a partition can look loose and still land on real populations.

The biological test is a "same RNA, different cluster assignments" experiment. PCA and
scGPT decide which cells are grouped together; the original RNA expression says what
those groups are. Because the marker evidence comes from the RNA assay and not from
either representation, it is independent of the thing being evaluated — neither arm
gets to score itself.

That independence only holds if the partitions stop moving. From here on the labels
are an input, never an output: no resolution retuning, no re-clustering, no dropping
of cells. This stage writes the labels to their own file so that every later stage
reads the same frozen partition, and records a checksum so a later stage can prove
it is reading what this stage wrote.

The cluster numbers carry no biological meaning. Cluster 3 in one arm has nothing to
do with cluster 3 in the other; Louvain numbers by descending size, so the numbering
is an artefact of how many cells happened to fall where. What is frozen is the
grouping, not the labels' names.

Four partitions are frozen. The main comparison is the authors' published pipeline
(Stage 03: per-sample SCTransform, anchor integration, 15 PCs) against zero-shot scGPT
with no correction at all (Stage 06's scgpt_raw). The two Harmony arms were the first
comparison and stay frozen so that the work built on them can still be re-run; the
columns they occupied are never changed when partitions are added.

This stage performs no clustering and makes no biological claim. It selects, checks,
and freezes.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd


# --- What is being frozen ----------------------------------------------------
# The Harmony arms are named as Stage 06 named them. The PCA arm is the 15-dimension
# one, the number the authors used, matching Stages 07 and 11. The 30-dimension arm
# stays in Stage 06 only as a sensitivity check.
PCA_ARM = "pca_harmony"
SCGPT_ARM = "scgpt_harmony"
# The main comparison. scgpt_raw is Stage 06's, not Stage 04's: both cluster the same
# embedding with the same settings, but Stage 04's run predates a change of linear
# algebra library on this machine and its partition differs slightly (ARI 0.98, one
# 344-cell group split off). Stage 06's is the current run and the one Stage 17 reads.
SCGPT_RAW_ARM = "scgpt_raw"

# The names the frozen file uses, which are what later stages will refer to. Stage 03's
# partition is not a Stage 06 arm and is read from Stage 03's own output.
FROZEN_NAMES = {
    PCA_ARM: "pca15_harmony_cluster",
    SCGPT_ARM: "scgpt_harmony_cluster",
    SCGPT_RAW_ARM: "scgpt_raw_cluster",
}
STAGE03_COLUMN = "seurat_cluster"
STAGE03_NAME = "authors_pipeline_cluster"
ALL_FROZEN = [*FROZEN_NAMES.values(), STAGE03_NAME]

EXPECTED_CELLS = 70391


@dataclass(frozen=True)
class PipelinePaths:
    """All Stage 12 input and output paths below one workspace directory."""

    workspace: Path

    @classmethod
    def from_workspace(cls, workspace: str | Path) -> "PipelinePaths":
        return cls(Path(workspace).expanduser().resolve())

    @property
    def benchmark_dir(self) -> Path:
        return self.workspace / "results" / "06_harmony_benchmark"

    @property
    def stage03_clusters_file(self) -> Path:
        return (
            self.workspace / "results" / "03_seurat_reproduction"
            / "GSE200164_umap_clusters.csv.gz"
        )

    @property
    def seurat_input_dir(self) -> Path:
        return self.workspace / "results" / "03_seurat_reproduction" / "seurat_input"

    @property
    def frozen_dir(self) -> Path:
        return self.workspace / "results" / "12_frozen_labels"

    def create_directories(self) -> None:
        self.frozen_dir.mkdir(parents=True, exist_ok=True)


def default_workspace() -> Path:
    """Use the repository's Age_trajectory directory unless --workspace is given."""

    return Path(__file__).resolve().parents[1] / "Age_trajectory"


def require_files(paths: Iterable[Path]) -> None:
    """Stop before doing any work if an input this stage depends on is absent."""

    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required inputs:\n  " + "\n  ".join(missing))


# --- Step 1: read the partitions out of Stages 03 and 06 ---------------------
# Stage 06 wrote all five arms into one table and only three are taken; taking exactly
# those is what makes the freeze a decision rather than a copy. Stage 03's partition
# is joined by barcode, never by position.
def load_partitions(paths: PipelinePaths) -> pd.DataFrame:
    """Take the frozen partitions from Stage 06's table and Stage 03's output."""

    assignments_file = paths.benchmark_dir / "benchmark_cluster_assignments.csv.gz"
    require_files([assignments_file, paths.stage03_clusters_file])

    assignments = pd.read_csv(assignments_file, dtype=str)
    for arm in FROZEN_NAMES:
        if arm not in assignments.columns:
            raise ValueError(f"Stage 06 output does not contain the arm {arm}")
    frozen = assignments[["cell_barcode", *FROZEN_NAMES]].rename(columns=FROZEN_NAMES)

    stage03 = pd.read_csv(paths.stage03_clusters_file, dtype=str)
    stage03 = stage03[["cell_barcode", STAGE03_COLUMN]].rename(
        columns={STAGE03_COLUMN: STAGE03_NAME}
    )
    if len(stage03) != len(frozen) or set(stage03["cell_barcode"]) != set(
        frozen["cell_barcode"]
    ):
        raise ValueError("Stage 03 and Stage 06 do not describe the same cells")
    frozen = frozen.merge(stage03, on="cell_barcode", how="left", validate="one_to_one")
    return frozen


# --- Step 2: confirm both partitions describe the same cells as the RNA ------
# The whole design rests on the two arms labelling one fixed set of cells, and on the
# RNA the markers will be read from being that same set. A silent mismatch here would
# not crash anything downstream; it would attach marker expression to the wrong cells
# and produce a clean, wrong annotation. So the barcodes are compared against the
# matrix the RNA will be read from, not merely counted.
def verify_cells_match_rna(frozen: pd.DataFrame, paths: PipelinePaths) -> None:
    """Check the frozen labels cover the RNA matrix's cells, exactly and in order."""

    barcodes_file = paths.seurat_input_dir / "barcodes.tsv"
    require_files([barcodes_file])
    rna_barcodes = barcodes_file.read_text(encoding="utf-8").split()

    if len(frozen) != EXPECTED_CELLS:
        raise ValueError(
            f"Expected {EXPECTED_CELLS} cells to freeze, found {len(frozen)}"
        )
    if frozen["cell_barcode"].duplicated().any():
        raise ValueError("Stage 06's assignments contain duplicate cell barcodes")
    if frozen[ALL_FROZEN].isna().any().any():
        raise ValueError("A cell is missing a cluster label in one of the arms")
    if len(rna_barcodes) != EXPECTED_CELLS:
        raise ValueError(
            f"The RNA matrix describes {len(rna_barcodes)} cells, "
            f"not the expected {EXPECTED_CELLS}"
        )
    if list(frozen["cell_barcode"]) != rna_barcodes:
        if set(frozen["cell_barcode"]) == set(rna_barcodes):
            raise ValueError(
                "The frozen labels and the RNA matrix hold the same cells in "
                "different orders; later stages must join on barcode, not position"
            )
        raise ValueError(
            "The frozen labels and the RNA matrix do not describe the same cells"
        )

    print(
        f"Cells: {len(frozen)}, matching the RNA matrix barcode for barcode, in order."
    )


# --- Step 3: report the partitions being frozen ------------------------------
# Printed so the freeze is inspectable at the moment it happens. The cross-tabulation
# is the substance of the main comparison: the same cells partitioned two ways, and
# where the two disagree is where the marker evidence will have to decide.
def describe_partitions(frozen: pd.DataFrame) -> dict[str, object]:
    """State each partition's cluster sizes and how far the main two disagree."""

    sources = {name: f"Stage 06 arm {arm}" for arm, name in FROZEN_NAMES.items()}
    sources[STAGE03_NAME] = f"Stage 03 {STAGE03_COLUMN}"
    report: dict[str, object] = {}
    for column in ALL_FROZEN:
        sizes = frozen[column].value_counts()
        print(
            f"\n{column} ({sources[column]}): {len(sizes)} clusters, "
            f"{sizes.min()}-{sizes.max()} cells per cluster"
        )
        report[column] = {
            "source": sources[column],
            "n_clusters": int(len(sizes)),
            "cluster_sizes": {str(k): int(v) for k, v in sizes.sort_index().items()},
        }

    overlap = pd.crosstab(frozen[STAGE03_NAME], frozen[FROZEN_NAMES[SCGPT_RAW_ARM]])
    print("\nCells shared between clusters (rows authors' pipeline, columns scGPT raw):")
    print(overlap.to_string())
    print(
        "\nCluster numbers are not comparable across arms. Louvain numbers by "
        "descending size, so this table shows which cells moved, not which "
        "clusters correspond."
    )
    return report


# --- Step 4: write the freeze and a checksum that proves it ------------------
# The checksum is over the label columns only, not the file bytes, so that it stays
# valid if the file is rewritten identically and breaks if any cell's assignment
# changes. Later stages record it alongside their results, which is what makes their
# claim "these markers describe this partition" checkable after the fact.
def write_frozen_labels(frozen: pd.DataFrame, paths: PipelinePaths) -> tuple[Path, str]:
    """Write the frozen labels and return the path and their checksum."""

    # A partition may be added to the freeze, but one already frozen may never change:
    # every column the existing file holds must come back cell for cell identical.
    frozen_file = paths.frozen_dir / "frozen_cluster_labels.csv.gz"
    if frozen_file.exists():
        existing = pd.read_csv(frozen_file, dtype=str)
        unknown = set(existing.columns) - set(frozen.columns)
        changed = [
            column for column in existing.columns
            if column in frozen.columns
            and not existing[column].equals(frozen[column].astype(str))
        ]
        if unknown or changed:
            raise RuntimeError(
                f"{frozen_file} already exists and holds different labels "
                f"({sorted(unknown | set(changed))}). The partitions are frozen; if "
                "they genuinely need to change, delete this file deliberately and "
                "re-run every stage that read it."
            )
        if list(existing.columns) == list(frozen.columns):
            print(f"\n{frozen_file.name} already exists and is identical; left in place.")
        else:
            added = [c for c in frozen.columns if c not in existing.columns]
            frozen.to_csv(frozen_file, index=False)
            print(f"\nAdded {added} to {frozen_file.name}; existing columns unchanged.")
    else:
        frozen.to_csv(frozen_file, index=False)

    columns = ["cell_barcode", *ALL_FROZEN]
    payload = frozen[columns].to_csv(index=False).encode("utf-8")
    checksum = hashlib.sha256(payload).hexdigest()
    return frozen_file, checksum


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--workspace",
        default=str(default_workspace()),
        help="Directory holding data/ and results/ (default: Age_trajectory)",
    )
    parser.add_argument(
        "--verify-only",
        action="store_true",
        help="Run the checks and print the partitions, then stop without writing",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    paths = PipelinePaths.from_workspace(args.workspace)
    paths.create_directories()

    frozen = load_partitions(paths)
    verify_cells_match_rna(frozen, paths)
    report = describe_partitions(frozen)

    if args.verify_only:
        print("\nVerification complete; stopping before writing as requested.")
        return

    frozen_file, checksum = write_frozen_labels(frozen, paths)

    audit = {
        "stage": "12_freeze_cluster_labels",
        "sources": [
            "results/06_harmony_benchmark/benchmark_cluster_assignments.csv.gz",
            "results/03_seurat_reproduction/GSE200164_umap_clusters.csv.gz",
        ],
        "arms_frozen": {**FROZEN_NAMES, STAGE03_COLUMN: STAGE03_NAME},
        "n_cells": int(len(frozen)),
        "label_checksum_sha256": checksum,
        "partitions": report,
    }
    audit_path = paths.frozen_dir / "frozen_labels_audit.json"
    audit_path.write_text(json.dumps(audit, indent=2), encoding="utf-8")

    print(
        f"\nStage 12 complete. The partitions are frozen.\n"
        f"  Labels: {frozen_file}\n"
        f"  Audit: {audit_path}\n"
        f"  Checksum: {checksum}"
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as error:  # noqa: BLE001
        print(f"Stage 12 failed: {error}", file=sys.stderr)
        raise SystemExit(1) from error
