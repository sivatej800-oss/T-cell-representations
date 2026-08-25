#!/usr/bin/env python3
"""Stage 13: ask the RNA what the two frozen partitions are.

Stage 12 froze two partitions of the same 70,391 cells -- one from Harmony-corrected
PCA at 30 dimensions, one from the Harmony-corrected scGPT embedding. Everything up to
here compared them geometrically: silhouette, WCSS, donor mixing. None of it says
whether either partition's communities are cell types.

This stage runs the comparison the way the paper's authors ran their annotation, with
one change that makes it a comparison rather than two separate analyses. The authors
clustered on the integrated assay and then switched to expression to interpret the
result. Here both partitions are interpreted against one shared expression object --
the same matrix, normalized once, with both label sets attached as metadata columns.
FindAllMarkers runs twice over that single object with identical assay, slot, and
thresholds. The only thing that differs between the two runs is which cells are called
cluster 0, cluster 1, and so on.

That is what makes the evidence independent. Marker expression comes from the RNA, not
from PCA and not from scGPT, so neither representation can score itself.

The cluster numbers mean nothing across arms, and nothing across papers. Louvain
numbers by descending size, so PCA cluster 3 and scGPT cluster 3 are unrelated, and
neither has anything to do with the paper's cluster 3.

This file checks inputs and hands them to R. It performs no analysis.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd


# presto is required, not optional: it supplies the Wilcoxon implementation, and both
# partitions must be tested by the same one for a difference between their marker
# lists to mean biology rather than code.
REQUIRED_R_PACKAGES = ("Seurat", "Matrix", "jsonlite", "presto")

# The frozen label columns, named as Stage 12 named them.
PARTITIONS = ("pca30_harmony_cluster", "scgpt_harmony_cluster")

EXPECTED_CELLS = 70391
EXPECTED_GENES = 19689
DEFAULT_SEED = 42

# The authors' own thresholds: FindAllMarkers(s, min.pct = 0.25, logfc.threshold = 0.25).
DEFAULT_MIN_PCT = 0.25
DEFAULT_LOGFC = 0.25


@dataclass(frozen=True)
class PipelinePaths:
    """All Stage 13 input and output paths below one workspace directory."""

    workspace: Path

    @classmethod
    def from_workspace(cls, workspace: str | Path) -> "PipelinePaths":
        return cls(Path(workspace).expanduser().resolve())

    @property
    def seurat_input_dir(self) -> Path:
        return self.workspace / "results" / "03_seurat_reproduction" / "seurat_input"

    @property
    def frozen_dir(self) -> Path:
        return self.workspace / "results" / "12_frozen_labels"

    @property
    def marker_dir(self) -> Path:
        return self.workspace / "results" / "13_marker_discovery"

    def create_directories(self) -> None:
        self.marker_dir.mkdir(parents=True, exist_ok=True)


def default_workspace() -> Path:
    """Use the repository's Age_trajectory directory unless --workspace is given."""

    return Path(__file__).resolve().parents[1] / "Age_trajectory"


def require_files(paths: Iterable[Path]) -> None:
    """Stop before doing any work if an input this stage depends on is absent."""

    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required inputs:\n  " + "\n  ".join(missing))


# --- Step 1: confirm the frozen labels are intact ----------------------------
# The freeze is the premise of this stage, so it is verified rather than assumed. The
# checksum is recomputed the way Stage 12 computed it; a mismatch means the labels
# moved after being frozen, which would invalidate the comparison silently.
def verify_frozen_labels(paths: PipelinePaths) -> tuple[pd.DataFrame, str]:
    """Check Stage 12's freeze is present, complete, and unchanged since it was written."""

    import hashlib

    labels_file = paths.frozen_dir / "frozen_cluster_labels.csv.gz"
    audit_file = paths.frozen_dir / "frozen_labels_audit.json"
    require_files([labels_file, audit_file])

    frozen = pd.read_csv(labels_file, dtype=str)
    for column in ("cell_barcode", *PARTITIONS):
        if column not in frozen.columns:
            raise ValueError(f"Frozen label file is missing the column {column}")
    if len(frozen) != EXPECTED_CELLS:
        raise ValueError(
            f"Frozen labels cover {len(frozen)} cells, expected {EXPECTED_CELLS}"
        )

    columns = ["cell_barcode", *PARTITIONS]
    payload = frozen[columns].to_csv(index=False).encode("utf-8")
    checksum = hashlib.sha256(payload).hexdigest()

    recorded = json.loads(audit_file.read_text(encoding="utf-8"))
    expected = recorded.get("label_checksum_sha256")
    if expected != checksum:
        raise ValueError(
            "The frozen labels have changed since Stage 12 wrote them.\n"
            f"  Stage 12 recorded: {expected}\n"
            f"  These labels are:  {checksum}\n"
            "The partitions are meant to be frozen; re-run Stage 12 deliberately "
            "if they genuinely need to change."
        )

    print(f"Frozen labels verified against Stage 12 (checksum {checksum[:16]}...).")
    for partition in PARTITIONS:
        sizes = frozen[partition].value_counts()
        print(
            f"  {partition}: {len(sizes)} clusters, "
            f"{sizes.min()}-{sizes.max()} cells per cluster"
        )
    return frozen, checksum


# --- Step 2: confirm the labels and the expression describe the same cells ---
# Both partitions must be judged against the same cells, and those cells must be the
# ones the RNA matrix holds. A mismatch here would not crash anything downstream; it
# would attach expression to the wrong cells and produce a plausible, wrong annotation.
def verify_expression_matches(frozen: pd.DataFrame, paths: PipelinePaths) -> None:
    """Check the RNA matrix holds exactly the cells and genes the labels expect."""

    barcodes_file = paths.seurat_input_dir / "barcodes.tsv"
    genes_file = paths.seurat_input_dir / "genes.tsv"
    counts_file = paths.seurat_input_dir / "counts.mtx"
    require_files([barcodes_file, genes_file, counts_file])

    barcodes = barcodes_file.read_text(encoding="utf-8").split()
    genes = genes_file.read_text(encoding="utf-8").split()

    if len(barcodes) != EXPECTED_CELLS:
        raise ValueError(
            f"The RNA matrix describes {len(barcodes)} cells, expected {EXPECTED_CELLS}"
        )
    if len(genes) != EXPECTED_GENES:
        raise ValueError(
            f"The RNA matrix holds {len(genes)} genes, expected {EXPECTED_GENES}"
        )
    if set(frozen["cell_barcode"]) != set(barcodes):
        raise ValueError(
            "The frozen labels and the RNA matrix do not describe the same cells"
        )

    print(
        f"Expression: {len(genes)} genes x {len(barcodes)} cells, "
        "covering every labelled cell."
    )


def audit_r_environment(rscript: str | None = None) -> dict[str, str]:
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
        raise RuntimeError("Missing required R packages: " + ", ".join(missing))
    print("R packages: " + ", ".join(f"{k} {v}" for k, v in versions.items()))
    return versions


def bundled_r_script() -> Path:
    """The R script that sits beside this one and does the actual work."""

    script = Path(__file__).resolve().with_suffix(".R")
    if not script.exists():
        raise FileNotFoundError(f"Companion R script not found at {script}")
    return script


def run_marker_discovery(
    paths: PipelinePaths,
    seed: int,
    min_pct: float,
    logfc: float,
    rscript: str | None = None,
) -> Path:
    """Start the R stage and confirm it produced a result for both partitions."""

    executable = rscript or shutil.which("Rscript")
    if executable is None:
        raise RuntimeError("Rscript was not found on PATH")

    command = [
        executable,
        "--vanilla",
        str(bundled_r_script()),
        str(paths.seurat_input_dir),
        str(paths.frozen_dir),
        str(paths.marker_dir),
        str(seed),
        str(min_pct),
        str(logfc),
    ]
    print("\nStarting R:\n  " + " ".join(command) + "\n")
    subprocess.run(command, check=True)

    for partition in PARTITIONS:
        result = paths.marker_dir / f"markers_{partition}.csv.gz"
        if not result.exists():
            raise RuntimeError(f"R finished without writing markers for {partition}")
    panel = paths.marker_dir / "canonical_panel_by_cluster.csv"
    if not panel.exists():
        raise RuntimeError("R finished without writing the canonical panel table")
    return panel


# --- Step 3: report what each partition produced -----------------------------
# A count of significant markers per cluster is the first thing worth seeing. A cluster
# with almost none is one the RNA cannot distinguish from its neighbours, whatever the
# geometry said about it.
def report_marker_yield(paths: PipelinePaths) -> dict[str, object]:
    """Summarize how many genes distinguish each cluster, per partition."""

    report: dict[str, object] = {}
    for partition in PARTITIONS:
        markers = pd.read_csv(paths.marker_dir / f"markers_{partition}.csv.gz")
        significant = markers[(markers["p_val_adj"] < 0.05) & (markers["avg_log2FC"] > 0)]
        per_cluster = significant.groupby("cluster").size().sort_index()
        print(f"\n{partition}: {len(significant)} significant positive markers")
        print("  per cluster: " + ", ".join(f"{c}:{n}" for c, n in per_cluster.items()))
        thin = per_cluster[per_cluster < 10]
        if len(thin) > 0:
            print(
                "  clusters the RNA barely separates (<10 markers): "
                + ", ".join(str(c) for c in thin.index)
            )
        report[partition] = {
            "significant_positive_markers": int(len(significant)),
            "per_cluster": {str(k): int(v) for k, v in per_cluster.items()},
        }
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--workspace",
        default=str(default_workspace()),
        help="Directory holding data/ and results/ (default: Age_trajectory)",
    )
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED, help="Random seed")
    parser.add_argument(
        "--min-pct",
        type=float,
        default=DEFAULT_MIN_PCT,
        help=f"FindAllMarkers min.pct (default: {DEFAULT_MIN_PCT}, the authors' value)",
    )
    parser.add_argument(
        "--logfc-threshold",
        type=float,
        default=DEFAULT_LOGFC,
        help=f"FindAllMarkers logfc.threshold (default: {DEFAULT_LOGFC}, the authors' value)",
    )
    parser.add_argument("--rscript", help="Optional explicit path to Rscript")
    parser.add_argument(
        "--verify-only",
        action="store_true",
        help="Run the input checks, then stop before the analysis",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    paths = PipelinePaths.from_workspace(args.workspace)
    paths.create_directories()

    frozen, checksum = verify_frozen_labels(paths)
    verify_expression_matches(frozen, paths)

    if args.verify_only:
        print("\nVerification complete; stopping before the analysis as requested.")
        return

    versions = audit_r_environment(rscript=args.rscript)
    panel = run_marker_discovery(
        paths,
        seed=args.seed,
        min_pct=args.min_pct,
        logfc=args.logfc_threshold,
        rscript=args.rscript,
    )
    yields = report_marker_yield(paths)

    audit = {
        "stage": "13_marker_discovery",
        "seed": args.seed,
        "min_pct": args.min_pct,
        "logfc_threshold": args.logfc_threshold,
        "assay": "RNA",
        "slot": "data",
        "partitions": list(PARTITIONS),
        "frozen_label_checksum_sha256": checksum,
        "r_packages": versions,
        "marker_yield": yields,
    }
    audit_path = paths.marker_dir / "marker_discovery_environment_audit.json"
    audit_path.write_text(json.dumps(audit, indent=2), encoding="utf-8")

    print(f"\nStage 13 complete.\n  Panel table: {panel}\n  Audit: {audit_path}")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:  # noqa: BLE001
        print(f"Stage 13 failed: {error}", file=sys.stderr)
        raise SystemExit(1) from error
