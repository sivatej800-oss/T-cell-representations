#!/usr/bin/env python3
"""Stage 14: make Stage 13's marker results legible, without interpreting them.

Stage 13 answered "which genes distinguish each cluster" and wrote the answer as two
large tables. This stage turns those into what a person reads: ranked marker lists per
cluster, dot plots of the paper's canonical panel, top-marker heatmaps, and one
evidence row per cluster carrying its size, donor spread, markers, and whether markers
of more than one lineage are present.

No identities are assigned. A cluster expressing CD3D is reported as expressing CD3D.
The separation matters because the evidence should be legible before anyone decides
what it means, and because the same evidence has to be visible for both partitions on
the same scale -- a dot plot scaled per panel makes every partition look equally well
separated, since the strongest cluster in each is always the darkest dot.

Nothing is re-clustered and nothing is re-normalized. The expression object is rebuilt
exactly as Stage 13 built it and checked against the fingerprint Stage 13 recorded, so
the two stages are known to be describing one matrix rather than assumed to be.

This file checks inputs and hands them to R. It performs no analysis.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd


REQUIRED_R_PACKAGES = ("Seurat", "Matrix", "ggplot2", "patchwork", "jsonlite")

PARTITIONS = ("pca30_harmony_cluster", "scgpt_harmony_cluster")

# The paper's panel, in the paper's order. CD44 is included here and was not in the
# 13-gene set Stage 13 scored, so Stage 14 recomputes the panel rather than reading
# Stage 13's table back.
MARKER_PANEL = (
    "CD3D", "CD4", "CD8A", "CD8B", "FOXP3",
    "NKG7", "GNLY", "FCER1A", "CD68",
    "CD14", "FCGR3A", "CD19", "JCHAIN", "CD44",
)

EXPECTED_CELLS = 70391
EXPECTED_GENES = 19689


@dataclass(frozen=True)
class PipelinePaths:
    workspace: Path

    @classmethod
    def from_workspace(cls, workspace: str | Path) -> "PipelinePaths":
        return cls(Path(workspace).expanduser().resolve())

    @property
    def seurat_input_dir(self) -> Path:
        return self.workspace / "results" / "03_seurat_reproduction" / "seurat_input"

    @property
    def benchmark_dir(self) -> Path:
        return self.workspace / "results" / "06_harmony_benchmark"

    @property
    def frozen_dir(self) -> Path:
        return self.workspace / "results" / "12_frozen_labels"

    @property
    def marker_dir(self) -> Path:
        return self.workspace / "results" / "13_marker_discovery"

    @property
    def evidence_dir(self) -> Path:
        return self.workspace / "results" / "14_marker_evidence"

    @property
    def figure_dir(self) -> Path:
        return self.workspace / "results" / "figures"

    def create_directories(self) -> None:
        self.evidence_dir.mkdir(parents=True, exist_ok=True)
        self.figure_dir.mkdir(parents=True, exist_ok=True)


def default_workspace() -> Path:
    return Path(__file__).resolve().parents[1] / "Age_trajectory"


def require_files(paths: Iterable[Path]) -> None:
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required inputs:\n  " + "\n  ".join(missing))


# --- Step 1: confirm Stage 13 finished and both arms ran the same way --------
# The two marker files are only comparable if they came from one run over one object.
# Stage 13's audit records the parameters and the expression fingerprint; both are
# checked here rather than assumed, and the fingerprint is checked again inside R
# against the object this stage rebuilds.
def verify_stage13(paths: PipelinePaths) -> dict:
    """Check Stage 13's output is present, complete, and internally consistent."""

    audit_file = paths.marker_dir / "marker_discovery_audit.json"
    marker_files = [
        paths.marker_dir / f"markers_{partition}.csv.gz" for partition in PARTITIONS
    ]
    require_files([audit_file, *marker_files])

    audit = json.loads(audit_file.read_text(encoding="utf-8"))

    if audit.get("assay") != "RNA" or audit.get("slot") != "data":
        raise ValueError(
            f"Stage 13 did not run on the RNA data layer "
            f"(assay={audit.get('assay')}, slot={audit.get('slot')})"
        )
    if sorted(audit.get("partitions", [])) != sorted(PARTITIONS):
        raise ValueError("Stage 13 did not test the two partitions this stage expects")

    fingerprint = audit.get("expression_fingerprint", {})
    if fingerprint.get("n_cells") != EXPECTED_CELLS or fingerprint.get("n_genes") != EXPECTED_GENES:
        raise ValueError("Stage 13's expression was not the full matrix")

    print(
        f"Stage 13: assay {audit['assay']}/{audit['slot']}, "
        f"min.pct {audit['min_pct']}, logfc {audit['logfc_threshold']}, "
        f"{audit.get('wilcox_implementation', 'wilcox')}"
    )
    print(f"  expression fingerprint: {fingerprint['data_sum']}")
    return audit


# --- Step 2: report what the marker tables contain ---------------------------
# Positive markers are what describes a cluster; the files also hold negative ones,
# because Stage 13 ran with only.pos = FALSE. Saying so here keeps the counts in this
# stage's tables from looking inconsistent with the raw files.
def report_marker_tables(paths: PipelinePaths) -> dict[str, object]:
    """Summarize each partition's marker table before it is ranked."""

    report: dict[str, object] = {}
    for partition in PARTITIONS:
        markers = pd.read_csv(paths.marker_dir / f"markers_{partition}.csv.gz")
        positive = markers[(markers["p_val_adj"] < 0.05) & (markers["avg_log2FC"] > 0)]
        per_cluster = positive.groupby("cluster").size()
        print(
            f"\n{partition}: {len(markers)} rows, "
            f"{len(positive)} significant positive across {markers['cluster'].nunique()} clusters"
        )
        print(
            f"  per cluster: {per_cluster.min()}-{per_cluster.max()} "
            f"(median {int(per_cluster.median())})"
        )
        thin = per_cluster[per_cluster < 20]
        if len(thin) > 0:
            print(
                "  clusters with fewer than 20 positive markers: "
                + ", ".join(f"{c} ({n})" for c, n in thin.items())
                + " — their top-20 lists will be short"
            )
        report[partition] = {
            "rows": int(len(markers)),
            "significant_positive": int(len(positive)),
            "min_per_cluster": int(per_cluster.min()),
            "max_per_cluster": int(per_cluster.max()),
        }
    return report


def verify_layouts_available(paths: PipelinePaths) -> None:
    """FeaturePlot needs a layout; Stage 06 computed one per corrected arm."""

    layouts = paths.benchmark_dir / "benchmark_cluster_assignments.csv.gz"
    require_files([layouts])
    header = pd.read_csv(layouts, nrows=0).columns
    needed = [
        f"{prefix}_UMAP_{axis}"
        for prefix in ("pca_harmony", "scgpt_harmony")
        for axis in (1, 2)
    ]
    missing = [column for column in needed if column not in header]
    if missing:
        raise ValueError("Stage 06 output is missing UMAP columns: " + ", ".join(missing))
    print("Layouts: both Harmony UMAPs present in Stage 06's output.")


def verify_panel_available(paths: PipelinePaths) -> None:
    """Check every panel gene is in the matrix before R spends time loading it."""

    genes_file = paths.seurat_input_dir / "genes.tsv"
    require_files([genes_file])
    genes = set(genes_file.read_text(encoding="utf-8").split())
    missing = [gene for gene in MARKER_PANEL if gene not in genes]
    if missing:
        raise ValueError("Panel genes absent from the matrix: " + ", ".join(missing))
    print(f"\nPanel: all {len(MARKER_PANEL)} genes present in the matrix.")


def audit_r_environment(rscript: str | None = None) -> dict[str, str]:
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
    script = Path(__file__).resolve().with_suffix(".R")
    if not script.exists():
        raise FileNotFoundError(f"Companion R script not found at {script}")
    return script


def run_evidence(paths: PipelinePaths, rscript: str | None = None) -> None:
    executable = rscript or shutil.which("Rscript")
    if executable is None:
        raise RuntimeError("Rscript was not found on PATH")

    command = [
        executable,
        "--vanilla",
        str(bundled_r_script()),
        str(paths.seurat_input_dir),
        str(paths.frozen_dir),
        str(paths.benchmark_dir),
        str(paths.marker_dir),
        str(paths.evidence_dir),
        str(paths.figure_dir),
    ]
    print("\nStarting R:\n  " + " ".join(command) + "\n")
    # A threaded BLAS plus fork() is what broke the first attempt at this stage: the
    # children inherit lock state they cannot resolve and never return. The table
    # steps still use mclapply, so the BLAS is held to one thread for the whole run.
    environment = {
        **os.environ,
        "OPENBLAS_NUM_THREADS": "1",
        "OMP_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
    }
    subprocess.run(command, check=True, env=environment)

    expected = [
        paths.evidence_dir / "cluster_evidence_table.csv",
        paths.evidence_dir / "canonical_panel_by_cluster.csv",
        *[paths.evidence_dir / f"top_markers_{p}.csv" for p in PARTITIONS],
        *[paths.figure_dir / f"14_featureplots_{p}.png" for p in PARTITIONS],
        *[paths.figure_dir / f"14_clusters_{p}.png" for p in PARTITIONS],
    ]
    missing = [str(path) for path in expected if not path.exists()]
    if missing:
        raise RuntimeError("R finished without writing:\n  " + "\n  ".join(missing))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", default=str(default_workspace()))
    parser.add_argument("--rscript", help="Optional explicit path to Rscript")
    parser.add_argument(
        "--verify-only",
        action="store_true",
        help="Run the input checks and the marker-table report, then stop",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    paths = PipelinePaths.from_workspace(args.workspace)
    paths.create_directories()

    audit = verify_stage13(paths)
    tables = report_marker_tables(paths)
    verify_panel_available(paths)
    verify_layouts_available(paths)

    if args.verify_only:
        print("\nVerification complete; stopping before the analysis as requested.")
        return

    versions = audit_r_environment(rscript=args.rscript)
    run_evidence(paths, rscript=args.rscript)

    evidence = pd.read_csv(paths.evidence_dir / "cluster_evidence_table.csv")
    print("\n=== Evidence table ===")
    for partition in PARTITIONS:
        block = evidence[evidence["partition"] == partition]
        print(f"\n{partition}:")
        print(
            block[
                ["cluster", "n_cells", "n_donors", "largest_donor_share",
                 "canonical_present"]
            ].to_string(index=False)
        )

    audit_path = paths.evidence_dir / "marker_evidence_environment_audit.json"
    audit_path.write_text(
        json.dumps(
            {
                "stage": "14_marker_evidence",
                "panel": list(MARKER_PANEL),
                "partitions": list(PARTITIONS),
                "stage13_fingerprint": audit["expression_fingerprint"]["data_sum"],
                "marker_tables": tables,
                "r_packages": versions,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    print(f"\nStage 14 complete.\n  Tables: {paths.evidence_dir}\n  Figures: {paths.figure_dir}")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:  # noqa: BLE001
        print(f"Stage 14 failed: {error}", file=sys.stderr)
        raise SystemExit(1) from error
