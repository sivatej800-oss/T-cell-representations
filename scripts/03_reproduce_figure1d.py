#!/usr/bin/env python3
"""Stage 03 command: run or resume the R/Seurat Figure 1D reproduction.

This is the final command users run. Python validates the environment and checkpoints;
the bundled R backend performs the Seurat operations.
"""

from __future__ import annotations

import argparse

from figure1d_pipeline.config import PipelinePaths, default_workspace


def parse_args() -> argparse.Namespace:
    """Define the Stage 03 command-line options and parse user input."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--workspace",
        default=default_workspace(),
        help="Generated-data workspace (default: checkout-local Age_trajectory)",
    )
    parser.add_argument("--rscript", help="Optional explicit path to Rscript")
    parser.add_argument("--seed", type=int, default=1234, help="Random seed")
    parser.add_argument(
        "--force",
        action="store_true",
        help="Recompute existing sample checkpoints and the final result",
    )
    return parser.parse_args()


def main() -> None:
    """Run the Seurat backend and report the canonical final RDS path."""

    args = parse_args()
    # The import is delayed so users can inspect ``--help`` without loading R helpers.
    from figure1d_pipeline.seurat_runner import run_seurat_figure1d

    result = run_seurat_figure1d(
        PipelinePaths.from_workspace(args.workspace),
        rscript=args.rscript,
        seed=args.seed,
        force=args.force,
    )
    print(f"Integrated Seurat result: {result}")


if __name__ == "__main__":
    # Importing this file is safe; only direct command-line execution starts Stage 03.
    main()
