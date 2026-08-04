#!/usr/bin/env python3
"""Stage 02 command: audit AnnData, calculate QC, and prepare Seurat inputs.

Run this file after Stage 01. The implementation is kept separately in
``figure1d_pipeline.figure1d_preparation`` so the command remains easy to read.
"""

from __future__ import annotations

import argparse

from figure1d_pipeline.config import PipelinePaths, default_workspace


def parse_args() -> argparse.Namespace:
    """Define the Stage 02 command-line options and parse user input."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--workspace",
        default=default_workspace(),
        help="Generated-data workspace (default: checkout-local Age_trajectory)",
    )
    parser.add_argument("--rscript", help="Optional explicit path to Rscript")
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite the Figure 1D checkpoint and Seurat export",
    )
    return parser.parse_args()


def main() -> None:
    """Run Stage 02 using the workspace and R executable selected by the user."""

    args = parse_args()
    # Delay scientific-library imports so ``--help`` works in a fresh environment.
    from figure1d_pipeline.figure1d_preparation import prepare_figure1d

    prepare_figure1d(
        PipelinePaths.from_workspace(args.workspace),
        force=args.force,
        rscript=args.rscript,
    )


if __name__ == "__main__":
    # Prevent accidental execution when another Python file imports this module.
    main()
