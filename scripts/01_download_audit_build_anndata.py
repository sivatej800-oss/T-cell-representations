#!/usr/bin/env python3
"""Stage 01 command: download GSE200164 and build the canonical AnnData file.

Run this file from the repository root. It contains only the command-line interface;
the data-processing functions live in ``figure1d_pipeline.data_preparation``.
"""

from __future__ import annotations

import argparse

from figure1d_pipeline.config import PipelinePaths, default_workspace


def parse_args() -> argparse.Namespace:
    """Define the options shown by ``python ... --help`` and parse user input."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--workspace",
        default=default_workspace(),
        help="Generated-data workspace (default: checkout-local Age_trajectory)",
    )
    parser.add_argument("--chunk-size", type=int, default=100, help="Gene rows per CSV chunk")
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-download and overwrite existing generated outputs",
    )
    return parser.parse_args()


def main() -> None:
    """Translate command-line options into one call to the Stage 01 implementation."""

    args = parse_args()
    # Import the heavy scientific packages only after parsing arguments. This keeps
    # ``--help`` fast and usable even before the Python dependencies are installed.
    from figure1d_pipeline.data_preparation import run_data_preparation

    # PipelinePaths turns one workspace folder into all input/output paths used below.
    run_data_preparation(
        PipelinePaths.from_workspace(args.workspace),
        chunk_size=args.chunk_size,
        force=args.force,
    )


if __name__ == "__main__":
    # This guard runs main() only when this file is used as a command. Importing the
    # file from a test or another module will not start a multi-gigabyte data build.
    main()
