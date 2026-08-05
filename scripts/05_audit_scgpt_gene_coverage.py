#!/usr/bin/env python3
"""Stage 05: measure what scGPT's fixed vocabulary excluded from this dataset.

scGPT does not read expression as a vector in the dataset's column order. It reads
genes as tokens looked up in a vocabulary frozen at pretraining, and a symbol with
no entry is dropped before the model sees any cell. For this dataset that is 3,729
of 19,689 genes — 18.9%, which sounds severe.

It is not, and this stage exists to show why rather than to assert it. Counting
genes weights a clone-based locus that is never expressed the same as a histone
transcribed in every cell. Summing the expression behind each side of the split
instead gives the number that matters: how much of the measured signal the model
actually received.

Two questions are answered here, and they are different questions:

  1. Which genes were dropped, and is any of them a marker used to identify a cell
     population? A lost lineage marker would be serious, because this pipeline has
     already had rare populations destroyed once by feature selection.
  2. How much expression did the dropped genes carry? This is the size of the loss,
     as opposed to its count.

Both were open when the Stage 04 comparison was first run, and both are settled by
running this. It reads only finished artifacts and writes only a report, so it can
be run at any time and changes nothing.
"""

from __future__ import annotations

import argparse
import json
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np


# --- The vocabulary this run used --------------------------------------------
# The model directory's own vocab.json is not in this repository; it lived with the
# GPU run. The vocabulary vendored under archive/ is not byte-identical to it — the
# model file prepends special tokens, so the checksums differ — but its gene content
# is what matters, and that is checkable: its 60,694 tokens reproduce the run's
# recorded 15,960 / 3,729 split exactly. A different vocabulary would not.
VENDORED_VOCABULARY = Path("archive/vendor/scGPT-main/scgpt/tokenizer/default_census_vocab.json")

# Symbols with no stable name. These are largely absent from curated references, so
# their absence from the vocabulary costs little, and separating them from ordinary
# gene symbols is what turns a list of 3,729 names into an interpretable result.
UNSTABLE_NAME_PATTERNS = {
    "clone-based locus (AC/AL/AP…)": r"^(AC|AL|AP|BX|CR|CT|FP|Z\d)\d+\.\d+$",
    "pseudogene (…P<n>)": r"P\d+$",
    "LINC lncRNA": r"^LINC\d+",
    "antisense (…-AS<n>)": r"-AS\d+$",
    "divergent transcript (…-DT)": r"-DT$",
    "uncharacterized ORF": r"orf\d+",
    "small RNA (MIR/SNOR/RNU)": r"^(MIR|SNOR|RNU|RNVU|SCARNA)",
    "novel transcript (….<n>)": r"\.\d+$",
}

# The markers used to identify populations in CSF: T-cell subsets, B and plasma
# cells, monocytes and macrophages, microglia, dendritic cells. If scGPT's
# vocabulary is missing any of these, a cluster difference between the two arms
# could be the missing marker rather than the representation, which is precisely
# the ambiguity this stage exists to rule out.
CANONICAL_MARKERS = (
    "CD3E", "CD8A", "CD4", "IL7R", "CCR7", "SELL", "GZMB", "GZMK", "NKG7", "GNLY",
    "KLRD1", "KLRB1", "TRAC", "TRBC1", "TRBC2", "FOXP3", "CTLA4", "PDCD1", "TCF7",
    "LEF1", "CCR6", "KLRG1", "TOX", "CD27", "CD38",
    "MS4A1", "CD79A", "CD79B", "CD19", "BANK1", "JCHAIN", "MZB1", "DERL3", "IGHM",
    "IGHG1", "XBP1", "PRDM1", "SDC1",
    "LYZ", "CD14", "FCGR3A", "ITGAX", "CST3", "AIF1", "C1QA", "C1QB", "APOE",
    "TREM2", "SPP1", "S100A8", "S100A9", "VCAN",
    "P2RY12", "TMEM119", "CX3CR1",
    "HLA-DRA", "HLA-DRB1", "CD74", "CLEC9A", "FCER1A", "IRF8",
    "PPBP", "HBB", "PF4", "NRGN",
)


# --- Workspace paths ---------------------------------------------------------
# The numbered stages repeat these helpers rather than importing them, so each can
# be read and run on its own. This stage follows that decision.
@dataclass(frozen=True)
class PipelinePaths:
    """All Stage 05 input and output paths below one workspace directory."""

    workspace: Path

    @classmethod
    def from_workspace(cls, workspace: str | Path) -> "PipelinePaths":
        return cls(Path(workspace).expanduser().resolve())

    # The matrix scGPT was actually given. Using this rather than any other copy is
    # what guarantees the gene column order matches the lookup the model performed.
    @property
    def scgpt_input(self) -> Path:
        return self.workspace / "data" / "processed" / "GSE200164_scgpt_input.h5ad"

    @property
    def run_metadata(self) -> Path:
        return self.workspace / "data" / "scgpt" / "scgpt_run_metadata.json"

    @property
    def scgpt_result_dir(self) -> Path:
        return self.workspace / "results" / "figure1d_scgpt"


def default_workspace() -> Path:
    """Use the repository's Age_trajectory directory unless --workspace is given."""

    return Path(__file__).resolve().parents[1] / "Age_trajectory"


def repository_root() -> Path:
    """The checkout this script lives in, which is where the vocabulary is vendored."""

    return Path(__file__).resolve().parents[1]


# --- Reading an AnnData column without loading the matrix --------------------
def read_column(group: h5py.Group, key: str) -> np.ndarray:
    """Read one obs/var column, which may be plain, categorical, or nullable."""

    node = group[key]
    if isinstance(node, h5py.Group):
        if "categories" in node:
            values = node["categories"][:][node["codes"][:]]
        else:
            values = node["values"][:]
    else:
        values = node[:]
    return np.array([v.decode() if isinstance(v, bytes) else v for v in values])


# --- Step 1: split the dataset's genes by whether scGPT could look them up ----
def split_genes(symbols: np.ndarray, vocabulary: set[str], run_metadata: dict) -> np.ndarray:
    """Return a mask of genes present in the vocabulary, checked against the run."""

    kept_mask = np.array([symbol in vocabulary for symbol in symbols])
    kept = int(kept_mask.sum())
    dropped = int((~kept_mask).sum())
    # The run recorded how many genes it matched. Reproducing that number exactly is
    # what identifies this vocabulary as the one the run used, despite the model's
    # own file carrying a different checksum. A mismatch means the two disagree and
    # nothing below can be trusted.
    if kept != run_metadata["matched_genes"]:
        raise ValueError(
            f"This vocabulary matches {kept} genes but the run recorded "
            f"{run_metadata['matched_genes']}; it is not the vocabulary that was used"
        )
    if kept + dropped != run_metadata["total_dataset_genes"]:
        raise ValueError("Gene count disagrees with the recorded run metadata")
    print(f"genes {kept + dropped}  kept {kept}  dropped {dropped}")
    return kept_mask


# --- Step 2: what kind of names were dropped ---------------------------------
# A dropped symbol is either a locus with no stable name, which costs little, or an
# ordinary gene symbol, which is the case worth looking at. Nearly all of the latter
# turn out to be nomenclature drift: the dataset's names come from a GRCh38
# CellRanger reference predating HGNC renames the vocabulary already reflects, so
# the gene is known to scGPT under a newer symbol.
def classify_dropped_names(dropped: list[str]) -> tuple[Counter, list[str]]:
    """Group dropped symbols by name shape, separating out conventional symbols."""

    category: Counter = Counter()
    conventional: list[str] = []
    for symbol in dropped:
        for label, pattern in UNSTABLE_NAME_PATTERNS.items():
            if re.search(pattern, symbol):
                category[label] += 1
                break
        else:
            category["conventional HGNC-style symbol"] += 1
            conventional.append(symbol)
    print("\ndropped genes by name class:")
    for label, count in category.most_common():
        print(f"  {count:6d}  {label}")
    print(f"\nconventional symbols dropped: {len(conventional)}")
    return category, sorted(conventional)


# --- Step 3: how much expression the dropped genes carried -------------------
# The measurement that reframes the whole question. The matrix holds 126.6 million
# recorded values, so it is summed in chunks rather than held in memory at once.
# Values are the SoupX-corrected fractional expression this pipeline preserves
# deliberately; nothing is rounded.
def expression_by_gene(path: Path, n_genes: int) -> np.ndarray:
    """Total expression per gene across every cell, read straight from the CSR data."""

    per_gene = np.zeros(n_genes, dtype=np.float64)
    with h5py.File(path, "r") as handle:
        indices = handle["X"]["indices"]
        data = handle["X"]["data"]
        n_entries = data.shape[0]
        chunk = 8_000_000
        for start in range(0, n_entries, chunk):
            stop = min(start + chunk, n_entries)
            np.add.at(per_gene, indices[start:stop], data[start:stop].astype(np.float64))
    return per_gene


# --- Stage 05 controller -----------------------------------------------------
def audit_gene_coverage(paths: PipelinePaths) -> dict[str, object]:
    """Answer both questions and write the report beside the Stage 04 results."""

    for required in (paths.scgpt_input, paths.run_metadata):
        if not required.is_file():
            raise FileNotFoundError(f"Required file was not found: {required}")
    vocabulary_path = repository_root() / VENDORED_VOCABULARY
    if not vocabulary_path.is_file():
        raise FileNotFoundError(f"Vendored scGPT vocabulary was not found: {vocabulary_path}")

    run_metadata = json.loads(paths.run_metadata.read_text(encoding="utf-8"))
    vocabulary = set(json.load(vocabulary_path.open(encoding="utf-8")))
    print(f"vocabulary: {len(vocabulary)} tokens from {VENDORED_VOCABULARY}")

    # The gene symbols in the order scGPT received them. The run recorded which
    # column it matched on, so that column is used rather than the frame's index.
    with h5py.File(paths.scgpt_input, "r") as handle:
        symbols = read_column(handle["var"], run_metadata["input_gene_column"])
        n_cells, n_genes = (int(value) for value in handle["X"].attrs["shape"])
    if [n_cells, n_genes] != run_metadata["input_shape"]:
        raise ValueError("Input matrix no longer matches the recorded run metadata")

    kept_mask = split_genes(symbols, vocabulary, run_metadata)
    dropped_symbols = [str(s) for s in symbols[~kept_mask]]
    category, conventional = classify_dropped_names(dropped_symbols)

    # The size of the loss, as opposed to its count.
    per_gene = expression_by_gene(paths.scgpt_input, n_genes)
    total = per_gene.sum()
    kept_total = per_gene[kept_mask].sum()
    dropped_total = per_gene[~kept_mask].sum()
    print(f"\ntotal expression across all cells and genes: {total:,.0f}")
    print(f"  kept    {kept_total:>18,.0f}   {kept_total / total:7.3%}")
    print(f"  dropped {dropped_total:>18,.0f}   {dropped_total / total:7.3%}")

    # Whether the loss is diffuse or concentrated changes how it should be read. A
    # handful of abundant genes carrying most of it is a narrow, nameable gap; the
    # same total spread evenly across 3,729 genes would not be.
    dropped_index = np.flatnonzero(~kept_mask)
    ranked = dropped_index[np.argsort(per_gene[dropped_index])[::-1]]
    print("\nhighest-expression dropped genes:")
    for index in ranked[:20]:
        print(f"  {symbols[index]:<18} {per_gene[index]:>14,.0f}   {per_gene[index] / total:8.4%}")
    top_20_share = float(per_gene[ranked[:20]].sum() / dropped_total)
    print(f"\ntop 20 dropped genes account for {top_20_share:.1%} of the dropped signal")

    # The question the whole stage is for: did a marker used to identify a cell
    # population fall out of the vocabulary?
    in_dataset = set(str(s) for s in symbols)
    dropped_set = set(dropped_symbols)
    absent_from_dataset = [m for m in CANONICAL_MARKERS if m not in in_dataset]
    dropped_markers = [m for m in CANONICAL_MARKERS if m in in_dataset and m in dropped_set]
    print(f"\ncanonical markers checked: {len(CANONICAL_MARKERS)}")
    print(f"  absent from the dataset itself: {absent_from_dataset or 'none'}")
    print(f"  present but dropped by scGPT:   {dropped_markers or 'none'}")

    report = {
        "vocabulary_file": str(VENDORED_VOCABULARY),
        "vocabulary_tokens": len(vocabulary),
        "vocabulary_identified_by": (
            "reproduces the run's recorded matched-gene count exactly; the model's own "
            "vocab.json is not in this repository and carries a different checksum "
            "because it prepends special tokens"
        ),
        "genes_in_dataset": int(n_genes),
        "genes_kept": int(kept_mask.sum()),
        "genes_dropped": int((~kept_mask).sum()),
        "gene_coverage_fraction": float(kept_mask.sum() / n_genes),
        "total_expression": float(total),
        "expression_share_kept": float(kept_total / total),
        "expression_share_dropped": float(dropped_total / total),
        "dropped_name_classes": dict(category),
        "dropped_conventional_symbols": conventional,
        "highest_expression_dropped_genes": [
            {
                "gene": str(symbols[index]),
                "expression": float(per_gene[index]),
                "share_of_total": float(per_gene[index] / total),
            }
            for index in ranked[:20]
        ],
        "top_20_share_of_dropped_signal": top_20_share,
        "dropped_genes_with_zero_expression": int((per_gene[dropped_index] == 0).sum()),
        "canonical_markers_checked": len(CANONICAL_MARKERS),
        "canonical_markers_absent_from_dataset": absent_from_dataset,
        "canonical_markers_dropped_by_scgpt": dropped_markers,
    }
    paths.scgpt_result_dir.mkdir(parents=True, exist_ok=True)
    output = paths.scgpt_result_dir / "scgpt_gene_coverage_audit.json"
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"\nReport written to {output}")
    return report


def parse_args() -> argparse.Namespace:
    """Read Stage 05 command-line options."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--workspace",
        default=default_workspace(),
        help="Generated-data workspace (default: checkout-local Age_trajectory)",
    )
    return parser.parse_args()


def main() -> None:
    """Run the complete fifth stage from this one Python file."""

    # Example: python scripts/05_audit_scgpt_gene_coverage.py
    args = parse_args()
    audit_gene_coverage(PipelinePaths.from_workspace(args.workspace))


if __name__ == "__main__":
    main()
