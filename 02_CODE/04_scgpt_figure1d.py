#!/usr/bin/env python3
"""Stage 04: cluster the same cells from a frozen scGPT embedding instead of PCA.

This is the second arm of the comparison the README describes. Stage 03 describes
each cell by 15 principal components distilled from SCTransform-normalized,
batch-integrated expression. This stage describes the same 70,391 cells by the
512-dimensional vector a pretrained scGPT assigned to each one, and then runs the
identical downstream analysis: nearest neighbors, Louvain at resolution 0.3, UMAP
at n.neighbors 30.

Only the representation differs. Everything after it is the same Seurat code with
the same settings and the same seed, which is the point — a difference in the
result can then be attributed to the representation rather than to the tooling.

The embedding was produced elsewhere, on a GPU, and arrives as a finished artifact.
This file's job is to prove that artifact describes the cells it claims to describe
before any conclusion is drawn from it, then hand it to R. It performs no analysis.

The verification matters more here than it would for a normal input. The embedding
is 70,391 anonymous rows of 512 numbers; nothing about it is legible to inspection.
If its rows were ordered differently from the cells, every downstream step would
still run and would still produce a plausible figure, and it would be wrong. The
run recorded a hash of the cell order it embedded, and that hash is what makes the
correspondence checkable rather than assumed.
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
import psutil


# --- Required R software -----------------------------------------------------
# Fewer packages than Stage 03: this arm skips SCTransform and integration
# entirely, so glmGamPoi is not needed. The representation already exists.
REQUIRED_R_PACKAGES = ("Seurat", "Matrix", "ggplot2")

# The embedding is unit-normalized by scGPT, so no dimension dominates by scale and
# all 512 carry equal weight. See the note on dimensionality in run_scgpt_clustering.
EXPECTED_EMBEDDING_DIMS = 512
EXPECTED_CELLS = 70391


# --- Workspace paths ---------------------------------------------------------
# Stage 03 repeats these helpers rather than importing them, so that each stage can
# be read and run on its own. This stage follows that decision.
@dataclass(frozen=True)
class PipelinePaths:
    """All Stage 04 input and output paths below one workspace directory."""

    workspace: Path

    @classmethod
    def from_workspace(cls, workspace: str | Path) -> "PipelinePaths":
        return cls(Path(workspace).expanduser().resolve())

    # The unpacked scGPT bundle: the embedding, its manifest, and its run metadata.
    @property
    def scgpt_dir(self) -> Path:
        return self.workspace / "data" / "scgpt"

    # Stage 03's output directory, read here only to confirm the two arms describe
    # the same cells in the same order. Nothing in it is written or modified.
    @property
    def pca_arm_dir(self) -> Path:
        return self.workspace / "results" / "03_seurat_reproduction"

    # This arm's results, kept separate from Stage 03's so the two never collide.
    @property
    def scgpt_result_dir(self) -> Path:
        return self.workspace / "results" / "04_scgpt_embedding"

    # What R reads: the embedding as a flat binary block, and the per-cell details.
    @property
    def r_input_dir(self) -> Path:
        return self.scgpt_result_dir / "r_input"

    def create_directories(self) -> None:
        """Create this stage's folders, leaving any previous result in place."""

        for directory in (self.scgpt_result_dir, self.r_input_dir):
            directory.mkdir(parents=True, exist_ok=True)


def default_workspace() -> Path:
    """Use the repository's Age_trajectory directory unless --workspace is given."""

    return Path(__file__).resolve().parents[1] / "Age_trajectory"


# --- Small validation helper ------------------------------------------------
def require_files(paths: Iterable[Path]) -> None:
    """Stop with a readable list when an expected input file is missing."""

    missing = [path for path in paths if not path.is_file()]
    if missing:
        formatted = "\n".join(f"  - {path}" for path in missing)
        raise FileNotFoundError(f"Required files were not found:\n{formatted}")


# --- Step 1: confirm the embedding is the file the run says it is ------------
# The bundle ships checksums for its own contents. Checking them first separates a
# transfer problem from an analysis problem, which are otherwise hard to tell apart
# once the numbers are loaded.
def verify_bundle_checksums(paths: PipelinePaths) -> None:
    """Recompute the bundle's shipped SHA-256 sums and insist they still match."""

    checksum_file = paths.scgpt_dir / "SHA256SUMS.txt"
    require_files([checksum_file])
    for line in checksum_file.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        expected, name = line.split()
        target = paths.scgpt_dir / name
        require_files([target])
        digest = hashlib.sha256()
        with target.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        if digest.hexdigest() != expected:
            raise ValueError(f"Checksum mismatch for {name}; the bundle is corrupt")
    print("scGPT bundle checksums: all match")


# --- Step 2: confirm the embedding describes this pipeline's cells -----------
# The load-bearing check of this stage. Three separate things have to agree: the
# embedding's row count, the manifest that names each row, and the cells Stage 03
# ran on. The scGPT run also recorded a hash of the exact cell order it was given,
# which is reproduced here — that is what rules out a silent reordering, the one
# failure this stage could not otherwise detect.
def verify_cell_correspondence(
    paths: PipelinePaths,
) -> tuple[np.ndarray, pd.DataFrame, dict[str, object]]:
    """Load the embedding and prove its rows are this pipeline's cells, in order."""

    metadata_file = paths.scgpt_dir / "scgpt_run_metadata.json"
    embedding_file = paths.scgpt_dir / "X_scgpt.npy"
    manifest_file = paths.scgpt_dir / "scgpt_cell_manifest.csv.gz"
    require_files([metadata_file, embedding_file, manifest_file])
    run_metadata = json.loads(metadata_file.read_text(encoding="utf-8"))

    # The embedding itself. A non-finite value would propagate silently through the
    # neighbor graph, and a zero-length row would have no direction at all.
    embedding = np.load(embedding_file)
    if list(embedding.shape) != run_metadata["output_shape"]:
        raise ValueError("Embedding shape disagrees with the recorded run metadata")
    if embedding.shape != (EXPECTED_CELLS, EXPECTED_EMBEDDING_DIMS):
        raise ValueError(f"Expected a {EXPECTED_CELLS} x {EXPECTED_EMBEDDING_DIMS} embedding")
    if not np.isfinite(embedding).all():
        raise ValueError("Embedding contains NaN or infinite values")
    norms = np.linalg.norm(embedding, axis=1)
    if norms.min() <= 0:
        raise ValueError("At least one cell has an all-zero embedding")

    # The manifest names the cell behind each row and must be in embedding order.
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

    # The cells Stage 03 ran on. Identical order means the two arms' results can be
    # compared row by row, with no join and no opportunity to misalign them.
    pca_barcodes_file = paths.pca_arm_dir / "seurat_input" / "barcodes.tsv"
    pca_metadata_file = paths.pca_arm_dir / "seurat_input" / "metadata.csv"
    require_files([pca_barcodes_file, pca_metadata_file])
    pca_barcodes = pca_barcodes_file.read_text(encoding="utf-8").split()
    if pca_barcodes != list(cell_ids):
        raise ValueError("Embedded cells do not match Stage 03's cells in the same order")

    # The per-cell details travelled with the embedding, so they are checked against
    # the copy Stage 03 used rather than trusted. A patient mislabelled here would
    # corrupt every donor-mixing comparison the two arms are built to support.
    pca_metadata = pd.read_csv(pca_metadata_file).set_index("cell_barcode").loc[cell_ids]
    for column in ("internal_sample_id", "patient_id", "diagnosis", "sex", "age"):
        if not (pca_metadata[column].to_numpy() == manifest[column].to_numpy()).all():
            raise ValueError(f"Manifest disagrees with Stage 03 metadata on {column}")
    if manifest["internal_sample_id"].nunique() != 59:
        raise ValueError("Expected 59 internal sample IDs")

    print(
        f"Cell correspondence verified: {embedding.shape[0]} cells, "
        f"{embedding.shape[1]} dimensions, order hash reproduced, "
        f"metadata agrees with Stage 03"
    )
    return embedding, manifest, run_metadata


# --- Step 3: hand the embedding to R -----------------------------------------
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

    # The same columns Stage 03's metadata.csv carries, under the same names, so the
    # two arms' output tables line up without translation.
    metadata = manifest.rename(columns={"cell_id": "cell_barcode"})
    metadata.to_csv(paths.r_input_dir / "metadata.csv", index=False)
    print(f"Wrote R inputs to {paths.r_input_dir}")


# --- R and computer-resource preflight ---------------------------------------
# The same check Stage 03 performs, repeated here so this stage can be run on its
# own. Recording package versions matters for the comparison specifically: both
# arms must run the same Seurat, or a difference between them could come from the
# software rather than from the representation.
def audit_r_environment(
    paths: PipelinePaths,
    run_metadata: dict[str, object],
    *,
    rscript: str | None = None,
) -> dict[str, object]:
    """Confirm R and its required packages are available, and record the setting."""

    executable = rscript or shutil.which("Rscript")
    if executable is None:
        raise RuntimeError("Rscript was not found on PATH")
    package_vector = ", ".join(json.dumps(package) for package in REQUIRED_R_PACKAGES)
    code = f"""
packages <- c({package_vector})
for (pkg in packages) {{
  installed <- requireNamespace(pkg, quietly = TRUE)
  version <- if (installed) as.character(packageVersion(pkg)) else NA_character_
  cat(pkg, installed, version, sep = "\\t")
  cat("\\n")
}}
"""
    result = subprocess.run(
        [executable, "-e", code], check=True, capture_output=True, text=True
    )
    packages: dict[str, dict[str, object]] = {}
    for line in result.stdout.splitlines():
        name, installed, version = line.split("\t")
        packages[name] = {
            "installed": installed == "TRUE",
            "version": None if version == "NA" else version,
        }
    missing = [name for name in REQUIRED_R_PACKAGES if not packages[name]["installed"]]
    memory = psutil.virtual_memory()
    disk = shutil.disk_usage(paths.scgpt_result_dir)

    # An open question recorded with the result rather than left in conversation.
    # scGPT reads genes as vocabulary tokens, and 3,729 of this dataset's 19,689
    # gene symbols had no token, so they were dropped before the model saw any cell.
    # That sounds severe and is not, because counting genes weights a never-expressed
    # clone-based locus the same as a histone transcribed in every cell. Summing the
    # expression behind each side of the split instead, the dropped genes account for
    # 1.138% of all measured signal: 18.9% of the genes carry 1.1% of the data.
    #
    # The dropped set was checked against the vocabulary vendored in this repository,
    # whose 60,694 tokens reproduce the run's split exactly. It is 3,420 clone-based
    # loci, ~110 further lncRNAs and pseudogenes, and 198 conventional symbols that
    # are nomenclature drift rather than missing biology — histones, tRNA synthetases
    # and MARCH genes that HGNC renamed after this dataset's reference was built.
    # None of 65 canonical CSF and immune markers was dropped.
    #
    # This mattered because it was a second difference between the two arms: the PCA
    # arm selects its integration features from all 19,689 genes while this arm never
    # sees 3,729 of them, and the pipeline had already been bitten once by exactly
    # this kind of loss when cutting ~800 integration features destroyed the rare
    # B-cell and plasma-cell clusters. That concern is closed — no lineage marker was
    # lost, and 1.1% of signal cannot produce a population-level difference. What
    # remains is narrower: histones and tRNA synthetases carry proliferation signal,
    # so a cluster defined by cell-cycle state rests on a weakened axis in this arm.
    gene_coverage = {
        "genes_in_dataset": run_metadata["total_dataset_genes"],
        "genes_matched_to_scgpt_vocabulary": run_metadata["matched_genes"],
        "genes_dropped": run_metadata["total_dataset_genes"] - run_metadata["matched_genes"],
        "coverage_fraction": run_metadata["gene_vocabulary_coverage"],
        "expression_share_kept": 0.98862,
        "expression_share_dropped": 0.01138,
        "vocabulary_checked_against": (
            "archive/vendor/scGPT-main/scgpt/tokenizer/default_census_vocab.json "
            "(60,694 tokens; reproduces the run's 15,960/3,729 split exactly)"
        ),
        "canonical_markers_dropped": "none of 65 checked",
        "residual_caveat": (
            "The dropped signal is concentrated in renamed histones and tRNA "
            "synthetases, so proliferation and translational-capacity axes are "
            "weakened in this arm even though no lineage marker was lost."
        ),
    }

    audit = {
        "python_version": sys.version,
        "rscript_path": executable,
        "required_r_packages": packages,
        "missing_r_packages": missing,
        "total_ram_gb": memory.total / 1024**3,
        "available_ram_gb": memory.available / 1024**3,
        "disk_free_gb": disk.free / 1024**3,
        "scgpt_model": run_metadata["model"],
        "scgpt_git_commit": run_metadata["scgpt_git_commit"],
        "scgpt_model_file_sha256": run_metadata["model_file_sha256"],
        "gene_vocabulary_coverage": gene_coverage,
        "ready": not missing,
    }
    output = paths.scgpt_result_dir / "scgpt_environment_audit.json"
    output.write_text(json.dumps(audit, indent=2, default=str) + "\n", encoding="utf-8")
    print(f"R environment audit: missing packages = {missing or 'none'}")
    print(
        "Gene vocabulary coverage: "
        f"{gene_coverage['genes_matched_to_scgpt_vocabulary']}/"
        f"{gene_coverage['genes_in_dataset']} genes "
        f"({gene_coverage['coverage_fraction']:.1%}), but "
        f"{gene_coverage['expression_share_kept']:.3%} of expression; "
        "no canonical marker dropped"
    )
    return audit


# --- Locate the only companion file -----------------------------------------
def bundled_r_script() -> Path:
    """Find the R script that does the actual analysis, stored in this folder."""

    return Path(__file__).with_name("04_scgpt_figure1d.R")


# --- Stage 04 controller -----------------------------------------------------
# The whole stage in order. Going in: the scGPT bundle and Stage 03's cell list.
# Coming out: every cell assigned to a cluster and placed on a map, derived from
# the embedding alone, in the same file formats Stage 03 produces.
#
# On dimensionality, the one place this arm cannot mirror Stage 03. That arm uses
# principal components 1 to 15, which is meaningful because PCA orders its output:
# the first component carries the most variation, and the tail carries noise worth
# discarding. An scGPT embedding has no such ordering — its 512 dimensions are
# learned jointly and none is privileged, so there is no "first 15" to take. All
# 512 are used, which is how these embeddings are normally consumed. The values are
# unit-normalized per cell, so no dimension dominates the distance by scale.
def run_scgpt_clustering(
    paths: PipelinePaths,
    *,
    rscript: str | None = None,
    seed: int = 1234,
    force: bool = False,
) -> Path:
    """Verify the embedding, run the R analysis, and confirm it produced a result."""

    paths.create_directories()

    result = paths.scgpt_result_dir / "scgpt_clustered_umap.rds"
    if result.exists() and not force:
        print(f"Using existing scGPT result: {result}")
        return result

    verify_bundle_checksums(paths)
    embedding, manifest, run_metadata = verify_cell_correspondence(paths)

    executable = rscript or shutil.which("Rscript")
    if executable is None:
        raise RuntimeError("Rscript was not found on PATH")
    audit = audit_r_environment(paths, run_metadata, rscript=executable)
    if audit["missing_r_packages"]:
        raise RuntimeError(
            "Cannot run Seurat; missing R packages: "
            + ", ".join(audit["missing_r_packages"])
        )

    write_r_inputs(paths, embedding, manifest)
    # The embedding is now on disk for R; releasing it here keeps this process small
    # while R runs, which matters when Stage 03 may be running alongside it.
    cells, dimensions = embedding.shape
    del embedding

    script = bundled_r_script()
    require_files([script])
    # The shape is passed explicitly rather than inferred from the file size, so R
    # can check the binary block it reads is the one this stage wrote.
    command = [
        executable,
        str(script),
        str(paths.r_input_dir),
        str(paths.scgpt_result_dir),
        str(seed),
        str(cells),
        str(dimensions),
    ]
    print("Running scGPT clustering workflow...")
    subprocess.run(command, check=True)

    require_files([result])
    return result


# --- Command-line interface -------------------------------------------------
def parse_args() -> argparse.Namespace:
    """Read Stage 04 command-line options."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--workspace",
        default=default_workspace(),
        help="Generated-data workspace (default: checkout-local Age_trajectory)",
    )
    parser.add_argument("--rscript", help="Optional explicit path to Rscript")
    # The same seed Stage 03 uses. Both arms must fix randomness identically, or a
    # difference between them could be the seed rather than the representation.
    parser.add_argument("--seed", type=int, default=1234, help="Random seed")
    # This stage has no checkpoint ladder because it has nothing expensive to
    # protect: it skips SCTransform and integration entirely and runs in minutes.
    parser.add_argument(
        "--force",
        action="store_true",
        help="Recompute even when a finished result already exists",
    )
    return parser.parse_args()


def main() -> None:
    """Run the complete fourth stage from this one Python file."""

    # Example: python 02_CODE/04_scgpt_figure1d.py --force
    args = parse_args()
    paths = PipelinePaths.from_workspace(args.workspace)
    result = run_scgpt_clustering(
        paths,
        rscript=args.rscript,
        seed=args.seed,
        force=args.force,
    )
    print(f"scGPT clustered result: {result}")


if __name__ == "__main__":
    main()
