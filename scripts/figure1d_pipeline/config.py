"""URLs and generated file paths for the maintained Figure 1D pipeline.

No stage constructs its own ad-hoc path. Given one workspace directory,
``PipelinePaths`` provides every location used by the pipeline.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


# Public source identifiers are constants because all three stages target one dataset.
GEO_ACCESSION = "GSE200164"
COUNTS_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE200nnn/"
    "GSE200164/suppl/GSE200164_counts.csv.gz"
)
SERIES_MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE200nnn/"
    "GSE200164/matrix/GSE200164_series_matrix.txt.gz"
)


@dataclass(frozen=True)
class PipelinePaths:
    """Derive all generated paths from one user-selectable workspace.

    ``frozen=True`` prevents a stage from accidentally changing the workspace halfway
    through a run. The properties below return paths only; they do not create files.
    """

    workspace: Path

    @classmethod
    def from_workspace(cls, workspace: str | Path) -> "PipelinePaths":
        """Expand ``~`` and make the workspace absolute for consistent logging."""

        return cls(Path(workspace).expanduser().resolve())

    # Stage 01 data directories.
    @property
    def raw_dir(self) -> Path:
        return self.workspace / "data" / "raw"

    @property
    def interim_dir(self) -> Path:
        return self.workspace / "data" / "interim"

    @property
    def processed_dir(self) -> Path:
        return self.workspace / "data" / "processed"

    @property
    def metadata_dir(self) -> Path:
        return self.workspace / "data" / "metadata"

    # Small audit records that explain how the large artifacts were created.
    @property
    def logs_dir(self) -> Path:
        return self.workspace / "logs"

    @property
    def tables_dir(self) -> Path:
        return self.workspace / "tables"

    # Stage 02/03 interchange files, checkpoints, and final results.
    @property
    def figure1d_dir(self) -> Path:
        # Preserve the established output directory so existing multi-gigabyte
        # Seurat checkpoints remain resumable.
        return self.workspace / "results" / "figure1d_python"

    @property
    def seurat_input_dir(self) -> Path:
        return self.figure1d_dir / "seurat_input"

    # Individual file properties keep filenames in one place. If a filename changes,
    # downstream stages automatically use the same new name.
    @property
    def counts_gzip(self) -> Path:
        return self.raw_dir / "GSE200164_counts.csv.gz"

    @property
    def series_matrix_gzip(self) -> Path:
        return self.raw_dir / "GSE200164_series_matrix.txt.gz"

    @property
    def genes_by_cells_npz(self) -> Path:
        return self.interim_dir / "GSE200164_expression_genes_by_cells_csr_float32.npz"

    @property
    def cells_by_genes_npz(self) -> Path:
        return self.interim_dir / "GSE200164_expression_cells_by_genes_csr_float32.npz"

    @property
    def genes_csv(self) -> Path:
        return self.interim_dir / "GSE200164_genes.csv"

    @property
    def barcodes_txt(self) -> Path:
        return self.interim_dir / "GSE200164_barcodes.txt"

    @property
    def unannotated_h5ad(self) -> Path:
        return self.processed_dir / "GSE200164_expression_unannotated.h5ad"

    @property
    def annotated_h5ad(self) -> Path:
        return self.processed_dir / "GSE200164_expression_with_clinical_metadata.h5ad"

    @property
    def canonical_h5ad(self) -> Path:
        return self.processed_dir / "GSE200164_canonical_input.h5ad"

    @property
    def figure1d_checkpoint_h5ad(self) -> Path:
        return self.figure1d_dir / "GSE200164_after_qc_audit.h5ad"

    def create_directories(self) -> None:
        """Create the standard workspace tree; existing folders are left untouched."""

        for directory in (
            self.workspace,
            self.raw_dir,
            self.interim_dir,
            self.processed_dir,
            self.metadata_dir,
            self.logs_dir,
            self.tables_dir,
            self.figure1d_dir,
            self.seurat_input_dir,
        ):
            directory.mkdir(parents=True, exist_ok=True)


def default_workspace() -> Path:
    """Return the checkout-local workspace instead of a machine-specific path."""

    # config.py -> figure1d_pipeline/ -> scripts/ -> repository root.
    return Path(__file__).resolve().parents[2] / "Age_trajectory"
