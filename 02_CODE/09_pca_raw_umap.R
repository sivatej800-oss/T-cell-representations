#!/usr/bin/env Rscript
# Stage 09 — lay out the uncorrected PCA control for viewing.
#
# Stage 06 deliberately drew only the two Harmony-corrected arms, on the grounds
# that drawing an uncorrected arm invites reading donor separation as biology.
# That reasoning holds for interpretation, not for display: the donor structure
# in `pca_raw` is the confound Stage 06 exists to remove, and showing it beside
# the corrected panel is what makes the removal legible. The caption carries the
# warning instead of the figure being withheld.
#
# Nothing is re-clustered here. The `pca_raw` labels are read from Stage 06's
# saved assignments, so this script only computes the layout Stage 06 skipped.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: 09_pca_raw_umap.R <benchmark_dir> <figures_dir> <seed>")
}
benchmark_dir <- args[[1]]
figures_dir <- args[[2]]
random_seed <- as.integer(args[[3]])

# Stage 06's settings, restated so this script cannot silently drift from them.
UMAP_NEIGHBORS <- 30
CLUSTER_RESOLUTION <- 0.3
PCA_RAW_DIMS <- 1:15

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading Stage 06 cluster assignments")
assignments <- read.csv(
  file.path(benchmark_dir, "benchmark_cluster_assignments.csv.gz"),
  colClasses = "character"
)

message("Reading Stage 06 representations (this is a large object)")
reference <- readRDS(file.path(benchmark_dir, "benchmark_representations.rds"))

# The layout must be drawn for the same cells, in the same order, that were
# clustered. Anything else silently pairs a cell with another cell's label.
if (!identical(colnames(reference), assignments$cell_barcode)) {
  stop("Cell order in the representations does not match the saved assignments")
}

message("Building UMAP for pca_raw")
set.seed(random_seed)
reference <- RunUMAP(
  reference,
  reduction = "pca",
  dims = PCA_RAW_DIMS,
  n.neighbors = UMAP_NEIGHBORS,
  reduction.name = "umap_pca_raw",
  seed.use = random_seed,
  verbose = FALSE
)
coordinates <- Embeddings(reference, reduction = "umap_pca_raw")
if (nrow(coordinates) != ncol(reference) || anyNA(coordinates)) {
  stop("UMAP failed validation for pca_raw")
}

plot_data <- data.frame(
  UMAP_1 = coordinates[, 1],
  UMAP_2 = coordinates[, 2],
  cluster = assignments$pca_raw,
  internal_sample_id = assignments$internal_sample_id
)

# Same geometry, palette, and point settings as Stage 06's two panels, so the
# three figures can be read against one another without style differences
# standing in for real ones.
figure <- ggplot(plot_data, aes(UMAP_1, UMAP_2, colour = cluster)) +
  geom_point(size = 0.1, alpha = 0.5) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme_classic() +
  ggtitle(paste0("pca_raw, Louvain resolution ", CLUSTER_RESOLUTION))
ggsave(
  file.path(figures_dir, "06_pca_raw_umap.png"),
  figure,
  width = 8,
  height = 6,
  dpi = 300
)

# The same layout coloured by patient rather than by cluster. This is the panel
# that shows what "uncorrected" means: if the clusters above are donor structure,
# the colours here separate rather than mix.
donor_figure <- ggplot(plot_data, aes(UMAP_1, UMAP_2, colour = internal_sample_id)) +
  geom_point(size = 0.1, alpha = 0.5) +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("pca_raw, coloured by patient (59 samples, no correction)")
ggsave(
  file.path(figures_dir, "06_pca_raw_umap_by_donor.png"),
  donor_figure,
  width = 8,
  height = 6,
  dpi = 300
)

# Keep the coordinates so the layout can be redrawn without re-running UMAP.
write.csv(
  data.frame(
    cell_barcode = assignments$cell_barcode,
    pca_raw_UMAP_1 = coordinates[, 1],
    pca_raw_UMAP_2 = coordinates[, 2]
  ),
  gzfile(file.path(benchmark_dir, "benchmark_pca_raw_umap_coordinates.csv.gz")),
  row.names = FALSE
)

message("Stage 09 complete. Figures written to ", figures_dir)
