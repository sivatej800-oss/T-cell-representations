#!/usr/bin/env Rscript

# A PCA arm with the same number of dimensions as the embedding, to separate
# representation from dimensionality.
#
# Stage 07 compared a 13-cluster PCA partition in 30 dimensions against an
# 11-cluster scGPT partition in 512, and found that every partition scores higher
# in scGPT coordinates than in PCA coordinates — including the PCA partition
# itself, which improved seventeenfold simply by being measured elsewhere. That
# result has two candidate explanations and Stage 07 cannot separate them:
#
#   1. the embedding genuinely organizes these cells more cleanly, or
#   2. silhouette behaves differently in 512 dimensions than in 30.
#
# Explanation 2 is not a minor worry. Distances in high dimensions concentrate
# toward a common value, and silhouette is built entirely from distances, so the
# number of dimensions moves the score on its own. Comparing 30 against 512 leaves
# representation and dimensionality varying together.
#
# This stage removes that confound the only way it can be removed: by building a
# PCA arm at 512 dimensions and clustering it fresh, so the comparison against
# scGPT holds dimensionality fixed and varies only the representation.
#
# The clustering is genuinely new. Reusing the 30-dimension labels and merely
# re-scoring them would answer a different question — how the same partition
# behaves in a larger space — and would not produce the PCA-512 partition this
# comparison needs. Everything downstream of the decomposition matches Stage 06
# exactly: same Harmony, same neighbour graph settings, same resolution, same seed.
#
# What this stage cannot fix: a 2,815-gene pool does not carry 512 dimensions of
# real structure, so the trailing components are close to noise. That is a true
# property of PCA at this width rather than a flaw in the setup, and it is the
# honest answer to "what does PCA look like at 512 dimensions". The variance
# actually carried past component 30 is reported so the result can be read with it.
#
# This script is started by 10_pca512_arm.py and is not run directly.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: 10_pca512_arm.R <benchmark_dir> <figures_dir> <seed> <n_dims> <n_cells>")
}

benchmark_dir <- normalizePath(args[[1]], mustWork = TRUE)
figures_dir <- normalizePath(args[[2]], mustWork = FALSE)
random_seed <- as.integer(args[[3]])
n_dims <- as.integer(args[[4]])
expected_cells <- as.integer(args[[5]])

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
})

# Stage 06's settings, restated so this arm cannot silently drift from the ones it
# is meant to be comparable to.
CLUSTER_RESOLUTION <- 0.3
UMAP_NEIGHBORS <- 30

checkpoint <- file.path(benchmark_dir, "benchmark_pca512_arm.rds")

if (file.exists(checkpoint)) {
  message("Resuming from ", checkpoint)
  reference <- readRDS(checkpoint)
} else {
  # The uncorrected checkpoint carries the merged, scaled expression Stage 06 built:
  # 59 samples SCTransformed separately then merged without anchors, so the only
  # donor correction applied anywhere is the Harmony below.
  message("Reading Stage 06 uncorrected checkpoint...")
  reference <- readRDS(file.path(benchmark_dir, "benchmark_pca_uncorrected.rds"))
  if (ncol(reference) != expected_cells) {
    stop("Checkpoint holds ", ncol(reference), " cells, expected ", expected_cells)
  }

  features <- VariableFeatures(reference)
  message("Feature pool: ", length(features), " genes")
  # PCA cannot return more components than the data has genes or cells.
  if (n_dims > length(features)) {
    stop("Cannot compute ", n_dims, " components from ", length(features), " genes")
  }

  # Recomputed rather than extended from the saved 50. Seurat's PCA is approximate,
  # and how many components are requested slightly perturbs the leading ones, so a
  # 512-component run is its own decomposition rather than the 50-component one with
  # more columns appended.
  #
  # approx = FALSE forces exact SVD instead of irlba. irlba is a truncated method
  # built for taking a few components from a large matrix, and it stops converging
  # when asked for a large fraction of them: at 512 of a possible 2,815 it returned
  # "did not converge--results might be invalid!" and produced components that
  # cannot be trusted. Exact SVD has no such limit. It costs more, but the whole
  # point of this arm is the trailing components, so they have to be right.
  message("Running PCA to ", n_dims, " components, exact SVD (this is the slow step)...")
  set.seed(random_seed)
  pca_warnings <- character(0)
  reference <- withCallingHandlers(
    RunPCA(
      reference,
      features = features,
      npcs = n_dims,
      approx = FALSE,
      reduction.name = "pca512",
      reduction.key = "PC512_",
      verbose = FALSE
    ),
    warning = function(w) {
      pca_warnings <<- c(pca_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  # A decomposition that did not converge still returns a full matrix of numbers,
  # and every check below it — dimensions, no missing values — passes on that
  # matrix. Nothing downstream would reveal the problem: Harmony would correct it,
  # Louvain would cluster it, and the figure would look entirely ordinary. The
  # warning is the only evidence, so it is promoted to a hard failure here.
  convergence_failures <- grep("converge", pca_warnings, value = TRUE, ignore.case = TRUE)
  if (length(convergence_failures) > 0) {
    stop(
      "PCA did not converge, so the components are not trustworthy:\n  ",
      paste(convergence_failures, collapse = "\n  ")
    )
  }
  for (warning_text in pca_warnings) {
    message("PCA warning: ", warning_text)
  }
  saveRDS(reference, checkpoint, compress = FALSE)
}

embeddings <- Embeddings(reference, reduction = "pca512")
if (ncol(embeddings) != n_dims || nrow(embeddings) != expected_cells) {
  stop("PCA did not return ", n_dims, " components for ", expected_cells, " cells")
}

# How much of the data's variation the trailing components actually carry. This is
# the number that says whether "PCA at 512 dimensions" is mostly signal or mostly
# noise, and it qualifies every score computed downstream.
variation <- apply(embeddings, 2, var)
total_variation <- sum(variation)
share_first_30 <- sum(variation[1:30]) / total_variation
share_past_30 <- sum(variation[31:n_dims]) / total_variation
message(
  "Variance captured: PCs 1-30 hold ", round(100 * share_first_30, 2),
  "%, PCs 31-", n_dims, " hold ", round(100 * share_past_30, 2), "%"
)
# A component with no variation at all contributes nothing to any distance and
# indicates the decomposition ran past the rank of the data.
degenerate <- sum(variation <= 0 | is.na(variation))
if (degenerate > 0) {
  message("WARNING: ", degenerate, " components carry no variation")
}

# --- Harmony, neighbours, clustering, UMAP: Stage 06's settings exactly -------
if (!"pca512_harmony" %in% Reductions(reference)) {
  message("Running Harmony on all ", n_dims, " components...")
  set.seed(random_seed)
  reference <- RunHarmony(
    reference,
    group.by.vars = "internal_sample_id",
    reduction.use = "pca512",
    reduction.save = "pca512_harmony",
    dims.use = 1:n_dims,
    verbose = TRUE
  )
  corrected <- Embeddings(reference, reduction = "pca512_harmony")
  if (nrow(corrected) != expected_cells || anyNA(corrected)) {
    stop("Harmony returned invalid coordinates for the PCA-512 arm")
  }
  saveRDS(reference, checkpoint, compress = FALSE)
}

message("Clustering on all ", n_dims, " corrected dimensions...")
set.seed(random_seed)
reference <- FindNeighbors(
  reference,
  reduction = "pca512_harmony",
  dims = 1:n_dims,
  verbose = FALSE
)
reference <- FindClusters(
  reference,
  resolution = CLUSTER_RESOLUTION,
  algorithm = 1,
  random.seed = random_seed,
  verbose = FALSE
)
labels <- Idents(reference)
sizes <- as.vector(table(labels))
message(
  "PCA-512 produced ", length(sizes), " clusters, smallest ", min(sizes),
  ", largest ", max(sizes)
)

message("Building UMAP...")
set.seed(random_seed)
reference <- RunUMAP(
  reference,
  reduction = "pca512_harmony",
  dims = 1:n_dims,
  n.neighbors = UMAP_NEIGHBORS,
  reduction.name = "umap_pca512",
  seed.use = random_seed,
  verbose = FALSE
)
coordinates <- Embeddings(reference, reduction = "umap_pca512")
if (nrow(coordinates) != expected_cells || anyNA(coordinates)) {
  stop("UMAP failed validation for the PCA-512 arm")
}

saveRDS(reference, checkpoint, compress = FALSE)

# --- Save the labels and the layout ------------------------------------------
# Written in the object's cell order, and read back against Stage 06's assignments
# by the silhouette stage rather than assumed to align.
output <- data.frame(
  cell_barcode = colnames(reference),
  pca512_harmony = as.character(labels),
  pca512_UMAP_1 = coordinates[, 1],
  pca512_UMAP_2 = coordinates[, 2],
  stringsAsFactors = FALSE
)
write.csv(
  output,
  gzfile(file.path(benchmark_dir, "benchmark_pca512_assignments.csv.gz")),
  row.names = FALSE
)

summary_table <- data.frame(
  arm = "pca512_harmony",
  dims_used = n_dims,
  n_clusters = length(sizes),
  smallest_cluster = min(sizes),
  largest_cluster = max(sizes),
  variance_share_first_30 = share_first_30,
  variance_share_past_30 = share_past_30,
  stringsAsFactors = FALSE
)
write.csv(
  summary_table,
  file.path(benchmark_dir, "benchmark_pca512_summary.csv"),
  row.names = FALSE
)

# Same geometry and point settings as Stage 06's panels, so the figures can be read
# against one another without style differences standing in for real ones.
plot_data <- data.frame(
  UMAP_1 = coordinates[, 1],
  UMAP_2 = coordinates[, 2],
  cluster = as.character(labels)
)
figure <- ggplot(plot_data, aes(UMAP_1, UMAP_2, colour = cluster)) +
  geom_point(size = 0.1, alpha = 0.5) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme_classic() +
  ggtitle(paste0("pca512_harmony, Louvain resolution ", CLUSTER_RESOLUTION))
ggsave(
  file.path(figures_dir, "10_pca512_harmony_umap.png"),
  figure,
  width = 8,
  height = 6,
  dpi = 300
)

message("Stage 10 complete.")
