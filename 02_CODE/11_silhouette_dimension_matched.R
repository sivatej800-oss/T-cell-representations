#!/usr/bin/env Rscript

# Silhouette with dimensionality held fixed, so the comparison isolates the
# representation.
#
# Stage 07's result was that every partition scores higher in scGPT coordinates
# than in PCA coordinates — including the PCA partition itself, which improved
# seventeenfold simply by being measured somewhere else. Read one way that says the
# embedding organizes these cells more cleanly. Read another it says nothing about
# the representations at all, because the two spaces differ in dimension as well as
# in kind, and silhouette responds to dimension on its own: distances in high
# dimensions concentrate toward a common value, which changes every score built
# from them.
#
# Stage 07 cannot separate those readings, because 30 against 512 varies
# representation and dimensionality together. This stage separates them by scoring
# PCA-512 against scGPT-512, where the only remaining difference is how the
# coordinates were produced.
#
#                        scored in PCA-512 space   scored in scGPT-512 space
#   PCA-512 labels               native                      cross
#   scGPT-512 labels             cross                       native
#
# Both spaces now have 512 dimensions, both partitions were clustered at the same
# resolution with the same seed on the same cells, and both were Harmony-corrected
# for donor. Comparisons still run down a column, for the same reason as in Stage
# 07: within a column the geometry is fixed, and the question "in this space, which
# labels form tighter groups" is well-posed. What is new is that reading across a
# row is no longer meaningless, because the two geometries now share a dimension.
#
# The PCA-30 partition is carried along as a third set of labels, scored in both
# 512-dimensional spaces. It is not part of the matched comparison; it is there to
# show what changes when only the dimensionality of the scoring space changes and
# the labels stay the same.
#
# The same two limits from Stage 07 remain. Cluster counts are not equalized and
# silhouette rewards coarser partitions. And silhouette measures geometric
# separation, not biological correctness, which needs marker annotation these cells
# do not have.
#
# This script is started by 11_silhouette_dimension_matched.py.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: 11_silhouette_dimension_matched.R <benchmark_dir> <output_dir> <seed> <n_sample>")
}

benchmark_dir <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- normalizePath(args[[2]], mustWork = FALSE)
random_seed <- as.integer(args[[3]])
sample_size <- as.integer(args[[4]])

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(cluster)
  library(ggplot2)
})
set.seed(random_seed)

N_REPEATS <- 5
PCA30_ARM <- "pca_harmony_30dims"
SCGPT_ARM <- "scgpt_harmony"


# --- Step 1: gather the three partitions and the two 512-dimensional spaces ---
message("Reading Stage 06 output...")
object <- readRDS(file.path(benchmark_dir, "benchmark_representations.rds"))
assignments <- read.csv(
  gzfile(file.path(benchmark_dir, "benchmark_cluster_assignments.csv.gz")),
  stringsAsFactors = FALSE
)
if (!identical(assignments$cell_barcode, colnames(object))) {
  stop("Stage 06 assignments are not in the same cell order as the saved object")
}

message("Reading Stage 10 PCA-512 arm...")
pca512_object <- readRDS(file.path(benchmark_dir, "benchmark_pca512_arm.rds"))
pca512_assignments <- read.csv(
  gzfile(file.path(benchmark_dir, "benchmark_pca512_assignments.csv.gz")),
  stringsAsFactors = FALSE
)
# Three files describing the same cells, each written by a different run. Any one of
# them silently out of order would still score and still produce a plausible figure,
# so the order is checked rather than assumed.
if (!identical(pca512_assignments$cell_barcode, colnames(pca512_object))) {
  stop("Stage 10 assignments are not in the same cell order as the Stage 10 object")
}
if (!identical(pca512_assignments$cell_barcode, assignments$cell_barcode)) {
  stop("Stage 10 and Stage 06 do not describe the same cells in the same order")
}

spaces <- list(
  pca512 = Embeddings(pca512_object, reduction = "pca512_harmony"),
  scgpt = Embeddings(object, reduction = "scgpt_harmony")
)
rm(pca512_object)
gc()

# The coordinates are used exactly as each arm clustered on them, with no rescaling.
# Stage 07 documents why: flattening PCA's variance ordering promotes trailing noise
# components to the weight of the leading ones and destroys the thing being measured.
spaces <- lapply(spaces, function(coordinates) {
  spread <- apply(coordinates, 2, sd)
  constant <- spread == 0 | is.na(spread)
  if (any(constant)) {
    message("Dropping ", sum(constant), " dimensions with no variation")
    coordinates <- coordinates[, !constant, drop = FALSE]
  }
  coordinates
})
for (name in names(spaces)) {
  if (nrow(spaces[[name]]) != ncol(object) || anyNA(spaces[[name]])) {
    stop("Coordinates failed validation for space ", name)
  }
  message("Space ", name, ": ", ncol(spaces[[name]]), " dimensions")
}
# The comparison's whole premise is that the two spaces have the same width. If they
# do not, the confound Stage 07 left open is still present and the result would be
# read as controlled when it is not.
if (ncol(spaces$pca512) != ncol(spaces$scgpt)) {
  stop(
    "Spaces are not dimensionality-matched: pca512 has ", ncol(spaces$pca512),
    " and scgpt has ", ncol(spaces$scgpt),
    ". The comparison this stage exists to make is not available."
  )
}

partitions <- list(
  pca512_labels = as.factor(pca512_assignments$pca512_harmony),
  scgpt_labels = as.factor(assignments[[SCGPT_ARM]]),
  pca30_labels = as.factor(assignments[[PCA30_ARM]])
)
for (name in names(partitions)) {
  message("Partition ", name, ": ", nlevels(partitions[[name]]), " clusters")
}


# --- Step 2: score one partition in one space --------------------------------
# R's dist() is a single-threaded C loop, which leaves fifteen of sixteen cores idle
# on a machine bought for the other fifteen. The same Euclidean distances come out of
# the identity ||a-b||^2 = ||a||^2 + ||b||^2 - 2(a.b), where the only expensive term
# is a matrix product that OpenBLAS threads across every core. Measured on a
# 4,000-cell draw in 512 dimensions this is about twenty times faster, and the two
# routes agree to around 1e-06 — floating-point noise, far below anything silhouette
# resolves.
#
# pmax(., 0) guards the one numerical trap: for two cells at nearly the same
# position the subtraction can land a hair below zero, and sqrt() of that is NaN.
# Clamping restores the true value of zero.
euclidean_distances <- function(coordinates) {
  squared_norms <- rowSums(coordinates^2)
  gram <- tcrossprod(coordinates)
  squared <- outer(squared_norms, squared_norms, "+") - 2 * gram
  as.dist(sqrt(pmax(squared, 0)))
}


# Drawn proportionally from each cluster, as Stage 07 does, so a 110-cell cluster is
# still represented in a draw that also holds a 30,000-cell one.
stratified_sample <- function(labels, size) {
  counts <- table(labels)
  per_cluster <- pmax(2, round(size * as.vector(counts) / length(labels)))
  names(per_cluster) <- names(counts)
  unlist(lapply(levels(labels), function(level) {
    available <- which(labels == level)
    sample(available, min(length(available), per_cluster[[level]]))
  }))
}

score_partition <- function(coordinates, labels, repeats, seed) {
  per_repeat <- lapply(seq_len(repeats), function(repeat_index) {
    set.seed(seed + repeat_index)
    chosen <- stratified_sample(labels, sample_size)
    drawn_labels <- droplevels(labels[chosen])
    if (nlevels(drawn_labels) < 2) {
      stop("A draw contained fewer than two clusters")
    }
    distances <- euclidean_distances(coordinates[chosen, , drop = FALSE])
    result <- silhouette(as.integer(drawn_labels), distances)
    per_cluster <- tapply(result[, "sil_width"], drawn_labels, mean)
    list(
      overall = mean(result[, "sil_width"]),
      negative_fraction = mean(result[, "sil_width"] < 0),
      per_cluster = per_cluster
    )
  })

  overall <- vapply(per_repeat, function(x) x$overall, numeric(1))
  negative <- vapply(per_repeat, function(x) x$negative_fraction, numeric(1))
  cluster_names <- levels(labels)
  per_cluster <- vapply(cluster_names, function(name) {
    mean(vapply(per_repeat, function(x) {
      value <- x$per_cluster[[name]]
      if (is.null(value)) NA_real_ else value
    }, numeric(1)), na.rm = TRUE)
  }, numeric(1))

  list(
    mean_silhouette = mean(overall),
    sd_silhouette = sd(overall),
    negative_fraction = mean(negative),
    per_cluster = per_cluster
  )
}


# --- Step 3: every partition in both matched spaces --------------------------
message("\nScoring ", N_REPEATS, " draws of ", sample_size, " cells per combination...")
combinations <- expand.grid(
  partition = names(partitions),
  space = names(spaces),
  stringsAsFactors = FALSE
)

native_space <- c(
  pca512_labels = "pca512",
  scgpt_labels = "scgpt",
  pca30_labels = NA_character_
)

results <- list()
per_cluster_rows <- list()
for (row in seq_len(nrow(combinations))) {
  partition_name <- combinations$partition[[row]]
  space_name <- combinations$space[[row]]
  message("  ", partition_name, " in ", space_name, " space")
  scored <- score_partition(
    spaces[[space_name]],
    partitions[[partition_name]],
    N_REPEATS,
    random_seed
  )
  results[[row]] <- data.frame(
    partition = partition_name,
    space = space_name,
    # The PCA-30 labels have no native space here; both 512-dimensional spaces are
    # foreign to them, which is the point of carrying them along.
    is_native = identical(native_space[[partition_name]], space_name),
    is_matched = partition_name %in% c("pca512_labels", "scgpt_labels"),
    n_clusters = nlevels(partitions[[partition_name]]),
    mean_silhouette = scored$mean_silhouette,
    sd_across_draws = scored$sd_silhouette,
    negative_fraction = scored$negative_fraction,
    stringsAsFactors = FALSE
  )
  per_cluster_rows[[row]] <- data.frame(
    partition = partition_name,
    space = space_name,
    cluster = names(scored$per_cluster),
    cluster_size = as.vector(table(partitions[[partition_name]])),
    mean_silhouette = as.vector(scored$per_cluster),
    stringsAsFactors = FALSE
  )
}

summary_table <- do.call(rbind, results)
per_cluster_table <- do.call(rbind, per_cluster_rows)
message("\n--- Silhouette at matched dimensionality ---")
print(summary_table, row.names = FALSE)


# --- Step 4: save ------------------------------------------------------------
write.csv(
  summary_table,
  file.path(output_dir, "silhouette_matched_summary.csv"),
  row.names = FALSE
)
write.csv(
  per_cluster_table,
  file.path(output_dir, "silhouette_matched_per_cluster.csv"),
  row.names = FALSE
)

# Figures for every stage are collected in results/figures/ so the report can link
# them from one place; each stage's bulky data stays in its own directory.
figures_dir <- file.path(dirname(output_dir), "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# The 2x2 that is the point of the stage: two matched partitions, two matched
# spaces. Columns are geometries, so bars are compared within a panel.
matched <- summary_table[summary_table$is_matched, ]
matched$space_label <- ifelse(
  matched$space == "pca512", "scored in PCA-512 space", "scored in scGPT-512 space"
)
matched$partition_label <- ifelse(
  matched$partition == "pca512_labels", "PCA-512 labels", "scGPT-512 labels"
)
figure <- ggplot(
  matched,
  aes(x = partition_label, y = mean_silhouette, fill = is_native)
) +
  geom_col(width = 0.6) +
  geom_errorbar(
    aes(ymin = mean_silhouette - sd_across_draws, ymax = mean_silhouette + sd_across_draws),
    width = 0.15
  ) +
  facet_wrap(~space_label) +
  scale_fill_manual(
    values = c("TRUE" = "#4a7a4a", "FALSE" = "#8a8a8a"),
    labels = c("TRUE" = "native", "FALSE" = "cross"),
    name = NULL
  ) +
  labs(
    x = NULL,
    y = "mean silhouette width",
    title = "Silhouette at matched dimensionality (512 vs 512)",
    subtitle = "compare bars within a panel; the geometry is fixed there"
  ) +
  theme_classic()
ggsave(
  file.path(figures_dir, "11_silhouette_matched.png"),
  figure,
  width = 9,
  height = 5,
  dpi = 300
)

# Per-cluster detail for the two matched partitions in their own spaces.
native_per_cluster <- per_cluster_table[
  (per_cluster_table$partition == "pca512_labels" & per_cluster_table$space == "pca512") |
    (per_cluster_table$partition == "scgpt_labels" & per_cluster_table$space == "scgpt"),
]
cluster_figure <- ggplot(
  native_per_cluster,
  aes(x = reorder(cluster, mean_silhouette), y = mean_silhouette, fill = cluster_size)
) +
  geom_col() +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  facet_wrap(~partition, scales = "free_x") +
  scale_fill_viridis_c(trans = "log10", name = "cells") +
  labs(x = "cluster", y = "mean silhouette width") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
ggsave(
  file.path(figures_dir, "11_silhouette_matched_per_cluster.png"),
  cluster_figure,
  width = 10,
  height = 5,
  dpi = 300
)

message("\nStage 11 complete. Outputs written to ", output_dir)
