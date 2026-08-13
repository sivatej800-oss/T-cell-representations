#!/usr/bin/env Rscript

# How quickly cluster compactness improves as Louvain is asked for finer groupings,
# in each of the two representations.
#
# This is the elbow plot's question — how much of a dataset's spread is accounted
# for by dividing it into k groups — asked of a method that has no k. Louvain takes a
# resolution rather than a cluster count, so the sweep runs resolution from 0.1 to
# 1.0 and records whatever number of clusters each setting happens to produce. The
# resulting curve is read against the number of clusters actually obtained, not
# against the resolution that produced it.
#
# Compactness is measured as the within-cluster sum of squared distances to each
# cluster's centroid, divided by the total sum of squared distances of every cell
# from the overall centroid. That ratio starts at 1 when every cell is in one
# cluster and falls toward 0 as clusters multiply, and it is what the elbow method
# plots. Its complement, 1 minus the ratio, is the fraction of the representation's
# total dispersion that the grouping accounts for.
#
# The comparison this supports is narrower than it looks. The two curves are computed
# in different spaces — 30 principal components against 512 learned dimensions — and
# squared distances in those spaces are not the same quantity. Dividing by the total
# spread puts both on a 0 to 1 scale, but it does not make the heights comparable:
# dispersion in 512 dimensions distributes differently, and a gap between the curves
# reflects that as much as it reflects compactness. What is comparable is the shape
# of each curve, meaning how steeply it falls as clusters are added, since that is a
# property of how the structure within one space is organized rather than of the
# space's dimension.
#
# Nothing is recomputed. The coordinates come from Stage 06's saved object, already
# corrected for patient, and the clustering is the same Seurat procedure at the same
# seed. Stage 06's own outputs are read only.
#
# This script is started by 08_wcss_resolution_sweep.py and is not run directly.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: 08_wcss_resolution_sweep.R <benchmark_dir> <output_dir> <seed>")
}

benchmark_dir <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- normalizePath(args[[2]], mustWork = FALSE)
random_seed <- as.integer(args[[3]])

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})
set.seed(random_seed)

# The settings Stage 06 clustered at, so that the point at resolution 0.3 reproduces
# the partition Stage 06 reported rather than merely resembling it.
RESOLUTIONS <- seq(0.1, 1.0, by = 0.1)
STAGE06_RESOLUTION <- 0.3
CLUSTER_ALGORITHM <- 1


# --- Step 1: read Stage 06's corrected coordinates ---------------------------
message("Reading Stage 06 output...")
object <- readRDS(file.path(benchmark_dir, "benchmark_representations.rds"))
reference_assignments <- read.csv(
  gzfile(file.path(benchmark_dir, "benchmark_cluster_assignments.csv.gz")),
  stringsAsFactors = FALSE
)
if (!identical(reference_assignments$cell_barcode, colnames(object))) {
  stop("Stage 06 assignments are not in the same cell order as its saved object")
}

# The two representations, each with the dimensions its arm was clustered on. Both
# are the Harmony-corrected forms, so patient differences have already been removed
# from each and the sweep measures structure rather than donor separation.
representations <- list(
  pca_harmony_30dims = list(
    reduction = "pca_harmony",
    dims = 1:30,
    reference_column = "pca_harmony_30dims"
  ),
  scgpt_harmony = list(
    reduction = "scgpt_harmony",
    dims = 1:ncol(Embeddings(object, reduction = "scgpt_harmony")),
    reference_column = "scgpt_harmony"
  )
)


# --- Step 2: the dispersion measures -----------------------------------------
# Total spread: how far every cell sits from the centre of the whole representation,
# squared and summed. This is the quantity a clustering is trying to account for, and
# it is fixed for a representation regardless of how the cells are grouped.
total_sum_of_squares <- function(coordinates) {
  overall_centroid <- colMeans(coordinates)
  sum((sweep(coordinates, 2, overall_centroid, "-"))^2)
}

# Remaining spread: how far every cell sits from the centre of its own cluster. A
# grouping that captures real structure leaves little of this; a grouping that cuts
# through a uniform cloud leaves nearly all of it.
within_cluster_sum_of_squares <- function(coordinates, labels) {
  sum(vapply(levels(labels), function(level) {
    members <- coordinates[labels == level, , drop = FALSE]
    cluster_centroid <- colMeans(members)
    sum((sweep(members, 2, cluster_centroid, "-"))^2)
  }, numeric(1)))
}


# --- Step 3: sweep the resolution in each representation ---------------------
# Each representation is clustered once per resolution. The neighbour graph depends
# only on the coordinates, not on the resolution, so it is built once per
# representation and reused across the sweep — the resolution affects how that fixed
# graph is divided, not how it is constructed.
rows <- list()
for (name in names(representations)) {
  spec <- representations[[name]]
  coordinates <- Embeddings(object, reduction = spec$reduction)[, spec$dims, drop = FALSE]
  if (nrow(coordinates) != ncol(object) || anyNA(coordinates)) {
    stop("Coordinates failed validation for ", name)
  }
  message(
    "\n", name, ": ", nrow(coordinates), " cells, ",
    ncol(coordinates), " dimensions"
  )

  total_spread <- total_sum_of_squares(coordinates)
  message("  total sum of squares: ", format(total_spread, digits = 6))

  set.seed(random_seed)
  object <- FindNeighbors(
    object,
    reduction = spec$reduction,
    dims = spec$dims,
    verbose = FALSE
  )

  for (resolution in RESOLUTIONS) {
    set.seed(random_seed)
    object <- FindClusters(
      object,
      resolution = resolution,
      algorithm = CLUSTER_ALGORITHM,
      random.seed = random_seed,
      verbose = FALSE
    )
    labels <- droplevels(Idents(object))
    within_spread <- within_cluster_sum_of_squares(coordinates, labels)
    normalized <- within_spread / total_spread
    # The ratio is bounded: cluster centroids minimize squared distance to their own
    # members, so the within-cluster total can never exceed the total about a single
    # overall centroid, and it cannot be negative.
    if (normalized < 0 || normalized > 1) {
      stop("Normalized WCSS outside [0, 1] for ", name, " at resolution ", resolution)
    }
    message(sprintf(
      "  resolution %.1f -> %2d clusters, normalized WCSS %.4f",
      resolution, nlevels(labels), normalized
    ))
    rows[[length(rows) + 1]] <- data.frame(
      representation = name,
      # Rounded so that the sweep's own step values match the setting the rest of
      # the project reports. Counting up in tenths does not land exactly on 0.3, and
      # the point marking Stage 06's resolution is located by matching that value.
      resolution = round(resolution, 1),
      n_clusters = nlevels(labels),
      normalized_wcss = normalized,
      stringsAsFactors = FALSE
    )
  }

  # The sweep passes through the setting Stage 06 used, so the cluster count there
  # must match what Stage 06 recorded. A mismatch would mean this sweep is not
  # clustering the same way Stage 06 did, and every point on the curve would be
  # describing a different procedure than the one the rest of the project reports.
  matching <- Filter(
    function(row) row$representation == name && row$resolution == STAGE06_RESOLUTION,
    rows
  )
  if (length(matching) != 1) {
    stop(
      "The sweep did not produce exactly one point at resolution ",
      STAGE06_RESOLUTION, " for ", name
    )
  }
  stage06_row <- matching[[1]]
  expected_clusters <- length(unique(reference_assignments[[spec$reference_column]]))
  if (stage06_row$n_clusters != expected_clusters) {
    stop(
      "At resolution ", STAGE06_RESOLUTION, " this sweep produced ",
      stage06_row$n_clusters, " clusters for ", name,
      " but Stage 06 recorded ", expected_clusters
    )
  }
  message(
    "  reproduces Stage 06 at resolution ", STAGE06_RESOLUTION,
    ": ", expected_clusters, " clusters"
  )
}

sweep_table <- do.call(rbind, rows)
sweep_table$dispersion_explained <- 1 - sweep_table$normalized_wcss
message("\n--- Resolution sweep ---")
print(sweep_table, row.names = FALSE)


# --- Step 4: save the table and both plots -----------------------------------
write.csv(
  sweep_table[, c("representation", "resolution", "n_clusters", "normalized_wcss")],
  file.path(output_dir, "wcss_resolution_sweep.csv"),
  row.names = FALSE
)

stage06_points <- sweep_table[sweep_table$resolution == STAGE06_RESOLUTION, ]
# Both curves must contribute one marked point. An empty selection would draw the
# plots without them and give no sign anything was missing.
if (nrow(stage06_points) != length(representations)) {
  stop(
    "Expected one point per representation at resolution ", STAGE06_RESOLUTION,
    ", found ", nrow(stage06_points)
  )
}

# Both curves against the number of clusters Louvain actually produced. The points
# from resolution 0.3 are marked because that is the setting every other stage of
# this project reports, and it locates the rest of the project on this curve.
compactness_plot <- ggplot(
  sweep_table,
  aes(n_clusters, normalized_wcss, colour = representation)
) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  geom_point(
    data = stage06_points,
    shape = 21, size = 5, stroke = 1.2, fill = NA, show.legend = FALSE
  ) +
  geom_text(
    data = stage06_points,
    aes(label = sprintf("resolution 0.3\n%d clusters", n_clusters)),
    vjust = -1.1, size = 3, show.legend = FALSE
  ) +
  scale_x_continuous(breaks = sort(unique(sweep_table$n_clusters))) +
  # Both curves lie between about 0.29 and 0.50, so the axis starts at 0.25 rather
  # than 0. Anchoring at zero would compress both into the top half of the panel and
  # hide the difference in how steeply each falls, which is what this plot is for.
  # The scale is a ratio of sums of squares, so the shortened axis exaggerates
  # apparent differences and the printed values remain the quantity to read.
  coord_cartesian(ylim = c(0.25, NA)) +
  labs(
    x = "number of clusters Louvain produced",
    y = "normalized WCSS  (within / total sum of squares)",
    title = "Cluster compactness against cluster count",
    subtitle = paste(
      "Lower is more compact. Curves are computed in different spaces",
      "(30 vs 512 dimensions);\ncompare their shape rather than their height."
    ),
    colour = NULL
  ) +
  theme_classic() +
  theme(legend.position = "bottom")
ggsave(
  file.path(output_dir, "wcss_compactness.png"),
  compactness_plot,
  width = 8,
  height = 6,
  dpi = 300
)

dispersion_plot <- ggplot(
  sweep_table,
  aes(n_clusters, dispersion_explained, colour = representation)
) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  geom_point(
    data = stage06_points,
    shape = 21, size = 5, stroke = 1.2, fill = NA, show.legend = FALSE
  ) +
  geom_text(
    data = stage06_points,
    aes(label = sprintf("resolution 0.3\n%d clusters", n_clusters)),
    vjust = -1.1, size = 3, show.legend = FALSE
  ) +
  scale_x_continuous(breaks = sort(unique(sweep_table$n_clusters))) +
  expand_limits(y = 0) +
  labs(
    x = "number of clusters Louvain produced",
    y = "fraction of total dispersion explained",
    title = "Fraction of total dispersion explained by the grouping",
    subtitle = paste(
      "Higher is more structure accounted for. Curves are computed in different",
      "spaces\n(30 vs 512 dimensions); compare their shape rather than their height."
    ),
    colour = NULL
  ) +
  theme_classic() +
  theme(legend.position = "bottom")
ggsave(
  file.path(output_dir, "wcss_dispersion_explained.png"),
  dispersion_plot,
  width = 8,
  height = 6,
  dpi = 300
)

message("\nStage 08 complete. Outputs written to ", output_dir)
