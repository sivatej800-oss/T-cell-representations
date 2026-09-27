#!/usr/bin/env Rscript

# How well-defined are the two arms' clusters, measured rather than judged by eye.
#
# Stage 06 established that the two descriptions find related but not identical
# structure, and left open which one draws better-separated groups. The Louvain
# clusters from that stage are kept exactly as they are; nothing is re-clustered and
# no k is chosen. What is added is a measurement of how well each cell sits inside
# the cluster it was given.
#
# Silhouette asks, for one cell: how far is it from the other cells in its own
# cluster, against how far it is from the cells of the nearest cluster it does not
# belong to. A cell deep inside a well-separated group scores near 1. A cell equally
# close to two groups scores near 0, meaning the boundary between them is arbitrary.
# A cell closer to another cluster than its own scores below 0, meaning it was
# probably assigned wrongly. Averaging over cells gives a number for the partition,
# and averaging within each cluster says which specific populations are solid.
#
# The comparison this enables is not as simple as scoring each arm in its own space
# and reading off the larger number. Silhouette is built from distances, and the two
# spaces are not comparable: 30 principal components and 512 learned dimensions have
# different scales, and distances in high dimensions concentrate toward a common
# value, which drags silhouette down for reasons that have nothing to do with the
# quality of the grouping. Comparing 30-dimensional scores against 512-dimensional
# ones would mostly measure the difference in dimension.
#
# So each partition is scored in both spaces, giving four numbers:
#
#                        scored in PCA space   scored in scGPT space
#   PCA labels                  native                 cross
#   scGPT labels                cross                  native
#
# Within a column the geometry is fixed, so the two partitions are directly
# comparable there and the question becomes well-posed: in this space, which set of
# labels forms the tighter groups? A partition that wins in both columns is better
# defined in a way that does not depend on where it is measured. A partition that
# wins only in its own column tells the opposite story — that each description simply
# favours the clusters it produced, which is what one would expect by construction
# and says nothing about which is right.
#
# Comparisons run down a column only. Reading across one compares numbers computed in
# different geometries, where 30 variance-ordered components and 512 learned
# dimensions produce distances that are not the same quantity.
#
# Two things this cannot answer. It does not say which clustering is biologically
# correct: a tight, well-separated group of cells is not necessarily one cell type,
# and two genuinely distinct types with similar expression will always score poorly.
# That requires marker annotation, which does not exist for these cells. And it does
# not account for the cluster counts, 13 against 11. Silhouette tends to reward
# coarser partitions, because merging two adjacent groups removes the boundary cells
# that scored badly, so some part of any advantage the 11-cluster arm shows may be
# the count rather than the separation. The partitions are scored exactly as Stage 06
# produced them and that effect is left uncontrolled; it qualifies the numbers rather
# than being removed from them.
#
# This script is started by 07_silhouette_comparison.py and is not run directly.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: 07_silhouette_comparison.R <benchmark_dir> <output_dir> <seed> <n_sample>")
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

# The two arms being compared, named as Stage 06 named them.
PCA_ARM <- "pca_harmony_30dims"
SCGPT_ARM <- "scgpt_harmony"
# Repeats of the subsample. Silhouette needs every pairwise distance within the
# sample, so the whole dataset is out of reach; drawing several times shows how much
# of any difference between the arms is just which cells were drawn.
N_REPEATS <- 5


# --- Step 1: read the finished representations and their clusters ------------
message("Reading Stage 06 output...")
object <- readRDS(file.path(benchmark_dir, "benchmark_representations.rds"))
assignments <- read.csv(
  gzfile(file.path(benchmark_dir, "benchmark_cluster_assignments.csv.gz")),
  stringsAsFactors = FALSE
)
if (!identical(assignments$cell_barcode, colnames(object))) {
  stop("Cluster assignments are not in the same cell order as the saved object")
}
for (column in c(PCA_ARM, SCGPT_ARM)) {
  if (!column %in% names(assignments)) {
    stop("Stage 06 output does not contain the arm ", column)
  }
}

# The coordinates each arm clustered on. These are the spaces the distances are
# measured in, and they are the corrected ones, so donor differences have already
# been removed from both.
spaces <- list(
  pca = Embeddings(object, reduction = "pca_harmony")[, 1:30, drop = FALSE],
  scgpt = Embeddings(object, reduction = "scgpt_harmony")
)
# The coordinates are used exactly as the clustering used them, with no rescaling.
# Stage 06 built its neighbour graphs from these same values, so measuring distances
# in anything else would score a partition in a geometry it was never formed in.
#
# Rescaling would also actively damage the PCA side. Principal components come out
# ordered by how much variation each carries, and here the first spans about ten
# times the spread of the thirtieth. Giving every component equal weight promotes the
# trailing components, which are mostly noise, to the same footing as the leading
# ones that carry the biology. Measured that way the PCA clusters score near zero
# with more than half their cells apparently misassigned; measured on the components
# as they are, the same clusters score fifteen times higher. The variance ordering is
# the information, and flattening it destroys the thing being measured.
#
# A dimension that never varies contributes nothing to any distance, so it is dropped
# rather than carried through every calculation.
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

partitions <- list(
  pca_labels = as.factor(assignments[[PCA_ARM]]),
  scgpt_labels = as.factor(assignments[[SCGPT_ARM]])
)
for (name in names(partitions)) {
  message("Partition ", name, ": ", nlevels(partitions[[name]]), " clusters")
}


# --- Step 2: score one partition in one space --------------------------------
# The sample is drawn proportionally from each cluster rather than uniformly, so that
# a cluster holding 110 cells is still represented when a cluster holding 30,413 is
# in the same draw. Uniform sampling would leave the rare populations with too few
# cells to have a meaningful within-cluster distance, and those are the populations
# the comparison most needs to see.
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
  # A cluster reduced to a single cell in the draw has no within-cluster distance,
  # and a draw containing one cluster has no other cluster to compare against.
  per_repeat <- lapply(seq_len(repeats), function(repeat_index) {
    set.seed(seed + repeat_index)
    chosen <- stratified_sample(labels, sample_size)
    drawn_labels <- droplevels(labels[chosen])
    if (nlevels(drawn_labels) < 2) {
      stop("A draw contained fewer than two clusters")
    }
    distances <- dist(coordinates[chosen, , drop = FALSE])
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


# --- Step 3: the four combinations -------------------------------------------
# Both partitions scored in both spaces. The two on the diagonal are each arm judged
# by its own geometry; the two off it are each arm judged by the other's.
message("\nScoring ", N_REPEATS, " draws of ", sample_size, " cells per combination...")
combinations <- expand.grid(
  partition = names(partitions),
  space = names(spaces),
  stringsAsFactors = FALSE
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
    is_native = (partition_name == "pca_labels") == (space_name == "pca"),
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
message("\n--- Silhouette by partition and space ---")
print(summary_table, row.names = FALSE)


# --- Step 4: save ------------------------------------------------------------
# The partitions are scored as Stage 06 produced them, 13 clusters against 11. No
# attempt is made to equalize the counts. That leaves one effect uncontrolled and it
# should be read alongside the numbers: silhouette tends to reward coarser
# partitions, because merging two adjacent groups removes the boundary cells that
# scored badly. Part of any advantage the 11-cluster arm shows may be the count
# rather than the separation, and the size of that part is not measured here.
write.csv(summary_table, file.path(output_dir, "silhouette_summary.csv"), row.names = FALSE)
write.csv(per_cluster_table, file.path(output_dir, "silhouette_per_cluster.csv"), row.names = FALSE)

# Each cluster's own score, in the space its arm was built from. A cluster near zero
# is one whose boundary the data does not really support; a large cluster scoring
# well alongside small ones scoring badly is the signature of rare populations that
# were separated by name but not by measurement.
native_per_cluster <- per_cluster_table[
  (per_cluster_table$partition == "pca_labels" & per_cluster_table$space == "pca") |
    (per_cluster_table$partition == "scgpt_labels" & per_cluster_table$space == "scgpt"),
]
native_per_cluster$arm <- ifelse(
  native_per_cluster$partition == "pca_labels", PCA_ARM, SCGPT_ARM
)
figure <- ggplot(
  native_per_cluster,
  aes(x = reorder(cluster, mean_silhouette), y = mean_silhouette, fill = cluster_size)
) +
  geom_col() +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  facet_wrap(~arm, scales = "free_x") +
  scale_fill_viridis_c(trans = "log10", name = "cells") +
  labs(x = "cluster", y = "mean silhouette width") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
# Figures for every stage are collected in results/figures/ so the report can link
# them from one place; each stage's bulky data stays in its own directory.
figures_dir <- file.path(dirname(output_dir), "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  file.path(figures_dir, "07_silhouette_per_cluster.png"),
  figure,
  width = 10,
  height = 5,
  dpi = 300
)

message("\nStage 07 complete. Outputs written to ", output_dir)
