#!/usr/bin/env Rscript

# The scGPT arm of the comparison: group the same 70,391 cells and lay them out,
# using a pretrained model's description of each cell instead of one distilled from
# the expression measurements.
#
# Stage 03 spends hours arriving at its description of a cell. It normalizes each of
# the 59 samples separately, aligns them so the same cell type from different
# patients lands in the same place, and condenses thousands of gene measurements
# down to 15 combined dimensions. That work is what makes its coordinates meaningful.
#
# This script starts after the equivalent work has already happened somewhere else.
# A model pretrained on tens of millions of cells was shown each of these cells once
# and returned a 512-number description of it. No normalization, no batch alignment,
# and no dimension reduction happen here, because the embedding is already the
# reduced description — whether it also handles the batch differences that Stage 03
# corrects explicitly is one of the things this comparison is meant to reveal.
#
# From that point the two arms are the same code with the same settings:
#
#   1. Link each cell to the cells it most resembles, forming a network.
#   2. Divide that network into clusters at resolution 0.3.
#   3. Lay the cells out in two dimensions with n.neighbors 30.
#
# Because only the description of a cell differs, a difference in the result can be
# attributed to that description rather than to the software or the settings.
#
# This script is started by 04_scgpt_figure1d.py and is not run directly.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(
    "Usage: 04_scgpt_figure1d.R <r_input_dir> <output_dir> <seed> <n_cells> <n_dims>"
  )
}

input_dir <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- normalizePath(args[[2]], mustWork = FALSE)
random_seed <- as.integer(args[[3]])
expected_cells <- as.integer(args[[4]])
expected_dims <- as.integer(args[[5]])

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
final_path <- file.path(output_dir, "scgpt_clustered_umap.rds")

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})
# Fix the randomness up front so the same input always gives the same figure. This
# is the seed Stage 03 uses, because a difference between the arms must not be one.
set.seed(random_seed)


# --- Read the embedding and the cells it describes ---------------------------
# The embedding arrives as a flat block of 32-bit numbers rather than a table,
# because it is dense and R has no reader for the format it was written in. The
# numbers are in row order: the first 512 belong to the first cell, and so on.
message("Reading scGPT embedding...")
embedding_path <- file.path(input_dir, "embedding_float32.bin")
metadata_path <- file.path(input_dir, "metadata.csv")
missing_files <- c(embedding_path, metadata_path)[
  !file.exists(c(embedding_path, metadata_path))
]
if (length(missing_files) > 0) {
  stop("Missing scGPT input files: ", paste(missing_files, collapse = ", "))
}

connection <- file(embedding_path, "rb")
values <- readBin(
  connection,
  what = "numeric",
  n = expected_cells * expected_dims,
  size = 4,
  endian = "little"
)
close(connection)
if (length(values) != expected_cells * expected_dims) {
  stop("Embedding binary is not the size Stage 04 said it would be")
}
embedding <- matrix(values, nrow = expected_cells, ncol = expected_dims, byrow = TRUE)
rm(values)
# Reassembling a flat block into rows and columns is the one place the values could
# silently transpose. scGPT returns unit-length vectors, so every row should measure
# 1: a transposed read would scramble that and is caught here rather than becoming
# a plausible, wrong figure.
row_norms <- sqrt(rowSums(embedding^2))
if (any(!is.finite(embedding)) || max(abs(row_norms - 1)) > 1e-4) {
  stop("Embedding rows are not unit length; the binary block was misread")
}

metadata <- utils::read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(metadata) != expected_cells) {
  stop("Metadata row count does not match the embedding")
}
if (anyDuplicated(metadata$cell_barcode)) {
  stop("Metadata cell identifiers must be unique")
}
if (length(unique(metadata$internal_sample_id)) != 59) {
  stop("Expected 59 internal sample IDs")
}
# Row order is the correspondence between a cell and its embedding, and Stage 04
# verified it against a hash before writing these files. Naming the rows here is
# what carries that correspondence into the analysis.
rownames(embedding) <- metadata$cell_barcode
colnames(embedding) <- paste0("scgpt_", seq_len(expected_dims))
rownames(metadata) <- metadata$cell_barcode


# --- Build the object the analysis works on ----------------------------------
# Seurat organizes an analysis around measurements of features in cells, with
# reduced coordinates attached alongside. Here there are no gene measurements to
# attach: the embedding is the only description of a cell this arm uses, so it
# serves as both. The gene-level data is deliberately left out — it is untouched in
# Stage 03's export and can be joined back by cell barcode when the clusters are
# later labelled from marker expression, which is where it is actually needed.
message("Building Seurat object from the embedding...")
object <- CreateSeuratObject(
  counts = t(embedding),
  assay = "embedding",
  project = "GSE200164_CSF_scgpt",
  meta.data = metadata,
  min.cells = 0,
  min.features = 0
)
object[["scgpt"]] <- CreateDimReducObject(
  embeddings = embedding,
  key = "scgpt_",
  assay = "embedding"
)
if (ncol(object) != expected_cells) {
  stop("Seurat object does not contain every cell")
}
if (!identical(rownames(Embeddings(object, "scgpt")), colnames(object))) {
  stop("Embedding rows are not aligned to the object's cells")
}
rm(embedding)
gc()


# --- Group the cells, then lay them out for drawing --------------------------
# Every setting below is Stage 03's, unchanged, and they are the authors' settings:
# resolution 0.3 for how finely the network is divided, n.neighbors 30 for the
# layout. Holding them fixed is what makes the two arms comparable.
#
# The exception is which dimensions are used, where the two arms cannot match. Stage
# 03 takes principal components 1 to 15, and taking a leading subset is meaningful
# there because PCA orders its output by how much variation each component carries,
# so the tail is discardable noise. An embedding has no such order — its 512
# dimensions are learned together and none is privileged, so there is no leading 15
# to take. All 512 are used, as these embeddings normally are. The vectors are unit
# length, so no dimension dominates the distances through scale alone.
#
# Cluster count is a diagnostic here, exactly as it is in Stage 03, and not a target.
# It should not be expected to match either that arm's count or the paper's 13. What
# is worth reading from it is the same thing: whether small clusters appear at all,
# since their absence is the signature of rare cell types having been lost.
message("Running neighbors, Louvain clustering, and UMAP on the scGPT embedding...")
set.seed(random_seed)
object <- FindNeighbors(
  object,
  reduction = "scgpt",
  dims = 1:expected_dims,
  verbose = TRUE
)
object <- FindClusters(
  object,
  resolution = 0.3,
  algorithm = 1,
  random.seed = random_seed,
  verbose = TRUE
)
object <- RunUMAP(
  object,
  reduction = "scgpt",
  dims = 1:expected_dims,
  n.neighbors = 30,
  seed.use = random_seed,
  verbose = TRUE
)
# Confirm every cell received a position on the map before anything is saved.
umap <- Embeddings(object, reduction = "umap")
if (nrow(umap) != ncol(object) || ncol(umap) != 2 || anyNA(umap)) {
  stop("UMAP result failed validation")
}
message("Clusters found: ", nlevels(object$seurat_clusters))
print(table(object$seurat_clusters))


# --- Save the results --------------------------------------------------------
# The same three outputs Stage 03 writes, in the same formats and under matching
# names, so the two arms can be compared without any conversion step.
#
# The object is written under a temporary name first, so an interruption cannot
# destroy a previous good result.
message("Writing final checkpoint atomically...")
temporary_path <- paste0(final_path, ".tmp")
if (file.exists(temporary_path)) {
  unlink(temporary_path)
}
saveRDS(object, temporary_path, compress = "gzip")
if (file.exists(final_path)) {
  unlink(final_path)
}
if (!file.rename(temporary_path, final_path)) {
  stop("Could not rename temporary checkpoint to final path")
}

# One row per cell giving its position on the map, its cluster, and its patient.
# The columns match Stage 03's table exactly, and because both arms hold the cells
# in the same order the two tables can be set side by side as they are.
result_table <- data.frame(
  cell_barcode = rownames(umap),
  UMAP_1 = umap[, 1],
  UMAP_2 = umap[, 2],
  scgpt_cluster = as.character(object$seurat_clusters),
  internal_sample_id = as.character(object$internal_sample_id),
  stringsAsFactors = FALSE
)
cluster_connection <- gzfile(
  file.path(output_dir, "GSE200164_scgpt_umap_clusters.csv.gz"),
  open = "wt"
)
utils::write.csv(result_table, cluster_connection, row.names = FALSE)
close(cluster_connection)

# The figure: the same drawing as Stage 03's, so the two can be set beside each
# other and read the same way. The cells are painted onto a grid at least as fine as
# the image they are saved into, since the default grid is far coarser and washes
# the map out to a pale smear for reasons unrelated to the analysis.
figure_width_inches <- 9
figure_height_inches <- 7
figure_dpi <- 300
plot <- DimPlot(
  object,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  raster = TRUE,
  raster.dpi = c(
    figure_width_inches * figure_dpi,
    figure_height_inches * figure_dpi
  ),
  # Two rather than Stage 03's three. This embedding separates its clusters into
  # compact, well-spaced islands, so it needs less ink than Stage 03's single
  # overlapping mass to read clearly, and a smaller size keeps the sparse cluster 12
  # filament and the cluster 11 arc visible rather than thickening them.
  pt.size = 2
) + ggtitle("scGPT embedding, Louvain resolution 0.3")
# Figures for every stage are collected in results/figures/ so the report can link
# them from one place; each stage's bulky data stays in its own directory.
figures_dir <- file.path(dirname(output_dir), "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = file.path(figures_dir, "04_scgpt_raw_umap_clusters.png"),
  plot = plot,
  width = figure_width_inches,
  height = figure_height_inches,
  dpi = figure_dpi
)

message("scGPT clustering workflow complete: ", final_path)
