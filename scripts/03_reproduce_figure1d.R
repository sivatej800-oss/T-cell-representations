#!/usr/bin/env Rscript

# The analysis itself: turn 70,391 individually measured cells from 59 patients
# into the clustered map that Figure 1D of the paper shows.
#
# The problem this solves is that the 59 samples cannot simply be pooled. Each was
# processed separately, so systematic differences between samples — batch effects —
# would dominate any comparison, and cells would group by which patient they came
# from rather than by what type of cell they are. The analysis therefore proceeds:
#
#   1. Normalize each sample on its own, so measurements become comparable within it
#      and technical effects like cell size and damage are accounted for.
#   2. Merge the 59 samples, aligning them so that the same cell type from different
#      patients lands in the same place.
#   3. Summarize the merged data down to a small number of dimensions capturing the
#      most meaningful variation.
#   4. Group cells into clusters by which other cells they most resemble.
#   5. Lay everything out in two dimensions so it can be drawn.
#
# Steps 1 and 2 take hours, so progress is saved along the way and a rerun picks up
# from the furthest point it can rather than starting over.
#
# This script is started by 03_reproduce_figure1d.py and is not run directly.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(
    "Usage: 03_reproduce_figure1d.R <seurat_input_dir> <output_dir> <seed> <restart_from>"
  )
}

input_dir <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- normalizePath(args[[2]], mustWork = FALSE)
random_seed <- as.integer(args[[3]])

# How much work to redo. Each level discards that step's checkpoint and every
# checkpoint after it, so a cheap re-run never triggers an expensive one.
#   none        resume from the most advanced checkpoint (default)
#   clustering  redo neighbors/Louvain/UMAP, keep the PCA checkpoint (minutes)
#   pca         redo PCA onward, keep the integrated object     (~30 min)
#   integration redo anchors/IntegrateData, keep the 59 SCT checkpoints (hours)
#   all         redo everything including per-sample SCTransform (many hours)
RESTART_LEVELS <- c("none", "clustering", "pca", "integration", "all")
restart_from <- tolower(args[[4]])
# Accept the previous true/false argument so an older invocation still works.
if (restart_from %in% c("true", "1", "yes")) {
  restart_from <- "all"
} else if (restart_from %in% c("false", "0", "no", "")) {
  restart_from <- "none"
}
if (!restart_from %in% RESTART_LEVELS) {
  stop(
    "restart_from must be one of: ",
    paste(RESTART_LEVELS, collapse = ", "),
    " (received: ", args[[4]], ")"
  )
}

# True when the requested restart reaches back far enough to include this step, so
# its saved result must be thrown away and recomputed. The levels are listed
# cheapest first, so asking for an expensive restart also covers the cheaper steps
# that follow it, while "none" recomputes nothing.
redo_from <- function(step) {
  match(restart_from, RESTART_LEVELS) >= match(step, RESTART_LEVELS)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
final_path <- file.path(output_dir, "integrated_clustered_umap.rds")
if (file.exists(final_path) && restart_from == "none") {
  message("Using existing integrated result: ", final_path)
  quit(save = "no", status = 0)
}
message("Restart level: ", restart_from)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})
# Work through the samples one at a time rather than several at once. Running them
# in parallel would multiply an already large memory requirement by however many
# ran together, and this analysis is limited by memory rather than by speed.
if (requireNamespace("future", quietly = TRUE)) {
  future::plan("sequential")
}
options(future.globals.maxSize = 60 * 1024^3)
# Fix the randomness up front so the same input always gives the same figure.
set.seed(random_seed)


# --- Steps 1 and 2: normalize each sample, then merge them -------------------
# Everything from reading the data to producing one merged object. This is the
# expensive part, and it only runs when there is no saved progress to resume from.
build_integrated_object <- function() {
# Read the four files Stage 02 wrote: the measurements, and the gene names, cell
# names, and per-cell details that go with them.
required_files <- file.path(
  input_dir,
  c("counts.mtx", "genes.tsv", "barcodes.tsv", "metadata.csv")
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing Seurat input files: ", paste(missing_files, collapse = ", "))
}

message("Reading sparse matrix and aligned metadata...")
counts <- Matrix::readMM(file.path(input_dir, "counts.mtx"))
counts <- as(counts, "CsparseMatrix")
genes <- readLines(file.path(input_dir, "genes.tsv"), warn = FALSE)
barcodes <- readLines(file.path(input_dir, "barcodes.tsv"), warn = FALSE)
metadata <- utils::read.csv(
  file.path(input_dir, "metadata.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Put the names back onto the measurements, and match each cell's details to the
# right column by name rather than trusting the file order. Getting this wrong
# would attach every cell to the wrong patient, so it is checked rather than assumed.
if (nrow(counts) != length(genes) || ncol(counts) != length(barcodes)) {
  stop("Matrix dimensions do not match the gene and barcode files")
}
if (anyDuplicated(genes) || anyDuplicated(barcodes) || anyDuplicated(metadata$cell_barcode)) {
  stop("Gene names, barcodes, and metadata cell identifiers must be unique")
}
rownames(counts) <- genes
colnames(counts) <- barcodes
rownames(metadata) <- metadata$cell_barcode
metadata <- metadata[barcodes, , drop = FALSE]
if (!identical(rownames(metadata), colnames(counts))) {
  stop("Metadata could not be aligned to matrix columns")
}
# The paper's patient identifiers must be present, since the merge below aligns
# everything against two specific patients named in that system.
if (!"internal_sample_id" %in% colnames(metadata)) {
  stop("metadata.csv is missing internal_sample_id; rerun stage 02")
}
if (length(unique(metadata$internal_sample_id)) != 59) {
  stop("Expected 59 internal sample IDs")
}

# Combine the measurements and the per-cell details into the object the rest of the
# analysis works on. Nothing is filtered out here: the thresholds are recalculated
# so they can be reported, but every cell continues through, matching the decision
# made in Stage 02 to leave the published dataset as its authors filtered it.
message("Creating Seurat object...")
csf <- CreateSeuratObject(
  counts = counts,
  assay = "RNA",
  project = "GSE200164_CSF",
  meta.data = metadata,
  min.cells = 0,
  min.features = 0
)
csf[["percent.mt"]] <- PercentageFeatureSet(csf, pattern = "^MT-")
if (any(csf$nFeature_RNA < 200) || any(csf$percent.mt >= 10)) {
  warning("Some cells do not meet the paper's measurable QC thresholds")
}
rm(counts, metadata)
gc()

# Separate the cells back out by patient, because the normalization that follows
# has to be done within each sample rather than across the pooled data.
message("Splitting object into 59 biological samples...")
csf_list <- SplitObject(csf, split.by = "internal_sample_id")
if (length(csf_list) != 59) {
  stop("SplitObject did not produce 59 sample objects")
}
if (sum(vapply(csf_list, ncol, numeric(1))) != ncol(csf)) {
  stop("Cell count changed while splitting samples")
}

# Step 1: normalize each patient's cells separately.
#
# Raw measurements are not comparable between cells. A cell where more RNA was
# captured shows higher numbers for every gene, regardless of biology. Normalization
# corrects for that, so a difference between two cells reflects what they were
# actually expressing rather than how much material happened to be recovered.
#
# Three known technical influences are removed at the same time: total RNA captured,
# how many distinct genes were detected, and mitochondrial share — the same three
# measurements Stage 02 calculated. Removing them prevents cells from clustering by
# how well they were captured rather than by what type of cell they are.
#
# Each sample also has its 1,000 most variable genes identified — the ones differing
# most between cells, and therefore the ones carrying information about cell type.
# Genes that behave identically everywhere cannot help distinguish anything.
#
# All settings here match the authors' published workflow. The one exception is
# processing samples one at a time and saving each result, which is a practical
# concession: this loop takes hours, and without saved progress any interruption
# would mean starting the whole thing again.
sct_dir <- file.path(output_dir, "sct_samples")
dir.create(sct_dir, recursive = TRUE, showWarnings = FALSE)
message("Running or resuming per-sample SCTransform...")
for (sample_name in names(csf_list)) {
  checkpoint <- file.path(sct_dir, paste0(sample_name, "_sct.rds"))
  # Only a full restart recomputes these; every cheaper level keeps them.
  if (file.exists(checkpoint) && !redo_from("all")) {
    message("Skipping ", sample_name, ": checkpoint exists")
    next
  }
  message("SCTransform ", sample_name, " (", ncol(csf_list[[sample_name]]), " cells)")
  transformed <- SCTransform(
    object = csf_list[[sample_name]],
    assay = "RNA",
    new.assay.name = "SCT",
    variable.features.n = 1000,
    ncells = 10000,
    vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
    # The original method rather than its later revision, matching what the authors
    # used. The two give noticeably different results, so this is not incidental.
    vst.flavor = "v1",
    conserve.memory = TRUE,
    verbose = FALSE,
    seed.use = random_seed
  )
  if (!"SCT" %in% Assays(transformed)) {
    stop("SCTransform failed to create an SCT assay for ", sample_name)
  }
  saveRDS(transformed, checkpoint, compress = FALSE)
  rm(transformed)
  gc()
}

# Gather the 59 normalized samples back together, confirming none went missing and
# that each still holds exactly the cells it started with. A resumed run reads these
# from disk, so this also guards against a checkpoint left behind by different data.
expected_samples <- names(csf_list)
sct_files <- file.path(sct_dir, paste0(expected_samples, "_sct.rds"))
if (!all(file.exists(sct_files))) {
  stop("One or more per-sample SCTransform checkpoints are missing")
}
sct_list <- lapply(sct_files, readRDS)
names(sct_list) <- expected_samples
if (sum(vapply(sct_list, ncol, numeric(1))) != ncol(csf)) {
  stop("SCTransform checkpoints do not contain all cells")
}
for (sample_name in expected_samples) {
  checkpoint_object <- sct_list[[sample_name]]
  source_object <- csf_list[[sample_name]]
  if (!"SCT" %in% Assays(checkpoint_object)) {
    stop("SCT assay missing from checkpoint for ", sample_name)
  }
  if (!identical(colnames(checkpoint_object), colnames(source_object))) {
    stop("Cell identities in the checkpoint do not match sample ", sample_name)
  }
}

# Step 2: merge the 59 samples into one dataset.
#
# The merge needs a shared vocabulary: genes that are informative across the samples
# generally, rather than within any one of them. Comparing samples on genes that are
# uninformative in most of them would align on noise.
#
# The pool of candidates is smaller than it appears. It begins as the union of all 59
# samples' 1,000 variable genes, then is cut to genes measured in every sample — and
# the smallest sample here holds 151 cells, so it contributes few. That intersection,
# not the number requested, is the real limit: asking for 3,000 returns about 2,816,
# which is the entire pool.
#
# Asking for the whole pool is the point. Candidates are ranked by how many samples
# each gene varies in, so any cut keeps genes that vary everywhere and discards genes
# that vary in only a few. Markers of rare cell types are necessarily the latter — a
# plasma-cell gene can only look variable in the samples that contained plasma cells.
# Requesting 2,000 therefore discards roughly 800 candidates and takes the rare-
# population markers with them, leaving those cells with nothing to separate on. The
# authors requested 2,000, by leaving the argument unwritten, but their 2021 software
# left them a larger pool to take 2,000 from than this version leaves.
message("Selecting integration features and preparing SCT residuals...")
integration_features <- SelectIntegrationFeatures(
  object.list = sct_list,
  nfeatures = 3000
)
message("Integration features available: ", length(integration_features))
# Receiving the full number asked for means the pool was larger than the request, so
# genes were dropped by the ranking rather than by availability — the situation this
# setting exists to avoid, and a sign the request should be raised again.
if (length(integration_features) >= 3000) {
  warning(
    "Integration features were limited by the request rather than by the data; ",
    "markers of rare cell types may have been ranked out"
  )
}
sct_list <- PrepSCTIntegration(
  object.list = sct_list,
  anchor.features = integration_features,
  verbose = TRUE
)
# Two patients, A4 and A6, serve as the reference everything else is aligned to.
# This is the authors' choice, reproduced here. It also makes the merge tractable:
# aligning every sample against these two is far less work than comparing all 59
# against each other, which grows unmanageably as samples are added.
reference_indices <- match(c("A4", "A6"), names(sct_list))
if (anyNA(reference_indices)) {
  stop("Could not locate the A4 and A6 integration references")
}

# Find anchors: pairs of cells in different samples that look enough alike to be
# taken as the same kind of cell. These pairings are what reveal how each sample is
# systematically offset from the reference, and therefore how to correct it.
#
# This step and the merge below both draw on random numbers but accept no setting to
# fix them, so they continue from wherever the preceding work happened to leave the
# sequence. That is a problem here specifically because this script resumes: a run
# that recomputed 59 normalizations arrives with a different position in the sequence
# than a run that read them from disk, so the same data could yield different anchors
# depending on how much was recomputed. Restarting the sequence immediately before
# each step makes the result a function of the data alone.
set.seed(random_seed)
message("Finding SCT integration anchors...")
integration_anchors <- FindIntegrationAnchors(
  object.list = sct_list,
  normalization.method = "SCT",
  anchor.features = integration_features,
  reference = reference_indices,
  dims = 1:30,
  verbose = TRUE
)

# Apply those corrections and produce one merged dataset. From here on, cells that
# group together should do so because they are the same type of cell, not because
# they came from the same patient. The cell count must be unchanged: merging adjusts
# the measurements, it never adds or discards cells.
set.seed(random_seed)
message("Integrating samples...")
integrated <- IntegrateData(
  anchorset = integration_anchors,
  normalization.method = "SCT",
  dims = 1:30,
  verbose = TRUE
)
if (!"integrated" %in% Assays(integrated) || ncol(integrated) != ncol(csf)) {
  stop("Integrated Seurat object failed validation")
}
rm(integration_anchors, sct_list, csf_list, csf)
gc()

integrated
}

# --- Step 3: summarize the merged data down to its main patterns -------------
# Start from the furthest point already reached. Saved progress from a previous run
# is preferred over redoing hours of work; the requested restart level decides how
# far back to go.
pca_checkpoint <- file.path(output_dir, "GSE200164_after_PCA.rds")
integrated_checkpoint <- file.path(
  output_dir,
  "GSE200164_seurat_integrated.rds"
)

if (file.exists(pca_checkpoint) && !redo_from("pca")) {
  message("Resuming from post-PCA checkpoint: ", pca_checkpoint)
  integrated <- readRDS(pca_checkpoint)
  if (!"pca" %in% Reductions(integrated)) {
    stop("Post-PCA checkpoint does not contain a PCA reduction")
  }
} else {
  if (file.exists(integrated_checkpoint) && !redo_from("integration")) {
    message("Resuming from integrated checkpoint: ", integrated_checkpoint)
    integrated <- readRDS(integrated_checkpoint)
  } else {
    integrated <- build_integrated_object()
  }
  if (!"integrated" %in% Assays(integrated)) {
    stop("Resumed object does not contain the integrated assay")
  }
  # Each cell is currently described by thousands of gene measurements, most of
  # which move together — genes involved in the same cell function rise and fall as
  # a group. This condenses that down to 30 combined dimensions, ordered so the
  # first captures the largest pattern of variation among cells, the second the next
  # largest, and so on. Later steps use only the strongest 15 of these, on the
  # reasoning that the leading patterns carry the biology while the trailing ones
  # carry noise.
  #
  # How many are calculated is not simply an upper bound on how many are used. The
  # decomposition is solved by successive approximation rather than exactly, so
  # asking for more dimensions slightly changes the ones that are kept. The effect on
  # the leading 15 is small, but it is not nothing, and 30 is the number behind the
  # closer reference run. The authors left this unwritten, which meant 50 for them.
  DefaultAssay(integrated) <- "integrated"
  set.seed(random_seed)
  integrated <- RunPCA(
    integrated,
    assay = "integrated",
    npcs = 30,
    verbose = TRUE
  )
  saveRDS(integrated, pca_checkpoint, compress = FALSE)
}


# --- Steps 4 and 5: group the cells, then lay them out for drawing -----------
# Both steps work from the 15 strongest patterns found above, the number the authors
# used. Using more would fold noise into the result; using fewer would discard real
# structure, and either choice changes how many clusters appear.
#
# Grouping happens in two parts. First each cell is linked to the cells it most
# resembles, forming a network of similarity. Then that network is divided into
# clusters by finding regions where cells are more connected to each other than to
# the rest. The resolution setting of 0.3 controls how finely it divides: higher
# values split the cells into more, smaller clusters. This value is the authors'.
#
# Finally the cells are arranged on a flat two-dimensional map for drawing, keeping
# cells that resemble each other close together. The map is a way of viewing the
# structure, not a measurement — distances across it should not be read as
# quantities, and the clustering above does not depend on it.
#
# One caveat on interpreting the output: the paper reported 13 clusters and this is
# unlikely to recover exactly that. Cluster count is a sensitive readout, and both the
# grouping and the layout depend on implementation details that have changed across
# Seurat versions, so the published settings do not guarantee a published answer on
# modern software.
#
# What the count does diagnose is whether rare populations survived. An earlier
# version of this script cut the integration features above to 2,000 and returned 11
# clusters, none smaller than 671 cells — the paper's B-cell and undetermined
# populations were absorbed rather than resolved. Taking the full pool of features
# recovered 12 clusters, the smallest holding 120 cells. A result with no small
# clusters means the rare cell types were lost upstream, and the place to look is the
# feature selection, not this resolution setting.
message("Running neighbors, Louvain clustering, and UMAP...")
DefaultAssay(integrated) <- "integrated"
set.seed(random_seed)
integrated <- FindNeighbors(integrated, reduction = "pca", dims = 1:15, verbose = TRUE)
integrated <- FindClusters(
  integrated,
  resolution = 0.3,
  algorithm = 1,
  random.seed = random_seed,
  verbose = TRUE
)
integrated <- RunUMAP(
  integrated,
  reduction = "pca",
  dims = 1:15,
  n.neighbors = 30,
  seed.use = random_seed,
  verbose = TRUE
)
# Confirm every cell received a position on the map before anything is saved.
umap <- Embeddings(integrated, reduction = "umap")
if (nrow(umap) != ncol(integrated) || ncol(umap) != 2 || anyNA(umap)) {
  stop("UMAP result failed validation")
}

# --- Save the results --------------------------------------------------------
# Three things are written: the complete analysis object, a plain table of the
# result, and the figure.
#
# The complete object is several gigabytes and holds everything — every normalized
# measurement, the merged data, the patterns, the clusters, and the map. It is what
# any further analysis would start from. As with earlier saves, it is written under
# a temporary name first so an interruption cannot destroy a previous good result.
message("Writing final checkpoint atomically...")
temporary_path <- paste0(final_path, ".tmp")
if (file.exists(temporary_path)) {
  unlink(temporary_path)
}
saveRDS(integrated, temporary_path, compress = "gzip")
# Reaching this point means a restart was requested: the early exit above already
# returned when a final result existed and restart_from was "none".
if (file.exists(final_path)) {
  unlink(final_path)
}
if (!file.rename(temporary_path, final_path)) {
  stop("Could not rename temporary checkpoint to final path")
}

# The table is the same result in a form anyone can open: one row per cell giving
# its position on the map, which cluster it was assigned to, and which patient it
# came from. A few megabytes rather than several gigabytes, and readable without R.
result_table <- data.frame(
  cell_barcode = rownames(umap),
  UMAP_1 = umap[, 1],
  UMAP_2 = umap[, 2],
  seurat_cluster = as.character(integrated$seurat_clusters),
  internal_sample_id = as.character(integrated$internal_sample_id),
  stringsAsFactors = FALSE
)
cluster_connection <- gzfile(
  file.path(output_dir, "GSE200164_umap_clusters.csv.gz"),
  open = "wt"
)
utils::write.csv(result_table, cluster_connection, row.names = FALSE)
close(cluster_connection)

# The figure: every cell drawn at its position on the map, coloured by cluster.
# This is the reproduction of Figure 1D, minus the cell-type names the authors
# assigned afterwards — the clusters here are numbered, not identified.
#
# The cells are painted onto a fixed grid of pixels rather than drawn as 70,391
# separate shapes, which keeps this fast. The grid has to be at least as fine as the
# image it is saved into, or each cell is stretched over several pixels and the map
# washes out to a pale smear. Its default is 512 by 512, far coarser than the image
# below, so the two are tied to the same numbers here rather than set independently.
figure_width_inches <- 9
figure_height_inches <- 7
figure_dpi <- 300
plot <- DimPlot(
  integrated,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  raster = TRUE,
  raster.dpi = c(
    figure_width_inches * figure_dpi,
    figure_height_inches * figure_dpi
  ),
  # Each cell is drawn three grid pixels across. The authors wrote 1, but they were
  # not painting onto a grid: their 2021 call had no raster argument, so each cell
  # became a real shape and size 1 was a visible dot. Rasterized at 2700 by 2100 the
  # same 1 covers about one pixel, so even the densest part of a cluster never
  # accumulates enough ink to read as a solid colour. Departing from their number is
  # what preserves their appearance here. Three is the largest size that darkens the
  # clusters while keeping the 120-cell cluster 10 and the thin cluster 11 wedge
  # separate from their neighbours.
  pt.size = 3
)
# Figures for every stage are collected in results/figures/ so the report can link
# them from one place; each stage's bulky data stays in its own directory.
figures_dir <- file.path(dirname(output_dir), "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = file.path(figures_dir, "03_seurat_umap_clusters.png"),
  plot = plot,
  width = figure_width_inches,
  height = figure_height_inches,
  dpi = figure_dpi
)

message("Figure 1D workflow complete: ", final_path)
