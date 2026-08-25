# Stage 13: what the two frozen partitions are, biologically.
#
# Stages 06 through 11 compared the two representations in their own terms, and could
# not say whether either one's communities are real cell types. This stage asks the
# RNA that question directly.
#
# The design is one expression object, two sets of labels. Both partitions are marker-
# tested against the identical matrix, the identical assay and layer, and the identical
# thresholds; the only thing that differs between the two runs is which cells are
# called cluster 0, cluster 1, and so on. Because the evidence comes from the RNA and
# not from either representation, neither arm gets to grade its own work.
#
# The expression is the RNA assay, log-normalized, not SCT. The authors set the
# default assay to SCT for their marker call, but their SCT ran on one merged object
# and so covered one gene set. Here SCTransform ran per sample, and the 59 checkpoints
# share only 9,013 of 19,689 genes -- the genes that drop are the ones measured in
# some samples but not all, which is exactly what a rare population's markers look
# like. Taking the intersection would repeat the failure the 2,000-feature cut caused
# in Stage 03. RNA carries every gene for every cell, is what the authors themselves
# used for their final marker heatmap, and is the same values for both arms, which is
# what the comparison actually requires.
#
# This script performs no clustering. The labels arrive frozen from Stage 12 and are
# used exactly as given.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

# presto replaces Seurat's own Wilcoxon implementation with a compiled one. It is
# required rather than optional because the two partitions must be tested by the same
# code: if one arm ran under presto and the other did not, a difference between their
# marker lists could be the implementation rather than the biology. Requiring it here
# makes that impossible to get wrong by accident.
if (!requireNamespace("presto", quietly = TRUE)) {
  stop("presto is not installed. Both partitions must be tested by the same ",
       "implementation, so this stage will not run without it.")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop("Usage: 13_marker_discovery.R <seurat_input_dir> <frozen_dir> <output_dir> <seed> <min_pct> <logfc_threshold>")
}
input_dir <- args[1]
frozen_dir <- args[2]
output_dir <- args[3]
seed <- as.integer(args[4])
min_pct <- as.numeric(args[5])
logfc_threshold <- as.numeric(args[6])

expected_cells <- 70391
expected_genes <- 19689

# The two frozen label columns, named as Stage 12 named them.
partitions <- c("pca30_harmony_cluster", "scgpt_harmony_cluster")

# The panel the paper used to name its populations, in the paper's order.
marker_panel <- c(
  "CD3D", "CD4", "CD8A", "CD8B", "FOXP3",
  "NKG7", "GNLY", "FCER1A", "CD68",
  "CD14", "FCGR3A", "CD19", "JCHAIN"
)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(seed)

# --- Step 1: build the one expression object both partitions will be judged on ----
# Read the same Matrix Market hand-off Stage 03 reads, so that the expression here is
# the published data and not a downstream product of either arm.
message("Reading the RNA matrix...")
counts <- Matrix::readMM(file.path(input_dir, "counts.mtx"))
counts <- as(counts, "CsparseMatrix")
genes <- readLines(file.path(input_dir, "genes.tsv"), warn = FALSE)
barcodes <- readLines(file.path(input_dir, "barcodes.tsv"), warn = FALSE)

if (nrow(counts) != length(genes) || ncol(counts) != length(barcodes)) {
  stop("Matrix dimensions do not match the gene and barcode files")
}
if (nrow(counts) != expected_genes || ncol(counts) != expected_cells) {
  stop("Expected ", expected_genes, " x ", expected_cells, ", found ",
       nrow(counts), " x ", ncol(counts))
}
if (anyDuplicated(genes) || anyDuplicated(barcodes)) {
  stop("Gene names and barcodes must be unique")
}
rownames(counts) <- genes
colnames(counts) <- barcodes

# --- Step 2: attach both frozen label sets to that one object --------------------
# Joined by barcode rather than by position. The files are known to agree in order,
# but a positional join that silently succeeded on mismatched files would attach
# marker expression to the wrong cells and produce a clean, wrong annotation.
message("Reading frozen labels from Stage 12...")
frozen <- utils::read.csv(
  file.path(frozen_dir, "frozen_cluster_labels.csv.gz"),
  stringsAsFactors = FALSE,
  colClasses = "character"
)
for (column in c("cell_barcode", partitions)) {
  if (!column %in% colnames(frozen)) {
    stop("Frozen label file is missing the column ", column)
  }
}
if (nrow(frozen) != expected_cells) {
  stop("Frozen labels cover ", nrow(frozen), " cells, expected ", expected_cells)
}
rownames(frozen) <- frozen$cell_barcode
if (!setequal(frozen$cell_barcode, barcodes)) {
  stop("The frozen labels and the RNA matrix do not describe the same cells")
}
frozen <- frozen[barcodes, , drop = FALSE]
if (!identical(rownames(frozen), colnames(counts))) {
  stop("Frozen labels could not be aligned to the matrix columns")
}
if (any(is.na(frozen[, partitions]))) {
  stop("A cell is missing a cluster label in one of the partitions")
}

message("Creating the shared Seurat object...")
csf <- CreateSeuratObject(
  counts = counts,
  assay = "RNA",
  project = "GSE200164_markers",
  meta.data = frozen[, partitions, drop = FALSE],
  min.cells = 0,
  min.features = 0
)
rm(counts)
gc(verbose = FALSE)

if (ncol(csf) != expected_cells || nrow(csf) != expected_genes) {
  stop("The shared object does not hold the full matrix")
}

# Log-normalize once. Both partitions read these same values; normalizing once rather
# than per arm is what makes "same RNA" literally true rather than approximately so.
message("Log-normalizing (once, shared by both partitions)...")
csf <- NormalizeData(csf, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
DefaultAssay(csf) <- "RNA"

# The expression both runs will see, captured before either run touches it, so that
# the claim that the two runs shared an assay can be checked rather than trusted.
expression_fingerprint <- list(
  assay = "RNA",
  layer = "data",
  n_genes = nrow(csf),
  n_cells = ncol(csf),
  data_sum = sum(GetAssayData(csf, assay = "RNA", layer = "data"))
)
message(sprintf(
  "Shared expression: %d genes x %d cells, data layer sum %.6f",
  expression_fingerprint$n_genes, expression_fingerprint$n_cells,
  expression_fingerprint$data_sum
))

# --- Step 3: run marker discovery once per partition -----------------------------
# The parameters are the authors' own: min.pct 0.25 and logfc.threshold 0.25, the
# Wilcoxon default. They are applied identically to both partitions. The seed is reset
# before each call so that neither run inherits state from the other.
marker_counts <- list()
for (partition in partitions) {
  message("\n=== FindAllMarkers on ", partition, " ===")
  Idents(csf) <- factor(csf@meta.data[[partition]],
                        levels = sort(unique(csf@meta.data[[partition]])))
  message("Clusters: ", nlevels(Idents(csf)),
          " (sizes ", paste(range(table(Idents(csf))), collapse = "-"), ")")

  set.seed(seed)
  markers <- FindAllMarkers(
    csf,
    assay = "RNA",
    slot = "data",
    only.pos = FALSE,
    min.pct = min_pct,
    logfc.threshold = logfc_threshold,
    verbose = TRUE
  )
  markers$partition <- partition

  out <- file.path(output_dir, paste0("markers_", partition, ".csv.gz"))
  utils::write.csv(markers, gzfile(out), row.names = FALSE)
  message("Wrote ", nrow(markers), " marker rows to ", basename(out))
  marker_counts[[partition]] <- nrow(markers)

  # Confirm the assay was not disturbed by this run before the next one uses it.
  if (sum(GetAssayData(csf, assay = "RNA", layer = "data")) != expression_fingerprint$data_sum) {
    stop("The expression changed during the run on ", partition)
  }
}

# --- Step 4: the paper's marker panel, per cluster, for both partitions -----------
# Mean log-normalized expression and the fraction of cells expressing each panel gene,
# which is what a dot plot draws. Written as a table so Step 4 of the comparison can
# be read as numbers and not only as a picture.
missing_panel <- setdiff(marker_panel, rownames(csf))
if (length(missing_panel) > 0) {
  message("Panel genes absent from the matrix: ", paste(missing_panel, collapse = ", "))
}
panel_present <- intersect(marker_panel, rownames(csf))
panel_data <- GetAssayData(csf, assay = "RNA", layer = "data")[panel_present, , drop = FALSE]

panel_rows <- list()
for (partition in partitions) {
  labels <- csf@meta.data[[partition]]
  for (cluster in sort(unique(labels))) {
    in_cluster <- labels == cluster
    block <- panel_data[, in_cluster, drop = FALSE]
    panel_rows[[length(panel_rows) + 1]] <- data.frame(
      partition = partition,
      cluster = cluster,
      n_cells = sum(in_cluster),
      gene = panel_present,
      mean_expression = as.numeric(Matrix::rowMeans(block)),
      pct_expressing = as.numeric(Matrix::rowSums(block > 0)) / sum(in_cluster),
      stringsAsFactors = FALSE
    )
  }
}
panel_table <- do.call(rbind, panel_rows)
panel_path <- file.path(output_dir, "canonical_panel_by_cluster.csv")
utils::write.csv(panel_table, panel_path, row.names = FALSE)
message("\nWrote the canonical panel table to ", basename(panel_path))

# --- Step 5: record what ran, so the two runs can be shown to have matched --------
audit <- list(
  stage = "13_marker_discovery",
  seed = seed,
  assay = "RNA",
  slot = "data",
  normalization = "LogNormalize, scale.factor 10000",
  min_pct = min_pct,
  logfc_threshold = logfc_threshold,
  test_used = "wilcox",
  wilcox_implementation = "presto",
  presto_version = as.character(packageVersion("presto")),
  partitions = partitions,
  marker_rows = marker_counts,
  expression_fingerprint = expression_fingerprint,
  panel_genes_used = panel_present,
  panel_genes_missing = missing_panel,
  seurat_version = as.character(packageVersion("Seurat")),
  r_version = paste(R.version$major, R.version$minor, sep = ".")
)
writeLines(
  jsonlite::toJSON(audit, auto_unbox = TRUE, pretty = TRUE, digits = 10),
  file.path(output_dir, "marker_discovery_audit.json")
)

message("\nStage 13 complete.")
