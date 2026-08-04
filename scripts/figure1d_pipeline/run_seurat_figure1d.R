#!/usr/bin/env Rscript

# Faithful, resumable implementation of maintained pipeline Stage 03.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(
    "Usage: run_seurat_figure1d.R <seurat_input_dir> <output_dir> <seed> <force>"
  )
}

input_dir <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- normalizePath(args[[2]], mustWork = FALSE)
random_seed <- as.integer(args[[3]])
force <- tolower(args[[4]]) %in% c("true", "1", "yes")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
final_path <- file.path(output_dir, "integrated_clustered_umap.rds")
if (file.exists(final_path) && !force) {
  message("Using existing integrated result: ", final_path)
  quit(save = "no", status = 0)
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})
if (requireNamespace("future", quietly = TRUE)) {
  future::plan("sequential")
}
options(future.globals.maxSize = 60 * 1024^3)
set.seed(random_seed)

build_integrated_object <- function() {
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
if (!"internal_sample_id" %in% colnames(metadata)) {
  stop("metadata.csv is missing internal_sample_id; rerun stage 02")
}
if (length(unique(metadata$internal_sample_id)) != 59) {
  stop("Expected 59 internal sample IDs")
}

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

message("Splitting object into 59 biological samples...")
csf_list <- SplitObject(csf, split.by = "internal_sample_id")
if (length(csf_list) != 59) {
  stop("SplitObject did not produce 59 sample objects")
}
if (sum(vapply(csf_list, ncol, numeric(1))) != ncol(csf)) {
  stop("Cell count changed while splitting samples")
}

sct_dir <- file.path(output_dir, "sct_samples")
dir.create(sct_dir, recursive = TRUE, showWarnings = FALSE)
message("Running or resuming per-sample SCTransform...")
for (sample_name in names(csf_list)) {
  checkpoint <- file.path(sct_dir, paste0(sample_name, "_sct.rds"))
  if (file.exists(checkpoint) && !force) {
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

message("Selecting integration features and preparing SCT residuals...")
integration_features <- SelectIntegrationFeatures(
  object.list = sct_list,
  nfeatures = 2000
)
sct_list <- PrepSCTIntegration(
  object.list = sct_list,
  anchor.features = integration_features,
  verbose = TRUE
)
reference_indices <- match(c("A4", "A6"), names(sct_list))
if (anyNA(reference_indices)) {
  stop("Could not locate the A4 and A6 integration references")
}

message("Finding SCT integration anchors...")
integration_anchors <- FindIntegrationAnchors(
  object.list = sct_list,
  normalization.method = "SCT",
  anchor.features = integration_features,
  reference = reference_indices,
  dims = 1:30,
  verbose = TRUE
)

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

pca_checkpoint <- file.path(output_dir, "GSE200164_after_PCA.rds")
integrated_checkpoint <- file.path(
  output_dir,
  "GSE200164_seurat_integrated.rds"
)

if (file.exists(pca_checkpoint) && !force) {
  message("Resuming from post-PCA checkpoint: ", pca_checkpoint)
  integrated <- readRDS(pca_checkpoint)
  if (!"pca" %in% Reductions(integrated)) {
    stop("Post-PCA checkpoint does not contain a PCA reduction")
  }
} else {
  if (file.exists(integrated_checkpoint) && !force) {
    message("Resuming from integrated checkpoint: ", integrated_checkpoint)
    integrated <- readRDS(integrated_checkpoint)
  } else {
    integrated <- build_integrated_object()
  }
  if (!"integrated" %in% Assays(integrated)) {
    stop("Resumed object does not contain the integrated assay")
  }
  DefaultAssay(integrated) <- "integrated"
  integrated <- RunPCA(
    integrated,
    assay = "integrated",
    npcs = 50,
    verbose = TRUE
  )
  saveRDS(integrated, pca_checkpoint, compress = FALSE)
}

message("Running neighbors, Louvain clustering, and UMAP...")
DefaultAssay(integrated) <- "integrated"
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
umap <- Embeddings(integrated, reduction = "umap")
if (nrow(umap) != ncol(integrated) || ncol(umap) != 2 || anyNA(umap)) {
  stop("UMAP result failed validation")
}

message("Writing final checkpoint atomically...")
temporary_path <- paste0(final_path, ".tmp")
if (file.exists(temporary_path)) {
  unlink(temporary_path)
}
saveRDS(integrated, temporary_path, compress = "gzip")
if (file.exists(final_path)) {
  if (!force) {
    stop("Final checkpoint exists and force is false")
  }
  unlink(final_path)
}
if (!file.rename(temporary_path, final_path)) {
  stop("Could not rename temporary checkpoint to final path")
}

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
plot <- DimPlot(
  integrated,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  raster = TRUE,
  pt.size = 0.15
)
ggsave(
  filename = file.path(output_dir, "figure1d_unsupervised_clusters.png"),
  plot = plot,
  width = 9,
  height = 7,
  dpi = 300
)

message("Figure 1D workflow complete: ", final_path)
