#!/usr/bin/env Rscript

# Figure 1D pipeline: make a UMAP plot of CSF immune cells.
#
# The analysis follows this order:
# 1. Load R packages and find the project folders.
# 2. Create/reuse a sparse version of the large input count table.
# 3. Perform single-cell quality control, normalization, clustering, and UMAP.
# 4. Label clusters from known marker genes.
# 5. Save the plot and analysis outputs.

# Load the packages used by later steps. Hiding startup messages makes the
# script's own progress messages easier to notice.
suppressPackageStartupMessages({
  # Matrix provides compact sparse matrices for mostly-zero count data.
  library(Matrix)
  # Seurat provides the single-cell analysis functions.
  library(Seurat)
  # ggplot2 formats and saves the UMAP plot.
  library(ggplot2)
  # dplyr summarizes marker scores by cluster.
  library(dplyr)
  # readr reads and writes CSV files.
  library(readr)
})

# Support both `Rscript file.R` and interactive `source("file.R")` execution.
# `--file` is only present in the former, while source() records the path in
# the calling frame's `ofile` field.
script_path <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[1])))
  }

  frame_files <- vapply(sys.frames(), function(frame) {
    if (!is.null(frame$ofile)) {
      return(frame$ofile)
    }
    NA_character_
  }, character(1))
  frame_files <- frame_files[!is.na(frame_files) & nzchar(frame_files)]
  if (length(frame_files) > 0) {
    return(normalizePath(frame_files[length(frame_files)]))
  }

  stop("Could not determine script path. Run with Rscript or source() from a file path.")
})
script_dir <- dirname(script_path)
# This script is in code/test_code/, so going up two folders reaches the root
# folder that contains data/ and results/.
project_root <- normalizePath(file.path(script_dir, "../.."))

# Read command-line options. An option can override a default without editing
# the script, for example: Rscript generate_fig1d_umap.R --resolution=0.5.
args <- commandArgs(trailingOnly = TRUE)
# A flag has no value, such as --rebuild-sparse.
has_arg <- function(flag) flag %in% args
# A value option looks like --name=value. Return its default if it is absent.
arg_value <- function(flag, default = NULL) {
  match <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(match) == 0) return(default)
  sub(paste0("^", flag, "="), "", match[1])
}

# Set the paths for source data, the cache, and final outputs. The defaults use
# the standard project layout, but each location can be overridden by an option.
counts_csv <- arg_value("--counts", file.path(project_root, "data/ad_datasets/GSE200164_counts.csv"))
sample_metadata_csv <- arg_value("--metadata", file.path(project_root, "data/ad_datasets/GSE200164_metadata.csv"))
output_dir <- arg_value("--out-dir", file.path(project_root, "results/figure1d"))
sparse_dir <- arg_value("--sparse-dir", file.path(output_dir, "sparse_cache"))
# The Python helper converts the wide CSV to a compact cache before R reads it.
converter <- file.path(script_dir, "prepare_gse200164_sparse.py")
# Resolution controls how finely Seurat divides cells into clusters.
resolution <- as.numeric(arg_value("--resolution", "0.3"))
# Convert a text range such as "1:15" into the vector 1, 2, ..., 15.
dims <- as.integer(strsplit(arg_value("--dims", "1:15"), ":", fixed = TRUE)[[1]])
dims <- seq(dims[1], dims[2])
# max_cells is optional. It enables a smaller, quicker exploratory run.
max_cells <- arg_value("--max-cells", NULL)
max_cells <- if (is.null(max_cells)) NULL else as.integer(max_cells)

# Create destination folders if this is the first run.
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(sparse_dir, showWarnings = FALSE, recursive = TRUE)

# The Python conversion is intentionally cached because repeatedly parsing the
# wide GEO CSV dominates the runtime of this figure-generation workflow.
required_sparse_files <- file.path(
  sparse_dir,
  c("GSE200164_counts.mtx", "GSE200164_genes.txt", "GSE200164_barcodes.txt", "GSE200164_cell_metadata.csv")
)

if (!all(file.exists(required_sparse_files)) || has_arg("--rebuild-sparse")) {
  message("Building sparse cache from the wide count CSV. This is the slow step, but it is reusable.")
  # Assemble the arguments passed to the converter script.
  converter_args <- c(
    converter,
    "--counts", normalizePath(counts_csv),
    "--metadata", normalizePath(sample_metadata_csv),
    "--out-dir", normalizePath(sparse_dir, mustWork = FALSE)
  )
  if (has_arg("--rebuild-sparse")) {
    converter_args <- c(converter_args, "--force")
  }
  # Run Python. Exit status 0 means the conversion completed successfully.
  status <- system2("python3", converter_args)
  if (!identical(status, 0L)) {
    stop("Sparse conversion failed.")
  }
}

# Step 1: read the sparse cache. The .mtx file contains counts; the two text
# files provide the matching gene names (rows) and cell barcodes (columns).
message("Reading sparse count matrix.")
counts <- readMM(file.path(sparse_dir, "GSE200164_counts.mtx"))
# Convert to the sparse class that Seurat expects.
counts <- as(counts, "dgCMatrix")
genes <- readLines(file.path(sparse_dir, "GSE200164_genes.txt"))
barcodes <- readLines(file.path(sparse_dir, "GSE200164_barcodes.txt"))
rownames(counts) <- genes
colnames(counts) <- barcodes

# Read per-cell sample metadata generated by the converter.
cell_metadata <- read_csv(file.path(sparse_dir, "GSE200164_cell_metadata.csv"), show_col_types = FALSE)
cell_metadata <- as.data.frame(cell_metadata)
rownames(cell_metadata) <- cell_metadata$cell
# Preserve the matrix column order so Seurat attaches each cell's metadata to
# the correct barcode.
cell_metadata <- cell_metadata[colnames(counts), , drop = FALSE]

if (!is.null(max_cells) && max_cells < ncol(counts)) {
  # Optional, reproducible downsampling for quicker exploratory runs.
  set.seed(2026)
  selected_cells <- sample(colnames(counts), max_cells)
  counts <- counts[, selected_cells, drop = FALSE]
  cell_metadata <- cell_metadata[selected_cells, , drop = FALSE]
}

# Step 2: combine counts and metadata into Seurat's main analysis object.
# min.cells removes genes seen in fewer than 3 cells; min.features removes cells
# where fewer than 200 genes were detected.
message("Creating Seurat object.")
s <- CreateSeuratObject(
  counts = counts,
  meta.data = cell_metadata,
  project = "GSE200164_CSF_Aging",
  min.cells = 3,
  min.features = 200
)
# The count matrix is now stored inside `s`, so free this extra copy before the
# analysis steps that need more memory. gc() asks R to reclaim unused memory.
rm(counts)
gc()

# Calculate the percentage of counts from mitochondrial genes in each cell.
# High mitochondrial percentages can indicate low-quality or damaged cells.
s[["percent.mt"]] <- PercentageFeatureSet(s, pattern = "^MT-")
# Keep cells with at least 200 detected genes and no more than 20% mitochondrial
# counts. This is the quality-control filtering step.
s <- subset(s, subset = nFeature_RNA >= 200 & percent.mt <= 20)

# Step 3: prepare the data for clustering and calculate a two-dimensional UMAP.
message("Running Seurat normalization, PCA, clustering, and UMAP.")
# Normalize cell counts so cells with different sequencing depths are comparable.
s <- NormalizeData(s, verbose = FALSE)
# Select the 2,000 genes whose variation is most informative across cells.
s <- FindVariableFeatures(s, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
# Center and scale those genes before dimension reduction.
s <- ScaleData(s, features = VariableFeatures(s), verbose = FALSE)
# PCA summarizes the gene-expression variation into a smaller set of axes.
s <- RunPCA(s, features = VariableFeatures(s), npcs = max(dims), verbose = FALSE)
# Build a graph linking cells that are close in the selected PCA dimensions.
s <- FindNeighbors(s, dims = dims, verbose = FALSE)
# FindClusters identifies groups of strongly connected, similar cells.
s <- FindClusters(s, resolution = resolution, verbose = FALSE)
# UMAP projects the selected PCA dimensions into two coordinates for the plot.
s <- RunUMAP(s, dims = dims, n.neighbors = 30, verbose = FALSE)

# Cluster annotations are assigned from the average normalized expression of
# these canonical marker sets, rather than from an external label transfer.
marker_sets <- list(
  "CD4+ T Cells" = c("CD3D", "CD3E", "CD4", "IL7R", "CCR7"),
  "CD8+ T Cells" = c("CD3D", "CD3E", "CD8A", "CD8B", "GZMK", "NKG7"),
  "T Regulatory Cells" = c("FOXP3", "IL2RA", "CTLA4", "TIGIT"),
  "NK Cells" = c("NKG7", "GNLY", "KLRD1", "KLRF1"),
  "DC" = c("FCER1A", "CLEC10A", "CST3", "LILRA4"),
  "CD14+/CD16-/CD68lo Monocytes" = c("LYZ", "LST1", "S100A8", "S100A9", "CD14"),
  "CD14+/CD16+/CD68mid Monocytes" = c("LYZ", "LST1", "FCGR3A", "MS4A7"),
  "CD14+/CD16+/CD68hi Monoctyes" = c("CD68", "LYZ", "LST1", "CTSS"),
  "B Cells" = c("MS4A1", "CD79A", "CD79B", "CD74"),
  "Plasma Cells" = c("JCHAIN", "MZB1", "IGKC", "IGHG1")
)

# Step 4: use marker genes to give each Seurat cluster an interpretable label.
message("Assigning marker-based cluster labels.")
# Retrieve normalized RNA expression. This is used for marker scoring rather
# than raw counts, so cells with different total counts are more comparable.
expr <- GetAssayData(s, assay = "RNA", slot = "data")
# For each marker set, calculate its average expression in every cell. Markers
# that do not exist in the dataset are skipped.
score_matrix <- vapply(marker_sets, function(markers) {
  present <- intersect(markers, rownames(expr))
  if (length(present) == 0) {
    return(rep(NA_real_, ncol(expr)))
  }
  Matrix::colMeans(expr[present, , drop = FALSE])
}, numeric(ncol(expr)))
rownames(score_matrix) <- colnames(expr)

# Convert the score matrix into a table and record each cell's cluster ID.
score_df <- as.data.frame(score_matrix)
score_df$seurat_clusters <- as.character(s$seurat_clusters)
# Score clusters rather than individual cells to avoid fragmented labels in the
# plotted UMAP.
cluster_scores <- score_df %>%
  group_by(seurat_clusters) %>%
  summarize(across(all_of(names(marker_sets)), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

# Start all clusters as Undetermined, then replace that label only when a marker
# set has a positive, non-missing average score.
cluster_labels <- setNames(rep("Undetermined", nrow(cluster_scores)), cluster_scores$seurat_clusters)
for (i in seq_len(nrow(cluster_scores))) {
  # Extract the marker scores for one cluster and find its strongest marker set.
  scores <- as.numeric(cluster_scores[i, names(marker_sets)])
  if (all(is.na(scores)) || max(scores, na.rm = TRUE) <= 0) {
    cluster_labels[[cluster_scores$seurat_clusters[i]]] <- "Undetermined"
  } else {
    cluster_labels[[cluster_scores$seurat_clusters[i]]] <- names(marker_sets)[which.max(scores)]
  }
}

# Copy each cluster label back onto every cell in that cluster for plotting.
s$marker_celltype <- unname(cluster_labels[as.character(s$seurat_clusters)])
s$marker_celltype[is.na(s$marker_celltype)] <- "Undetermined"

# Fix legend order so it stays consistent across repeated runs.
celltype_levels <- c(
  "CD4+ T Cells", "CD8+ T Cells", "CD4+/CD8+ T Cells", "T Regulatory Cells",
  "NK Cells", "DC", "CD14+/CD16-/CD68lo Monocytes", "CD14+/CD16+/CD68mid Monocytes",
  "CD14+/CD16+/CD68hi Monoctyes", "B Cells", "Plasma Cells", "Undetermined"
)
s$marker_celltype <- factor(s$marker_celltype, levels = celltype_levels)

# Provide one color for every legend label, in the same order as above.
color_vector <- c(
  "darkturquoise", "lawngreen", "dodgerblue", "red",
  "lightpink", "hotpink", "gold1", "navy", "violet",
  "darkorange", "darkviolet", "gray"
)
names(color_vector) <- celltype_levels

# Step 5: create the UMAP. Each dot represents one cell; color shows its
# marker-based cell-type label. shuffle draws cells in random order so a large
# group does not always hide smaller groups behind it.
figure <- DimPlot(
  s,
  group.by = "marker_celltype",
  cols = color_vector,
  shuffle = TRUE,
  reduction = "umap",
  pt.size = 0.15,
  label = FALSE
) +
  labs(title = "Figure 1D recreation: marker-annotated CSF immune cell UMAP", color = NULL) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

# Step 6: choose output filenames. PDF preserves vector quality, PNG is a
# high-resolution image, CSV records the marker evidence, and RDS saves the full
# Seurat object so it can be reopened without rerunning the analysis.
pdf_path <- file.path(output_dir, "figure1d_marker_annotated_umap.pdf")
png_path <- file.path(output_dir, "figure1d_marker_annotated_umap.png")
scores_path <- file.path(output_dir, "cluster_marker_scores.csv")
object_path <- file.path(output_dir, "gse200164_fig1d_seurat.rds")

# Write the plot in two formats, then write the score table and Seurat object.
ggsave(pdf_path, figure, width = 10, height = 8)
ggsave(png_path, figure, width = 10, height = 8, dpi = 300)
write_csv(cluster_scores, scores_path)
saveRDS(s, object_path)

message("Wrote:")
message("  ", pdf_path)
message("  ", png_path)
message("  ", scores_path)
message("  ", object_path)
