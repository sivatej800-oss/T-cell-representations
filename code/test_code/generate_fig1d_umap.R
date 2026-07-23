#!/usr/bin/env Rscript

# Reproduce Figure 1D from Piehl et al., Cell 2022, using GSE200164.
#
# This workflow follows the authors' released code and STAR Methods:
#   1. Load their SoupX-corrected, DoubletFinder/QC-filtered GEO matrix.
#   2. Enforce the published feature and <10% mitochondrial QC thresholds.
#   3. Run SCTransform separately for every sample and perform SCT integration.
#   4. Run PCA, Louvain clustering (resolution 0.3), and UMAP (PCs 1:15,
#      30 neighbors) on the integrated assay.
#   5. Annotate broad canonical immune lineages from the authors' cluster map
#      and verify each assignment with expression-based marker scores.
#   6. Recluster T cells and monocytes, as the authors did, and transfer their
#      intended subtype annotations back to the all-cell UMAP.
#   7. Export the data-generated Figure 1D and complete audit diagnostics.

# Use Seurat v3-style assay objects so Seurat 5 behaves like the authors'
# Seurat 4.1.0 workflow during SplitObject and SCT integration.
options(Seurat.object.assay.version = "v3")

# Readable inverse of base R's membership operator for parameter validation.
`%!in%` <- function(x, table) !(x %in% table)

# Load only packages directly used by this reproducible workflow.
suppressPackageStartupMessages({
  # Matrix stores gene-by-cell expression without expanding millions of zeros.
  library(Matrix)
  # Seurat implements the authors' single-cell normalization and embedding steps.
  library(Seurat)
  # ggplot2 writes the publication panel and diagnostic plots.
  library(ggplot2)
  # dplyr creates transparent QC and annotation summary tables.
  library(dplyr)
  # readr reads metadata and writes auditable CSV outputs.
  library(readr)
  # reticulate imports the lossless SciPy sparse cache created by Python.
  library(reticulate)
})

# Keep Seurat's future-backed integration helpers local and deterministic.
future::plan(future::sequential)
# Permit the 59-sample SCT object list to pass future's serialization guard.
# This changes only an in-process object-size limit, not the analysis or memory allocation.
options(future.globals.maxSize = 64 * 1024^3)

# Fix all stochastic Seurat, UMAP, and plotting operations to a documented seed.
DEFAULT_SEED <- 42L
# Match SCTransform's Seurat 4.1.0 default seed from the authors' environment.
SCT_SEED <- 1448145L
# Match FindClusters' Seurat 4.1.0 default Louvain seed from the authors' code.
LOUVAIN_SEED <- 0L

# Recover the script's absolute path under both Rscript and interactive source().
script_path <- local({
  # Rscript supplies the current file as a --file command-line argument.
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  # Use the Rscript path when it is available.
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[1])))
  }
  # source() stores the active file in one of the evaluation frames.
  frame_files <- vapply(sys.frames(), function(frame) {
    # Return the frame's source file or an explicit missing value.
    if (!is.null(frame$ofile)) return(frame$ofile)
    NA_character_
  }, character(1))
  # Discard frames that do not correspond to a sourced file.
  frame_files <- frame_files[!is.na(frame_files) & nzchar(frame_files)]
  # Use the most recently sourced file when running interactively.
  if (length(frame_files) > 0) {
    return(normalizePath(frame_files[length(frame_files)]))
  }
  # Stop rather than guess a project root from the user's working directory.
  stop("Could not determine this script's path. Run it with Rscript or source().")
})

# This file lives in code/test_code, two levels below the project root.
script_dir <- dirname(script_path)
# Normalize the project root once so every default path is independent of cwd.
project_root <- normalizePath(file.path(script_dir, "../.."))

# Collect user-supplied command-line options after R has handled its own flags.
args <- commandArgs(trailingOnly = TRUE)

# Return TRUE when a valueless flag such as --rebuild-sparse was supplied.
has_arg <- function(flag) {
  # Exact matching avoids treating similarly named options as the same flag.
  flag %in% args
}

# Read an option written as --name=value and otherwise return its default.
arg_value <- function(flag, default = NULL) {
  # Escape regex metacharacters because option names are treated literally.
  pattern <- paste0("^", flag, "=")
  # Find every matching command-line token.
  matches <- grep(pattern, args, value = TRUE)
  # Use the supplied default when the option is absent.
  if (length(matches) == 0) return(default)
  # Reject duplicate settings because silently choosing one is not reproducible.
  if (length(matches) > 1) stop("Option supplied more than once: ", flag)
  # Remove the --name= prefix and return only the option value.
  sub(pattern, "", matches[1])
}

# Return the first existing path from a list of compatible project layouts.
first_existing_path <- function(paths) {
  # Test paths in preference order so the modern layout wins when both exist.
  existing <- paths[file.exists(paths)]
  # Return the first match when at least one input is present.
  if (length(existing) > 0) return(existing[1])
  # Return the preferred path to produce a clear downstream missing-file error.
  paths[1]
}

# Parse a dimension interval such as 1:15 into the exact integer PC sequence.
parse_dims <- function(value) {
  # Split on the literal colon expected by the command-line interface.
  endpoints <- strsplit(value, ":", fixed = TRUE)[[1]]
  # Require exactly one start and one end value.
  if (length(endpoints) != 2) stop("--dims must look like 1:15")
  # Convert both endpoints to integers.
  endpoints <- suppressWarnings(as.integer(endpoints))
  # Reject missing, reversed, or non-positive PC intervals.
  if (anyNA(endpoints) || endpoints[1] < 1 || endpoints[2] < endpoints[1]) {
    stop("--dims must be a positive increasing interval such as 1:15")
  }
  # Return every principal-component index in the inclusive interval.
  seq.int(endpoints[1], endpoints[2])
}

# Locate the count file under either the uploaded or documented data folder name.
default_counts <- first_existing_path(c(
  file.path(project_root, "data", "AD datasets", "GSE200164_counts.csv"),
  file.path(project_root, "data", "ad_datasets", "GSE200164_counts.csv")
))
# Locate matching GEO-derived sample metadata under the same compatible layouts.
default_metadata <- first_existing_path(c(
  file.path(project_root, "data", "AD datasets", "GSE200164_metadata.csv"),
  file.path(project_root, "data", "ad_datasets", "GSE200164_metadata.csv")
))

# Allow explicit input paths while making the uploaded data work without edits.
counts_csv <- arg_value("--counts", default_counts)
# Use sample metadata to split integration and retain diagnosis/age/sex covariates.
sample_metadata_csv <- arg_value("--metadata", default_metadata)
# Keep all generated artifacts in one ignored, reproducible results directory.
output_dir <- arg_value("--out-dir", file.path(project_root, "results", "figure1d"))
# Reuse a lossless sparse cache because reparsing the 4+ GB CSV is expensive.
sparse_dir <- arg_value("--sparse-dir", file.path(output_dir, "sparse_cache"))
# Allow an exact continuation from a validated sparse cache when the raw GEO CSV
# remains archived elsewhere and is intentionally absent from the compute VM.
cache_only <- has_arg("--cache-only")
# Point to the streaming Python converter maintained beside this script.
converter <- file.path(script_dir, "prepare_gse200164_sparse.py")
# Match the authors' Louvain resolution unless explicitly overridden.
resolution <- as.numeric(arg_value("--resolution", "0.3"))
# Match the authors' first 15 principal components unless explicitly overridden.
dims <- parse_dims(arg_value("--dims", "1:15"))
# Use the fixed default seed unless a valid integer override is supplied.
seed <- as.integer(arg_value("--seed", as.character(DEFAULT_SEED)))
# Optionally retain at most this many cells per sample for a smoke test.
max_cells_per_sample <- arg_value("--max-cells-per-sample", NULL)
# Convert the optional sampling limit only when the user supplied it.
max_cells_per_sample <- if (is.null(max_cells_per_sample)) {
  NULL
} else {
  as.integer(max_cells_per_sample)
}
# Use the authors' default 30 integration dimensions for the full recreation.
# A 50-cell/sample smoke test uses 15 to avoid an ill-conditioned truncated SVD.
integration_dims <- if (is.null(max_cells_per_sample)) 1:30 else 1:15
# Keep the authors' default integration weighting for the full data.
# The deliberately tiny smoke samples need fewer anchor neighbors than 100.
integration_k_weight <- if (is.null(max_cells_per_sample)) 100L else 20L

# Reject malformed numeric parameters before creating any partial output.
if (is.na(resolution) || resolution <= 0) stop("--resolution must be positive")
# Reject a missing or negative random seed.
if (is.na(seed) || seed < 0) stop("--seed must be a non-negative integer")
# Reject sampling limits too small to support the default 30 integration PCs.
if (!is.null(max_cells_per_sample) &&
    (is.na(max_cells_per_sample) || max_cells_per_sample < 50)) {
  stop("--max-cells-per-sample must be at least 50")
}

# Cache-only execution still needs the converter script for provenance, but it
# deliberately does not require raw counts or source metadata on the VM.
required_input_paths <- if (cache_only) converter else c(
  counts_csv, sample_metadata_csv, converter
)
# Require the inputs needed for the selected execution mode before integration.
for (path in required_input_paths) {
  # Name the exact missing path in the error shown to the user.
  if (!file.exists(path)) stop("Required input does not exist: ", path)
}

# Create result and cache directories recursively on the first run.
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
# Create the sparse cache independently when the user points it elsewhere.
dir.create(sparse_dir, showWarnings = FALSE, recursive = TRUE)
# Keep diagnostic plots separate from the final requested panel.
diagnostics_dir <- file.path(output_dir, "diagnostics")
# Create the diagnostic directory before any plot-writing functions run.
dir.create(diagnostics_dir, showWarnings = FALSE, recursive = TRUE)

# Define every cache file required for a complete lossless sparse import.
required_sparse_files <- file.path(sparse_dir, c(
  "GSE200164_counts_csr.npz",
  "GSE200164_genes.txt",
  "GSE200164_barcodes.txt",
  "GSE200164_cell_metadata.csv",
  "GSE200164_sparse_manifest.json",
  ".sparse_cache_complete"
))

# A cache-only run may never regenerate a cache because doing so requires the
# deliberately omitted raw source files.
if (cache_only && has_arg("--rebuild-sparse")) {
  stop("--cache-only cannot be combined with --rebuild-sparse")
}

# Stop clearly when a VM transfer omitted any cache artifact needed to resume.
if (cache_only && !all(file.exists(required_sparse_files))) {
  stop(
    "--cache-only requires a complete sparse cache; missing: ",
    paste(basename(required_sparse_files[!file.exists(required_sparse_files)]),
          collapse = ", ")
  )
}

# Rebuild when any component is missing or when the user requests fresh input.
if (!cache_only && (!all(file.exists(required_sparse_files)) ||
                    has_arg("--rebuild-sparse"))) {
  # Tell the user why this one-time conversion can take several minutes.
  message("Building a lossless sparse cache from the 4+ GB GEO count table.")
  # Quote path arguments because this workspace and the uploaded data folders contain spaces.
  converter_args <- c(
    shQuote(converter),
    "--counts", shQuote(normalizePath(counts_csv)),
    "--metadata", shQuote(normalizePath(sample_metadata_csv)),
    "--out-dir", shQuote(normalizePath(sparse_dir, mustWork = FALSE))
  )
  # Forward the explicit rebuild request to the cache validator.
  if (has_arg("--rebuild-sparse")) converter_args <- c(converter_args, "--force")
  # Run the converter with the same Python executable discoverable by R.
  conversion_status <- system2("python3", converter_args)
  # Any nonzero process status means the matrix cache is unsafe to use.
  if (!identical(conversion_status, 0L)) stop("Sparse conversion failed")
}

# Recheck all outputs because an interrupted converter must never be accepted.
if (!all(file.exists(required_sparse_files))) {
  # List missing pieces to distinguish a partial cache from a path typo.
  stop(
    "Sparse cache is incomplete; missing: ",
    paste(basename(required_sparse_files[!file.exists(required_sparse_files)]),
          collapse = ", ")
  )
}

# Read gene symbols exactly in sparse row order.
genes <- readLines(file.path(sparse_dir, "GSE200164_genes.txt"))
# Read cell barcodes exactly in sparse column order.
barcodes <- readLines(file.path(sparse_dir, "GSE200164_barcodes.txt"))
# Read the prevalidated one-row-per-cell metadata table.
cell_metadata <- read_csv(
  file.path(sparse_dir, "GSE200164_cell_metadata.csv"),
  show_col_types = FALSE
)
# Convert tibble semantics to the base data.frame expected by Seurat metadata.
cell_metadata <- as.data.frame(cell_metadata)
# Require a unique row key before aligning metadata to matrix columns.
if (anyDuplicated(cell_metadata$cell)) stop("Cell metadata contains duplicate barcodes")
# Use cell barcodes as row names for Seurat's alignment checks.
rownames(cell_metadata) <- cell_metadata$cell
# Reorder metadata explicitly to match every sparse matrix column.
cell_metadata <- cell_metadata[barcodes, , drop = FALSE]
# Stop on missing joins instead of letting Seurat silently insert NA metadata.
if (anyNA(cell_metadata$sample_id)) stop("Cell metadata did not align to all barcodes")

# Force reticulate to use the same system Python that provides SciPy and pandas.
Sys.setenv(RETICULATE_PYTHON = Sys.which("python3"))
# Import SciPy without automatic conversion so loading stays sparse in Python.
scipy_sparse <- reticulate::import("scipy.sparse", convert = FALSE)
# Load the standard CSR archive without constructing a dense gene-by-cell array.
counts_python <- scipy_sparse$load_npz(
  normalizePath(file.path(sparse_dir, "GSE200164_counts_csr.npz"))
)
# Convert SciPy CSR directly to a Matrix sparse row matrix in R.
counts <- reticulate::py_to_r(counts_python)
# Release the Python proxy before later Seurat stages reach peak memory.
rm(counts_python)
# Convert once to the column-compressed representation Seurat operates on best.
# Matrix 1.7 exposes the safe dgRMatrix -> CsparseMatrix coercion via this parent class.
counts <- as(counts, "CsparseMatrix")
# Verify independent name vectors against the serialized sparse dimensions.
if (!identical(dim(counts), c(length(genes), length(barcodes)))) {
  stop("Sparse matrix dimensions do not match gene and barcode files")
}
# Attach feature names to sparse rows without altering their order.
rownames(counts) <- genes
# Attach cell barcodes to sparse columns without altering their order.
colnames(counts) <- barcodes

# Optionally downsample each sample independently for a representative smoke test.
if (!is.null(max_cells_per_sample)) {
  # Fix sampling so repeated test runs retain the same cells.
  set.seed(seed)
  # Split barcodes by sample while preserving the GEO metadata sample order.
  cells_by_sample <- split(
    barcodes,
    factor(cell_metadata[barcodes, "sample_id"],
           levels = unique(cell_metadata[barcodes, "sample_id"]))
  )
  # Sample no more than the requested number from every biological sample.
  selected_cells <- unlist(lapply(cells_by_sample, function(sample_cells) {
    # Keep all cells from small samples and randomly sample only large samples.
    if (length(sample_cells) <= max_cells_per_sample) return(sample_cells)
    sample(sample_cells, max_cells_per_sample)
  }), use.names = FALSE)
  # Subset matrix columns to the selected test cells.
  counts <- counts[, selected_cells, drop = FALSE]
  # Subset and reorder metadata to the exact same selected cell sequence.
  cell_metadata <- cell_metadata[selected_cells, , drop = FALSE]
  # Report that this is a method smoke test rather than the final full recreation.
  message("Smoke-test subset retained ", length(selected_cells), " cells.")
}

# Build a parameter record used to detect incompatible checkpoint reuse.
# Do not resolve absent raw-data paths on a cache-only VM.  The existing
# checkpoint and cache manifest retain the original source provenance instead.
source_parameter_values <- if (cache_only) {
  c("<cache-only>", NA_character_, NA_character_, "<cache-only>")
} else {
  c(
    normalizePath(counts_csv),
    as.character(file.info(counts_csv)$size),
    as.character(file.info(counts_csv)$mtime),
    normalizePath(sample_metadata_csv)
  )
}

# Build the run settings that determine whether a checkpoint is compatible.
analysis_parameters <- data.frame(
  parameter = c(
    "counts_path", "counts_size", "counts_mtime", "metadata_path",
    "resolution", "dims", "integration_dims", "integration_k_weight",
    "pca_umap_seed", "sct_seed", "louvain_seed",
    "max_cells_per_sample", "seurat_version"
  ),
  value = c(
    source_parameter_values,
    as.character(resolution),
    paste(range(dims), collapse = ":"),
    paste(range(integration_dims), collapse = ":"),
    as.character(integration_k_weight),
    as.character(seed),
    as.character(SCT_SEED),
    as.character(LOUVAIN_SEED),
    if (is.null(max_cells_per_sample)) "all" else as.character(max_cells_per_sample),
    as.character(packageVersion("Seurat"))
  ),
  stringsAsFactors = FALSE
)
# Store parameters beside checkpoints so future resumes are auditable.
parameter_path <- file.path(output_dir, "analysis_parameters.csv")
# A parameter change invalidates analysis checkpoints but not the source cache.
if (file.exists(parameter_path) && !has_arg("--rebuild-analysis")) {
  # Read the exact parameter table used by the previous analysis run.
  previous_parameters <- read_csv(parameter_path, show_col_types = FALSE)
  # A cache-only VM has a different filesystem path for archived raw inputs;
  # compare biological and computational settings while retaining the original
  # raw-source provenance in the cache manifest and existing parameter record.
  source_parameter_names <- c(
    "counts_path", "counts_size", "counts_mtime", "metadata_path"
  )
  compared_parameters <- if (cache_only) {
    analysis_parameters$parameter %!in% source_parameter_names
  } else {
    rep(TRUE, nrow(analysis_parameters))
  }
  # Compare values after normalizing both inputs to base data frames.
  # CSV reads receive fresh sequential row names, whereas subsetting the
  # in-memory table retains its original indices.  Compare the actual fields,
  # not this incidental row-name attribute.
  previous_compared <- as.data.frame(
    previous_parameters[compared_parameters, , drop = FALSE]
  )
  current_compared <- analysis_parameters[compared_parameters, , drop = FALSE]
  rownames(previous_compared) <- NULL
  rownames(current_compared) <- NULL
  if (!identical(previous_compared, current_compared)) {
    stop(
      "Analysis parameters changed. Re-run with --rebuild-analysis to replace ",
      "incompatible checkpoints."
    )
  }
}
# Preserve the original raw-source parameter record during a cache-only resume;
# it remains the audit trail linking this cache to the archived GEO files.
if (!cache_only || !file.exists(parameter_path)) {
  write_csv(analysis_parameters, parameter_path)
}

# Name checkpoints for the QC, integrated, and final annotated stages.
qc_checkpoint <- file.path(output_dir, "checkpoint_qc_seurat.rds")
# Save anchors separately so an IntegrateData compatibility error is inexpensive to retry.
anchors_checkpoint <- file.path(output_dir, "checkpoint_sct_anchors.rds")
# Save integration separately because it is the most expensive reproducible step.
integrated_checkpoint <- file.path(output_dir, "checkpoint_sct_integrated.rds")
# Save the final object after marker-verified hierarchical annotation.
final_object_path <- file.path(output_dir, "gse200164_figure1d_seurat.rds")

# Rebuild the QC object unless a compatible checkpoint is available.
if (file.exists(qc_checkpoint) && !has_arg("--rebuild-analysis")) {
  # Resume from the validated QC checkpoint to avoid repeated object creation.
  message("Loading QC checkpoint: ", qc_checkpoint)
  s <- readRDS(qc_checkpoint)
  # Release the duplicate raw sparse matrix after loading the saved Seurat object.
  rm(counts)
  # Ask R to reclaim the released matrix memory before SCT integration.
  invisible(gc())
} else {
  # Build Seurat's central analysis object from the lossless fractional counts.
  message("Creating Seurat object from SoupX-adjusted counts.")
  s <- CreateSeuratObject(
    counts = counts,
    meta.data = cell_metadata,
    project = "GSE200164_CSF_Aging",
    min.cells = 3,
    min.features = 200
  )
  # Reassert original sample identities in the GEO sample order.
  s$orig.ident <- factor(
    s$sample_id,
    levels = unique(cell_metadata$sample_id)
  )
  # Calculate mitochondrial abundance using human HGNC mitochondrial prefixes.
  s[["percent.mt"]] <- PercentageFeatureSet(s, pattern = "^MT-")
  # Record the number of cells before applying the paper's mitochondrial cutoff.
  cells_before_mt_qc <- ncol(s)
  # Match the authors' strict threshold: retain percent.mt below, not equal to, 10.
  s <- subset(s, subset = nFeature_RNA >= 200 & percent.mt < 10)
  # Report whether the published GEO matrix was already fully QC filtered.
  message(
    "Published QC retained ", ncol(s), " of ", cells_before_mt_qc,
    " input cells (<10% mitochondrial reads)."
  )
  # Normalize the RNA assay only for annotation evidence and diagnostic plots.
  # Clustering below still uses the authors' independently SCTransformed assays.
  s <- NormalizeData(s, assay = "RNA", normalization.method = "LogNormalize",
                     scale.factor = 10000, verbose = FALSE)
  # Release the duplicate input matrix before SCT integration peaks in memory.
  rm(counts)
  # Ask R to reclaim released matrix blocks before the next stage.
  invisible(gc())
  # Save an uncompressed local checkpoint quickly; it remains a sparse R object.
  saveRDS(s, qc_checkpoint, compress = FALSE)
}

# Summarize per-sample cell retention as a formal QC artifact.
qc_summary <- s@meta.data %>%
  # Convert row names to an explicit cell column for portable CSV output.
  mutate(cell = rownames(s@meta.data)) %>%
  # Group QC statistics at the biological-sample level.
  group_by(sample_id, diagnosis, sex, age) %>%
  # Report retained cells and the distributions that governed filtering.
  summarize(
    retained_cells = n(),
    median_features = median(nFeature_RNA),
    median_counts = median(nCount_RNA),
    median_percent_mt = median(percent.mt),
    maximum_percent_mt = max(percent.mt),
    .groups = "drop"
  )
# Write QC evidence before performing integration.
write_csv(qc_summary, file.path(output_dir, "qc_summary_by_sample.csv"))

# Resume from SCT integration when a compatible checkpoint already exists.
integrated_checkpoint_loaded <- FALSE
if (file.exists(integrated_checkpoint) && !has_arg("--rebuild-analysis")) {
  # An interrupted serialization can leave a truncated checkpoint.  Retain the
  # valid anchor checkpoint and recompute integration instead of blocking resume.
  message("Loading SCT integration checkpoint: ", integrated_checkpoint)
  s <- tryCatch(
    readRDS(integrated_checkpoint),
    error = function(error) {
      warning(
        "Ignoring unreadable SCT integration checkpoint: ", conditionMessage(error),
        ". Recomputing integration from SCT anchors."
      )
      NULL
    }
  )
  integrated_checkpoint_loaded <- !is.null(s)
}
if (!integrated_checkpoint_loaded) {
  # Reuse a completed anchor set when only a later Seurat compatibility step failed.
  if (file.exists(anchors_checkpoint) && !has_arg("--rebuild-analysis")) {
    # Loading the anchors also restores the prepared per-sample object list they contain.
    message("Loading SCT anchor checkpoint: ", anchors_checkpoint)
    anchors <- readRDS(anchors_checkpoint)
  } else {
    # Split cells so SCTransform models every biological sample independently.
    message("Running per-sample SCTransform for ", length(unique(s$orig.ident)), " samples.")
    s_split <- SplitObject(s, split.by = "orig.ident")
    # Process samples serially to keep memory use predictable and reproducible.
    for (sample_name in names(s_split)) {
      # Make progress visible during the longest normalization stage.
      message("  SCTransform sample ", sample_name, " (", ncol(s_split[[sample_name]]), " cells)")
      # Match Seurat 4.1/sctransform 0.3.3 behavior with the original v1 model.
      s_split[[sample_name]] <- SCTransform(
        s_split[[sample_name]],
        assay = "RNA",
        new.assay.name = "SCT",
        variable.features.n = 1000,
        ncells = min(10000, ncol(s_split[[sample_name]])),
        vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
        vst.flavor = "v1",
        seed.use = SCT_SEED,
        verbose = FALSE
      )
      # Make the sample's SCT assay explicit for integration feature selection.
      DefaultAssay(s_split[[sample_name]]) <- "SCT"
    }
    # Select shared integration features exactly as in the authors' released code.
    message("Selecting SCT integration features.")
    integration_features <- SelectIntegrationFeatures(
      object.list = s_split,
      assay = rep("SCT", length(s_split))
    )
    # Calculate residuals for the selected features across every sample model.
    s_split <- PrepSCTIntegration(
      object.list = s_split,
      anchor.features = integration_features,
      verbose = FALSE
    )
    # Locate the authors' stated A4 healthy and A6 impaired reference samples.
    reference_indices <- match(c("A4", "A6"), names(s_split))
    # Fall back to all samples only if a smoke-test input omitted a reference.
    reference_indices <- reference_indices[!is.na(reference_indices)]
    # Find cross-sample anchors while retaining the authors' memory-saving references.
    message("Finding SCT integration anchors.")
    anchors <- FindIntegrationAnchors(
      object.list = s_split,
      normalization.method = "SCT",
      anchor.features = integration_features,
      reference = if (length(reference_indices) > 0) reference_indices else NULL,
      dims = integration_dims,
      verbose = FALSE
    )
    # Checkpoint the expensive anchor set before beginning version-sensitive integration.
    saveRDS(anchors, anchors_checkpoint, compress = FALSE)
    # Release the separate split list; the anchor set retains what IntegrateData needs.
    rm(s_split)
    # Reclaim the duplicate list before integrated matrix construction.
    invisible(gc())
  }
  # Harmonize batch/sample effects into the integrated expression assay.
  message("Integrating all samples on the SCT anchors.")
  s <- IntegrateData(
    anchorset = anchors,
    normalization.method = "SCT",
    dims = integration_dims,
    k.weight = integration_k_weight,
    verbose = FALSE
  )
  # Remove anchors before PCA and UMAP allocate more memory.
  rm(anchors)
  # Reclaim their memory before serializing the integrated checkpoint.
  invisible(gc())
  # Serialize to a temporary file, then replace the checkpoint only when the
  # write succeeded, so an interruption cannot invalidate a completed resume.
  integrated_checkpoint_tmp <- paste0(integrated_checkpoint, ".tmp")
  saveRDS(s, integrated_checkpoint_tmp, compress = FALSE)
  if (!file.rename(integrated_checkpoint_tmp, integrated_checkpoint)) {
    stop("Could not promote completed SCT integration checkpoint into place.")
  }
}

# Use the harmonized assay for the global Figure 1D embedding and clusters.
DefaultAssay(s) <- "integrated"
# Generate enough PCs to cover the requested interval while matching the authors.
message("Running integrated PCA, Louvain clustering, and UMAP.")
s <- RunPCA(s, npcs = max(50L, max(dims)), seed.use = seed, verbose = FALSE)
# Construct the shared-nearest-neighbor graph from the first 15 PCs by default.
s <- FindNeighbors(s, dims = dims, verbose = FALSE)
# Match the authors' Louvain algorithm and resolution 0.3 with a fixed seed.
s <- FindClusters(
  s,
  resolution = resolution,
  algorithm = 1,
  random.seed = LOUVAIN_SEED,
  verbose = FALSE
)
# Preserve global cluster IDs before compartment-specific reclustering.
s$seurat_clusters_allcells <- as.character(s$seurat_clusters)
# Project the same PCs with the authors' 30-neighbor UMAP setting.
s <- RunUMAP(
  s,
  dims = dims,
  n.neighbors = 30,
  seed.use = seed,
  verbose = FALSE
)

# Define canonical lineage signatures that verify, rather than replace, the
# authors' expression-guided manual cluster annotations.
broad_signatures <- list(
  # Pan-T markers distinguish T cells from NKG7-positive natural killer cells.
  "T Cells" = c("CD3D", "CD3E", "TRAC", "TRBC1", "LTB"),
  # NK cytotoxic genes are evaluated without shared pan-T markers.
  "NK Cells" = c("GNLY", "KLRD1", "KLRF1", "PRF1", "NKG7"),
  # Dendritic markers reproduce the FCER1A evidence shown in Figure 1D/E.
  "DC" = c("FCER1A", "CD1C", "CLEC10A", "CST3"),
  # Myeloid lineage markers encompass both classical and non-classical monocytes.
  "Monocytes" = c("LST1", "TYROBP", "FCER1G", "CTSS", "LILRB1"),
  # B-cell markers include the paper's CD19 plus robust B-lineage transcripts.
  "B Cells" = c("CD19", "MS4A1", "CD79A", "CD74", "CD37"),
  # Antibody-secreting-cell markers include the paper's JCHAIN marker.
  "Plasma Cells" = c("JCHAIN", "MZB1", "SDC1", "IGHG1", "IGKC")
)

# Score cluster-level marker evidence without modifying expression or the UMAP.
score_cluster_signatures <- function(object, group_column, signatures) {
  # Pull sparse log-normalized RNA values used only for cell-type evidence.
  expression <- GetAssayData(object, assay = "RNA", layer = "data")
  # Read group IDs in exact cell order from Seurat metadata.
  groups <- as.character(object[[group_column]][, 1])
  # Sort numeric-looking cluster IDs numerically for readable output.
  unique_groups <- unique(groups)
  # Use numeric ordering only when every group can be parsed as a number.
  numeric_groups <- suppressWarnings(as.numeric(unique_groups))
  # Choose stable numeric or lexical ordering as appropriate.
  group_levels <- if (!anyNA(numeric_groups)) {
    unique_groups[order(numeric_groups)]
  } else {
    sort(unique_groups)
  }
  # Build a sparse cell-by-cluster indicator matrix.
  group_design <- sparse.model.matrix(
    ~ 0 + factor(groups, levels = group_levels)
  )
  # Give indicator columns the original cluster labels.
  colnames(group_design) <- group_levels
  # Restrict calculations to the union of canonical markers present in the data.
  marker_union <- intersect(unique(unlist(signatures)), rownames(expression))
  # Require marker evidence rather than producing arbitrary empty scores.
  if (length(marker_union) == 0) stop("None of the annotation markers are present")
  # Sum sparse expression for every marker in every cluster.
  expression_sums <- expression[marker_union, , drop = FALSE] %*% group_design
  # Divide each cluster column by its number of cells to obtain mean expression.
  expression_means <- sweep(
    as.matrix(expression_sums),
    MARGIN = 2,
    STATS = Matrix::colSums(group_design),
    FUN = "/"
  )
  # Center each marker across clusters so lineage-specific enrichment is comparable.
  marker_z <- t(scale(t(expression_means)))
  # Constant markers provide no cluster evidence and therefore receive score zero.
  marker_z[is.na(marker_z)] <- 0
  # Create one long-form record for every cluster/signature pair.
  score_rows <- lapply(names(signatures), function(signature_name) {
    # Retain only signature markers actually measured in GSE200164.
    present <- intersect(signatures[[signature_name]], rownames(expression_means))
    # Calculate mean z-score evidence across the signature's present markers.
    z_score <- colMeans(marker_z[present, , drop = FALSE])
    # Also report unscaled mean log-normalized expression for interpretation.
    mean_expression <- colMeans(expression_means[present, , drop = FALSE])
    # Return one row per cluster for this biological signature.
    data.frame(
      cluster = group_levels,
      signature = signature_name,
      z_score = unname(z_score),
      mean_expression = unname(mean_expression),
      present_markers = paste(present, collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  # Combine signature-specific tables into one auditable data frame.
  bind_rows(score_rows)
}

# Return the highest marker-signature label for each cluster in a score table.
best_signature_by_cluster <- function(score_table) {
  # Split evidence by cluster and retain the largest standardized score.
  best <- score_table %>%
    group_by(cluster) %>%
    slice_max(order_by = z_score, n = 1, with_ties = FALSE) %>%
    ungroup()
  # Return a named vector convenient for cell-level mapping.
  setNames(best$signature, best$cluster)
}

# Calculate broad-lineage evidence for the global Seurat clusters.
broad_scores <- score_cluster_signatures(
  s, "seurat_clusters_allcells", broad_signatures
)
# Predict the strongest canonical lineage independently for audit comparison.
broad_predictions <- best_signature_by_cluster(broad_scores)
# Encode the broad cluster identities published in the authors' Figure 1D script.
author_broad_map <- c(
  "0" = "T Cells", "1" = "T Cells", "2" = "T Cells",
  "3" = "Monocytes", "4" = "T Cells", "5" = "DC",
  "6" = "T Cells", "7" = "NK Cells", "8" = "Monocytes",
  "9" = "T Cells", "10" = "Plasma Cells", "11" = "B Cells",
  "12" = "Undetermined"
)
# Read the cluster IDs actually produced by this software environment.
observed_global_clusters <- sort(unique(as.character(s$seurat_clusters_allcells)))
# Use the exact author map only when the expected 13-cluster solution was reproduced.
if (setequal(observed_global_clusters, names(author_broad_map))) {
  # Initialize labels from the authors' expression-guided assignments.
  broad_map <- author_broad_map
  # Record marker agreement for every non-undetermined author cluster.
  broad_validation <- data.frame(
    cluster = names(broad_map),
    author_label = unname(broad_map),
    strongest_marker_signature = unname(broad_predictions[names(broad_map)]),
    stringsAsFactors = FALSE
  ) %>%
    # Undetermined has no positive signature and is intentionally exempt.
    mutate(
      marker_agrees = author_label == strongest_marker_signature |
        author_label == "Undetermined"
    )
  # Warn about disagreements while preserving the exact published mapping.
  if (any(!broad_validation$marker_agrees)) {
    warning(
      "Some author broad-cluster labels were not the highest marker score; ",
      "inspect broad_cluster_annotation_validation.csv."
    )
  }
  # Record that labels follow the exact published cluster mapping.
  s$broad_annotation_source <- "authors_cluster_map_marker_checked"
} else {
  # Fall back to reproducible canonical-marker labels when cluster numbering differs.
  warning(
    "Global cluster IDs differ from the authors' 0:12 solution; using the ",
    "strongest canonical lineage signature for each cluster."
  )
  # Use expression-derived predictions as a functionally similar valid annotation.
  broad_map <- broad_predictions
  # Preserve a complete validation table in the fallback case as well.
  broad_validation <- data.frame(
    cluster = names(broad_map),
    author_label = NA_character_,
    strongest_marker_signature = unname(broad_map),
    marker_agrees = NA,
    stringsAsFactors = FALSE
  )
  # Record the fallback provenance on every cell.
  s$broad_annotation_source <- "canonical_marker_score_fallback"
}
# Map each global cluster to one broad immune compartment.
s$cluster_ident <- unname(broad_map[as.character(s$seurat_clusters_allcells)])
# Replace any unmapped value with an explicit non-biological placeholder.
s$cluster_ident[is.na(s$cluster_ident)] <- "Undetermined"
# Write broad marker evidence and author-map validation before subclustering.
write_csv(broad_scores, file.path(output_dir, "broad_cluster_marker_scores.csv"))
# Write the direct comparison that makes manual marker checking reproducible.
write_csv(
  broad_validation,
  file.path(output_dir, "broad_cluster_annotation_validation.csv")
)

# Define a reusable version of the authors' compartment-specific reclustering.
recluster_compartment <- function(object, broad_label, dims, resolution, seed) {
  # Retain only cells assigned to the requested broad immune compartment.
  compartment <- subset(object, subset = cluster_ident == broad_label)
  # Cluster on the harmonized assay, exactly as in the authors' code.
  DefaultAssay(compartment) <- "integrated"
  # Recalculate variable features within the restricted biological compartment.
  compartment <- FindVariableFeatures(
    compartment,
    assay = "integrated",
    selection.method = "vst",
    nfeatures = 2000,
    verbose = FALSE
  )
  # Recompute PCA so T-cell or monocyte variation is not dominated by other lineages.
  compartment <- RunPCA(
    compartment,
    npcs = max(50L, max(dims)),
    seed.use = seed,
    verbose = FALSE
  )
  # Construct a compartment-specific shared-nearest-neighbor graph.
  compartment <- FindNeighbors(compartment, dims = dims, verbose = FALSE)
  # Apply the same Louvain algorithm and resolution used for all cells.
  compartment <- FindClusters(
    compartment,
    resolution = resolution,
    algorithm = 1,
    random.seed = LOUVAIN_SEED,
    verbose = FALSE
  )
  # Generate a diagnostic compartment UMAP using the same 30-neighbor setting.
  compartment <- RunUMAP(
    compartment,
    dims = dims,
    n.neighbors = 30,
    seed.use = seed,
    verbose = FALSE
  )
  # Return the reclustered compartment without modifying the global UMAP.
  compartment
}

# Recluster the broad T-cell compartment to resolve CD4, CD8, DP, and Treg cells.
message("Reclustering the broad T-cell compartment.")
t_cells <- recluster_compartment(s, "T Cells", dims, resolution, seed)
# Preserve subcluster IDs as characters to avoid the authors' factor +1 coercion bug.
t_cluster_ids <- as.character(t_cells$seurat_clusters)
# Add a named T-cell subcluster column aligned by cell barcode.
s$seurat_clusters_tcells <- NA_character_
# Transfer character IDs by names, never by factor integer codes.
s@meta.data[colnames(t_cells), "seurat_clusters_tcells"] <- t_cluster_ids

# Define marker panels for independent T-subcluster annotation diagnostics.
t_signatures <- list(
  # Conventional helper/naive T-cell evidence.
  "CD4+ T Cells" = c("CD4", "IL7R", "CCR7", "LTB", "MAL"),
  # Cytotoxic T-cell evidence separate from NK broad-lineage assignment.
  "CD8+ T Cells" = c("CD8A", "CD8B", "CCL5", "GZMK", "GZMA"),
  # Coexpression evidence used for the authors' double-positive clusters.
  "CD4+/CD8+ T Cells" = c("CD4", "CD8A", "CD8B"),
  # Regulatory T-cell genes include the FOXP3 marker shown by the authors.
  "T Regulatory Cells" = c("FOXP3", "IL2RA", "CTLA4", "TIGIT")
)
# Score T-cell marker evidence on the original log-normalized RNA assay.
t_scores <- score_cluster_signatures(t_cells, "seurat_clusters", t_signatures)
# Calculate a reproducible marker-only fallback for changed cluster solutions.
t_predictions <- best_signature_by_cluster(t_scores)
# Correct the authors' intended 1-based factor-code mapping to actual 0-based IDs.
author_t_map <- c(
  "0" = "CD4+ T Cells",
  "1" = "CD4+ T Cells",
  "2" = "CD8+ T Cells",
  "3" = "CD4+/CD8+ T Cells",
  "4" = "CD4+ T Cells",
  "5" = "T Regulatory Cells",
  "6" = "CD4+/CD8+ T Cells"
)
# Use the corrected exact map when the expected seven T clusters are reproduced.
if (setequal(unique(t_cluster_ids), names(author_t_map))) {
  # Apply the intended author labels with explicit character keys.
  t_map <- author_t_map
  # Record exact-map versus strongest-signature evidence for every T subcluster.
  t_validation <- data.frame(
    cluster = names(t_map),
    author_label = unname(t_map),
    strongest_marker_signature = unname(t_predictions[names(t_map)]),
    marker_agrees = unname(t_map) == unname(t_predictions[names(t_map)]),
    stringsAsFactors = FALSE
  )
} else {
  # Use canonical expression signatures if software changes cluster numbering/count.
  warning("T-cell subclusters differ from the authors' 0:6 solution; using marker scores.")
  # Retain the functionally similar marker-derived subtype assignments.
  t_map <- t_predictions
  # Record fallback assignments with no inapplicable author ID comparison.
  t_validation <- data.frame(
    cluster = names(t_map), author_label = NA_character_,
    strongest_marker_signature = unname(t_map), marker_agrees = NA,
    stringsAsFactors = FALSE
  )
}
# Transfer final T subtype labels to the corresponding all-cell metadata rows.
s@meta.data[colnames(t_cells), "cluster_ident"] <- unname(t_map[t_cluster_ids])
# Write the evidence used to validate or derive every T-cell subtype.
write_csv(t_scores, file.path(output_dir, "tcell_subcluster_marker_scores.csv"))
# Write the corrected author-map comparison as a separate audit artifact.
write_csv(t_validation, file.path(output_dir, "tcell_annotation_validation.csv"))

# Recluster the broad monocyte compartment to reproduce its three expression states.
message("Reclustering the broad monocyte compartment.")
if (any(s$cluster_ident == "Monocytes")) {
  monocytes <- recluster_compartment(s, "Monocytes", dims, resolution, seed)
# Preserve actual Seurat IDs rather than their one-based R factor codes.
monocyte_cluster_ids <- as.character(monocytes$seurat_clusters)
# Add an aligned monocyte-subcluster column to the all-cell object.
s$seurat_clusters_monocytes <- NA_character_
# Transfer character IDs by cell names to prevent order or coercion errors.
s@meta.data[colnames(monocytes), "seurat_clusters_monocytes"] <- monocyte_cluster_ids

# Define monocyte marker panels matching the names used in published Figure 1D.
monocyte_signatures <- list(
  # The classical CD14-positive, FCGR3A-negative, lower-CD68 state.
  "CD14+/CD16-/CD68lo Monocytes" = c("CD14", "S100A8", "S100A9", "VCAN"),
  # The intermediate/non-classical FCGR3A/MS4A7 state.
  "CD14+/CD16+/CD68mid Monocytes" = c("FCGR3A", "MS4A7", "LST1", "IFITM3"),
  # The higher-CD68 phagolysosomal state annotated by the authors.
  "CD14+/CD16+/CD68hi Monocytes" = c("CD68", "CTSS", "FCER1G", "TYROBP")
)
# Score each monocyte subcluster from its RNA expression evidence.
monocyte_scores <- score_cluster_signatures(
  monocytes, "seurat_clusters", monocyte_signatures
)
# Produce a deterministic marker-only fallback for changed cluster solutions.
monocyte_predictions <- best_signature_by_cluster(monocyte_scores)
# Correct the authors' factor-code assignments to the actual 0-based cluster IDs.
author_monocyte_map <- c(
  "0" = "CD14+/CD16+/CD68hi Monocytes",
  "1" = "CD14+/CD16-/CD68lo Monocytes",
  "2" = "CD14+/CD16+/CD68hi Monocytes",
  "3" = "CD14+/CD16+/CD68mid Monocytes",
  "4" = "CD14+/CD16+/CD68mid Monocytes",
  "5" = "CD14+/CD16+/CD68mid Monocytes"
)
# Use the corrected exact author map when all six expected monocyte IDs exist.
if (setequal(unique(monocyte_cluster_ids), names(author_monocyte_map))) {
  # Apply the authors' intended biological labels without factor conversion.
  monocyte_map <- author_monocyte_map
  # Compare the intended labels to standardized marker evidence for auditing.
  monocyte_validation <- data.frame(
    cluster = names(monocyte_map),
    author_label = unname(monocyte_map),
    strongest_marker_signature = unname(
      monocyte_predictions[names(monocyte_map)]
    ),
    marker_agrees = unname(monocyte_map) == unname(
      monocyte_predictions[names(monocyte_map)]
    ),
    stringsAsFactors = FALSE
  )
} else {
  # Fall back only when the software environment changes the subcluster solution.
  warning("Monocyte subclusters differ from the authors' 0:5 solution; using marker scores.")
  # Use the strongest valid canonical marker signature for each new cluster.
  monocyte_map <- monocyte_predictions
  # Record that no direct published numeric-ID comparison was applicable.
  monocyte_validation <- data.frame(
    cluster = names(monocyte_map), author_label = NA_character_,
    strongest_marker_signature = unname(monocyte_map), marker_agrees = NA,
    stringsAsFactors = FALSE
  )
}
# Transfer final monocyte-state labels to the all-cell metadata by cell barcode.
s@meta.data[colnames(monocytes), "cluster_ident"] <- unname(
  monocyte_map[monocyte_cluster_ids]
)
# Save monocyte signature strengths for biological review.
write_csv(
  monocyte_scores,
  file.path(output_dir, "monocyte_subcluster_marker_scores.csv")
)
# Save the corrected author-map validation table.
write_csv(
  monocyte_validation,
  file.path(output_dir, "monocyte_annotation_validation.csv")
)
} else {
  # Tiny smoke subsets can omit a broad lineage after marker-score fallback.
  # This is not a valid reason to discard the completed global embedding.
  warning("No cells were assigned to Monocytes; skipping monocyte reclustering for this run.")
  s$seurat_clusters_monocytes <- NA_character_
  write_csv(
    data.frame(
      cluster = character(), signature = character(), z_score = numeric(),
      mean_expression = numeric(), present_markers = character()
    ),
    file.path(output_dir, "monocyte_subcluster_marker_scores.csv")
  )
  write_csv(
    data.frame(
      cluster = character(), author_label = character(),
      strongest_marker_signature = character(), marker_agrees = logical()
    ),
    file.path(output_dir, "monocyte_annotation_validation.csv")
  )
}

# Correct the typographical "Monoctyes" label in the authors' public code/figure key.
celltype_levels <- c(
  "CD4+ T Cells",
  "CD8+ T Cells",
  "CD4+/CD8+ T Cells",
  "T Regulatory Cells",
  "NK Cells",
  "DC",
  "CD14+/CD16-/CD68lo Monocytes",
  "CD14+/CD16+/CD68mid Monocytes",
  "CD14+/CD16+/CD68hi Monocytes",
  "B Cells",
  "Plasma Cells",
  "Undetermined"
)
# Replace any unexpected missing subtype assignment with Undetermined.
s$cluster_ident[is.na(s$cluster_ident)] <- "Undetermined"
# Use a factor to preserve the authors' intended legend order.
s$cluster_ident <- factor(s$cluster_ident, levels = celltype_levels)
# Stop if an unrecognized label was accidentally converted to NA by factorization.
if (anyNA(s$cluster_ident)) stop("One or more cells have an invalid final cell-type label")

# Define the original Figure 1D colors before they are first used.
# This fixes the authors' clean-session error where color_vector was defined later.
color_vector <- c(
  "CD4+ T Cells" = "darkturquoise",
  "CD8+ T Cells" = "lawngreen",
  "CD4+/CD8+ T Cells" = "dodgerblue",
  "T Regulatory Cells" = "red",
  "NK Cells" = "lightpink",
  "DC" = "hotpink",
  "CD14+/CD16-/CD68lo Monocytes" = "gold1",
  "CD14+/CD16+/CD68mid Monocytes" = "navy",
  "CD14+/CD16+/CD68hi Monocytes" = "violet",
  "B Cells" = "darkorange",
  "Plasma Cells" = "darkviolet",
  "Undetermined" = "gray"
)

# Restore the global all-cell identity and UMAP for the requested Figure 1D panel.
Idents(s) <- s$cluster_ident
# Set the random draw order so overplotting is reproducible across repeated runs.
set.seed(seed)
# Generate Figure 1D directly from the integrated UMAP and data-derived annotations.
figure_1d <- DimPlot(
  s,
  group.by = "cluster_ident",
  cols = color_vector,
  shuffle = TRUE,
  reduction = "umap",
  pt.size = 1,
  label = FALSE
)

# Name the vector and raster versions of the requested panel.
figure_pdf <- file.path(output_dir, "figure1d_recreated.pdf")
# Keep a high-resolution PNG for convenient visual comparison outside R.
figure_png <- file.path(output_dir, "figure1d_recreated.png")
# Save at the 10-by-8-inch panel size used in the authors' source code.
ggsave(figure_pdf, figure_1d, width = 10, height = 8, useDingbats = FALSE)
# Save a 300-dpi raster counterpart without manually editing any plotted element.
ggsave(figure_png, figure_1d, width = 10, height = 8, dpi = 300)

# Create a marker dot plot that makes every final label biologically reviewable.
paper_markers <- c(
  "CD3D", "CD4", "CD8A", "CD8B", "FOXP3", "NKG7", "GNLY",
  "FCER1A", "CD68", "CD14", "FCGR3A", "CD19", "JCHAIN", "CD44"
)
# Remove any paper marker absent from the deposited expression matrix.
paper_markers <- intersect(paper_markers, rownames(s))
# Use the log-normalized RNA assay rather than integrated correction values.
DefaultAssay(s) <- "RNA"
# Plot expression prevalence and average abundance by final cell type.
marker_dotplot <- DotPlot(
  s,
  features = paper_markers,
  group.by = "cluster_ident",
  assay = "RNA"
) +
  # Rotate feature labels so all canonical genes remain legible.
  RotatedAxis() +
  # Label axes explicitly because this plot is an audit artifact, not Figure 1D.
  labs(x = "Canonical marker", y = NULL, title = "Figure 1D annotation checks")
# Save marker evidence separately from the recreated panel.
ggsave(
  file.path(diagnostics_dir, "figure1d_annotation_marker_dotplot.pdf"),
  marker_dotplot,
  width = 10,
  height = 6,
  useDingbats = FALSE
)

# Summarize final cell counts and proportions for every legend entry.
celltype_composition <- as.data.frame(table(s$cluster_ident), stringsAsFactors = FALSE)
# Give the output columns explicit, analysis-friendly names.
names(celltype_composition) <- c("cell_type", "n_cells")
# Correct the authors' self-reference bug by dividing by the actual total count.
celltype_composition$proportion <- celltype_composition$n_cells /
  sum(celltype_composition$n_cells)
# Save the composition table used to sanity-check the recreated panel.
write_csv(celltype_composition, file.path(output_dir, "celltype_composition.csv"))

# Export cell-level UMAP coordinates and annotations for independent verification.
umap_export <- as.data.frame(Embeddings(s, reduction = "umap"))
# Add cell barcodes as a normal column rather than relying on CSV row names.
umap_export$cell <- rownames(umap_export)
# Copy only metadata needed to trace every plotted point.
umap_export$sample_id <- as.character(s$sample_id[rownames(umap_export)])
# Record the all-cell cluster that determined the broad compartment.
umap_export$seurat_cluster_allcells <- as.character(
  s$seurat_clusters_allcells[rownames(umap_export)]
)
# Record T subclusters only for cells that entered T-cell reclustering.
umap_export$seurat_cluster_tcells <- s$seurat_clusters_tcells[rownames(umap_export)]
# Record monocyte subclusters only for cells that entered myeloid reclustering.
umap_export$seurat_cluster_monocytes <- s$seurat_clusters_monocytes[rownames(umap_export)]
# Record the final plotted cell-type identity.
umap_export$cell_type <- as.character(s$cluster_ident[rownames(umap_export)])
# Write a compressed text table that remains easy to inspect and share.
write_csv(umap_export, file.path(output_dir, "figure1d_umap_coordinates.csv.gz"))

# Save compartment diagnostic UMAPs without substituting them for the all-cell panel.
set.seed(seed)
# Draw T-cell subcluster IDs to inspect the first hierarchical annotation stage.
t_diagnostic <- DimPlot(
  t_cells, group.by = "seurat_clusters", reduction = "umap",
  label = TRUE, shuffle = TRUE, pt.size = 0.4
)
# Save the T-cell diagnostic plot at a compact review size.
ggsave(
  file.path(diagnostics_dir, "tcell_reclustering_umap.pdf"),
  t_diagnostic, width = 6, height = 5, useDingbats = FALSE
)
# Draw monocyte subcluster IDs to inspect the second hierarchical stage when
# the current input contains a monocyte compartment.
if (exists("monocytes")) {
  monocyte_diagnostic <- DimPlot(
    monocytes, group.by = "seurat_clusters", reduction = "umap",
    label = TRUE, shuffle = TRUE, pt.size = 0.4
  )
  # Save the monocyte diagnostic plot at the same review size.
  ggsave(
    file.path(diagnostics_dir, "monocyte_reclustering_umap.pdf"),
    monocyte_diagnostic, width = 6, height = 5, useDingbats = FALSE
  )
}

# Remove compartment copies before serializing the already complete all-cell object.
rm(t_cells)
if (exists("monocytes")) rm(monocytes)
# Reclaim their assays and reductions to reduce peak final-save memory.
invisible(gc())
# Save the full annotated object so no downstream user must rerun integration.
saveRDS(s, final_object_path, compress = FALSE)

# Capture exact R and package versions for reproducibility across environments.
session_info_path <- file.path(output_dir, "sessionInfo.txt")
# Divert printed session information into the audit text file.
capture.output(sessionInfo(), file = session_info_path)

# Report the main deliverables after every write has completed successfully.
message("Figure 1D recreation completed.")
# Print the vector panel path suitable for publication workflows.
message("  PDF: ", figure_pdf)
# Print the high-resolution raster panel path suitable for quick inspection.
message("  PNG: ", figure_png)
# Print the final Seurat object path for downstream analysis.
message("  Seurat object: ", final_object_path)
# Point to the marker evidence that validates the biological annotations.
message("  Diagnostics: ", diagnostics_dir)
