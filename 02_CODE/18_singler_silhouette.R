# Stage 18 (supplementary): how well each space separates the reference cell types.
#
# Stage 17 asks whether each arm's clusters follow SingleR's cell types. This asks the
# question the other way round, without any clustering in between: take SingleR's
# broad labels as they are, and measure how cleanly each representation keeps those
# cell types apart. The spaces are the two in the main comparison -- the authors'
# integrated PCA (PCs 1-15, Stage 03) and the uncorrected scGPT embedding (all 512
# dimensions, the coordinates Stage 06's scgpt_raw arm clustered on).
#
# The two SingleR numbers cannot be compared with each other directly. Stage 07 sets
# out why: silhouette is built from distances, and 15 variance-ordered components and
# 512 learned dimensions produce distances that are not the same quantity; distances
# in high dimensions also concentrate toward a common value, which lowers silhouette
# for reasons unrelated to the grouping. So each space also scores both arms' own
# cluster partitions, on the same cells. Within one space the geometry is fixed and
# every row is comparable: the useful readings are how SingleR's cell types score
# against that space's own clusters, and which arm's clusters score higher there.
#
# Every labelling is scored on the same cells: the 69,861 that SingleR assigned. Its
# 530 unassigned cells are not a cell type and are left out of every row, including
# the cluster rows, so that no row is scored on a different population. SingleR's
# labels are scored twice, as assigned and with DC and Monocytes merged into Myeloid,
# matching Stage 17: the two are adjacent in expression and SingleR's split between
# them is the least certain call it makes.
#
# The method is Stage 07's unchanged: Euclidean distances on the coordinates exactly
# as the clustering used them, 12,000 cells per draw sampled in proportion to each
# group with at least two from every group, five draws, seed 42. A label with only a
# handful of cells (SingleR calls 17 Progenitors and 5 Plasma Cells) contributes
# almost nothing to the mean, and its per-label value is not interpretable.
#
# Usage: Rscript 02_CODE/18_singler_silhouette.R <workspace>

suppressPackageStartupMessages({
  library(Seurat)
  library(cluster)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: 18_singler_silhouette.R <workspace>")
results_dir <- file.path(normalizePath(args[1]), "results")
output_dir <- file.path(results_dir, "18_singler_silhouette")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

expected_cells <- 70391
scgpt_dims <- 512
random_seed <- 42
sample_size <- 12000
n_repeats <- 5


# --- Step 1: the labels ------------------------------------------------------------
frozen <- utils::read.csv(file.path(results_dir, "12_frozen_labels", "frozen_cluster_labels.csv.gz"),
                          colClasses = "character")
singler <- utils::read.csv(file.path(results_dir, "16_singler_labels", "singler_cell_labels.csv.gz"),
                           colClasses = "character")
if (nrow(frozen) != expected_cells || nrow(singler) != expected_cells) {
  stop("Expected ", expected_cells, " cells in both the frozen labels and SingleR's")
}
if (!setequal(frozen$cell_barcode, singler$cell_barcode)) {
  stop("The frozen labels and SingleR describe different cells")
}
barcodes <- frozen$cell_barcode
singler <- singler[match(barcodes, singler$cell_barcode), ]
if (!identical(singler$cell_barcode, barcodes)) stop("SingleR labels could not be aligned")

labelings <- list(
  singler_broad = singler$broad_label,
  singler_myeloid_merged = ifelse(singler$broad_label %in% c("DC", "Monocytes"),
                                  "Myeloid", singler$broad_label),
  authors_pipeline_cluster = frozen$authors_pipeline_cluster,
  scgpt_raw_cluster = frozen$scgpt_raw_cluster
)
assigned <- singler$broad_label != "unassigned"
message("Cells SingleR assigned: ", sum(assigned), " of ", length(assigned))


# --- Step 2: the authors' PCA space ------------------------------------------------
# Read from Stage 03's final object, once, then kept as a small file: the object is
# 6 GB and only 15 columns of it are needed. The object's clusters are checked against
# the frozen partition, which proves these are the coordinates that partition came from.
pca_cache <- file.path(output_dir, "authors_pca15_coordinates.rds")
if (file.exists(pca_cache)) {
  authors_pca <- readRDS(pca_cache)
} else {
  message("Reading Stage 03's integrated object (once)...")
  integrated <- readRDS(file.path(results_dir, "03_seurat_reproduction", "integrated_clustered_umap.rds"))
  authors_pca <- Embeddings(integrated, reduction = "pca")[, 1:15, drop = FALSE]
  object_clusters <- as.character(integrated$seurat_clusters)
  names(object_clusters) <- colnames(integrated)
  rm(integrated); invisible(gc())
  if (!setequal(rownames(authors_pca), barcodes)) {
    stop("Stage 03's object does not describe the same cells as the frozen labels")
  }
  if (!identical(unname(object_clusters[barcodes]), frozen$authors_pipeline_cluster)) {
    stop("Stage 03's object clusters differ from the frozen authors' partition")
  }
  saveRDS(authors_pca, pca_cache)
}
authors_pca <- authors_pca[barcodes, , drop = FALSE]


# --- Step 3: the uncorrected scGPT space -------------------------------------------
# Read exactly as Stage 04 reads it, including its check that every row is unit
# length, which is what catches a transposed read of the flat binary block.
input_dir <- file.path(results_dir, "06_harmony_benchmark", "r_input")
connection <- file(file.path(input_dir, "embedding_float32.bin"), "rb")
values <- readBin(connection, what = "numeric", n = expected_cells * scgpt_dims,
                  size = 4, endian = "little")
close(connection)
if (length(values) != expected_cells * scgpt_dims) stop("Embedding binary is the wrong size")
scgpt <- matrix(values, nrow = expected_cells, ncol = scgpt_dims, byrow = TRUE)
rm(values)
if (any(!is.finite(scgpt)) || max(abs(sqrt(rowSums(scgpt^2)) - 1)) > 1e-4) {
  stop("Embedding rows are not unit length; the binary block was misread")
}
embedding_cells <- utils::read.csv(file.path(input_dir, "metadata.csv"),
                                   colClasses = "character")$cell_barcode
if (!setequal(embedding_cells, barcodes)) stop("The embedding describes different cells")
rownames(scgpt) <- embedding_cells
scgpt <- scgpt[barcodes, , drop = FALSE]

spaces <- list(authors_pca15 = authors_pca, scgpt_raw_512 = scgpt)


# --- Step 4: score every labelling in every space ----------------------------------
# Stage 07's procedure. The draws are seeded per repeat and not per space, so within
# one labelling every space sees the same cells.
stratified_sample <- function(labels, size) {
  counts <- table(labels)
  per_group <- pmax(2, round(size * as.vector(counts) / length(labels)))
  names(per_group) <- names(counts)
  unlist(lapply(levels(labels), function(level) {
    available <- which(labels == level)
    sample(available, min(length(available), per_group[[level]]))
  }))
}

# The draws are fixed here, in the parent, before any scoring starts. Each silhouette
# (labelling x space x draw) is then independent, and they are scored in parallel:
# one 12,000-cell silhouette in 512 dimensions takes minutes, and forty in sequence
# kept one core busy for hours. Because the draws are fixed first, the parallel run
# gives exactly the numbers a sequential one would.
labels_by <- lapply(labelings, function(x) factor(x[assigned]))
draws_by <- lapply(labels_by, function(labels) {
  lapply(seq_len(n_repeats), function(i) {
    set.seed(random_seed + i)
    stratified_sample(labels, sample_size)
  })
})
spaces <- lapply(spaces, function(coordinates) coordinates[assigned, , drop = FALSE])
jobs <- expand.grid(draw = seq_len(n_repeats), space = names(spaces),
                    labelling = names(labelings), stringsAsFactors = FALSE)
# Every core is used. Each worker peaks near 1.2 GB (one 12,000-cell distance matrix
# plus working copies), so sixteen fit comfortably in this machine's memory.
n_workers <- max(1, min(parallel::detectCores(), nrow(jobs)))
message("Scoring ", nrow(jobs), " silhouettes on ", n_workers, " workers...")
scored <- parallel::mclapply(seq_len(nrow(jobs)), function(j) {
  labels <- labels_by[[jobs$labelling[j]]]
  chosen <- draws_by[[jobs$labelling[j]]][[jobs$draw[j]]]
  drawn <- droplevels(labels[chosen])
  widths <- silhouette(as.integer(drawn),
                       dist(spaces[[jobs$space[j]]][chosen, , drop = FALSE]))[, "sil_width"]
  list(overall = mean(widths), negative = mean(widths < 0),
       per_label = tapply(widths, drawn, mean))
}, mc.cores = n_workers, mc.preschedule = FALSE)
if (any(vapply(scored, inherits, logical(1), "try-error"))) stop("A silhouette job failed")

summary_rows <- list()
per_label_rows <- list()
for (labelling in names(labelings)) {
  labels <- labels_by[[labelling]]
  for (space in names(spaces)) {
    per_draw <- scored[jobs$labelling == labelling & jobs$space == space]
    overall <- vapply(per_draw, `[[`, numeric(1), "overall")
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      labelling = labelling, space = space, n_groups = nlevels(labels),
      n_cells = length(labels), mean_silhouette = mean(overall),
      sd_across_draws = sd(overall),
      negative_fraction = mean(vapply(per_draw, `[[`, numeric(1), "negative")))
    per_label_rows[[length(per_label_rows) + 1]] <- data.frame(
      labelling = labelling, space = space, label = levels(labels),
      n_cells = as.vector(table(labels)),
      mean_silhouette = vapply(levels(labels), function(level) {
        mean(vapply(per_draw, function(d) d$per_label[[level]], numeric(1)))
      }, numeric(1)))
  }
}

summary_table <- do.call(rbind, summary_rows)
per_label_table <- do.call(rbind, per_label_rows)
utils::write.csv(summary_table, file.path(output_dir, "singler_silhouette_summary.csv"), row.names = FALSE)
utils::write.csv(per_label_table, file.path(output_dir, "singler_silhouette_per_label.csv"), row.names = FALSE)

writeLines(toJSON(list(
  stage = "18_singler_silhouette",
  spaces = list(authors_pca15 = "Stage 03 integrated PCA, PCs 1-15",
                scgpt_raw_512 = "scGPT embedding, all 512 dims, uncorrected"),
  cells_scored = sum(assigned), cells_excluded_singler_unassigned = sum(!assigned),
  sample_size = sample_size, repeats = n_repeats, seed = random_seed,
  distance = "euclidean", method = "Stage 07 stratified subsample",
  comparable = "within a space only; never across spaces",
  seurat_version = as.character(packageVersion("Seurat")),
  r_version = paste(R.version$major, R.version$minor, sep = ".")
), auto_unbox = TRUE, pretty = TRUE), file.path(output_dir, "singler_silhouette_audit.json"))

print(summary_table, row.names = FALSE)
message("\nStage 18 complete: ", output_dir)
