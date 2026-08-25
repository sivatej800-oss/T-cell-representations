# Stage 14: the authors' own marker visualisation, applied to both partitions.
#
# The authors interpreted their clusters by setting the default assay to expression,
# plotting each canonical marker over the UMAP with FeaturePlot, and reading the
# result by eye before typing cluster numbers into a manual annotation block. This
# stage reproduces that procedure. It is deliberately not a rule engine: the point is
# to put the same evidence in front of a person that the authors put in front of
# themselves, for both partitions, drawn the same way.
#
# The only intended divergences from their script are the two this project exists to
# study -- the scGPT representation, and Harmony in place of anchor integration. The
# marker panel, the plot type, the colours, and the point size are theirs.
#
# One further divergence is forced rather than chosen: the authors plotted the SCT
# assay, whose single merged model covered one gene set. Here SCTransform ran per
# sample, and the 59 checkpoints share only 9,013 of 19,689 genes, so the RNA assay is
# used instead. It carries every gene for every cell and is identical for both arms.
#
# Nothing is clustered here and no identity is assigned. Stage 15 does the annotation,
# by hand, the way the authors did.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop("Usage: 14_marker_evidence.R <seurat_input_dir> <frozen_dir> <benchmark_dir> <marker_dir> <output_dir> <figure_dir>")
}
input_dir <- args[1]
frozen_dir <- args[2]
benchmark_dir <- args[3]
marker_dir <- args[4]
output_dir <- args[5]
figure_dir <- args[6]

expected_cells <- 70391
expected_genes <- 19689

partitions <- c("pca30_harmony_cluster", "scgpt_harmony_cluster")
partition_labels <- c(
  pca30_harmony_cluster = "PCA-30 Harmony",
  scgpt_harmony_cluster = "scGPT Harmony"
)
# Which Stage 06 UMAP belongs to which partition. Each arm is drawn in its own layout,
# because a partition plotted over the other arm's coordinates would look incoherent
# for reasons that have nothing to do with the biology.
umap_prefix <- c(
  pca30_harmony_cluster = "pca_harmony",
  scgpt_harmony_cluster = "scgpt_harmony"
)

# The authors' panel, from their line 121, in their order.
marker_panel <- c(
  "CD3D", "CD4", "CD8A", "CD8B", "FOXP3",
  "NKG7", "GNLY", "FCER1A", "CD68",
  "CD14", "FCGR3A", "CD19", "JCHAIN", "CD44"
)

# The authors' FeaturePlot settings, from their line 133.
FEATURE_COLS <- c("grey88", "deeppink")
FEATURE_PT_SIZE <- 1

TOP_N_TABLE <- 20

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
marker_figure_dir <- file.path(figure_dir, "14_feature_plots")
dir.create(marker_figure_dir, showWarnings = FALSE, recursive = TRUE)

available_gb <- tryCatch({
  info <- readLines("/proc/meminfo")
  as.numeric(gsub("[^0-9]", "", grep("^MemAvailable:", info, value = TRUE))) / 1024 / 1024
}, error = function(e) NA_real_)
n_workers <- if (is.na(available_gb)) 2 else {
  max(1, min(detectCores() - 2, 8, floor((available_gb * 0.5) / 3)))
}
message("Memory available: ", round(available_gb, 1), " GB; using ", n_workers, " workers")

# --- Step 1: the shared expression object, as Stage 13 built it ------------------
message("Reading the RNA matrix (peak memory is here)...")
counts <- Matrix::readMM(file.path(input_dir, "counts.mtx"))
counts <- as(counts, "CsparseMatrix")
genes <- readLines(file.path(input_dir, "genes.tsv"), warn = FALSE)
barcodes <- readLines(file.path(input_dir, "barcodes.tsv"), warn = FALSE)
if (nrow(counts) != expected_genes || ncol(counts) != expected_cells) {
  stop("Expected ", expected_genes, " x ", expected_cells, ", found ",
       nrow(counts), " x ", ncol(counts))
}
rownames(counts) <- genes
colnames(counts) <- barcodes
invisible(gc(full = TRUE))

frozen <- utils::read.csv(file.path(frozen_dir, "frozen_cluster_labels.csv.gz"),
                          stringsAsFactors = FALSE, colClasses = "character")
rownames(frozen) <- frozen$cell_barcode
if (!setequal(frozen$cell_barcode, barcodes)) {
  stop("The frozen labels and the RNA matrix do not describe the same cells")
}
frozen <- frozen[barcodes, , drop = FALSE]

metadata <- utils::read.csv(file.path(input_dir, "metadata.csv"), stringsAsFactors = FALSE)
rownames(metadata) <- metadata$cell_barcode
metadata <- metadata[barcodes, , drop = FALSE]
frozen$internal_sample_id <- metadata$internal_sample_id

message("Creating the shared Seurat object...")
csf <- CreateSeuratObject(
  counts = counts, assay = "RNA", project = "GSE200164_evidence",
  meta.data = frozen[, c(partitions, "internal_sample_id")],
  min.cells = 0, min.features = 0
)
rm(counts)
invisible(gc(full = TRUE))

csf <- NormalizeData(csf, normalization.method = "LogNormalize",
                     scale.factor = 10000, verbose = FALSE)
DefaultAssay(csf) <- "RNA"

audit <- jsonlite::fromJSON(file.path(marker_dir, "marker_discovery_audit.json"))
observed_sum <- sum(GetAssayData(csf, assay = "RNA", layer = "data"))
if (!isTRUE(all.equal(observed_sum, audit$expression_fingerprint$data_sum, tolerance = 1e-8))) {
  stop("This expression does not match the one Stage 13 tested against")
}
message("Expression matches Stage 13's fingerprint: ", format(observed_sum, digits = 16))

missing_panel <- setdiff(marker_panel, rownames(csf))
if (length(missing_panel) > 0) {
  stop("Panel genes absent from the matrix: ", paste(missing_panel, collapse = ", "))
}

# --- Step 2: attach each arm's UMAP, so FeaturePlot has a layout to draw on -------
# The coordinates come from Stage 06, which computed them on the corrected
# representation each arm is defined by. They are read, not recomputed.
message("Attaching Stage 06 UMAP layouts...")
layouts <- utils::read.csv(
  file.path(benchmark_dir, "benchmark_cluster_assignments.csv.gz"),
  stringsAsFactors = FALSE
)
rownames(layouts) <- layouts$cell_barcode
if (!setequal(layouts$cell_barcode, barcodes)) {
  stop("Stage 06's layouts do not describe the same cells as the expression")
}
layouts <- layouts[barcodes, , drop = FALSE]

for (partition in partitions) {
  prefix <- umap_prefix[[partition]]
  coordinates <- as.matrix(layouts[, paste0(prefix, c("_UMAP_1", "_UMAP_2"))])
  colnames(coordinates) <- paste0(prefix, "UMAP_", 1:2)
  rownames(coordinates) <- colnames(csf)
  if (any(!is.finite(coordinates))) {
    stop("Stage 06's UMAP for ", partition, " contains non-finite coordinates")
  }
  csf[[paste0("umap_", prefix)]] <- CreateDimReducObject(
    embeddings = coordinates, key = paste0(prefix, "UMAP_"), assay = "RNA"
  )
}

# --- Step 3: FeaturePlot per marker, the authors' way ----------------------------
# One plot per marker per partition, same colours and point size as their script.
# Written individually as they wrote them, and also assembled into one panel per
# partition, because fourteen separate files are hard to read side by side and the
# comparison between arms is the entire purpose here.
message("Drawing feature plots...")

# The authors' call and nothing more (their lines 132-135). In particular there is
# no `order` argument: Seurat defaults it to FALSE, and setting it TRUE redraws the
# expressing cells last, on top of the grey ones, which makes a marker read as
# stronger and more contiguous than it is. The method here is reading these panels
# by eye, so that one argument would change the answer.
feature_plot <- function(partition, marker) {
  FeaturePlot(
    csf, features = marker, pt.size = FEATURE_PT_SIZE,
    cols = FEATURE_COLS, reduction = paste0("umap_", umap_prefix[[partition]])
  )
}

# Smaller type is only needed by the fourteen-panel grid, which is this project's
# addition rather than theirs. The per-marker files stay as the authors drew them.
compact_for_grid <- function(plot) {
  plot + theme(plot.title = element_text(size = 11, face = "bold"),
               axis.title = element_text(size = 7), axis.text = element_text(size = 6),
               legend.key.width = unit(0.25, "cm"), legend.text = element_text(size = 6))
}

# Drawn in a plain loop, the way the authors drew them. An earlier version forked
# this with mclapply and deadlocked: a process that has already driven a threaded
# BLAS and a graphics device cannot reliably fork children that finish. It showed
# up as four zero-byte PNGs -- exactly the first worker's block -- and a parent that
# waited five minutes for results that never arrived. Drawing is a few minutes
# either way, so the loop is not worth the failure mode.
for (partition in partitions) {
  message("  ", partition)
  plots <- vector("list", length(marker_panel))
  names(plots) <- marker_panel
  for (marker in marker_panel) {
    plot <- feature_plot(partition, marker)
    target <- file.path(marker_figure_dir, paste0(partition, "_", marker, ".png"))
    ggsave(target, plot, width = 5, height = 4, dpi = 150)
    if (!file.exists(target) || file.info(target)$size == 0) {
      stop("FeaturePlot wrote nothing for ", partition, " / ", marker)
    }
    plots[[marker]] <- plot
  }

  combined <- wrap_plots(lapply(plots, compact_for_grid), ncol = 4) +
    plot_annotation(
      title = paste0(partition_labels[[partition]], " — canonical markers on RNA expression"),
      theme = theme(plot.title = element_text(size = 14, face = "bold"))
    )
  ggsave(file.path(figure_dir, paste0("14_featureplots_", partition, ".png")),
         combined, width = 18, height = 16, dpi = 130)
}

# --- Step 4: the cluster layout itself, as the authors also plotted -------------
message("Drawing cluster layouts...")
for (partition in partitions) {
  Idents(csf) <- factor(csf@meta.data[[partition]],
                        levels = sort(unique(as.integer(csf@meta.data[[partition]]))))
  plot <- DimPlot(csf, reduction = paste0("umap_", umap_prefix[[partition]]),
                  label = TRUE, pt.size = 1, shuffle = TRUE) +
    ggtitle(partition_labels[[partition]]) + NoLegend()
  target <- file.path(figure_dir, paste0("14_clusters_", partition, ".png"))
  ggsave(target, plot, width = 6, height = 5.5, dpi = 150)
  if (!file.exists(target) || file.info(target)$size == 0) {
    stop("DimPlot wrote nothing for ", partition)
  }
}

# --- Step 5: the numbers behind the plots ----------------------------------------
# The authors read their annotation off the plots alone. The same numbers are written
# here as a table, because a percentage is easier to check than a colour, and because
# the manual annotation in Stage 15 should be auditable after the fact.
message("Scoring the canonical panel...")
panel_data <- GetAssayData(csf, assay = "RNA", layer = "data")[marker_panel, , drop = FALSE]
label_lookup <- lapply(partitions, function(p) csf@meta.data[[p]])
names(label_lookup) <- partitions

panel_jobs <- do.call(rbind, lapply(partitions, function(partition) {
  data.frame(partition = partition,
             cluster = sort(unique(label_lookup[[partition]]), method = "radix"),
             stringsAsFactors = FALSE)
}))

panel_rows <- mclapply(seq_len(nrow(panel_jobs)), function(i) {
  partition <- panel_jobs$partition[i]; cluster <- panel_jobs$cluster[i]
  in_cluster <- label_lookup[[partition]] == cluster
  block <- panel_data[, in_cluster, drop = FALSE]
  data.frame(
    partition = partition, cluster = cluster, n_cells = sum(in_cluster),
    gene = factor(marker_panel, levels = marker_panel),
    mean_expression = as.numeric(Matrix::rowMeans(block)),
    pct_expressing = 100 * as.numeric(Matrix::rowSums(block > 0)) / sum(in_cluster),
    stringsAsFactors = FALSE
  )
}, mc.cores = n_workers)
if (any(vapply(panel_rows, inherits, logical(1), "try-error"))) stop("Panel scoring failed")
panel_table <- do.call(rbind, panel_rows)
utils::write.csv(panel_table, file.path(output_dir, "canonical_panel_by_cluster.csv"),
                 row.names = FALSE)

for (partition in partitions) {
  block <- panel_table[panel_table$partition == partition, ]
  wide <- reshape(block[, c("cluster", "gene", "pct_expressing")],
                  idvar = "cluster", timevar = "gene", direction = "wide")
  colnames(wide) <- sub("^pct_expressing\\.", "", colnames(wide))
  wide <- wide[order(as.integer(wide$cluster)), ]
  utils::write.csv(wide, file.path(output_dir,
                   paste0("canonical_pct_expressing_", partition, ".csv")), row.names = FALSE)
}

# --- Step 6: ranked markers, and one evidence row per cluster --------------------
message("Ranking positive markers...")
top_tables <- list()
for (partition in partitions) {
  markers <- utils::read.csv(file.path(marker_dir, paste0("markers_", partition, ".csv.gz")),
                             stringsAsFactors = FALSE)
  markers$cluster <- as.character(markers$cluster)
  positive <- markers[markers$p_val_adj < 0.05 & markers$avg_log2FC > 0, ]
  ranked <- do.call(rbind, lapply(sort(unique(positive$cluster), method = "radix"), function(cluster) {
    block <- positive[positive$cluster == cluster, ]
    block <- utils::head(block[order(-block$avg_log2FC), ], TOP_N_TABLE)
    block$rank <- seq_len(nrow(block))
    block
  }))
  ranked$partition <- partition
  utils::write.csv(
    ranked[, c("partition", "cluster", "rank", "gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")],
    file.path(output_dir, paste0("top_markers_", partition, ".csv")), row.names = FALSE
  )
  top_tables[[partition]] <- ranked
}

message("Building the evidence table...")
PRESENT <- 25
evidence_rows <- list()
for (partition in partitions) {
  labels <- label_lookup[[partition]]
  donors <- csf@meta.data$internal_sample_id
  ranked <- top_tables[[partition]]
  block_all <- panel_table[panel_table$partition == partition, ]
  rows <- mclapply(sort(unique(labels), method = "radix"), function(cluster) {
    in_cluster <- labels == cluster
    panel_block <- block_all[block_all$cluster == cluster, ]
    present <- as.character(panel_block$gene[panel_block$pct_expressing >= PRESENT])
    donor_counts <- table(donors[in_cluster])
    donor_counts <- sort(donor_counts[donor_counts > 0], decreasing = TRUE)
    data.frame(
      partition = partition, cluster = cluster, n_cells = sum(in_cluster),
      n_donors = length(donor_counts),
      largest_donor_share = round(100 * max(donor_counts) / sum(in_cluster), 1),
      top_markers = paste(utils::head(ranked$gene[ranked$cluster == cluster], 20), collapse = ", "),
      canonical_present = if (length(present)) paste(present, collapse = ", ") else "none",
      stringsAsFactors = FALSE
    )
  }, mc.cores = n_workers)
  if (any(vapply(rows, inherits, logical(1), "try-error"))) stop("Evidence table failed")
  evidence_rows <- c(evidence_rows, rows)
}
utils::write.csv(do.call(rbind, evidence_rows),
                 file.path(output_dir, "cluster_evidence_table.csv"), row.names = FALSE)

writeLines(jsonlite::toJSON(list(
  stage = "14_marker_evidence",
  follows = "author_reference/figure1d_author_workflow.R lines 119-140",
  partitions = partitions, panel = marker_panel,
  plot_type = "FeaturePlot", feature_cols = FEATURE_COLS, feature_pt_size = FEATURE_PT_SIZE,
  assay = "RNA", assay_divergence = "authors used SCT; per-sample SCT here shares only 9013 of 19689 genes",
  presence_threshold_pct = PRESENT, expression_fingerprint = observed_sum,
  seurat_version = as.character(packageVersion("Seurat")),
  r_version = paste(R.version$major, R.version$minor, sep = ".")
), auto_unbox = TRUE, pretty = TRUE, digits = 10),
file.path(output_dir, "marker_evidence_audit.json"))

message("\nStage 14 complete.")
