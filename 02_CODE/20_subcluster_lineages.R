# Stage 20: recluster the T cells and the monocytes, in each arm, the authors' way.
#
# The broad annotation (Stage 15) is the authors' first pass. Their Figure 1D labels --
# CD4+, CD8+, CD4+/CD8+ and regulatory T cells, and three monocyte subtypes -- come from
# a second pass: they took the cells labelled T Cells, and separately Monocytes, and
# clustered each set again on its own (05_REFERENCES/author_code/
# figure1d_author_workflow.R, lines 211-271). This stage does that second pass for both
# arms, on each arm's own Stage 15 labels.
#
#   Authors' pipeline: their recipe exactly. Subset the integrated object; on the
#     integrated assay re-find variable features and re-run PCA; neighbours on PCs
#     1-15, Louvain at resolution 0.3, UMAP with 30 neighbours. There is no ScaleData
#     call, because they made none: the integrated assay's stored scale.data is used.
#   scGPT raw: the same cells reclustered on the 512-dimension embedding with the same
#     neighbour, Louvain and UMAP settings. This is the one departure, and it is forced:
#     an embedding has no assay to re-find variable features or re-run PCA on. It is
#     how Stage 06 clustered the full set.
#
# Each subcluster is then described from the RNA, as Stages 13-14 described the broad
# clusters: one LogNormalized RNA matrix shared by both arms (checked against Stage
# 13's fingerprint), Wilcoxon markers with the authors' thresholds, and the authors'
# reclustering panel (their lines 274-275) as the share of cells with a raw count above
# zero and as mean expression. Mean expression is kept because the authors' monocyte
# labels are levels (CD68 hi, mid, lo), not presence. RNA is used rather than SCT for
# the reason Stage 13 gives: the per-sample SCT assays share only 9,013 of 19,689 genes.
#
# No subcluster is named here. The labels are assigned by hand afterwards, as the
# authors did, and the arms are compared only after that.
#
# Usage: Rscript 02_CODE/20_subcluster_lineages.R <workspace> <authors_pipeline|scgpt_raw>

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(jsonlite)
})
if (!requireNamespace("presto", quietly = TRUE)) {
  stop("presto is not installed; Stage 13's marker test used it and this stage must match")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2 || !args[2] %in% c("authors_pipeline", "scgpt_raw")) {
  stop("Usage: 20_subcluster_lineages.R <workspace> <authors_pipeline|scgpt_raw>")
}
results_dir <- file.path(normalizePath(args[1]), "results")
arm <- args[2]
output_dir <- file.path(results_dir, "20_subclusters")
figure_dir <- file.path(results_dir, "figures")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

expected_cells <- 70391
expected_genes <- 19689
stage13_data_sum <- 176503092.94553709
scgpt_dims <- 512
random_seed <- 42
lineages <- c("T Cells", "Monocytes")
# The authors' reclustering panel, their lines 274-275, in their order.
panel <- c("CD3D", "CD4", "CD8A", "CD8B", "FOXP3", "CD14", "FCGR3A", "CD68")
TOP_N <- 20


# --- Step 1: which cells each lineage holds, from Stage 15 --------------------------
annotation <- utils::read.csv(file.path(results_dir, "15_annotation", "annotated_cluster_labels.csv.gz"),
                              colClasses = "character")
identity <- annotation[[paste0(arm, "_cluster_identity")]]
broad_cluster <- annotation[[paste0(arm, "_cluster")]]
if (nrow(annotation) != expected_cells || is.null(identity) || anyNA(identity)) {
  stop("Stage 15's file does not hold a complete annotation for ", arm)
}
names(identity) <- annotation$cell_barcode


# --- Step 2: the shared RNA expression ----------------------------------------------
message("Reading the RNA matrix...")
input_dir <- file.path(results_dir, "03_seurat_reproduction", "seurat_input")
counts <- as(Matrix::readMM(file.path(input_dir, "counts.mtx")), "CsparseMatrix")
genes <- readLines(file.path(input_dir, "genes.tsv"), warn = FALSE)
barcodes <- readLines(file.path(input_dir, "barcodes.tsv"), warn = FALSE)
if (nrow(counts) != expected_genes || ncol(counts) != expected_cells) {
  stop("Expected ", expected_genes, " x ", expected_cells, ", found ", nrow(counts), " x ", ncol(counts))
}
rownames(counts) <- genes
colnames(counts) <- barcodes
if (!setequal(barcodes, names(identity))) stop("Stage 15 and the RNA matrix describe different cells")
missing_panel <- setdiff(panel, genes)
if (length(missing_panel) > 0) stop("Panel genes absent: ", paste(missing_panel, collapse = ", "))

metadata <- utils::read.csv(file.path(input_dir, "metadata.csv"), colClasses = "character")
patients <- metadata$internal_sample_id[match(barcodes, metadata$cell_barcode)]
if (anyNA(patients)) stop("A cell has no patient")

rna <- CreateSeuratObject(counts = counts, assay = "RNA", project = "GSE200164_subclusters",
                          min.cells = 0, min.features = 0)
rna$patient <- patients
rna <- NormalizeData(rna, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
observed_sum <- sum(GetAssayData(rna, assay = "RNA", layer = "data"))
if (abs(observed_sum - stage13_data_sum) / stage13_data_sum > 1e-9) {
  stop("The RNA expression does not match Stage 13's fingerprint")
}
rm(counts); invisible(gc())


# --- Step 3: the representation each arm reclusters in ------------------------------
if (arm == "authors_pipeline") {
  message("Reading Stage 03's integrated object...")
  integrated <- readRDS(file.path(results_dir, "03_seurat_reproduction", "integrated_clustered_umap.rds"))
  if (!setequal(colnames(integrated), barcodes)) stop("Stage 03's object holds different cells")
  object_clusters <- as.character(integrated$seurat_clusters)[match(barcodes, colnames(integrated))]
  if (!identical(object_clusters, broad_cluster[match(barcodes, annotation$cell_barcode)])) {
    stop("Stage 03's object clusters differ from the frozen authors' partition")
  }
  if (!"integrated" %in% Assays(integrated)) stop("Stage 03's object has no integrated assay")
} else {
  input <- file.path(results_dir, "06_harmony_benchmark", "r_input")
  connection <- file(file.path(input, "embedding_float32.bin"), "rb")
  values <- readBin(connection, what = "numeric", n = expected_cells * scgpt_dims, size = 4, endian = "little")
  close(connection)
  if (length(values) != expected_cells * scgpt_dims) stop("Embedding binary is the wrong size")
  embedding <- matrix(values, nrow = expected_cells, ncol = scgpt_dims, byrow = TRUE)
  rm(values)
  if (any(!is.finite(embedding)) || max(abs(sqrt(rowSums(embedding^2)) - 1)) > 1e-4) {
    stop("Embedding rows are not unit length; the binary block was misread")
  }
  rownames(embedding) <- utils::read.csv(file.path(input, "metadata.csv"), colClasses = "character")$cell_barcode
  if (!setequal(rownames(embedding), barcodes)) stop("The embedding describes different cells")
  colnames(embedding) <- paste0("scgpt_", seq_len(scgpt_dims))
}


# --- Step 4: recluster each lineage --------------------------------------------------
recluster <- function(cells) {
  if (arm == "authors_pipeline") {
    # The authors' lines 233-271, unchanged.
    sub <- subset(integrated, cells = cells)
    DefaultAssay(sub) <- "integrated"
    sub <- FindVariableFeatures(sub, assay = "integrated", verbose = FALSE)
    set.seed(random_seed)
    sub <- RunPCA(sub, verbose = FALSE, seed.use = random_seed)
    reduction <- "pca"; dims <- 1:15
  } else {
    sub <- CreateSeuratObject(counts = GetAssayData(rna, assay = "RNA", layer = "counts")[, cells],
                              assay = "RNA", min.cells = 0, min.features = 0)
    sub[["scgpt"]] <- CreateDimReducObject(embeddings = embedding[cells, , drop = FALSE],
                                           key = "scgpt_", assay = "RNA")
    reduction <- "scgpt"; dims <- seq_len(scgpt_dims)
  }
  set.seed(random_seed)
  sub <- FindNeighbors(sub, reduction = reduction, dims = dims, verbose = FALSE)
  sub <- FindClusters(sub, resolution = 0.3, algorithm = 1, random.seed = random_seed, verbose = FALSE)
  set.seed(random_seed)
  sub <- RunUMAP(sub, reduction = reduction, dims = dims, n.neighbors = 30,
                 seed.use = random_seed, verbose = FALSE)
  coordinates <- Embeddings(sub, reduction = "umap")
  data.frame(cell_barcode = colnames(sub),
             subcluster = as.character(sub$seurat_clusters),
             UMAP_1 = coordinates[, 1], UMAP_2 = coordinates[, 2],
             stringsAsFactors = FALSE)
}

summary <- list()
for (lineage in lineages) {
  tag <- paste0(arm, "_", sub(" ", "_", tolower(lineage)))
  cells <- names(identity)[identity == lineage]
  message(sprintf("\n=== %s: %s, %d cells ===", arm, lineage, length(cells)))
  assignments <- recluster(cells)
  if (!setequal(assignments$cell_barcode, cells)) stop(tag, ": reclustering lost or gained cells")
  assignments$patient <- patients[match(assignments$cell_barcode, barcodes)]
  assignments$broad_cluster <- broad_cluster[match(assignments$cell_barcode, annotation$cell_barcode)]
  utils::write.csv(assignments, gzfile(file.path(output_dir, paste0("subclusters_", tag, ".csv.gz"))),
                   row.names = FALSE)

  # --- Step 5: describe each subcluster from the RNA ---------------------------------
  sub_rna <- subset(rna, cells = assignments$cell_barcode)
  sub_rna$subcluster <- assignments$subcluster[match(colnames(sub_rna), assignments$cell_barcode)]
  sub_rna[["umap"]] <- CreateDimReducObject(
    embeddings = as.matrix(data.frame(assignments[match(colnames(sub_rna), assignments$cell_barcode),
                                                  c("UMAP_1", "UMAP_2")], row.names = colnames(sub_rna))),
    key = "UMAP_", assay = "RNA")
  levels_sorted <- as.character(sort(unique(as.integer(sub_rna$subcluster))))
  Idents(sub_rna) <- factor(sub_rna$subcluster, levels = levels_sorted)

  set.seed(random_seed)
  markers <- FindAllMarkers(sub_rna, assay = "RNA", slot = "data", only.pos = TRUE,
                            min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
  utils::write.csv(markers, gzfile(file.path(output_dir, paste0("markers_", tag, ".csv.gz"))), row.names = FALSE)

  raw <- GetAssayData(sub_rna, assay = "RNA", layer = "counts")[panel, , drop = FALSE]
  lognorm <- GetAssayData(sub_rna, assay = "RNA", layer = "data")[panel, , drop = FALSE]
  evidence <- do.call(rbind, lapply(levels_sorted, function(cluster) {
    in_cluster <- sub_rna$subcluster == cluster
    donor_counts <- sort(table(sub_rna$patient[in_cluster]), decreasing = TRUE)
    top <- markers[markers$cluster == cluster & markers$p_val_adj < 0.05, ]
    top <- utils::head(top[order(-top$avg_log2FC), "gene"], TOP_N)
    row <- data.frame(arm = arm, lineage = lineage, subcluster = cluster,
                      n_cells = sum(in_cluster), n_patients = length(donor_counts),
                      largest_patient_share = round(100 * donor_counts[[1]] / sum(in_cluster), 1),
                      top_markers = paste(top, collapse = ", "), stringsAsFactors = FALSE)
    pct <- 100 * Matrix::rowSums(raw[, in_cluster, drop = FALSE] > 0) / sum(in_cluster)
    mean_expression <- Matrix::rowMeans(lognorm[, in_cluster, drop = FALSE])
    cbind(row, t(setNames(round(pct, 1), paste0("pct_", panel))),
          t(setNames(round(mean_expression, 3), paste0("mean_", panel))))
  }))
  utils::write.csv(evidence, file.path(output_dir, paste0("evidence_", tag, ".csv")), row.names = FALSE)

  # --- Step 6: the authors' plots for this lineage -----------------------------------
  arm_label <- if (arm == "authors_pipeline") "Authors' pipeline" else "scGPT raw"
  cluster_plot <- DimPlot(sub_rna, reduction = "umap", label = TRUE, pt.size = 0.6, shuffle = TRUE) +
    ggtitle(sprintf("%s — %s reclustered (%d subclusters)", arm_label, lineage, length(levels_sorted))) +
    NoLegend()
  ggsave(file.path(figure_dir, paste0("20_subclusters_", tag, ".png")), cluster_plot,
         width = 7, height = 6, dpi = 150)
  feature_plots <- lapply(panel, function(gene) {
    FeaturePlot(sub_rna, features = gene, reduction = "umap", pt.size = 0.6,
                cols = c("grey88", "deeppink")) +
      theme(plot.title = element_text(size = 11, face = "bold"),
            axis.title = element_text(size = 7), axis.text = element_text(size = 6))
  })
  ggsave(file.path(figure_dir, paste0("20_featureplots_", tag, ".png")),
         wrap_plots(feature_plots, ncol = 4) +
           plot_annotation(title = sprintf("%s — %s: the authors' reclustering panel", arm_label, lineage)),
         width = 18, height = 8.5, dpi = 130)

  summary[[lineage]] <- list(n_cells = length(cells), n_subclusters = length(levels_sorted),
                             sizes = as.list(table(assignments$subcluster)))
  message(sprintf("%s: %d subclusters", lineage, length(levels_sorted)))
  print(evidence[, c("subcluster", "n_cells", "n_patients", "largest_patient_share",
                     paste0("pct_", panel))], row.names = FALSE)
}

writeLines(toJSON(list(
  stage = "20_subcluster_lineages", arm = arm,
  follows = "05_REFERENCES/author_code/figure1d_author_workflow.R lines 211-275",
  representation = if (arm == "authors_pipeline") "integrated assay: variable features + PCA recomputed per lineage, PCs 1-15"
                   else "scGPT embedding, all 512 dims (no assay to recompute features on)",
  clustering = list(resolution = 0.3, algorithm = "Louvain", umap_neighbors = 30, seed = random_seed),
  expression = "RNA, LogNormalize 10000, Stage 13 fingerprint checked",
  panel = panel, lineages = summary,
  seurat_version = as.character(packageVersion("Seurat"))
), auto_unbox = TRUE, pretty = TRUE), file.path(output_dir, paste0("subcluster_audit_", arm, ".json")))

message("\nStage 20 complete for ", arm)
