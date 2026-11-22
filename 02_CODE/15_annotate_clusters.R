# Stage 15: assign broad cell type identities, by hand, the way the authors did.
#
# The authors did not compute their annotation. They plotted each canonical marker
# over the UMAP, looked at the result, and typed cluster numbers into a block of
# assignments (their lines 181-202). This stage is that block, written once per
# partition, and nothing more. The main comparison is the authors' pipeline against
# uncorrected scGPT; the two Harmony blocks are the earlier comparison, kept as they
# were so its outputs can still be reproduced.
#
# It is deliberately not a rule engine. An earlier version of this stage inferred
# identities from marker thresholds and invented finer categories than the panel can
# support; it was deleted. The categories here are exactly the authors' seven, and the
# assignments below are entered by a person reading Stage 14's figures.
#
# The cluster numbers are not transferable. Louvain numbers by descending size, so the
# authors' "cluster 7 = NK Cells" says nothing about cluster 7 in either partition
# here, and PCA cluster 3 has nothing to do with scGPT cluster 3. Each block below was
# read off that partition's own plots.

suppressPackageStartupMessages({
  library(jsonlite)
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: 15_annotate_clusters.R <frozen_dir> <evidence_dir> <benchmark_dir> <output_dir> <figure_dir>")
}
frozen_dir <- args[1]
evidence_dir <- args[2]
benchmark_dir <- args[3]
output_dir <- args[4]
figure_dir <- args[5]

expected_cells <- 70391

# The authors' seven categories, from their lines 181-202. Nothing finer: CD4 against
# CD8 against Treg, and monocyte subtypes, come from reclustering each lineage, which
# has not happened yet. "Undetermined" is theirs too -- they used it rather than
# forcing a label onto a cluster the panel could not resolve.
VALID_IDENTITIES <- c(
  "T Cells", "NK Cells", "Monocytes", "DC",
  "B Cells", "Plasma Cells", "Undetermined"
)

# The authors' Figure 1D colours (their lines 412-416), taken for the categories this
# pass can actually assign. Their vector is built for the twelve post-reclustering
# labels, so the four T subtypes and three monocyte subtypes each collapse here; the
# first colour of each run is used, and NK, DC and B keep theirs exactly. Reusing
# their palette now means the broad plot and the eventual Figure 1D read as one family.
IDENTITY_LEVELS <- c("T Cells", "NK Cells", "DC", "Monocytes", "B Cells",
                     "Plasma Cells", "Undetermined")
IDENTITY_COLOURS <- c(
  "T Cells"      = "darkturquoise",
  "NK Cells"     = "lightpink",
  "DC"           = "hotpink",
  "Monocytes"    = "gold1",
  "B Cells"      = "darkorange",
  "Plasma Cells" = "darkviolet",
  "Undetermined" = "gray"
)

# --- The manual annotation -------------------------------------------------------
# Read off Stage 14's feature plots and panel tables. Each entry is a decision, not a
# computation; the evidence for each is recorded in annotation_rationale.csv.

# PCA-15 Harmony: 12 clusters, built from SCT residuals. Cluster 3 (APOE, C1Q genes,
# LYVE1, CD14 and CD68 high) is macrophage-like; it goes under Monocytes because the
# authors' seven categories have nothing closer. Cluster 8 is dominated by XIST and
# long-transcript genes (SYNE2, MACF1, UTRN) with no canonical lineage marker above
# 50%; it spans 58 patients, so it is not one donor's cells, and the panel cannot
# resolve it. The authors used Undetermined for exactly that case.
pca15_annotation <- list(
  "T Cells"      = c("0", "1", "2", "4", "10"),
  "Monocytes"    = c("3", "7"),
  "DC"           = c("5", "9"),
  "NK Cells"     = c("6"),
  "B Cells"      = c("11"),
  "Undetermined" = c("8")
)

# scGPT Harmony: 11 clusters. Cluster 7 holds the same cells as PCA cluster 8 (85% of
# that cluster lands here) with the same XIST and long-gene markers, so it gets the
# same label.
scgpt_annotation <- list(
  "T Cells"      = c("0", "1", "4", "10"),
  "Monocytes"    = c("2", "6"),
  "DC"           = c("3", "8"),
  "NK Cells"     = c("5"),
  "B Cells"      = c("9"),
  "Undetermined" = c("7")
)

# Authors' pipeline (Stage 03): 12 clusters, read from Stage 14's evidence for this
# partition alone, without SingleR, which Stage 17 scores the annotation against.
# - Cluster 3 (APOE, C1Q genes, LYVE1; CD14 and CD68 high) is macrophage-like and goes
#   under Monocytes, as PCA cluster 3 did above.
# - Cluster 10 is plasmacytoid DC (LRRC26, PTCRA, SCT, PTPRS, TCL1A). Its JCHAIN is
#   high, which the authors' panel reads as plasma, but pDCs express JCHAIN too, and
#   there is no plasma programme among its top markers. Stage 16 also files pDCs
#   under DC. About 11% of its cells express CD19: this partition has no separate B
#   cell cluster, and its B cells appear to sit here. That is left visible as a real
#   difference between the arms, not corrected.
# - Cluster 11 is the XIST and long-transcript cluster seen in every arm, with no
#   canonical marker above 50% across 57 patients: Undetermined.
authors_pipeline_annotation <- list(
  "T Cells"      = c("0", "1", "2", "4", "8", "9"),
  "Monocytes"    = c("3", "7"),
  "DC"           = c("5", "10"),
  "NK Cells"     = c("6"),
  "Undetermined" = c("11")
)

# scGPT, uncorrected: 12 clusters, read the same way.
# - Cluster 2 is the same macrophage-like population: Monocytes.
# - Cluster 9 is plasmacytoid DC, on the same markers as the authors' cluster 10.
# - Cluster 10 (CD79A, CD19, VPREB3, FCRLA) is B cells.
# - Cluster 7 is naive/central memory T (LEF1, IL7R, SELL), but 47% of it comes from
#   one patient: uncorrected scGPT kept part of that patient's cells together. The
#   label follows the biology; the patient concentration is reported with Stage 17.
# - Cluster 11 (103 cells) is cycling (RRM2, UBE2C, BIRC5, TYMS) with CD3D 71%: T Cells.
# - Cluster 8 is the XIST and long-transcript cluster: Undetermined.
scgpt_raw_annotation <- list(
  "T Cells"      = c("0", "1", "4", "7", "11"),
  "Monocytes"    = c("2", "6"),
  "DC"           = c("3", "9"),
  "NK Cells"     = c("5"),
  "B Cells"      = c("10"),
  "Undetermined" = c("8")
)

annotations <- list(
  authors_pipeline_cluster = authors_pipeline_annotation,
  scgpt_raw_cluster = scgpt_raw_annotation,
  pca15_harmony_cluster = pca15_annotation,
  scgpt_harmony_cluster = scgpt_annotation
)

# --- Apply and check --------------------------------------------------------------
# The checks below are what stops a typed cluster number from silently going missing:
# every cluster must receive exactly one identity, and every identity must be one of
# the authors' seven.
frozen <- utils::read.csv(file.path(frozen_dir, "frozen_cluster_labels.csv.gz"),
                          stringsAsFactors = FALSE, colClasses = "character")
if (nrow(frozen) != expected_cells) {
  stop("Frozen labels cover ", nrow(frozen), " cells, expected ", expected_cells)
}

summaries <- list()
for (partition in names(annotations)) {
  assignment <- annotations[[partition]]

  bad <- setdiff(names(assignment), VALID_IDENTITIES)
  if (length(bad) > 0) {
    stop("Not one of the authors' categories: ", paste(bad, collapse = ", "))
  }

  clusters <- sort(unique(frozen[[partition]]), method = "radix")
  assigned <- unlist(assignment, use.names = FALSE)

  duplicated_clusters <- assigned[duplicated(assigned)]
  if (length(duplicated_clusters) > 0) {
    stop(partition, ": cluster(s) assigned twice: ",
         paste(unique(duplicated_clusters), collapse = ", "))
  }
  unassigned <- setdiff(clusters, assigned)
  if (length(unassigned) > 0) {
    stop(partition, ": no identity given for cluster(s): ",
         paste(unassigned, collapse = ", "),
         "\n  Every cluster needs a decision; use \"Undetermined\" if the panel cannot resolve it.")
  }
  unknown <- setdiff(assigned, clusters)
  if (length(unknown) > 0) {
    stop(partition, ": annotation names cluster(s) that do not exist: ",
         paste(unknown, collapse = ", "))
  }

  lookup <- stats::setNames(
    rep(names(assignment), lengths(assignment)),
    unlist(assignment, use.names = FALSE)
  )
  frozen[[paste0(partition, "_identity")]] <- unname(lookup[frozen[[partition]]])

  counts <- table(frozen[[paste0(partition, "_identity")]])
  summaries[[partition]] <- list(
    n_clusters = length(clusters),
    identities = as.list(stats::setNames(as.integer(counts), names(counts)))
  )
  message("\n", partition, " (", length(clusters), " clusters):")
  for (identity in names(counts)) {
    which_clusters <- assignment[[identity]]
    message(sprintf("  %-13s %6d cells   clusters %s",
                    identity, counts[[identity]], paste(which_clusters, collapse = ", ")))
  }
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
utils::write.csv(frozen, gzfile(file.path(output_dir, "annotated_cluster_labels.csv.gz")),
                 row.names = FALSE)

# The evidence each decision rests on, carried alongside the decision itself.
evidence <- utils::read.csv(file.path(evidence_dir, "cluster_evidence_table.csv"),
                            stringsAsFactors = FALSE, colClasses = c(cluster = "character"))
evidence$assigned_identity <- vapply(seq_len(nrow(evidence)), function(i) {
  assignment <- annotations[[evidence$partition[i]]]
  hit <- names(assignment)[vapply(assignment, function(x) evidence$cluster[i] %in% x, logical(1))]
  if (length(hit) == 1) hit else NA_character_
}, character(1))
utils::write.csv(evidence, file.path(output_dir, "annotation_rationale.csv"), row.names = FALSE)

# --- The labelled UMAP, one per arm ----------------------------------------------
# The authors' identity DimPlot (their lines 146-153) with cluster_ident in place of
# seurat_clusters. Only coordinates and metadata are needed to draw it, so rather than
# rebuild the 19,689 x 70,391 expression object this holds a one-row placeholder assay:
# DimPlot reads the reduction and the metadata and never touches expression.
message("\nDrawing labelled UMAPs...")
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

# Each partition is drawn on the layout Stage 14 used for it: the Harmony arms on
# Stage 06's, the authors' pipeline on Stage 03's integrated UMAP, and uncorrected
# scGPT on Stage 04's UMAP of the same embedding. Every file is joined by barcode.
read_layout <- function(path, columns, prefix) {
  table <- utils::read.csv(path, stringsAsFactors = FALSE)
  if (!setequal(table$cell_barcode, frozen$cell_barcode)) {
    stop(basename(path), " does not describe the same cells as the frozen labels")
  }
  rownames(table) <- table$cell_barcode
  table <- table[frozen$cell_barcode, columns, drop = FALSE]
  colnames(table) <- paste0(prefix, c("_UMAP_1", "_UMAP_2"))
  table
}
results_dir <- dirname(benchmark_dir)
benchmark_file <- file.path(benchmark_dir, "benchmark_cluster_assignments.csv.gz")
layouts <- cbind(
  read_layout(benchmark_file, c("pca_harmony_UMAP_1", "pca_harmony_UMAP_2"), "pca_harmony"),
  read_layout(benchmark_file, c("scgpt_harmony_UMAP_1", "scgpt_harmony_UMAP_2"), "scgpt_harmony"),
  read_layout(file.path(results_dir, "03_seurat_reproduction", "GSE200164_umap_clusters.csv.gz"),
              c("UMAP_1", "UMAP_2"), "authors_pipeline"),
  read_layout(file.path(results_dir, "04_scgpt_embedding", "GSE200164_scgpt_umap_clusters.csv.gz"),
              c("UMAP_1", "UMAP_2"), "scgpt_raw")
)

placeholder <- Matrix::Matrix(0, nrow = 1, ncol = nrow(frozen), sparse = TRUE)
rownames(placeholder) <- "placeholder"
colnames(placeholder) <- frozen$cell_barcode
csf <- CreateSeuratObject(counts = placeholder, assay = "RNA")
umap_prefix <- c(authors_pipeline_cluster = "authors_pipeline",
                 scgpt_raw_cluster = "scgpt_raw",
                 pca15_harmony_cluster = "pca_harmony",
                 scgpt_harmony_cluster = "scgpt_harmony")
partition_labels <- c(authors_pipeline_cluster = "Authors' pipeline (Stage 03)",
                      scgpt_raw_cluster = "scGPT, uncorrected",
                      pca15_harmony_cluster = "PCA-15 Harmony",
                      scgpt_harmony_cluster = "scGPT Harmony")

for (partition in names(annotations)) {
  identity_column <- paste0(partition, "_identity")
  csf@meta.data[[identity_column]] <- factor(
    frozen[[identity_column]],
    levels = intersect(IDENTITY_LEVELS, unique(frozen[[identity_column]]))
  )

  prefix <- umap_prefix[[partition]]
  coordinates <- as.matrix(layouts[, paste0(prefix, c("_UMAP_1", "_UMAP_2"))])
  colnames(coordinates) <- paste0(prefix, "UMAP_", 1:2)
  rownames(coordinates) <- colnames(csf)
  csf[[paste0("umap_", prefix)]] <- CreateDimReducObject(
    embeddings = coordinates, key = paste0(gsub("_", "", prefix), "UMAP_"), assay = "RNA"
  )

  Idents(csf) <- csf@meta.data[[identity_column]]
  plot <- DimPlot(csf, reduction = paste0("umap_", prefix),
                  cols = IDENTITY_COLOURS[levels(Idents(csf))],
                  label = TRUE, label.size = 4, repel = TRUE,
                  pt.size = 1, shuffle = TRUE) +
    ggtitle(paste0(partition_labels[[partition]], " — broad cell types"))
  target <- file.path(figure_dir, paste0("15_celltypes_", partition, ".png"))
  ggsave(target, plot, width = 7.5, height = 6, dpi = 150)
  if (!file.exists(target) || file.info(target)$size == 0) {
    stop("Labelled UMAP wrote nothing for ", partition)
  }
  message("  ", target)
}

writeLines(jsonlite::toJSON(list(
  stage = "15_annotate_clusters",
  method = "manual, following 05_REFERENCES/author_code/figure1d_author_workflow.R lines 181-202",
  categories = VALID_IDENTITIES,
  annotations = annotations,
  summaries = summaries
), auto_unbox = TRUE, pretty = TRUE),
file.path(output_dir, "annotation_audit.json"))

message("\nStage 15 complete.")
