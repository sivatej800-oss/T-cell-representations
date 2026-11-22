# Stage 17: does cluster structure follow biology or the patient?
#
# Stage 06 measured how strongly each clustering tracks which patient a cell came from
# (donor NMI). That is only half of the question: a clustering can be mostly biology and
# still carry some patient signal, or the reverse. This stage measures the other half --
# how strongly each clustering tracks cell type -- using labels that belong to neither
# arm: SingleR's per-cell calls against purified immune profiles (Stage 16). It covers the
# PCA (15 PCs) and scGPT arms, each uncorrected and Harmony-corrected.
#
# Three things are measured.
#   A. Patient NMI and biology NMI side by side, per arm.
#   B. How closely each Harmony arm's hand annotation (Stage 15) agrees with SingleR, cell
#      by cell, with the monocyte/DC boundary kept apart and then merged into one group.
#   C. For the Monocyte/DC disagreements specifically, what marker genes say about
#      which side is right.
#
# SingleR's "unassigned" cells (pruned label NA) are excluded from every biology
# comparison. Patient NMI uses all cells, because patient identity is known for every
# cell and the check against Stage 06's value is on the full set.
#
# The expression is Stage 13/16's: SoupX-corrected counts, fractional on purpose,
# log-normalized once. Nothing is rounded, cast or filtered.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(Seurat)
  library(Matrix)
  library(SingleR)
  library(celldex)
  library(BiocParallel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9) {
  stop("Usage: 17_cluster_structure_vs_singler.R <seurat_input_dir> <cluster_csv> ",
       "<donor_leakage_csv> <singler_csv> <old_stage15_csv> <annotation_csv> <output_dir> ",
       "<workers> <figure_dir>")
}
input_dir <- args[1]
cluster_csv <- args[2]
leakage_csv <- args[3]
singler_csv <- args[4]
old_stage15_csv <- args[5]
annotation_csv <- args[6]
output_dir <- args[7]
n_workers <- as.integer(args[8])
figure_dir <- args[9]

expected_cells <- 70391
expected_genes <- 19689

# The RNA data-layer sum Stage 13 recorded (and Stage 16 reproduced). Stage 17 reads the
# same unchanged seurat_input/, so it must reproduce this number.
stage13_data_sum <- 176503092.94553709
sum_tolerance <- 1e-9

# --- The arms ---------------------------------------------------------------------
# One entry per clustering. `cluster_column` is the arm's column in Stage 06's cluster file.
# An arm with an `annotation_column` has a hand annotation in Stage 15's file (and the
# `stage15_cluster_column` that annotation was made on, checked against Stage 06's); an arm
# without one gets the NMI measurements only. `corrected_from` names the uncorrected arm a
# Harmony arm was derived from (it draws the arrow); `family` gives each representation one
# colour; `umap_prefix` names the arm's own UMAP columns. To add an arm, add an entry.
#
# pca_harmony is the main 15-PC arm, as in the authors' dims = 1:15. pca_harmony_30dims is
# deliberately not here.
arms <- list(
  list(name = "pca_raw", label = "PCA raw", family = "PCA", corrected_from = NA_character_,
       cluster_column = "pca_raw", annotation_column = NULL),
  list(name = "scgpt_raw", label = "scGPT raw", family = "scGPT", corrected_from = NA_character_,
       cluster_column = "scgpt_raw", annotation_column = NULL),
  list(name = "pca_harmony", label = "PCA + Harmony", family = "PCA", corrected_from = "pca_raw",
       cluster_column = "pca_harmony", annotation_column = "pca15_harmony_cluster_identity",
       stage15_cluster_column = "pca15_harmony_cluster", umap_prefix = "pca_harmony"),
  list(name = "scgpt_harmony", label = "scGPT + Harmony", family = "scGPT",
       corrected_from = "scgpt_raw",
       cluster_column = "scgpt_harmony", annotation_column = "scgpt_harmony_cluster_identity",
       stage15_cluster_column = "scgpt_harmony_cluster", umap_prefix = "scgpt_harmony")
)

# Marker genes for telling monocytes from dendritic cells: classical/non-classical
# monocyte (CD14, VCAN, S100A8, FCGR3A), cDC2 and general DC (FCER1A, CD1C, CLEC10A),
# plasmacytoid DC (LILRA4).
marker_genes <- c("CD14", "VCAN", "S100A8", "FCGR3A", "FCER1A", "CD1C", "CLEC10A", "LILRA4")

# Copied from 02_CODE/06_benchmark_representations.R so this stage runs on its own.
# NMI: how much knowing a cell's cluster tells you about its other label, scaled so 1 is
# complete determination and 0 is none.
normalized_mutual_information <- function(labels_a, labels_b) {
  joint <- table(labels_a, labels_b) / length(labels_a)
  margin_a <- rowSums(joint)
  margin_b <- colSums(joint)
  entropy <- function(p) {
    p <- p[p > 0]
    -sum(p * log(p))
  }
  outer_product <- outer(margin_a, margin_b)
  positive <- joint > 0
  mutual <- sum(joint[positive] * log(joint[positive] / outer_product[positive]))
  denominator <- sqrt(entropy(margin_a) * entropy(margin_b))
  if (denominator == 0) return(NA_real_)
  mutual / denominator
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
checks <- list()

# --- Step 1: rebuild Stage 13/16's expression exactly -------------------------------
message("Reading the RNA matrix...")
counts <- as(Matrix::readMM(file.path(input_dir, "counts.mtx")), "CsparseMatrix")
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

# The marker genes are read from the SoupX-corrected counts themselves, not the
# normalized layer, so "expressing" means a corrected count above zero.
missing_markers <- setdiff(marker_genes, genes)
if (length(missing_markers) > 0) {
  stop("Marker genes absent from the matrix: ", paste(missing_markers, collapse = ", "))
}
marker_counts <- counts[marker_genes, , drop = FALSE]
checks$marker_genes_present <- length(marker_genes)

so <- CreateSeuratObject(counts = counts, assay = "RNA", project = "GSE200164_stage17",
                         min.cells = 0, min.features = 0)
rm(counts)
gc(verbose = FALSE)
so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 10000,
                    verbose = FALSE)
lognorm <- GetAssayData(so, assay = "RNA", layer = "data")
observed_sum <- sum(lognorm)
relative_error <- abs(observed_sum - stage13_data_sum) / stage13_data_sum
if (relative_error > sum_tolerance) {
  stop(sprintf("Data-layer sum %.8f differs from Stage 13's %.8f; relative error %.3g",
               observed_sum, stage13_data_sum, relative_error))
}
message(sprintf("Expression: %d genes x %d cells, data sum %.6f (relative error %.2e)",
                nrow(lognorm), ncol(lognorm), observed_sum, relative_error))
checks$data_sum_fingerprint <- list(expected = stage13_data_sum, observed = observed_sum,
                                    relative_error = relative_error,
                                    tolerance = sum_tolerance)
seurat_version <- as.character(packageVersion("Seurat"))
rm(so)
gc(verbose = FALSE)

# --- Step 2: load every label source and align it to the matrix by barcode -----------
# Joined by barcode, never position. A positional join on mismatched files would attach
# one cell's cluster to another cell's expression and produce a clean, wrong result.
align_to_matrix <- function(table, barcode_column, source) {
  if (nrow(table) != expected_cells || anyDuplicated(table[[barcode_column]])) {
    stop(source, ": expected ", expected_cells, " unique barcodes")
  }
  if (!setequal(table[[barcode_column]], barcodes)) {
    stop(source, ": barcodes differ from seurat_input/barcodes.tsv")
  }
  table[match(barcodes, table[[barcode_column]]), , drop = FALSE]
}

clusters <- align_to_matrix(
  utils::read.csv(cluster_csv, colClasses = "character"), "cell_barcode", "cluster assignments")
singler <- align_to_matrix(
  utils::read.csv(singler_csv, colClasses = "character"), "cell_barcode", "Stage 16 labels")
old15 <- align_to_matrix(
  utils::read.csv(old_stage15_csv, colClasses = "character"), "cell_barcode", "old Stage 15")
annotation <- align_to_matrix(
  utils::read.csv(annotation_csv, colClasses = "character"), "cell_barcode", "Stage 15 annotation")
stopifnot(identical(clusters$cell_barcode, barcodes), identical(singler$cell_barcode, barcodes),
          identical(old15$cell_barcode, barcodes), identical(annotation$cell_barcode, barcodes))
patients <- clusters$internal_sample_id
if (anyNA(patients)) stop("A cell has no patient identifier")

assigned <- singler$broad_label != "unassigned"
n_excluded <- sum(!assigned)
if (n_excluded != sum(singler$is_unassigned == "TRUE")) {
  stop("broad_label 'unassigned' and is_unassigned disagree")
}
message(sprintf("SingleR-unassigned cells excluded from biology comparisons: %d", n_excluded))
checks$cells <- list(total = length(barcodes), singler_unassigned_excluded = n_excluded,
                     assigned = sum(assigned), n_patients = length(unique(patients)))

# --- Step 3: broad-level SingleR confidence ------------------------------------------
# Stage 16's fine-level run supplies the label. This second run, against Monaco's ten
# main types, is used only for a confidence score: delta.next is the gap between a cell's
# best-scoring type and its runner-up, and at the main level a small gap means the cell
# sits between two lineages rather than between two subsets of one.
message("Loading celldex::MonacoImmuneData()...")
ref <- celldex::MonacoImmuneData()
shared_genes <- intersect(rownames(lognorm), rownames(ref))
message(sprintf("Genes shared with the reference: %d of %d", length(shared_genes), nrow(lognorm)))
test_matrix <- lognorm[shared_genes, , drop = FALSE]
ref_matrix <- ref[shared_genes, ]
rm(lognorm)
gc(verbose = FALSE)

message(sprintf("Running main-level SingleR on %d cells with %d workers...",
                ncol(test_matrix), n_workers))
set.seed(42)
timing <- system.time(
  main_run <- SingleR(test = test_matrix, ref = ref_matrix, labels = ref_matrix$label.main,
                      BPPARAM = BiocParallel::MulticoreParam(workers = n_workers))
)
message(sprintf("Main-level SingleR wall time %.1f s", unname(timing["elapsed"])))
if (!identical(rownames(main_run), barcodes)) {
  stop("SingleR output is not in the matrix's cell order")
}
rm(test_matrix)
gc(verbose = FALSE)

confidence <- data.frame(
  cell_barcode = barcodes,
  broad_label = singler$broad_label,
  stage16_label_main = singler$label_main,
  main_run_pruned_label = main_run$pruned.labels,
  main_run_unpruned_label = main_run$labels,
  delta_next_main = main_run$delta.next,
  stringsAsFactors = FALSE
)
utils::write.csv(confidence, gzfile(file.path(output_dir, "broad_level_confidence.csv.gz")),
                 row.names = FALSE)

# Median confidence per Stage 16 broad label, assigned cells only.
conf_assigned <- confidence[assigned, ]
confidence_by_label <- do.call(rbind, lapply(sort(unique(conf_assigned$broad_label)), function(label) {
  values <- conf_assigned$delta_next_main[conf_assigned$broad_label == label]
  data.frame(broad_label = label, n_cells = length(values),
             median_delta_next_main = median(values), stringsAsFactors = FALSE)
}))
utils::write.csv(confidence_by_label, file.path(output_dir, "confidence_by_broad_label.csv"),
                 row.names = FALSE)
message("\nMedian main-level delta_next per broad label:"); print(confidence_by_label, row.names = FALSE)

# Where the main-level run and Stage 16's fine-derived main label disagree. Cells
# unassigned in either run cannot disagree on a label they do not have, so they are counted
# apart rather than folded in.
both_assigned <- !is.na(main_run$pruned.labels) & assigned
disagree <- both_assigned & main_run$pruned.labels != singler$label_main
checks$main_vs_fine_derived_main <- list(
  cells_assigned_in_both = sum(both_assigned),
  disagree = sum(disagree),
  pct_disagree_of_both_assigned = 100 * sum(disagree) / sum(both_assigned),
  main_run_unassigned = sum(is.na(main_run$pruned.labels)),
  stage16_unassigned = n_excluded,
  main_run_unassigned_among_stage16_assigned = sum(is.na(main_run$pruned.labels) & assigned)
)
message(sprintf("\nMain-run vs fine-derived main label: %d of %d cells assigned in both disagree (%.2f%%)",
                sum(disagree), sum(both_assigned), 100 * sum(disagree) / sum(both_assigned)))

# --- Step 4: each annotated arm's per-cell identity, from Stage 15 --------------------
# Stage 15's file is the only source of the annotation: this stage no longer carries its own
# copy of the cluster-to-identity tables. It is checked, not trusted: Stage 15 must have
# annotated the same partitions Stage 06 produced, and every cluster must carry exactly one
# identity. "Undetermined" is Stage 15's label for a cluster the marker panel could not
# resolve (XIST and long-transcript genes, no lineage marker).
VALID_IDENTITIES <- c("T Cells", "NK Cells", "DC", "Monocytes", "B Cells", "Undetermined")

identity_by_arm <- list()
annotation_rows <- list()
for (arm in arms) {
  if (is.null(arm$annotation_column)) next
  if (!identical(annotation[[arm$stage15_cluster_column]], clusters[[arm$cluster_column]])) {
    stop(arm$name, ": Stage 15's clusters differ from Stage 06's, cell by cell")
  }
  identity <- annotation[[arm$annotation_column]]
  if (anyNA(identity) || any(identity == "")) stop(arm$name, ": a cell has no identity")
  if (!all(identity %in% VALID_IDENTITIES)) {
    stop(arm$name, ": unexpected identities: ",
         paste(setdiff(identity, VALID_IDENTITIES), collapse = ", "))
  }
  identities_per_cluster <- tapply(identity, clusters[[arm$cluster_column]],
                                   function(x) length(unique(x)))
  if (any(identities_per_cluster != 1)) stop(arm$name, ": a cluster carries two identities")
  identity_by_arm[[arm$name]] <- identity
  pairs <- unique(data.frame(cluster = clusters[[arm$cluster_column]], identity = identity,
                             stringsAsFactors = FALSE))
  annotation_rows[[arm$name]] <- cbind(arm = arm$name, pairs[order(as.integer(pairs$cluster)), ])
  message(sprintf("%s identities: %s", arm$name,
                  paste(names(table(identity)), table(identity), sep = " = ", collapse = ", ")))
}
annotation_table <- do.call(rbind, annotation_rows)
utils::write.csv(annotation_table, file.path(output_dir, "annotation_table_used.csv"),
                 row.names = FALSE)

# Merged-myeloid resolution. The monocyte/DC boundary is where SingleR is least stable
# (Step 3), so every comparison is repeated with DC and Monocytes folded into one group,
# in both SingleR's labels and the annotation, to see whether that boundary alone decides
# the result. Every label either side can carry appears here, so none can silently miss.
merge_table <- data.frame(
  label  = c("T Cells", "NK Cells", "DC",      "Monocytes", "B Cells", "Plasma Cells",
             "Progenitors", "Undetermined", "unassigned"),
  merged = c("T Cells", "NK Cells", "Myeloid", "Myeloid",   "B Cells", "Plasma Cells",
             "Progenitors", "Undetermined", "unassigned"),
  stringsAsFactors = FALSE
)
merge_labels <- function(labels) {
  merged <- merge_table$merged[match(labels, merge_table$label)]
  if (anyNA(merged)) {
    stop("Labels missing from the merge table: ", paste(unique(labels[is.na(merged)]), collapse = ", "))
  }
  merged
}
utils::write.csv(merge_table, file.path(output_dir, "myeloid_merge_table.csv"), row.names = FALSE)

# The rebuild's annotation against the superseded Stage 15, where the old file has the
# same column. Reported, not asserted: the rebuilt Stage 15 added the Undetermined group.
old_match <- list()
for (arm in arms) {
  if (is.null(arm$annotation_column)) next
  old_identity <- old15[[arm$annotation_column]]
  if (is.null(old_identity)) next
  old_match[[arm$name]] <- list(pct_cells_same_identity_as_old_stage15 =
                                  100 * mean(identity_by_arm[[arm$name]] == old_identity))
  message(sprintf("%s: %.2f%% of cells keep their superseded Stage 15 identity", arm$name,
                  old_match[[arm$name]]$pct_cells_same_identity_as_old_stage15))
}
checks$identity_vs_old_stage15 <- old_match
checks$identity_counts <- lapply(identity_by_arm, function(x) as.list(table(x)))

# --- Step 5: A -- patient NMI and biology NMI per arm ---------------------------------
# The patient NMI is recomputed here on all cells and must equal the value Stage 06 wrote;
# agreement shows this stage reads the same clusters and the same patients as Stage 06.
# Biology NMI uses assigned cells only, at full resolution and with Myeloid merged.
leakage <- utils::read.csv(leakage_csv, stringsAsFactors = FALSE)
expected_patient_nmi <- c(scgpt_raw = 0.0721503826775315, scgpt_harmony = 0.0506991838687104,
                          pca_raw = 0.0489403923224483, pca_harmony = 0.0380828878753405)

nmi_rows <- list()
for (arm in arms) {
  cluster_ids <- clusters[[arm$cluster_column]]
  if (anyNA(cluster_ids)) stop(arm$name, ": a cell has no cluster")
  patient_nmi <- normalized_mutual_information(cluster_ids, patients)
  recorded <- leakage$donor_nmi[leakage$arm == arm$name]
  if (length(recorded) != 1 || abs(patient_nmi - recorded) > 1e-12) {
    stop(arm$name, ": patient NMI ", patient_nmi, " differs from Stage 06's ", recorded)
  }
  if (!arm$name %in% names(expected_patient_nmi) ||
      abs(patient_nmi - expected_patient_nmi[[arm$name]]) > 1e-12) {
    stop(arm$name, ": patient NMI differs from the expected value")
  }
  nmi_rows[[arm$name]] <- data.frame(
    arm = arm$name,
    n_clusters = length(unique(cluster_ids)),
    patient_nmi = patient_nmi,
    biology_nmi = normalized_mutual_information(cluster_ids[assigned], singler$broad_label[assigned]),
    n_cells_scored = sum(assigned),
    biology_nmi_merged_myeloid = normalized_mutual_information(
      cluster_ids[assigned], merge_labels(singler$broad_label[assigned])),
    stringsAsFactors = FALSE
  )
}
nmi_table <- do.call(rbind, nmi_rows)
utils::write.csv(nmi_table, file.path(output_dir, "nmi_by_arm.csv"), row.names = FALSE)
message("\n--- A: patient vs biology NMI ---"); print(nmi_table, row.names = FALSE)
checks$patient_nmi_matches_stage06 <- TRUE

# --- Step 6: B -- hand annotation vs SingleR ------------------------------------------
# Agreement is an exact label match on cells SingleR assigned. Cells in an "Undetermined"
# cluster have no annotation to agree or disagree with, so they are left out of the
# agreement figures and counted on their own; their SingleR labels still appear in the
# confusion table, as a row of their own. SingleR's Plasma Cells and Progenitors have no
# counterpart in the annotation, so those cells can only disagree; they are listed apart.
label_order <- c("T Cells", "NK Cells", "DC", "Monocytes", "Myeloid", "B Cells", "Undetermined",
                 "Plasma Cells", "Progenitors")
agreement_for <- function(arm_name, resolution) {
  mine <- identity_by_arm[[arm_name]][assigned]
  theirs <- singler$broad_label[assigned]
  if (resolution == "merged_myeloid") {
    mine <- merge_labels(mine)
    theirs <- merge_labels(theirs)
  }
  scored <- mine != "Undetermined"
  same <- mine[scored] == theirs[scored]

  confusion <- as.data.frame(table(my_identity = mine, singler_broad_label = theirs),
                             stringsAsFactors = FALSE)
  confusion <- confusion[confusion$Freq > 0, ]
  names(confusion)[3] <- "n_cells"
  confusion$row_pct <- 100 * confusion$n_cells / as.integer(table(mine)[confusion$my_identity])
  confusion$arm <- arm_name
  confusion$resolution <- resolution
  confusion <- confusion[order(match(confusion$my_identity, label_order), -confusion$n_cells),
                         c("arm", "resolution", "my_identity", "singler_broad_label", "n_cells", "row_pct")]

  types <- label_order[label_order %in% setdiff(unique(mine), "Undetermined")]
  per_type <- do.call(rbind, lapply(types, function(type) {
    in_type <- scored & mine == type
    data.frame(arm = arm_name, resolution = resolution, my_identity = type,
               n_cells = sum(in_type), n_agree = sum(in_type & theirs == type),
               pct_agree = 100 * sum(in_type & theirs == type) / sum(in_type),
               stringsAsFactors = FALSE)
  }))
  overall <- data.frame(
    arm = arm_name, resolution = resolution,
    n_cells_scored = sum(scored), n_agree = sum(same), pct_agree = 100 * mean(same),
    n_undetermined_excluded = sum(!scored),
    pct_agree_if_undetermined_counted_as_disagreement = 100 * sum(same) / length(mine),
    singler_plasma_cells = sum(theirs[scored] == "Plasma Cells"),
    singler_progenitors = sum(theirs[scored] == "Progenitors"),
    stringsAsFactors = FALSE
  )
  list(overall = overall, per_type = per_type, confusion = confusion)
}

agreement_results <- list()
overall_rows <- list()
for (arm in arms) {
  if (is.null(arm$annotation_column)) next
  for (resolution in c("full", "merged_myeloid")) {
    result <- agreement_for(arm$name, resolution)
    agreement_results[[arm$name]][[resolution]] <- result
    overall_rows[[length(overall_rows) + 1]] <- result$overall
  }
  utils::write.csv(agreement_results[[arm$name]]$full$confusion,
                   file.path(output_dir, paste0("confusion_", arm$name, ".csv")), row.names = FALSE)
  utils::write.csv(agreement_results[[arm$name]]$merged_myeloid$confusion,
                   file.path(output_dir, paste0("confusion_merged_myeloid_", arm$name, ".csv")),
                   row.names = FALSE)
  utils::write.csv(rbind(agreement_results[[arm$name]]$full$per_type,
                         agreement_results[[arm$name]]$merged_myeloid$per_type),
                   file.path(output_dir, paste0("agreement_per_type_", arm$name, ".csv")),
                   row.names = FALSE)
  message("\n--- B: ", arm$name, " annotation vs SingleR (full labels) ---")
  print(agreement_results[[arm$name]]$full$overall, row.names = FALSE)
  print(agreement_results[[arm$name]]$full$per_type, row.names = FALSE)
  print(agreement_results[[arm$name]]$full$confusion, row.names = FALSE)
}
agreement_table <- do.call(rbind, overall_rows)
utils::write.csv(agreement_table, file.path(output_dir, "annotation_vs_singler_agreement.csv"),
                 row.names = FALSE)
message("\n--- B: overall agreement ---"); print(agreement_table, row.names = FALSE)

# --- Step 7: C -- Monocyte/DC marker evidence -----------------------------------------
# The two labellings disagree most between Monocytes and DC. Four groups separate the
# cells both sides call the same from the ones they call differently; the marker table
# shows which side the expression sides with.
marker_rows <- list()
for (arm in arms) {
  if (is.null(arm$annotation_column)) next
  mine_all <- identity_by_arm[[arm$name]]
  group_defs <- list(
    "mine Monocytes / SingleR DC" = assigned & mine_all == "Monocytes" & singler$broad_label == "DC",
    "mine DC / SingleR Monocytes" = assigned & mine_all == "DC" & singler$broad_label == "Monocytes",
    "both Monocytes"              = assigned & mine_all == "Monocytes" & singler$broad_label == "Monocytes",
    "both DC"                     = assigned & mine_all == "DC" & singler$broad_label == "DC"
  )
  for (group in names(group_defs)) {
    cells <- which(group_defs[[group]])
    expressing <- Matrix::rowSums(marker_counts[, cells, drop = FALSE] > 0)
    marker_rows[[length(marker_rows) + 1]] <- data.frame(
      arm = arm$name,
      group = group,
      n_cells = length(cells),
      median_delta_next_main = median(main_run$delta.next[cells]),
      t(setNames(100 * as.numeric(expressing) / length(cells), paste0("pct_", marker_genes))),
      check.names = FALSE, stringsAsFactors = FALSE
    )
  }
}
marker_table <- do.call(rbind, marker_rows)
utils::write.csv(marker_table, file.path(output_dir, "monocyte_dc_marker_groups.csv"),
                 row.names = FALSE)
message("\n--- C: Monocyte/DC marker groups (percent expressing, counts > 0) ---")
print(marker_table, row.names = FALSE)

# --- Step 8: record what ran ------------------------------------------------------------
audit <- list(
  stage = "17_cluster_structure_vs_singler",
  packages = list(SingleR = as.character(packageVersion("SingleR")),
                  celldex = as.character(packageVersion("celldex")),
                  BiocParallel = as.character(packageVersion("BiocParallel")),
                  Seurat = seurat_version),
  r_version = paste(R.version$major, R.version$minor, sep = "."),
  arms = vapply(arms, function(a) a$name, character(1)),
  excluded_arms = list(pca_harmony_30dims =
    "30-PC variant of the PCA arm; the 15-PC pca_harmony is the main arm, as in the authors' dims = 1:15"),
  reference = "celldex::MonacoImmuneData(), labels = label.main",
  myeloid_merge = merge_table,
  n_shared_genes = length(shared_genes),
  workers = n_workers,
  main_run_wall_seconds = unname(timing["elapsed"]),
  checks = checks
)
writeLines(jsonlite::toJSON(audit, auto_unbox = TRUE, pretty = TRUE, digits = 12, na = "null"),
           file.path(output_dir, "stage17_audit_r.json"))

# --- Step 9: figures ----------------------------------------------------------------------
# One palette and one function per figure kind, called once per arm, so PCA and scGPT are
# drawn identically. Figures go to the tracked figures folder.
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

# Stage 15's palette (the authors' Figure 1D colours). Progenitors and "unassigned" have no
# Stage 15 colour and keep the ones chosen earlier in this stage. Stage 15 gives
# "Undetermined" the same gray as "unassigned" here; it is darkened so the two stay
# distinguishable when they share a legend.
CELLTYPE_COLOURS <- c(
  "T Cells"      = "darkturquoise",
  "NK Cells"     = "lightpink",
  "DC"           = "hotpink",
  "Monocytes"    = "gold1",
  "B Cells"      = "darkorange",
  "Plasma Cells" = "darkviolet",
  "Progenitors"  = "firebrick",
  "Undetermined" = "gray35",
  "unassigned"   = "gray"
)
stopifnot(all(singler$broad_label %in% names(CELLTYPE_COLOURS)))
family_colours <- c("scGPT" = "steelblue4", "PCA" = "darkorange3")

# Figure 1. Patient NMI against biology NMI, one point per arm. An arrow runs from each
# uncorrected arm to the Harmony arm derived from it.
points <- do.call(rbind, lapply(arms, function(a) {
  row <- nmi_table[nmi_table$arm == a$name, ]
  data.frame(arm = a$name, label = a$label, family = a$family,
             corrected_from = a$corrected_from,
             patient_nmi = row$patient_nmi, biology_nmi = row$biology_nmi,
             stringsAsFactors = FALSE)
}))
arrows <- points[!is.na(points$corrected_from), ]
arrows$x_start <- points$patient_nmi[match(arrows$corrected_from, points$arm)]
arrows$y_start <- points$biology_nmi[match(arrows$corrected_from, points$arm)]
x_limits <- c(0, max(0.1, 1.3 * max(points$patient_nmi)))
y_pad <- max(0.1, 0.5 * diff(range(points$biology_nmi)))
y_limits <- c(max(0, min(points$biology_nmi) - y_pad), min(1, max(points$biology_nmi) + y_pad))
# Pull each arrow back from both points so the head lands beside the corrected point.
shrink <- 0.12
arrows$x_from <- arrows$x_start + shrink * (arrows$patient_nmi - arrows$x_start)
arrows$y_from <- arrows$y_start + shrink * (arrows$biology_nmi - arrows$y_start)
arrows$x_to <- arrows$patient_nmi - shrink * (arrows$patient_nmi - arrows$x_start)
arrows$y_to <- arrows$biology_nmi - shrink * (arrows$biology_nmi - arrows$y_start)
# Uncorrected arms are labelled to the right of their point, corrected arms to the left.
points$corrected <- !is.na(points$corrected_from)
label_offset <- 0.012 * diff(x_limits) / 0.1

figure_scatter <- ggplot(points, aes(x = patient_nmi, y = biology_nmi, colour = family)) +
  geom_segment(data = arrows, aes(x = x_from, y = y_from, xend = x_to, yend = y_to),
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
               linewidth = 0.6, show.legend = FALSE) +
  geom_point(size = 3.5) +
  geom_text(aes(x = patient_nmi + ifelse(corrected, -label_offset, label_offset),
                label = label, hjust = ifelse(corrected, 1, 0)), size = 4, show.legend = FALSE) +
  scale_colour_manual(values = family_colours, guide = "none") +
  coord_cartesian(xlim = x_limits, ylim = y_limits) +
  labs(x = "Patient NMI (lower is better)", y = "Biology NMI (higher is better)") +
  theme_classic(base_size = 12)
ggsave(file.path(figure_dir, "17_patient_vs_biology_scatter.png"), figure_scatter,
       width = 6.5, height = 5, dpi = 150)

# Figure 2. Annotation against SingleR, one heatmap per Harmony arm. Rows, columns and the
# colour scale are fixed so the two arms can be read against each other; each row is a
# percentage of that annotated type, with its cell count beside the name.
heatmap_rows <- c("T Cells", "NK Cells", "DC", "Monocytes", "B Cells", "Undetermined")
heatmap_columns <- c("T Cells", "NK Cells", "DC", "Monocytes", "B Cells", "Plasma Cells", "Progenitors")
draw_heatmap <- function(arm) {
  confusion <- agreement_results[[arm$name]]$full$confusion
  if (!all(confusion$my_identity %in% heatmap_rows) ||
      !all(confusion$singler_broad_label %in% heatmap_columns)) {
    stop(arm$name, ": the confusion table has labels outside the fixed heatmap axes")
  }
  grid <- expand.grid(my_identity = heatmap_rows, singler_broad_label = heatmap_columns,
                      stringsAsFactors = FALSE)
  grid$row_pct <- confusion$row_pct[match(paste(grid$my_identity, grid$singler_broad_label),
                                          paste(confusion$my_identity, confusion$singler_broad_label))]
  grid$row_pct[is.na(grid$row_pct)] <- 0
  row_n <- tapply(confusion$n_cells, confusion$my_identity, sum)[heatmap_rows]
  if (any(abs(tapply(grid$row_pct, grid$my_identity, sum)[heatmap_rows] - 100) > 1e-9)) {
    stop(arm$name, ": a row of the heatmap does not sum to 100%")
  }
  grid$my_identity <- factor(grid$my_identity, levels = rev(heatmap_rows),
                             labels = rev(paste0(heatmap_rows, " (n = ",
                                                 prettyNum(row_n, big.mark = ","), ")")))
  grid$singler_broad_label <- factor(grid$singler_broad_label, levels = heatmap_columns)
  grid$text <- ifelse(grid$row_pct == 0, "0",
                      ifelse(grid$row_pct < 0.1, "<0.1", sprintf("%.1f", grid$row_pct)))
  figure <- ggplot(grid, aes(x = singler_broad_label, y = my_identity, fill = row_pct)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = text, colour = row_pct > 50), size = 3.8, show.legend = FALSE) +
    scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "white")) +
    scale_fill_gradient(low = "white", high = "steelblue4", limits = c(0, 100), name = "% of row") +
    scale_x_discrete(position = "top") +
    labs(x = "SingleR label", y = "Annotated cell type", title = arm$label) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(), axis.ticks = element_blank(),
          plot.title = element_text(hjust = 0))
  ggsave(file.path(figure_dir, paste0("17_agreement_heatmap_", arm$name, ".png")), figure,
         width = 10, height = 4.8, dpi = 150)
}

# Figure 3. The arm's own UMAP three times with identical axes: coloured by SingleR, by the
# annotation, and the cells where the two disagree. Cells are drawn in a seeded shuffled
# order so no type sits on top of another. Only the SingleR panel carries the legend.
draw_umap <- function(arm) {
  umap <- data.frame(
    x = as.numeric(clusters[[paste0(arm$umap_prefix, "_UMAP_1")]]),
    y = as.numeric(clusters[[paste0(arm$umap_prefix, "_UMAP_2")]]),
    singler = singler$broad_label,
    annotation = identity_by_arm[[arm$name]],
    stringsAsFactors = FALSE
  )
  if (anyNA(umap$x) || anyNA(umap$y)) stop(arm$name, ": a cell has no UMAP coordinates")
  # Undetermined cells have no annotation to disagree with, so they are not drawn as
  # disagreements; this is the same rule the agreement table uses.
  umap$disagree <- assigned & umap$annotation != "Undetermined" & umap$annotation != umap$singler
  full <- agreement_results[[arm$name]]$full$overall
  if (sum(umap$disagree) != full$n_cells_scored - full$n_agree) {
    stop(arm$name, ": the disagreement count does not match the agreement table")
  }
  set.seed(42)
  umap <- umap[sample(nrow(umap)), ]

  shown_types <- names(CELLTYPE_COLOURS)[names(CELLTYPE_COLOURS) %in% c(umap$singler, umap$annotation)]
  axes <- list(
    coord_cartesian(xlim = range(umap$x), ylim = range(umap$y)),
    labs(x = "UMAP 1", y = "UMAP 2"),
    theme_classic(base_size = 12)
  )
  # A type absent from the panel that carries the legend would get a blank key. An
  # invisible point per type (alpha 0 in the plot, full in the legend) keeps every key.
  legend_keys <- data.frame(x = umap$x[1], y = umap$y[1], type = shown_types)
  colour_panel <- function(column, title, legend) {
    ggplot(umap, aes(x = x, y = y, colour = .data[[column]])) +
      geom_point(size = 0.3, stroke = 0) +
      geom_point(data = legend_keys, aes(colour = type), alpha = 0, size = 0.3, stroke = 0) +
      scale_colour_manual(values = CELLTYPE_COLOURS, limits = shown_types, name = NULL,
                          guide = if (legend) guide_legend(override.aes = list(size = 5, alpha = 1)) else "none") +
      ggtitle(title) + axes +
      theme(legend.key.size = unit(0.7, "cm"), legend.text = element_text(size = 12))
  }
  disagree_panel <- ggplot(umap, aes(x = x, y = y)) +
    geom_point(colour = "gray90", size = 0.3, stroke = 0) +
    geom_point(data = umap[umap$disagree, ], colour = "blue3", size = 0.3, stroke = 0) +
    ggtitle(sprintf("Disagree (n = %s)", prettyNum(sum(umap$disagree), big.mark = ","))) + axes
  figure <- (colour_panel("singler", paste(arm$label, "- SingleR"), TRUE) |
               colour_panel("annotation", paste(arm$label, "- Annotation"), FALSE) |
               disagree_panel) +
    plot_layout(guides = "collect")
  ggsave(file.path(figure_dir, paste0("17_umap_singler_vs_annotation_", arm$name, ".png")),
         figure, width = 17, height = 5.5, dpi = 150)
}

for (arm in arms) {
  if (is.null(arm$annotation_column)) next
  draw_heatmap(arm)
  draw_umap(arm)
}

message("\nStage 17 (R) complete.")
