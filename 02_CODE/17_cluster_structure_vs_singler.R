# Stage 17: does cluster structure follow biology or the patient? (scGPT arms)
#
# Stage 06 measured how strongly each clustering tracks which patient a cell came from
# (donor NMI). That is only half of the question: a clustering can be mostly biology and
# still carry some patient signal, or the reverse. This stage measures the other half --
# how strongly each clustering tracks cell type -- using labels that belong to neither
# arm: SingleR's per-cell calls against purified immune profiles (Stage 16).
#
# Three things are measured for the scGPT arms.
#   A. Patient NMI and biology NMI side by side, per arm.
#   B. How closely the hand annotation of scGPT-Harmony clusters (Stage 15's table)
#      agrees with SingleR, cell by cell.
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
if (length(args) != 8) {
  stop("Usage: 17_cluster_structure_vs_singler.R <seurat_input_dir> <cluster_csv> ",
       "<donor_leakage_csv> <singler_csv> <old_stage15_csv> <output_dir> <workers> <figure_dir>")
}
input_dir <- args[1]
cluster_csv <- args[2]
leakage_csv <- args[3]
singler_csv <- args[4]
old_stage15_csv <- args[5]
output_dir <- args[6]
n_workers <- as.integer(args[7])
figure_dir <- args[8]

expected_cells <- 70391
expected_genes <- 19689

# The RNA data-layer sum Stage 13 recorded (and Stage 16 reproduced). Stage 17 reads the
# same unchanged seurat_input/, so it must reproduce this number.
stage13_data_sum <- 176503092.94553709
sum_tolerance <- 1e-9

# --- The arms ---------------------------------------------------------------------
# One entry per clustering. `annotation` is the hand annotation (cluster -> identity)
# where one exists; an arm without one gets the NMI measurements only. To add the PCA
# arm later, add an entry here and rerun.
#
# scGPT-Harmony annotation: copied from 02_CODE/15_annotate_clusters.R. It is the author's
# reading of Stage 14's marker evidence, applied to the rebuilt Stage 06 clusters. Louvain
# numbers clusters by size, so the numbers carry no meaning outside this clustering.
arms <- list(
  list(
    name = "scgpt_raw",
    cluster_column = "scgpt_raw",
    annotation = NULL
  ),
  list(
    name = "scgpt_harmony",
    cluster_column = "scgpt_harmony",
    annotation = list(
      "T Cells"   = c("0", "1", "4", "7", "10"),
      "Monocytes" = c("2", "6"),
      "DC"        = c("3", "8"),
      "NK Cells"  = c("5"),
      "B Cells"   = c("9")
    )
  )
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
stopifnot(identical(clusters$cell_barcode, barcodes), identical(singler$cell_barcode, barcodes),
          identical(old15$cell_barcode, barcodes))
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

# --- Step 4: A -- patient NMI and biology NMI per arm --------------------------------
# The patient NMI is recomputed here on all cells and must equal the value Stage 06 wrote;
# agreement shows this stage reads the same clusters and the same patients as Stage 06.
leakage <- utils::read.csv(leakage_csv, stringsAsFactors = FALSE)
expected_patient_nmi <- c(scgpt_raw = 0.0721503826775315, scgpt_harmony = 0.0506991838687104)

nmi_rows <- list()
for (arm in arms) {
  cluster_ids <- clusters[[arm$cluster_column]]
  if (anyNA(cluster_ids)) stop(arm$name, ": a cell has no cluster")
  patient_nmi <- normalized_mutual_information(cluster_ids, patients)
  recorded <- leakage$donor_nmi[leakage$arm == arm$name]
  if (length(recorded) != 1 || abs(patient_nmi - recorded) > 1e-12) {
    stop(arm$name, ": patient NMI ", patient_nmi, " differs from Stage 06's ", recorded)
  }
  if (arm$name %in% names(expected_patient_nmi) &&
      abs(patient_nmi - expected_patient_nmi[[arm$name]]) > 1e-12) {
    stop(arm$name, ": patient NMI differs from the expected value")
  }
  biology_nmi <- normalized_mutual_information(cluster_ids[assigned], singler$broad_label[assigned])
  nmi_rows[[arm$name]] <- data.frame(
    arm = arm$name,
    n_clusters = length(unique(cluster_ids)),
    patient_nmi = patient_nmi,
    biology_nmi = biology_nmi,
    n_cells_scored = sum(assigned),
    stringsAsFactors = FALSE
  )
}
nmi_table <- do.call(rbind, nmi_rows)
utils::write.csv(nmi_table, file.path(output_dir, "nmi_by_arm.csv"), row.names = FALSE)
message("\n--- A: patient vs biology NMI ---"); print(nmi_table, row.names = FALSE)
checks$patient_nmi_matches_stage06 <- TRUE

# --- Step 5: build each annotated arm's per-cell identity ----------------------------
# Every cluster must receive exactly one identity. A typed cluster number that silently
# missed would leave cells without an identity and still produce a clean-looking table.
identity_by_arm <- list()
annotation_rows <- list()
for (arm in arms) {
  if (is.null(arm$annotation)) next
  cluster_ids <- clusters[[arm$cluster_column]]
  all_clusters <- sort(unique(cluster_ids))
  listed <- unlist(arm$annotation, use.names = FALSE)
  if (anyDuplicated(listed)) stop(arm$name, ": a cluster is listed under two identities")
  if (!setequal(listed, all_clusters)) {
    stop(arm$name, ": annotation and clusters differ. Unannotated: ",
         paste(setdiff(all_clusters, listed), collapse = ","), "; nonexistent: ",
         paste(setdiff(listed, all_clusters), collapse = ","))
  }
  lookup <- setNames(rep(names(arm$annotation), lengths(arm$annotation)), listed)
  identity_by_arm[[arm$name]] <- unname(lookup[cluster_ids])
  if (anyNA(identity_by_arm[[arm$name]])) stop(arm$name, ": a cell has no identity")
  annotation_rows[[arm$name]] <- data.frame(arm = arm$name, cluster = names(lookup),
                                            identity = unname(lookup), stringsAsFactors = FALSE)
}
annotation_table <- do.call(rbind, annotation_rows)
annotation_table <- annotation_table[order(annotation_table$arm,
                                           as.integer(annotation_table$cluster)), ]
utils::write.csv(annotation_table, file.path(output_dir, "annotation_table_used.csv"),
                 row.names = FALSE)

# The rebuild's annotation should reproduce the old Stage 15 identities for most cells.
old_match <- list()
for (arm in arms) {
  if (is.null(arm$annotation)) next
  old_identity <- old15[[paste0(arm$cluster_column, "_cluster_identity")]]
  if (is.null(old_identity)) next
  old_match[[arm$name]] <- list(pct_cells_same_identity_as_old_stage15 =
                                  100 * mean(identity_by_arm[[arm$name]] == old_identity))
  message(sprintf("%s: %.2f%% of cells keep their old Stage 15 identity", arm$name,
                  old_match[[arm$name]]$pct_cells_same_identity_as_old_stage15))
}
checks$identity_vs_old_stage15 <- old_match

# --- Step 6: B -- hand annotation vs SingleR ------------------------------------------
# Agreement is an exact label match on assigned cells. SingleR's Plasma Cells and
# Progenitors have no counterpart in the hand annotation, so those cells can only
# disagree; they are counted separately so they are not mistaken for ordinary mislabels.
agreement_rows <- list()
for (arm in arms) {
  if (is.null(arm$annotation)) next
  mine <- identity_by_arm[[arm$name]][assigned]
  theirs <- singler$broad_label[assigned]
  same <- mine == theirs

  confusion <- as.data.frame(table(my_identity = mine, singler_broad_label = theirs),
                             stringsAsFactors = FALSE)
  confusion <- confusion[confusion$Freq > 0, ]
  names(confusion)[3] <- "n_cells"
  row_totals <- table(mine)
  confusion$row_pct <- 100 * confusion$n_cells / as.integer(row_totals[confusion$my_identity])
  confusion$arm <- arm$name
  confusion <- confusion[order(confusion$my_identity, -confusion$n_cells), c("arm", "my_identity",
                                                         "singler_broad_label", "n_cells", "row_pct")]
  utils::write.csv(confusion, file.path(output_dir, paste0("confusion_", arm$name, ".csv")),
                   row.names = FALSE)

  per_type <- do.call(rbind, lapply(names(row_totals), function(type) {
    data.frame(arm = arm$name, my_identity = type, n_cells = as.integer(row_totals[type]),
               n_agree = sum(same & mine == type),
               pct_agree = 100 * sum(same & mine == type) / as.integer(row_totals[type]),
               stringsAsFactors = FALSE)
  }))
  utils::write.csv(per_type, file.path(output_dir, paste0("agreement_per_type_", arm$name, ".csv")),
                   row.names = FALSE)

  no_counterpart <- table(theirs[!theirs %in% names(arm$annotation)])
  agreement_rows[[arm$name]] <- data.frame(
    arm = arm$name,
    n_cells_scored = length(mine),
    n_agree = sum(same),
    pct_agree = 100 * mean(same),
    singler_plasma_cells = as.integer(sum(theirs == "Plasma Cells")),
    singler_progenitors = as.integer(sum(theirs == "Progenitors")),
    singler_other_without_counterpart = as.integer(sum(no_counterpart)) -
      sum(theirs == "Plasma Cells") - sum(theirs == "Progenitors"),
    stringsAsFactors = FALSE
  )
  message("\n--- B: ", arm$name, " annotation vs SingleR ---")
  print(agreement_rows[[arm$name]], row.names = FALSE)
  print(per_type, row.names = FALSE)
  print(confusion, row.names = FALSE)
}
agreement_table <- do.call(rbind, agreement_rows)
utils::write.csv(agreement_table, file.path(output_dir, "annotation_vs_singler_agreement.csv"),
                 row.names = FALSE)

# --- Step 7: C -- Monocyte/DC marker evidence -----------------------------------------
# The two labellings disagree most between Monocytes and DC. Four groups separate the
# cells both sides call the same from the ones they call differently; the marker table
# shows which side the expression sides with.
marker_rows <- list()
for (arm in arms) {
  if (is.null(arm$annotation)) next
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
  reference = "celldex::MonacoImmuneData(), labels = label.main",
  n_shared_genes = length(shared_genes),
  workers = n_workers,
  main_run_wall_seconds = unname(timing["elapsed"]),
  checks = checks
)
writeLines(jsonlite::toJSON(audit, auto_unbox = TRUE, pretty = TRUE, digits = 12, na = "null"),
           file.path(output_dir, "stage17_audit_r.json"))

# --- Step 9: figures ----------------------------------------------------------------------
# Three plain figures for the scGPT results, written to the tracked figures folder. Cell
# types keep one colour across all of them.
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

# Stage 15's palette (the authors' Figure 1D colours), copied here so this stage stands
# alone. "unassigned" and "Progenitors" have no Stage 15 colour; they get gray and
# firebrick, which collide with nothing else here.
CELLTYPE_COLOURS <- c(
  "T Cells"      = "darkturquoise",
  "NK Cells"     = "lightpink",
  "DC"           = "hotpink",
  "Monocytes"    = "gold1",
  "B Cells"      = "darkorange",
  "Plasma Cells" = "darkviolet",
  "Progenitors"  = "firebrick",
  "unassigned"   = "gray"
)
stopifnot(all(singler$broad_label %in% names(CELLTYPE_COLOURS)))

# Figure 1. Patient vs biology NMI, one group per arm in the order of the arm list, so a
# new arm adds its group with no other change.
nmi_long <- rbind(
  data.frame(arm = nmi_table$arm, measure = "Patient", nmi = nmi_table$patient_nmi),
  data.frame(arm = nmi_table$arm, measure = "Cell type (SingleR)", nmi = nmi_table$biology_nmi)
)
nmi_long$arm <- factor(nmi_long$arm, levels = vapply(arms, function(a) a$name, character(1)))
nmi_long$measure <- factor(nmi_long$measure, levels = c("Patient", "Cell type (SingleR)"))
figure_nmi <- ggplot(nmi_long, aes(x = arm, y = nmi, fill = measure)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", nmi)), position = position_dodge(width = 0.8),
            vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = c("Patient" = "gray55", "Cell type (SingleR)" = "steelblue4")) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(x = NULL, y = "Normalized mutual information with clusters", fill = NULL) +
  theme_classic(base_size = 12)
ggsave(file.path(figure_dir, "17_nmi_patient_vs_biology.png"), figure_nmi,
       width = 6, height = 5, dpi = 150)

# Figures 2 and 3 are for the scGPT-Harmony arm, the one with a hand annotation.
figure_arm <- "scgpt_harmony"
mine_identity <- identity_by_arm[[figure_arm]]

# Figure 2. Agreement with SingleR per annotated cell type, read back from the table the
# run wrote, with the number of cells behind each bar.
per_type_table <- utils::read.csv(
  file.path(output_dir, paste0("agreement_per_type_", figure_arm, ".csv")),
  stringsAsFactors = FALSE)
per_type_table$my_identity <- factor(
  per_type_table$my_identity,
  levels = intersect(names(CELLTYPE_COLOURS), per_type_table$my_identity))
figure_agreement <- ggplot(per_type_table, aes(x = my_identity, y = pct_agree, fill = my_identity)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = format(n_cells, big.mark = ",")), vjust = 1.5, size = 3.5) +
  scale_fill_manual(values = CELLTYPE_COLOURS, guide = "none") +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  labs(x = "Annotated cell type (cells on bar)", y = "Agreement with SingleR (%)") +
  theme_classic(base_size = 12)
ggsave(file.path(figure_dir, "17_agreement_per_celltype_scgpt_harmony.png"), figure_agreement,
       width = 6.5, height = 5, dpi = 150)

# Figure 3. The same UMAP twice, coloured by SingleR and by the hand annotation. Cells are
# drawn commonest type first so rare types are not buried. The legend is shared, so a
# type that appears in only one panel still has a colour; only the SingleR panel, which has
# every type, carries the legend.
umap <- data.frame(
  x = as.numeric(clusters[[paste0(figure_arm, "_UMAP_1")]]),
  y = as.numeric(clusters[[paste0(figure_arm, "_UMAP_2")]]),
  singler = singler$broad_label,
  annotation = mine_identity,
  stringsAsFactors = FALSE
)
if (anyNA(umap$x) || anyNA(umap$y)) stop("A cell has no UMAP coordinates")
shown_types <- names(CELLTYPE_COLOURS)[names(CELLTYPE_COLOURS) %in% c(umap$singler, umap$annotation)]
umap_panel <- function(column, label, legend = TRUE) {
  frequency <- table(umap[[column]])
  drawn <- umap[order(-as.integer(frequency[umap[[column]]])), ]
  ggplot(drawn, aes(x = x, y = y, colour = .data[[column]])) +
    geom_point(size = 0.15, stroke = 0) +
    scale_colour_manual(values = CELLTYPE_COLOURS, limits = shown_types, name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    coord_cartesian(xlim = range(umap$x), ylim = range(umap$y)) +
    labs(x = "UMAP 1", y = "UMAP 2", title = label) +
    theme_classic(base_size = 12) +
    theme(legend.position = if (legend) "right" else "none")
}
# The SingleR panel holds every type, so its legend serves both panels.
figure_umap <- umap_panel("singler", "SingleR") | umap_panel("annotation", "Annotation", legend = FALSE)
ggsave(file.path(figure_dir, "17_umap_singler_vs_annotation_scgpt_harmony.png"), figure_umap,
       width = 12, height = 5.5, dpi = 150)

message("\nStage 17 (R) complete.")
