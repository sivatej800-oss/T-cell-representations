# Stage 16: an independent, reference-based cell-type label for every cell.
#
# Stages 13-15 named clusters from marker genes read by eye. That judgement is only as
# good as the reader, and it is tied to the partitions it was applied to. This stage
# asks a different question that involves neither partition: which known immune cell
# type does each individual cell's expression most resemble?
#
# SingleR answers it by correlating each cell against purified bulk-RNA profiles of
# sorted immune populations (the Monaco reference, 29 fine types in 10 main types). It
# is run per cell, with no clusters, so the label for one cell never depends on which
# cluster a representation put it in. That independence is the point: the labels can
# later grade both the PCA arm and the scGPT arm without belonging to either.
#
# Cells whose best label does not clearly beat the alternatives are "unassigned"
# (SingleR's pruned label is NA). They stay unassigned in every label column. Forcing a
# type onto an ambiguous cell would manufacture agreement.
#
# The expression is the same RNA Stage 13 used: SoupX-corrected counts from the GEO
# hand-off, fractional on purpose, log-normalized once. Nothing is rounded or filtered.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(SingleR)
  library(celldex)
  library(BiocParallel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: 16_singler_reference_labels.R <seurat_input_dir> <output_dir> <stage13_audit_json>")
}
input_dir <- args[1]
output_dir <- args[2]
stage13_audit_path <- args[3]

expected_cells <- 70391
expected_genes <- 19689

# The RNA data-layer sum Stage 13 recorded in marker_discovery_audit.json. Stage 16 reads
# the same unchanged seurat_input/ and normalizes it the same way, so it must reproduce
# this number. It is fixed here because that audit file may not exist when this runs.
stage13_data_sum <- 176503092.94553709
sum_tolerance <- 1e-9

# One core stays free for other work on the machine.
n_workers <- parallel::detectCores() - 1

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(42)

peak_rss_gb <- function() {
  # VmHWM is the parent process's peak resident memory; forked workers are not included.
  status <- readLines("/proc/self/status")
  kb <- as.numeric(gsub("[^0-9]", "", grep("^VmHWM:", status, value = TRUE)))
  kb / 1024^2
}

# --- Step 1: rebuild Stage 13's expression exactly ---------------------------------
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

so <- CreateSeuratObject(counts = counts, assay = "RNA", project = "GSE200164_singler",
                         min.cells = 0, min.features = 0)
rm(counts)
gc(verbose = FALSE)
so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 10000,
                    verbose = FALSE)

# Confirm this is the same expression Stage 13 judged, not merely the same shape. A
# changed matrix would give plausible labels for the wrong data.
lognorm <- GetAssayData(so, assay = "RNA", layer = "data")
observed_sum <- sum(lognorm)
check_sum <- function(expected, source) {
  relative_error <- abs(observed_sum - expected) / expected
  if (relative_error > sum_tolerance) {
    stop(sprintf("Data-layer sum %.8f differs from %s (%.8f); relative error %.3g",
                 observed_sum, source, expected, relative_error))
  }
  relative_error
}
relative_error <- check_sum(stage13_data_sum, "Stage 13's recorded value")
stage13_file_sum <- NA_real_
if (file.exists(stage13_audit_path)) {
  stage13_file_sum <- jsonlite::fromJSON(stage13_audit_path)$expression_fingerprint$data_sum
  check_sum(stage13_file_sum, "marker_discovery_audit.json")
  message("Also matches the live Stage 13 audit file.")
} else {
  message("Stage 13 audit file not present; relying on the recorded value.")
}
message(sprintf("Expression: %d genes x %d cells, data sum %.6f (relative error %.2e)",
                nrow(lognorm), ncol(lognorm), observed_sum, relative_error))

seurat_version <- as.character(packageVersion("Seurat"))
rm(so)
gc(verbose = FALSE)

# --- Step 2: load the reference and verify the label vocabulary --------------------
message("Loading celldex::MonacoImmuneData()...")
ref <- celldex::MonacoImmuneData()
message("label.main values:\n  ", paste(sort(unique(ref$label.main)), collapse = "\n  "))
message("label.fine values:\n  ", paste(sort(unique(ref$label.fine)), collapse = "\n  "))

# The main label is taken from the reference's own fine-to-main correspondence, so the
# two levels cannot disagree. That only works if each fine type belongs to one main type.
pairs <- unique(data.frame(label_fine = ref$label.fine, label_main = ref$label.main,
                           stringsAsFactors = FALSE))
if (anyDuplicated(pairs$label_fine)) {
  stop("A fine label maps to more than one main label in the reference")
}
fine_to_main <- setNames(pairs$label_main, pairs$label_fine)

# --- Step 3: the broad-label mapping, written out as a table -----------------------
# The broad vocabulary matches the paper's coarse populations. Rules apply in order:
#   1. Plasmablasts become Plasma Cells before any main-label rule, because Monaco files
#      them under "B cells" but they are a distinct population here.
#   2-6. Monaco's main labels fold into the paper's names. Plasmacytoid DCs stay in DC.
#   7. Anything else (Neutrophils, Basophils, Progenitors) keeps its Monaco name and is
#      never merged into a neighbouring type.
#   8. Unassigned cells stay "unassigned".
fine_rule <- data.frame(label_fine = "Plasmablasts", broad_label = "Plasma Cells",
                        rule = 1, stringsAsFactors = FALSE)
main_rules <- data.frame(
  label_main = c("CD4+ T cells", "CD8+ T cells", "T cells", "Monocytes",
                 "Dendritic cells", "NK cells", "B cells",
                 "Neutrophils", "Basophils", "Progenitors"),
  broad_label = c("T Cells", "T Cells", "T Cells", "Monocytes",
                  "DC", "NK Cells", "B Cells",
                  "Neutrophils", "Basophils", "Progenitors"),
  rule = c(2, 2, 2, 3, 4, 5, 6, 7, 7, 7),
  stringsAsFactors = FALSE
)

# Every name in the table must exist in the reference, and every reference main label
# must be covered. A rule that silently matches nothing would leave cells unmapped.
if (!all(fine_rule$label_fine %in% ref$label.fine)) {
  stop("A fine label in the mapping table is not in the reference")
}
if (!all(main_rules$label_main %in% ref$label.main)) {
  stop("A main label in the mapping table is not in the reference: ",
       paste(setdiff(main_rules$label_main, ref$label.main), collapse = ", "))
}
if (!setequal(main_rules$label_main, ref$label.main)) {
  stop("The reference has main labels the mapping table does not cover: ",
       paste(setdiff(ref$label.main, main_rules$label_main), collapse = ", "))
}

mapping <- data.frame(label_fine = names(fine_to_main), label_main = unname(fine_to_main),
                      stringsAsFactors = FALSE)
mapping$broad_label <- main_rules$broad_label[match(mapping$label_main, main_rules$label_main)]
mapping$rule <- main_rules$rule[match(mapping$label_main, main_rules$label_main)]
is_plasmablast <- mapping$label_fine %in% fine_rule$label_fine
mapping$broad_label[is_plasmablast] <- fine_rule$broad_label
mapping$rule[is_plasmablast] <- fine_rule$rule
if (anyNA(mapping$broad_label)) stop("A fine label was left without a broad label")
mapping <- rbind(mapping, data.frame(label_fine = "unassigned", label_main = "unassigned",
                                     broad_label = "unassigned", rule = 8))
utils::write.csv(mapping, file.path(output_dir, "mapping_table.csv"), row.names = FALSE)

# --- Step 4: SingleR, once, per cell ------------------------------------------------
shared_genes <- intersect(rownames(lognorm), rownames(ref))
message(sprintf("Genes shared with the reference: %d of %d", length(shared_genes),
                nrow(lognorm)))
test_matrix <- lognorm[shared_genes, , drop = FALSE]
ref_matrix <- ref[shared_genes, ]
rm(lognorm)
gc(verbose = FALSE)

message(sprintf("Running SingleR on %d cells with %d workers...", ncol(test_matrix), n_workers))
gc(reset = TRUE, verbose = FALSE)
set.seed(42)
timing <- system.time(
  singler <- SingleR(test = test_matrix, ref = ref_matrix, labels = ref_matrix$label.fine,
                     BPPARAM = BiocParallel::MulticoreParam(workers = n_workers))
)
wall_seconds <- unname(timing["elapsed"])
gc_after <- gc(verbose = FALSE)
peak_gb <- peak_rss_gb()
message(sprintf("SingleR wall time %.1f s; parent peak memory %.2f GB", wall_seconds, peak_gb))

if (!identical(rownames(singler), barcodes)) {
  stop("SingleR output is not in the matrix's cell order; barcodes cannot be trusted")
}

# --- Step 5: assemble the per-cell labels --------------------------------------------
# Pruned labels are the assignments; the unpruned best label is kept only for later
# inspection of what the unassigned cells looked like.
pruned <- singler$pruned.labels
unassigned <- is.na(pruned)
label_fine <- ifelse(unassigned, "unassigned", pruned)
label_main <- ifelse(unassigned, "unassigned", unname(fine_to_main[pruned]))
broad_label <- mapping$broad_label[match(label_fine, mapping$label_fine)]
if (anyNA(label_main) || anyNA(broad_label)) stop("A cell was left without a label")

cells <- data.frame(
  cell_barcode = rownames(singler),
  label_fine = label_fine,
  label_main = label_main,
  broad_label = broad_label,
  is_unassigned = unassigned,
  unpruned_label_fine = singler$labels,
  delta_next = singler$delta.next,
  best_score = apply(singler$scores, 1, max),
  stringsAsFactors = FALSE
)

if (nrow(cells) != expected_cells || anyDuplicated(cells$cell_barcode) ||
    !setequal(cells$cell_barcode, barcodes)) {
  stop("The label table does not cover exactly the cells in seurat_input")
}
if (any(cells$is_unassigned & cells$broad_label != "unassigned")) {
  stop("An unassigned cell received a type")
}

utils::write.csv(cells, gzfile(file.path(output_dir, "singler_cell_labels.csv.gz")),
                 row.names = FALSE)

count_table <- function(values) {
  counts <- sort(table(values), decreasing = TRUE)
  data.frame(label = names(counts), n_cells = as.integer(counts),
             pct = 100 * as.integer(counts) / length(values), stringsAsFactors = FALSE)
}
broad_counts <- count_table(cells$broad_label)
fine_counts <- count_table(cells$label_fine)
utils::write.csv(broad_counts, file.path(output_dir, "label_counts_broad.csv"), row.names = FALSE)
utils::write.csv(fine_counts, file.path(output_dir, "label_counts_fine.csv"), row.names = FALSE)
message("\nBroad labels:"); print(broad_counts, row.names = FALSE)
message("\nFine labels:"); print(fine_counts, row.names = FALSE)

# --- Step 6: record what ran -----------------------------------------------------------
audit <- list(
  stage = "16_singler_reference_labels",
  packages = list(SingleR = as.character(packageVersion("SingleR")),
                  celldex = as.character(packageVersion("celldex")),
                  BiocParallel = as.character(packageVersion("BiocParallel")),
                  Seurat = seurat_version),
  r_version = paste(R.version$major, R.version$minor, sep = "."),
  reference = "celldex::MonacoImmuneData()",
  reference_labels = list(n_main = length(unique(ref$label.main)),
                          n_fine = length(unique(ref$label.fine))),
  n_shared_genes = length(shared_genes),
  fingerprint = list(
    expected_data_sum = stage13_data_sum,
    expected_source = "Stage 13 marker_discovery_audit.json, recorded value",
    observed_data_sum = observed_sum,
    relative_error = relative_error,
    relative_tolerance = sum_tolerance,
    checked_against_live_audit_file = file.exists(stage13_audit_path)
  ),
  n_cells = nrow(cells),
  n_unassigned = sum(cells$is_unassigned),
  pct_unassigned = 100 * mean(cells$is_unassigned),
  workers = n_workers,
  singler_wall_seconds = wall_seconds,
  peak_memory_gb_parent_process = peak_gb,
  gc_max_used_mb_after_singler = as.numeric(sum(gc_after[, 6]))
)
writeLines(jsonlite::toJSON(audit, auto_unbox = TRUE, pretty = TRUE, digits = 12),
           file.path(output_dir, "singler_audit.json"))
message("\nStage 16 complete.")
