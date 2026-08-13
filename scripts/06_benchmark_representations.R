#!/usr/bin/env Rscript

# A like-for-like comparison of the two cell descriptions, with the batch confound
# removed from both.
#
# Stages 03 and 04 cannot be compared directly. Stage 03 aligns the 59 samples with
# integration anchors before reducing dimensions; Stage 04 embeds each cell with no
# alignment at all, because scGPT expects count-scale input and correcting expression
# beforehand would put its input outside the range the model was trained on. So one
# arm is batch-corrected and the other is not, and any difference between their
# clusterings could be donor structure surviving in the second rather than anything
# about the descriptions themselves.
#
# The fix is to correct after the description exists rather than before it. Harmony
# works on a finished low-dimensional representation: it repeatedly finds groups of
# cells that ought to be one cell type, measures how far each donor sits from that
# group's centre, and subtracts the offset. It needs no assumption that the
# dimensions are ordered by variance or independently meaningful, which is what makes
# it valid on a learned embedding where neither holds.
#
# That gives a 2x2. Both descriptions appear corrected and uncorrected, so the effect
# of the correction is visible separately from the effect of the description:
#
#                    uncorrected          Harmony-corrected
#   PCA 50 -> 15      pca_raw              pca_harmony
#   scGPT 512         scgpt_raw            scgpt_harmony
#
# A fifth arm clusters the corrected PCA on 30 dimensions instead of 15, to show
# whether the conclusion depends on that cut. Everything downstream of the
# representation is identical across arms: same neighbours, same resolution 0.3,
# same seed.
#
# Correcting for donor equally is the confound this script removes. Three others
# remain, upstream of anything here, and they qualify any conclusion drawn from it.
#
# The first is not a difference at all, and is recorded because it looks like one.
# SCTransform normalizes each cell for sequencing depth and stabilizes variance
# before the PCA arm reduces dimensions; the scGPT arm appears to skip that. It does
# not need it. scGPT sorts each cell's values into 51 bins cut at that same cell's
# own quantiles, which is a ranking within the cell, so scaling a cell by its library
# size or taking a logarithm leaves every bin unchanged. The normalization the PCA
# arm performs explicitly is absorbed by the binning rather than omitted.
#
# The second is a real asymmetry that cannot be removed without retraining. The
# per-sample normalization also regressed out mitochondrial fraction, so that
# technical signal is gone from the PCA arm. Binning removes depth but not
# mitochondrial content: a stressed cell keeps its mitochondrial genes ranked high
# and the bins preserve it. The scGPT arm therefore still carries a technical
# gradient the PCA arm does not.
#
# The third is the largest and is a limit on precision rather than a bias. scGPT
# reads at most 1,200 genes per cell and samples them at random when a cell has more.
# The median cell here has 1,691 genes detected and 83% exceed the limit, so for five
# cells in six roughly a third of the measured genes were discarded before the model
# saw them, chosen by an unseeded permutation. Re-embedding would produce a different
# result. How different is not known: it would take embedding the data two or three
# times and comparing, which is the measurement that would say how large a difference
# between the arms has to be before it means anything. Until that exists, small
# differences in the numbers below should not be interpreted.
#
# What this script does NOT do is decide which description is better. There are no
# cell-type labels for these cells — the scGPT manifest carries a cell_type column
# and every value in it is "unknown" — so there is no external truth to score
# against. What is measurable is how well each arm mixes donors, whether rare
# populations survive as small clusters, and how far the two descriptions agree with
# each other. Deciding which is right requires annotating the clusters by marker
# expression, which is a separate piece of work.
#
# One consequence of the donor correction must be understood before any of these
# outputs are reused. Each of the 59 samples is one patient, and age and diagnosis
# are properties of the patient, constant across all of that patient's cells. Donor
# therefore determines age exactly. Removing between-donor variation removes all
# age and all diagnosis variation with it, necessarily and by construction, not as a
# side effect that could be tuned away. That is correct for grouping cells by type,
# which is all this script does. It makes the corrected representations unusable for
# the age trajectory question: that analysis must return to the uncorrected
# expression and test within clusters, never on these coordinates.
#
# This script is started by 06_benchmark_representations.py and is not run directly.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop(
    "Usage: 06_benchmark_representations.R <scgpt_input_dir> <sct_dir> ",
    "<output_dir> <seed> <n_cells> <n_dims>"
  )
}

scgpt_input_dir <- normalizePath(args[[1]], mustWork = TRUE)
sct_dir <- normalizePath(args[[2]], mustWork = TRUE)
output_dir <- normalizePath(args[[3]], mustWork = FALSE)
random_seed <- as.integer(args[[4]])
expected_cells <- as.integer(args[[5]])
expected_dims <- as.integer(args[[6]])

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(harmony)
  library(ggplot2)
})
set.seed(random_seed)

# Settings shared by every arm, so that no arm differs from another by anything
# except the coordinates it is handed.
CLUSTER_RESOLUTION <- 0.3
UMAP_NEIGHBORS <- 30
MIXING_NEIGHBORS <- 30


# --- Agreement and mixing measures -------------------------------------------
# Two clusterings of the same cells are compared by how often they place a pair of
# cells the same way — both together, or both apart. Raw agreement is misleading
# because two unrelated clusterings already agree on most pairs by chance, so the
# chance level is subtracted and the result scaled: 0 means no better than chance,
# 1 means identical groupings.
adjusted_rand_index <- function(labels_a, labels_b) {
  contingency <- table(labels_a, labels_b)
  pairs <- function(x) sum(x * (x - 1) / 2)
  index <- pairs(as.vector(contingency))
  a_pairs <- pairs(rowSums(contingency))
  b_pairs <- pairs(colSums(contingency))
  total_pairs <- pairs(sum(contingency))
  expected <- a_pairs * b_pairs / total_pairs
  maximum <- (a_pairs + b_pairs) / 2
  if (isTRUE(all.equal(maximum, expected))) return(NA_real_)
  (index - expected) / (maximum - expected)
}

# A second view of the same question, in terms of information rather than pairs: how
# much knowing a cell's cluster in one arm tells you about its cluster in the other,
# scaled so that 1 is complete determination. It is reported alongside the pair
# measure because the two disagree when cluster sizes are very uneven, and these
# clusterings contain populations differing by more than an order of magnitude.
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

# How well the donors are mixed in a representation, judged locally. For each cell,
# look at the cells nearest it and count how many came from a different patient. If
# the representation still separates patients, a cell's neighbours are mostly its own
# patient's cells and the count is low.
#
# The raw fraction cannot be read on its own, because how high it can go depends on
# how the cells are distributed across patients: with 59 unequal samples, a randomly
# chosen neighbour comes from a different patient about 97% of the time. Dividing by
# that ceiling puts the result on a fixed scale where 1 is as mixed as the sample
# sizes allow and 0 is complete separation by patient.
donor_mixing_score <- function(object, reduction, dims, donors) {
  coordinates <- Embeddings(object, reduction = reduction)[, dims, drop = FALSE]
  neighbors <- FindNeighbors(
    coordinates,
    k.param = MIXING_NEIGHBORS + 1,
    return.neighbor = TRUE,
    verbose = FALSE
  )
  indices <- Indices(neighbors)
  # The first neighbour of a cell is itself, and a cell always shares its own donor,
  # so including it would drag every score down by a fixed amount.
  indices <- indices[, -1, drop = FALSE]
  neighbor_donors <- matrix(
    donors[as.vector(indices)],
    nrow = nrow(indices)
  )
  observed <- mean(neighbor_donors != donors)
  proportions <- as.vector(table(donors)) / length(donors)
  ceiling_value <- 1 - sum(proportions^2)
  c(observed = observed, ceiling = ceiling_value, score = observed / ceiling_value)
}


# --- Step 1: rebuild the PCA description without anchor integration ----------
# Stage 03's principal components cannot be reused here. They were computed on
# expression that anchor integration had already corrected, so they carry a batch
# correction of a different kind and would make the PCA arm doubly corrected while
# the scGPT arm is corrected once. This recomputes them from the same per-sample
# normalizations Stage 03 starts from, merged but left uncorrected, which is the
# input Harmony expects.
#
# The gene set is chosen exactly as Stage 03 chooses it, so that the two arms differ
# in how batch is handled and in nothing else. Requesting 3,000 returns the whole
# available pool of about 2,816; the reason to ask for more than is available is that
# any cut ranks out genes that vary in only a few samples, which is precisely what
# markers of rare cell types look like.
pca_checkpoint <- file.path(output_dir, "benchmark_pca_uncorrected.rds")

if (file.exists(pca_checkpoint)) {
  message("Resuming from uncorrected PCA checkpoint: ", pca_checkpoint)
  reference <- readRDS(pca_checkpoint)
} else {
  message("Reading per-sample SCTransform checkpoints...")
  sct_files <- list.files(sct_dir, pattern = "_sct\\.rds$", full.names = TRUE)
  sample_names <- sub("_sct\\.rds$", "", basename(sct_files))
  if (length(sct_files) != 59) {
    stop("Expected 59 per-sample SCTransform checkpoints, found ", length(sct_files))
  }
  sct_list <- lapply(sct_files, readRDS)
  names(sct_list) <- sample_names
  if (sum(vapply(sct_list, ncol, numeric(1))) != expected_cells) {
    stop("SCTransform checkpoints do not contain all ", expected_cells, " cells")
  }

  message("Selecting features...")
  features <- SelectIntegrationFeatures(object.list = sct_list, nfeatures = 3000)
  # A gene is only usable if every sample measured it; one sample here holds 151
  # cells and contributes few genes, so this intersection is the real constraint.
  for (sample_name in sample_names) {
    features <- intersect(
      features,
      rownames(GetAssayData(sct_list[[sample_name]], assay = "SCT", layer = "data"))
    )
  }
  message("Features shared by all samples: ", length(features))
  if (length(features) < 1000) {
    stop("Too few shared features to describe cell types: ", length(features))
  }

  # Assemble the normalized values into one matrix. The samples are stacked as they
  # are, with no adjustment between them: leaving the donor differences in place is
  # what gives Harmony something to remove, and what makes the uncorrected arm a
  # genuine control rather than a partially corrected one.
  message("Assembling merged expression matrix...")
  expression <- do.call(cbind, lapply(sample_names, function(sample_name) {
    GetAssayData(sct_list[[sample_name]], assay = "SCT", layer = "data")[features, ]
  }))
  cell_donors <- rep(sample_names, vapply(sct_list, ncol, numeric(1))[sample_names])
  rm(sct_list)
  gc()
  if (ncol(expression) != expected_cells) {
    stop("Merged matrix holds ", ncol(expression), " cells, expected ", expected_cells)
  }

  reference <- CreateSeuratObject(
    counts = expression,
    assay = "SCTmerged",
    meta.data = data.frame(
      internal_sample_id = cell_donors,
      row.names = colnames(expression)
    )
  )
  reference <- SetAssayData(reference, assay = "SCTmerged", layer = "data", new.data = expression)
  VariableFeatures(reference) <- features
  rm(expression)
  gc()

  # Centre and rescale each gene so that highly expressed genes do not dominate the
  # decomposition purely through their magnitude, then condense to 50 dimensions.
  # Fifty is the field default and more than any arm uses; computing beyond what is
  # needed costs little and keeps the decomposition's approximation error out of the
  # comparison, since how many are requested slightly perturbs the leading ones.
  message("Scaling and running PCA...")
  set.seed(random_seed)
  reference <- ScaleData(reference, features = features, verbose = FALSE)
  reference <- RunPCA(reference, features = features, npcs = 50, verbose = FALSE)
  saveRDS(reference, pca_checkpoint, compress = FALSE)
}

donors <- reference$internal_sample_id
if (length(unique(donors)) != 59) {
  stop("Expected 59 donors, found ", length(unique(donors)))
}


# --- Step 2: attach the scGPT description to the same cells ------------------
# The embedding arrives as a flat block of 32-bit numbers in row order, the same
# artifact Stage 04 reads. It describes the same cells but was written in its own
# order, so it is reordered to match this object rather than assumed to align. An
# embedding whose rows were silently misassigned would still cluster and still
# produce a plausible figure, and every conclusion drawn from it would be wrong.
message("Reading scGPT embedding...")
embedding_path <- file.path(scgpt_input_dir, "embedding_float32.bin")
metadata_path <- file.path(scgpt_input_dir, "metadata.csv")
connection <- file(embedding_path, "rb")
values <- readBin(
  connection,
  what = "numeric",
  n = expected_cells * expected_dims,
  size = 4,
  endian = "little"
)
close(connection)
if (length(values) != expected_cells * expected_dims) {
  stop("Embedding binary is not the expected size")
}
embedding <- matrix(values, nrow = expected_cells, ncol = expected_dims, byrow = TRUE)
rm(values)

scgpt_metadata <- read.csv(metadata_path, stringsAsFactors = FALSE)
if (!"cell_barcode" %in% names(scgpt_metadata)) {
  stop("scGPT metadata is missing the cell_barcode column")
}
rownames(embedding) <- scgpt_metadata$cell_barcode
colnames(embedding) <- paste0("scgpt_", seq_len(expected_dims))
if (!setequal(rownames(embedding), colnames(reference))) {
  stop("The embedding does not describe the same cells as the PCA arm")
}
embedding <- embedding[colnames(reference), , drop = FALSE]
if (!identical(rownames(embedding), colnames(reference))) {
  stop("Failed to align the embedding to the cell order")
}
reference[["scgpt"]] <- CreateDimReducObject(
  embeddings = embedding,
  key = "scgpt_",
  assay = "SCTmerged"
)
rm(embedding)
gc()


# --- Step 3: correct each description for donor ------------------------------
# Harmony is given the patient each cell came from and nothing else. It is run
# separately on each description, so neither inherits the other's correction.
message("Running Harmony on the PCA description...")
set.seed(random_seed)
reference <- RunHarmony(
  reference,
  group.by.vars = "internal_sample_id",
  reduction.use = "pca",
  reduction.save = "pca_harmony",
  dims.use = 1:50,
  verbose = TRUE
)

message("Running Harmony on the scGPT description...")
set.seed(random_seed)
reference <- RunHarmony(
  reference,
  group.by.vars = "internal_sample_id",
  reduction.use = "scgpt",
  reduction.save = "scgpt_harmony",
  dims.use = 1:expected_dims,
  verbose = TRUE
)
# Correction shifts cells within the space; it does not add or remove any. The scGPT
# vectors arrive unit length and do not stay that way, which changes the geometry
# neighbours are measured in from angular to ordinary distance.
for (name in c("pca_harmony", "scgpt_harmony")) {
  if (nrow(Embeddings(reference, reduction = name)) != expected_cells) {
    stop("Harmony returned the wrong number of cells for ", name)
  }
  if (anyNA(Embeddings(reference, reduction = name))) {
    stop("Harmony produced missing values for ", name)
  }
}


# --- Step 4: cluster every arm the same way ----------------------------------
# Each arm supplies coordinates; everything after that is fixed. The PCA arms use the
# leading 15 dimensions, the number the authors used and the number Stage 03 uses.
# The scGPT arms use all 512, because a learned embedding has no ordering that would
# make a leading subset meaningful. The fifth arm repeats the corrected PCA on 30
# dimensions, so that a difference between arms can be checked against how much the
# result moves when only that cut changes.
arms <- list(
  list(name = "pca_raw", reduction = "pca", dims = 1:15),
  list(name = "pca_harmony", reduction = "pca_harmony", dims = 1:15),
  list(name = "pca_harmony_30dims", reduction = "pca_harmony", dims = 1:30),
  list(name = "scgpt_raw", reduction = "scgpt", dims = 1:expected_dims),
  list(name = "scgpt_harmony", reduction = "scgpt_harmony", dims = 1:expected_dims)
)

cluster_table <- data.frame(cell_barcode = colnames(reference), stringsAsFactors = FALSE)
summaries <- list()

for (arm in arms) {
  message("Clustering arm: ", arm$name)
  set.seed(random_seed)
  reference <- FindNeighbors(
    reference,
    reduction = arm$reduction,
    dims = arm$dims,
    verbose = FALSE
  )
  reference <- FindClusters(
    reference,
    resolution = CLUSTER_RESOLUTION,
    algorithm = 1,
    random.seed = random_seed,
    verbose = FALSE
  )
  assignments <- Idents(reference)
  cluster_table[[arm$name]] <- as.character(assignments)

  # How well donors mix in this arm's own coordinates, and whether small groups
  # survived. A smallest cluster in the thousands means rare cell types were
  # absorbed into larger ones rather than resolved, which is the failure mode
  # Stage 03 documents from its own earlier run.
  mixing <- donor_mixing_score(reference, arm$reduction, arm$dims, donors)
  sizes <- as.vector(table(assignments))
  summaries[[arm$name]] <- data.frame(
    arm = arm$name,
    reduction = arm$reduction,
    dims_used = length(arm$dims),
    n_clusters = length(sizes),
    smallest_cluster = min(sizes),
    largest_cluster = max(sizes),
    donor_mixing_observed = unname(mixing["observed"]),
    donor_mixing_score = unname(mixing["score"]),
    stringsAsFactors = FALSE
  )
}

summary_table <- do.call(rbind, summaries)
rownames(summary_table) <- NULL
message("\n--- Arm summary ---")
print(summary_table, row.names = FALSE)


# --- Step 5: measure how far the arms agree ----------------------------------
# With no cell-type labels to score against, the remaining question the data can
# answer is whether the two descriptions find the same structure. High agreement
# between the corrected arms would mean the embedding recovers what the expression
# already showed; low agreement means they disagree, and which one is right is not
# decidable here.
arm_names <- vapply(arms, function(a) a$name, character(1))
pairs <- t(combn(arm_names, 2))
agreement_table <- data.frame(
  arm_a = pairs[, 1],
  arm_b = pairs[, 2],
  adjusted_rand = apply(pairs, 1, function(p) {
    adjusted_rand_index(cluster_table[[p[1]]], cluster_table[[p[2]]])
  }),
  normalized_mutual_info = apply(pairs, 1, function(p) {
    normalized_mutual_information(cluster_table[[p[1]]], cluster_table[[p[2]]])
  }),
  stringsAsFactors = FALSE
)
message("\n--- Agreement between arms ---")
print(agreement_table, row.names = FALSE)

# How much of each arm's clustering is explained by which patient a cell came from.
# An uncorrected arm that splits largely along patients will score high here, and
# that is the confound the correction exists to remove.
donor_leakage <- data.frame(
  arm = arm_names,
  donor_nmi = vapply(
    arm_names,
    function(name) normalized_mutual_information(cluster_table[[name]], donors),
    numeric(1)
  ),
  stringsAsFactors = FALSE
)
message("\n--- Cluster agreement with donor identity ---")
print(donor_leakage, row.names = FALSE)


# --- Step 6: lay out the corrected arms for viewing --------------------------
# Only the two corrected arms are drawn. The uncorrected ones exist as controls for
# the numbers above; drawing them invites reading donor separation as biology.
for (arm in arms[c(2, 5)]) {
  message("Building UMAP for ", arm$name)
  set.seed(random_seed)
  reference <- RunUMAP(
    reference,
    reduction = arm$reduction,
    dims = arm$dims,
    n.neighbors = UMAP_NEIGHBORS,
    reduction.name = paste0("umap_", arm$name),
    seed.use = random_seed,
    verbose = FALSE
  )
  coordinates <- Embeddings(reference, reduction = paste0("umap_", arm$name))
  if (nrow(coordinates) != expected_cells || anyNA(coordinates)) {
    stop("UMAP failed validation for ", arm$name)
  }
  cluster_table[[paste0(arm$name, "_UMAP_1")]] <- coordinates[, 1]
  cluster_table[[paste0(arm$name, "_UMAP_2")]] <- coordinates[, 2]
}


# --- Step 7: save ------------------------------------------------------------
cluster_table$internal_sample_id <- donors
write.csv(
  cluster_table,
  gzfile(file.path(output_dir, "benchmark_cluster_assignments.csv.gz")),
  row.names = FALSE
)
write.csv(
  summary_table,
  file.path(output_dir, "benchmark_arm_summary.csv"),
  row.names = FALSE
)
write.csv(
  agreement_table,
  file.path(output_dir, "benchmark_arm_agreement.csv"),
  row.names = FALSE
)
write.csv(
  donor_leakage,
  file.path(output_dir, "benchmark_donor_leakage.csv"),
  row.names = FALSE
)
saveRDS(reference, file.path(output_dir, "benchmark_representations.rds"), compress = FALSE)

# The two corrected arms side by side, coloured by cluster. The panels are separate
# layouts of separate coordinates, so a cluster's position in one has no relationship
# to a position in the other; only the grouping structure is comparable.
for (arm_name in c("pca_harmony", "scgpt_harmony")) {
  plot_data <- data.frame(
    UMAP_1 = cluster_table[[paste0(arm_name, "_UMAP_1")]],
    UMAP_2 = cluster_table[[paste0(arm_name, "_UMAP_2")]],
    cluster = cluster_table[[arm_name]]
  )
  figure <- ggplot(plot_data, aes(UMAP_1, UMAP_2, colour = cluster)) +
    geom_point(size = 0.1, alpha = 0.5) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    theme_classic() +
    ggtitle(paste0(arm_name, ", Louvain resolution ", CLUSTER_RESOLUTION))
  ggsave(
    file.path(output_dir, paste0("benchmark_", arm_name, "_umap.png")),
    figure,
    width = 8,
    height = 6,
    dpi = 300
  )
}

message("\nStage 06 complete. Outputs written to ", output_dir)
