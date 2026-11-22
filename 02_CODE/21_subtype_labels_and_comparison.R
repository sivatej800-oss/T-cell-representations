# Stage 21: name the subclusters, draw Figure 1D for both arms, and compare the arms.
#
# The authors named their T-cell and monocyte subclusters by hand (their lines 309-357)
# and wrote the names back onto the full dataset, which is Figure 1D: twelve labels on
# the full UMAP. This stage does both for each arm, from Stage 20's subclusters, and
# then compares the two arms' subtypes without any reference labels.
#
# The names below were entered by a person reading Stage 20's evidence tables and
# feature plots, as the authors read theirs, and approved before this stage ran. The
# categories are exactly the authors' Figure 1D ones. Their "Monoctyes" typo in the
# CD68hi label is corrected; nothing else is changed.
#
# The comparison has three parts, each within a lineage (T cells, monocytes).
#   1. Cross-tabulation of the two arms' subtype labels on the cells both arms placed
#      in the lineage: for every cell, whether the arms agree, and where they don't.
#   2. Marker purity from raw counts: the share of each label's cells expressing the
#      authors' panel genes, and the same for each disagreement pair, which says which
#      arm the expression sides with where they disagree.
#   3. Patient NMI of each arm's subclustering, on its own lineage cells: how much of
#      the subclustering is explained by which patient a cell came from.
#
# Usage: Rscript 02_CODE/21_subtype_labels_and_comparison.R <workspace>

suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: 21_subtype_labels_and_comparison.R <workspace>")
results_dir <- file.path(normalizePath(args[1]), "results")
subcluster_dir <- file.path(results_dir, "20_subclusters")
output_dir <- file.path(results_dir, "21_subtypes")
figure_dir <- file.path(results_dir, "figures")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

expected_cells <- 70391
arms <- c("authors_pipeline", "scgpt_raw")
arm_labels <- c(authors_pipeline = "Authors' pipeline", scgpt_raw = "scGPT raw")
lineages <- c(t_cells = "T Cells", monocytes = "Monocytes")

CD4 <- "CD4+ T Cells"; CD8 <- "CD8+ T Cells"; DP <- "CD4+/CD8+ T Cells"; TREG <- "T Regulatory Cells"
LO <- "CD14+/CD16-/CD68lo Monocytes"; MID <- "CD14+/CD16+/CD68mid Monocytes"
HI <- "CD14+/CD16+/CD68hi Monocytes"; UND <- "Undetermined"


# --- The manual subtype annotation ---------------------------------------------------
# Read off Stage 20's evidence. CD8A/CD8B are the sharper CD4-vs-CD8 split: CD4 mRNA is
# captured in only about 40% of cells even in clear CD4 subclusters.
subtype_annotation <- list(
  # Authors' pipeline T cells. 4 is the only FOXP3-enriched subcluster (17%, every
  # other at 2% or less) and carries naive/central-memory genes beside TIGIT and ICOS,
  # so it is Tregs mixed with naive CD4 cells. 3 is cytotoxic and HLA-DR+ with CD8A 37%
  # and CD4 16%: the authors' mixed CD4+/CD8+ category.
  authors_pipeline_t_cells = list("0" = CD4, "1" = CD4, "2" = CD8, "3" = DP, "4" = TREG, "5" = CD4),
  # scGPT T cells. No subcluster is FOXP3-enriched (at most 5%), so no Treg label is
  # given. 2 carries an immediate-early signature (CD69, JUN, NFKBIA) and 4 is naive
  # CD4 with 42% of its cells from one patient; both are CD4 on the panel. 3 (LTB,
  # ribosomal; CD4 14%, CD8B 17%) and 5 (cycling) are mixed.
  scgpt_raw_t_cells = list("0" = CD8, "1" = CD4, "2" = CD4, "3" = DP, "4" = CD4, "5" = DP),
  # Monocytes, both arms: CD16 (FCGR3A) presence plus CD68 level. The LYVE1/TIMD4 and
  # the C1Q/APOE/TREM2 macrophage-like subclusters are CD68hi, the second on its
  # macrophage programme as much as its CD68 level.
  authors_pipeline_monocytes = list("0" = LO, "1" = HI, "2" = MID, "3" = MID, "4" = HI),
  # scGPT monocyte 4 expresses CD3D in 85% of cells alongside T-cell and monocyte
  # genes, with 53% from one patient: T-cell/monocyte doublets or complexes.
  scgpt_raw_monocytes = list("0" = HI, "1" = LO, "2" = MID, "3" = HI, "4" = UND, "5" = MID)
)

# The authors' Figure 1D order and colours (their lines 381-386 and 412-416).
FIGURE1D_LEVELS <- c(CD4, CD8, DP, TREG, "NK Cells", "DC", LO, MID, HI,
                     "B Cells", "Plasma Cells", UND)
FIGURE1D_COLOURS <- setNames(c("darkturquoise", "lawngreen", "dodgerblue", "red",
                               "lightpink", "hotpink", "gold1", "navy", "violet",
                               "darkorange", "darkviolet", "gray"), FIGURE1D_LEVELS)

# Copied from 02_CODE/06_benchmark_representations.R so this stage runs on its own.
normalized_mutual_information <- function(labels_a, labels_b) {
  joint <- table(labels_a, labels_b) / length(labels_a)
  margin_a <- rowSums(joint); margin_b <- colSums(joint)
  entropy <- function(p) { p <- p[p > 0]; -sum(p * log(p)) }
  outer_product <- outer(margin_a, margin_b)
  positive <- joint > 0
  mutual <- sum(joint[positive] * log(joint[positive] / outer_product[positive]))
  denominator <- sqrt(entropy(margin_a) * entropy(margin_b))
  if (denominator == 0) return(NA_real_)
  mutual / denominator
}


# --- Step 1: apply the names, checked ------------------------------------------------
annotation <- utils::read.csv(file.path(results_dir, "15_annotation", "annotated_cluster_labels.csv.gz"),
                              colClasses = "character")
if (nrow(annotation) != expected_cells) stop("Stage 15's file should hold ", expected_cells, " cells")
barcodes <- annotation$cell_barcode

subclusters <- list()
for (arm in arms) for (tag in names(lineages)) {
  key <- paste0(arm, "_", tag)
  table <- utils::read.csv(file.path(subcluster_dir, paste0("subclusters_", key, ".csv.gz")),
                           colClasses = c(cell_barcode = "character", subcluster = "character"))
  expected <- barcodes[annotation[[paste0(arm, "_cluster_identity")]] == lineages[[tag]]]
  if (!setequal(table$cell_barcode, expected)) stop(key, ": Stage 20 and Stage 15 disagree on the cells")
  assignment <- subtype_annotation[[key]]
  if (!setequal(names(assignment), unique(table$subcluster))) {
    stop(key, ": every subcluster needs exactly one name, and only existing ones may be named")
  }
  if (!all(unlist(assignment) %in% FIGURE1D_LEVELS)) stop(key, ": a name is not a Figure 1D category")
  table$subtype <- unname(unlist(assignment)[table$subcluster])
  subclusters[[key]] <- table
}

# Figure 1D labels per cell: the subtype where the cell was reclustered, otherwise its
# Stage 15 broad identity.
figure1d <- data.frame(cell_barcode = barcodes, stringsAsFactors = FALSE)
for (arm in arms) {
  labels <- annotation[[paste0(arm, "_cluster_identity")]]
  for (tag in names(lineages)) {
    table <- subclusters[[paste0(arm, "_", tag)]]
    labels[match(table$cell_barcode, barcodes)] <- table$subtype
  }
  if (!all(labels %in% FIGURE1D_LEVELS)) stop(arm, ": a cell carries a label outside Figure 1D")
  figure1d[[arm]] <- labels
}
utils::write.csv(figure1d, gzfile(file.path(output_dir, "figure1d_labels.csv.gz")), row.names = FALSE)
counts <- do.call(rbind, lapply(arms, function(arm) {
  data.frame(arm = arm, label = FIGURE1D_LEVELS,
             n_cells = as.integer(table(factor(figure1d[[arm]], levels = FIGURE1D_LEVELS))))
}))
utils::write.csv(counts, file.path(output_dir, "figure1d_label_counts.csv"), row.names = FALSE)
message("Figure 1D label counts:"); print(reshape(counts, idvar = "label", timevar = "arm", direction = "wide"),
                                         row.names = FALSE)


# --- Step 2: raw counts for the marker panel -----------------------------------------
message("\nReading the RNA matrix for marker purity...")
input_dir <- file.path(results_dir, "03_seurat_reproduction", "seurat_input")
raw_all <- as(Matrix::readMM(file.path(input_dir, "counts.mtx")), "CsparseMatrix")
genes <- readLines(file.path(input_dir, "genes.tsv"), warn = FALSE)
matrix_barcodes <- readLines(file.path(input_dir, "barcodes.tsv"), warn = FALSE)
if (ncol(raw_all) != expected_cells || !setequal(matrix_barcodes, barcodes)) stop("RNA matrix cells differ")
rownames(raw_all) <- genes; colnames(raw_all) <- matrix_barcodes
panel <- c("CD3D", "CD4", "CD8A", "CD8B", "FOXP3", "CD14", "FCGR3A", "CD68")
library_size <- Matrix::colSums(raw_all)
raw <- raw_all[panel, barcodes, drop = FALSE]
library_size <- library_size[barcodes]
rm(raw_all); invisible(gc())
# LogNormalize, as Stage 13: log1p of counts per 10,000.
lognorm <- log1p(t(t(as.matrix(raw)) / library_size) * 1e4)
metadata <- utils::read.csv(file.path(input_dir, "metadata.csv"), colClasses = "character")
patients <- metadata$internal_sample_id[match(barcodes, metadata$cell_barcode)]

marker_profile <- function(cells) {
  index <- match(cells, barcodes)
  c(n_cells = length(index),
    setNames(round(100 * Matrix::rowSums(raw[, index, drop = FALSE] > 0) / length(index), 1),
             paste0("pct_", panel)),
    setNames(round(rowMeans(lognorm[, index, drop = FALSE]), 3), paste0("mean_", panel)))
}


# --- Step 3: compare the arms within each lineage ------------------------------------
purity_rows <- list(); pair_rows <- list(); agreement_rows <- list(); nmi_rows <- list()
crosstab_plots <- list()
for (tag in names(lineages)) {
  a <- subclusters[[paste0("authors_pipeline_", tag)]]
  s <- subclusters[[paste0("scgpt_raw_", tag)]]

  # Marker purity per arm and label, on that arm's own lineage cells.
  for (arm in arms) {
    table <- subclusters[[paste0(arm, "_", tag)]]
    for (label in intersect(FIGURE1D_LEVELS, unique(table$subtype))) {
      purity_rows[[length(purity_rows) + 1]] <- data.frame(
        lineage = lineages[[tag]], arm = arm, subtype = label,
        t(marker_profile(table$cell_barcode[table$subtype == label])), check.names = FALSE)
    }
    nmi_rows[[length(nmi_rows) + 1]] <- data.frame(
      lineage = lineages[[tag]], arm = arm, n_cells = nrow(table),
      n_subclusters = length(unique(table$subcluster)),
      patient_nmi_subclusters = normalized_mutual_information(table$subcluster, table$patient),
      patient_nmi_subtypes = normalized_mutual_information(table$subtype, table$patient))
  }

  # Cells both arms placed in this lineage, and the cross-tabulation of their labels.
  shared <- intersect(a$cell_barcode, s$cell_barcode)
  label_a <- a$subtype[match(shared, a$cell_barcode)]
  label_s <- s$subtype[match(shared, s$cell_barcode)]
  levels_here <- FIGURE1D_LEVELS[FIGURE1D_LEVELS %in% c(label_a, label_s)]
  crosstab <- table(authors = factor(label_a, levels = levels_here),
                    scgpt = factor(label_s, levels = levels_here))
  utils::write.csv(as.data.frame.matrix(crosstab),
                   file.path(output_dir, paste0("crosstab_", tag, ".csv")))
  agreement_rows[[length(agreement_rows) + 1]] <- data.frame(
    lineage = lineages[[tag]],
    authors_cells = nrow(a), scgpt_cells = nrow(s), shared_cells = length(shared),
    only_authors = nrow(a) - length(shared), only_scgpt = nrow(s) - length(shared),
    same_label = sum(label_a == label_s),
    pct_same_label = round(100 * mean(label_a == label_s), 2))

  # The marker profile of every (authors' label, scGPT label) pair with cells in it.
  pairs <- unique(data.frame(authors = label_a, scgpt = label_s, stringsAsFactors = FALSE))
  for (i in seq_len(nrow(pairs))) {
    cells <- shared[label_a == pairs$authors[i] & label_s == pairs$scgpt[i]]
    pair_rows[[length(pair_rows) + 1]] <- data.frame(
      lineage = lineages[[tag]], authors_label = pairs$authors[i], scgpt_label = pairs$scgpt[i],
      agree = pairs$authors[i] == pairs$scgpt[i], t(marker_profile(cells)), check.names = FALSE)
  }

  # Heatmap: rows are the authors' labels, each cell the share of that row's cells
  # carrying each scGPT label; counts beside the names.
  grid <- as.data.frame(as.table(100 * sweep(crosstab, 1, pmax(rowSums(crosstab), 1), "/")),
                        stringsAsFactors = FALSE)
  colnames(grid) <- c("authors", "scgpt", "pct")
  row_n <- rowSums(crosstab); col_n <- colSums(crosstab)
  grid$authors <- factor(grid$authors, levels = rev(levels_here),
                         labels = rev(sprintf("%s (n = %s)", levels_here, prettyNum(row_n, big.mark = ","))))
  grid$scgpt <- factor(grid$scgpt, levels = levels_here,
                       labels = sprintf("%s\n(n = %s)", levels_here, prettyNum(col_n, big.mark = ",")))
  grid$text <- ifelse(grid$pct < 0.05, "0", sprintf("%.1f", grid$pct))
  crosstab_plots[[tag]] <- ggplot(grid, aes(x = scgpt, y = authors, fill = pct)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = text, colour = pct > 50), size = 3.6, show.legend = FALSE) +
    scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "white")) +
    scale_fill_gradient(low = "white", high = "steelblue4", limits = c(0, 100), name = "% of row") +
    scale_x_discrete(position = "top") +
    labs(x = "scGPT raw", y = "Authors' pipeline",
         title = sprintf("%s: same label for %.1f%% of %s shared cells", lineages[[tag]],
                         100 * mean(label_a == label_s), prettyNum(length(shared), big.mark = ","))) +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(), axis.ticks = element_blank(),
          axis.text.x.top = element_text(size = 9), plot.title = element_text(hjust = 0))
  ggsave(file.path(figure_dir, paste0("21_crosstab_", tag, ".png")), crosstab_plots[[tag]],
         width = 11, height = 4.5 + 0.35 * length(levels_here), dpi = 150)
}
purity <- do.call(rbind, purity_rows); pairs_table <- do.call(rbind, pair_rows)
agreement <- do.call(rbind, agreement_rows); nmi <- do.call(rbind, nmi_rows)
utils::write.csv(purity, file.path(output_dir, "marker_purity_by_subtype.csv"), row.names = FALSE)
utils::write.csv(pairs_table, file.path(output_dir, "marker_purity_by_label_pair.csv"), row.names = FALSE)
utils::write.csv(agreement, file.path(output_dir, "subtype_agreement.csv"), row.names = FALSE)
utils::write.csv(nmi, file.path(output_dir, "subcluster_patient_nmi.csv"), row.names = FALSE)
message("\nAgreement:"); print(agreement, row.names = FALSE)
message("\nPatient NMI:"); print(nmi, row.names = FALSE)


# --- Step 4: the UMAPs ---------------------------------------------------------------
# Cells are drawn in a seeded shuffled order so no label sits on top of another. Both
# panels of a figure share one colour key (`shown`, the labels either arm uses), drawn
# once, from the first panel. A label that panel lacks would get a blank key, so one
# invisible point per label (alpha 0 in the plot, full in the legend) keeps every key
# coloured. Points grow as the cell count falls, so 7,000 monocytes are as legible as
# 70,000 cells.
umap_panel <- function(frame, title, shown, show_legend) {
  set.seed(42)
  frame <- frame[sample(nrow(frame)), ]
  frame$label <- factor(frame$label, levels = shown)
  keys <- data.frame(x = frame$x[1], y = frame$y[1], label = factor(shown, levels = shown))
  point_size <- min(1.2, max(0.3, 0.35 * sqrt(50000 / nrow(frame))))
  ggplot(frame, aes(x = x, y = y, colour = label)) +
    geom_point(size = point_size, stroke = 0) +
    geom_point(data = keys, alpha = 0, size = point_size, stroke = 0) +
    scale_colour_manual(values = FIGURE1D_COLOURS[shown], limits = shown, name = NULL,
                        guide = if (show_legend) guide_legend(override.aes = list(size = 4, alpha = 1))
                                else "none") +
    labs(x = "UMAP 1", y = "UMAP 2", title = title) +
    theme_classic(base_size = 12) +
    theme(legend.text = element_text(size = 10))
}

# Each lineage on its own reclustered UMAP, the two arms side by side.
shared_key <- function(labels) FIGURE1D_LEVELS[FIGURE1D_LEVELS %in% labels]
for (tag in names(lineages)) {
  shown <- shared_key(unlist(lapply(arms, function(arm) subclusters[[paste0(arm, "_", tag)]]$subtype)))
  panels <- lapply(arms, function(arm) {
    table <- subclusters[[paste0(arm, "_", tag)]]
    umap_panel(data.frame(x = table$UMAP_1, y = table$UMAP_2, label = table$subtype),
               sprintf("%s — %s (%s cells)", arm_labels[[arm]], lineages[[tag]],
                       prettyNum(nrow(table), big.mark = ",")), shown, arm == arms[1])
  })
  ggsave(file.path(figure_dir, paste0("21_subtypes_", tag, ".png")),
         wrap_plots(panels, ncol = 2) + plot_layout(guides = "collect"),
         width = 15, height = 6, dpi = 150)
}

# Figure 1D: the full dataset with all twelve labels, each arm on its own full UMAP
# (the authors' on Stage 03's, scGPT's on Stage 04's of the same embedding).
layouts <- list(
  authors_pipeline = file.path(results_dir, "03_seurat_reproduction", "GSE200164_umap_clusters.csv.gz"),
  scgpt_raw = file.path(results_dir, "04_scgpt_embedding", "GSE200164_scgpt_umap_clusters.csv.gz"))
shown_1d <- shared_key(unlist(figure1d[arms]))
figure1d_panels <- lapply(arms, function(arm) {
  layout <- utils::read.csv(layouts[[arm]], colClasses = c(cell_barcode = "character"))
  if (!setequal(layout$cell_barcode, barcodes)) stop(arm, ": layout cells differ")
  layout <- layout[match(barcodes, layout$cell_barcode), ]
  umap_panel(data.frame(x = layout$UMAP_1, y = layout$UMAP_2, label = figure1d[[arm]]),
             paste0(arm_labels[[arm]], " — Figure 1D labels"), shown_1d, arm == arms[1])
})
ggsave(file.path(figure_dir, "21_figure1d.png"),
       wrap_plots(figure1d_panels, ncol = 2) + plot_layout(guides = "collect"),
       width = 17, height = 7, dpi = 150)

writeLines(toJSON(list(
  stage = "21_subtype_labels_and_comparison",
  follows = "05_REFERENCES/author_code/figure1d_author_workflow.R lines 309-386",
  annotation = subtype_annotation, panel = panel,
  expression = "raw counts > 0 for presence; log1p(counts / library size * 1e4) for means",
  agreement = agreement, patient_nmi = nmi
), auto_unbox = TRUE, pretty = TRUE, digits = 10), file.path(output_dir, "subtype_audit.json"))
message("\nStage 21 complete: ", output_dir)
