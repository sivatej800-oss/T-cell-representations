# Stage 17 figures: two plain pictures drawn from Stage 17's finished tables.
#
# Nothing is recomputed here. The numbers come straight from nmi_by_arm.csv and
# confusion_<arm>.csv, written by 17_cluster_structure_vs_singler.py, so the figures can be
# redrawn in seconds without rerunning SingleR.
#
#   1. Patient vs biology: one point per arm. A clustering that tracks cell type and not
#      the patient sits top-left. An arrow runs from each uncorrected arm to its
#      Harmony-corrected version, showing what the correction did.
#   2. Agreement with SingleR: for each hand-annotated arm, what SingleR calls the cells of
#      each annotated type, as a percentage of that type.
#
# Usage: Rscript 17_make_figures.R <results/17_biology_vs_donor directory>
# Figures are written to figures/ inside that directory.

suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: 17_make_figures.R <17_biology_vs_donor directory>")
results_dir <- args[1]
figure_dir <- file.path(results_dir, "figures")
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

# --- The arms -------------------------------------------------------------------------
# One entry per clustering, in plotting order. `corrected_from` names the uncorrected arm a
# Harmony arm was derived from (NA for an uncorrected arm); it draws the arrow. `annotated`
# says whether a hand annotation exists, which is what the heatmap needs. `family` groups
# arms by representation so each gets one colour. To add the PCA arm, add entries here.
arms <- list(
  list(name = "scgpt_raw",     label = "scGPT raw",       family = "scGPT",
       corrected_from = NA_character_, annotated = FALSE),
  list(name = "scgpt_harmony", label = "scGPT + Harmony", family = "scGPT",
       corrected_from = "scgpt_raw",   annotated = TRUE)
)
family_colours <- c("scGPT" = "steelblue4", "PCA" = "darkorange3")
arm_names <- vapply(arms, function(a) a$name, character(1))

# The cell-type colours used elsewhere in Stage 17 are not needed here: the heatmap is one
# colour scale on purpose, so the eye reads percentages rather than categories.

# --- Figure 1: patient NMI vs biology NMI ----------------------------------------------
nmi <- utils::read.csv(file.path(results_dir, "nmi_by_arm.csv"), stringsAsFactors = FALSE)
if (!all(arm_names %in% nmi$arm)) {
  stop("nmi_by_arm.csv lacks arms: ", paste(setdiff(arm_names, nmi$arm), collapse = ", "))
}
points <- do.call(rbind, lapply(arms, function(a) {
  row <- nmi[nmi$arm == a$name, ]
  data.frame(arm = a$name, label = a$label, family = a$family,
             corrected_from = a$corrected_from,
             patient_nmi = row$patient_nmi, biology_nmi = row$biology_nmi,
             stringsAsFactors = FALSE)
}))

arrows <- points[!is.na(points$corrected_from), ]
arrows$x_start <- points$patient_nmi[match(arrows$corrected_from, points$arm)]
arrows$y_start <- points$biology_nmi[match(arrows$corrected_from, points$arm)]

# Axis limits follow the data, with room around it, so an added arm is not clipped.
x_limits <- c(0, max(0.1, 1.3 * max(points$patient_nmi)))
y_pad <- max(0.1, 0.5 * diff(range(points$biology_nmi)))
y_limits <- c(max(0, min(points$biology_nmi) - y_pad), min(1, max(points$biology_nmi) + y_pad))

# Pull each arrow back from both points, by a fraction of its length measured in axis
# units, so the head lands beside the corrected point instead of on top of it.
shrink <- 0.12
arrows$x_from <- arrows$x_start + shrink * (arrows$patient_nmi - arrows$x_start)
arrows$y_from <- arrows$y_start + shrink * (arrows$biology_nmi - arrows$y_start)
arrows$x_to <- arrows$patient_nmi - shrink * (arrows$patient_nmi - arrows$x_start)
arrows$y_to <- arrows$biology_nmi - shrink * (arrows$biology_nmi - arrows$y_start)

# Uncorrected arms are labelled to the right of their point, corrected arms to the left,
# which keeps both labels clear of the arrow between them.
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
       width = 5.5, height = 5, dpi = 150)

# --- Figure 2: annotation vs SingleR heatmap -------------------------------------------
# Rows are the hand-annotated types, columns SingleR's labels, each cell the percentage of
# its row, so every row sums to 100. The row's cell count is printed beside the type name,
# since a percentage of a few cells is not the same evidence as a percentage of thousands.
label_order <- c("T Cells", "NK Cells", "DC", "Monocytes", "B Cells", "Plasma Cells", "Progenitors")
for (a in arms) {
  if (!a$annotated) next
  confusion <- utils::read.csv(file.path(results_dir, paste0("confusion_", a$name, ".csv")),
                               stringsAsFactors = FALSE)
  rows <- label_order[label_order %in% confusion$my_identity]
  columns <- label_order[label_order %in% confusion$singler_broad_label]
  grid <- expand.grid(my_identity = rows, singler_broad_label = columns,
                      stringsAsFactors = FALSE)
  grid$row_pct <- confusion$row_pct[match(paste(grid$my_identity, grid$singler_broad_label),
                                          paste(confusion$my_identity, confusion$singler_broad_label))]
  grid$row_pct[is.na(grid$row_pct)] <- 0

  row_n <- tapply(confusion$n_cells, confusion$my_identity, sum)[rows]
  if (any(abs(tapply(grid$row_pct, grid$my_identity, sum)[rows] - 100) > 1e-9)) {
    stop(a$name, ": a row of the heatmap does not sum to 100%")
  }
  row_labels <- paste0(rows, " (n = ", prettyNum(row_n, big.mark = ","), ")")
  grid$my_identity <- factor(grid$my_identity, levels = rev(rows), labels = rev(row_labels))
  grid$singler_broad_label <- factor(grid$singler_broad_label, levels = columns)
  grid$text <- ifelse(grid$row_pct == 0, "0",
                      ifelse(grid$row_pct < 0.1, "<0.1", sprintf("%.1f", grid$row_pct)))

  figure_heatmap <- ggplot(grid, aes(x = singler_broad_label, y = my_identity, fill = row_pct)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = text, colour = row_pct > 50), size = 3.8, show.legend = FALSE) +
    scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "white")) +
    scale_fill_gradient(low = "white", high = "steelblue4", limits = c(0, 100),
                        name = "% of row") +
    scale_x_discrete(position = "top") +
    labs(x = "SingleR label", y = "Annotated cell type") +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(), axis.ticks = element_blank())
  ggsave(file.path(figure_dir, paste0("17_agreement_heatmap_", a$name, ".png")),
         figure_heatmap, width = 10, height = 4.5, dpi = 150)
}

message("Wrote figures to ", figure_dir)
