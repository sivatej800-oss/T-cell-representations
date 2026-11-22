# Stage 19: how far the two arms' clusterings agree with each other.
#
# Stage 17 compares each arm with SingleR. This compares the two arms directly, with
# no reference labels at all: the authors' published pipeline (Stage 03) against scGPT
# used out of the box (Stage 06's scgpt_raw), on the same 70,391 cells.
#
# The measure is the cross-tabulation of the frozen partitions from Stage 12: where
# each cluster of one arm lands in the other, in cells and as percentages in both
# directions, plus a one-line best match for every cluster.
#
# Single summary scores (ARI, NMI between the arms) were computed in an earlier version
# and dropped. They added nothing the table does not show, and misled: both are
# dominated by the largest groups, so the T-cell boundaries pulled ARI down to 0.40
# while every non-T population matched almost one to one.
#
# Each cluster is labelled with its Stage 15 broad cell type so the table can be read,
# but the measures use the cluster numbers alone. Cluster numbers carry no meaning
# across arms: Louvain numbers by size, so "cluster 3" in one arm has nothing to do
# with "cluster 3" in the other.
#
# Usage: Rscript 02_CODE/19_arm_agreement.R <workspace>

suppressPackageStartupMessages({
  library(ggplot2)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: 19_arm_agreement.R <workspace>")
results_dir <- file.path(normalizePath(args[1]), "results")
output_dir <- file.path(results_dir, "19_arm_agreement")
figure_dir <- file.path(results_dir, "figures")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

expected_cells <- 70391

# --- Step 1: the partitions and their labels ---------------------------------------
annotation <- utils::read.csv(file.path(results_dir, "15_annotation", "annotated_cluster_labels.csv.gz"),
                              colClasses = "character")
if (nrow(annotation) != expected_cells || anyDuplicated(annotation$cell_barcode)) {
  stop("Stage 15's file should hold ", expected_cells, " unique cells")
}
needed <- c("authors_pipeline_cluster", "authors_pipeline_cluster_identity",
            "scgpt_raw_cluster", "scgpt_raw_cluster_identity")
missing <- setdiff(needed, colnames(annotation))
if (length(missing) > 0) stop("Stage 15's file is missing: ", paste(missing, collapse = ", "))
if (anyNA(annotation[, needed]) || any(annotation[, needed] == "")) {
  stop("A cell is missing a cluster or an identity")
}


authors <- annotation$authors_pipeline_cluster
scgpt <- annotation$scgpt_raw_cluster


# --- Step 2: the cross-tabulation --------------------------------------------------
# Clusters are ordered by broad cell type and then by size, so that matching
# populations sit near the diagonal and a split or a merge is easy to see.
identity_order <- c("T Cells", "NK Cells", "DC", "Monocytes", "B Cells", "Plasma Cells",
                    "Undetermined")
cluster_info <- function(clusters, identities) {
  sizes <- table(clusters)
  info <- unique(data.frame(cluster = clusters, identity = identities, stringsAsFactors = FALSE))
  if (anyDuplicated(info$cluster)) stop("A cluster carries two identities")
  info$n_cells <- as.integer(sizes[info$cluster])
  info <- info[order(match(info$identity, identity_order), -info$n_cells), ]
  info$label <- sprintf("%s · %s (%s)", info$cluster, info$identity,
                        prettyNum(info$n_cells, big.mark = ","))
  info
}
authors_info <- cluster_info(authors, annotation$authors_pipeline_cluster_identity)
scgpt_info <- cluster_info(scgpt, annotation$scgpt_raw_cluster_identity)

counts <- table(factor(authors, levels = authors_info$cluster),
                factor(scgpt, levels = scgpt_info$cluster))
if (sum(counts) != expected_cells) stop("The cross-tabulation does not cover every cell")
row_pct <- 100 * sweep(counts, 1, rowSums(counts), "/")
col_pct <- 100 * sweep(counts, 2, colSums(counts), "/")

write_table <- function(matrix, file) {
  table <- as.data.frame.matrix(matrix)
  colnames(table) <- paste0("scgpt_", scgpt_info$cluster, " ", scgpt_info$identity)
  table <- cbind(authors_cluster = authors_info$cluster,
                 authors_identity = authors_info$identity, table)
  utils::write.csv(table, file.path(output_dir, file), row.names = FALSE)
}
write_table(counts, "crosstab_counts.csv")
write_table(round(row_pct, 2), "crosstab_pct_of_authors_cluster.csv")
write_table(round(col_pct, 2), "crosstab_pct_of_scgpt_cluster.csv")

# For every cluster in either arm: where most of its cells go in the other arm, and
# what share that is. This is the one-line version of the table.
best_match <- function(pct, info_from, info_to, arm_from, by_row) {
  do.call(rbind, lapply(seq_len(nrow(info_from)), function(i) {
    shares <- if (by_row) pct[i, ] else pct[, i]
    top <- order(-shares)[1:2]
    data.frame(arm = arm_from, cluster = info_from$cluster[i],
               identity = info_from$identity[i], n_cells = info_from$n_cells[i],
               best_match_cluster = info_to$cluster[top[1]],
               best_match_identity = info_to$identity[top[1]],
               best_match_pct = round(shares[top[1]], 1),
               second_match_cluster = info_to$cluster[top[2]],
               second_match_identity = info_to$identity[top[2]],
               second_match_pct = round(shares[top[2]], 1),
               stringsAsFactors = FALSE)
  }))
}
correspondence <- rbind(
  best_match(row_pct, authors_info, scgpt_info, "authors_pipeline", TRUE),
  best_match(col_pct, scgpt_info, authors_info, "scgpt_raw", FALSE)
)
utils::write.csv(correspondence, file.path(output_dir, "cluster_correspondence.csv"), row.names = FALSE)
message("\nWhere each cluster lands in the other arm:")
print(correspondence, row.names = FALSE)


# --- Step 3: the figure ------------------------------------------------------------
# Rows are the authors' clusters, columns scGPT's; each cell is the share of the row's
# cells that land in that column, so every row sums to 100%.
grid <- as.data.frame(as.table(row_pct), stringsAsFactors = FALSE)
colnames(grid) <- c("authors", "scgpt", "pct")
grid$authors <- factor(authors_info$label[match(grid$authors, authors_info$cluster)],
                       levels = rev(authors_info$label))
grid$scgpt <- factor(scgpt_info$label[match(grid$scgpt, scgpt_info$cluster)],
                     levels = scgpt_info$label)
grid$text <- ifelse(grid$pct < 0.5, "", sprintf("%.0f", grid$pct))
figure <- ggplot(grid, aes(x = scgpt, y = authors, fill = pct)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = text, colour = pct > 50), size = 3.4, show.legend = FALSE) +
  scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "white")) +
  scale_fill_gradient(low = "white", high = "steelblue4", limits = c(0, 100),
                      name = "% of authors'\ncluster") +
  scale_x_discrete(position = "top") +
  labs(x = "scGPT raw cluster", y = "Authors' pipeline cluster",
       title = "Where each authors' pipeline cluster lands in scGPT raw") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), axis.ticks = element_blank(),
        axis.text.x.top = element_text(angle = 40, hjust = 0, vjust = 0),
        plot.title = element_text(hjust = 0))
target <- file.path(figure_dir, "19_arm_crosstab_heatmap.png")
ggsave(target, figure, width = 12, height = 7.5, dpi = 150)
if (!file.exists(target) || file.info(target)$size == 0) stop("The heatmap was not written")

writeLines(toJSON(list(
  stage = "19_arm_agreement",
  arms = c("authors_pipeline_cluster", "scgpt_raw_cluster"),
  source = "results/15_annotation/annotated_cluster_labels.csv.gz",
  n_cells = length(authors),
  clusters = list(authors_pipeline = length(unique(authors)), scgpt_raw = length(unique(scgpt))),
  r_version = paste(R.version$major, R.version$minor, sep = ".")
), auto_unbox = TRUE, pretty = TRUE, digits = 12), file.path(output_dir, "arm_agreement_audit.json"))

message("\nStage 19 complete: ", output_dir)
