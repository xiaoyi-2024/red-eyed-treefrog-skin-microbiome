suppressPackageStartupMessages({
  library(phyloseq)
  library(cluster)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(svglite)
  library(ragg)
})

set.seed(20260814)
input_file <- Sys.getenv("PHYLOSEQ_INPUT",
                         unset = "data/processed/ps_clean_47_rebuilt_from_01_filtered.rds")
result_dir <- Sys.getenv("ANALYSIS_RESULT_DIR", unset = "results/pam_clustering")
community_label <- Sys.getenv("COMMUNITY_LABEL", unset = "whole bacterial community")
figure_dir <- file.path(result_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(input_file))

ps <- readRDS(input_file)
otu <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) otu <- t(otu)
tax <- as(tax_table(ps), "matrix")[colnames(otu), , drop = FALSE]
meta <- data.frame(sample_data(ps), check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(nrow(otu) == 47L, identical(rownames(otu), rownames(meta)))

# Aggregate ASVs to genus. Missing genera are grouped by the lowest available
# family label instead of being silently deleted.
genus <- as.character(tax[, "Genus"])
missing_genus <- is.na(genus) | trimws(genus) == ""
family <- as.character(tax[, "Family"])
family[is.na(family) | trimws(family) == ""] <- "Bacteria"
genus[missing_genus] <- paste0("Unclassified_", family[missing_genus])
genus_counts <- t(rowsum(t(otu), group = genus, reorder = FALSE))

# Main solution retains genera occurring in >=20% of samples. This threshold
# was selected before choosing k to reduce sparse-ASV singleton clusters.
prevalence_min <- ceiling(0.20 * nrow(genus_counts))
keep <- colSums(genus_counts > 0) >= prevalence_min
filtered <- genus_counts[, keep, drop = FALSE]
clr <- log(filtered + 1)
clr <- clr - rowMeans(clr)
aitchison <- dist(clr, method = "euclidean")

site_key <- c(
  "Lost Iguana Resort, La Fortuna, San Carlos" = "Lost Iguana",
  "Veragua Rainforest, Las Brisas de Veragua" = "Veragua",
  "PN Altos de Campana" = "Altos de Campana",
  "PN Soberanía, Las Cumbres" = "Soberanía"
)
site <- unname(site_key[meta$Locality])
stopifnot(!anyNA(site))

# Subsampling stability: repeatedly cluster 80% of samples and, for each
# original cluster, retain the best Jaccard overlap with a resampled cluster.
cluster_stability <- function(original, k, B = 500L, fraction = 0.80) {
  out <- matrix(NA_real_, nrow = B, ncol = k)
  for (b in seq_len(B)) {
    idx <- sort(sample(seq_along(original), ceiling(fraction * length(original))))
    db <- as.dist(as.matrix(aitchison)[idx, idx, drop = FALSE])
    cb <- pam(db, k = k, diss = TRUE)$clustering
    for (g in seq_len(k)) {
      original_members <- idx[original[idx] == g]
      if (!length(original_members)) next
      out[b, g] <- max(vapply(seq_len(k), function(h) {
        boot_members <- idx[cb == h]
        length(intersect(original_members, boot_members)) /
          length(union(original_members, boot_members))
      }, numeric(1)))
    }
  }
  c(mean = mean(out, na.rm = TRUE), minimum_cluster = min(colMeans(out, na.rm = TRUE)))
}

k_values <- 2:8
fits <- lapply(k_values, function(k) pam(aitchison, k = k, diss = TRUE))
evaluation <- do.call(rbind, Map(function(k, fit) {
  sizes <- as.integer(table(factor(fit$clustering, levels = seq_len(k))))
  stab <- cluster_stability(fit$clustering, k)
  data.frame(k = k, mean_silhouette = fit$silinfo$avg.width,
             mean_bootstrap_Jaccard = stab["mean"],
             minimum_cluster_Jaccard = stab["minimum_cluster"],
             minimum_cluster_size = min(sizes),
             cluster_sizes = paste(sizes, collapse = ";"))
}, k_values, fits))

# Avoid solutions containing clusters smaller than five samples. Among the
# eligible solutions, choose the largest mean silhouette width.
eligible <- evaluation$minimum_cluster_size >= 5
if (!any(eligible)) stop("No PAM solution has a minimum cluster size >= 5")
selected_row <- which.max(ifelse(eligible, evaluation$mean_silhouette, -Inf))
selected_k <- evaluation$k[selected_row]
fit <- fits[[match(selected_k, k_values)]]
cluster_id <- factor(fit$clustering, levels = seq_len(selected_k),
                     labels = paste0("PAM", seq_len(selected_k)))

sil <- silhouette(fit$clustering, aitchison)
pca <- prcomp(clr, center = FALSE, scale. = FALSE)
scores <- data.frame(Sample_ID = rownames(clr), PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                     Site = site, Cluster = cluster_id,
                     Silhouette_width = sil[, "sil_width"],
                     stringsAsFactors = FALSE)
scores$PC1_percent <- 100 * summary(pca)$importance[2, 1]
scores$PC2_percent <- 100 * summary(pca)$importance[2, 2]

tab <- table(scores$Site, scores$Cluster)
set.seed(20260814)
fisher <- fisher.test(tab, simulate.p.value = TRUE, B = 99999)
association <- data.frame(
  Test = "Fisher exact test with Monte Carlo simulation",
  Site_n = length(unique(site)), Sample_n = nrow(scores), Selected_k = selected_k,
  p = fisher$p.value, Replicates = 99999L,
  Note = "Sample-level association; clustering strength must be evaluated separately"
)

assignment <- cbind(scores, Locality_original = meta$Locality,
                    Is_medoid = rownames(clr) %in% fit$medoids)
write.csv(evaluation, file.path(result_dir, "pam_k_evaluation.csv"), row.names = FALSE)
write.csv(assignment, file.path(result_dir, "pam_sample_assignments.csv"), row.names = FALSE)
write.csv(as.data.frame.matrix(tab), file.path(result_dir, "pam_cluster_by_site_table.csv"))
write.csv(association, file.path(result_dir, "pam_site_association.csv"), row.names = FALSE)
write.csv(data.frame(Parameter = c("Input_file", "Community", "Samples", "Original_ASVs", "Aggregated_genera",
                                   "Retained_genera", "Prevalence_min_samples", "Pseudocount",
                                   "Selected_k", "Bootstrap_replicates"),
                     Value = c(input_file, community_label, nrow(otu), ncol(otu), ncol(genus_counts), ncol(filtered),
                               prevalence_min, 1, selected_k, 500)),
          file.path(result_dir, "pam_analysis_parameters.csv"), row.names = FALSE)

cluster_pal <- c("PAM1" = "#0072B2", "PAM2" = "#D55E00", "PAM3" = "#009E73",
                 "PAM4" = "#CC79A7", "PAM5" = "#E69F00", "PAM6" = "#56B4E9",
                 "PAM7" = "#000000", "PAM8" = "#999999")
site_shapes <- c("Lost Iguana" = 21, "Veragua" = 22,
                 "Altos de Campana" = 24, "Soberanía" = 23)
theme_pub <- theme_classic(base_size = 8, base_family = "sans") +
  theme(axis.line = element_line(linewidth = 0.4), axis.ticks = element_line(linewidth = 0.4),
        plot.title = element_text(size = 9, face = "bold"),
        legend.title = element_text(size = 7), legend.text = element_text(size = 6.5))

p_a <- ggplot(scores, aes(PC1, PC2, fill = Cluster, shape = Site)) +
  stat_ellipse(data = scores,
               aes(x = PC1, y = PC2, colour = Cluster, group = Cluster),
               inherit.aes = FALSE, type = "norm", linewidth = 0.55,
               linetype = "dashed", show.legend = FALSE) +
  geom_point(size = 3.0, stroke = 0.55, colour = "black") +
  scale_fill_manual(values = cluster_pal) + scale_colour_manual(values = cluster_pal) +
  scale_shape_manual(values = site_shapes) +
  guides(fill = guide_legend(override.aes = list(shape = 21, colour = "black",
                                                  fill = unname(cluster_pal[levels(cluster_id)])))) +
  labs(title = "a  CLR-PCA representation",
       subtitle = "Colour: PAM cluster; shape: sampling locality",
       x = sprintf("PC1 (%.1f%%)", scores$PC1_percent[1]),
       y = sprintf("PC2 (%.1f%%)", scores$PC2_percent[1])) + theme_pub

eval_long <- rbind(
  data.frame(k = evaluation$k, Metric = "Mean silhouette", Value = evaluation$mean_silhouette),
  data.frame(k = evaluation$k, Metric = "Bootstrap Jaccard", Value = evaluation$mean_bootstrap_Jaccard)
)
p_b <- ggplot(eval_long, aes(k, Value, colour = Metric, shape = Metric)) +
  geom_line(linewidth = 0.65) + geom_point(size = 2.5) +
  geom_vline(xintercept = selected_k, linetype = "dotted", colour = "grey30") +
  annotate("text", x = selected_k + 0.12, y = 0.98, label = paste0("selected k = ", selected_k),
           hjust = 0, vjust = 1, size = 2.2) +
  scale_colour_manual(values = c("Mean silhouette" = "#7A1F5C",
                                 "Bootstrap Jaccard" = "#0072B2")) +
  scale_shape_manual(values = c("Mean silhouette" = 16, "Bootstrap Jaccard" = 17)) +
  scale_x_continuous(breaks = k_values) + scale_y_continuous(limits = c(0, 1)) +
  labs(title = "b  Cluster separation and stability", x = "Number of clusters (k)",
       y = "Metric value", colour = NULL, shape = NULL) + theme_pub +
  theme(legend.position = "bottom", legend.direction = "horizontal")

# Distance heat map ordered by PAM cluster and then locality.
ord <- order(scores$Cluster, scores$Site)
dm <- as.matrix(aitchison)[ord, ord]
heat <- as.data.frame(as.table(dm), stringsAsFactors = FALSE)
names(heat) <- c("Sample_y", "Sample_x", "Distance")
levels_order <- rownames(dm)
heat$Sample_x <- factor(heat$Sample_x, levels = levels_order)
heat$Sample_y <- factor(heat$Sample_y, levels = rev(levels_order))
p_c <- ggplot(heat, aes(Sample_x, Sample_y, fill = Distance)) +
  geom_raster() +
  scale_fill_viridis_c(option = "C", name = "Aitchison\ndistance") +
  labs(title = "c  Sample-to-sample distance matrix", x = "Samples", y = "Samples") +
  theme_minimal(base_size = 8, base_family = "sans") +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        plot.title = element_text(size = 9, face = "bold"),
        legend.title = element_text(size = 7), legend.text = element_text(size = 6.5))

figure <- ((p_a | p_b) / p_c) + plot_layout(heights = c(1, 1.12)) +
  plot_annotation(
    title = paste0("PAM clustering of the ", community_label),
    subtitle = sprintf(paste0("Genus-level CLR-Aitchison distance; genera present in >=20%% of samples; pseudocount = 1. ",
                              "Selected k = %d; mean silhouette = %.3f; Fisher p %s."),
                       selected_k, evaluation$mean_silhouette[selected_row],
                       ifelse(fisher$p.value <= 1e-4, "< 0.0001",
                              paste0("= ", formatC(fisher$p.value, digits = 4, format = "f")))),
    theme = theme(plot.title = element_text(size = 11, face = "bold"),
                  plot.subtitle = element_text(size = 7.2))
  )

save_plot <- function(plot, stem, width_mm = 183, height_mm = 175) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  ggsave(paste0(stem, ".png"), plot, width = w, height = h, dpi = 600, bg = "white")
  ragg::agg_tiff(paste0(stem, ".tiff"), width = w, height = h, units = "in",
                 res = 600, background = "white"); print(plot); dev.off()
  svglite::svglite(paste0(stem, ".svg"), width = w, height = h); print(plot); dev.off()
  grDevices::pdf(paste0(stem, ".pdf"), width = w, height = h,
                 family = "Helvetica", useDingbats = FALSE); print(plot); dev.off()
}
save_plot(figure, file.path(figure_dir, "figure_pam_clustering_all_samples"))

cat(sprintf("Selected k = %d\n", selected_k))
print(evaluation)
print(tab)
print(association)
