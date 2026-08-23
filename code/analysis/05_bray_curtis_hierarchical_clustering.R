suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(ggplot2)
  library(ggdendro)
  library(patchwork)
  library(svglite)
  library(ragg)
})

input_file <- Sys.getenv("PHYLOSEQ_INPUT",
                         unset = "data/processed/ps_clean_47_rebuilt_from_01_filtered.rds")
result_dir <- Sys.getenv("ANALYSIS_RESULT_DIR", unset = "results/bray_curtis_clustering")
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

# Match the retained feature set used in the primary PAM analysis: aggregate
# ASVs to genus and retain genera occurring in at least 20% (10/47) samples.
genus <- as.character(tax[, "Genus"])
missing <- is.na(genus) | trimws(genus) == ""
family <- as.character(tax[, "Family"])
family[is.na(family) | trimws(family) == ""] <- "Bacteria"
genus[missing] <- paste0("Unclassified_", family[missing])
genus_counts <- t(rowsum(t(otu), group = genus, reorder = FALSE))
prevalence_min <- ceiling(0.20 * nrow(genus_counts))
keep <- colSums(genus_counts > 0) >= prevalence_min
filtered <- genus_counts[, keep, drop = FALSE]
relative <- filtered / rowSums(filtered)

bray <- vegan::vegdist(relative, method = "bray")
tree <- hclust(bray, method = "average")
cophen <- cophenetic(tree)
pearson <- cor(as.vector(bray), as.vector(cophen), method = "pearson")
spearman <- cor(as.vector(bray), as.vector(cophen), method = "spearman")

site_key <- c(
  "Lost Iguana Resort, La Fortuna, San Carlos" = "Lost Iguana",
  "Veragua Rainforest, Las Brisas de Veragua" = "Veragua",
  "PN Altos de Campana" = "Altos de Campana",
  "PN Soberanía, Las Cumbres" = "Soberanía"
)
site <- unname(site_key[meta$Locality])
names(site) <- rownames(meta)
stopifnot(!anyNA(site))

ordered_samples <- tree$labels[tree$order]
sample_order <- data.frame(
  Tree_order = seq_along(ordered_samples), Sample_ID = ordered_samples,
  Site = unname(site[ordered_samples]), stringsAsFactors = FALSE
)
metrics <- data.frame(
  Input_file = input_file, Community = community_label,
  Samples = nrow(relative), Original_ASVs = ncol(otu),
  Aggregated_genera = ncol(genus_counts), Retained_genera = ncol(relative),
  Prevalence_threshold_percent = 20, Prevalence_min_samples = prevalence_min,
  Distance = "Bray-Curtis on genus relative abundance",
  Linkage = "average (UPGMA)",
  Cophenetic_Pearson = pearson, Cophenetic_Spearman = spearman,
  Fixed_k_or_tree_cut = "none"
)
write.csv(metrics, file.path(result_dir, "bray_curtis_tree_metrics.csv"), row.names = FALSE)
write.csv(sample_order, file.path(result_dir, "bray_curtis_tree_sample_order.csv"), row.names = FALSE)
write.csv(as.matrix(bray), file.path(result_dir, "bray_curtis_distance_matrix.csv"))

dd <- ggdendro::dendro_data(as.dendrogram(tree), type = "rectangle")
segments <- ggdendro::segment(dd)
labels <- ggdendro::label(dd)
labels$Site <- unname(site[labels$label])

site_pal <- c("Lost Iguana" = "#CC79A7", "Veragua" = "#0072B2",
              "Altos de Campana" = "#009E73", "Soberanía" = "#E69F00")
theme_pub <- theme_classic(base_size = 8, base_family = "sans") +
  theme(axis.line = element_line(linewidth = 0.4), axis.ticks = element_line(linewidth = 0.4),
        plot.title = element_text(size = 9, face = "bold"),
        legend.title = element_text(size = 7), legend.text = element_text(size = 6.5))

p_tree <- ggplot(segments) +
  geom_segment(aes(x, y, xend = xend, yend = yend), colour = "#3F4850",
               linewidth = 0.38, lineend = "round") +
  geom_point(data = labels, aes(x, y = 0, fill = Site), shape = 21,
             size = 2.2, stroke = 0.35, colour = "white") +
  scale_fill_manual(values = site_pal) +
  scale_x_continuous(breaks = labels$x, labels = labels$label,
                     expand = expansion(mult = c(0.012, 0.012))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.06))) +
  labs(title = "a  Bray-Curtis hierarchical tree",
       subtitle = "Average linkage; coloured tips indicate sampling locality",
       x = NULL, y = "Bray-Curtis fusion height", fill = "Locality") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.0),
        axis.ticks.x = element_blank(), legend.position = "right")

dm <- as.matrix(bray)[ordered_samples, ordered_samples]
heat <- as.data.frame(as.table(dm), stringsAsFactors = FALSE)
names(heat) <- c("Sample_y", "Sample_x", "Distance")
heat$Sample_x <- factor(heat$Sample_x, levels = ordered_samples)
heat$Sample_y <- factor(heat$Sample_y, levels = rev(ordered_samples))

# Two thin locality annotation strips make the sample ordering readable without
# imposing or displaying any arbitrary cluster cut.
strip_x <- data.frame(Sample = factor(ordered_samples, levels = ordered_samples),
                      Site = unname(site[ordered_samples]), y = 1)
p_strip <- ggplot(strip_x, aes(Sample, y, fill = Site)) +
  geom_tile() + scale_fill_manual(values = site_pal, guide = "none") +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Locality") +
  theme_void(base_family = "sans") +
  theme(axis.title.y = element_text(size = 7, angle = 90, margin = margin(r = 5)))

p_heat <- ggplot(heat, aes(Sample_x, Sample_y, fill = Distance)) +
  geom_raster() +
  scale_fill_viridis_c(option = "C", limits = c(0, max(heat$Distance)),
                       name = "Bray-Curtis\ndissimilarity") +
  labs(title = "b  Tree-ordered sample dissimilarity matrix",
       subtitle = sprintf("Cophenetic Pearson r = %.3f; Spearman rho = %.3f",
                          pearson, spearman),
       x = "Samples ordered by hierarchical tree", y = "Samples") +
  theme_minimal(base_size = 8, base_family = "sans") +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        plot.title = element_text(size = 9, face = "bold"),
        plot.subtitle = element_text(size = 7),
        legend.title = element_text(size = 7), legend.text = element_text(size = 6.5))

lower <- p_strip / p_heat + plot_layout(heights = c(0.045, 1))
figure <- p_tree / lower + plot_layout(heights = c(0.92, 1.08)) +
  plot_annotation(
    title = paste0("Bray-Curtis clustering of the ", community_label),
    subtitle = paste0("Genus relative abundance; genera present in >=20% of samples. ",
                      "The complete hierarchy is shown without selecting k or cutting the tree."),
    theme = theme(plot.title = element_text(size = 11, face = "bold"),
                  plot.subtitle = element_text(size = 7.2))
  )

save_plot <- function(plot, stem, width_mm = 183, height_mm = 205) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  ggsave(paste0(stem, ".png"), plot, width = w, height = h, dpi = 600, bg = "white")
  ragg::agg_tiff(paste0(stem, ".tiff"), width = w, height = h, units = "in",
                 res = 600, background = "white"); print(plot); dev.off()
  svglite::svglite(paste0(stem, ".svg"), width = w, height = h); print(plot); dev.off()
  grDevices::pdf(paste0(stem, ".pdf"), width = w, height = h,
                 family = "Helvetica", useDingbats = FALSE); print(plot); dev.off()
}
save_plot(figure, file.path(figure_dir, "figure_bray_curtis_hierarchical_clustering"))
print(metrics)
