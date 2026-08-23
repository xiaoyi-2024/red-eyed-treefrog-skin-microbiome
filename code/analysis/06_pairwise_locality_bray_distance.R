suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
})

set.seed(20260814)
input_file <- Sys.getenv("PHYLOSEQ_INPUT",
                         unset = "data/processed/ps_clean_47_rebuilt_from_01_filtered.rds")
result_dir <- Sys.getenv("ANALYSIS_RESULT_DIR", unset = "results/pairwise_bray_distance")
community_label <- Sys.getenv("COMMUNITY_LABEL", unset = "whole community")
figure_dir <- file.path(result_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

ps <- readRDS(input_file)
otu <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) otu <- t(otu)
tax <- as(tax_table(ps), "matrix")[colnames(otu), , drop = FALSE]
meta <- data.frame(sample_data(ps), check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(nrow(otu) == 47L, identical(rownames(otu), rownames(meta)))

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

site_key <- c(
  "Lost Iguana Resort, La Fortuna, San Carlos" = "Lost Iguana",
  "Veragua Rainforest, Las Brisas de Veragua" = "Veragua",
  "PN Altos de Campana" = "Altos de Campana",
  "PN Soberanía, Las Cumbres" = "Soberanía"
)
site <- unname(site_key[meta$Locality])
site_levels <- c("Lost Iguana", "Veragua", "Altos de Campana", "Soberanía")
site <- factor(site, levels = site_levels)
stopifnot(!anyNA(site))

bray <- vegan::vegdist(relative, method = "bray")
dm <- as.matrix(bray)

# Long table of unique sample pairs. Distances sharing a sample are not
# independent; this table is used for descriptive summaries and plotting only.
pair_idx <- which(upper.tri(dm), arr.ind = TRUE)
sample_pairs <- data.frame(
  Sample_1 = rownames(dm)[pair_idx[, 1]], Sample_2 = colnames(dm)[pair_idx[, 2]],
  Site_1 = as.character(site[pair_idx[, 1]]),
  Site_2 = as.character(site[pair_idx[, 2]]),
  Bray_Curtis = dm[pair_idx], stringsAsFactors = FALSE
)
sample_pairs$Pair <- ifelse(
  match(sample_pairs$Site_1, site_levels) <= match(sample_pairs$Site_2, site_levels),
  paste(sample_pairs$Site_1, sample_pairs$Site_2, sep = " vs "),
  paste(sample_pairs$Site_2, sample_pairs$Site_1, sep = " vs ")
)

site_pairs <- combn(site_levels, 2, simplify = FALSE)
pairwise_results <- do.call(rbind, lapply(seq_along(site_pairs), function(i) {
  z <- site_pairs[[i]]
  idx <- which(site %in% z)
  dat <- data.frame(Site = droplevels(site[idx]))
  dsub <- vegan::vegdist(relative[idx, , drop = FALSE], method = "bray")
  set.seed(20260814 + i)
  per <- vegan::adonis2(dsub ~ Site, data = dat, permutations = 9999)
  bd <- vegan::betadisper(dsub, dat$Site, type = "centroid", bias.adjust = TRUE)
  set.seed(20260814 + 100 + i)
  bdp <- permutest(bd, permutations = 9999)
  cross <- sample_pairs[(sample_pairs$Site_1 == z[1] & sample_pairs$Site_2 == z[2]) |
                        (sample_pairs$Site_1 == z[2] & sample_pairs$Site_2 == z[1]), ]
  data.frame(
    Site_1 = z[1], Site_2 = z[2], Pair = paste(z, collapse = " vs "),
    N_site_1 = sum(site == z[1]), N_site_2 = sum(site == z[2]),
    Cross_sample_pairs = nrow(cross),
    Mean_Bray_Curtis_distance = mean(cross$Bray_Curtis),
    Median_Bray_Curtis_distance = median(cross$Bray_Curtis),
    SD_Bray_Curtis_distance = sd(cross$Bray_Curtis),
    PERMANOVA_pseudo_F = per$F[1], PERMANOVA_R2 = per$R2[1],
    PERMANOVA_p = per$`Pr(>F)`[1],
    Betadisper_F = bdp$tab$F[1], Betadisper_p = bdp$tab$`Pr(>F)`[1]
  )
}))
pairwise_results$PERMANOVA_FDR <- p.adjust(pairwise_results$PERMANOVA_p, method = "BH")
pairwise_results$Betadisper_FDR <- p.adjust(pairwise_results$Betadisper_p, method = "BH")
pairwise_results$Interpretation <- ifelse(
  pairwise_results$PERMANOVA_FDR < 0.05 & pairwise_results$Betadisper_FDR >= 0.05,
  "Centroid difference supported without detectable dispersion difference",
  ifelse(pairwise_results$PERMANOVA_FDR < 0.05 & pairwise_results$Betadisper_FDR < 0.05,
         "PERMANOVA significant; dispersion difference may contribute",
         "No FDR-significant pairwise centroid difference")
)

within_results <- do.call(rbind, lapply(site_levels, function(s) {
  z <- sample_pairs[sample_pairs$Site_1 == s & sample_pairs$Site_2 == s, ]
  data.frame(Site = s, Sample_n = sum(site == s), Within_sample_pairs = nrow(z),
             Mean_within_Bray_Curtis_distance = mean(z$Bray_Curtis),
             Median_within_Bray_Curtis_distance = median(z$Bray_Curtis),
             SD_within_Bray_Curtis_distance = sd(z$Bray_Curtis))
}))

write.csv(pairwise_results, file.path(result_dir, "pairwise_locality_bray_permanova.csv"), row.names = FALSE)
write.csv(within_results, file.path(result_dir, "within_locality_bray_distances.csv"), row.names = FALSE)
write.csv(sample_pairs, file.path(result_dir, "descriptive_sample_pair_distances.csv"), row.names = FALSE)
write.csv(data.frame(Parameter = c("Input_file", "Community", "Samples", "Original_ASVs", "Aggregated_genera",
                                   "Retained_genera", "Prevalence_min_samples",
                                   "PERMANOVA_permutations", "Betadisper_permutations",
                                   "Multiple_testing"),
                     Value = c(input_file, community_label, nrow(otu), ncol(otu), ncol(genus_counts), ncol(filtered),
                               prevalence_min, 9999, 9999, "BH-FDR across six locality pairs")),
          file.path(result_dir, "pairwise_bray_distance_parameters.csv"), row.names = FALSE)

# Symmetric heat-map table: diagonal is mean within-locality distance;
# off-diagonals are mean cross-locality distances.
heat <- expand.grid(Site_y = site_levels, Site_x = site_levels, stringsAsFactors = FALSE)
heat$Distance <- NA_real_; heat$FDR <- NA_real_
for (i in seq_len(nrow(heat))) {
  a <- heat$Site_y[i]; b <- heat$Site_x[i]
  if (a == b) {
    heat$Distance[i] <- within_results$Mean_within_Bray_Curtis_distance[within_results$Site == a]
  } else {
    hit <- pairwise_results[(pairwise_results$Site_1 == a & pairwise_results$Site_2 == b) |
                            (pairwise_results$Site_1 == b & pairwise_results$Site_2 == a), ]
    heat$Distance[i] <- hit$Mean_Bray_Curtis_distance
    heat$FDR[i] <- hit$PERMANOVA_FDR
  }
}
heat$Site_x <- factor(heat$Site_x, levels = site_levels)
heat$Site_y <- factor(heat$Site_y, levels = rev(site_levels))
heat$Label <- ifelse(is.na(heat$FDR), sprintf("%.3f\nwithin", heat$Distance),
                     sprintf("%.3f%s", heat$Distance,
                             ifelse(heat$FDR < 0.001, "***",
                                    ifelse(heat$FDR < 0.01, "**",
                                           ifelse(heat$FDR < 0.05, "*", "")))))

pair_order <- pairwise_results$Pair[order(pairwise_results$Mean_Bray_Curtis_distance, decreasing = TRUE)]
between_pairs <- sample_pairs[sample_pairs$Site_1 != sample_pairs$Site_2, ]
between_pairs$Pair <- factor(between_pairs$Pair, levels = pair_order)
pairwise_results$Pair <- factor(pairwise_results$Pair, levels = pair_order)
short_site <- c("Lost Iguana" = "Lost", "Veragua" = "Veragua",
                "Altos de Campana" = "Campana", "Soberanía" = "Soberanía")
short_pair <- function(x) {
  parts <- strsplit(as.character(x), " vs ", fixed = TRUE)
  vapply(parts, function(z) paste(short_site[z], collapse = " vs "), character(1))
}
short_pair_levels <- short_pair(pair_order)
between_pairs$Pair_short <- factor(short_pair(between_pairs$Pair), levels = short_pair_levels)
pairwise_results$Pair_short <- factor(short_pair(pairwise_results$Pair), levels = short_pair_levels)

theme_pub <- theme_classic(base_size = 8, base_family = "sans") +
  theme(axis.line = element_line(linewidth = 0.4), axis.ticks = element_line(linewidth = 0.4),
        plot.title = element_text(size = 8.5, face = "bold"),
        plot.subtitle = element_text(size = 7),
        legend.title = element_text(size = 7), legend.text = element_text(size = 6.5))

p_a <- ggplot(heat, aes(Site_x, Site_y, fill = Distance)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = Label), size = 2.25, lineheight = 0.9) +
  scale_fill_viridis_c(option = "C", direction = -1,
                       name = "Mean distance") +
  labs(title = "a  Mean Bray-Curtis distance",
       subtitle = "Diagonal: within-locality; off-diagonal: between-locality",
       x = NULL, y = NULL) +
  scale_x_discrete(labels = short_site) + scale_y_discrete(labels = short_site) +
  theme_minimal(base_size = 8, base_family = "sans") +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1),
        plot.title = element_text(size = 8.5, face = "bold"), plot.subtitle = element_text(size = 6.5),
        legend.title = element_text(size = 7), legend.text = element_text(size = 6.5))

p_b <- ggplot(between_pairs, aes(Pair_short, Bray_Curtis)) +
  geom_violin(fill = "#D8E6F3", colour = "#56758F", linewidth = 0.4,
              scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", linewidth = 0.35) +
  geom_point(data = pairwise_results,
             aes(Pair_short, Mean_Bray_Curtis_distance), inherit.aes = FALSE,
             shape = 23, size = 2.3, fill = "#D55E00", colour = "black") +
  coord_flip() +
  labs(title = "c  Between-locality distance distributions",
       subtitle = "Diamonds are means; sample-pair values are descriptive and non-independent",
       x = NULL, y = "Bray-Curtis distance") + theme_pub

pairwise_results$Sig <- ifelse(pairwise_results$PERMANOVA_FDR < 0.001, "q < 0.001",
                               ifelse(pairwise_results$PERMANOVA_FDR < 0.01, "q < 0.01",
                                      ifelse(pairwise_results$PERMANOVA_FDR < 0.05, "q < 0.05", "q >= 0.05")))
p_c <- ggplot(pairwise_results, aes(Pair_short, PERMANOVA_R2, colour = Sig)) +
  geom_segment(aes(xend = Pair_short, y = 0, yend = PERMANOVA_R2), colour = "grey75", linewidth = 0.45) +
  geom_point(size = 2.8) + coord_flip() +
  scale_colour_manual(values = c("q < 0.001" = "#7A1F5C", "q < 0.01" = "#0072B2",
                                 "q < 0.05" = "#009E73", "q >= 0.05" = "grey55"), drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "b  PERMANOVA effect sizes",
       subtitle = "BH-FDR-adjusted significance",
       x = NULL, y = expression(PERMANOVA~R^2), colour = NULL) + theme_pub +
  theme(legend.position = "bottom")

figure <- ((p_a | p_c) / p_b) + plot_layout(heights = c(1.0, 1.05)) +
  plot_annotation(
    title = paste0("Pairwise locality comparison: ", community_label),
    subtitle = paste0("Genus relative abundance; 20% prevalence filter; Bray-Curtis dissimilarity. ",
                      "PERMANOVA and betadisper use 9,999 permutations."),
    caption = "Heat-map symbols denote pairwise PERMANOVA after BH-FDR correction: *q<0.05, **q<0.01, ***q<0.001.",
    theme = theme(plot.title = element_text(size = 10, face = "bold"),
                  plot.subtitle = element_text(size = 6.8),
                  plot.caption = element_text(size = 6.2, hjust = 0))
  )

save_plot <- function(plot, stem, width_mm = 183, height_mm = 155) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  ggsave(paste0(stem, ".png"), plot, width = w, height = h, dpi = 600, bg = "white")
  ragg::agg_tiff(paste0(stem, ".tiff"), width = w, height = h, units = "in",
                 res = 600, background = "white"); print(plot); dev.off()
  svglite::svglite(paste0(stem, ".svg"), width = w, height = h); print(plot); dev.off()
  grDevices::pdf(paste0(stem, ".pdf"), width = w, height = h,
                 family = "Helvetica", useDingbats = FALSE); print(plot); dev.off()
}
save_plot(figure, file.path(figure_dir, "figure_pairwise_bray_distance"))
print(pairwise_results)
print(within_results)
