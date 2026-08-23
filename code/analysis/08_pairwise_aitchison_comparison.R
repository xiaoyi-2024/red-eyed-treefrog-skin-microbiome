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
result_dir <- Sys.getenv("ANALYSIS_RESULT_DIR", unset = "results/pairwise_aitchison_comparison")
community_label <- Sys.getenv("COMMUNITY_LABEL", unset = "whole bacterial community")
figure_dir <- file.path(result_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

ps <- readRDS(input_file)
otu <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) otu <- t(otu)
tax <- as(tax_table(ps), "matrix")[colnames(otu), , drop = FALSE]
meta <- data.frame(sample_data(ps), check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(nrow(otu) == 47L, identical(rownames(otu), rownames(meta)))

genus <- as.character(tax[, "Genus"])
missing_genus <- is.na(genus) | trimws(genus) == ""
family <- as.character(tax[, "Family"])
family[is.na(family) | trimws(family) == ""] <- "Bacteria"
genus[missing_genus] <- paste0("Unclassified_", family[missing_genus])
genus_counts <- t(rowsum(t(otu), group = genus, reorder = FALSE))
prevalence_min <- ceiling(0.20 * nrow(genus_counts))
keep <- colSums(genus_counts > 0) >= prevalence_min
filtered <- genus_counts[, keep, drop = FALSE]

site_key <- c(
  "Lost Iguana Resort, La Fortuna, San Carlos" = "Lost",
  "Veragua Rainforest, Las Brisas de Veragua" = "Veragua",
  "PN Altos de Campana" = "Campana",
  "PN Soberanía, Las Cumbres" = "Soberanía"
)
site_levels <- c("Lost", "Veragua", "Campana", "Soberanía")
site <- factor(unname(site_key[meta$Locality]), levels = site_levels)
stopifnot(!anyNA(site))

# Ordinary Aitchison: pseudocount 1, CLR, Euclidean distance.
clr <- log(filtered + 1)
clr <- clr - rowMeans(clr)
# Robust Aitchison: vegan rCLR with matrix completion, then Euclidean distance.
rclr <- vegan::decostand(filtered, method = "rclr")
distances <- list(
  `CLR-Aitchison` = as.matrix(dist(clr)),
  `Robust Aitchison` = as.matrix(dist(rclr))
)

site_pairs <- combn(site_levels, 2, simplify = FALSE)
nperm <- 9999L

extract_cross_pairs <- function(D, metric) {
  do.call(rbind, lapply(site_pairs, function(pair) {
    i <- which(site == pair[1]); j <- which(site == pair[2])
    z <- expand.grid(Sample_1 = rownames(D)[i], Sample_2 = rownames(D)[j],
                     stringsAsFactors = FALSE)
    z$Site_1 <- pair[1]; z$Site_2 <- pair[2]
    z$Pair <- paste(pair, collapse = " vs ")
    z$Metric <- metric
    z$Distance <- D[cbind(match(z$Sample_1, rownames(D)),
                          match(z$Sample_2, colnames(D)))]
    z
  }))
}

pairwise_distances <- do.call(rbind, Map(extract_cross_pairs, distances, names(distances)))
write.csv(pairwise_distances, file.path(result_dir, "cross_locality_sample_pair_distances.csv"),
          row.names = FALSE)

distance_summary <- do.call(rbind, lapply(split(pairwise_distances,
                                                interaction(pairwise_distances$Metric,
                                                            pairwise_distances$Pair,
                                                            drop = TRUE)), function(z) {
  data.frame(Metric = z$Metric[1], Pair = z$Pair[1],
             Site_1 = z$Site_1[1], Site_2 = z$Site_2[1],
             N_cross_sample_pairs = nrow(z), Mean_distance = mean(z$Distance),
             Median_distance = median(z$Distance), SD_distance = sd(z$Distance),
             Q1 = unname(quantile(z$Distance, 0.25)),
             Q3 = unname(quantile(z$Distance, 0.75)),
             Note = "Cross-sample pairs are descriptive and non-independent")
}))
rownames(distance_summary) <- NULL
write.csv(distance_summary, file.path(result_dir, "pairwise_distance_summary.csv"), row.names = FALSE)

mean_distance_matrix <- function(D, metric) {
  out <- matrix(NA_real_, length(site_levels), length(site_levels),
                dimnames = list(site_levels, site_levels))
  for (a in site_levels) for (b in site_levels) {
    ia <- which(site == a); ib <- which(site == b)
    if (a == b) {
      block <- D[ia, ib, drop = FALSE]
      out[a, b] <- mean(block[upper.tri(block)])
    } else {
      out[a, b] <- mean(D[ia, ib, drop = FALSE])
    }
  }
  data.frame(Metric = metric, Site_1 = rep(site_levels, each = length(site_levels)),
             Site_2 = rep(site_levels, times = length(site_levels)),
             Mean_distance = as.vector(t(out)), row.names = NULL)
}
heatmap_data <- do.call(rbind, Map(mean_distance_matrix, distances, names(distances)))
write.csv(heatmap_data, file.path(result_dir, "mean_distance_matrices_long.csv"), row.names = FALSE)

pairwise_tests <- do.call(rbind, lapply(seq_along(distances), function(m) {
  D <- distances[[m]]; metric <- names(distances)[m]
  do.call(rbind, lapply(seq_along(site_pairs), function(i) {
    pair <- site_pairs[[i]]; take <- which(site %in% pair)
    group <- droplevels(site[take])
    dsub <- as.dist(D[take, take, drop = FALSE])
    set.seed(20262000 + m * 100 + i)
    ad <- vegan::adonis2(dsub ~ group, permutations = nperm)
    bd <- vegan::betadisper(dsub, group, type = "centroid", bias.adjust = TRUE)
    set.seed(20263000 + m * 100 + i)
    bp <- vegan::permutest(bd, permutations = nperm)
    data.frame(
      Metric = metric, Site_1 = pair[1], Site_2 = pair[2],
      Pair = paste(pair, collapse = " vs "),
      N_site_1 = sum(group == pair[1]), N_site_2 = sum(group == pair[2]),
      PERMANOVA_R2 = ad$R2[1], PERMANOVA_pseudo_F = ad$F[1],
      PERMANOVA_p = ad$`Pr(>F)`[1],
      Betadisper_F = bp$tab$F[1], Betadisper_p = bp$tab$`Pr(>F)`[1],
      Permutations = nperm, stringsAsFactors = FALSE
    )
  }))
}))
pairwise_tests$PERMANOVA_p_FDR_BH <- ave(pairwise_tests$PERMANOVA_p,
                                         pairwise_tests$Metric,
                                         FUN = function(x) p.adjust(x, "BH"))
pairwise_tests$Betadisper_p_FDR_BH <- ave(pairwise_tests$Betadisper_p,
                                          pairwise_tests$Metric,
                                          FUN = function(x) p.adjust(x, "BH"))
pairwise_tests$Interpretation_flag <- ifelse(
  pairwise_tests$PERMANOVA_p_FDR_BH < 0.05 & pairwise_tests$Betadisper_p_FDR_BH < 0.05,
  "Location effect significant; dispersion also differs",
  ifelse(pairwise_tests$PERMANOVA_p_FDR_BH < 0.05,
         "Location effect significant; no FDR-significant dispersion difference",
         "No FDR-significant location effect")
)
write.csv(pairwise_tests, file.path(result_dir, "pairwise_permanova_betadisper.csv"),
          row.names = FALSE)

write.csv(data.frame(
  Parameter = c("Input_file", "Community", "Samples", "Original_ASVs", "Aggregated_genera", "Retained_genera",
                "Prevalence_threshold", "Prevalence_min_samples", "CLR_pseudocount",
                "PERMANOVA_permutations", "Betadisper_permutations", "vegan_version"),
  Value = c(input_file, community_label, nrow(otu), ncol(otu), ncol(genus_counts), ncol(filtered), "20%", prevalence_min,
            1, nperm, nperm, as.character(packageVersion("vegan")))
), file.path(result_dir, "analysis_parameters.csv"), row.names = FALSE)

palette_metric <- c("CLR-Aitchison" = "#8E1B61", "Robust Aitchison" = "#0072B2")
theme_pub <- theme_classic(base_size = 8, base_family = "sans") +
  theme(axis.line = element_line(linewidth = 0.35),
        axis.ticks = element_line(linewidth = 0.35),
        plot.title = element_text(size = 8.5, face = "bold"),
        legend.position = "none")

make_metric_row <- function(metric) {
  h <- heatmap_data[heatmap_data$Metric == metric, ]
  h$Site_1 <- factor(h$Site_1, levels = rev(site_levels))
  h$Site_2 <- factor(h$Site_2, levels = site_levels)
  p_a <- ggplot(h, aes(Site_2, Site_1, fill = Mean_distance)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.1f", Mean_distance)), size = 2.25) +
    scale_fill_viridis_c(option = "C", direction = 1) +
    coord_equal() +
    labs(title = paste0("a  Mean ", metric, " distance"), x = NULL, y = NULL,
         fill = "Mean\ndistance") +
    theme_minimal(base_size = 8) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1),
          plot.title = element_text(size = 8.5, face = "bold"), legend.position = "right",
          legend.key.height = grid::unit(8, "mm"))

  z <- pairwise_distances[pairwise_distances$Metric == metric, ]
  order_pairs <- distance_summary$Pair[distance_summary$Metric == metric][
    order(distance_summary$Mean_distance[distance_summary$Metric == metric])]
  z$Pair <- factor(z$Pair, levels = order_pairs)
  p_b <- ggplot(z, aes(Pair, Distance)) +
    geom_violin(fill = palette_metric[metric], colour = NA, alpha = 0.20,
                width = 0.85, trim = TRUE) +
    geom_boxplot(width = 0.17, outlier.shape = NA, fill = "white",
                 colour = palette_metric[metric], linewidth = 0.4) +
    geom_jitter(width = 0.12, height = 0, size = 0.42, alpha = 0.20,
                colour = palette_metric[metric]) +
    coord_flip() +
    labs(title = "b  Cross-locality sample-pair distances",
         subtitle = "Points are descriptive and non-independent", x = NULL, y = "Distance") +
    theme_pub

  t <- pairwise_tests[pairwise_tests$Metric == metric, ]
  t$Pair <- factor(t$Pair, levels = order_pairs)
  t$Dispersion <- ifelse(t$Betadisper_p_FDR_BH < 0.05, "Dispersion FDR < 0.05",
                         "Dispersion FDR >= 0.05")
  xmax <- max(t$PERMANOVA_R2) * 1.70
  p_c <- ggplot(t, aes(Pair, PERMANOVA_R2, fill = Dispersion, shape = Dispersion)) +
    geom_point(size = 2.8, stroke = 0.55) + coord_flip() +
    geom_text(aes(y = xmax * 0.98,
                  label = paste0("FDR ", formatC(PERMANOVA_p_FDR_BH,
                                                  format = "g", digits = 2))),
              hjust = 1, size = 1.85, show.legend = FALSE) +
    scale_fill_manual(values = c("Dispersion FDR < 0.05" = "#D55E00",
                                 "Dispersion FDR >= 0.05" = "white")) +
    scale_shape_manual(values = c("Dispersion FDR < 0.05" = 21,
                                  "Dispersion FDR >= 0.05" = 21)) +
    scale_y_continuous(limits = c(0, xmax),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(title = expression(bold("c  PERMANOVA " * R^2)), x = NULL,
         y = expression(PERMANOVA~R^2), fill = NULL, shape = NULL) +
    theme_pub
  p_a + p_b + p_c + plot_layout(widths = c(0.92, 1.22, 1.08))
}

fig <- make_metric_row("CLR-Aitchison") / make_metric_row("Robust Aitchison") +
  plot_annotation(
    title = paste0("Pairwise locality comparison of the ", community_label),
    subtitle = "Genus level; 20% prevalence filter; 47 independent frog samples.",
    caption = paste0("Diagonal cells: mean within-locality distance. Orange points: betadisper FDR < 0.05. ",
                     "PERMANOVA and betadisper: 9,999 permutations; BH correction within each metric."),
    theme = theme(plot.title = element_text(size = 11, face = "bold"),
                  plot.subtitle = element_text(size = 7.5),
                  plot.caption = element_text(size = 6.2, hjust = 0))
  )

save_plot <- function(plot, stem, width_mm = 183, height_mm = 190) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  ggsave(paste0(stem, ".png"), plot, width = w, height = h, dpi = 600, bg = "white")
  ragg::agg_tiff(paste0(stem, ".tiff"), width = w, height = h, units = "in",
                 res = 600, background = "white"); print(plot); dev.off()
  svglite::svglite(paste0(stem, ".svg"), width = w, height = h); print(plot); dev.off()
  grDevices::pdf(paste0(stem, ".pdf"), width = w, height = h,
                 family = "Helvetica", useDingbats = FALSE); print(plot); dev.off()
}
save_plot(fig, file.path(figure_dir, "figure_pairwise_aitchison_three_panel"))

print(distance_summary)
print(pairwise_tests)
