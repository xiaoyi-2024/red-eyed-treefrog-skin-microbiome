suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
})

whole_dir <- "results/pairwise_aitchison_comparison"
anti_dir <- "results/anti_bd_pairwise_aitchison_comparison"
result_dir <- "results/whole_vs_antibd_aitchison"
figure_dir <- file.path(result_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

whole_heat <- read.csv(file.path(whole_dir, "mean_distance_matrices_long.csv"),
                       stringsAsFactors = FALSE)
anti_heat <- read.csv(file.path(anti_dir, "mean_distance_matrices_long.csv"),
                      stringsAsFactors = FALSE)
whole_pair <- read.csv(file.path(whole_dir, "pairwise_distance_summary.csv"),
                       stringsAsFactors = FALSE)
anti_pair <- read.csv(file.path(anti_dir, "pairwise_distance_summary.csv"),
                      stringsAsFactors = FALSE)

whole_heat$Community <- "Whole community"
anti_heat$Community <- "Anti-Bd community"
heat <- rbind(whole_heat, anti_heat)

short <- c("Lost" = "Lost", "Veragua" = "Veragua", "Campana" = "Campana",
           "Soberanía" = "Soberanía")
site_levels <- names(short)

comparison <- merge(
  transform(whole_pair, Whole_mean_distance = Mean_distance)[,
    c("Metric", "Pair", "Site_1", "Site_2", "N_cross_sample_pairs", "Whole_mean_distance")],
  transform(anti_pair, Anti_Bd_mean_distance = Mean_distance)[,
    c("Metric", "Pair", "Anti_Bd_mean_distance")],
  by = c("Metric", "Pair"), sort = FALSE
)
comparison$Difference_AntiBd_minus_whole <-
  comparison$Anti_Bd_mean_distance - comparison$Whole_mean_distance
within_mean <- function(metric, community, site_1, site_2) {
  z <- heat[heat$Metric == metric & heat$Community == community, ]
  mean(c(z$Mean_distance[z$Site_1 == site_1 & z$Site_2 == site_1],
         z$Mean_distance[z$Site_1 == site_2 & z$Site_2 == site_2]))
}
comparison$Whole_reference_within_distance <- mapply(
  within_mean, metric = comparison$Metric, site_1 = comparison$Site_1,
  site_2 = comparison$Site_2, MoreArgs = list(community = "Whole community")
)
comparison$Anti_Bd_reference_within_distance <- mapply(
  within_mean, metric = comparison$Metric, site_1 = comparison$Site_1,
  site_2 = comparison$Site_2, MoreArgs = list(community = "Anti-Bd community")
)
comparison$Whole_between_within_ratio <- comparison$Whole_mean_distance /
  comparison$Whole_reference_within_distance
comparison$Anti_Bd_between_within_ratio <- comparison$Anti_Bd_mean_distance /
  comparison$Anti_Bd_reference_within_distance
comparison$Ratio_difference_AntiBd_minus_whole <-
  comparison$Anti_Bd_between_within_ratio - comparison$Whole_between_within_ratio
comparison$Inference_note <- paste0(
  "Descriptive cross-sample-pair means; pairs sharing a frog are non-independent; ",
  "distance scales are metric-specific")

write.csv(heat, file.path(result_dir, "whole_vs_antibd_mean_distance_matrices.csv"),
          row.names = FALSE)
write.csv(comparison, file.path(result_dir, "whole_vs_antibd_pairwise_distances.csv"),
          row.names = FALSE)

theme_pub <- theme_classic(base_size = 8, base_family = "sans") +
  theme(axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4),
        plot.title = element_text(size = 8.5, face = "bold"),
        plot.subtitle = element_text(size = 6.8))

heat_panel <- function(metric, community, tag, limits, breaks, show_legend) {
  z <- heat[heat$Metric == metric & heat$Community == community, ]
  z$Site_1 <- factor(z$Site_1, levels = rev(site_levels))
  z$Site_2 <- factor(z$Site_2, levels = site_levels)
  ggplot(z, aes(Site_2, Site_1, fill = Mean_distance)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.1f", Mean_distance)), size = 2.15) +
    scale_fill_stepsn(colors = viridisLite::viridis(8, option = "D"),
                      limits = limits, breaks = breaks, oob = scales::squish) +
    coord_equal() +
    labs(title = paste0(tag, "  ", community), x = NULL, y = NULL,
         fill = "Mean distance") +
    theme_minimal(base_size = 7.5) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 35, hjust = 1),
          plot.title = element_text(size = 8.5, face = "bold"),
          legend.position = if (show_legend) "bottom" else "none",
          legend.key.width = grid::unit(13, "mm"))
}

dumbbell_panel <- function(metric, tag, xmax) {
  z <- comparison[comparison$Metric == metric, ]
  pair_order <- z$Pair[order(z$Whole_mean_distance)]
  z$Pair <- factor(z$Pair, levels = pair_order)
  long <- rbind(
    data.frame(Pair = z$Pair, Community = "Whole community",
               Ratio = z$Whole_between_within_ratio),
    data.frame(Pair = z$Pair, Community = "Anti-Bd community",
               Ratio = z$Anti_Bd_between_within_ratio)
  )
  ggplot(z, aes(Pair)) +
    geom_segment(aes(y = Whole_between_within_ratio, yend = Anti_Bd_between_within_ratio,
                     xend = Pair), colour = "#7A7A7A", linewidth = 0.8,
                 alpha = 0.75) +
    geom_point(data = long, aes(y = Ratio, fill = Community, shape = Community),
               size = 2.6, stroke = 0.6) +
    coord_flip() +
    scale_fill_manual(values = c("Whole community" = "white",
                                 "Anti-Bd community" = "#0072B2")) +
    scale_shape_manual(values = c("Whole community" = 21,
                                  "Anti-Bd community" = 22)) +
    scale_y_continuous(limits = c(0, xmax), expand = expansion(mult = c(0, 0.03))) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey45", linewidth = 0.4) +
    labs(title = paste0(tag, "  Standardized separation"),
         subtitle = "Between-locality / mean within-locality distance",
         x = NULL, y = "Between/within distance ratio", fill = NULL, shape = NULL) +
    theme_pub + theme(legend.position = "top", legend.text = element_text(size = 6.3))
}

# Within a metric, both communities use the same scale. CLR and robust
# Aitchison distances are not placed on a shared numerical scale.
clr_limits <- c(10, 48)
robust_limits <- c(3, 21)

row_clr <- heat_panel("CLR-Aitchison", "Whole community", "a", clr_limits,
                      c(11, 13, 15, 20, 30, 40, 47), FALSE) |
  heat_panel("CLR-Aitchison", "Anti-Bd community", "b", clr_limits,
             c(11, 13, 15, 20, 30, 40, 47), TRUE) |
  dumbbell_panel("CLR-Aitchison", "c", 1.55)
row_robust <- heat_panel("Robust Aitchison", "Whole community", "a", robust_limits,
                         c(3.5, 5, 7, 10, 14, 18, 20.5), FALSE) |
  heat_panel("Robust Aitchison", "Anti-Bd community", "b", robust_limits,
             c(3.5, 5, 7, 10, 14, 18, 20.5), TRUE) |
  dumbbell_panel("Robust Aitchison", "c", 3.25)

figure_clr <- row_clr +
  plot_layout(widths = c(0.92, 0.92, 1.25)) +
  plot_annotation(
    title = "Whole-community and anti-Bd CLR-Aitchison distances",
    subtitle = "Genus level; 20% prevalence filtering within each community; pseudocount = 1.",
    caption = paste0("Heat maps: raw metric-specific distances. Right: between-locality / mean within-locality distance; ",
                     "dashed line = 1. Sample pairs are non-independent."),
    theme = theme(plot.title = element_text(size = 11, face = "bold"),
                  plot.subtitle = element_text(size = 7.5),
                  plot.caption = element_text(size = 6.2, hjust = 0))
  )

figure_robust <- row_robust +
  plot_layout(widths = c(0.92, 0.92, 1.25)) +
  plot_annotation(
    title = "Whole-community and anti-Bd robust Aitchison distances",
    subtitle = "Genus level; 20% prevalence filtering within each community; rCLR transformation.",
    caption = paste0("Heat maps: raw metric-specific distances. Right: between-locality / mean within-locality distance; ",
                     "dashed line = 1. Sample pairs are non-independent."),
    theme = theme(plot.title = element_text(size = 11, face = "bold"),
                  plot.subtitle = element_text(size = 7.5),
                  plot.caption = element_text(size = 6.2, hjust = 0))
  )

save_plot <- function(plot, stem, width_mm = 183, height_mm = 100) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  ggsave(paste0(stem, ".png"), plot, width = w, height = h, dpi = 600, bg = "white")
  ragg::agg_tiff(paste0(stem, ".tiff"), width = w, height = h, units = "in",
                 res = 600, background = "white"); print(plot); dev.off()
  svglite::svglite(paste0(stem, ".svg"), width = w, height = h); print(plot); dev.off()
  grDevices::pdf(paste0(stem, ".pdf"), width = w, height = h,
                 family = "Helvetica", useDingbats = FALSE); print(plot); dev.off()
}
save_plot(figure_clr, file.path(figure_dir, "figure_whole_vs_antibd_clr_aitchison"))
save_plot(figure_robust, file.path(figure_dir, "figure_whole_vs_antibd_robust_aitchison"))
print(comparison)
