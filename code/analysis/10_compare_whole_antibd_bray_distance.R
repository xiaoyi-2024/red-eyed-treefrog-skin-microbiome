suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
  library(scales)
})

whole_pair_file <- "results/pairwise_bray_distance/pairwise_locality_bray_permanova.csv"
whole_within_file <- "results/pairwise_bray_distance/within_locality_bray_distances.csv"
anti_pair_file <- "results/anti_bd_pairwise_bray_distance/pairwise_locality_bray_permanova.csv"
anti_within_file <- "results/anti_bd_pairwise_bray_distance/within_locality_bray_distances.csv"
result_dir <- "results/whole_vs_antibd_bray_distance"
figure_dir <- file.path(result_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

whole_pair <- read.csv(whole_pair_file, check.names = FALSE, stringsAsFactors = FALSE)
whole_within <- read.csv(whole_within_file, check.names = FALSE, stringsAsFactors = FALSE)
anti_pair <- read.csv(anti_pair_file, check.names = FALSE, stringsAsFactors = FALSE)
anti_within <- read.csv(anti_within_file, check.names = FALSE, stringsAsFactors = FALSE)

short <- c("Lost Iguana" = "Lost", "Veragua" = "Veragua",
           "Altos de Campana" = "Campana", "Soberanía" = "Soberanía")
site_levels <- unname(short)

make_heat <- function(pair, within, community) {
  out <- expand.grid(Site_1 = site_levels, Site_2 = site_levels,
                     stringsAsFactors = FALSE)
  out$Distance <- NA_real_
  for (i in seq_len(nrow(out))) {
    a <- out$Site_1[i]; b <- out$Site_2[i]
    if (a == b) {
      z <- within[short[within$Site] == a, ]
      out$Distance[i] <- z$Mean_within_Bray_Curtis_distance
    } else {
      z <- pair[(short[pair$Site_1] == a & short[pair$Site_2] == b) |
                (short[pair$Site_1] == b & short[pair$Site_2] == a), ]
      out$Distance[i] <- z$Mean_Bray_Curtis_distance
    }
  }
  out$Community <- community
  out
}

heat <- rbind(make_heat(whole_pair, whole_within, "Whole community"),
              make_heat(anti_pair, anti_within, "Anti-Bd community"))
write.csv(heat, file.path(result_dir, "whole_vs_antibd_bray_distance_matrices.csv"),
          row.names = FALSE)

comparison <- merge(
  transform(whole_pair,
            Pair_short = paste(short[Site_1], short[Site_2], sep = " vs "),
            Whole_mean_distance = Mean_Bray_Curtis_distance)[,
              c("Pair_short", "Cross_sample_pairs", "Whole_mean_distance")],
  transform(anti_pair,
            Pair_short = paste(short[Site_1], short[Site_2], sep = " vs "),
            Anti_Bd_mean_distance = Mean_Bray_Curtis_distance)[,
              c("Pair_short", "Anti_Bd_mean_distance")],
  by = "Pair_short", sort = FALSE
)
comparison$Difference_AntiBd_minus_whole <-
  comparison$Anti_Bd_mean_distance - comparison$Whole_mean_distance
comparison$Inference_note <- paste0(
  "Descriptive cross-sample-pair means; pairs sharing a frog are non-independent; ",
  "no inferential test is performed on the difference")
write.csv(comparison, file.path(result_dir, "whole_vs_antibd_pairwise_bray_distances.csv"),
          row.names = FALSE)

theme_pub <- theme_classic(base_size = 8, base_family = "sans") +
  theme(axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4),
        plot.title = element_text(size = 9, face = "bold"),
        plot.subtitle = element_text(size = 7))

heat_panel <- function(data, title, tag, show_legend = TRUE) {
  data$Site_1 <- factor(data$Site_1, levels = rev(site_levels))
  data$Site_2 <- factor(data$Site_2, levels = site_levels)
  ggplot(data, aes(Site_2, Site_1, fill = Distance)) +
    geom_tile(colour = "white", linewidth = 0.7) +
    geom_text(aes(label = sprintf("%.3f", Distance)), size = 2.35) +
    scale_fill_stepsn(colors = viridisLite::viridis(8, option = "D", direction = -1),
                      limits = c(0.4, 0.8),
                      breaks = seq(0.4, 0.8, 0.1),
                      oob = scales::squish) +
    coord_equal() +
    labs(title = paste0(tag, "  ", title), x = NULL, y = NULL,
         fill = "Mean distance") +
    theme_minimal(base_size = 8) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 35, hjust = 1),
          plot.title = element_text(size = 9, face = "bold"),
          legend.position = if (show_legend) "bottom" else "none",
          legend.key.width = grid::unit(9, "mm"))
}

p_a <- heat_panel(heat[heat$Community == "Whole community", ],
                  "Whole bacterial community", "a", FALSE)
p_b <- heat_panel(heat[heat$Community == "Anti-Bd community", ],
                  "Anti-Bd bacterial community", "b")

comparison$Pair_short <- factor(
  comparison$Pair_short,
  levels = comparison$Pair_short[order(comparison$Whole_mean_distance)]
)
long <- rbind(
  data.frame(Pair_short = comparison$Pair_short, Community = "Whole community",
             Distance = comparison$Whole_mean_distance),
  data.frame(Pair_short = comparison$Pair_short, Community = "Anti-Bd community",
             Distance = comparison$Anti_Bd_mean_distance)
)
comparison$Direction <- ifelse(comparison$Difference_AntiBd_minus_whole >= 0,
                               "Anti-Bd farther", "Anti-Bd closer")
p_c <- ggplot(comparison, aes(Pair_short)) +
  geom_segment(aes(y = Whole_mean_distance, yend = Anti_Bd_mean_distance,
                   xend = Pair_short, colour = Direction), linewidth = 1.0,
               alpha = 0.65) +
  geom_point(data = long,
             aes(y = Distance, fill = Community, shape = Community),
             size = 3.0, stroke = 0.65) +
  coord_flip() +
  scale_colour_manual(values = c("Anti-Bd farther" = "#D55E00",
                                 "Anti-Bd closer" = "#009E73"), guide = "none") +
  scale_fill_manual(values = c("Whole community" = "white",
                               "Anti-Bd community" = "#0072B2")) +
  scale_shape_manual(values = c("Whole community" = 21, "Anti-Bd community" = 22)) +
  scale_y_continuous(limits = c(0, 0.85), breaks = seq(0, 0.8, 0.2),
                     expand = expansion(mult = c(0, 0.03))) +
  labs(title = "c  Direct comparison of between-locality distance",
       subtitle = "Lines connect estimates for the same locality pair",
       x = NULL, y = "Mean Bray-Curtis distance",
       fill = NULL, shape = NULL, colour = NULL) +
  theme_pub + theme(legend.position = "top")

figure <- (p_a | p_b) / p_c +
  plot_layout(heights = c(1.0, 1.05)) +
  plot_annotation(
    title = "Whole-community and anti-Bd Bray-Curtis distances",
    subtitle = "Genus relative abundance; 20% prevalence filter applied separately to each community.",
    caption = paste0("Heat-map diagonals are within-locality means; off-diagonals are between-locality means.\n",
                     "Orange line: anti-Bd distance larger; green: smaller. Cross-sample pairs are descriptive and non-independent."),
    theme = theme(plot.title = element_text(size = 11, face = "bold"),
                  plot.subtitle = element_text(size = 7.5),
                  plot.caption = element_text(size = 6.3, hjust = 0))
  )

save_plot <- function(plot, stem, width_mm = 183, height_mm = 180) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  ggsave(paste0(stem, ".png"), plot, width = w, height = h, dpi = 600, bg = "white")
  ragg::agg_tiff(paste0(stem, ".tiff"), width = w, height = h, units = "in",
                 res = 600, background = "white"); print(plot); dev.off()
  svglite::svglite(paste0(stem, ".svg"), width = w, height = h); print(plot); dev.off()
  grDevices::pdf(paste0(stem, ".pdf"), width = w, height = h,
                 family = "Helvetica", useDingbats = FALSE); print(plot); dev.off()
}
save_plot(figure, file.path(figure_dir, "figure_whole_vs_antibd_bray_distance"))
print(comparison)
