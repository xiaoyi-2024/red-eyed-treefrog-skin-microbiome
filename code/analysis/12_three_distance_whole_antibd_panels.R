#!/usr/bin/env Rscript

# Publication panels comparing whole-community (top row) and anti-Bd
# community (bottom row) patterns for Bray-Curtis, CLR-Aitchison, and robust
# Aitchison distances. Sample-pair distances are descriptive/non-independent;
# inference is based on PERMANOVA with betadisper reported alongside it.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
})

out_dir <- "results/whole_vs_antibd_three_distances"
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

sites <- c("Lost", "Veragua", "Campana", "Soberanía")
short <- c("Lost Iguana" = "Lost", "Veragua" = "Veragua",
           "Altos de Campana" = "Campana", "Soberanía" = "Soberanía")

read_bray <- function(dir, community) {
  pair <- read.csv(file.path(dir, "pairwise_locality_bray_permanova.csv"),
                   check.names = FALSE, stringsAsFactors = FALSE)
  within <- read.csv(file.path(dir, "within_locality_bray_distances.csv"),
                     check.names = FALSE, stringsAsFactors = FALSE)
  raw <- read.csv(file.path(dir, "descriptive_sample_pair_distances.csv"),
                  check.names = FALSE, stringsAsFactors = FALSE)
  heat <- expand.grid(Site_1 = sites, Site_2 = sites, stringsAsFactors = FALSE)
  heat$Mean_distance <- NA_real_
  for (i in seq_len(nrow(heat))) {
    a <- heat$Site_1[i]; b <- heat$Site_2[i]
    if (a == b) {
      heat$Mean_distance[i] <- within$Mean_within_Bray_Curtis_distance[
        match(a, unname(short[within$Site]))]
    } else {
      hit <- pair[(short[pair$Site_1] == a & short[pair$Site_2] == b) |
                  (short[pair$Site_1] == b & short[pair$Site_2] == a), ]
      heat$Mean_distance[i] <- hit$Mean_Bray_Curtis_distance
    }
  }
  raw <- raw[raw$Site_1 != raw$Site_2, ]
  raw$Pair <- vapply(strsplit(raw$Pair, " vs ", fixed = TRUE), function(z)
    paste(unname(short[z]), collapse = " vs "), character(1))
  tests <- transform(pair,
    Site_1 = unname(short[Site_1]), Site_2 = unname(short[Site_2]),
    Pair = paste(unname(short[Site_1]), unname(short[Site_2]), sep = " vs "),
    PERMANOVA_p_FDR_BH = PERMANOVA_FDR,
    Betadisper_p_FDR_BH = Betadisper_FDR)
  list(heat = heat, raw = transform(raw, Distance = Bray_Curtis),
       tests = tests, community = community)
}

read_aitch <- function(dir, metric, community) {
  heat <- read.csv(file.path(dir, "mean_distance_matrices_long.csv"),
                   check.names = FALSE, stringsAsFactors = FALSE)
  raw <- read.csv(file.path(dir, "cross_locality_sample_pair_distances.csv"),
                  check.names = FALSE, stringsAsFactors = FALSE)
  tests <- read.csv(file.path(dir, "pairwise_permanova_betadisper.csv"),
                    check.names = FALSE, stringsAsFactors = FALSE)
  list(heat = heat[heat$Metric == metric, ],
       raw = raw[raw$Metric == metric, ],
       tests = tests[tests$Metric == metric, ], community = community)
}

whole_bray <- read_bray("results/pairwise_bray_distance", "Whole community")
anti_bray <- read_bray("results/anti_bd_pairwise_bray_distance", "Anti-Bd community")
whole_clr <- read_aitch("results/pairwise_aitchison_comparison",
                        "CLR-Aitchison", "Whole community")
anti_clr <- read_aitch("results/anti_bd_pairwise_aitchison_comparison",
                       "CLR-Aitchison", "Anti-Bd community")
whole_robust <- read_aitch("results/pairwise_aitchison_comparison",
                           "Robust Aitchison", "Whole community")
anti_robust <- read_aitch("results/anti_bd_pairwise_aitchison_comparison",
                          "Robust Aitchison", "Anti-Bd community")

pretty_fdr <- function(x) {
  ifelse(x < 0.001, paste0("FDR ", formatC(x, format = "fg", digits = 2)),
         paste0("FDR ", formatC(x, format = "f", digits = 3)))
}

make_row <- function(dat, metric_label, tags, palette, distribution_colour,
                     common_limits = NULL) {
  h <- dat$heat
  h$Site_1 <- factor(h$Site_1, levels = rev(sites))
  h$Site_2 <- factor(h$Site_2, levels = sites)
  limits <- if (is.null(common_limits)) range(h$Mean_distance, finite = TRUE) else common_limits
  p_heat <- ggplot(h, aes(Site_2, Site_1, fill = Mean_distance)) +
    geom_tile(colour = "white", linewidth = 0.7) +
    geom_text(aes(label = ifelse(Mean_distance < 1,
                                 sprintf("%.3f", Mean_distance),
                                 sprintf("%.1f", Mean_distance))), size = 2.05) +
    scale_fill_gradientn(colours = palette, limits = limits,
                         oob = scales::squish, name = "Mean distance",
                         guide = guide_colorbar(direction = "horizontal",
                                                title.position = "top",
                                                barwidth = grid::unit(23, "mm"),
                                                barheight = grid::unit(2.2, "mm"))) +
    labs(title = paste0(tags[1], "  ", dat$community),
         subtitle = paste0("Mean ", metric_label),
         x = NULL, y = NULL) +
    theme_minimal(base_size = 7, base_family = "sans") +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 35, hjust = 1, size = 5.8),
          axis.text.y = element_text(size = 5.8),
          plot.title = element_text(face = "bold", size = 7.5),
          plot.subtitle = element_text(size = 5.8),
          legend.position = "bottom", legend.justification = "center",
          legend.title = element_text(size = 6), legend.text = element_text(size = 5.5),
          legend.box.margin = margin(t = -2, r = 0, b = 0, l = 0))

  order_pairs <- dat$tests$Pair[order(dat$tests$PERMANOVA_R2, decreasing = TRUE)]
  raw <- dat$raw
  raw$Pair <- factor(raw$Pair, levels = rev(order_pairs))
  p_dist <- ggplot(raw, aes(Distance, Pair)) +
    geom_violin(fill = scales::alpha(distribution_colour, 0.24),
                colour = distribution_colour, linewidth = 0.35,
                scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white", linewidth = 0.3) +
    geom_point(position = position_jitter(height = 0.08, width = 0),
               size = 0.35, alpha = 0.34, colour = distribution_colour) +
    labs(title = paste0(tags[2], "  Cross-locality sample-pair distances"),
         subtitle = "Descriptive, non-independent sample pairs",
         x = "Distance", y = NULL) +
    theme_classic(base_size = 7, base_family = "sans") +
    theme(axis.text = element_text(size = 5.5),
          plot.title = element_text(face = "bold", size = 7.5),
          plot.subtitle = element_text(size = 5.8))

  tt <- dat$tests
  tt$Pair <- factor(tt$Pair, levels = rev(order_pairs))
  tt$Dispersion <- tt$Betadisper_p_FDR_BH < 0.05
  tt$FDR_label <- pretty_fdr(tt$PERMANOVA_p_FDR_BH)
  xmax <- max(tt$PERMANOVA_R2) * 1.55
  p_r2 <- ggplot(tt, aes(PERMANOVA_R2, Pair)) +
    geom_segment(aes(x = 0, xend = PERMANOVA_R2, yend = Pair),
                 colour = "grey78", linewidth = 0.35) +
    geom_point(aes(fill = Dispersion), shape = 21, size = 2,
               colour = "black", stroke = 0.35) +
    geom_text(aes(x = xmax * 0.98, label = FDR_label), hjust = 1,
              size = 1.8, colour = "grey20") +
    scale_fill_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "white"),
                      guide = "none") +
    scale_x_continuous(limits = c(0, xmax), expand = expansion(mult = c(0, 0))) +
    labs(title = paste0(tags[3], "  PERMANOVA effect size"),
         x = expression(PERMANOVA~R^2), y = NULL) +
    theme_classic(base_size = 7, base_family = "sans") +
    theme(axis.text = element_text(size = 5.3),
          plot.title = element_text(face = "bold", size = 7.5),
          axis.title.x = element_text(size = 6.5))
  p_heat | p_dist | p_r2
}

save_pub <- function(plot, stem) {
  w <- 183 / 25.4; h <- 190 / 25.4
  ggsave(paste0(stem, ".png"), plot, width = w, height = h, dpi = 600, bg = "white")
  svglite::svglite(paste0(stem, ".svg"), width = w, height = h); print(plot); dev.off()
  grDevices::pdf(paste0(stem, ".pdf"), width = w, height = h,
                 family = "Helvetica", useDingbats = FALSE); print(plot); dev.off()
  ragg::agg_tiff(paste0(stem, ".tiff"), width = w, height = h, units = "in",
                 res = 600, compression = "lzw", background = "white")
  print(plot); dev.off()
}

make_figure <- function(whole, anti, metric_label, palette, distribution_colour,
                        stem, common_limits = NULL) {
  top <- make_row(whole, metric_label, c("a", "b", "c"), palette,
                  distribution_colour, common_limits)
  bottom <- make_row(anti, metric_label, c("d", "e", "f"), palette,
                     distribution_colour, common_limits)
  fig <- (top / bottom) +
    plot_layout(widths = c(1.0, 1.25, 0.9), heights = c(1, 1), guides = "keep") +
    plot_annotation(
      title = paste0(metric_label, ": whole community and anti-Bd community"),
      subtitle = "Genus level; 20% prevalence filter; 47 independent frog samples",
      caption = paste0("Top: whole community; bottom: anti-Bd community. Diagonal cells show mean within-locality distances.\n",
                       "Orange points: betadisper FDR < 0.05; PERMANOVA and betadisper use 9,999 permutations; BH correction within metric."),
      theme = theme(plot.title = element_text(face = "bold", size = 10),
                    plot.subtitle = element_text(size = 6.8),
                    plot.caption = element_text(size = 5.6, hjust = 0)))
  save_pub(fig, file.path(fig_dir, stem))
}

# Bray is bounded on the same 0-1 scale for both subcommunities.
make_figure(whole_bray, anti_bray, "Bray-Curtis distance",
            c("#FFF200", "#F59E42", "#C43C8C", "#5A00A5", "#160078"),
            "#7A0177",
            "figure_bray_whole_above_antibd", c(0.4, 0.8))
# Absolute Aitchison scales depend on the retained feature set, so each row has
# its own scale even though the layouts and inferential procedures are matched.
make_figure(whole_clr, anti_clr, "CLR-Aitchison distance",
            c("#3B0F70", "#B73779", "#F1605D", "#FBB32B", "#FDE725"),
            "#A11A67",
            "figure_clr_aitchison_whole_above_antibd")
make_figure(whole_robust, anti_robust, "Robust Aitchison distance",
            c("#08306B", "#2171B5", "#6BAED6", "#BDD7E7", "#F7FBFF"),
            "#B23A48",
            "figure_robust_aitchison_whole_above_antibd")

write.csv(data.frame(
  Figure = c("Bray-Curtis", "CLR-Aitchison", "Robust Aitchison"),
  Whole_row = "top", Anti_Bd_row = "bottom",
  Cross_community_absolute_scale_comparable = c(TRUE, FALSE, FALSE),
  stringsAsFactors = FALSE),
  file.path(out_dir, "three_distance_figure_notes.csv"), row.names = FALSE)
