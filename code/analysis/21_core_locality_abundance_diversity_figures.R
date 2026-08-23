#!/usr/bin/env Rscript

# Publication figures that complete the core analysis narrative:
# 1) overall locality differentiation and within-locality heterogeneity;
# 2) anti-Bd relative abundance plus whole/anti-Bd Hill D0-D2 profiles;
# 3) whole/anti-Bd community composition along estimated Bd-arrival year.

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(iNEXT)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(svglite)
  library(ragg)
})

set.seed(20260817)

result_dir <- "results/core_locality_abundance_diversity"
figure_dir <- file.path(result_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

full_file <- "data/processed/ps_clean_47_rebuilt_from_01_filtered.rds"
anti_file <- "data/processed/blast_97id_90cov/ps_anti_bd_47_rebuilt_97id_90cov.rds"
wave_file <- "results/wave_exploratory/site_metrics_single_wave.csv"
bb_file <- "results/wave_exploratory/beta_binomial_site_estimates.csv"
stopifnot(file.exists(full_file), file.exists(anti_file),
          file.exists(wave_file), file.exists(bb_file))

site_key <- c(
  "Lost Iguana Resort, La Fortuna, San Carlos" = "Lost Iguana",
  "Veragua Rainforest, Las Brisas de Veragua" = "Veragua",
  "PN Altos de Campana" = "Altos de Campana",
  "PN Soberanía, Las Cumbres" = "Soberanía"
)
site_levels <- c("Lost Iguana", "Veragua", "Altos de Campana", "Soberanía")
site_palette <- c("Lost Iguana" = "#B44682", "Veragua" = "#216DA3",
                  "Altos de Campana" = "#21865F", "Soberanía" = "#D17A00")
community_palette <- c("Whole community" = "#365E8D", "Anti-Bd community" = "#B54B88")

otu_samples_rows <- function(ps) {
  x <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) x <- t(x)
  x
}

aggregate_genus <- function(ps, counts) {
  tx <- as(tax_table(ps), "matrix")[colnames(counts), , drop = FALSE]
  genus <- as.character(tx[, "Genus"])
  family <- as.character(tx[, "Family"])
  family[is.na(family) | trimws(family) == ""] <- "Bacteria"
  miss <- is.na(genus) | trimws(genus) == ""
  genus[miss] <- paste0("Unclassified_", family[miss])
  t(rowsum(t(counts), group = genus, reorder = FALSE))
}

ps_full <- readRDS(full_file)
ps_anti <- readRDS(anti_file)
full_asv <- otu_samples_rows(ps_full)
anti_asv <- otu_samples_rows(ps_anti)[rownames(full_asv), , drop = FALSE]
md <- data.frame(sample_data(ps_full), check.names = FALSE,
                 stringsAsFactors = FALSE)[rownames(full_asv), , drop = FALSE]
site <- factor(unname(site_key[md$Locality]), levels = site_levels)
stopifnot(nrow(full_asv) == 47L, !anyNA(site), all(rowSums(anti_asv) > 0))

full_genus <- aggregate_genus(ps_full, full_asv)
anti_genus <- aggregate_genus(ps_anti, anti_asv)

# Prevalence filtering is used for beta-diversity, while Hill numbers retain
# all observed genera and standardize sample completeness by coverage.
prevalence_filter <- function(x, fraction = 0.20) {
  keep <- colSums(x > 0) >= ceiling(nrow(x) * fraction)
  x[, keep, drop = FALSE]
}

full_beta <- prevalence_filter(full_genus)
anti_beta <- prevalence_filter(anti_genus)
full_rclr <- vegan::decostand(full_beta, method = "rclr")
anti_rclr <- vegan::decostand(anti_beta, method = "rclr")
stopifnot(all(is.finite(full_rclr)), all(is.finite(anti_rclr)))

wave <- read.csv(wave_file, stringsAsFactors = FALSE, check.names = FALSE)
site_year <- setNames(wave$pred_wave_year, wave$Site)
year <- unname(site_year[as.character(site)])
stopifnot(!anyNA(year))

permutations_all <- function(x) {
  if (length(x) == 1L) return(matrix(x, nrow = 1L))
  do.call(rbind, lapply(seq_along(x), function(i) {
    cbind(x[i], permutations_all(x[-i]))
  }))
}

theme_pub <- theme_classic(base_size = 8, base_family = "sans") +
  theme(
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.4),
    plot.title = element_text(size = 9, face = "bold"),
    plot.subtitle = element_text(size = 6.8),
    strip.text = element_text(size = 7.5, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 6.5),
    plot.caption = element_text(size = 5.8, hjust = 0)
  )

save_pub <- function(plot, stem, width_mm = 183, height_mm = 145, dpi = 600) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  ggsave(paste0(stem, ".png"), plot, width = w, height = h, dpi = dpi,
         bg = "white")
  svglite::svglite(paste0(stem, ".svg"), width = w, height = h)
  print(plot); dev.off()
  grDevices::pdf(paste0(stem, ".pdf"), width = w, height = h,
                 family = "Helvetica", useDingbats = FALSE)
  print(plot); dev.off()
  ragg::agg_tiff(paste0(stem, ".tiff"), width = w, height = h,
                 units = "in", res = dpi, background = "white")
  print(plot); dev.off()
}

community_objects <- list(
  "Whole community" = full_rclr,
  "Anti-Bd community" = anti_rclr
)

# -------------------------------------------------------------------------
# Figure 1: overall locality differentiation + within-locality heterogeneity
# -------------------------------------------------------------------------

locality_tests <- do.call(rbind, lapply(names(community_objects), function(nm) {
  coord <- community_objects[[nm]]
  D <- dist(coord)
  ad <- vegan::adonis2(D ~ site, permutations = 9999)
  bd <- vegan::betadisper(D, site, type = "centroid", bias.adjust = TRUE)
  pt <- vegan::permutest(bd, permutations = 9999)
  data.frame(
    Community = nm, Sample_n = nrow(coord), Genera_n = ncol(coord),
    PERMANOVA_R2 = unname(ad$R2[1]), PERMANOVA_pseudo_F = unname(ad$F[1]),
    PERMANOVA_p = unname(ad$`Pr(>F)`[1]),
    Betadisper_F = unname(pt$tab[1, "F"]),
    Betadisper_p = unname(pt$tab[1, "Pr(>F)"]),
    Permutations = 9999L, stringsAsFactors = FALSE
  )
}))
locality_tests$PERMANOVA_Holm_p <- p.adjust(locality_tests$PERMANOVA_p, "holm")
locality_tests$Betadisper_Holm_p <- p.adjust(locality_tests$Betadisper_p, "holm")
write.csv(locality_tests, file.path(result_dir, "overall_locality_tests.csv"),
          row.names = FALSE)

ordination_panel <- function(coord, nm, tag) {
  pc <- prcomp(coord, center = TRUE, scale. = FALSE)
  pct <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  z <- data.frame(Sample_ID = rownames(coord), PC1 = pc$x[, 1], PC2 = pc$x[, 2],
                  Site = site)
  hit <- locality_tests[locality_tests$Community == nm, ]
  ggplot(z, aes(PC1, PC2, colour = Site)) +
    stat_ellipse(aes(group = Site), type = "norm", level = 0.68,
                 linewidth = 0.45, alpha = 0.75, show.legend = FALSE) +
    geom_point(size = 1.8, alpha = 0.72) +
    scale_colour_manual(values = site_palette) +
    labs(
      title = paste0(tag, "  ", nm, " composition"),
      subtitle = sprintf("Robust Aitchison PERMANOVA: R² = %.3f; p = %.4f",
                         hit$PERMANOVA_R2, hit$PERMANOVA_p),
      x = sprintf("rCLR-PC1 (%.1f%%)", pct[1]),
      y = sprintf("rCLR-PC2 (%.1f%%)", pct[2])
    ) + theme_pub + theme(legend.position = "none")
}

dispersion_values <- do.call(rbind, lapply(names(community_objects), function(nm) {
  bd <- vegan::betadisper(dist(community_objects[[nm]]), site,
                          type = "centroid", bias.adjust = TRUE)
  data.frame(Sample_ID = names(bd$distances), Site = site,
             Distance = unname(bd$distances), Community = nm)
}))
write.csv(dispersion_values,
          file.path(result_dir, "sample_robust_aitchison_centroid_distances.csv"),
          row.names = FALSE)

dispersion_panel <- function(nm, tag) {
  z <- dispersion_values[dispersion_values$Community == nm, ]
  hit <- locality_tests[locality_tests$Community == nm, ]
  ggplot(z, aes(Site, Distance, fill = Site, colour = Site)) +
    geom_violin(width = 0.78, alpha = 0.18, linewidth = 0.35,
                trim = FALSE, show.legend = FALSE) +
    geom_boxplot(width = 0.28, alpha = 0.60, outlier.shape = NA,
                 linewidth = 0.4, show.legend = FALSE) +
    geom_point(position = position_jitter(width = 0.09, height = 0),
               size = 1.25, alpha = 0.72, show.legend = FALSE) +
    scale_fill_manual(values = site_palette) +
    scale_colour_manual(values = site_palette) +
    labs(
      title = paste0(tag, "  ", nm, " heterogeneity"),
      subtitle = sprintf("PERMDISP: F = %.2f; p = %.3f",
                         hit$Betadisper_F, hit$Betadisper_p),
      x = NULL, y = "Distance to locality centroid"
    ) + theme_pub +
    theme(axis.text.x = element_text(angle = 24, hjust = 1),
          legend.position = "none")
}

p_locality <- (
  ordination_panel(full_rclr, "Whole community", "a") |
    ordination_panel(anti_rclr, "Anti-Bd community", "b")
) / (
  dispersion_panel("Whole community", "c") |
    dispersion_panel("Anti-Bd community", "d")
) +
  plot_annotation(
    title = "Locality differentiation and within-locality heterogeneity",
    subtitle = "Whole and anti-Bd genus communities; 20% prevalence filter; robust Aitchison geometry.",
    caption = paste0(
      "Points are individual frogs (n = 47). Ellipses show 68% normal-data regions for visualization.\n",
      "PERMANOVA and PERMDISP use 9,999 permutations; boxes show medians and interquartile ranges."
    ),
    theme = theme(plot.title = element_text(size = 11, face = "bold"),
                  plot.subtitle = element_text(size = 7.2),
                  plot.caption = element_text(size = 6, hjust = 0))
  )
save_pub(p_locality, file.path(figure_dir, "figure_overall_locality_composition_heterogeneity"))

# -------------------------------------------------------------------------
# Figure 2: anti-Bd relative abundance + Hill D0/D1/D2 diversity spectra
# -------------------------------------------------------------------------

hill_by_coverage <- function(x, community) {
  assemblages <- lapply(seq_len(nrow(x)), function(i) {
    z <- x[i, ]
    unname(z[z > 0])
  })
  names(assemblages) <- rownames(x)
  info <- iNEXT::DataInfo(assemblages, datatype = "abundance")
  coverage <- min(info$SC)
  est <- iNEXT::estimateD(assemblages, datatype = "abundance", base = "coverage",
                          level = coverage, q = c(0, 1, 2), nboot = 0)
  out <- data.frame(
    Sample_ID = est$Assemblage,
    Order_q = est$Order.q,
    Hill_D = est$qD,
    Community = community,
    Common_coverage = coverage,
    stringsAsFactors = FALSE
  )
  out$Site <- factor(as.character(site[match(out$Sample_ID, rownames(x))]),
                     levels = site_levels)
  list(values = out, coverage = info)
}

hill_full <- hill_by_coverage(full_genus, "Whole community")
hill_anti <- hill_by_coverage(anti_genus, "Anti-Bd community")
hill_values <- rbind(hill_full$values, hill_anti$values)
write.csv(hill_values, file.path(result_dir, "sample_hill_D0_D1_D2.csv"),
          row.names = FALSE)
write.csv(rbind(transform(hill_full$coverage, Community = "Whole community"),
                transform(hill_anti$coverage, Community = "Anti-Bd community")),
          file.path(result_dir, "hill_sample_coverage.csv"), row.names = FALSE)

bootstrap_hill <- function(z, B = 4999L) {
  set.seed(20260817 + z$Order_q[1] + ifelse(z$Community[1] == "Anti-Bd community", 10, 0))
  b <- replicate(B, mean(sample(z$Hill_D, nrow(z), replace = TRUE)))
  data.frame(
    Community = z$Community[1], Site = z$Site[1], Order_q = z$Order_q[1],
    Sample_n = nrow(z), Mean = mean(z$Hill_D), SE = sd(z$Hill_D) / sqrt(nrow(z)),
    CI_low = unname(quantile(b, 0.025)), CI_high = unname(quantile(b, 0.975)),
    Common_coverage = z$Common_coverage[1]
  )
}
hill_summary <- do.call(rbind, lapply(
  split(hill_values, interaction(hill_values$Community, hill_values$Site,
                                 hill_values$Order_q, drop = TRUE)),
  bootstrap_hill
))
hill_summary$Site <- factor(as.character(hill_summary$Site), levels = site_levels)
hill_summary$Order_label <- factor(
  paste0("D", hill_summary$Order_q), levels = c("D0", "D1", "D2")
)
write.csv(hill_summary, file.path(result_dir, "locality_hill_D0_D1_D2_summary.csv"),
          row.names = FALSE)

sample_abundance <- data.frame(
  Sample_ID = rownames(full_asv), Site = site,
  Anti_Bd_percent = 100 * rowSums(anti_asv) / rowSums(full_asv)
)
bb <- read.csv(bb_file, stringsAsFactors = FALSE, check.names = FALSE)
bb$Site <- factor(bb$Site, levels = site_levels)
bb$CI_low <- pmax(0, bb$beta_binomial_mean_percent - 1.96 * bb$beta_binomial_mean_SE_percent)
bb$CI_high <- bb$beta_binomial_mean_percent + 1.96 * bb$beta_binomial_mean_SE_percent
write.csv(sample_abundance, file.path(result_dir, "sample_anti_bd_relative_abundance.csv"),
          row.names = FALSE)

p_abundance <- ggplot(sample_abundance,
                      aes(Site, Anti_Bd_percent, fill = Site, colour = Site)) +
  geom_violin(width = 0.80, alpha = 0.17, linewidth = 0.35,
              trim = FALSE, show.legend = FALSE) +
  geom_boxplot(width = 0.26, alpha = 0.48, outlier.shape = NA,
               linewidth = 0.4, show.legend = FALSE) +
  geom_point(position = position_jitter(width = 0.09, height = 0),
             size = 1.2, alpha = 0.65, show.legend = FALSE) +
  geom_errorbar(data = bb, aes(Site, beta_binomial_mean_percent,
                               ymin = CI_low, ymax = CI_high),
                inherit.aes = FALSE, width = 0.12, linewidth = 0.7,
                colour = "black") +
  geom_point(data = bb, aes(Site, beta_binomial_mean_percent),
             inherit.aes = FALSE, shape = 21, fill = "white", colour = "black",
             size = 2.7, stroke = 0.8) +
  scale_fill_manual(values = site_palette) +
  scale_colour_manual(values = site_palette) +
  labs(title = "a  Anti-Bd relative abundance",
       subtitle = "White circles: beta-binomial locality means and approximate 95% CIs",
       x = NULL, y = "Anti-Bd reads (% of total reads)") +
  theme_pub + theme(axis.text.x = element_text(angle = 24, hjust = 1),
                    legend.position = "none")

hill_panel <- function(community, tag) {
  z <- hill_summary[hill_summary$Community == community, ]
  ggplot(z, aes(Order_label, Mean, colour = Site, group = Site)) +
    geom_line(linewidth = 0.55, alpha = 0.75) +
    geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.08,
                  linewidth = 0.45) +
    geom_point(size = 2.2) +
    scale_colour_manual(values = site_palette) +
    labs(title = paste0(tag, "  ", community, " Hill diversity profile"),
         subtitle = sprintf("Coverage-standardized; common coverage = %.3f",
                            unique(round(z$Common_coverage, 3))),
         x = "Hill order", y = "Effective number of genera") +
    theme_pub + theme(legend.position = "bottom")
}

p_diversity <- p_abundance / (
  hill_panel("Anti-Bd community", "b") |
    hill_panel("Whole community", "c")
) +
  plot_layout(heights = c(1.05, 1)) +
  plot_annotation(
    title = "Anti-Bd relative abundance and genus diversity across localities",
    subtitle = "Hill D0 emphasizes richness, D1 common genera and evenness, and D2 dominant genera.",
    caption = paste0(
      "All 47 frogs are retained. Hill numbers are standardized to the minimum observed sample coverage within each community.\n",
      "Error bars in Hill panels are locality bootstrap 95% CIs (4,999 resamples)."
    ),
    theme = theme(plot.title = element_text(size = 11, face = "bold"),
                  plot.subtitle = element_text(size = 7.2),
                  plot.caption = element_text(size = 6, hjust = 0))
  )
save_pub(p_diversity, file.path(figure_dir, "figure_locality_abundance_hill_D0_D1_D2"),
         height_mm = 158)

# -------------------------------------------------------------------------
# Figure 3: whole/anti-Bd composition along estimated Bd-arrival year
# -------------------------------------------------------------------------

perm_years <- permutations_all(unname(site_year[site_levels]))
year_tests <- do.call(rbind, lapply(names(community_objects), function(nm) {
  coord <- community_objects[[nm]]
  D <- dist(coord)
  fit_stat <- function(y) {
    ad <- vegan::adonis2(D ~ y, permutations = 0)
    c(R2 = ad$R2[1], F = ad$F[1])
  }
  observed <- fit_stat(year)
  perm_F <- apply(perm_years, 1, function(py) {
    yp <- py[match(as.character(site), site_levels)]
    fit_stat(yp)["F"]
  })
  data.frame(
    Community = nm, R2 = unname(observed["R2"]),
    pseudo_F = unname(observed["F"]),
    Exact_p = mean(perm_F >= observed["F"] - 1e-12),
    Locality_n = 4L, Frog_n = nrow(coord), Exact_permutations = nrow(perm_years)
  )
}))
year_tests$Holm_p <- p.adjust(year_tests$Exact_p, method = "holm")
year_tests$Community_factor <- factor(
  year_tests$Community, levels = rev(names(community_objects))
)
write.csv(year_tests, file.path(result_dir, "robust_aitchison_bd_year_tests.csv"),
          row.names = FALSE)

trajectory_panel <- function(coord, nm, tag) {
  pc <- prcomp(coord, center = TRUE, scale. = FALSE)
  pct <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  z <- data.frame(Sample_ID = rownames(coord), PC1 = pc$x[, 1], PC2 = pc$x[, 2],
                  Site = site, Bd_arrival_year = year)
  cent <- aggregate(cbind(PC1, PC2, Bd_arrival_year) ~ Site, z, mean)
  cent <- cent[order(cent$Bd_arrival_year), ]
  hit <- year_tests[year_tests$Community == nm, ]
  ggplot(z, aes(PC1, PC2, colour = Site)) +
    geom_point(size = 1.45, alpha = 0.42) +
    geom_path(data = cent, aes(PC1, PC2), inherit.aes = FALSE,
              linewidth = 0.85, colour = "#4A4A4A",
              arrow = arrow(length = grid::unit(1.8, "mm"), type = "closed")) +
    geom_point(data = cent, aes(PC1, PC2, fill = Site), shape = 21,
               colour = "white", size = 3.6, stroke = 0.7) +
    geom_text_repel(data = cent,
                    aes(PC1, PC2,
                        label = paste0(as.character(Site), " ", round(Bd_arrival_year, 1))),
                    size = 1.75, seed = 20260817, min.segment.length = 0,
                    show.legend = FALSE) +
    scale_colour_manual(values = site_palette) +
    scale_fill_manual(values = site_palette) +
    labs(
      title = paste0(tag, "  ", nm),
      subtitle = sprintf("R² = %.3f; exact p = %.3f",
                         hit$R2, hit$Exact_p),
      x = sprintf("rCLR-PC1 (%.1f%%)", pct[1]),
      y = sprintf("rCLR-PC2 (%.1f%%)", pct[2])
    ) + theme_pub + theme(legend.position = "none")
}

p_effect <- ggplot(year_tests,
                   aes(R2, Community_factor,
                       colour = Community)) +
  geom_segment(aes(x = 0, xend = R2, yend = Community_factor),
    linewidth = 0.65, colour = "grey75") +
  geom_point(size = 3) +
  geom_text(aes(label = sprintf("R² %.3f\np %.3f", R2, Exact_p)),
            nudge_x = 0.025, hjust = 0, size = 2.4, colour = "black") +
  scale_colour_manual(values = community_palette) +
  scale_x_continuous(limits = c(0, max(year_tests$R2) + 0.14),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(title = "c  Bd-year effect sizes", x = "PERMANOVA R²", y = NULL) +
  theme_pub + theme(legend.position = "none")

p_year <- (
  trajectory_panel(full_rclr, "Whole community", "a") |
    trajectory_panel(anti_rclr, "Anti-Bd community", "b") |
    p_effect
) +
  plot_layout(widths = c(1.05, 1.05, 0.90)) +
  plot_annotation(
    title = "Community composition along estimated Bd-arrival year",
    subtitle = "Arrows connect locality centroids from earlier to later estimated arrival; they show trajectory, not causation.",
    caption = paste0(
      "Faint points are individual frogs; large points are locality centroids.\n",
      "Exact p values enumerate all 24 assignments of arrival years among four independent localities."
    ),
    theme = theme(plot.title = element_text(size = 11, face = "bold"),
                  plot.subtitle = element_text(size = 7.2),
                  plot.caption = element_text(size = 6, hjust = 0))
  )
# Keep this earlier PCA-style diagnostic under a distinct name.  The canonical
# two-panel CAP/dbRDA figure is generated by 22_wdstar_cap_bd_year.R.
save_pub(p_year, file.path(figure_dir, "figure_legacy_pca_composition_bd_arrival_year"),
         height_mm = 105)

cat("Completed core locality, abundance/diversity, and Bd-year composition figures.\n")
print(locality_tests)
print(year_tests)
