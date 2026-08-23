suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(svglite)
  library(ragg)
})

set.seed(20260816)
result_dir <- "results/bd_arrival_community"
figure_dir <- file.path(result_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

full_file <- "data/processed/ps_clean_47_rebuilt_from_01_filtered.rds"
anti_file <- "data/processed/blast_97id_90cov/ps_anti_bd_47_rebuilt_97id_90cov.rds"
wave_file <- "results/wave_exploratory/site_metrics_single_wave.csv"
stopifnot(file.exists(full_file), file.exists(anti_file), file.exists(wave_file))

otu_samples_rows <- function(ps) {
  x <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) x <- t(x)
  x
}

ps_full <- readRDS(full_file)
ps_anti <- readRDS(anti_file)
full_asv <- otu_samples_rows(ps_full)
anti_asv <- otu_samples_rows(ps_anti)[rownames(full_asv), , drop = FALSE]
md <- data.frame(sample_data(ps_full), check.names = FALSE,
                 stringsAsFactors = FALSE)[rownames(full_asv), , drop = FALSE]

site_key <- c(
  "Lost Iguana Resort, La Fortuna, San Carlos" = "Lost Iguana",
  "Veragua Rainforest, Las Brisas de Veragua" = "Veragua",
  "PN Altos de Campana" = "Altos de Campana",
  "PN Soberanía, Las Cumbres" = "Soberanía"
)
site_levels <- c("Lost Iguana", "Veragua", "Altos de Campana", "Soberanía")
site <- factor(unname(site_key[md$Locality]), levels = site_levels)
stopifnot(!anyNA(site), nrow(full_asv) == 47L)

wave <- read.csv(wave_file, stringsAsFactors = FALSE, check.names = FALSE)
site_year <- setNames(wave$pred_wave_year, wave$Site)
stopifnot(all(site_levels %in% names(site_year)))
year <- unname(site_year[as.character(site)])

# The anti-Bd hypothesis is genus-level. Unclassified genera are retained as
# family-specific groups rather than discarded.
anti_tax <- as(tax_table(ps_anti), "matrix")[colnames(anti_asv), , drop = FALSE]
genus <- as.character(anti_tax[, "Genus"])
family <- as.character(anti_tax[, "Family"])
family[is.na(family) | trimws(family) == ""] <- "Bacteria"
missing_genus <- is.na(genus) | trimws(genus) == ""
genus[missing_genus] <- paste0("Unclassified_", family[missing_genus])
anti_genus <- t(rowsum(t(anti_asv), group = genus, reorder = FALSE))
stopifnot(all(rowSums(anti_genus) > 0))

clr_transform <- function(x, pseudocount = 1) {
  z <- log(x + pseudocount)
  z - rowMeans(z)
}

coordinates <- list(
  `Whole community|CLR-Aitchison` = clr_transform(full_asv),
  `Whole community|Robust Aitchison` = vegan::decostand(full_asv, method = "rclr"),
  `Anti-Bd genera|CLR-Aitchison` = clr_transform(anti_genus),
  `Anti-Bd genera|Robust Aitchison` = vegan::decostand(anti_genus, method = "rclr")
)

permutations <- function(x) {
  if (length(x) == 1L) return(matrix(x, nrow = 1L))
  do.call(rbind, lapply(seq_along(x), function(i) cbind(x[i], permutations(x[-i]))))
}
perm_years <- permutations(unname(site_year[site_levels]))

# Continuous-year PERMANOVA with the locality as the exchangeable unit. All
# frogs at one locality receive the same permuted year, avoiding pseudoreplication.
exact_year_permanova <- function(coord, community, metric) {
  D <- dist(coord)
  fit_one <- function(y) {
    ad <- vegan::adonis2(D ~ y, permutations = 0)
    c(R2 = ad$R2[1], pseudo_F = ad$F[1])
  }
  observed <- fit_one(year)
  perm_stats <- t(apply(perm_years, 1, function(py) {
    y_perm <- py[match(as.character(site), site_levels)]
    fit_one(y_perm)
  }))
  exact_p <- mean(perm_stats[, "pseudo_F"] >= observed["pseudo_F"] - 1e-12)
  list(
    result = data.frame(
      Community = community, Metric = metric,
      R2 = unname(observed["R2"]), pseudo_F = unname(observed["pseudo_F"]),
      Exact_p = exact_p, Locality_n = 4L, Frog_n = nrow(coord),
      Permutations = nrow(perm_years),
      Permutation_unit = "Locality (all frogs at a locality moved together)",
      stringsAsFactors = FALSE
    ),
    permutations = data.frame(
      Community = community, Metric = metric,
      Permutation = seq_len(nrow(perm_stats)),
      R2 = perm_stats[, "R2"], pseudo_F = perm_stats[, "pseudo_F"],
      stringsAsFactors = FALSE
    )
  )
}

composition_fits <- lapply(names(coordinates), function(nm) {
  bits <- strsplit(nm, "\\|", fixed = FALSE)[[1]]
  exact_year_permanova(coordinates[[nm]], bits[1], bits[2])
})
composition_tests <- do.call(rbind, lapply(composition_fits, `[[`, "result"))
composition_tests$Analysis_role <- ifelse(composition_tests$Metric == "CLR-Aitchison",
                                           "Primary", "Sensitivity")
composition_tests$Holm_p_within_metric <- ave(
  composition_tests$Exact_p, composition_tests$Metric,
  FUN = function(p) p.adjust(p, method = "holm")
)
composition_tests$Conclusion_at_0_05 <- ifelse(
  composition_tests$Holm_p_within_metric < 0.05,
  "Evidence of association with estimated Bd arrival year",
  "Insufficient evidence to reject no association"
)
write.csv(composition_tests,
          file.path(result_dir, "bd_arrival_composition_exact_permanova.csv"),
          row.names = FALSE)
write.csv(do.call(rbind, lapply(composition_fits, `[[`, "permutations")),
          file.path(result_dir, "bd_arrival_composition_permutation_distribution.csv"),
          row.names = FALSE)

# Bias-adjusted distance to locality centroid. For Euclidean coordinate spaces,
# sqrt(n/(n-1)) is the small-sample correction used by betadisper.
centroid_distances <- function(coord) {
  out <- numeric(nrow(coord))
  for (s in site_levels) {
    idx <- which(site == s)
    ctr <- colMeans(coord[idx, , drop = FALSE])
    raw <- sqrt(rowSums((coord[idx, , drop = FALSE] - ctr)^2))
    out[idx] <- raw * sqrt(length(idx) / (length(idx) - 1))
  }
  out
}

bootstrap_mean_dispersion <- function(coord, B = 4999L, seed = 1L) {
  set.seed(seed)
  do.call(rbind, lapply(site_levels, function(s) {
    idx <- which(site == s)
    observed <- mean(centroid_distances(coord)[idx])
    boot <- replicate(B, {
      take <- sample(idx, length(idx), replace = TRUE)
      z <- coord[take, , drop = FALSE]
      ctr <- colMeans(z)
      mean(sqrt(rowSums((z - ctr)^2)) * sqrt(length(idx) / (length(idx) - 1)))
    })
    data.frame(Site = s, Sample_n = length(idx), Mean_distance = observed,
               CI_low = unname(quantile(boot, 0.025)),
               CI_high = unname(quantile(boot, 0.975)))
  }))
}

dispersion_tables <- Map(function(coord, nm, seed) {
  bits <- strsplit(nm, "\\|", fixed = FALSE)[[1]]
  z <- bootstrap_mean_dispersion(coord, seed = seed)
  z$Community <- bits[1]
  z$Metric <- bits[2]
  z
}, coordinates, names(coordinates), 20260816 + seq_along(coordinates))
dispersion_summary <- do.call(rbind, dispersion_tables)
dispersion_summary$Bd_arrival_year <- unname(site_year[dispersion_summary$Site])
write.csv(dispersion_summary,
          file.path(result_dir, "within_locality_dispersion_bootstrap95.csv"),
          row.names = FALSE)

dispersion_tests <- do.call(rbind, lapply(names(coordinates), function(nm) {
  bits <- strsplit(nm, "\\|", fixed = FALSE)[[1]]
  D <- dist(coordinates[[nm]])
  bd <- vegan::betadisper(D, site, type = "centroid", bias.adjust = TRUE)
  pt <- vegan::permutest(bd, permutations = 9999)
  data.frame(Community = bits[1], Metric = bits[2], F = unname(pt$tab[1, "F"]),
             p_value = unname(pt$tab[1, "Pr(>F)"]), Permutations = 9999L)
}))
dispersion_tests$BH_FDR_within_metric <- ave(
  dispersion_tests$p_value, dispersion_tests$Metric,
  FUN = function(p) p.adjust(p, method = "BH")
)
write.csv(dispersion_tests,
          file.path(result_dir, "within_locality_dispersion_betadisper.csv"),
          row.names = FALSE)

# Ordination figure for the two prespecified CLR-Aitchison composition tests.
palette_site <- c("Lost Iguana" = "#B44682", "Veragua" = "#166AA5",
                  "Altos de Campana" = "#18835C", "Soberanía" = "#D07A00")
theme_pub <- theme_classic(base_size = 8, base_family = "sans") +
  theme(axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4),
        plot.title = element_text(size = 9, face = "bold"),
        plot.subtitle = element_text(size = 6.8),
        legend.position = "bottom", legend.title = element_blank())

ordination_panel <- function(coord, community, tag) {
  pc <- prcomp(coord, center = FALSE, scale. = FALSE)
  variance <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  scores <- data.frame(Sample_ID = rownames(coord), PC1 = pc$x[, 1], PC2 = pc$x[, 2],
                       Site = site, Bd_arrival_year = year)
  cent <- aggregate(cbind(PC1, PC2, Bd_arrival_year) ~ Site, scores, mean)
  cent <- cent[order(cent$Bd_arrival_year), ]
  hit <- composition_tests[composition_tests$Community == community &
                             composition_tests$Metric == "CLR-Aitchison", ]
  ggplot(scores, aes(PC1, PC2, colour = Site)) +
    stat_ellipse(aes(group = Site), type = "norm", level = 0.68,
                 linewidth = 0.45, alpha = 0.75, show.legend = FALSE) +
    geom_point(size = 1.7, alpha = 0.65) +
    geom_path(data = cent, aes(PC1, PC2), inherit.aes = FALSE,
              colour = "#4D4D4D", linewidth = 0.65,
              arrow = arrow(length = grid::unit(1.8, "mm"), type = "closed")) +
    geom_point(data = cent, aes(PC1, PC2, colour = Site), size = 3.2,
               shape = 21, fill = "white", stroke = 0.9) +
    geom_text_repel(data = cent,
                    aes(PC1, PC2, label = sprintf("%.1f", Bd_arrival_year), colour = Site),
                    size = 2.15, seed = 20260816, show.legend = FALSE,
                    min.segment.length = 0) +
    scale_colour_manual(values = palette_site) +
    labs(title = paste0(tag, "  ", community),
         subtitle = sprintf("CLR-Aitchison PERMANOVA: R² = %.3f; exact p = %.3f; Holm p = %.3f",
                            hit$R2, hit$Exact_p, hit$Holm_p_within_metric),
         x = sprintf("CLR-PC1 (%.1f%%)", variance[1]),
         y = sprintf("CLR-PC2 (%.1f%%)", variance[2])) + theme_pub
}

p_comp <- ordination_panel(coordinates[["Whole community|CLR-Aitchison"]],
                           "Whole community", "a") |
  ordination_panel(coordinates[["Anti-Bd genera|CLR-Aitchison"]],
                   "Anti-Bd genera", "b")
p_comp <- p_comp + plot_layout(guides = "collect") +
  plot_annotation(
    title = "Community composition along estimated Bd-arrival year",
    subtitle = "Labels give estimated arrival year; arrows connect locality centroids from earlier to later arrival.",
    caption = paste0("Each point is one frog (n = 47); independent exposure units are four localities. ",
                     "Exact p values enumerate all 24 locality-level assignments of arrival years."),
    theme = theme(plot.title = element_text(size = 11, face = "bold"),
                  plot.subtitle = element_text(size = 7.2),
                  plot.caption = element_text(size = 6.1, hjust = 0))) &
  theme(legend.position = "bottom")

effect_plot_data <- composition_tests
effect_plot_data$Label <- sprintf("exact p = %.3f; Holm p = %.3f",
                                  effect_plot_data$Exact_p,
                                  effect_plot_data$Holm_p_within_metric)
effect_plot_data$Endpoint <- factor(
  paste(effect_plot_data$Community, effect_plot_data$Metric, sep = ": "),
  levels = rev(c("Whole community: CLR-Aitchison",
                 "Anti-Bd genera: CLR-Aitchison",
                 "Whole community: Robust Aitchison",
                 "Anti-Bd genera: Robust Aitchison")))
p_effect <- ggplot(effect_plot_data,
                   aes(R2, Endpoint, colour = Community, shape = Analysis_role)) +
  geom_segment(aes(x = 0, xend = R2, yend = Endpoint),
               colour = "grey78", linewidth = 0.6) +
  geom_point(size = 3.1) +
  geom_text(aes(label = Label), hjust = -0.08, size = 2.35,
            colour = "#333333", show.legend = FALSE) +
  scale_colour_manual(values = c("Whole community" = "#166AA5",
                                 "Anti-Bd genera" = "#B44682")) +
  scale_shape_manual(values = c("Primary" = 16, "Sensitivity" = 17)) +
  scale_x_continuous(limits = c(0, 0.56), breaks = seq(0, 0.5, 0.1),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Association of community composition with estimated Bd-arrival year",
       subtitle = "Continuous-year PERMANOVA with all 24 locality-level permutations",
       x = expression(PERMANOVA~R^2), y = NULL, colour = NULL, shape = NULL,
       caption = paste0("CLR-Aitchison analyses are prespecified primary tests; robust Aitchison analyses test sensitivity to sparsity and zero handling.\n",
                        "Holm correction is applied to the two community endpoints within each distance method.")) +
  theme_pub + theme(legend.position = "bottom",
                    axis.text.y = element_text(size = 7.2))

dispersion_summary$Community <- factor(
  dispersion_summary$Community, levels = c("Whole community", "Anti-Bd genera"))
dispersion_summary$Metric <- factor(
  dispersion_summary$Metric, levels = c("CLR-Aitchison", "Robust Aitchison"))
dispersion_summary$Site <- factor(
  dispersion_summary$Site,
  levels = names(sort(site_year[site_levels])))
p_disp <- ggplot(dispersion_summary,
                 aes(Site, Mean_distance, colour = Site)) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                width = 0.12, linewidth = 0.45) +
  geom_point(size = 2.5) +
  facet_grid(Metric ~ Community, scales = "free_y") +
  scale_colour_manual(values = palette_site) +
  labs(title = "Within-locality community heterogeneity",
       subtitle = "Bias-adjusted mean distance to locality centroid",
       x = NULL, y = NULL,
       caption = "Error bars are locality-wise nonparametric bootstrap 95% confidence intervals (4,999 resamples).") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        strip.background = element_blank(), strip.text = element_text(face = "bold"),
        legend.position = "none")

save_pub <- function(plot, stem, width_mm = 183, height_mm = 105) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  ggsave(paste0(stem, ".png"), plot, width = w, height = h, dpi = 600, bg = "white")
  svglite::svglite(paste0(stem, ".svg"), width = w, height = h); print(plot); dev.off()
  grDevices::pdf(paste0(stem, ".pdf"), width = w, height = h,
                 family = "Helvetica", useDingbats = FALSE); print(plot); dev.off()
  ragg::agg_tiff(paste0(stem, ".tiff"), width = w, height = h,
                 units = "in", res = 600, background = "white"); print(plot); dev.off()
}

save_pub(p_effect, file.path(figure_dir, "figure_bd_arrival_composition_permanova"),
         height_mm = 92)
save_pub(p_disp, file.path(figure_dir, "figure_within_locality_dispersion_improved"),
         height_mm = 132)

methods <- data.frame(
  Question = c("Within-locality consistency", "Bd-arrival association with composition"),
  Primary_analysis = c(
    "Bias-adjusted CLR-Aitchison distance to locality centroid with bootstrap 95% CI and PERMDISP",
    "Continuous-year CLR-Aitchison PERMANOVA with exact locality-level permutation"
  ),
  Sensitivity_analysis = c("Robust Aitchison dispersion", "Robust Aitchison PERMANOVA"),
  Independent_unit = c("Frog for dispersion distribution; locality for Bd-year association",
                       "Locality for permutation of Bd arrival year"),
  stringsAsFactors = FALSE
)
write.csv(methods, file.path(result_dir, "analysis_design_summary.csv"), row.names = FALSE)

cat("Composition tests\n")
print(composition_tests)
cat("\nDispersion tests\n")
print(dispersion_tests)
