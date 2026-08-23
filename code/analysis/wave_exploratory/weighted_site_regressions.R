suppressPackageStartupMessages({
  library(phyloseq)
  library(iNEXT)
  library(glmmTMB)
  library(vegan)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(svglite)
  library(ragg)
})

set.seed(20260815)
result_dir <- "results/wave_exploratory"
figure_dir <- file.path(result_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

full_file <- "data/processed/ps_clean_47_rebuilt_from_01_filtered.rds"
anti_file <- "data/processed/blast_97id_90cov/ps_anti_bd_47_rebuilt_97id_90cov.rds"
wave_file <- file.path(result_dir, "site_metrics_single_wave.csv")
stopifnot(file.exists(full_file), file.exists(anti_file), file.exists(wave_file))

otu_matrix <- function(ps) {
  x <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) x <- t(x)
  x
}

ps_full <- readRDS(full_file)
ps_anti <- readRDS(anti_file)
full <- otu_matrix(ps_full)
anti <- otu_matrix(ps_anti)[rownames(full), , drop = FALSE]
md <- data.frame(sample_data(ps_full), check.names = FALSE, stringsAsFactors = FALSE)
md <- md[rownames(full), , drop = FALSE]

site_key <- c(
  "Lost Iguana Resort, La Fortuna, San Carlos" = "Lost Iguana",
  "Veragua Rainforest, Las Brisas de Veragua" = "Veragua",
  "PN Altos de Campana" = "Altos de Campana",
  "PN Soberanía, Las Cumbres" = "Soberanía"
)
site <- unname(site_key[md$Locality])
stopifnot(!anyNA(site), nrow(full) == 47L)

# Anti-Bd genus Hill D1. Coverage standardisation is requested at the minimum
# observed sample coverage. In this data set every genus-level assemblage has
# estimated coverage 1 (there are no genus singletons), so the standardised
# estimate is numerically identical to exp(Shannon entropy).
tax <- as(tax_table(ps_anti), "matrix")[colnames(anti), , drop = FALSE]
genus <- as.character(tax[, "Genus"])
family <- as.character(tax[, "Family"])
family[is.na(family) | trimws(family) == ""] <- "Bacteria"
missing_genus <- is.na(genus) | trimws(genus) == ""
genus[missing_genus] <- paste0("Unclassified_", family[missing_genus])
anti_genus <- t(rowsum(t(anti), group = genus, reorder = FALSE))
assemblages <- lapply(seq_len(nrow(anti_genus)), function(i) {
  z <- anti_genus[i, ]
  unname(z[z > 0])
})
names(assemblages) <- rownames(anti_genus)
coverage_info <- iNEXT::DataInfo(assemblages, datatype = "abundance")
common_coverage <- min(coverage_info$SC)
d1_est <- iNEXT::estimateD(
  assemblages, datatype = "abundance", base = "coverage",
  level = common_coverage, q = 1, nboot = 0
)
d1 <- setNames(d1_est$qD, d1_est$Assemblage)[rownames(full)]

# Whole-community CLR-Aitchison distance from each frog to its locality centroid.
clr <- log(full + 1)
clr <- clr - rowMeans(clr)
distance_to_centroid <- numeric(nrow(clr))
anti_clr <- log(anti + 1)
anti_clr <- anti_clr - rowMeans(anti_clr)
anti_distance_to_centroid <- numeric(nrow(anti_clr))
# Robust CLR (rclr) avoids adding a pseudocount to zeros and uses low-rank
# matrix completion. Euclidean distances in this coordinate space are robust
# Aitchison distances. These are retained as sensitivity endpoints.
whole_rclr <- vegan::decostand(full, method = "rclr")
anti_rclr <- vegan::decostand(anti, method = "rclr")
stopifnot(!anyNA(whole_rclr), !anyNA(anti_rclr),
          all(is.finite(whole_rclr)), all(is.finite(anti_rclr)))
whole_robust_distance_to_centroid <- numeric(nrow(whole_rclr))
anti_robust_distance_to_centroid <- numeric(nrow(anti_rclr))
for (s in unique(site)) {
  idx <- which(site == s)
  centroid <- colMeans(clr[idx, , drop = FALSE])
  distance_to_centroid[idx] <- sqrt(rowSums((clr[idx, , drop = FALSE] - centroid)^2))
  anti_centroid <- colMeans(anti_clr[idx, , drop = FALSE])
  anti_distance_to_centroid[idx] <-
    sqrt(rowSums((anti_clr[idx, , drop = FALSE] - anti_centroid)^2))
  whole_robust_centroid <- colMeans(whole_rclr[idx, , drop = FALSE])
  whole_robust_distance_to_centroid[idx] <- sqrt(rowSums(
    (whole_rclr[idx, , drop = FALSE] - whole_robust_centroid)^2
  ))
  anti_robust_centroid <- colMeans(anti_rclr[idx, , drop = FALSE])
  anti_robust_distance_to_centroid[idx] <- sqrt(rowSums(
    (anti_rclr[idx, , drop = FALSE] - anti_robust_centroid)^2
  ))
}

sample_metrics <- data.frame(
  Sample_ID = rownames(full), Site = factor(site),
  anti_reads = rowSums(anti), non_anti_reads = rowSums(full) - rowSums(anti),
  total_reads = rowSums(full), anti_bd_percent = 100 * rowSums(anti) / rowSums(full),
  anti_bd_genus_Hill_D1 = unname(d1),
  whole_community_distance_to_centroid = distance_to_centroid,
  anti_bd_distance_to_centroid = anti_distance_to_centroid,
  whole_community_robust_distance_to_centroid = whole_robust_distance_to_centroid,
  anti_bd_robust_distance_to_centroid = anti_robust_distance_to_centroid,
  stringsAsFactors = FALSE
)
stopifnot(all(sample_metrics$non_anti_reads >= 0), all(is.finite(d1)))
write.csv(sample_metrics, file.path(result_dir, "sample_three_metrics.csv"), row.names = FALSE)
write.csv(coverage_info, file.path(result_dir, "anti_bd_genus_sample_coverage.csv"), row.names = FALSE)

# Beta-binomial models separate locality effects on the expected proportion
# (conditional mean) from locality effects on beta-binomial overdispersion.
bb_full <- glmmTMB(
  cbind(anti_reads, non_anti_reads) ~ Site,
  dispformula = ~ Site, family = betabinomial(link = "logit"), data = sample_metrics
)
bb_no_mean_site <- update(bb_full, . ~ 1)
bb_no_disp_site <- update(bb_full, dispformula = ~ 1)

mean_lrt <- anova(bb_no_mean_site, bb_full)
disp_lrt <- anova(bb_no_disp_site, bb_full)
bb_tests <- data.frame(
  Test = c("Locality effect on mean anti-Bd relative abundance",
           "Locality effect on beta-binomial dispersion"),
  Null_model = c("constant mean; locality-specific dispersion",
                 "locality-specific mean; constant dispersion"),
  Alternative_model = c("locality-specific mean and dispersion",
                        "locality-specific mean and dispersion"),
  Chisq = c(mean_lrt$Chisq[2], disp_lrt$Chisq[2]),
  df = c(mean_lrt$`Chi Df`[2], disp_lrt$`Chi Df`[2]),
  p_value = c(mean_lrt$`Pr(>Chisq)`[2], disp_lrt$`Pr(>Chisq)`[2]),
  Sample_n = nrow(sample_metrics), Locality_n = nlevels(sample_metrics$Site),
  stringsAsFactors = FALSE
)
write.csv(bb_tests, file.path(result_dir, "beta_binomial_mean_and_dispersion_tests.csv"), row.names = FALSE)

site_levels <- levels(sample_metrics$Site)
new_sites <- data.frame(Site = factor(site_levels, levels = site_levels))
mean_pred <- predict(bb_full, newdata = new_sites, type = "response", se.fit = TRUE)
disp_pred <- predict(bb_full, newdata = new_sites, type = "disp", se.fit = TRUE)
bb_site <- data.frame(
  Site = site_levels,
  beta_binomial_mean_percent = 100 * as.numeric(mean_pred$fit),
  beta_binomial_mean_SE_percent = 100 * as.numeric(mean_pred$se.fit),
  beta_binomial_dispersion = as.numeric(disp_pred$fit),
  beta_binomial_dispersion_SE = as.numeric(disp_pred$se.fit),
  # glmmTMB reports beta-binomial precision phi. Convert it to the
  # intraclass-correlation/overdispersion scale rho = 1 / (phi + 1),
  # for which larger values mean stronger extra-binomial variation.
  beta_binomial_overdispersion_rho = 1 / (as.numeric(disp_pred$fit) + 1),
  beta_binomial_overdispersion_rho_SE = as.numeric(disp_pred$se.fit) /
    (as.numeric(disp_pred$fit) + 1)^2,
  stringsAsFactors = FALSE
)
write.csv(bb_site, file.path(result_dir, "beta_binomial_site_estimates.csv"), row.names = FALSE)

summarise_site <- function(z) {
  data.frame(
    Site = as.character(z$Site[1]), Sample_n = nrow(z),
    anti_raw_mean_percent = mean(z$anti_bd_percent),
    anti_raw_se_percent = sd(z$anti_bd_percent) / sqrt(nrow(z)),
    d1_mean = mean(z$anti_bd_genus_Hill_D1),
    d1_se = sd(z$anti_bd_genus_Hill_D1) / sqrt(nrow(z)),
    dispersion_mean = mean(z$whole_community_distance_to_centroid),
    dispersion_se = sd(z$whole_community_distance_to_centroid) / sqrt(nrow(z)),
    anti_dispersion_mean = mean(z$anti_bd_distance_to_centroid),
    anti_dispersion_se = sd(z$anti_bd_distance_to_centroid) / sqrt(nrow(z)),
    robust_dispersion_mean = mean(z$whole_community_robust_distance_to_centroid),
    robust_dispersion_se = sd(z$whole_community_robust_distance_to_centroid) / sqrt(nrow(z)),
    anti_robust_dispersion_mean = mean(z$anti_bd_robust_distance_to_centroid),
    anti_robust_dispersion_se = sd(z$anti_bd_robust_distance_to_centroid) / sqrt(nrow(z))
  )
}
site_summary <- do.call(rbind, lapply(split(sample_metrics, sample_metrics$Site), summarise_site))
rownames(site_summary) <- NULL
wave <- read.csv(wave_file, check.names = FALSE, stringsAsFactors = FALSE)
wave <- wave[, c("Site", "pred_wave_year")]
site_summary <- merge(site_summary, wave, by = "Site", all.x = TRUE, sort = FALSE)
site_summary <- merge(site_summary, bb_site, by = "Site", all.x = TRUE, sort = FALSE)
stopifnot(nrow(site_summary) == 4L, !anyNA(site_summary))
write.csv(site_summary, file.path(result_dir, "site_four_metrics.csv"), row.names = FALSE)

endpoint_definitions <- data.frame(
  Endpoint = c("Anti-Bd relative abundance", "Anti-Bd genus Hill D1 diversity",
               "Anti-Bd relative-abundance variation",
               "Whole-community within-locality heterogeneity",
               "Anti-Bd within-locality heterogeneity",
               "Whole-community robust within-locality heterogeneity",
               "Anti-Bd robust within-locality heterogeneity"),
  Role = c("Primary", "Primary", "Secondary", "Secondary", "Secondary",
           "Sensitivity", "Sensitivity"),
  Estimator = c("Beta-binomial conditional mean",
                "Coverage-standardized Hill number q = 1",
                "Beta-binomial extra-variation rho = 1 / (phi + 1)",
                "Mean CLR-Aitchison distance to locality centroid",
                "Mean anti-Bd CLR-Aitchison distance to locality centroid",
                "Mean robust Aitchison distance to locality centroid",
                "Mean anti-Bd robust Aitchison distance to locality centroid"),
  Biological_scope = c("Anti-Bd ASVs versus all bacterial reads",
                       "Anti-Bd genera", "Anti-Bd ASVs versus all bacterial reads",
                       "All bacterial ASVs", "Anti-Bd ASVs",
                       "All bacterial ASVs", "Anti-Bd ASVs"),
  Hypothesis_role = c("Tests the stated hypothesis", "Tests the stated hypothesis",
                      "Supporting heterogeneity analysis",
                      "Supporting community-heterogeneity analysis",
                      "Supporting anti-Bd community-heterogeneity analysis",
                      "Zero-aware sensitivity analysis of community heterogeneity",
                      "Zero-aware sensitivity analysis of anti-Bd heterogeneity"),
  stringsAsFactors = FALSE
)
write.csv(endpoint_definitions,
          file.path(result_dir, "analysis_endpoint_definitions.csv"), row.names = FALSE)

# The former data-adaptive smoothing-spline section is deliberately disabled.
# The confirmatory hypothesis-aligned trend models and figure are produced by
# test_bd_arrival_hypothesis.R.
print(bb_tests)
quit(save = "no", status = 0)

# Data-adaptive smoothing splines: no linear, exponential, or monotonic shape is
# selected in advance. With only four independent localities these are
# descriptive fits; smoothing is selected by leave-one-out cross-validation.
fit_spline <- function(y, label, transform = identity, inverse = identity,
                       scale_note = "identity") {
  ord <- order(site_summary$pred_wave_year)
  x <- site_summary$pred_wave_year[ord]
  yy <- site_summary[[y]][ord]
  yy_model <- transform(yy)
  w <- site_summary$Sample_n[ord]
  fit <- smooth.spline(x, yy_model, w = w, cv = TRUE, all.knots = TRUE)
  fitted_y <- inverse(predict(fit, x)$y)
  sse <- sum(w * (yy - fitted_y)^2)
  null_sse <- sum(w * (yy - weighted.mean(yy, w))^2)
  list(
    fit = fit, x = x, y = yy, w = w, label = label,
    inverse = inverse, scale_note = scale_note,
    pseudo_R2 = 1 - sse / null_sse,
    weighted_RMSE = sqrt(sse / sum(w)),
    effective_df = fit$df, spar = fit$spar, cv_criterion = fit$cv.crit
  )
}

fits <- list(
  fit_spline("beta_binomial_mean_percent",
             "Beta-binomial estimated anti-Bd relative abundance (%)",
             transform = function(x) qlogis(x / 100),
             inverse = function(x) 100 * plogis(x), scale_note = "logit proportion"),
  fit_spline("beta_binomial_overdispersion_rho",
             "Beta-binomial extra-variation parameter (rho)",
             transform = qlogis, inverse = plogis, scale_note = "logit rho"),
  fit_spline("d1_mean", "Coverage-standardized anti-Bd genus Hill D1",
             transform = log, inverse = exp, scale_note = "log D1"),
  fit_spline("dispersion_mean",
             "Mean whole-community Aitchison distance to locality centroid",
             transform = log, inverse = exp, scale_note = "log distance")
)
fit_results <- do.call(rbind, lapply(fits, function(z) data.frame(
  Outcome = z$label,
  Model = "weighted cubic smoothing spline; smoothing selected by leave-one-out CV",
  Fitting_scale = z$scale_note,
  Effective_df = z$effective_df, Smoothing_parameter = z$spar,
  CV_criterion = z$cv_criterion,
  Weighted_pseudo_R2 = z$pseudo_R2, Weighted_RMSE = z$weighted_RMSE,
  Independent_units = 4L, Weight = "frog sample count per locality",
  Inference_note = "descriptive exploratory fit; no curve p value with four localities",
  stringsAsFactors = FALSE
)))
write.csv(fit_results, file.path(result_dir, "four_metric_flexible_spline_results.csv"), row.names = FALSE)

pal <- c("Lost Iguana" = "#B44682", "Veragua" = "#166AA5",
         "Altos de Campana" = "#18835C", "Soberanía" = "#D07A00")
base_theme <- theme_classic(base_size = 8, base_family = "sans") +
  theme(axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4),
        plot.title = element_text(size = 9, face = "bold"),
        plot.subtitle = element_text(size = 6.8), legend.position = "none")

make_panel <- function(y, error, fit, title, ylab, tag, role) {
  grid <- data.frame(pred_wave_year = seq(min(fit$x), max(fit$x), length.out = 300))
  grid$fitted <- fit$inverse(predict(fit$fit, grid$pred_wave_year)$y)
  ggplot(site_summary, aes(pred_wave_year, .data[[y]], colour = Site)) +
    geom_errorbar(aes(ymin = pmax(0, .data[[y]] - .data[[error]]),
                      ymax = .data[[y]] + .data[[error]]),
                  width = 0.25, linewidth = 0.4) +
    geom_line(data = grid, aes(pred_wave_year, fitted), inherit.aes = FALSE,
              colour = "#6A3D9A", linewidth = 0.9) +
    geom_point(aes(size = Sample_n), alpha = 0.95) +
    geom_text_repel(aes(label = Site), size = 2.05, seed = 20260815,
                    min.segment.length = 0, box.padding = 0.25) +
    scale_colour_manual(values = pal) +
    scale_size_continuous(range = c(2.7, 4.8)) +
    scale_x_continuous(breaks = seq(1985, 2010, 5)) +
    labs(title = paste0(tag, "  ", title),
         subtitle = sprintf("%s; CV-selected spline; effective df = %.2f",
                            role, fit$effective_df),
         x = "Estimated Bd arrival year", y = ylab) + base_theme
}

p1 <- make_panel("beta_binomial_mean_percent", "beta_binomial_mean_SE_percent", fits[[1]],
                 "Anti-Bd relative abundance", "Beta-binomial estimate (%)", "a",
                 "Primary endpoint") +
  expand_limits(y = 0)
p2 <- make_panel("d1_mean", "d1_se", fits[[3]], "Anti-Bd Hill D1 diversity",
                 "Effective number of common genera", "b", "Primary endpoint") +
  expand_limits(y = 0)
p3 <- make_panel("beta_binomial_overdispersion_rho",
                 "beta_binomial_overdispersion_rho_SE", fits[[2]],
                 "Anti-Bd relative-abundance variation",
                 expression("Extra-binomial variation " * rho), "c",
                 "Secondary endpoint") + expand_limits(y = 0)
p4 <- make_panel("dispersion_mean", "dispersion_se", fits[[4]],
                 "Whole-community heterogeneity", "Mean Aitchison distance", "d",
                 "Secondary endpoint")

figure <- (p1 | p2) / (p3 | p4) +
  plot_annotation(
    title = "Bd arrival time, anti-Bd abundance and diversity",
    subtitle = paste0("Hypothesis: anti-Bd relative abundance and genus diversity vary systematically ",
                      "with Bd arrival year; four independent localities."),
    caption = paste0("Curves are descriptive and have no inferential p value. Error bars are s.e.m.\n",
                     "Panels a-b use beta-binomial locality estimates; larger rho indicates greater ",
                     "extra-binomial variation. Hill D1 uses anti-Bd genera."),
    theme = theme(plot.title = element_text(size = 11, face = "bold"),
                  plot.subtitle = element_text(size = 7.2),
                  plot.caption = element_text(size = 6.1, hjust = 0))
  )

save_plot <- function(plot, stem, width_mm = 183, height_mm = 155) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  ggsave(paste0(stem, ".png"), plot, width = w, height = h, dpi = 600, bg = "white")
  ragg::agg_tiff(paste0(stem, ".tiff"), width = w, height = h,
                 units = "in", res = 600, background = "white")
  print(plot); dev.off()
  svglite::svglite(paste0(stem, ".svg"), width = w, height = h)
  print(plot); dev.off()
  # Base PDF is used because this macOS R installation reports Cairo support
  # but lacks the external XQuartz libraries required to load cairo.so.
  grDevices::pdf(paste0(stem, ".pdf"), width = w, height = h,
                 family = "Helvetica", useDingbats = FALSE)
  print(plot); dev.off()
}

save_plot(figure, file.path(figure_dir, "figure_four_metrics_flexible_fits"))
print(bb_tests)
print(fit_results)
