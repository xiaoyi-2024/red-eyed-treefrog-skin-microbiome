suppressPackageStartupMessages({
  library(glmmTMB)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(svglite)
  library(ragg)
})

set.seed(20260815)
result_dir <- "results/wave_exploratory"
sample_file <- file.path(result_dir, "sample_three_metrics.csv")
site_file <- file.path(result_dir, "site_four_metrics.csv")
stopifnot(file.exists(sample_file), file.exists(site_file))

sample_data <- read.csv(sample_file, stringsAsFactors = FALSE, check.names = FALSE)
site_data <- read.csv(site_file, stringsAsFactors = FALSE, check.names = FALSE)
site_data <- site_data[order(site_data$Site), , drop = FALSE]
sample_data$Site <- factor(sample_data$Site, levels = site_data$Site)
sample_data$wave_year <- site_data$pred_wave_year[match(sample_data$Site, site_data$Site)]
sample_data$year10 <- (sample_data$wave_year - mean(site_data$pred_wave_year)) / 10
stopifnot(nrow(site_data) == 4L, nrow(sample_data) == 47L, !anyNA(sample_data$wave_year))

permutations <- function(x) {
  if (length(x) == 1L) return(matrix(x, nrow = 1L))
  do.call(rbind, lapply(seq_along(x), function(i) {
    cbind(x[i], permutations(x[-i]))
  }))
}

lrt_stat <- function(null, alt) {
  max(0, 2 * (as.numeric(logLik(alt)) - as.numeric(logLik(null))))
}

# Primary endpoint 1: anti-Bd relative abundance. The beta-binomial response
# retains each frog's anti-Bd and non-anti-Bd read counts. Locality-specific
# dispersion is retained in both models. The year effect is tested by all 4!
# assignments of the four observed arrival years to localities.
bb_mean_null <- glmmTMB(
  cbind(anti_reads, non_anti_reads) ~ 1,
  dispformula = ~ Site, family = betabinomial(link = "logit"), data = sample_data
)
bb_mean_alt <- update(bb_mean_null, . ~ year10)
mean_obs <- lrt_stat(bb_mean_null, bb_mean_alt)

site_year <- setNames(site_data$pred_wave_year, site_data$Site)
perm_mat <- permutations(unname(site_year))
mean_perm_stats <- numeric(nrow(perm_mat))
mean_perm_beta <- numeric(nrow(perm_mat))
for (i in seq_len(nrow(perm_mat))) {
  perm_map <- setNames(perm_mat[i, ], names(site_year))
  z <- sample_data
  z$year10 <- (unname(perm_map[as.character(z$Site)]) - mean(site_year)) / 10
  alt <- update(bb_mean_null, . ~ year10, data = z)
  mean_perm_stats[i] <- lrt_stat(bb_mean_null, alt)
  mean_perm_beta[i] <- fixef(alt)$cond["year10"]
}
mean_exact_p <- mean(mean_perm_stats >= mean_obs - 1e-10)
mean_coef <- summary(bb_mean_alt)$coefficients$cond["year10", ]

# Secondary endpoint 1: extra-binomial variation. Conditional means are free
# to differ among localities in both models; the beta-binomial precision phi is
# constant under the null and changes with arrival year under the alternative.
# Because rho = 1/(phi+1), a positive phi slope means declining extra variation.
bb_disp_null <- glmmTMB(
  cbind(anti_reads, non_anti_reads) ~ Site,
  dispformula = ~ 1, family = betabinomial(link = "logit"), data = sample_data
)
bb_disp_alt <- update(bb_disp_null, dispformula = ~ year10)
disp_obs <- lrt_stat(bb_disp_null, bb_disp_alt)
disp_perm_stats <- numeric(nrow(perm_mat))
disp_perm_beta <- numeric(nrow(perm_mat))
for (i in seq_len(nrow(perm_mat))) {
  perm_map <- setNames(perm_mat[i, ], names(site_year))
  z <- sample_data
  z$year10 <- (unname(perm_map[as.character(z$Site)]) - mean(site_year)) / 10
  alt <- update(bb_disp_null, dispformula = ~ year10, data = z)
  disp_perm_stats[i] <- lrt_stat(bb_disp_null, alt)
  disp_perm_beta[i] <- fixef(alt)$disp["year10"]
}
disp_exact_p <- mean(disp_perm_stats >= disp_obs - 1e-10)
disp_coef <- summary(bb_disp_alt)$coefficients$disp["year10", ]

# For positive continuous endpoints, fit one multiplicative trend parameter on
# the log scale. Inverse-variance weights retain the precision of each locality
# mean. Exact p values use all 24 assignments of years to localities.
fit_log_trend <- function(y, se) {
  w <- 1 / se^2
  x <- (site_data$pred_wave_year - mean(site_data$pred_wave_year)) / 10
  null <- lm(log(y) ~ 1, weights = w)
  alt <- lm(log(y) ~ x, weights = w)
  observed <- max(0, deviance(null) - deviance(alt))
  perm_stats <- numeric(nrow(perm_mat))
  perm_beta <- numeric(nrow(perm_mat))
  for (i in seq_len(nrow(perm_mat))) {
    xp <- (perm_mat[i, ] - mean(site_data$pred_wave_year)) / 10
    fp <- lm(log(y) ~ xp, weights = w)
    perm_stats[i] <- max(0, deviance(null) - deviance(fp))
    perm_beta[i] <- coef(fp)["xp"]
  }
  list(
    fit = alt, beta = unname(coef(alt)["x"]), ratio = exp(unname(coef(alt)["x"])),
    statistic = observed, exact_p = mean(perm_stats >= observed - 1e-12),
    perm_stats = perm_stats, perm_beta = perm_beta
  )
}

d1_trend <- fit_log_trend(site_data$d1_mean, site_data$d1_se)
aitch_trend <- fit_log_trend(site_data$dispersion_mean, site_data$dispersion_se)
anti_aitch_trend <- fit_log_trend(site_data$anti_dispersion_mean,
                                  site_data$anti_dispersion_se)
robust_aitch_trend <- fit_log_trend(site_data$robust_dispersion_mean,
                                    site_data$robust_dispersion_se)
anti_robust_aitch_trend <- fit_log_trend(site_data$anti_robust_dispersion_mean,
                                         site_data$anti_robust_dispersion_se)

# Exact rank tests are retained as nonparametric sensitivity analyses.
d1_test <- cor.test(site_data$pred_wave_year, site_data$d1_mean,
                    method = "spearman", exact = TRUE, alternative = "two.sided")
aitch_test <- cor.test(site_data$pred_wave_year, site_data$dispersion_mean,
                       method = "spearman", exact = TRUE, alternative = "two.sided")
mean_rank <- cor.test(site_data$pred_wave_year, site_data$beta_binomial_mean_percent,
                      method = "spearman", exact = TRUE, alternative = "two.sided")
rho_rank <- cor.test(site_data$pred_wave_year,
                     site_data$beta_binomial_overdispersion_rho,
                     method = "spearman", exact = TRUE, alternative = "two.sided")
anti_aitch_test <- cor.test(site_data$pred_wave_year, site_data$anti_dispersion_mean,
                            method = "spearman", exact = TRUE,
                            alternative = "two.sided")
robust_aitch_test <- cor.test(site_data$pred_wave_year,
                             site_data$robust_dispersion_mean,
                             method = "spearman", exact = TRUE,
                             alternative = "two.sided")
anti_robust_aitch_test <- cor.test(site_data$pred_wave_year,
                                  site_data$anti_robust_dispersion_mean,
                                  method = "spearman", exact = TRUE,
                                  alternative = "two.sided")

tests <- data.frame(
  Endpoint = c(
    "Anti-Bd relative abundance",
    "Anti-Bd genus Hill D1 diversity",
    "Anti-Bd relative-abundance extra variation",
    "Whole-community within-locality Aitchison dispersion",
    "Anti-Bd within-locality Aitchison dispersion",
    "Whole-community within-locality robust Aitchison dispersion",
    "Anti-Bd within-locality robust Aitchison dispersion"
  ),
  Role = c("Primary", "Primary", "Secondary", "Secondary", "Secondary",
           "Sensitivity", "Sensitivity"),
  Test = c(
    "Beta-binomial year effect; exact permutation of years among four localities",
    "Inverse-variance weighted log-linear year effect; exact permutation among four localities",
    "Beta-binomial dispersion year effect; exact permutation of years among four localities",
    "Inverse-variance weighted log-linear year effect; exact permutation among four localities",
    "Inverse-variance weighted log-linear year effect; exact permutation among four localities",
    "Inverse-variance weighted log-linear year effect; exact permutation among four localities",
    "Inverse-variance weighted log-linear year effect; exact permutation among four localities"
  ),
  Effect = c(
    unname(mean_coef["Estimate"]), d1_trend$beta,
    unname(disp_coef["Estimate"]), aitch_trend$beta, anti_aitch_trend$beta,
    robust_aitch_trend$beta, anti_robust_aitch_trend$beta
  ),
  Effect_definition = c(
    "log odds ratio per 10-year later Bd arrival",
    "log D1 ratio per 10-year later Bd arrival",
    "log precision ratio per 10-year later Bd arrival; positive means lower rho",
    "log whole-community Aitchison-distance ratio per 10-year later Bd arrival",
    "log anti-Bd Aitchison-distance ratio per 10-year later Bd arrival",
    "log whole-community robust-Aitchison-distance ratio per 10-year later Bd arrival",
    "log anti-Bd robust-Aitchison-distance ratio per 10-year later Bd arrival"
  ),
  Effect_on_ratio_scale = c(exp(mean_coef["Estimate"]), d1_trend$ratio,
                            exp(disp_coef["Estimate"]), aitch_trend$ratio,
                            anti_aitch_trend$ratio, robust_aitch_trend$ratio,
                            anti_robust_aitch_trend$ratio),
  Statistic = c(mean_obs, d1_trend$statistic, disp_obs, aitch_trend$statistic,
                anti_aitch_trend$statistic, robust_aitch_trend$statistic,
                anti_robust_aitch_trend$statistic),
  Exact_p = c(mean_exact_p, d1_trend$exact_p, disp_exact_p, aitch_trend$exact_p,
              anti_aitch_trend$exact_p, robust_aitch_trend$exact_p,
              anti_robust_aitch_trend$exact_p),
  Independent_localities = 4L,
  Frog_samples = rep(47L, 7L),
  stringsAsFactors = FALSE
)
tests$Multiplicity_family <- c("Two co-primary endpoints", "Two co-primary endpoints",
                               rep("Three secondary endpoints", 3L),
                               rep("Two robust-Aitchison sensitivity endpoints", 2L))
tests$Holm_p <- NA_real_
tests$Holm_p[tests$Role == "Primary"] <- p.adjust(tests$Exact_p[tests$Role == "Primary"],
                                                   method = "holm")
tests$Holm_p[tests$Role == "Secondary"] <- p.adjust(tests$Exact_p[tests$Role == "Secondary"],
                                                     method = "holm")
tests$Holm_p[tests$Role == "Sensitivity"] <- p.adjust(
  tests$Exact_p[tests$Role == "Sensitivity"], method = "holm"
)
tests$Conclusion_at_0_05 <- ifelse(
  tests$Holm_p < 0.05, "Evidence of an association",
  "Insufficient evidence to reject no association"
)
write.csv(tests, file.path(result_dir, "bd_arrival_hypothesis_tests.csv"), row.names = FALSE)

permutation_results <- rbind(
  data.frame(Test = "Beta-binomial mean", Permutation = seq_len(nrow(perm_mat)),
             LRT = mean_perm_stats, Coefficient_per_10_years = mean_perm_beta),
  data.frame(Test = "Beta-binomial dispersion", Permutation = seq_len(nrow(perm_mat)),
             LRT = disp_perm_stats, Coefficient_per_10_years = disp_perm_beta),
  data.frame(Test = "Hill D1 log-linear trend", Permutation = seq_len(nrow(perm_mat)),
             LRT = d1_trend$perm_stats, Coefficient_per_10_years = d1_trend$perm_beta),
  data.frame(Test = "Aitchison dispersion log-linear trend", Permutation = seq_len(nrow(perm_mat)),
             LRT = aitch_trend$perm_stats, Coefficient_per_10_years = aitch_trend$perm_beta),
  data.frame(Test = "Anti-Bd Aitchison dispersion log-linear trend",
             Permutation = seq_len(nrow(perm_mat)),
             LRT = anti_aitch_trend$perm_stats,
             Coefficient_per_10_years = anti_aitch_trend$perm_beta),
  data.frame(Test = "Robust Aitchison dispersion log-linear trend",
             Permutation = seq_len(nrow(perm_mat)),
             LRT = robust_aitch_trend$perm_stats,
             Coefficient_per_10_years = robust_aitch_trend$perm_beta),
  data.frame(Test = "Anti-Bd robust Aitchison dispersion log-linear trend",
             Permutation = seq_len(nrow(perm_mat)),
             LRT = anti_robust_aitch_trend$perm_stats,
             Coefficient_per_10_years = anti_robust_aitch_trend$perm_beta)
)
write.csv(permutation_results,
          file.path(result_dir, "bd_arrival_exact_permutation_distributions.csv"),
          row.names = FALSE)

sensitivity <- data.frame(
  Endpoint = c("Anti-Bd relative abundance", "Anti-Bd genus Hill D1 diversity",
               "Anti-Bd extra variation rho", "Whole-community Aitchison dispersion",
               "Anti-Bd Aitchison dispersion",
               "Whole-community robust Aitchison dispersion",
               "Anti-Bd robust Aitchison dispersion"),
  Spearman_rho = c(unname(mean_rank$estimate), unname(d1_test$estimate),
                   unname(rho_rank$estimate), unname(aitch_test$estimate),
                   unname(anti_aitch_test$estimate),
                   unname(robust_aitch_test$estimate),
                   unname(anti_robust_aitch_test$estimate)),
  Exact_two_sided_p = c(mean_rank$p.value, d1_test$p.value,
                        rho_rank$p.value, aitch_test$p.value,
                        anti_aitch_test$p.value, robust_aitch_test$p.value,
                        anti_robust_aitch_test$p.value),
  Independent_localities = 4L
)
write.csv(sensitivity,
          file.path(result_dir, "bd_arrival_rank_sensitivity_tests.csv"), row.names = FALSE)

# Hypothesis-aligned figure: one prespecified year coefficient per endpoint.
# No smoothing spline or data-selected curve complexity is used.
figure_dir <- file.path(result_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
grid <- data.frame(
  pred_wave_year = seq(min(site_data$pred_wave_year), max(site_data$pred_wave_year), length.out = 300)
)
grid$year10 <- (grid$pred_wave_year - mean(site_data$pred_wave_year)) / 10
grid$Site <- factor(levels(sample_data$Site)[1], levels = levels(sample_data$Site))
grid$anti_fit <- 100 * predict(bb_mean_alt, newdata = grid, type = "response")
grid$phi_fit <- predict(bb_disp_alt, newdata = grid, type = "disp")
grid$rho_fit <- 1 / (grid$phi_fit + 1)
grid$d1_fit <- exp(predict(d1_trend$fit, newdata = data.frame(x = grid$year10)))
grid$aitch_fit <- exp(predict(aitch_trend$fit, newdata = data.frame(x = grid$year10)))
grid$anti_aitch_fit <- exp(predict(anti_aitch_trend$fit,
                                   newdata = data.frame(x = grid$year10)))
grid$robust_aitch_fit <- exp(predict(robust_aitch_trend$fit,
                                     newdata = data.frame(x = grid$year10)))
grid$anti_robust_aitch_fit <- exp(predict(anti_robust_aitch_trend$fit,
                                          newdata = data.frame(x = grid$year10)))

pal <- c("Lost Iguana" = "#B44682", "Veragua" = "#166AA5",
         "Altos de Campana" = "#18835C", "Soberanía" = "#D07A00")
base_theme <- theme_classic(base_size = 8, base_family = "sans") +
  theme(axis.line = element_line(linewidth = 0.4), axis.ticks = element_line(linewidth = 0.4),
        plot.title = element_text(size = 9, face = "bold"),
        plot.subtitle = element_text(size = 6.8), legend.position = "none")

make_panel <- function(y, error, fit_y, title, ylab, tag, role, row_index) {
  p_adj <- tests$Holm_p[row_index]
  ggplot(site_data, aes(pred_wave_year, .data[[y]], colour = Site)) +
    geom_errorbar(aes(ymin = pmax(0, .data[[y]] - .data[[error]]),
                      ymax = .data[[y]] + .data[[error]]),
                  width = 0.25, linewidth = 0.4) +
    geom_line(data = grid, aes(pred_wave_year, .data[[fit_y]]), inherit.aes = FALSE,
              colour = "#4D4D4D", linewidth = 0.85) +
    geom_point(aes(size = Sample_n), alpha = 0.95) +
    geom_text_repel(aes(label = Site), size = 2.05, seed = 20260815,
                    min.segment.length = 0, box.padding = 0.25) +
    scale_colour_manual(values = pal) + scale_size_continuous(range = c(2.7, 4.8)) +
    scale_x_continuous(breaks = seq(1985, 2010, 5)) +
    labs(title = paste0(tag, "  ", title),
         subtitle = sprintf("%s; exact p = %.3f; Holm p = %.3f",
                            role, tests$Exact_p[row_index], p_adj),
         x = "Estimated Bd arrival year", y = ylab) + base_theme
}

p_a <- make_panel("beta_binomial_mean_percent", "beta_binomial_mean_SE_percent",
                  "anti_fit", "Anti-Bd relative abundance",
                  "Beta-binomial estimate (%)", "a", "Primary endpoint", 1) +
  expand_limits(y = 0)
p_b <- make_panel("d1_mean", "d1_se", "d1_fit", "Anti-Bd Hill D1 diversity",
                  "Effective number of common genera", "b", "Primary endpoint", 2) +
  expand_limits(y = 0)
p_c <- make_panel("beta_binomial_overdispersion_rho",
                  "beta_binomial_overdispersion_rho_SE", "rho_fit",
                  "Anti-Bd relative-abundance variation",
                  expression("Extra-binomial variation " * rho), "c", "Secondary endpoint", 3) +
  expand_limits(y = 0)
p_disp_whole <- make_panel("dispersion_mean", "dispersion_se", "aitch_fit",
                           "Whole-community heterogeneity", "Mean Aitchison distance", "a",
                           "Secondary endpoint", 4)
p_disp_anti <- make_panel("anti_dispersion_mean", "anti_dispersion_se", "anti_aitch_fit",
                          "Anti-Bd community heterogeneity", "Mean Aitchison distance", "b",
                          "Secondary endpoint", 5)
p_disp_whole_robust <- make_panel(
  "robust_dispersion_mean", "robust_dispersion_se", "robust_aitch_fit",
  "Whole-community heterogeneity", "Mean robust Aitchison distance", "a",
  "Sensitivity endpoint", 6
)
p_disp_anti_robust <- make_panel(
  "anti_robust_dispersion_mean", "anti_robust_dispersion_se",
  "anti_robust_aitch_fit", "Anti-Bd community heterogeneity",
  "Mean robust Aitchison distance", "b", "Sensitivity endpoint", 7
)

figure_primary <- (p_a | p_b) +
  plot_annotation(
    title = "Bd-arrival trends in anti-Bd abundance and diversity",
    subtitle = "One prespecified temporal trend parameter per endpoint; four independent localities.",
    caption = paste0("Lines show hypothesis-aligned beta-binomial or log-linear trends, not smoothing splines. ",
                     "Exact p values use all 24 permutations of Bd years among localities; error bars are s.e.m."),
    theme = theme(plot.title = element_text(size = 11, face = "bold"),
                  plot.subtitle = element_text(size = 7.2),
                  plot.caption = element_text(size = 6.1, hjust = 0))
  )

figure_aitchison <- (p_disp_whole | p_disp_anti) +
  plot_annotation(
    title = "Bd-arrival trends in within-locality Aitchison dispersion",
    subtitle = "Whole and anti-Bd bacterial communities; four independent localities.",
    caption = paste0("Lines show prespecified log-linear trends, not smoothing splines. ",
                     "Exact p values use all 24 permutations of Bd years among localities; error bars are s.e.m."),
    theme = theme(plot.title = element_text(size = 11, face = "bold"),
                  plot.subtitle = element_text(size = 7.2),
                  plot.caption = element_text(size = 6.1, hjust = 0))
  )

figure_robust_aitchison <- (p_disp_whole_robust | p_disp_anti_robust) +
  plot_annotation(
    title = "Bd-arrival trends in within-locality robust Aitchison dispersion",
    subtitle = paste0("Whole and anti-Bd bacterial communities; zero-aware sensitivity analysis; ",
                      "four independent localities."),
    caption = paste0("Lines show the same prespecified log-linear trend used for CLR-Aitchison. ",
                     "Exact p values use all 24 permutations of Bd years among localities; ",
                     "error bars are s.e.m."),
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
  grDevices::pdf(paste0(stem, ".pdf"), width = w, height = h,
                 family = "Helvetica", useDingbats = FALSE)
  print(plot); dev.off()
}
save_plot(figure_primary, file.path(figure_dir, "figure_primary_ab_bd_arrival_trends"),
          width_mm = 183, height_mm = 92)
save_plot(figure_aitchison, file.path(figure_dir, "figure_aitchison_ab_whole_vs_anti_bd"),
          width_mm = 183, height_mm = 92)
save_plot(figure_robust_aitchison,
          file.path(figure_dir, "figure_robust_aitchison_ab_whole_vs_anti_bd"),
          width_mm = 183, height_mm = 92)

cat("Hypothesis tests\n")
print(tests)
