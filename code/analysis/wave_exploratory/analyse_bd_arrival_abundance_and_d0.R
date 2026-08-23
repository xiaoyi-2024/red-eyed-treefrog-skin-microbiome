suppressPackageStartupMessages({
  library(glmmTMB)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(svglite)
  library(ragg)
})

set.seed(20260819)

# Run from submission2/. All inputs are previously generated, reproducible
# analysis tables; no absolute paths are used.
result_dir <- "results/wave_exploratory"
figure_dir <- file.path(result_dir, "figures")
hill_file <- "results/core_locality_abundance_diversity/sample_hill_D0_D1_D2.csv"
sample_file <- file.path(result_dir, "sample_three_metrics.csv")
bb_file <- file.path(result_dir, "beta_binomial_site_estimates.csv")
year_file <- file.path(result_dir, "bd_wave_spatial_model_site_predictions.csv")
stopifnot(file.exists(hill_file), file.exists(sample_file), file.exists(bb_file),
          file.exists(year_file))
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

hill <- read.csv(hill_file, stringsAsFactors = FALSE, check.names = FALSE)
samples <- read.csv(sample_file, stringsAsFactors = FALSE, check.names = FALSE)
bb_site <- read.csv(bb_file, stringsAsFactors = FALSE, check.names = FALSE)
years <- read.csv(year_file, stringsAsFactors = FALSE, check.names = FALSE)
names(years)[names(years) == "Site_short"] <- "Site"
names(years)[names(years) == "Predicted_DOD_year"] <- "bd_year"

site_order <- c("Lost Iguana", "Veragua", "Altos de Campana", "Soberanía")
site_colours <- c(
  "Lost Iguana" = "#B33B7D", "Veragua" = "#176AA6",
  "Altos de Campana" = "#13865E", "Soberanía" = "#D17800"
)

# D0 is genus richness (the number of detected genera), not abundance. Here it
# is calculated per frog and then summarized at the independent locality level.
d0 <- hill[hill$Order_q == 0, c("Sample_ID", "Hill_D", "Community", "Site")]
d0$Community <- sub(" community$", "", d0$Community)
d0$Community[d0$Community == "Anti-Bd"] <- "Anti-Bd"
d0$Community[d0$Community == "Whole"] <- "Whole community"

summarise_d0 <- function(z) {
  data.frame(
    Site = z$Site[1], Community = z$Community[1], n = nrow(z),
    D0_mean = mean(z$Hill_D), D0_SD = sd(z$Hill_D),
    D0_SE = sd(z$Hill_D) / sqrt(nrow(z)), stringsAsFactors = FALSE
  )
}
d0_site <- do.call(rbind, lapply(split(d0, interaction(d0$Site, d0$Community,
                                                       drop = TRUE)), summarise_d0))
rownames(d0_site) <- NULL
d0_site <- merge(d0_site, years[, c("Site", "bd_year")], by = "Site", all.x = TRUE)
d0_site$Site <- factor(d0_site$Site, levels = site_order)
d0_site <- d0_site[order(d0_site$Community, d0_site$Site), ]
stopifnot(nrow(d0_site) == 8L, !anyNA(d0_site$bd_year))
write.csv(d0_site, file.path(result_dir, "bd_arrival_D0_site_summary.csv"), row.names = FALSE)

permutations <- function(x) {
  if (length(x) == 1L) return(matrix(x, nrow = 1L))
  do.call(rbind, lapply(seq_along(x), function(i) cbind(x[i], permutations(x[-i]))))
}

# A single prespecified log-linear slope is fitted to the four locality means.
# Inverse-SE-squared weighting gives more weight to more precisely estimated
# locality means. The exact P value enumerates all 4! assignments of years.
fit_d0 <- function(z) {
  z <- z[order(z$Site), ]
  x <- (z$bd_year - mean(z$bd_year)) / 10
  w <- 1 / z$D0_SE^2
  null <- lm(log(D0_mean) ~ 1, weights = w, data = z)
  alt <- lm(log(D0_mean) ~ x, weights = w, data = z)
  obs <- max(0, deviance(null) - deviance(alt))
  pm <- permutations(z$bd_year)
  stat <- beta <- numeric(nrow(pm))
  for (i in seq_len(nrow(pm))) {
    xp <- (pm[i, ] - mean(z$bd_year)) / 10
    fp <- lm(log(z$D0_mean) ~ xp, weights = w)
    stat[i] <- max(0, deviance(null) - deviance(fp))
    beta[i] <- coef(fp)["xp"]
  }
  list(data = z, fit = alt, slope = unname(coef(alt)["x"]),
       ratio = exp(unname(coef(alt)["x"])), statistic = obs,
       exact_p = mean(stat >= obs - 1e-12), perm_stat = stat, perm_beta = beta)
}

d0_fits <- lapply(split(d0_site, d0_site$Community), fit_d0)
d0_tests <- do.call(rbind, lapply(names(d0_fits), function(nm) {
  x <- d0_fits[[nm]]
  data.frame(
    Metric = paste(nm, "genus Hill D0"),
    Effect = x$slope,
    Effect_definition = "log D0 ratio per 10-year later estimated Bd arrival",
    Ratio_per_10_years = x$ratio,
    Percent_change_per_10_years = 100 * (x$ratio - 1),
    Exact_p = x$exact_p,
    Localities = 4L,
    stringsAsFactors = FALSE
  )
}))
d0_tests$Holm_p <- p.adjust(d0_tests$Exact_p, method = "holm")
write.csv(d0_tests, file.path(result_dir, "bd_arrival_D0_exact_trend_tests.csv"), row.names = FALSE)

d0_perm <- do.call(rbind, lapply(names(d0_fits), function(nm) {
  x <- d0_fits[[nm]]
  data.frame(Metric = paste(nm, "genus Hill D0"),
             Permutation = seq_along(x$perm_stat), Statistic = x$perm_stat,
             Coefficient_per_10_years = x$perm_beta)
}))
write.csv(d0_perm, file.path(result_dir, "bd_arrival_D0_exact_permutations.csv"), row.names = FALSE)

# Summarize and test the complete Hill profile (q = 0, 1, 2) so the figure can
# distinguish changes in richness from changes in common and dominant genera.
summarise_hill <- function(z) {
  data.frame(
    Site = z$Site[1], Community = z$Community[1], Order_q = z$Order_q[1],
    n = nrow(z), Hill_mean = mean(z$Hill_D), Hill_SD = sd(z$Hill_D),
    Hill_SE = sd(z$Hill_D) / sqrt(nrow(z)), stringsAsFactors = FALSE
  )
}
hill$Community <- sub(" community$", "", hill$Community)
hill$Community[hill$Community == "Whole"] <- "Whole community"
hill_profile_site <- do.call(rbind, lapply(
  split(hill, interaction(hill$Site, hill$Community, hill$Order_q, drop = TRUE)),
  summarise_hill
))
rownames(hill_profile_site) <- NULL
hill_profile_site <- merge(hill_profile_site, years[, c("Site", "bd_year")],
                           by = "Site", all.x = TRUE)
hill_profile_site$Site <- factor(hill_profile_site$Site, levels = site_order)
hill_profile_site <- hill_profile_site[order(hill_profile_site$Community,
                                             hill_profile_site$Order_q,
                                             hill_profile_site$Site), ]
stopifnot(nrow(hill_profile_site) == 24L, !anyNA(hill_profile_site$bd_year))
write.csv(hill_profile_site,
          file.path(result_dir, "bd_arrival_Hill_D0_D1_D2_site_summary.csv"),
          row.names = FALSE)

fit_hill_profile <- function(z) {
  x <- z
  x$D0_mean <- x$Hill_mean
  x$D0_SE <- x$Hill_SE
  fit_d0(x)
}
hill_split <- split(hill_profile_site,
                    interaction(hill_profile_site$Community,
                                hill_profile_site$Order_q, drop = TRUE))
hill_fits <- lapply(hill_split, fit_hill_profile)
hill_tests <- do.call(rbind, lapply(names(hill_fits), function(nm) {
  z <- hill_split[[nm]]
  x <- hill_fits[[nm]]
  data.frame(
    Community = z$Community[1], Order_q = z$Order_q[1],
    Metric = paste0("Hill D", z$Order_q[1]), Effect = x$slope,
    Ratio_per_10_years = x$ratio,
    Percent_change_per_10_years = 100 * (x$ratio - 1),
    Exact_p = x$exact_p, Localities = 4L, stringsAsFactors = FALSE
  )
}))
hill_tests$Holm_p_all_six <- p.adjust(hill_tests$Exact_p, method = "holm")
write.csv(hill_tests,
          file.path(result_dir, "bd_arrival_Hill_D0_D1_D2_exact_trend_tests.csv"),
          row.names = FALSE)

# Retain the previously specified beta-binomial model for anti-Bd relative
# abundance because it uses anti-Bd and non-anti-Bd read counts from every frog.
samples$Site <- factor(samples$Site, levels = site_order)
samples$bd_year <- years$bd_year[match(samples$Site, years$Site)]
samples$year10 <- (samples$bd_year - mean(years$bd_year)) / 10
bb_null <- glmmTMB(cbind(anti_reads, non_anti_reads) ~ 1,
                   dispformula = ~ Site, family = betabinomial(link = "logit"),
                   data = samples)
bb_alt <- update(bb_null, . ~ year10)
bb_obs <- max(0, 2 * (as.numeric(logLik(bb_alt)) - as.numeric(logLik(bb_null))))
pm <- permutations(years$bd_year[match(site_order, years$Site)])
bb_stat <- bb_beta <- numeric(nrow(pm))
for (i in seq_len(nrow(pm))) {
  z <- samples
  map <- setNames(pm[i, ], site_order)
  z$year10 <- (map[as.character(z$Site)] - mean(years$bd_year)) / 10
  fp <- update(bb_null, . ~ year10, data = z)
  bb_stat[i] <- max(0, 2 * (as.numeric(logLik(fp)) - as.numeric(logLik(bb_null))))
  bb_beta[i] <- fixef(fp)$cond["year10"]
}
bb_coef <- fixef(bb_alt)$cond["year10"]
bb_test <- data.frame(
  Metric = "Anti-Bd relative abundance",
  Effect = unname(bb_coef), Effect_definition = "log odds ratio per 10-year later estimated Bd arrival",
  Ratio_per_10_years = exp(unname(bb_coef)),
  Percent_change_per_10_years = NA_real_,
  Exact_p = mean(bb_stat >= bb_obs - 1e-10), Holm_p = NA_real_, Localities = 4L
)
write.csv(bb_test, file.path(result_dir, "bd_arrival_anti_bd_abundance_exact_trend_test.csv"),
          row.names = FALSE)

# Publication figures ------------------------------------------------------
theme_paper <- theme_classic(base_size = 10, base_family = "Helvetica") +
  theme(plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8.5, colour = "grey25"),
        axis.title = element_text(face = "bold"),
        legend.position = "none", plot.margin = margin(6, 8, 6, 6))

save_all <- function(p, stem, width, height) {
  ggsave(paste0(stem, ".pdf"), p, width = width, height = height, units = "in",
         device = grDevices::pdf)
  ggsave(paste0(stem, ".svg"), p, width = width, height = height, units = "in",
         device = svglite::svglite)
  ggsave(paste0(stem, ".png"), p, width = width, height = height, units = "in",
         dpi = 400, device = ragg::agg_png)
  ggsave(paste0(stem, ".tiff"), p, width = width, height = height, units = "in",
         dpi = 600, compression = "lzw", device = ragg::agg_tiff)
}

ab_plot <- merge(bb_site, years[, c("Site", "bd_year")], by = "Site")
ab_plot$Site <- factor(ab_plot$Site, levels = site_order)
grid <- data.frame(bd_year = seq(min(years$bd_year), max(years$bd_year), length.out = 250))
grid$year10 <- (grid$bd_year - mean(years$bd_year)) / 10
pred <- predict(bb_alt, newdata = transform(grid, Site = factor("Lost Iguana",
                                                                levels = site_order)),
                type = "response", se.fit = TRUE)
grid$fit <- 100 * pred$fit
grid$low <- 100 * pmax(0, pred$fit - 1.96 * pred$se.fit)
grid$high <- 100 * pmin(1, pred$fit + 1.96 * pred$se.fit)

p_ab <- ggplot(ab_plot, aes(bd_year, beta_binomial_mean_percent)) +
  geom_ribbon(data = grid, aes(x = bd_year, ymin = low, ymax = high), inherit.aes = FALSE,
              fill = "#8E8E8E", alpha = 0.18) +
  geom_line(data = grid, aes(x = bd_year, y = fit), inherit.aes = FALSE,
            colour = "#3B3B3B", linewidth = 0.9) +
  geom_errorbar(aes(ymin = pmax(0, beta_binomial_mean_percent -
                                  beta_binomial_mean_SE_percent),
                    ymax = beta_binomial_mean_percent + beta_binomial_mean_SE_percent,
                    colour = Site), width = 0.25, linewidth = 0.65) +
  geom_point(aes(fill = Site, colour = Site), shape = 21, size = 3.4, stroke = 0.8) +
  ggrepel::geom_text_repel(aes(label = Site, colour = Site), size = 3,
                           min.segment.length = 0, seed = 12, show.legend = FALSE) +
  scale_colour_manual(values = site_colours) + scale_fill_manual(values = site_colours) +
  labs(title = "Anti-Bd relative abundance and estimated Bd-arrival year",
       subtitle = sprintf("Beta-binomial trend: OR per 10 years = %.2f; exact P = %.3f",
                          exp(bb_coef), bb_test$Exact_p),
       x = "Estimated Bd-arrival year", y = "Estimated anti-Bd relative abundance (%)") +
  theme_paper
save_all(p_ab, file.path(figure_dir, "figure_anti_bd_abundance_bd_arrival"), 6.4, 4.2)

q_labels <- c(`0` = "D0 · richness", `1` = "D1 · common genera",
              `2` = "D2 · dominant genera")
q_colours <- c(`0` = "#2B6F9C", `1` = "#D07A25", `2` = "#7A4E9D")
q_linetypes <- c(`0` = "solid", `1` = "longdash", `2` = "dotted")

make_profile_panel <- function(community, panel_title) {
  z <- hill_profile_site[hill_profile_site$Community == community, ]
  z$Order_q <- factor(z$Order_q, levels = 0:2)
  curves <- do.call(rbind, lapply(0:2, function(q) {
    key <- names(hill_split)[vapply(hill_split, function(x) {
      x$Community[1] == community && x$Order_q[1] == q
    }, logical(1))]
    stopifnot(length(key) == 1L)
    gx <- data.frame(bd_year = seq(min(z$bd_year), max(z$bd_year), length.out = 250))
    gx$x <- (gx$bd_year - mean(z$bd_year)) / 10
    gx$fit <- exp(predict(hill_fits[[key]]$fit, newdata = gx))
    gx$Order_q <- factor(q, levels = 0:2)
    gx
  }))
  stat <- hill_tests[hill_tests$Community == community, ]
  stat <- stat[order(stat$Order_q), ]
  stat_label <- paste0("Exact P: D0 = ", sprintf("%.3f", stat$Exact_p[1]),
                       "; D1 = ", sprintf("%.3f", stat$Exact_p[2]),
                       "; D2 = ", sprintf("%.3f", stat$Exact_p[3]))

  ggplot(z, aes(bd_year, Hill_mean, colour = Order_q, linetype = Order_q)) +
    geom_line(data = curves, aes(y = fit), linewidth = 0.9) +
    geom_errorbar(aes(ymin = pmax(0.01, Hill_mean - Hill_SE),
                      ymax = Hill_mean + Hill_SE), width = 0.24,
                  linewidth = 0.55, alpha = 0.8) +
    geom_point(aes(shape = Order_q), size = 2.8, stroke = 0.8, fill = "white") +
    scale_y_log10() +
    scale_colour_manual(values = q_colours, labels = q_labels) +
    scale_linetype_manual(values = q_linetypes, labels = q_labels) +
    scale_shape_manual(values = c(`0` = 16, `1` = 17, `2` = 15), labels = q_labels) +
    scale_x_continuous(breaks = round(years$bd_year[match(site_order, years$Site)]),
                       labels = c("1986", "1993", "2008", "2010")) +
    labs(title = panel_title, subtitle = stat_label,
         x = "Estimated Bd-arrival year",
         y = "Effective number of genera (log scale)",
         colour = NULL, linetype = NULL, shape = NULL) +
    theme_paper +
    theme(plot.title = element_text(size = 10), plot.subtitle = element_text(size = 7.5),
          legend.position = "top", legend.direction = "horizontal",
          legend.text = element_text(size = 7.3),
          axis.title = element_text(size = 8.5), axis.text = element_text(size = 7.5))
}

p_whole_profile <- make_profile_panel("Whole community", "a  Whole community")
p_anti_profile <- make_profile_panel("Anti-Bd", "b  Anti-Bd community")
p_ab_profile <- p_ab +
  scale_x_continuous(breaks = round(years$bd_year[match(site_order, years$Site)]),
                     labels = c("1986", "1993", "2008", "2010")) +
  labs(title = "c  Anti-Bd relative abundance",
       subtitle = sprintf("Beta-binomial OR per 10 years = %.2f; exact P = %.3f",
                          exp(bb_coef), bb_test$Exact_p),
       x = "Estimated Bd-arrival year",
       y = "Estimated relative abundance (%)") +
  guides(colour = "none", fill = "none") +
  theme(plot.title = element_text(size = 10), plot.subtitle = element_text(size = 7.5),
        axis.title = element_text(size = 8.5), axis.text = element_text(size = 7.5))

p_hill <- (p_whole_profile | p_anti_profile | p_ab_profile) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Microbiome diversity, anti-Bd abundance and estimated Bd-arrival year",
    subtitle = paste0("Panels a-b overlay D0 richness, D1 common genera and D2 dominant genera; ",
                      "the abundance panel uses beta-binomial locality estimates."),
    caption = paste0("Lines are prespecified inverse-variance weighted log-linear fits. ",
                     "Exact P values enumerate all 24 assignments of years to four localities; ",
                     "error bars are s.e.m. for Hill diversity and model s.e. for relative abundance.")
  ) & theme(legend.position = "top")
save_all(p_hill, file.path(figure_dir, "figure_whole_antibd_D0_bd_arrival"), 14.2, 5.4)
save_all(p_hill, file.path(figure_dir, "figure_whole_antibd_hill_D0_D1_D2_bd_arrival"),
         14.2, 5.4)

message("Created anti-Bd abundance and whole/anti-Bd D0 analyses and figures.")
