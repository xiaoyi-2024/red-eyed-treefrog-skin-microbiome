#!/usr/bin/env Rscript

# Bd-arrival-year trends for the seven anti-Bd genera that passed the robust
# ANCOM-BC2 global locality test. Beta-binomial models retain read depth;
# exact inference permutes the four locality-level Bd years (24 assignments).

suppressPackageStartupMessages({
  library(phyloseq)
  library(glmmTMB)
  library(ggplot2)
  library(svglite)
  library(ragg)
})

set.seed(20260819)

result_dir <- "results/genus_ancombc2_aldex2/robust_anti_bd_followup"
figure_dir <- file.path(result_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

whole_file <- "data/processed/ps_clean_47_rebuilt_from_01_filtered.rds"
anti_file <- "data/processed/blast_97id_90cov/ps_anti_bd_47_rebuilt_97id_90cov.rds"
test_file <- file.path(result_dir, "robust_genera_global_ancombc2.csv")
wave_file <- "results/wave_exploratory/site_metrics_single_wave.csv"
stopifnot(file.exists(whole_file), file.exists(anti_file),
          file.exists(test_file), file.exists(wave_file))

site_key <- c(
  "Lost Iguana Resort, La Fortuna, San Carlos" = "Lost Iguana",
  "Veragua Rainforest, Las Brisas de Veragua" = "Veragua",
  "PN Altos de Campana" = "Altos de Campana",
  "PN Soberanía, Las Cumbres" = "Soberanía"
)
site_levels <- c("Lost Iguana", "Veragua", "Altos de Campana", "Soberanía")
site_palette <- c(
  "Lost Iguana" = "#B44682", "Veragua" = "#216DA3",
  "Altos de Campana" = "#21865F", "Soberanía" = "#D17A00"
)

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

ps_whole <- readRDS(whole_file)
ps_anti <- readRDS(anti_file)
whole_asv <- otu_samples_rows(ps_whole)
anti_asv <- otu_samples_rows(ps_anti)[rownames(whole_asv), , drop = FALSE]
anti_genus <- aggregate_genus(ps_anti, anti_asv)
md <- data.frame(sample_data(ps_whole), check.names = FALSE,
                 stringsAsFactors = FALSE)[rownames(whole_asv), , drop = FALSE]
site <- factor(unname(site_key[md$Locality]), levels = site_levels)
total_reads <- rowSums(whole_asv)
stopifnot(nrow(whole_asv) == 47L, !anyNA(site), all(total_reads > 0))

robust <- read.csv(test_file, stringsAsFactors = FALSE, check.names = FALSE)
genera <- robust$taxon[robust$diff_robust_abn %in% TRUE & robust$q_val < 0.05]
stopifnot(length(genera) == 7L, all(genera %in% colnames(anti_genus)))

wave <- read.csv(wave_file, stringsAsFactors = FALSE, check.names = FALSE)
site_year <- setNames(wave$pred_wave_year, wave$Site)
year <- unname(site_year[as.character(site)])
year_center <- mean(unname(site_year[site_levels]))

permutations_all <- function(x) {
  if (length(x) == 1L) return(matrix(x, nrow = 1L))
  do.call(rbind, lapply(seq_along(x), function(i) {
    cbind(x[i], permutations_all(x[-i]))
  }))
}
perm_years <- permutations_all(unname(site_year[site_levels]))

fit_bb <- function(success, assigned_year) {
  d <- data.frame(Success = success, Failure = total_reads - success,
                  Year10 = (assigned_year - year_center) / 10)
  suppressWarnings(glmmTMB(cbind(Success, Failure) ~ Year10,
                           family = betabinomial(link = "logit"), data = d))
}

model_rows <- list()
plot_rows <- list()
pred_rows <- list()

for (g in genera) {
  success <- as.numeric(anti_genus[, g])
  observed <- fit_bb(success, year)
  beta <- unname(fixef(observed)$cond["Year10"])
  se <- sqrt(diag(vcov(observed)$cond))["Year10"]
  perm_beta <- apply(perm_years, 1, function(py) {
    yp <- py[match(as.character(site), site_levels)]
    fit <- try(fit_bb(success, yp), silent = TRUE)
    if (inherits(fit, "try-error")) return(NA_real_)
    unname(fixef(fit)$cond["Year10"])
  })
  exact_p <- mean(abs(perm_beta) >= abs(beta) - 1e-12, na.rm = TRUE)

  model_rows[[g]] <- data.frame(
    Genus = g, Beta_log_odds_per_10y = beta, SE = se,
    Odds_ratio_per_10y = exp(beta), Exact_p = exact_p,
    Successful_permutations = sum(is.finite(perm_beta)),
    Localities = 4L, Frog_samples = length(success), stringsAsFactors = FALSE
  )

  sample_df <- data.frame(
    Genus = g, Sample_ID = rownames(whole_asv), Site = site,
    Bd_arrival_year = year,
    Relative_abundance_percent = 100 * success / total_reads
  )
  site_df <- aggregate(Relative_abundance_percent ~ Genus + Site + Bd_arrival_year,
                       sample_df, function(x) c(Mean = mean(x),
                                                SE = sd(x) / sqrt(length(x)),
                                                N = length(x)))
  site_df <- cbind(site_df[, 1:3], as.data.frame(site_df[, 4]))
  names(site_df)[4:6] <- c("Mean_percent", "SE_percent", "N")
  plot_rows[[g]] <- site_df

  grid_year <- seq(min(year), max(year), length.out = 160)
  nd <- data.frame(Year10 = (grid_year - year_center) / 10)
  pr <- predict(observed, newdata = nd, type = "link", se.fit = TRUE)
  pred_rows[[g]] <- data.frame(
    Genus = g, Bd_arrival_year = grid_year,
    Fit_percent = 100 * plogis(pr$fit),
    Low_percent = 100 * plogis(pr$fit - 1.96 * pr$se.fit),
    High_percent = 100 * plogis(pr$fit + 1.96 * pr$se.fit)
  )
}

models <- do.call(rbind, model_rows)
models$Holm_p <- p.adjust(models$Exact_p, method = "holm")
models$Direction <- ifelse(models$Beta_log_odds_per_10y > 0, "Increase", "Decrease")
models <- models[order(models$Holm_p, models$Exact_p), ]
write.csv(models, file.path(result_dir, "robust_genera_bd_year_beta_binomial_tests.csv"),
          row.names = FALSE)

site_summary <- do.call(rbind, plot_rows)
predictions <- do.call(rbind, pred_rows)
write.csv(site_summary, file.path(result_dir, "robust_genera_bd_year_site_summary.csv"),
          row.names = FALSE)
write.csv(predictions, file.path(result_dir, "robust_genera_bd_year_predictions.csv"),
          row.names = FALSE)

display_genus <- models$Genus
display_genus[display_genus == "Burkholderia-Caballeronia-Paraburkholderia"] <-
  "Burkholderia complex"
## Keep facet strips short enough for a four-column journal layout. Full
## effect sizes and multiplicity-adjusted tests are reported in the companion
## CSV and the Results text rather than squeezed into the strip labels.
label_lookup <- setNames(display_genus, models$Genus)
site_summary$Panel <- factor(label_lookup[site_summary$Genus],
                             levels = label_lookup[models$Genus])
predictions$Panel <- factor(label_lookup[predictions$Genus],
                            levels = label_lookup[models$Genus])

p <- ggplot() +
  geom_ribbon(data = predictions,
              aes(Bd_arrival_year, ymin = Low_percent, ymax = High_percent),
              fill = "#8E8E8E", alpha = 0.18) +
  geom_line(data = predictions, aes(Bd_arrival_year, Fit_percent),
            colour = "#333333", linewidth = 0.75) +
  geom_errorbar(data = site_summary,
                aes(Bd_arrival_year,
                    ymin = pmax(0, Mean_percent - SE_percent),
                    ymax = Mean_percent + SE_percent, colour = Site),
                width = 0.22, linewidth = 0.55) +
  geom_point(data = site_summary,
             aes(Bd_arrival_year, Mean_percent, fill = Site, size = N),
             shape = 21, colour = "white", stroke = 0.55) +
  facet_wrap(~Panel, scales = "free_y", ncol = 4) +
  scale_colour_manual(values = site_palette) +
  scale_fill_manual(values = site_palette) +
  scale_size_continuous(range = c(2.2, 3.4), breaks = c(9, 10, 18)) +
  scale_x_continuous(breaks = c(1983, 1990, 2000, 2009)) +
  labs(
    title = "Bd-arrival-year trends in seven locality-associated anti-Bd genera",
    subtitle = paste0("Beta-binomial fits retain read depth; exact inference permutes ",
                      "Bd years at the locality level."),
    x = "Estimated Bd-arrival year", y = "Genus reads (% of total bacterial reads)",
    fill = "Locality", colour = "Locality", size = "Frogs",
    caption = paste0(
      "Points are locality means; error bars are s.e.m. Curves and 95% model intervals are descriptive.\n",
      "Holm correction covers the seven genera; locality and Bd year are completely linked (n = 4 independent localities)."
    )
  ) +
  theme_classic(base_size = 7.5, base_family = "sans") +
  theme(
    axis.line = element_line(linewidth = 0.4),
    strip.background = element_rect(fill = "grey95", colour = "grey75", linewidth = 0.3),
    strip.text = element_text(size = 6.4, face = "bold", lineheight = 0.95),
    plot.title = element_text(size = 10.5, face = "bold"),
    plot.subtitle = element_text(size = 7),
    plot.caption = element_text(size = 5.8, hjust = 0),
    legend.position = "top", legend.title = element_text(size = 6.4),
    legend.text = element_text(size = 6.2),
    panel.spacing = grid::unit(5, "mm")
  )

save_pub <- function(plot, stem, width_mm = 183, height_mm = 132, dpi = 600) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  ggsave(paste0(stem, ".png"), plot, width = w, height = h, dpi = dpi,
         bg = "white")
  svglite(paste0(stem, ".svg"), width = w, height = h); print(plot); dev.off()
  pdf(paste0(stem, ".pdf"), width = w, height = h, family = "Helvetica",
      useDingbats = FALSE); print(plot); dev.off()
  tryCatch({
    agg_tiff(paste0(stem, ".tiff"), width = w, height = h, units = "in",
             res = dpi, background = "white")
    print(plot); dev.off()
  }, error = function(e) {
    if (dev.cur() > 1L) try(dev.off(), silent = TRUE)
    warning("TIFF export failed after PNG/PDF/SVG succeeded: ", conditionMessage(e))
  })
}

save_pub(p, file.path(figure_dir, "figure_seven_robust_genera_bd_year_fits"))

cat("Completed seven robust anti-Bd genus Bd-year models.\n")
print(models)
