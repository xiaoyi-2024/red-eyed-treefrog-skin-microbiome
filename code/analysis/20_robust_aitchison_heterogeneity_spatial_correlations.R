#!/usr/bin/env Rscript

# Locality-level correlations between coordinates/elevation and within-locality
# robust Aitchison heterogeneity. The independent unit is locality (n = 4).

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(svglite)
  library(ragg)
})

result_dir <- "results/robust_aitchison_heterogeneity_spatial"
figure_dir <- file.path(result_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

input_file <- "results/spatial_elevation_abundance_diversity/site_spatial_metrics.csv"
z <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(nrow(z) == 4L)

long <- rbind(
  data.frame(Site = z$Site, Longitude = z$Longitude, Latitude = z$Latitude,
             Elevation_m = z$Elevation_m, Community = "Whole community",
             Mean_robust_heterogeneity = z$robust_dispersion_mean,
             SE_robust_heterogeneity = z$robust_dispersion_se),
  data.frame(Site = z$Site, Longitude = z$Longitude, Latitude = z$Latitude,
             Elevation_m = z$Elevation_m, Community = "Anti-Bd",
             Mean_robust_heterogeneity = z$anti_robust_dispersion_mean,
             SE_robust_heterogeneity = z$anti_robust_dispersion_se)
)
long$Community <- factor(long$Community, levels = c("Whole community", "Anti-Bd"))
write.csv(long, file.path(result_dir, "locality_robust_heterogeneity_spatial_data.csv"),
          row.names = FALSE)

n <- nrow(z)
perm_grid <- as.matrix(expand.grid(rep(list(seq_len(n)), n)))
perms <- perm_grid[apply(perm_grid, 1, function(v) length(unique(v)) == n), , drop = FALSE]
stopifnot(nrow(perms) == factorial(n))

exact_spearman <- function(x, y) {
  obs <- suppressWarnings(cor(x, y, method = "spearman"))
  null <- apply(perms, 1, function(idx)
    suppressWarnings(cor(x, y[idx], method = "spearman")))
  c(rho = obs, Exact_p = mean(abs(null) >= abs(obs) - 1e-12))
}

predictors <- c(Longitude = "Longitude", Latitude = "Latitude", Elevation = "Elevation_m")
tests <- do.call(rbind, lapply(levels(long$Community), function(comm) {
  d <- long[long$Community == comm, ]
  do.call(rbind, lapply(names(predictors), function(pred) {
    ans <- exact_spearman(d[[predictors[pred]]], d$Mean_robust_heterogeneity)
    fit <- lm(Mean_robust_heterogeneity ~ d[[predictors[pred]]], data = d)
    data.frame(
      Community = comm, Predictor = pred,
      Spearman_rho = ans["rho"], Exact_p = ans["Exact_p"],
      Descriptive_linear_R2 = summary(fit)$r.squared,
      Independent_localities = n, Permutations = nrow(perms),
      stringsAsFactors = FALSE
    )
  }))
}))
tests$Holm_p <- p.adjust(tests$Exact_p, method = "holm")
tests$Community <- factor(tests$Community, levels = levels(long$Community))
tests$Predictor <- factor(tests$Predictor, levels = names(predictors))
tests <- tests[order(tests$Community, tests$Predictor), ]
write.csv(tests, file.path(result_dir, "exact_spatial_correlation_tests.csv"), row.names = FALSE)

plot_data <- rbind(
  transform(long, Predictor = "Longitude", X = Longitude),
  transform(long, Predictor = "Latitude", X = Latitude),
  transform(long, Predictor = "Elevation", X = Elevation_m)
)
plot_data$Predictor <- factor(plot_data$Predictor, levels = names(predictors))

ann <- tests
ann$Label <- sprintf("rho = %.2f; exact p = %.3f\nHolm p = %.3f",
                     ann$Spearman_rho, ann$Exact_p, ann$Holm_p)

p <- ggplot(plot_data, aes(X, Mean_robust_heterogeneity)) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.6, colour = "grey38") +
  geom_errorbar(aes(ymin = Mean_robust_heterogeneity - SE_robust_heterogeneity,
                    ymax = Mean_robust_heterogeneity + SE_robust_heterogeneity),
                width = 0, linewidth = 0.35, colour = "grey35") +
  geom_point(aes(fill = Site), shape = 21, size = 2.8, colour = "black", stroke = 0.3) +
  ggrepel::geom_text_repel(aes(label = Site), size = 2.0, max.overlaps = Inf,
                           box.padding = 0.25, min.segment.length = 0) +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = Label), inherit.aes = FALSE,
            hjust = -0.05, vjust = 1.1, size = 2.25) +
  facet_wrap(vars(Community, Predictor), ncol = 3, scales = "free") +
  scale_fill_brewer(palette = "Dark2", guide = "none") +
  labs(
    title = "Spatial correlates of within-locality robust Aitchison heterogeneity",
    subtitle = "Top: whole community; bottom: anti-Bd community; mean distance to locality centroid ± SE",
    x = NULL, y = "Mean robust Aitchison distance to locality centroid",
    caption = paste0("Each point is one independent locality; lines are descriptive. Exact two-sided Spearman tests enumerate all 24 locality-label permutations.\n",
                     "Holm correction covers six tests. Whole and anti-Bd panels use separate y scales because feature dimensions differ.")
  ) +
  theme_classic(base_size = 8) +
  theme(plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 7), strip.background = element_rect(fill = "grey94"),
        strip.text = element_text(face = "bold", size = 7),
        plot.caption = element_text(size = 5.8, hjust = 0),
        panel.spacing = grid::unit(7, "mm"))

save_plot <- function(plot, stem, width_mm = 183, height_mm = 142) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  ggsave(paste0(stem, ".png"), plot, width = w, height = h, dpi = 600, bg = "white")
  svglite::svglite(paste0(stem, ".svg"), width = w, height = h); print(plot); dev.off()
  pdf(paste0(stem, ".pdf"), width = w, height = h, useDingbats = FALSE); print(plot); dev.off()
  ragg::agg_tiff(paste0(stem, ".tiff"), width = w, height = h, units = "in",
                 res = 600, compression = "lzw", background = "white")
  print(plot); dev.off()
}
save_plot(p, file.path(figure_dir, "figure_robust_heterogeneity_vs_coordinates_elevation"))

write.csv(data.frame(
  Item = c("Heterogeneity estimator", "Independent unit", "Predictors",
           "Inference", "Multiplicity", "Error bars"),
  Definition = c("Mean sample-to-locality-centroid distance in rCLR space",
                 "Locality (n = 4)", "Longitude, latitude and WorldClim elevation tested separately",
                 "Exact two-sided Spearman correlation; all 4! label permutations",
                 "Holm correction across 2 communities x 3 predictors",
                 "SE among frog-level distances within locality; descriptive")
), file.path(result_dir, "analysis_methods.csv"), row.names = FALSE)

print(long)
print(tests)
