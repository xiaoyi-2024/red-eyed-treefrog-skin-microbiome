#!/usr/bin/env Rscript

# Exploratory locality-level associations of anti-Bd relative abundance and
# diversity with geographic position and elevation. The four localities, not
# the 47 frogs, are the independent geographic units.

suppressPackageStartupMessages({
  library(terra)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
})

set.seed(20260815)
result_dir <- "results/spatial_elevation_abundance_diversity"
figure_dir <- file.path(result_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

metrics_file <- "results/wave_exploratory/site_four_metrics.csv"
sample_file <- "data/processed/ps_clean_47_rebuilt_from_01_filtered.rds"
stopifnot(file.exists(metrics_file), file.exists(sample_file))

suppressPackageStartupMessages(library(phyloseq))
ps <- readRDS(sample_file)
md <- data.frame(sample_data(ps), check.names = FALSE, stringsAsFactors = FALSE)

site_key <- c(
  "Lost Iguana Resort, La Fortuna, San Carlos" = "Lost Iguana",
  "Veragua Rainforest, Las Brisas de Veragua" = "Veragua",
  "PN Altos de Campana" = "Altos de Campana",
  "PN Soberanía, Las Cumbres" = "Soberanía"
)
md$Site <- unname(site_key[md$Locality])
stopifnot(!anyNA(md$Site))

site_coordinates <- aggregate(cbind(Latitude = md$Lat, Longitude = md$Long),
                              list(Site = md$Site), mean)

# Extract elevation at each frog coordinate from the locally archived WorldClim
# 30 arc-second DEMs, then average within locality. CRI and PAN rasters overlap
# the country in which each sample occurs; the non-missing extraction is used.
dem_cr <- rast("data/reference/elevation/CRI_elv_msk.tif")
dem_pa <- rast("data/reference/elevation/PAN_elv_msk.tif")
pts <- vect(md[, c("Long", "Lat")], geom = c("Long", "Lat"), crs = "EPSG:4326")
e_cr <- terra::extract(dem_cr, pts)[, 2]
e_pa <- terra::extract(dem_pa, pts)[, 2]
md$Elevation_m <- ifelse(is.na(e_cr), e_pa, e_cr)
stopifnot(!anyNA(md$Elevation_m))
site_elevation <- aggregate(Elevation_m ~ Site, md, mean)

z <- read.csv(metrics_file, check.names = FALSE, stringsAsFactors = FALSE)
whole_d1_file <- "results/bd_arrival_community/whole_community_hill_d1_sample_values.csv"
stopifnot(file.exists(whole_d1_file))
whole_d1_sample <- read.csv(whole_d1_file, check.names = FALSE, stringsAsFactors = FALSE)
whole_d1_site <- aggregate(Whole_community_genus_Hill_D1 ~ Locality_internal,
                           whole_d1_sample, mean)
names(whole_d1_site) <- c("Site", "whole_d1_mean")
z <- merge(z, site_coordinates, by = "Site", all.x = TRUE, sort = FALSE)
z <- merge(z, site_elevation, by = "Site", all.x = TRUE, sort = FALSE)
z <- merge(z, whole_d1_site, by = "Site", all.x = TRUE, sort = FALSE)
stopifnot(nrow(z) == 4L, !anyNA(z))

# The beta-binomial conditional mean is a proportion, so use the logit scale.
# Hill D1 is positive and multiplicative, so use the natural-log scale.
p <- pmin(pmax(z$beta_binomial_mean_percent / 100, 1e-8), 1 - 1e-8)
z$logit_anti_bd_mean <- qlogis(p)
z$log_hill_d1 <- log(z$d1_mean)
z$log_whole_hill_d1 <- log(z$whole_d1_mean)
write.csv(z, file.path(result_dir, "site_spatial_metrics.csv"), row.names = FALSE)

haversine_km <- function(lon1, lat1, lon2, lat2) {
  rad <- pi / 180
  dlon <- (lon2 - lon1) * rad; dlat <- (lat2 - lat1) * rad
  a <- sin(dlat / 2)^2 + cos(lat1 * rad) * cos(lat2 * rad) * sin(dlon / 2)^2
  6371.0088 * 2 * atan2(sqrt(a), sqrt(1 - a))
}

n <- nrow(z)
geo <- matrix(0, n, n, dimnames = list(z$Site, z$Site))
for (i in seq_len(n - 1L)) for (j in (i + 1L):n) {
  geo[i, j] <- geo[j, i] <- haversine_km(z$Longitude[i], z$Latitude[i],
                                         z$Longitude[j], z$Latitude[j])
}
unlink(file.path(result_dir, "locality_geographic_distance_km.csv"))

# All 4! locality-label permutations; exact tests include the observed order.
perm_grid <- as.matrix(expand.grid(rep(list(seq_len(n)), n)))
perms <- perm_grid[apply(perm_grid, 1, function(v) length(unique(v)) == n), , drop = FALSE]
stopifnot(nrow(perms) == factorial(n))
lower <- lower.tri(geo)

exact_spearman <- function(x, y) {
  obs <- suppressWarnings(cor(x, y, method = "spearman"))
  null <- apply(perms, 1, function(idx)
    suppressWarnings(cor(x, y[idx], method = "spearman")))
  c(Effect = obs, Exact_p = mean(abs(null) >= abs(obs) - 1e-12))
}

outcomes <- list(
  `Anti-Bd beta-binomial mean` = z$logit_anti_bd_mean,
  `Anti-Bd genus Hill D1` = z$log_hill_d1,
  `Whole-community genus Hill D1` = z$log_whole_hill_d1
)
tests <- do.call(rbind, lapply(names(outcomes), function(nm) {
  y <- outcomes[[nm]]
  e <- exact_spearman(z$Elevation_m, y)
  la <- exact_spearman(z$Latitude, y)
  lo <- exact_spearman(z$Longitude, y)
  data.frame(
    Endpoint = nm,
    Predictor = c("Longitude", "Latitude", "Elevation"),
    Test = "Exact Spearman",
    Effect = c(lo["Effect"], la["Effect"], e["Effect"]),
    Exact_p = c(lo["Exact_p"], la["Exact_p"], e["Exact_p"]),
    Independent_localities = n,
    Permutations = nrow(perms), stringsAsFactors = FALSE
  )
}))
tests$Holm_p <- p.adjust(tests$Exact_p, method = "holm")
write.csv(tests, file.path(result_dir, "exact_spatial_elevation_tests.csv"), row.names = FALSE)
unlink(file.path(result_dir, "descriptive_pairwise_spatial_differences.csv"))

theme_pub <- theme_classic(base_size = 8) +
  theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7), legend.position = "none")

make_panel <- function(xvar, yvar, colour, title, xlab, ylab) {
  ggplot(z, aes(x = .data[[xvar]], y = .data[[yvar]])) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.55, colour = "grey40") +
    geom_point(size = 2.3, colour = colour) +
    ggrepel::geom_text_repel(aes(label = Site), size = 2.0, max.overlaps = Inf) +
    labs(title = title, subtitle = "Four independent localities; line is descriptive",
         x = xlab, y = ylab) + theme_pub
}

p_a <- make_panel("Longitude", "logit_anti_bd_mean", "#8E3B73",
                  "a  Relative abundance vs longitude", "Longitude (°E)",
                  "Beta-binomial mean (logit)")
p_b <- make_panel("Latitude", "logit_anti_bd_mean", "#8E3B73",
                  "b  Relative abundance vs latitude", "Latitude (°N)",
                  "Beta-binomial mean (logit)")
p_c <- make_panel("Elevation_m", "logit_anti_bd_mean", "#8E3B73",
                  "c  Relative abundance vs elevation", "Elevation (m)",
                  "Beta-binomial mean (logit)")

div_long <- rbind(
  transform(z[, c("Site", "Longitude", "Latitude", "Elevation_m", "log_hill_d1")],
            Community = "Anti-Bd", log_D1 = log_hill_d1)[,
              c("Site", "Longitude", "Latitude", "Elevation_m", "Community", "log_D1")],
  transform(z[, c("Site", "Longitude", "Latitude", "Elevation_m", "log_whole_hill_d1")],
            Community = "Whole community", log_D1 = log_whole_hill_d1)[,
              c("Site", "Longitude", "Latitude", "Elevation_m", "Community", "log_D1")]
)
div_long$Community <- factor(div_long$Community,
                             levels = c("Whole community", "Anti-Bd"))

make_div_panel <- function(xvar, title, xlab) {
  ggplot(div_long, aes(x = .data[[xvar]], y = log_D1,
                       colour = Community, shape = Community)) +
    geom_smooth(aes(group = Community), method = "lm", se = FALSE,
                linewidth = 0.6, show.legend = FALSE) +
    geom_point(size = 2.4) +
    scale_colour_manual(values = c("Whole community" = "#0072B2", "Anti-Bd" = "#D55E00")) +
    scale_shape_manual(values = c("Whole community" = 16, "Anti-Bd" = 17)) +
    labs(title = title, subtitle = "Whole and anti-Bd shown on the same log scale",
         x = xlab, y = expression("log(Hill "^1*D*")"), colour = NULL, shape = NULL) +
    theme_pub + theme(legend.position = "top", legend.text = element_text(size = 6),
                      legend.key.width = grid::unit(8, "mm"))
}

p_d <- make_div_panel("Longitude", expression("d  Hill "^1*D*" vs longitude"),
                      "Longitude (°E)")
p_e <- make_div_panel("Latitude", expression("e  Hill "^1*D*" vs latitude"),
                      "Latitude (°N)")
p_f <- make_div_panel("Elevation_m", expression("f  Hill "^1*D*" vs elevation"),
                      "Elevation (m)")

fig <- ((p_a | p_b | p_c) / (p_d | p_e | p_f)) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Coordinate and elevational associations of anti-Bd abundance and bacterial diversity",
    subtitle = "Inference uses all 24 exact permutations of four locality labels",
    caption = paste0("Longitude, latitude and elevation are tested separately using exact locality-level Spearman correlations.\n",
                     "Holm correction is applied across all nine tests; fitted lines are descriptive."),
    theme = theme(plot.caption = element_text(size = 6.2, hjust = 0),
                  plot.margin = margin(5, 8, 7, 8))
  ) & theme(legend.position = "bottom")

save_plot <- function(plot, stem, width_mm = 183, height_mm = 150) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  ggsave(paste0(stem, ".png"), plot, width = w, height = h, dpi = 600, bg = "white")
  svglite::svglite(paste0(stem, ".svg"), width = w, height = h); print(plot); dev.off()
  pdf(paste0(stem, ".pdf"), width = w, height = h, useDingbats = FALSE); print(plot); dev.off()
  ragg::agg_tiff(paste0(stem, ".tiff"), width = w, height = h, units = "in",
                 res = 600, compression = "lzw", background = "white")
  print(plot); dev.off()
}
save_plot(fig, file.path(figure_dir, "figure_spatial_elevation_abundance_diversity"),
          width_mm = 183, height_mm = 132)

write.csv(data.frame(
  Item = c("Independent unit", "Coordinate method", "Elevation method",
           "Relative-abundance scale", "Diversity scale", "Multiplicity"),
  Definition = c("Locality (n = 4)",
                 "Separate exact locality-level Spearman correlations for longitude and latitude",
                 "Exact locality-level Spearman correlation",
                 "Logit of beta-binomial conditional mean anti-Bd proportion",
                 "Log of coverage-standardized genus-level Hill q = 1",
                 "Holm correction across 3 endpoints x 3 predictors")
), file.path(result_dir, "analysis_methods.csv"), row.names = FALSE)

print(z[, c("Site", "Latitude", "Longitude", "Elevation_m",
            "beta_binomial_mean_percent", "d1_mean", "whole_d1_mean")])
print(tests)
