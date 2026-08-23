suppressPackageStartupMessages({
  library(phyloseq)
  library(ggplot2)
  library(svglite)
  library(ragg)
  library(ggrepel)
})

# Reproduce panel pA from the original v3_dow workflow, using the current
# submission2 reference table and rebuilt 47-sample phyloseq object.
result_dir <- "results/wave_exploratory"
figure_dir <- file.path(result_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

reference <- read.csv("data/reference/bd_wave_reference_table.csv",
                      stringsAsFactors = FALSE, check.names = FALSE,
                      na.strings = c("", "NA"))
anchor_names <- c("A. senex", "A. varius Monteverde", "A. chiriquiensis",
                  "A. varius Fortuna", "A. varius Santa Fe",
                  "A. varius El Cope", "A. zeteki El Valle")
anchors <- reference[match(anchor_names, reference$species_event), ]
anchors$dow <- rowMeans(anchors[, c("DOD_low", "DOD_high")], na.rm = TRUE)
anchors$lon <- anchors$longitude
anchors$lat <- anchors$latitude
stopifnot(nrow(anchors) == 7L, !anyNA(anchors[, c("dow", "lon", "lat")]))

ps <- readRDS("data/processed/ps_clean_47_rebuilt_from_01_filtered.rds")
md <- data.frame(sample_data(ps), check.names = FALSE, stringsAsFactors = FALSE)
site_key <- c(
  "Lost Iguana Resort, La Fortuna, San Carlos" = "Lost Iguana",
  "Veragua Rainforest, Las Brisas de Veragua" = "Veragua",
  "PN Altos de Campana" = "Altos de Campana",
  "PN Soberanía, Las Cumbres" = "Soberanía"
)
md$Site_short <- unname(site_key[md$Locality])
sites <- aggregate(cbind(Latitude = md$Lat, Longitude = md$Long),
                   by = list(Site_short = md$Site_short), FUN = mean)
stopifnot(nrow(sites) == 4L, !anyNA(sites))

fit <- lm(dow ~ lon + lat, data = anchors)
lon_seq <- seq(min(c(anchors$lon, sites$Longitude)) - 0.2,
               max(c(anchors$lon, sites$Longitude)) + 0.2, length.out = 160)
lat_seq <- seq(min(c(anchors$lat, sites$Latitude)) - 0.2,
               max(c(anchors$lat, sites$Latitude)) + 0.2, length.out = 160)
grid <- expand.grid(lon = lon_seq, lat = lat_seq)
grid$pred <- as.numeric(predict(fit, newdata = grid))

# For a planar model, each fitted arrival-year contour is a straight line:
# year = beta0 + beta_lon * longitude + beta_lat * latitude.
fit_years <- seq(1988, 2008, 4)
coef_fit <- coef(fit)
fitted_lines <- do.call(rbind, lapply(fit_years, function(year) {
  latitude <- seq(min(lat_seq), max(lat_seq), length.out = 120)
  longitude <- (year - coef_fit[1] - coef_fit["lat"] * latitude) / coef_fit["lon"]
  data.frame(year = year, lon = longitude, lat = latitude)
}))
fitted_lines <- fitted_lines[
  fitted_lines$lon >= min(lon_seq) & fitted_lines$lon <= max(lon_seq), ]
line_labels <- data.frame(year = fit_years, lat = min(lat_seq) + 0.10)
line_labels$lon <- (line_labels$year - coef_fit[1] -
                      coef_fit["lat"] * line_labels$lat) / coef_fit["lon"]

pal <- c("Lost Iguana" = "#CC79A7", "Veragua" = "#0072B2",
         "Altos de Campana" = "#009E73", "Soberanía" = "#E69F00")

pA <- ggplot() +
  geom_raster(data = grid, aes(lon, lat, fill = pred), alpha = 0.82) +
  geom_line(data = fitted_lines, aes(lon, lat, group = year),
            colour = "white", linewidth = 0.42, alpha = 0.95) +
  geom_label(data = line_labels, aes(lon, lat, label = year),
             size = 1.9, linewidth = 0.12, label.padding = grid::unit(0.10, "lines"),
             fill = scales::alpha("white", 0.82), colour = "#6A3030") +
  geom_point(data = anchors, aes(lon, lat), shape = 21, size = 2.2,
             stroke = 0.5, fill = "white", colour = "black") +
  geom_text(data = anchors, aes(lon, lat, label = round(dow)),
            nudge_y = 0.10, size = 2.2) +
  geom_point(data = sites, aes(Longitude, Latitude, colour = Site_short), size = 3) +
  geom_text_repel(data = sites,
                  aes(Longitude, Latitude, label = Site_short, colour = Site_short),
                  seed = 20260814, size = 2.2, min.segment.length = 0,
                  show.legend = FALSE) +
  scale_fill_gradient(low = "#E6F2F2", high = "#A84949",
                      name = "Predicted\narrival year") +
  scale_colour_manual(values = pal, guide = "none") +
  coord_equal() +
  labs(title = "Approximate Bd wave-arrival surface",
       subtitle = "Linear spatial model fitted to seven published DOD anchors",
       x = "Longitude", y = "Latitude") +
  theme_classic(base_size = 8, base_family = "sans") +
  theme(axis.line = element_line(linewidth = 0.35),
        axis.ticks = element_line(linewidth = 0.35),
        plot.title = element_text(size = 9, face = "bold"),
        plot.subtitle = element_text(size = 7),
        legend.position = "right", legend.title = element_text(size = 7),
        legend.text = element_text(size = 6))

stem <- file.path(figure_dir, "figure_v3_dow_first_panel")
width_in <- 150 / 25.4
height_in <- 105 / 25.4
ggsave(paste0(stem, ".png"), pA, width = width_in, height = height_in,
       dpi = 600, bg = "white")
ragg::agg_tiff(paste0(stem, ".tiff"), width = width_in, height = height_in,
               units = "in", res = 600, background = "white")
print(pA)
dev.off()
svglite::svglite(paste0(stem, ".svg"), width = width_in, height = height_in)
print(pA)
dev.off()
grDevices::pdf(paste0(stem, ".pdf"), width = width_in, height = height_in,
               family = "Helvetica", useDingbats = FALSE)
print(pA)
dev.off()

model_summary <- data.frame(
  Formula = "DOD_year ~ longitude + latitude",
  Anchor_n = nrow(anchors),
  R_squared = summary(fit)$r.squared,
  Adjusted_R_squared = summary(fit)$adj.r.squared,
  RMSE_years = sqrt(mean(residuals(fit)^2)),
  LOOCV_RMSE_years = sqrt(mean(vapply(seq_len(nrow(anchors)), function(i) {
    fit_i <- lm(dow ~ lon + lat, data = anchors[-i, , drop = FALSE])
    anchors$dow[i] - predict(fit_i, newdata = anchors[i, , drop = FALSE])
  }, numeric(1))^2))
)
coef_table <- data.frame(
  Term = names(coef(fit)), Estimate = unname(coef(fit)),
  SE = coef(summary(fit))[, "Std. Error"],
  t_value = coef(summary(fit))[, "t value"],
  p_value = coef(summary(fit))[, "Pr(>|t|)"], row.names = NULL
)
site_prediction <- cbind(
  sites,
  as.data.frame(predict(
    fit,
    newdata = data.frame(lon = sites$Longitude, lat = sites$Latitude),
    interval = "prediction", level = 0.95
  ))
)
names(site_prediction)[names(site_prediction) == "fit"] <- "Predicted_DOD_year"
names(site_prediction)[names(site_prediction) == "lwr"] <- "Prediction_low_95"
names(site_prediction)[names(site_prediction) == "upr"] <- "Prediction_high_95"
anchor_output <- anchors[, c("species_event", "country", "dow", "lon", "lat")]
names(anchor_output) <- c("Species_event", "Country", "DOD_midpoint_year",
                          "Longitude", "Latitude")
write.csv(anchor_output,
          file.path(result_dir, "bd_wave_spatial_model_anchors.csv"), row.names = FALSE)
write.csv(model_summary,
          file.path(result_dir, "bd_wave_spatial_model_summary.csv"), row.names = FALSE)
write.csv(coef_table,
          file.path(result_dir, "bd_wave_spatial_model_coefficients.csv"), row.names = FALSE)
write.csv(site_prediction,
          file.path(result_dir, "bd_wave_spatial_model_site_predictions.csv"), row.names = FALSE)
print(model_summary)
print(site_prediction)
