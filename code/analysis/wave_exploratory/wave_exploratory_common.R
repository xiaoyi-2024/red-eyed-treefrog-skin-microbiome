suppressPackageStartupMessages({
  library(phyloseq)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
  library(ggrepel)
  library(maps)
  library(terra)
})

result_dir <- "results/wave_exploratory"
figure_dir <- file.path(result_dir, "figures")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

full_file <- "data/processed/ps_clean_47_rebuilt_from_01_filtered.rds"
anti_file <- "data/processed/blast_97id_90cov/ps_anti_bd_47_rebuilt_97id_90cov.rds"
reference_file <- "data/reference/bd_wave_reference_table.csv"
stopifnot(file.exists(full_file), file.exists(anti_file), file.exists(reference_file))

otu_matrix <- function(ps) {
  x <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) x <- t(x)
  x
}

# Recompute the microbiome metrics from the rebuilt 47-sample objects.
ps_full <- readRDS(full_file)
ps_anti <- readRDS(anti_file)
stopifnot(nsamples(ps_full) == 47L, nsamples(ps_anti) == 47L)
md <- data.frame(sample_data(ps_full), check.names = FALSE, stringsAsFactors = FALSE)
md$Sample.ID <- rownames(md)
site_key <- c(
  "Lost Iguana Resort, La Fortuna, San Carlos" = "Lost Iguana",
  "Veragua Rainforest, Las Brisas de Veragua" = "Veragua",
  "PN Altos de Campana" = "Altos de Campana",
  "PN Soberanía, Las Cumbres" = "Soberanía"
)
md$Site <- unname(site_key[md$Locality])
stopifnot(!anyNA(md$Site), length(unique(md$Site)) == 4L)

sites <- aggregate(cbind(latitude = md$Lat, longitude = md$Long),
                   by = list(Site = md$Site), FUN = mean)
sites$Sample_n <- as.integer(table(factor(md$Site, levels = sites$Site)))

full_counts <- otu_matrix(ps_full)
anti_counts <- otu_matrix(ps_anti)[rownames(full_counts), , drop = FALSE]
sample_metric <- data.frame(
  Sample.ID = rownames(full_counts), Site = md[rownames(full_counts), "Site"],
  anti_bd_percent = 100 * rowSums(anti_counts) / rowSums(full_counts),
  stringsAsFactors = FALSE
)

# Aitchison distance after adding a pseudocount of 1 to every ASV count.
clr <- log(full_counts + 1)
clr <- clr - rowMeans(clr)
sample_metric$distance_to_site_centroid <- NA_real_
for (s in unique(sample_metric$Site)) {
  idx <- which(sample_metric$Site == s)
  centroid <- colMeans(clr[idx, , drop = FALSE])
  sample_metric$distance_to_site_centroid[idx] <-
    sqrt(rowSums((clr[idx, , drop = FALSE] - centroid)^2))
}
site_summary <- do.call(rbind, lapply(split(sample_metric, sample_metric$Site), function(z) {
  data.frame(
    Site = z$Site[1], Sample_n = nrow(z),
    anti_mean_percent = mean(z$anti_bd_percent),
    anti_se_percent = sd(z$anti_bd_percent) / sqrt(nrow(z)),
    dispersion_mean = mean(z$distance_to_site_centroid),
    dispersion_se = sd(z$distance_to_site_centroid) / sqrt(nrow(z))
  )
}))
rownames(site_summary) <- NULL

# Build one directed Central American Bd wave from the published decline records.
reference <- read.csv(reference_file, stringsAsFactors = FALSE, check.names = FALSE,
                      na.strings = c("", "NA"))
stopifnot(nrow(reference) == 59L)
reference$DOD_year <- rowMeans(reference[, c("DOD_low", "DOD_high")], na.rm = TRUE)
reference$DOD_year[!is.finite(reference$DOD_year)] <- NA_real_
path_names <- c("A. senex", "A. chiriquiensis", "A. varius Fortuna",
                "A. varius Santa Fe", "A. varius El Cope", "A. zeteki El Valle")
map_point_names <- c(path_names, "A. varius Monteverde")
map_points <- reference[match(map_point_names, reference$species_event), , drop = FALSE]
stopifnot(nrow(map_points) == 7L, !anyNA(map_points$species_event),
          all(map_points$status == "Decline"), all(is.finite(map_points$DOD_year)))
# Preserve the Wave column only in the complete transcription; it is removed
# from all model-specific data because this reconstruction has one wave.
map_points$wave <- NULL

path <- map_points[match(path_names, map_points$species_event), , drop = FALSE]
stopifnot(!anyNA(path$species_event), all(diff(path$DOD_year) > 0))
path$path_order <- seq_len(nrow(path))
path$year_label_lon <- path$longitude + c(-0.18, 0.18, -0.18, 0.18, -0.18, -0.18)
path$year_label_lat <- path$latitude + c(-0.08, 0.08, -0.08, 0.08, 0.08, -0.15)
path_match <- match(map_points$species_event, path$species_event)
map_points$year_label_lon <- ifelse(
  is.na(path_match), map_points$longitude + 0.08, path$year_label_lon[path_match]
)
map_points$year_label_lat <- ifelse(
  is.na(path_match), map_points$latitude + 0.06, path$year_label_lat[path_match]
)

earth_radius_km <- 6371.0088
haversine_km <- function(lon1, lat1, lon2, lat2) {
  rad <- pi / 180
  dlon <- (lon2 - lon1) * rad; dlat <- (lat2 - lat1) * rad
  a <- sin(dlat / 2)^2 + cos(lat1 * rad) * cos(lat2 * rad) * sin(dlon / 2)^2
  2 * earth_radius_km * atan2(sqrt(a), sqrt(1 - a))
}

segments <- data.frame(
  segment = seq_len(nrow(path) - 1L),
  from_event = head(path$species_event, -1), to_event = tail(path$species_event, -1),
  lon1 = head(path$longitude, -1), lat1 = head(path$latitude, -1),
  lon2 = tail(path$longitude, -1), lat2 = tail(path$latitude, -1),
  year1 = head(path$DOD_year, -1), year2 = tail(path$DOD_year, -1),
  stringsAsFactors = FALSE
)
segments$distance_km <- with(segments, haversine_km(lon1, lat1, lon2, lat2))
segments$duration_years <- segments$year2 - segments$year1
segments$speed_km_per_year <- segments$distance_km / segments$duration_years
segments$cumulative_start_km <- c(0, head(cumsum(segments$distance_km), -1))
segments$cumulative_end_km <- cumsum(segments$distance_km)
segments$mid_lon <- (segments$lon1 + segments$lon2) / 2
segments$mid_lat <- (segments$lat1 + segments$lat2) / 2 + c(0.05, -0.08, 0.08, -0.08, 0.11)
segments$speed_label <- sprintf("%.0f km/yr", segments$speed_km_per_year)

# Project each sampling site to the nearest segment in a local equirectangular
# coordinate system, then linearly interpolate the wave year along that segment.
lat0 <- mean(path$latitude) * pi / 180
to_xy <- function(lon, lat) {
  cbind(x = earth_radius_km * cos(lat0) * lon * pi / 180,
        y = earth_radius_km * lat * pi / 180)
}
project_site <- function(lon, lat) {
  p <- as.numeric(to_xy(lon, lat))
  candidates <- lapply(seq_len(nrow(segments)), function(i) {
    a <- as.numeric(to_xy(segments$lon1[i], segments$lat1[i]))
    b <- as.numeric(to_xy(segments$lon2[i], segments$lat2[i]))
    ab <- b - a
    raw_t <- sum((p - a) * ab) / sum(ab^2)
    # The conceptual source shows the wave continuing beyond both terminal
    # anchors. Permit directional extrapolation only on the first/last segment.
    lower <- if (i == 1L) -Inf else 0
    upper <- if (i == nrow(segments)) Inf else 1
    t <- max(lower, min(upper, raw_t))
    q <- a + t * ab
    data.frame(segment = i, fraction = t, distance_to_path_km = sqrt(sum((p - q)^2)),
               projected_x = q[1], projected_y = q[2])
  })
  out <- do.call(rbind, candidates)
  out[which.min(out$distance_to_path_km), , drop = FALSE]
}
projection <- do.call(rbind, Map(project_site, sites$longitude, sites$latitude))
rownames(projection) <- NULL
sites <- cbind(sites, projection)
sites$pred_wave_year <- segments$year1[sites$segment] + sites$fraction *
  (segments$year2[sites$segment] - segments$year1[sites$segment])
sites$along_path_km <- segments$cumulative_start_km[sites$segment] +
  sites$fraction * segments$distance_km[sites$segment]
sites$years_before_sampling <- 2022 - sites$pred_wave_year
sites$Prediction_status <- ifelse(
  sites$fraction < 0, "Extrapolated before the first anchor along segment 1",
  ifelse(sites$fraction > 1, "Extrapolated beyond the last anchor along segment 5",
         "Interpolated on the nearest segment of the single wave path")
)
sites$label_lon <- sites$longitude + c(-0.08, 0.08, -0.08, 0.08)
sites$label_lat <- sites$latitude + c(-0.12, 0.07, -0.10, 0.07)
sites$label_hjust <- c(1, 0, 1, 0)

metrics <- merge(sites, site_summary, by = c("Site", "Sample_n"), sort = FALSE)
metrics <- metrics[match(sites$Site, metrics$Site), ]
test_anti <- cor.test(metrics$pred_wave_year, metrics$anti_mean_percent,
                      method = "spearman", exact = TRUE)
test_disp <- cor.test(metrics$pred_wave_year, metrics$dispersion_mean,
                      method = "spearman", exact = TRUE)
correlations <- data.frame(
  Model = "Single directed Central American wave",
  Outcome = c("Mean anti-Bd relative abundance (%)",
              "Mean Aitchison distance to site centroid"),
  Site_n = 4L,
  Spearman_rho = c(unname(test_anti$estimate), unname(test_disp$estimate)),
  Exact_p = c(test_anti$p.value, test_disp$p.value)
)

write.csv(segments, file.path(result_dir, "single_wave_path_segments_and_speeds.csv"), row.names = FALSE)
write.csv(metrics, file.path(result_dir, "site_metrics_single_wave.csv"), row.names = FALSE)
write.csv(correlations, file.path(result_dir, "site_correlations_single_wave.csv"), row.names = FALSE)

pal <- c("Lost Iguana" = "#CC79A7", "Veragua" = "#0072B2",
         "Altos de Campana" = "#009E73", "Soberanía" = "#E69F00")
theme_set(theme_classic(base_size = 7, base_family = "sans") +
  theme(axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35),
        plot.title = element_text(size = 8, face = "bold"),
        plot.subtitle = element_text(size = 6.5), legend.text = element_text(size = 6),
        legend.title = element_text(size = 6.5)))

coast <- map_data("world", region = c("Costa Rica", "Panama"))

# WorldClim 30 arc-second elevation, retained locally for reproducibility.
dem_cr <- terra::rast("data/reference/elevation/CRI_elv_msk.tif")
dem_pa <- terra::rast("data/reference/elevation/PAN_elv_msk.tif")
dem <- terra::merge(dem_cr, dem_pa)
dem <- terra::crop(dem, terra::ext(-85.35, -79.30, 8.25, 10.65))
names(dem) <- "Elevation_m"
dem_plot <- terra::aggregate(dem, fact = 3, fun = mean, na.rm = TRUE)
slope <- terra::terrain(dem_plot, v = "slope", unit = "radians")
aspect <- terra::terrain(dem_plot, v = "aspect", unit = "radians")
hill <- terra::shade(slope, aspect, angle = 40, direction = 315)
dem_df <- as.data.frame(dem_plot, xy = TRUE, na.rm = TRUE)
hill_df <- as.data.frame(hill, xy = TRUE, na.rm = TRUE)
names(hill_df)[3] <- "Hillshade"
dem_df <- merge(dem_df, hill_df, by = c("x", "y"))

p1 <- ggplot() +
  geom_raster(data = dem_df, aes(x, y, fill = Elevation_m)) +
  geom_raster(data = dem_df, aes(x, y, alpha = Hillshade), fill = "black",
              inherit.aes = FALSE) +
  geom_polygon(data = coast, aes(long, lat, group = group), fill = NA,
               colour = "#43484B", linewidth = 0.35) +
  # Broad corridor plus a high-contrast centreline with one terminal arrow.
  geom_path(data = path, aes(longitude, latitude), colour = "#F0B37E",
            linewidth = 4.2, alpha = 0.48, lineend = "round", linejoin = "round") +
  geom_path(data = path, aes(longitude, latitude), colour = "#7A1F5C",
            linewidth = 0.85, lineend = "round", linejoin = "round") +
  # A separate summary arrow communicates the overall west-to-east direction.
  annotate("segment", x = -82.85, y = 10.37, xend = -80.35, yend = 10.37,
           colour = "#7A1F5C", linewidth = 1.05,
           arrow = arrow(length = grid::unit(3.2, "mm"), type = "closed")) +
  annotate("text", x = -81.60, y = 10.47, label = "Bd wave direction",
           colour = "#7A1F5C", size = 2.25, fontface = "bold") +
  geom_point(data = map_points, aes(longitude, latitude), shape = 21,
             size = 2.3, stroke = 0.65, fill = "white", colour = "black") +
  geom_text(data = map_points,
            aes(year_label_lon, year_label_lat, label = round(DOD_year)), size = 2.2) +
  geom_label(data = segments, aes(mid_lon, mid_lat, label = speed_label), size = 1.9,
             linewidth = 0.15, fill = scales::alpha("white", 0.75)) +
  geom_point(data = sites, aes(longitude, latitude, colour = Site), size = 2.8) +
  geom_text(data = sites[sites$Site != "Lost Iguana", ],
            aes(label_lon, label_lat, label = Site, colour = Site,
                              hjust = label_hjust), size = 2.1, show.legend = FALSE) +
  scale_fill_gradientn(
    colours = c("#DDEBC7", "#B9D39A", "#91B873", "#C6A56B", "#92704D", "#F2EEE7"),
    values = scales::rescale(c(0, 100, 500, 1200, 2200, 3700)),
    limits = c(0, 3700), oob = scales::squish, name = "Elevation (m)"
  ) +
  scale_alpha(range = c(0.01, 0.22), guide = "none") +
  scale_colour_manual(values = pal, guide = "none") +
  coord_quickmap(xlim = c(-85.35, -79.30), ylim = c(8.25, 10.65), expand = FALSE) +
  labs(title = "Bd wave through Costa Rica and Panama",
       subtitle = NULL, x = "Longitude", y = "Latitude", caption = NULL) +
  theme(
    panel.background = element_rect(fill = "#CFE8F3", colour = NA),
    panel.grid.major = element_line(colour = scales::alpha("white", 0.65), linewidth = 0.25),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.key.height = grid::unit(14, "mm"),
    legend.title = element_text(face = "bold")
  )

p2 <- ggplot(metrics, aes(pred_wave_year, anti_mean_percent, colour = Site)) +
  geom_errorbar(aes(ymin = pmax(0, anti_mean_percent - anti_se_percent),
                    ymax = anti_mean_percent + anti_se_percent), width = 0.25, linewidth = 0.4) +
  geom_point(size = 2.8) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, colour = "grey35",
              linewidth = 0.5, linetype = "dashed") +
  geom_text_repel(aes(label = Site), size = 2, seed = 20260814, show.legend = FALSE) +
  scale_colour_manual(values = pal) +
  labs(title = "a  Wave year and anti-Bd abundance",
       subtitle = sprintf("Four sites; Spearman rho = %.2f, exact p = %.3f",
                          correlations$Spearman_rho[1], correlations$Exact_p[1]),
       x = "Estimated Bd wave year", y = "Mean anti-Bd relative abundance (%)") +
  theme(legend.position = "none")

p3 <- ggplot(metrics, aes(pred_wave_year, dispersion_mean, colour = Site)) +
  geom_errorbar(aes(ymin = pmax(0, dispersion_mean - dispersion_se),
                    ymax = dispersion_mean + dispersion_se), width = 0.25, linewidth = 0.4) +
  geom_point(size = 2.8) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, colour = "grey35",
              linewidth = 0.5, linetype = "dashed") +
  geom_text_repel(aes(label = Site), size = 2, seed = 20260814, show.legend = FALSE) +
  scale_colour_manual(values = pal) +
  labs(title = "b  Wave year and community dispersion",
       subtitle = sprintf("Four sites; Spearman rho = %.2f, exact p = %.3f",
                          correlations$Spearman_rho[2], correlations$Exact_p[2]),
       x = "Estimated Bd wave year", y = "Mean Aitchison distance to site centroid") +
  theme(legend.position = "none")

p23 <- (p2 | p3) + plot_annotation(
  title = "Exploratory microbiome associations with the Bd wave",
  subtitle = "Site-level summaries are descriptive (n = 4 sites) and do not establish temporal effects or causality."
)

# Standalone visualization of the piecewise arrival-time model.
path$cumulative_distance_km <- c(0, cumsum(segments$distance_km))
model_min_x <- min(c(0, metrics$along_path_km))
model_max_x <- max(c(max(path$cumulative_distance_km), metrics$along_path_km))
first_speed <- segments$speed_km_per_year[1]
last_speed <- tail(segments$speed_km_per_year, 1)
first_year <- path$DOD_year[1]
last_year <- tail(path$DOD_year, 1)
path_end_x <- max(path$cumulative_distance_km)
arrival_label_offsets <- data.frame(
  Site = c("Altos de Campana", "Lost Iguana", "Soberanía", "Veragua"),
  dx = c(5, 12, -4, -16), dy = c(-1.35, 0, 1.25, 1.25),
  hjust = c(0.5, 0, 0.5, 1), stringsAsFactors = FALSE
)
metrics <- merge(metrics, arrival_label_offsets, by = "Site", all.x = TRUE, sort = FALSE)
metrics <- metrics[match(sites$Site, metrics$Site), ]
metrics$arrival_label_x <- metrics$along_path_km + metrics$dx
metrics$arrival_label_y <- metrics$pred_wave_year + metrics$dy

p_arrival <- ggplot() +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf,
           fill = "#F2F2F2", alpha = 0.85) +
  annotate("rect", xmin = path_end_x, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#F2F2F2", alpha = 0.85) +
  annotate("segment", x = model_min_x, xend = 0,
           y = first_year + model_min_x / first_speed, yend = first_year,
           colour = "#7A1F5C", linewidth = 0.7, linetype = "dashed") +
  geom_segment(data = segments,
               aes(x = cumulative_start_km, xend = cumulative_end_km,
                   y = year1, yend = year2),
               colour = "#7A1F5C", linewidth = 1.0, lineend = "round") +
  annotate("segment", x = path_end_x, xend = model_max_x,
           y = last_year, yend = last_year + (model_max_x - path_end_x) / last_speed,
           colour = "#7A1F5C", linewidth = 0.7, linetype = "dashed") +
  geom_point(data = path, aes(cumulative_distance_km, DOD_year),
             shape = 21, size = 2.6, stroke = 0.7, fill = "white", colour = "black") +
  geom_text(data = path, aes(cumulative_distance_km, DOD_year, label = round(DOD_year)),
            nudge_y = 0.8, size = 2.1) +
  geom_point(data = metrics, aes(along_path_km, pred_wave_year, colour = Site),
             size = 3.0) +
  geom_text(data = metrics,
            aes(arrival_label_x, arrival_label_y,
                label = sprintf("%s\n%.1f", Site, pred_wave_year),
                colour = Site, hjust = hjust), size = 2.1, show.legend = FALSE) +
  annotate("text", x = model_min_x / 2, y = 2010.5, label = "Directional\nextrapolation",
           size = 2.0, colour = "grey35") +
  annotate("text", x = (path_end_x + model_max_x) / 2, y = 1981.5,
           label = "Directional\nextrapolation", size = 2.0, colour = "grey35") +
  scale_colour_manual(values = pal, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.07, 0.07))) +
  scale_y_continuous(breaks = seq(1985, 2010, 5),
                     expand = expansion(mult = c(0.08, 0.10))) +
  labs(title = "Bd arrival-time model",
       subtitle = "Piecewise interpolation between published DOD anchors; dashed lines denote endpoint extrapolation",
       x = "Distance along reconstructed wave path (km)",
       y = "Bd arrival year") +
  theme(panel.grid.major = element_line(colour = "#E6E6E6", linewidth = 0.3),
        panel.grid.minor = element_blank())

save_plot <- function(plot, stem, width_mm, height_mm) {
  width_in <- width_mm / 25.4; height_in <- height_mm / 25.4
  ggsave(paste0(stem, ".png"), plot, width = width_in, height = height_in, dpi = 600, bg = "white")
  ragg::agg_tiff(paste0(stem, ".tiff"), width = width_in, height = height_in,
                 units = "in", res = 600, background = "white"); print(plot); dev.off()
  svglite::svglite(paste0(stem, ".svg"), width = width_in, height = height_in); print(plot); dev.off()
  grDevices::pdf(paste0(stem, ".pdf"), width = width_in, height = height_in,
                 family = "Helvetica", useDingbats = FALSE); print(plot); dev.off()
}
save_plot(p1, file.path(figure_dir, "figure1_single_bd_wave"), 150, 90)
save_plot(p23, file.path(figure_dir, "figure2_3_single_wave_microbiome"), 183, 82)
save_plot(p_arrival, file.path(figure_dir, "figure_arrival_time_model"), 140, 92)

print(segments[, c("from_event", "to_event", "distance_km", "duration_years", "speed_km_per_year")])
print(metrics[, c("Site", "pred_wave_year", "distance_to_path_km")])
print(correlations)
