library(sf)
library(rnaturalearth)
library(ggplot2)
library(ggrepel)
library(ggspatial)
library(terra)
library(geodata)
source("code/analysis/common.R")
figures_root <- "results/maps/figures"
dir.create(figures_root, recursive = TRUE, showWarnings = FALSE)
input_file <- "data/sampling_localities.csv"
output_dir <- "results/maps/tables"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

sites <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)
sites$ShortName <- c(
  "Lost Iguana Resort",
  "Veragua Rainforest",
  "PN Altos de Campana",
  "PN Soberanía"
)
sites$BiogeographicRegion <- c(
  "Northeastern Costa Rica",
  "Southeastern Costa Rica",
  "Central Panama",
  "Central Panama"
)
# Derive sample counts from the final 47-sample object, not a manually entered
# table, so labels remain synchronized with the analysis population.
final_counts <- table(site_factor(ps_clean))
sites$SampleCount <- as.integer(final_counts[sites$ShortName])
stopifnot(!anyNA(sites$SampleCount),sum(sites$SampleCount)==47L)
elevation_dir <- "data/reference/elevation"
dir.create(elevation_dir, recursive = TRUE, showWarnings = FALSE)
dem_cr <- elevation_30s(country = "CRI", path = elevation_dir)
dem_pa <- elevation_30s(country = "PAN", path = elevation_dir)
dem <- merge(dem_cr, dem_pa)
dem <- crop(dem, ext(-86.2, -78.9, 7.7, 11.5))
names(dem) <- "Elevation_m"
site_elevation <- terra::extract(dem, sites[, c("Longitude", "Latitude")], ID = FALSE)[, 1]
sites$Elevation_m_DEM <- round(site_elevation)
sites$Elevation_source <- "GADM/WorldClim 30 arc-second DEM via geodata; coordinate extraction"
sites$MapLabel <- sprintf(
  "%d  %s\n%s\n%d m  |  n=%d",
  sites$SiteID, sites$ShortName, sites$BiogeographicRegion,
  sites$Elevation_m_DEM, sites$SampleCount
)
write.csv(sites,file.path(output_dir,"sampling_localities_47.csv"),row.names=FALSE)

# Downsample only the display raster; point elevation extraction above uses the
# full 30 arc-second DEM. Hillshade preserves mountain relief without a rainbow palette.
dem_plot <- aggregate(dem, fact = 4, fun = mean, na.rm = TRUE)
slope <- terrain(dem_plot, v = "slope", unit = "radians")
aspect <- terrain(dem_plot, v = "aspect", unit = "radians")
hill <- shade(slope, aspect, angle = 40, direction = 315)
dem_df <- as.data.frame(dem_plot, xy = TRUE, na.rm = TRUE)
hill_df <- as.data.frame(hill, xy = TRUE, na.rm = TRUE)
names(hill_df)[3] <- "Hillshade"
dem_df <- merge(dem_df, hill_df, by = c("x", "y"))

mountain_labels <- data.frame(
  Range = c("Cordillera de Tilarán", "Cordillera Central",
            "Cordillera de Talamanca", "Central Cordillera"),
  Longitude = c(-84.82, -84.05, -83.25, -81.10),
  Latitude = c(10.32, 10.02, 9.20, 8.80),
  Angle = c(-35, -35, -42, -15)
)

site_colors <- c(
  "1" = "#0072B2",
  "2" = "#009E73",
  "3" = "#E69F00",
  "4" = "#CC79A7"
)

world <- ne_countries(scale = "medium", returnclass = "sf")
region <- suppressWarnings(st_crop(world, xmin = -86.2, xmax = -78.9, ymin = 7.7, ymax = 11.5))
country_labels <- data.frame(
  Country = c("COSTA RICA", "PANAMA"),
  Longitude = c(-84.20, -80.75),
  Latitude = c(9.55, 8.15)
)

p <- ggplot() +
  geom_raster(data = dem_df, aes(x, y, fill = Elevation_m)) +
  geom_raster(data = dem_df, aes(x, y, alpha = Hillshade), fill = "black",
              inherit.aes = FALSE) +
  geom_sf(data = region, fill = NA, color = "#4D4D4D", linewidth = 0.42) +
  geom_text(data = mountain_labels,
            aes(Longitude, Latitude, label = Range, angle = Angle),
            color = "#493C2C", size = 2.35, fontface = "italic", alpha = 0.85) +
  geom_text(
    data = country_labels,
    aes(Longitude, Latitude, label = Country),
    color = "#3F3F3F",
    size = 3.0,
    fontface = "bold"
  ) +
  geom_point(
    data = sites,
    aes(Longitude, Latitude, color = factor(SiteID)),
    shape = 21,
    size = 5.5,
    stroke = 1.2,
    fill = "white"
  ) +
  geom_text(
    data = sites,
    aes(Longitude, Latitude, label = SiteID),
    color = "black",
    fontface = "bold",
    size = 3.0
  ) +
  geom_label_repel(
    data = sites,
    aes(Longitude, Latitude, label = MapLabel, color = factor(SiteID)),
    seed = 20260720,
    size = 2.75,
    fontface = "plain",
    box.padding = 0.45,
    point.padding = 0.55,
    min.segment.length = 0,
    segment.color = "#666666",
    segment.size = 0.35,
    label.padding = grid::unit(0.18, "lines"),
    label.r = grid::unit(0.08, "lines"),
    label.size = 0.25,
    fill = "white",
    show.legend = FALSE,
    max.overlaps = Inf
  ) +
  scale_fill_gradientn(
    colors = c("#DCEAF2", "#D9E4C5", "#AFC58B", "#C5AA73", "#92704D", "#F1EEE8"),
    values = scales::rescale(c(0, 100, 500, 1200, 2200, 3700)),
    limits = c(0, 3700), oob = scales::squish,
    name = "Elevation (m)"
  ) +
  scale_alpha(range = c(0.02, 0.25), guide = "none") +
  scale_color_manual(values = site_colors, guide = "none") +
  annotation_scale(
    location = "bl",
    width_hint = 0.17,
    text_cex = 0.72,
    line_width = 0.55,
    pad_x = grid::unit(0.35, "cm"),
    pad_y = grid::unit(0.35, "cm")
  ) +
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    height = grid::unit(0.8, "cm"),
    width = grid::unit(0.8, "cm"),
    pad_x = grid::unit(0.35, "cm"),
    pad_y = grid::unit(0.35, "cm"),
    style = north_arrow_minimal(text_size = 7)
  ) +
  coord_sf(
    xlim = c(-86.0, -79.1),
    ylim = c(7.9, 11.35),
    expand = FALSE,
    datum = st_crs(4326)
  ) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw(base_size = 8.5, base_family = "sans") +
  theme(
    panel.background = element_rect(fill = "#DCEAF2", color = NA),
    panel.grid.major = element_line(color = "white", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold", size = 8),
    axis.text = element_text(color = "black", size = 7),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.border = element_rect(color = "black", linewidth = 0.55),
    plot.margin = margin(6, 8, 5, 6)
  )

png_file <- file.path(figures_root, "sampling_localities_costa_rica_panama_map.png")
pdf_file <- file.path(figures_root, "sampling_localities_costa_rica_panama_map.pdf")
svg_file <- file.path(figures_root, "sampling_localities_costa_rica_panama_map.svg")

ggsave(png_file, p, width = 7.2, height = 4.8, dpi = 600, bg = "white")
ggsave(pdf_file, p, width = 7.2, height = 4.8, device = "pdf")
ggsave(svg_file, p, width = 7.2, height = 4.8, device = svglite::svglite)

print(sites[, c("SiteID", "ShortName", "BiogeographicRegion", "Country",
                "Latitude", "Longitude", "SampleCount")])
