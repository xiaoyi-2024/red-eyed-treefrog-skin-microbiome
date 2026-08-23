library(sf)
library(rnaturalearth)
library(ggplot2)
library(ggrepel)
library(ggspatial)
library(png)
source("code/analysis/common.R")
figures_root <- "results/maps/figures"
dir.create(figures_root, recursive = TRUE, showWarnings = FALSE)

output_dir <- "results/maps/tables"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

original_sites <- read.csv("data/sampling_localities.csv", stringsAsFactors = FALSE)

# Counts for the four focal microbiome sites are derived from the final object.
original_sites$ShortName <- c(
  "Lost Iguana Resort", "Veragua Rainforest",
  "PN Altos de Campana", "PN Soberanía"
)
original_sites$Country_position <- c(
  "Northeastern Costa Rica", "Southeastern Costa Rica",
  "Central Panama", "Central Panama"
)
original_sites$Phenotype <- factor(
  c("Northeastern Costa Rica", "Southeastern Costa Rica",
    "Central Panama", "Central Panama"),
  levels = c("Northeastern Costa Rica", "Southeastern Costa Rica", "Central Panama")
)
final_counts <- table(site_factor(ps_clean))
original_sites$SampleCount <- as.integer(final_counts[original_sites$ShortName])
stopifnot(sum(original_sites$SampleCount) == 47L)
original_sites$Label <- sprintf("%d  %s (n=%d)\n%s", original_sites$SiteID,
                                original_sites$ShortName, original_sites$SampleCount,
                                original_sites$Country_position)

phenotype_colors <- c(
  "Northeastern Costa Rica" = "#0072B2",
  "Southeastern Costa Rica" = "#CC79A7",
  "Central Panama" = "#E69F00"
)

world <- ne_countries(scale = "medium", returnclass = "sf")
region_map <- suppressWarnings(st_crop(world, xmin = -86.2, xmax = -78.9,
                                       ymin = 7.7, ymax = 11.2))
country_labels <- data.frame(
  label = c("COSTA RICA (CR)", "PANAMA (PA)"),
  Longitude = c(-85.0, -80.2), Latitude = c(10.85, 7.98)
)

# Leg-colour phenotype images are assigned from the supplied reference map.
phenotype_dir <- "data/reference/leg_colour_phenotypes"
phenotype_boxes <- data.frame(
  Phenotype = c("Northeastern Costa Rica", "Southeastern Costa Rica", "Central Panama"),
  File = file.path(phenotype_dir, c(
    "northeastern_cr.png", "southeastern_cr_pa.png", "central_panama.png"
  )),
  xmin = c(-85.45, -83.30, -80.90),
  xmax = c(-84.85, -82.70, -80.30),
  ymin = c(10.25, 10.25, 9.75),
  ymax = c(10.73, 10.73, 10.23),
  stringsAsFactors = FALSE
)
stopifnot(all(file.exists(phenotype_boxes$File)))

phenotype_layers <- lapply(seq_len(nrow(phenotype_boxes)), function(i) {
  b <- phenotype_boxes[i, ]
  annotation_custom(
    grid::rasterGrob(readPNG(b$File), interpolate = TRUE),
    xmin = b$xmin, xmax = b$xmax, ymin = b$ymin, ymax = b$ymax
  )
})

p <- ggplot() +
  geom_sf(data = region_map, fill = "#F3F0E9", color = "#696969", linewidth = 0.35) +
  phenotype_layers +
  geom_text(data = country_labels, aes(Longitude, Latitude, label = label),
            size = 3.0, fontface = "bold", color = "#666666") +
  geom_point(data = original_sites,
             aes(Longitude, Latitude, fill = Phenotype), shape = 23, size = 5.2,
             stroke = 0.9, color = "black") +
  geom_text(data = original_sites,
            aes(Longitude, Latitude, label = SiteID),
            size = 2.7, fontface = "bold", color = "black") +
  geom_label_repel(data = original_sites,
                   aes(Longitude, Latitude, label = Label),
                   seed = analysis_seed + 1L, size = 2.65,
                   box.padding = 0.45, point.padding = 0.55,
                   min.segment.length = 0, segment.color = "black",
                   label.size = 0.25, fill = "white", max.overlaps = Inf) +
  scale_fill_manual(values = phenotype_colors, drop = FALSE) +
  annotation_scale(location = "bl", width_hint = 0.16, text_cex = 0.7,
                   pad_x = grid::unit(0.35, "cm"), pad_y = grid::unit(0.35, "cm")) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         height = grid::unit(0.75, "cm"), width = grid::unit(0.75, "cm"),
                         pad_x = grid::unit(0.35, "cm"), pad_y = grid::unit(0.35, "cm"),
                         style = north_arrow_minimal(text_size = 7)) +
  coord_sf(xlim = c(-86.0, -79.1), ylim = c(7.85, 11.05),
           expand = FALSE, datum = st_crs(4326)) +
  labs(x = "Longitude", y = "Latitude", fill = "Biogeographic region",
       shape = NULL) +
  guides(fill = guide_legend(override.aes = list(shape = 23, size = 3.5, color = "black"))) +
  theme_bw(base_size = 8.5, base_family = "sans") +
  theme(
    panel.background = element_rect(fill = "#DCEAF2", color = NA),
    panel.grid.major = element_line(color = "white", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold", size = 8),
    axis.text = element_text(color = "black", size = 7),
    legend.position = "right", legend.title = element_text(face = "bold"),
    panel.border = element_rect(color = "black", linewidth = 0.55),
    plot.margin = margin(6, 8, 5, 6)
  )

stem <- "sampling_map_four_focal_sites_and_five_regions"
ggsave(file.path(figures_root, paste0(stem, ".png")), p,
       width = 9.2, height = 6.0, dpi = 600, bg = "white")
svg_file <- file.path(figures_root, paste0(stem, ".svg"))
pdf_file <- file.path(figures_root, paste0(stem, ".pdf"))
ggsave(svg_file, p,
       width = 9.2, height = 6.0, device = svglite::svglite)
# Converting the completed SVG avoids the XQuartz/Cairo dependency that can
# make raster-containing PDF devices fail on a clean macOS R installation.
if (requireNamespace("rsvg", quietly = TRUE)) {
  rsvg::rsvg_pdf(svg_file, pdf_file)
} else if (nzchar(Sys.which("sips"))) {
  status <- system2("sips", c("-s", "format", "pdf",
                               file.path(figures_root, paste0(stem, ".png")),
                               "--out", pdf_file), stdout = FALSE, stderr = FALSE)
  if (status != 0L) warning("PDF conversion with sips failed; PNG and SVG were still created.")
} else {
  warning("PDF conversion skipped: install rsvg; PNG and SVG were still created.")
}

write.csv(original_sites,
          file.path(output_dir, "four_focal_microbiome_sampling_sites.csv"),
          row.names = FALSE)
write.csv(phenotype_boxes,
          file.path(output_dir, "five_region_leg_colour_phenotype_mapping.csv"),
          row.names = FALSE)
print(table(original_sites$Phenotype))
