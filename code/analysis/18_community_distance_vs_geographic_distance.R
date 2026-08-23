#!/usr/bin/env Rscript

# Test whether between-locality microbiome distances increase with geographic
# separation. Six locality pairs are descriptive/non-independent; inference
# uses exact Mantel permutations of the four locality labels (4! = 24).

suppressPackageStartupMessages({
  library(ggplot2)
  library(svglite)
  library(ragg)
})

result_dir <- "results/community_distance_vs_geography"
figure_dir <- file.path(result_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

site_file <- "results/spatial_elevation_abundance_diversity/site_spatial_metrics.csv"
sites_df <- read.csv(site_file, check.names = FALSE, stringsAsFactors = FALSE)
sites <- c("Lost", "Veragua", "Campana", "Soberanía")
site_alias <- c("Lost Iguana" = "Lost", "Altos de Campana" = "Campana",
                "Veragua" = "Veragua", "Soberanía" = "Soberanía",
                "Lost" = "Lost", "Campana" = "Campana")
sites_df$Site_short <- unname(site_alias[sites_df$Site])
sites_df <- sites_df[match(sites, sites_df$Site_short), ]
stopifnot(!anyNA(sites_df$Site_short))

haversine_km <- function(lon1, lat1, lon2, lat2) {
  rad <- pi / 180
  dlon <- (lon2 - lon1) * rad; dlat <- (lat2 - lat1) * rad
  a <- sin(dlat / 2)^2 + cos(lat1 * rad) * cos(lat2 * rad) * sin(dlon / 2)^2
  6371.0088 * 2 * atan2(sqrt(a), sqrt(1 - a))
}

n <- length(sites)
geo <- matrix(0, n, n, dimnames = list(sites, sites))
for (i in seq_len(n - 1L)) for (j in (i + 1L):n) {
  geo[i, j] <- geo[j, i] <- haversine_km(
    sites_df$Longitude[i], sites_df$Latitude[i],
    sites_df$Longitude[j], sites_df$Latitude[j]
  )
}

read_bray <- function(path, community) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  data.frame(
    Community = community, Metric = "Bray-Curtis",
    Site_1 = unname(site_alias[x$Site_1]),
    Site_2 = unname(site_alias[x$Site_2]),
    Community_distance = x$Mean_Bray_Curtis_distance,
    stringsAsFactors = FALSE
  )
}

read_aitch <- function(path, community) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  data.frame(
    Community = community,
    Metric = ifelse(x$Metric == "CLR-Aitchison", "CLR-Aitchison", "Robust Aitchison"),
    Site_1 = unname(site_alias[x$Site_1]), Site_2 = unname(site_alias[x$Site_2]),
    Community_distance = x$Mean_distance,
    stringsAsFactors = FALSE
  )
}

dat <- rbind(
  read_bray("results/pairwise_bray_distance/pairwise_locality_bray_permanova.csv",
            "Whole community"),
  read_bray("results/anti_bd_pairwise_bray_distance/pairwise_locality_bray_permanova.csv",
            "Anti-Bd"),
  read_aitch("results/pairwise_aitchison_comparison/pairwise_distance_summary.csv",
             "Whole community"),
  read_aitch("results/anti_bd_pairwise_aitchison_comparison/pairwise_distance_summary.csv",
             "Anti-Bd")
)
dat$Geographic_km <- mapply(function(a, b) geo[a, b], dat$Site_1, dat$Site_2)
dat$Pair <- paste(dat$Site_1, dat$Site_2, sep = " vs ")
dat$Community <- factor(dat$Community, levels = c("Whole community", "Anti-Bd"))
dat$Metric <- factor(dat$Metric,
                     levels = c("Bray-Curtis", "CLR-Aitchison", "Robust Aitchison"))
write.csv(dat, file.path(result_dir, "community_and_geographic_pair_distances.csv"),
          row.names = FALSE)

perm_grid <- as.matrix(expand.grid(rep(list(seq_len(n)), n)))
perms <- perm_grid[apply(perm_grid, 1, function(v) length(unique(v)) == n), , drop = FALSE]
lower <- lower.tri(geo)

test_one <- function(d) {
  comm <- matrix(0, n, n, dimnames = list(sites, sites))
  for (i in seq_len(nrow(d))) {
    comm[d$Site_1[i], d$Site_2[i]] <- comm[d$Site_2[i], d$Site_1[i]] <-
      d$Community_distance[i]
  }
  obs_s <- suppressWarnings(cor(geo[lower], comm[lower], method = "spearman"))
  obs_p <- suppressWarnings(cor(geo[lower], comm[lower], method = "pearson"))
  null_s <- apply(perms, 1, function(idx)
    suppressWarnings(cor(geo[lower], comm[idx, idx][lower], method = "spearman")))
  null_p <- apply(perms, 1, function(idx)
    suppressWarnings(cor(geo[lower], comm[idx, idx][lower], method = "pearson")))
  lm_fit <- lm(Community_distance ~ I(Geographic_km / 100), data = d)
  data.frame(
    Mantel_Spearman_rho = obs_s,
    Exact_p_Spearman = mean(abs(null_s) >= abs(obs_s) - 1e-12),
    Mantel_Pearson_r = obs_p,
    Exact_p_Pearson = mean(abs(null_p) >= abs(obs_p) - 1e-12),
    Descriptive_slope_per_100km = unname(coef(lm_fit)[2]),
    Descriptive_linear_R2 = summary(lm_fit)$r.squared,
    Localities = n, Locality_pairs = sum(lower), Permutations = nrow(perms)
  )
}

split_dat <- split(dat, interaction(dat$Community, dat$Metric, drop = TRUE))
tests <- do.call(rbind, lapply(split_dat, function(d) {
  cbind(data.frame(Community = as.character(d$Community[1]),
                   Metric = as.character(d$Metric[1])), test_one(d))
}))
rownames(tests) <- NULL
tests$Holm_p_Spearman <- p.adjust(tests$Exact_p_Spearman, method = "holm")
tests$Holm_p_Pearson_sensitivity <- p.adjust(tests$Exact_p_Pearson, method = "holm")
tests <- tests[order(match(tests$Community, c("Whole community", "Anti-Bd")),
                     match(tests$Metric, levels(dat$Metric))), ]
write.csv(tests, file.path(result_dir, "exact_mantel_tests.csv"), row.names = FALSE)

ann <- tests
ann$Community <- factor(ann$Community, levels = levels(dat$Community))
ann$Metric <- factor(ann$Metric, levels = levels(dat$Metric))
ann$Label <- sprintf("rho = %.2f; exact p = %.3f\nHolm p = %.3f",
                     ann$Mantel_Spearman_rho, ann$Exact_p_Spearman,
                     ann$Holm_p_Spearman)

p <- ggplot(dat, aes(Geographic_km, Community_distance)) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.65, colour = "grey38") +
  geom_point(aes(fill = Pair), shape = 21, size = 2.7, colour = "black", stroke = 0.3) +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = Label),
            inherit.aes = FALSE, hjust = -0.06, vjust = 1.15, size = 2.3) +
  facet_wrap(vars(Community, Metric), ncol = 3, scales = "free_y") +
  scale_fill_brewer(palette = "Dark2", guide = guide_legend(nrow = 2)) +
  labs(
    title = "Microbiome compositional distance versus geographic separation",
    subtitle = "Top: whole community; bottom: anti-Bd community; six locality pairs per panel",
    x = "Geographic distance between localities (km)", y = "Mean community distance",
    fill = "Locality pair",
    caption = paste0("Six locality-pair points and fitted lines are descriptive/non-independent. Exact Mantel-Spearman tests use all 24 locality-label permutations;\n",
                     "Holm correction covers six tests. Aitchison magnitudes are not compared directly between whole and anti-Bd because feature dimensions differ.")
  ) +
  theme_classic(base_size = 8) +
  theme(plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 7), strip.background = element_rect(fill = "grey94"),
        strip.text = element_text(face = "bold", size = 7), legend.position = "bottom",
        legend.text = element_text(size = 6), plot.caption = element_text(size = 5.8, hjust = 0),
        panel.spacing = grid::unit(7, "mm"))

save_plot <- function(plot, stem, width_mm = 183, height_mm = 145) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  ggsave(paste0(stem, ".png"), plot, width = w, height = h, dpi = 600, bg = "white")
  svglite::svglite(paste0(stem, ".svg"), width = w, height = h); print(plot); dev.off()
  pdf(paste0(stem, ".pdf"), width = w, height = h, useDingbats = FALSE); print(plot); dev.off()
  ragg::agg_tiff(paste0(stem, ".tiff"), width = w, height = h, units = "in",
                 res = 600, compression = "lzw", background = "white")
  print(plot); dev.off()
}
save_plot(p, file.path(figure_dir, "figure_three_community_distances_vs_geography"))

write.csv(data.frame(
  Item = c("Independent unit", "Geographic separation", "Primary test",
           "Exact null", "Multiplicity", "Plot interpretation"),
  Definition = c("Locality (n = 4)", "Great-circle distance from mean coordinates",
                 "Mantel-Spearman correlation", "All 4! = 24 locality-label permutations",
                 "Holm correction across 3 distances x 2 communities",
                 "Six pair points and linear fit are descriptive/non-independent")
), file.path(result_dir, "analysis_methods.csv"), row.names = FALSE)

print(tests)
