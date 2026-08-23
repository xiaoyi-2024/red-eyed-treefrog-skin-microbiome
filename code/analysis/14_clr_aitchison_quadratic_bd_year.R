suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(svglite)
  library(ragg)
})

set.seed(20260816)
result_dir <- "results/bd_arrival_community"
figure_dir <- file.path(result_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

full_file <- "data/processed/ps_clean_47_rebuilt_from_01_filtered.rds"
anti_file <- "data/processed/blast_97id_90cov/ps_anti_bd_47_rebuilt_97id_90cov.rds"
wave_file <- "results/wave_exploratory/site_metrics_single_wave.csv"
stopifnot(file.exists(full_file), file.exists(anti_file), file.exists(wave_file))

otu_samples_rows <- function(ps) {
  x <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) x <- t(x)
  x
}

ps_full <- readRDS(full_file)
ps_anti <- readRDS(anti_file)
full_asv <- otu_samples_rows(ps_full)
anti_asv <- otu_samples_rows(ps_anti)[rownames(full_asv), , drop = FALSE]
md <- data.frame(sample_data(ps_full), check.names = FALSE,
                 stringsAsFactors = FALSE)[rownames(full_asv), , drop = FALSE]

site_key <- c(
  "Lost Iguana Resort, La Fortuna, San Carlos" = "Lost Iguana",
  "Veragua Rainforest, Las Brisas de Veragua" = "Veragua",
  "PN Altos de Campana" = "Altos de Campana",
  "PN Soberanía, Las Cumbres" = "Soberanía"
)
site_levels <- c("Lost Iguana", "Veragua", "Altos de Campana", "Soberanía")
site <- factor(unname(site_key[md$Locality]), levels = site_levels)
stopifnot(!anyNA(site))

wave <- read.csv(wave_file, stringsAsFactors = FALSE, check.names = FALSE)
site_year <- setNames(wave$pred_wave_year, wave$Site)
year <- unname(site_year[as.character(site)])
year10 <- (year - mean(unname(site_year[site_levels]))) / 10

tax <- as(tax_table(ps_anti), "matrix")[colnames(anti_asv), , drop = FALSE]
genus <- as.character(tax[, "Genus"])
family <- as.character(tax[, "Family"])
family[is.na(family) | trimws(family) == ""] <- "Bacteria"
missing_genus <- is.na(genus) | trimws(genus) == ""
genus[missing_genus] <- paste0("Unclassified_", family[missing_genus])
anti_genus <- t(rowsum(t(anti_asv), group = genus, reorder = FALSE))

clr <- function(x) {
  z <- log(x + 1)
  z - rowMeans(z)
}
coordinates <- list(
  "Whole community" = clr(full_asv),
  "Anti-Bd genera" = clr(anti_genus)
)

permutations <- function(x) {
  if (length(x) == 1L) return(matrix(x, nrow = 1L))
  do.call(rbind, lapply(seq_along(x), function(i) cbind(x[i], permutations(x[-i]))))
}
perm_years <- permutations(unname(site_year[site_levels]))

fit_stats <- function(D, y10) {
  overall <- vegan::adonis2(D ~ y10 + I(y10^2), permutations = 0)
  terms <- vegan::adonis2(D ~ y10 + I(y10^2), permutations = 0, by = "terms")
  linear <- vegan::adonis2(D ~ y10, permutations = 0)
  c(
    linear_R2 = linear$R2[1], linear_F = linear$F[1],
    quadratic_model_R2 = overall$R2[1], quadratic_model_F = overall$F[1],
    quadratic_increment_R2 = terms$R2[2], quadratic_increment_F = terms$F[2]
  )
}

quadratic_results <- list()
permutation_results <- list()
for (community in names(coordinates)) {
  D <- dist(coordinates[[community]])
  observed <- fit_stats(D, year10)
  perm_stats <- t(apply(perm_years, 1, function(py) {
    yp <- py[match(as.character(site), site_levels)]
    yp10 <- (yp - mean(unname(site_year[site_levels]))) / 10
    fit_stats(D, yp10)
  }))
  p_overall <- mean(perm_stats[, "quadratic_model_F"] >=
                      observed["quadratic_model_F"] - 1e-12)
  p_increment <- mean(perm_stats[, "quadratic_increment_F"] >=
                        observed["quadratic_increment_F"] - 1e-12)
  quadratic_results[[community]] <- data.frame(
    Community = community,
    Linear_R2 = unname(observed["linear_R2"]),
    Linear_pseudo_F = unname(observed["linear_F"]),
    Quadratic_model_R2 = unname(observed["quadratic_model_R2"]),
    Quadratic_model_pseudo_F = unname(observed["quadratic_model_F"]),
    Quadratic_model_exact_p = p_overall,
    Quadratic_increment_R2 = unname(observed["quadratic_increment_R2"]),
    Quadratic_increment_pseudo_F = unname(observed["quadratic_increment_F"]),
    Quadratic_increment_exact_p = p_increment,
    Locality_n = 4L, Frog_n = nrow(full_asv), Exact_permutations = 24L,
    Analysis_role = "Exploratory quadratic sensitivity analysis",
    stringsAsFactors = FALSE
  )
  permutation_results[[community]] <- data.frame(
    Community = community, Permutation = seq_len(nrow(perm_stats)), perm_stats,
    stringsAsFactors = FALSE
  )
}
quadratic_tests <- do.call(rbind, quadratic_results)
quadratic_tests$Holm_p_overall <- p.adjust(quadratic_tests$Quadratic_model_exact_p,
                                           method = "holm")
quadratic_tests$Holm_p_quadratic_increment <- p.adjust(
  quadratic_tests$Quadratic_increment_exact_p, method = "holm")
quadratic_tests$Conclusion <- ifelse(
  quadratic_tests$Holm_p_overall < 0.05 &
    quadratic_tests$Holm_p_quadratic_increment < 0.05,
  "Evidence supporting a curved association",
  "No multiplicity-adjusted evidence supporting a curved association"
)
write.csv(quadratic_tests,
          file.path(result_dir, "clr_aitchison_quadratic_bd_year_tests.csv"),
          row.names = FALSE)
write.csv(do.call(rbind, permutation_results),
          file.path(result_dir, "clr_aitchison_quadratic_exact_permutations.csv"),
          row.names = FALSE)

# Descriptive projection of the fitted multivariate trajectory. The curve is
# fitted to four locality centroids on CLR-PC1 and CLR-PC2 separately and is not
# an additional inferential test.
palette_site <- c("Lost Iguana" = "#B44682", "Veragua" = "#166AA5",
                  "Altos de Campana" = "#18835C", "Soberanía" = "#D07A00")
trajectory_panel <- function(coord, community, tag) {
  pc <- prcomp(coord, center = FALSE, scale. = FALSE)
  variance <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  scores <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2], Site = site, Year = year)
  cent <- aggregate(cbind(PC1, PC2, Year) ~ Site, scores, mean)
  cent <- cent[order(cent$Year), ]
  fit1 <- lm(PC1 ~ Year + I(Year^2), data = cent)
  fit2 <- lm(PC2 ~ Year + I(Year^2), data = cent)
  curve <- data.frame(Year = seq(min(cent$Year), max(cent$Year), length.out = 300))
  curve$PC1 <- predict(fit1, curve)
  curve$PC2 <- predict(fit2, curve)
  hit <- quadratic_tests[quadratic_tests$Community == community, ]
  ggplot(cent, aes(PC1, PC2, colour = Site)) +
    geom_path(data = curve, aes(PC1, PC2), inherit.aes = FALSE,
              colour = "#343434", linewidth = 0.9,
              arrow = arrow(length = grid::unit(1.8, "mm"), type = "closed")) +
    geom_point(size = 3.2) +
    geom_text_repel(aes(label = sprintf("%.1f", Year)), size = 2.2,
                    seed = 20260816, min.segment.length = 0,
                    show.legend = FALSE) +
    scale_colour_manual(values = palette_site) +
    labs(
      title = paste0(tag, "  ", community),
      subtitle = sprintf("Full curve: R² = %.3f, p = %.3f\nQuadratic addition: R² = %.3f, p = %.3f",
                         hit$Quadratic_model_R2, hit$Quadratic_model_exact_p,
                         hit$Quadratic_increment_R2, hit$Quadratic_increment_exact_p),
      x = sprintf("CLR-PC1 (%.1f%%)", variance[1]),
      y = sprintf("CLR-PC2 (%.1f%%)", variance[2])
    ) +
    theme_classic(base_size = 8, base_family = "sans") +
    theme(axis.line = element_line(linewidth = 0.4),
          axis.ticks = element_line(linewidth = 0.4),
          plot.title = element_text(size = 9, face = "bold"),
          plot.subtitle = element_text(size = 6.3, lineheight = 1.05),
          legend.position = "none")
}

p <- trajectory_panel(coordinates[["Whole community"]], "Whole community", "a") |
  trajectory_panel(coordinates[["Anti-Bd genera"]], "Anti-Bd genera", "b")
p <- p + plot_annotation(
  title = "Exploratory curved CLR-Aitchison trajectories across Bd-arrival years",
  subtitle = "Numbers label estimated arrival years; arrows run from earlier to later fitted values.",
  caption = paste0("Curves are quadratic fits to four locality centroids projected onto CLR-PC1 and CLR-PC2.\n",
                   "Inference uses full CLR-Aitchison distances and all 24 locality-level year permutations; curves are descriptive."),
  theme = theme(plot.title = element_text(size = 11, face = "bold"),
                plot.subtitle = element_text(size = 7.2),
                plot.caption = element_text(size = 6.1, hjust = 0)))

stem <- file.path(figure_dir, "figure_clr_aitchison_quadratic_bd_year")
w <- 183 / 25.4; h <- 105 / 25.4
ggsave(paste0(stem, ".png"), p, width = w, height = h, dpi = 600, bg = "white")
svglite::svglite(paste0(stem, ".svg"), width = w, height = h); print(p); dev.off()
grDevices::pdf(paste0(stem, ".pdf"), width = w, height = h,
               family = "Helvetica", useDingbats = FALSE); print(p); dev.off()
ragg::agg_tiff(paste0(stem, ".tiff"), width = w, height = h,
               units = "in", res = 600, background = "white"); print(p); dev.off()

cat("Exploratory CLR-Aitchison quadratic tests\n")
print(quadratic_tests)
