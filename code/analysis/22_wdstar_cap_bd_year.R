#!/usr/bin/env Rscript

# Heteroscedasticity-robust Wd* sensitivity analysis and Bd-year constrained
# ordination. This script uses exactly the same genus-level 20% prevalence
# filter and robust Aitchison (rCLR Euclidean) geometry as analysis 21.

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(WdStar)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(svglite)
  library(ragg)
})

set.seed(20260818)

result_dir <- "results/core_locality_abundance_diversity"
figure_dir <- file.path(result_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

full_file <- "data/processed/ps_clean_47_rebuilt_from_01_filtered.rds"
anti_file <- "data/processed/blast_97id_90cov/ps_anti_bd_47_rebuilt_97id_90cov.rds"
wave_file <- "results/wave_exploratory/bd_wave_spatial_model_site_predictions.csv"
stopifnot(file.exists(full_file), file.exists(anti_file), file.exists(wave_file))

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
community_palette <- c("Whole community" = "#365E8D",
                       "Anti-Bd community" = "#B54B88")

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

prevalence_filter <- function(x, fraction = 0.20) {
  x[, colSums(x > 0) >= ceiling(nrow(x) * fraction), drop = FALSE]
}

ps_full <- readRDS(full_file)
ps_anti <- readRDS(anti_file)
full_asv <- otu_samples_rows(ps_full)
anti_asv <- otu_samples_rows(ps_anti)[rownames(full_asv), , drop = FALSE]
md <- data.frame(sample_data(ps_full), check.names = FALSE,
                 stringsAsFactors = FALSE)[rownames(full_asv), , drop = FALSE]
site <- factor(unname(site_key[md$Locality]), levels = site_levels)
stopifnot(nrow(full_asv) == 47L, !anyNA(site))

full_rclr <- decostand(prevalence_filter(aggregate_genus(ps_full, full_asv)),
                       method = "rclr")
anti_rclr <- decostand(prevalence_filter(aggregate_genus(ps_anti, anti_asv)),
                       method = "rclr")
communities <- list("Whole community" = full_rclr,
                    "Anti-Bd community" = anti_rclr)

wave <- read.csv(wave_file, stringsAsFactors = FALSE, check.names = FALSE)
site_year <- setNames(wave$Predicted_DOD_year, wave$Site_short)
bd_year <- unname(site_year[as.character(site)])
stopifnot(!anyNA(bd_year))

# Wd* is an omnibus test of equality of multivariate locations that estimates
# a separate distance variance for every locality. Its statistic/omega2 is not
# numerically interchangeable with PERMANOVA R2.
wd_results <- do.call(rbind, lapply(names(communities), function(nm) {
  D <- dist(communities[[nm]])
  wd <- WdS.test(dm = D, f = site, nrep = 9999)
  ad <- adonis2(D ~ site, permutations = 9999)
  bd <- betadisper(D, site, type = "centroid", bias.adjust = TRUE)
  pt <- permutest(bd, permutations = 9999)
  data.frame(
    Community = nm, Frog_n = attr(D, "Size"), Locality_n = nlevels(site),
    Wd_star = unname(wd$statistic), Wd_omega2 = unname(wd$estimate),
    Wd_p = wd$p.value, Wd_permutations = 9999L,
    PERMANOVA_R2 = unname(ad$R2[1]), PERMANOVA_p = unname(ad$`Pr(>F)`[1]),
    PERMDISP_F = unname(pt$tab[1, "F"]),
    PERMDISP_p = unname(pt$tab[1, "Pr(>F)"]), stringsAsFactors = FALSE
  )
}))
wd_results$Wd_Holm_p <- p.adjust(wd_results$Wd_p, method = "holm")
write.csv(wd_results, file.path(result_dir, "wdstar_locality_sensitivity.csv"),
          row.names = FALSE)

# CAP/dbRDA with a continuous Bd-year constraint. The constrained axis is the
# compositional direction fitted to year; MDS1 is the leading residual axis.
cap_objects <- lapply(communities, function(coord) {
  dbrda(dist(coord) ~ bd_year)
})

cap_tests <- do.call(rbind, lapply(names(cap_objects), function(nm) {
  fit <- cap_objects[[nm]]
  eig <- fit$CCA$eig
  total <- fit$tot.chi
  data.frame(
    Community = nm,
    Constrained_fraction = fit$CCA$tot.chi / total,
    CAP1_fraction_total = eig[1] / total,
    CAP1_fraction_constrained = eig[1] / sum(eig),
    stringsAsFactors = FALSE
  )
}))
write.csv(cap_tests, file.path(result_dir, "cap_bd_year_axis_summary.csv"),
          row.names = FALSE)

# Exact inference for the continuous Bd-year term uses the locality as the
# exchangeable unit. All frogs from one locality receive the same permuted year.
permutations_all <- function(x) {
  if (length(x) == 1L) return(matrix(x, nrow = 1L))
  do.call(rbind, lapply(seq_along(x), function(i) {
    cbind(x[i], permutations_all(x[-i]))
  }))
}
perm_years <- permutations_all(unname(site_year[site_levels]))
year_tests <- do.call(rbind, lapply(names(communities), function(nm) {
  D <- dist(communities[[nm]])
  fit_stat <- function(y) {
    ad <- adonis2(D ~ y, permutations = 0)
    c(R2 = ad$R2[1], pseudo_F = ad$F[1])
  }
  observed <- fit_stat(bd_year)
  perm_F <- apply(perm_years, 1, function(py) {
    y_perm <- py[match(as.character(site), site_levels)]
    fit_stat(y_perm)["pseudo_F"]
  })
  data.frame(
    Community = nm, R2 = unname(observed["R2"]),
    pseudo_F = unname(observed["pseudo_F"]),
    Exact_p = mean(perm_F >= observed["pseudo_F"] - 1e-12),
    Locality_n = 4L, Frog_n = nrow(communities[[nm]]),
    Exact_permutations = nrow(perm_years), stringsAsFactors = FALSE
  )
}))
year_tests$Holm_p <- p.adjust(year_tests$Exact_p, method = "holm")
write.csv(year_tests,
          file.path(result_dir, "robust_aitchison_bd_year_tests.csv"),
          row.names = FALSE)

cap_scores <- do.call(rbind, lapply(names(cap_objects), function(nm) {
  sc <- scores(cap_objects[[nm]], display = "sites", choices = c(1, 2),
               scaling = 1)
  # vegan names the first constrained and first unconstrained axes CAP1/MDS1.
  z <- data.frame(Sample_ID = rownames(sc), CAP1 = sc[, 1], MDS1 = sc[, 2],
                  Site = site, Bd_arrival_year = bd_year, Community = nm)
  z
}))
write.csv(cap_scores, file.path(result_dir, "cap_bd_year_sample_scores.csv"),
          row.names = FALSE)

cap_centroids <- do.call(rbind, lapply(split(cap_scores, cap_scores$Community),
                                      function(z) {
  out <- aggregate(cbind(CAP1, MDS1, Bd_arrival_year) ~ Site, z, mean)
  out$Community <- z$Community[1]
  out[order(out$Bd_arrival_year), ]
}))
write.csv(cap_centroids, file.path(result_dir, "cap_bd_year_locality_centroids.csv"),
          row.names = FALSE)

theme_pub <- theme_classic(base_size = 8, base_family = "sans") +
  theme(axis.line = element_line(linewidth = 0.4),
        plot.title = element_text(size = 9, face = "bold"),
        plot.subtitle = element_text(size = 6.8),
        plot.caption = element_text(size = 5.8, hjust = 0),
        legend.title = element_blank(), legend.text = element_text(size = 6.3))

cap_panel <- function(nm, tag) {
  z <- cap_scores[cap_scores$Community == nm, ]
  cent <- cap_centroids[cap_centroids$Community == nm, ]
  stat <- cap_tests[cap_tests$Community == nm, ]
  ggplot(z, aes(CAP1, MDS1, colour = Site)) +
    geom_point(size = 1.45, alpha = 0.42) +
    geom_path(data = cent, aes(CAP1, MDS1), inherit.aes = FALSE,
              linewidth = 0.85, colour = "#4A4A4A",
              arrow = arrow(length = grid::unit(1.8, "mm"), type = "closed")) +
    geom_point(data = cent, aes(CAP1, MDS1, fill = Site), shape = 21,
               size = 3.4, colour = "white", stroke = 0.7) +
    geom_text_repel(data = cent,
                    aes(CAP1, MDS1,
                        label = paste0(as.character(Site), " ",
                                       round(Bd_arrival_year, 1))),
                    size = 1.85, seed = 20260818, min.segment.length = 0,
                    show.legend = FALSE) +
    scale_colour_manual(values = site_palette) +
    scale_fill_manual(values = site_palette) +
    labs(title = paste0(tag, "  ", nm),
         subtitle = sprintf("Bd year constrains %.1f%% of robust Aitchison variation",
                            100 * stat$Constrained_fraction),
         x = "CAP1: composition constrained by Bd-arrival year",
         y = "MDS1: leading residual composition axis") +
    theme_pub + theme(legend.position = "none")
}

wd_plot_data <- reshape(
  wd_results[, c("Community", "PERMANOVA_p", "PERMDISP_p", "Wd_p")],
  varying = c("PERMANOVA_p", "PERMDISP_p", "Wd_p"),
  v.names = "p_value", timevar = "Method", times = c("PERMANOVA", "PERMDISP", "Wd*"),
  direction = "long"
)
wd_plot_data$Method <- factor(wd_plot_data$Method,
                              levels = c("PERMANOVA", "PERMDISP", "Wd*"))
p_sensitivity <- ggplot(wd_plot_data,
                        aes(Method, -log10(p_value), colour = Community,
                            group = Community)) +
  geom_hline(yintercept = -log10(0.05), linewidth = 0.5, linetype = 2,
             colour = "grey55") +
  geom_point(size = 2.7, position = position_dodge(width = 0.12)) +
  geom_text(aes(label = format.pval(p_value, digits = 2, eps = 1e-4)),
            size = 2.1, vjust = -0.8, position = position_dodge(width = 0.12),
            show.legend = FALSE) +
  scale_colour_manual(values = community_palette) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.18))) +
  labs(title = "c  Sensitivity tests",
       subtitle = "Wd* is robust to unequal dispersion",
       x = NULL, y = expression(-log[10](italic(p)))) +
  theme_pub + theme(axis.text.x = element_text(angle = 20, hjust = 1),
                    legend.position = "none")

p_cap <- (cap_panel("Whole community", "a") |
            cap_panel("Anti-Bd community", "b")) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(
    title = "Bd-year-constrained community composition",
    subtitle = "CAP/dbRDA uses robust Aitchison distances for whole and anti-Bd communities.",
    caption = paste0(
      "Faint points are 47 frogs; large points are locality centroids ordered by the two-dimensional spatial-model estimate of Bd-arrival year.\n",
      "CAP trajectories are descriptive because Bd-arrival year is completely linked to four localities."
    ),
    theme = theme(plot.title = element_text(size = 11, face = "bold"),
                  plot.subtitle = element_text(size = 7.2),
                  plot.caption = element_text(size = 6, hjust = 0))
  )

save_pub <- function(plot, stem, width_mm = 183, height_mm = 108, dpi = 600) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  ggsave(paste0(stem, ".png"), plot, width = w, height = h, dpi = dpi,
         bg = "white")
  svglite(paste0(stem, ".svg"), width = w, height = h); print(plot); dev.off()
  cairo_pdf(paste0(stem, ".pdf"), width = w, height = h,
            family = "Helvetica"); print(plot); dev.off()
  agg_tiff(paste0(stem, ".tiff"), width = w, height = h, units = "in",
           res = dpi, background = "white"); print(plot); dev.off()
}

save_pub(p_cap, file.path(figure_dir, "figure_wdstar_cap_bd_year"),
         height_mm = 100)
# This is the canonical replacement for the former unconstrained PCA trajectory
# figure. Keeping the former filename prevents manuscript links from breaking.
save_pub(p_cap, file.path(figure_dir, "figure_robust_composition_bd_arrival_year"),
         height_mm = 100)

cat("Completed Wd* sensitivity analysis and CAP/dbRDA figure.\n")
print(wd_results)
print(cap_tests)
