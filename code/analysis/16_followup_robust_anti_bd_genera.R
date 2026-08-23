#!/usr/bin/env Rscript

# Follow-up visualization for anti-Bd genera with a robust ANCOM-BC2 global
# locality effect. This script intentionally does not run post-hoc pairwise
# tests; it displays all samples so the global pattern remains auditable.

suppressPackageStartupMessages({
  library(phyloseq)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
})

set.seed(20260815)
result_dir <- "results/genus_ancombc2_aldex2/robust_anti_bd_followup"
figure_dir <- file.path(result_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

whole_file <- "data/processed/ps_clean_47_rebuilt_from_01_filtered.rds"
anti_file <- "data/processed/blast_97id_90cov/ps_anti_bd_47_rebuilt_97id_90cov.rds"
global_file <- "results/genus_ancombc2_aldex2/anti_bd_ancombc2_global_test.csv"
stopifnot(all(file.exists(c(whole_file, anti_file, global_file))))

site_levels <- c("Lost Iguana", "Veragua", "Altos de Campana", "Soberanía")
site_colors <- c("Lost Iguana" = "#B44682", "Veragua" = "#166AA5",
                 "Altos de Campana" = "#18835C", "Soberanía" = "#D07A00")
site_key <- c(
  "Lost Iguana Resort, La Fortuna, San Carlos" = "Lost Iguana",
  "Veragua Rainforest, Las Brisas de Veragua" = "Veragua",
  "PN Altos de Campana" = "Altos de Campana",
  "PN Soberanía, Las Cumbres" = "Soberanía"
)

otu_taxa_rows <- function(ps) {
  x <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) x <- t(x)
  x
}
aggregate_genus <- function(ps) {
  x <- otu_taxa_rows(ps)
  tx <- as(tax_table(ps), "matrix")[rownames(x), , drop = FALSE]
  genus <- as.character(tx[, "Genus"])
  family <- as.character(tx[, "Family"])
  family[is.na(family) | trimws(family) == ""] <- "Bacteria"
  missing <- is.na(genus) | trimws(genus) == ""
  genus[missing] <- paste0("Unclassified_", family[missing])
  rowsum(x, group = genus, reorder = FALSE)
}

ps_whole <- readRDS(whole_file)
ps_anti <- readRDS(anti_file)
whole <- aggregate_genus(ps_whole)
anti <- aggregate_genus(ps_anti)
sample_ids <- colnames(whole)
anti <- anti[, sample_ids, drop = FALSE]
metadata <- data.frame(sample_data(ps_whole), check.names = FALSE,
                       stringsAsFactors = FALSE)[sample_ids, , drop = FALSE]
metadata$Site <- factor(unname(site_key[metadata$Locality]), levels = site_levels)

global <- read.csv(global_file, check.names = FALSE, stringsAsFactors = FALSE)
robust_genera <- global$taxon[global$diff_robust_abn %in% TRUE]
stopifnot(length(robust_genera) == 7L, all(robust_genera %in% rownames(anti)))

display_name <- function(x) {
  x <- gsub("^Unclassified_", "Uncl. ", x)
  x[x == "Burkholderia-Caballeronia-Paraburkholderia"] <- "Burkholderia complex"
  x[x == "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"] <- "Rhizobium complex"
  x
}

# Relative abundance is the inhibitory reads assigned to a genus divided by
# the complete-library read count, not divided by anti-Bd reads alone.
anti_percent <- 100 * sweep(anti, 2, colSums(whole), "/")
distribution <- do.call(rbind, lapply(robust_genera, function(g) {
  data.frame(Genus_full = g, Genus = display_name(g), Sample_ID = sample_ids,
             Site = metadata$Site,
             Relative_abundance_percent = anti_percent[g, sample_ids],
             Detected = anti[g, sample_ids] > 0,
             stringsAsFactors = FALSE)
}))
distribution$Site <- factor(distribution$Site, levels = site_levels)
distribution$Genus <- factor(distribution$Genus,
                              levels = display_name(robust_genera))
write.csv(distribution, file.path(result_dir, "robust_genera_sample_abundance.csv"),
          row.names = FALSE)

site_summary <- do.call(rbind, lapply(split(distribution,
                                            list(distribution$Genus_full,
                                                 distribution$Site), drop = TRUE),
  function(z) data.frame(
    Genus = z$Genus_full[1], Site = as.character(z$Site[1]), N = nrow(z),
    Detected_n = sum(z$Detected), Prevalence_percent = 100 * mean(z$Detected),
    Mean_percent = mean(z$Relative_abundance_percent),
    Median_percent = median(z$Relative_abundance_percent),
    Q1_percent = unname(quantile(z$Relative_abundance_percent, 0.25)),
    Q3_percent = unname(quantile(z$Relative_abundance_percent, 0.75))
  )))
rownames(site_summary) <- NULL
write.csv(site_summary, file.path(result_dir, "robust_genera_site_summary.csv"),
          row.names = FALSE)

p_distribution <- ggplot(distribution,
                          aes(Site, Relative_abundance_percent, fill = Site)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.35,
              colour = "grey35", linewidth = 0.25) +
  geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white", linewidth = 0.3) +
  geom_jitter(aes(shape = Detected), width = 0.09, height = 0,
              size = 1.05, alpha = 0.78, colour = "#202020") +
  facet_wrap(~Genus, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = site_colors, guide = "none") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), name = "Detected") +
  labs(x = NULL, y = "Relative abundance (% of all reads)",
       title = "Anti-Bd genera with robust global locality effects",
       subtitle = "Boxes show median and IQR; every frog is displayed",
       caption = "Open points denote non-detection. Note the free y-axis scale among genera.") +
  theme_classic(base_size = 8, base_family = "sans") +
  theme(strip.text = element_text(face = "italic", size = 7.2),
        axis.text.x = element_text(angle = 50, hjust = 1, size = 6.2),
        axis.title.y = element_text(face = "bold"), legend.position = "top",
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8),
        plot.caption = element_text(size = 6.5, hjust = 0),
        panel.spacing = grid::unit(0.8, "lines"))

save_pub <- function(plot, stem, width_mm, height_mm) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  ggsave(paste0(stem, ".png"), plot, width = w, height = h, dpi = 600, bg = "white")
  svglite::svglite(paste0(stem, ".svg"), width = w, height = h)
  print(plot); dev.off()
  grDevices::pdf(paste0(stem, ".pdf"), width = w, height = h,
                 family = "Helvetica", useDingbats = FALSE)
  print(plot); dev.off()
  ragg::agg_tiff(paste0(stem, ".tiff"), width = w, height = h,
                 units = "in", res = 600, background = "white", compression = "lzw")
  print(plot); dev.off()
}

save_pub(p_distribution,
         file.path(figure_dir, "figure_robust_anti_bd_genera_abundance"), 183, 132)
global_followup <- global[match(robust_genera, global$taxon),
                          c("taxon", "W", "p_val", "q_val", "passed_ss",
                            "diff_robust_abn")]
write.csv(global_followup,
          file.path(result_dir, "robust_genera_global_ancombc2.csv"), row.names = FALSE)
cat("Robust global genera:", length(robust_genera), "\n")
