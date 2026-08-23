#!/usr/bin/env Rscript

# Integrated visualization of the three complementary anti-Bd genus analyses:
# KSM distributional screening, ANCOM-BC2 primary differential abundance, and
# ALDEx2 Monte Carlo CLR sensitivity analysis. Also creates a nine-genus
# sample-level violin figure and its manuscript-facing source table.

suppressPackageStartupMessages({
  library(phyloseq)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
})

set.seed(20260816)
result_dir <- Sys.getenv(
  "ANTI_BD_DRIVER_RESULT_DIR",
  unset = "results/ksm_locality_driver_screening"
)
figure_dir <- file.path(result_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

is_10pct <- grepl("10pct", basename(result_dir), fixed = TRUE)
rank_file <- file.path(
  result_dir,
  if (is_10pct) "anti_bd_10pct_three_method_ranked.csv" else
    "anti_bd_genus_ksm_ranked.csv"
)
prevalence_threshold <- if (is_10pct) 0.10 else 0.20
minimum_detected_samples <- ceiling(47 * prevalence_threshold)
whole_file <- "data/processed/ps_clean_47_rebuilt_from_01_filtered.rds"
anti_file <- "data/processed/blast_97id_90cov/ps_anti_bd_47_rebuilt_97id_90cov.rds"
stopifnot(all(file.exists(c(rank_file, whole_file, anti_file))))

site_levels <- c("Lost Iguana", "Veragua", "Altos de Campana", "Soberanía")
site_colors <- c("Lost Iguana" = "#B44682", "Veragua" = "#166AA5",
                 "Altos de Campana" = "#18835C", "Soberanía" = "#D07A00")
site_key <- c(
  "Lost Iguana Resort, La Fortuna, San Carlos" = "Lost Iguana",
  "Veragua Rainforest, Las Brisas de Veragua" = "Veragua",
  "PN Altos de Campana" = "Altos de Campana",
  "PN Soberanía, Las Cumbres" = "Soberanía"
)

display_name <- function(x) {
  x <- gsub("^Unclassified_", "Uncl. ", x)
  x[x == "Burkholderia-Caballeronia-Paraburkholderia"] <- "Burkholderia complex"
  x
}

save_pub <- function(plot, stem, width_mm, height_mm) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  ggsave(paste0(stem, ".png"), plot, width = w, height = h,
         dpi = 600, bg = "white")
  svglite::svglite(paste0(stem, ".svg"), width = w, height = h)
  print(plot); dev.off()
  grDevices::pdf(paste0(stem, ".pdf"), width = w, height = h,
                 family = "Helvetica", useDingbats = FALSE)
  print(plot); dev.off()
  ragg::agg_tiff(paste0(stem, ".tiff"), width = w, height = h,
                 units = "in", res = 600, background = "white",
                 compression = "lzw")
  print(plot); dev.off()
}

ranked <- read.csv(rank_file, stringsAsFactors = FALSE, check.names = FALSE)
if (!"Cross_method_support" %in% names(ranked)) {
  ranked$Cross_method_support <- ranked$Evidence_class
}
ranked <- ranked[order(ranked$KSM_rank), ]
top20 <- head(ranked, 20L)
top20$Display_genus <- display_name(top20$Genus)
top20$Display_genus <- factor(top20$Display_genus,
                              levels = rev(top20$Display_genus))

# Panel a: KSM is a ranking measure. Bootstrap intervals are uncertainty in the
# observed ranking score and are not null-hypothesis confidence intervals.
p_ksm <- ggplot(top20, aes(KSM, Display_genus)) +
  geom_errorbar(aes(xmin = Bootstrap_CI_low, xmax = Bootstrap_CI_high),
                orientation = "y", width = 0, colour = "grey68", linewidth = 0.45) +
  geom_point(aes(size = Top10_frequency), shape = 21, fill = "#0072B2",
             colour = "black", stroke = 0.35) +
  scale_size_continuous(range = c(1.8, 4.8), limits = c(0, 1),
                        labels = scales::percent_format(accuracy = 1)) +
  labs(title = "a  KSM distributional screening", x = "KSM", y = NULL,
       size = "Bootstrap\ntop-10 frequency") +
  theme_classic(base_size = 7.6) +
  theme(plot.title = element_text(face = "bold", size = 9),
        axis.text.y = element_text(face = "italic", colour = "black", size = 6.4),
        legend.position = "none")

# Panel b: ANCOM-BC2 robust result incorporates its pseudo-count sensitivity
# assessment. Taxa without an estimable global q-value are shown explicitly.
top20$ANCOM_plot <- ifelse(is.na(top20$ANCOMBC2_q), 0,
                           -log10(pmax(top20$ANCOMBC2_q, 1e-300)))
top20$ANCOM_status <- ifelse(
  is.na(top20$ANCOMBC2_q), "Not estimated",
  ifelse(top20$ANCOMBC2_robust_FDR05, "Robust FDR < 0.05",
         ifelse(top20$ANCOMBC2_q < 0.05, "q < 0.05, not robust", "q >= 0.05"))
)
ancom_colors <- c("Robust FDR < 0.05" = "#D55E00",
                  "q < 0.05, not robust" = "#E69F00",
                  "q >= 0.05" = "#8A8A8A", "Not estimated" = "#FFFFFF")
ancom_shapes <- c("Robust FDR < 0.05" = 21, "q < 0.05, not robust" = 21,
                  "q >= 0.05" = 21, "Not estimated" = 4)
p_ancom <- ggplot(top20, aes(ANCOM_plot, Display_genus)) +
  geom_vline(xintercept = -log10(0.05), linetype = 2, colour = "grey45",
             linewidth = 0.4) +
  geom_point(aes(fill = ANCOM_status, shape = ANCOM_status), size = 2.7,
             colour = "black", stroke = 0.55) +
  scale_fill_manual(values = ancom_colors, drop = FALSE) +
  scale_shape_manual(values = ancom_shapes, drop = FALSE) +
  labs(title = "b  ANCOM-BC2", x = expression(-log[10](q)), y = NULL,
       fill = NULL, shape = NULL) +
  theme_classic(base_size = 7.6) +
  theme(plot.title = element_text(face = "bold", size = 9),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "none")

# Panel c: ALDEx2 global Kruskal-Wallis q-values from Monte Carlo CLR instances.
top20$ALDEx_plot <- -log10(pmax(top20$ALDEx2_q, 1e-300))
top20$ALDEx_status <- ifelse(top20$ALDEx2_FDR05, "FDR < 0.05", "q >= 0.05")
p_aldex <- ggplot(top20, aes(ALDEx_plot, Display_genus)) +
  geom_vline(xintercept = -log10(0.05), linetype = 2, colour = "grey45",
             linewidth = 0.4) +
  geom_point(aes(fill = ALDEx_status), shape = 21, size = 2.7,
             colour = "black", stroke = 0.45) +
  scale_fill_manual(values = c("FDR < 0.05" = "#009E73",
                               "q >= 0.05" = "#8A8A8A")) +
  labs(title = "c  ALDEx2", x = expression(-log[10](q)), y = NULL,
       fill = NULL) +
  theme_classic(base_size = 7.6) +
  theme(plot.title = element_text(face = "bold", size = 9),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "none")

three_panel <- p_ksm + p_ancom + p_aldex +
  plot_layout(widths = c(1.45, 1, 1), guides = "keep") +
  plot_annotation(
    title = "Three-method identification of locality-associated anti-Bd genera",
    subtitle = if (is_10pct) paste0(
      "Genus-level anti-Bd community; four localities; 47 frogs.\n",
      "KSM ranks distributional separation; ANCOM-BC2 is primary inference; ",
      "ALDEx2 is sensitivity analysis."
    ) else paste0(
      sprintf("Genus level; %.0f%% prevalence filter (>=%d/47 frogs).\n",
              100 * prevalence_threshold, minimum_detected_samples),
      "KSM ranks distributional separation; ANCOM-BC2 is primary inference; ",
      "ALDEx2 is sensitivity analysis."
    ),
    caption = paste0(
      "KSM: bootstrap 95% intervals; point size = top-10 stability. ",
      "Dashed: q = 0.05; orange/green: significant; grey: not significant; crosses: not estimated."
    ),
    theme = theme(plot.title = element_text(face = "bold", size = 11),
                  plot.subtitle = element_text(size = 7.5),
                  plot.caption = element_text(size = 6.2, hjust = 0))
  )

write.csv(top20, file.path(result_dir,
                           "figure_three_method_anti_bd_source_data.csv"),
          row.names = FALSE)
save_pub(three_panel,
         file.path(figure_dir, "figure_anti_bd_three_method_driver_comparison"),
         183, 165)

# Nine-genera rule: within the selected-prevalence anti-Bd feature set, retain genera
# supported by robust ANCOM-BC2 and/or ALDEx2 FDR < 0.05, rank by observed KSM,
# and take the first nine. This rule is fixed before plotting sample abundances.
eligible <- ranked[ranked$ANCOMBC2_robust_FDR05 | ranked$ALDEx2_FDR05, ]
top9 <- head(eligible[order(eligible$KSM_rank), ], 9L)
stopifnot(nrow(top9) == 9L)

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
anti_percent <- 100 * sweep(anti, 2, colSums(whole), "/")

distribution <- do.call(rbind, lapply(top9$Genus, function(g) {
  data.frame(Genus = g, Display_genus = display_name(g),
             Sample_ID = sample_ids, Site = metadata$Site,
             Relative_abundance_percent = anti_percent[g, sample_ids],
             Detected = anti[g, sample_ids] > 0, stringsAsFactors = FALSE)
}))
distribution$Site <- factor(distribution$Site, levels = site_levels)
distribution$Display_genus <- factor(distribution$Display_genus,
                                      levels = display_name(top9$Genus))
write.csv(distribution,
          file.path(result_dir, "top9_anti_bd_genera_sample_abundance.csv"),
          row.names = FALSE)

top9_table <- top9[, c("KSM_rank", "Genus", "KSM", "Bootstrap_CI_low",
                       "Bootstrap_CI_high", "Top10_frequency", "Prevalence",
                       "Mean_relative_abundance_percent", "ANCOMBC2_q",
                       "ANCOMBC2_robust_FDR05", "ALDEx2_q", "ALDEx2_FDR05",
                       "Cross_method_support")]
write.csv(top9_table,
          file.path(result_dir, "table_top9_anti_bd_driver_genera.csv"),
          row.names = FALSE)

p_violin <- ggplot(distribution,
                   aes(Site, Relative_abundance_percent, fill = Site)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.34,
              colour = "grey35", linewidth = 0.25) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", linewidth = 0.3) +
  geom_jitter(aes(shape = Detected), width = 0.085, height = 0,
              size = 0.9, alpha = 0.78, colour = "#202020") +
  facet_wrap(~Display_genus, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = site_colors, guide = "none") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), name = "Detection",
                     labels = c(`TRUE` = "Detected", `FALSE` = "Not detected")) +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.01),
                     expand = expansion(mult = c(0.02, 0.08))) +
  labs(
    title = "Sample-level distributions of nine principal anti-Bd driver genera",
    subtitle = paste0(
      "Selected by KSM rank among genera supported by robust ANCOM-BC2 and/or ",
      "ALDEx2 FDR < 0.05"
    ),
    x = NULL, y = "Anti-Bd relative abundance (% of all reads)",
    caption = paste0(
      "Every frog is shown (n = 47). Boxes show median and IQR; open points are ",
      "non-detections. Facets use independent y ranges."
    )
  ) +
  theme_classic(base_size = 7.8) +
  theme(strip.text = element_text(face = "italic", size = 7.4),
        axis.text.x = element_text(angle = 48, hjust = 1, size = 6.0),
        axis.title.y = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 7.5),
        plot.caption = element_text(size = 6.2, hjust = 0),
        legend.position = "top", panel.spacing = grid::unit(0.8, "lines"))

save_pub(p_violin,
         file.path(figure_dir, "figure_top9_anti_bd_driver_genera_violin"),
         183, 165)

cat("Three-method and top-nine anti-Bd driver figures completed.\n")
print(top9_table, row.names = FALSE)
