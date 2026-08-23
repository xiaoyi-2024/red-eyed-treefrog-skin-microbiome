#!/usr/bin/env Rscript

# Genus-level locality analysis for the 47 quality-controlled frog samples.
# Primary differential-abundance method: ANCOM-BC2.
# Sensitivity method: ALDEx2 Monte Carlo CLR with a four-group Kruskal-Wallis test.
# Robust Aitchison (rCLR) is used for the sample-resolved heat map.

suppressPackageStartupMessages({
  library(phyloseq)
  library(ANCOMBC)
  library(ALDEx2)
  library(vegan)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(svglite)
  library(ragg)
})

set.seed(20260815)
result_dir <- "results/genus_ancombc2_aldex2"
figure_dir <- file.path(result_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

whole_file <- "data/processed/ps_clean_47_rebuilt_from_01_filtered.rds"
anti_file <- "data/processed/blast_97id_90cov/ps_anti_bd_47_rebuilt_97id_90cov.rds"
annotation_file <- "data/processed/blast_97id_90cov/integrated_asv_taxonomy_function_rebuilt.csv"
stopifnot(file.exists(whole_file), file.exists(anti_file), file.exists(annotation_file))

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

genus_labels <- function(ps, taxa_ids) {
  tx <- as(tax_table(ps), "matrix")[taxa_ids, , drop = FALSE]
  genus <- as.character(tx[, "Genus"])
  family <- as.character(tx[, "Family"])
  family[is.na(family) | trimws(family) == ""] <- "Bacteria"
  missing <- is.na(genus) | trimws(genus) == ""
  genus[missing] <- paste0("Unclassified_", family[missing])
  genus
}

aggregate_genus <- function(ps) {
  x <- otu_taxa_rows(ps)
  g <- genus_labels(ps, rownames(x))
  rowsum(x, group = g, reorder = FALSE)
}

ps_whole <- readRDS(whole_file)
ps_anti <- readRDS(anti_file)
stopifnot(nsamples(ps_whole) == 47L,
          identical(sort(sample_names(ps_whole)), sort(sample_names(ps_anti))))

whole_genus <- aggregate_genus(ps_whole)
anti_genus <- aggregate_genus(ps_anti)
sample_ids <- colnames(whole_genus)
anti_genus <- anti_genus[, sample_ids, drop = FALSE]

metadata <- data.frame(sample_data(ps_whole), check.names = FALSE,
                       stringsAsFactors = FALSE)[sample_ids, , drop = FALSE]
metadata$Site <- factor(unname(site_key[metadata$Locality]), levels = site_levels)
stopifnot(!anyNA(metadata$Site))

# The same prespecified filter is used by both inferential methods and figures.
prevalence_cutoff <- 0.20
minimum_samples <- ceiling(ncol(whole_genus) * prevalence_cutoff)
whole_prevalence <- rowMeans(whole_genus > 0)
keep <- whole_prevalence >= prevalence_cutoff
filtered <- whole_genus[keep, , drop = FALSE]
stopifnot(nrow(filtered) > 0L, all(colSums(filtered) > 0))

filter_table <- data.frame(
  Genus = rownames(whole_genus),
  Prevalence = whole_prevalence,
  Prevalence_samples = rowSums(whole_genus > 0),
  Mean_relative_abundance = rowMeans(sweep(whole_genus, 2, colSums(whole_genus), "/")),
  Retained = keep,
  stringsAsFactors = FALSE
)
write.csv(filter_table, file.path(result_dir, "genus_filter_statistics.csv"), row.names = FALSE)
write.csv(data.frame(Genus = rownames(filtered), filtered, check.names = FALSE),
          file.path(result_dir, "filtered_genus_count_matrix.csv"), row.names = FALSE)

# Matrix + metadata input avoids requiring optional import helpers and makes the
# exact feature table supplied to ANCOM-BC2 explicit.
ancom_checkpoint <- file.path(result_dir, "ancombc2_complete_result.rds")
if (file.exists(ancom_checkpoint)) {
  ancom <- readRDS(ancom_checkpoint)
} else {
  ancom <- ANCOMBC::ancombc2(
    data = filtered, meta_data = metadata, taxa_are_rows = TRUE,
    fix_formula = "Site", group = "Site", prv_cut = 0, lib_cut = 0,
    p_adj_method = "BH", struc_zero = TRUE, neg_lb = TRUE,
    global = TRUE, pairwise = TRUE, pseudo_sens = TRUE,
    alpha = 0.05, n_cl = 1, verbose = FALSE
  )
  saveRDS(ancom, ancom_checkpoint)
}
write.csv(ancom$res_global, file.path(result_dir, "ancombc2_global_test.csv"), row.names = FALSE)
write.csv(ancom$res_pair, file.path(result_dir, "ancombc2_pairwise_tests.csv"), row.names = FALSE)
write.csv(ancom$res, file.path(result_dir, "ancombc2_model_coefficients.csv"), row.names = FALSE)
write.csv(ancom$zero_ind, file.path(result_dir, "ancombc2_structural_zeros.csv"), row.names = FALSE)
write.csv(ancom$ss_tab, file.path(result_dir, "ancombc2_pseudocount_sensitivity.csv"), row.names = FALSE)

# ALDEx2 sensitivity analysis: Monte Carlo CLR and a four-group KW test.
aldex_file <- file.path(result_dir, "aldex2_global_kw_sensitivity.csv")
if (file.exists(aldex_file)) {
  aldex_kw <- read.csv(aldex_file, stringsAsFactors = FALSE, check.names = FALSE)
} else {
  aldex_clr <- ALDEx2::aldex.clr(filtered, as.character(metadata$Site), mc.samples = 128,
                                denom = "all", verbose = FALSE)
  aldex_kw <- ALDEx2::aldex.kw(aldex_clr, verbose = FALSE)
  aldex_kw$taxon <- rownames(aldex_kw)
  rownames(aldex_kw) <- NULL
  write.csv(aldex_kw, aldex_file, row.names = FALSE)
}

# Join the two global tests without forcing agreement into a binary conclusion.
global <- merge(ancom$res_global, aldex_kw, by = "taxon", all.x = TRUE, sort = FALSE)
global$ANCOMBC2_robust_FDR05 <- global$diff_robust_abn %in% TRUE
global$ALDEx2_FDR05 <- !is.na(global$kw.eBH) & global$kw.eBH < 0.05
global$Method_agreement <- ifelse(global$ANCOMBC2_robust_FDR05 & global$ALDEx2_FDR05,
                                  "Both", ifelse(global$ANCOMBC2_robust_FDR05,
                                                 "ANCOM-BC2 only",
                                                 ifelse(global$ALDEx2_FDR05,
                                                        "ALDEx2 only", "Neither")))
global <- global[order(global$q_val, global$kw.eBH), ]
write.csv(global, file.path(result_dir, "ancombc2_aldex2_global_comparison.csv"), row.names = FALSE)

# Separate anti-Bd genus analysis. This prevents locality effects from
# non-inhibitory ASVs in the same named genus from driving the anti-Bd ranking.
anti_prevalence <- rowMeans(anti_genus > 0)
anti_keep <- anti_prevalence >= prevalence_cutoff
anti_filtered <- anti_genus[anti_keep, , drop = FALSE]
stopifnot(nrow(anti_filtered) > 1L, all(colSums(anti_filtered) > 0))
write.csv(data.frame(Genus = rownames(anti_filtered), anti_filtered, check.names = FALSE),
          file.path(result_dir, "filtered_anti_bd_genus_count_matrix.csv"), row.names = FALSE)

anti_ancom_checkpoint <- file.path(result_dir, "anti_bd_ancombc2_complete_result.rds")
if (file.exists(anti_ancom_checkpoint)) {
  anti_ancom <- readRDS(anti_ancom_checkpoint)
} else {
  anti_ancom <- ANCOMBC::ancombc2(
    data = anti_filtered, meta_data = metadata, taxa_are_rows = TRUE,
    fix_formula = "Site", group = "Site", prv_cut = 0, lib_cut = 0,
    p_adj_method = "BH", struc_zero = TRUE, neg_lb = TRUE,
    global = TRUE, pairwise = TRUE, pseudo_sens = TRUE,
    alpha = 0.05, n_cl = 1, verbose = FALSE
  )
  saveRDS(anti_ancom, anti_ancom_checkpoint)
}
write.csv(anti_ancom$res_global,
          file.path(result_dir, "anti_bd_ancombc2_global_test.csv"), row.names = FALSE)
write.csv(anti_ancom$res_pair,
          file.path(result_dir, "anti_bd_ancombc2_pairwise_tests.csv"), row.names = FALSE)
write.csv(anti_ancom$ss_tab,
          file.path(result_dir, "anti_bd_ancombc2_pseudocount_sensitivity.csv"),
          row.names = FALSE)

anti_aldex_file <- file.path(result_dir, "anti_bd_aldex2_global_kw_sensitivity.csv")
if (file.exists(anti_aldex_file)) {
  anti_aldex_kw <- read.csv(anti_aldex_file, stringsAsFactors = FALSE,
                            check.names = FALSE)
} else {
  anti_aldex_clr <- ALDEx2::aldex.clr(
    anti_filtered, as.character(metadata$Site), mc.samples = 128,
    denom = "all", verbose = FALSE
  )
  anti_aldex_kw <- ALDEx2::aldex.kw(anti_aldex_clr, verbose = FALSE)
  anti_aldex_kw$taxon <- rownames(anti_aldex_kw)
  rownames(anti_aldex_kw) <- NULL
  write.csv(anti_aldex_kw, anti_aldex_file, row.names = FALSE)
}
anti_global <- merge(anti_ancom$res_global, anti_aldex_kw,
                     by = "taxon", all.x = TRUE, sort = FALSE)
anti_global$ANCOMBC2_robust_FDR05 <- anti_global$diff_robust_abn %in% TRUE
anti_global$ALDEx2_FDR05 <- !is.na(anti_global$kw.eBH) & anti_global$kw.eBH < 0.05
anti_global$Method_agreement <- ifelse(
  anti_global$ANCOMBC2_robust_FDR05 & anti_global$ALDEx2_FDR05, "Both",
  ifelse(anti_global$ANCOMBC2_robust_FDR05, "ANCOM-BC2 only",
         ifelse(anti_global$ALDEx2_FDR05, "ALDEx2 only", "Neither"))
)
anti_global <- anti_global[order(anti_global$q_val, anti_global$kw.eBH), ]
write.csv(anti_global,
          file.path(result_dir, "anti_bd_ancombc2_aldex2_global_comparison.csv"),
          row.names = FALSE)

# Robust Aitchison heat map. Genera are selected by the global ANCOM-BC2 test;
# all 47 samples are retained. Samples are grouped by locality and clustered
# only within locality so the annotation remains easy to audit.
rclr <- vegan::decostand(t(filtered), method = "rclr")
stopifnot(!anyNA(rclr), all(is.finite(rclr)))
top20 <- head(global$taxon, min(20L, nrow(global)))
sample_order <- unlist(lapply(site_levels, function(s) {
  ids <- sample_ids[metadata$Site == s]
  if (length(ids) < 3L) return(ids)
  ids[hclust(dist(rclr[ids, top20, drop = FALSE]), method = "average")$order]
}), use.names = FALSE)
heat <- t(rclr[sample_order, top20, drop = FALSE])
heat_z <- t(scale(t(heat)))
heat_z[!is.finite(heat_z)] <- 0
heat_long <- as.data.frame(as.table(heat_z), stringsAsFactors = FALSE)
names(heat_long) <- c("Genus", "Sample_ID", "Row_scaled_rCLR")
heat_long$Genus <- factor(heat_long$Genus, levels = rev(top20))
heat_long$Sample_ID <- factor(heat_long$Sample_ID, levels = sample_order)
annotation <- data.frame(
  Sample_ID = factor(sample_order, levels = sample_order),
  Site = factor(metadata[sample_order, "Site"], levels = site_levels),
  Annotation = "Site"
)
write.csv(transform(heat_long, Genus = as.character(Genus),
                    Sample_ID = as.character(Sample_ID)),
          file.path(result_dir, "top20_rclr_heatmap_data.csv"), row.names = FALSE)
write.csv(transform(annotation, Sample_ID = as.character(Sample_ID)),
          file.path(result_dir, "top20_rclr_heatmap_annotation.csv"), row.names = FALSE)

p_ann <- ggplot(annotation, aes(Sample_ID, Annotation, fill = Site)) +
  geom_tile() + scale_fill_manual(values = site_colors) +
  labs(x = NULL, y = NULL, fill = "Locality") +
  theme_void(base_size = 7) +
  theme(legend.position = "top", legend.title = element_text(face = "bold"))
p_heat <- ggplot(heat_long, aes(Sample_ID, Genus, fill = Row_scaled_rCLR)) +
  geom_tile(colour = "white", linewidth = 0.16) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 0, limits = c(-2.5, 2.5), oob = scales::squish) +
  scale_y_discrete(drop = FALSE, expand = expansion(add = 0.5)) +
  labs(x = "Individual samples grouped by locality", y = NULL,
       fill = "Row-scaled\nrCLR") +
  theme_minimal(base_size = 7, base_family = "sans") +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.2),
        axis.text.y = element_text(face = "italic", colour = "black", size = 6.3),
        axis.title.x = element_text(face = "bold"))
figure_heat <- p_ann / p_heat + plot_layout(heights = c(0.55, 8.5), guides = "collect") +
  plot_annotation(
    title = "Locality-associated genus profiles across individual frogs",
    subtitle = paste0("Top 20 genera ranked by the ANCOM-BC2 global test; ",
                      "rCLR values are standardized within genus."),
    caption = paste0("All 47 frogs are shown. Features were retained at prevalence >=20% (>=",
                     minimum_samples, " samples)."),
    theme = theme(plot.title = element_text(size = 10, face = "bold"),
                  plot.subtitle = element_text(size = 7),
                  plot.caption = element_text(size = 6, hjust = 0))
  ) & theme(legend.position = "top")

# Top anti-Bd genera: all reads come from high-confidence Woodhams inhibitory
# matches. Ranking is transparent and hierarchical, not an arbitrary weighted score.
ann <- read.csv(annotation_file, stringsAsFactors = FALSE, check.names = FALSE)
ann_inhib <- ann[ann$Anti_Bd %in% TRUE | tolower(ann$Functional_category) == "inhibitory", ]
ann_genus <- as.character(ann_inhib$Genus)
ann_family <- as.character(ann_inhib$Family)
ann_family[is.na(ann_family) | trimws(ann_family) == ""] <- "Bacteria"
miss <- is.na(ann_genus) | trimws(ann_genus) == ""
ann_genus[miss] <- paste0("Unclassified_", ann_family[miss])
ann_inhib$Genus_label <- ann_genus

total_reads <- colSums(whole_genus)
anti_rel_whole <- sweep(anti_genus, 2, total_reads, "/")
anti_stats <- data.frame(
  Genus = rownames(anti_genus),
  Mean_relative_abundance_percent = 100 * rowMeans(anti_rel_whole),
  Prevalence = rowMeans(anti_genus > 0),
  Prevalence_samples = rowSums(anti_genus > 0),
  stringsAsFactors = FALSE
)
blast_summary <- do.call(rbind, lapply(split(ann_inhib, ann_inhib$Genus_label), function(z) {
  data.frame(Genus = z$Genus_label[1], Inhibitory_ASVs = nrow(z),
             Max_identity = max(z$Percent_identity, na.rm = TRUE),
             Max_query_coverage = max(z$Query_coverage, na.rm = TRUE))
}))
anti_rank <- merge(anti_stats, blast_summary, by = "Genus", all.x = TRUE, sort = FALSE)
whole_evidence <- global[, c("taxon", "q_val", "diff_robust_abn", "kw.eBH")]
names(whole_evidence) <- c("Genus", "Whole_ANCOMBC2_global_q",
                           "Whole_ANCOMBC2_robust_global",
                           "Whole_ALDEx2_global_q")
anti_evidence <- anti_global[, c("taxon", "q_val", "diff_robust_abn", "kw.eBH")]
names(anti_evidence) <- c("Genus", "Anti_only_ANCOMBC2_global_q",
                          "Anti_only_ANCOMBC2_robust_global",
                          "Anti_only_ALDEx2_global_q")
anti_rank <- merge(anti_rank, whole_evidence, by = "Genus", all.x = TRUE, sort = FALSE)
anti_rank <- merge(anti_rank, anti_evidence, by = "Genus", all.x = TRUE, sort = FALSE)
anti_rank$Woodhams_function <- "inhibitory"
anti_rank$Passes_20pct_anti_prevalence <- anti_rank$Prevalence >= prevalence_cutoff
anti_rank$Site_association_supported <-
  (anti_rank$Whole_ANCOMBC2_robust_global %in% TRUE) &
  (!is.na(anti_rank$Anti_only_ALDEx2_global_q) &
     anti_rank$Anti_only_ALDEx2_global_q < 0.05)
# "Top" means abundant and recurrent within the inhibitory reads. Locality
# evidence is reported alongside the ranking rather than allowed to promote a
# rare inhibitory genus merely because non-inhibitory ASVs of that genus vary.
anti_rank <- anti_rank[order(!anti_rank$Passes_20pct_anti_prevalence,
                             -anti_rank$Mean_relative_abundance_percent,
                             -anti_rank$Prevalence,
                             ifelse(is.na(anti_rank$Whole_ANCOMBC2_global_q), Inf,
                                    anti_rank$Whole_ANCOMBC2_global_q)), ]
anti_rank$Priority_rank <- seq_len(nrow(anti_rank))
write.csv(anti_rank, file.path(result_dir, "top_anti_bd_genus_priority_table.csv"),
          row.names = FALSE)

top15_anti <- head(anti_rank$Genus, min(15L, nrow(anti_rank)))
violin <- do.call(rbind, lapply(top15_anti, function(g) {
  data.frame(Genus = g, Sample_ID = sample_ids,
             Site = metadata[sample_ids, "Site"],
             Relative_abundance_percent = 100 * anti_rel_whole[g, sample_ids],
             Detected = anti_genus[g, sample_ids] > 0,
             stringsAsFactors = FALSE)
}))
display_genus <- top15_anti
display_genus[display_genus == "Unclassified_Alcaligenaceae"] <- "Uncl. Alcaligenaceae"
display_genus[display_genus == "Unclassified_Bacillaceae"] <- "Uncl. Bacillaceae"
display_genus[display_genus == "Burkholderia-Caballeronia-Paraburkholderia"] <-
  "Burkholderia complex"
display_genus[display_genus ==
                "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"] <-
  "Rhizobium complex"
prev_lab <- setNames(sprintf("%s\nprev. %.0f%%", display_genus,
                             100 * anti_rank$Prevalence[match(top15_anti, anti_rank$Genus)]),
                     top15_anti)
violin$Genus <- factor(violin$Genus, levels = top15_anti)
violin$Site <- factor(violin$Site, levels = site_levels)
write.csv(transform(violin, Genus = as.character(Genus), Site = as.character(Site)),
          file.path(result_dir, "top15_anti_bd_distribution_data.csv"), row.names = FALSE)

p_violin <- ggplot(violin, aes(Site, Relative_abundance_percent, fill = Site)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.45,
              colour = "grey35", linewidth = 0.25) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", linewidth = 0.3) +
  geom_jitter(aes(shape = Detected), width = 0.09, height = 0,
              size = 0.7, alpha = 0.72, colour = "#202020") +
  facet_wrap(~Genus, scales = "free_y", ncol = 5,
             labeller = as_labeller(prev_lab)) +
  scale_fill_manual(values = site_colors, guide = "none") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), name = "Detected") +
  labs(x = NULL, y = "Relative abundance (% of all reads)",
       caption = paste0("Boxes show medians and interquartile ranges; all individual frogs are shown. ",
                        "Facet labels report overall prevalence.")) +
  theme_classic(base_size = 7, base_family = "sans") +
  theme(strip.text = element_text(face = "italic", size = 6.2),
        axis.text.x = element_text(angle = 55, hjust = 1, size = 5.2),
        axis.title.y = element_text(face = "bold"),
        legend.position = "top", panel.spacing = grid::unit(0.75, "lines"),
        plot.caption = element_text(size = 5.8, hjust = 0))

# Exploratory anti-Bd-only ANCOM-BC2 heat map. Only 13 genera were estimable by
# ANCOM-BC2 after the prespecified prevalence filter, so all 13 are displayed
# rather than padding the figure to an arbitrary Top 15. Inclusion is ordered by
# the global test; rows are subsequently clustered by their anti-Bd rCLR profile.
anti_rclr <- vegan::decostand(t(anti_filtered), method = "rclr")
top_anti_ancom <- anti_global$taxon[order(anti_global$q_val,
                                          -abs(anti_global$W), na.last = TRUE)]
stopifnot(all(top_anti_ancom %in% colnames(anti_rclr)),
          !anyNA(anti_rclr), all(is.finite(anti_rclr)))
anti_sample_order <- unlist(lapply(site_levels, function(s) {
  ids <- sample_ids[metadata$Site == s]
  if (length(ids) < 3L) return(ids)
  ids[hclust(dist(anti_rclr[ids, top_anti_ancom, drop = FALSE]),
             method = "average")$order]
}), use.names = FALSE)
anti_heat <- t(anti_rclr[anti_sample_order, top_anti_ancom, drop = FALSE])
anti_genus_order <- rownames(anti_heat)[
  hclust(dist(anti_heat), method = "average")$order
]
anti_heat <- anti_heat[anti_genus_order, , drop = FALSE]
anti_heat_z <- t(scale(t(anti_heat)))
anti_heat_z[!is.finite(anti_heat_z)] <- 0
anti_heat_long <- as.data.frame(as.table(anti_heat_z), stringsAsFactors = FALSE)
names(anti_heat_long) <- c("Genus", "Sample_ID", "Row_scaled_rCLR")
anti_display <- gsub("^Unclassified_", "Uncl. ", anti_genus_order)
anti_display[anti_display ==
  "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"] <- "Rhizobium complex"
anti_display[anti_display ==
  "Burkholderia-Caballeronia-Paraburkholderia"] <- "Burkholderia complex"
anti_display_lookup <- setNames(anti_display, anti_genus_order)
anti_heat_long$Genus_full <- anti_heat_long$Genus
anti_heat_long$Genus <- factor(anti_display_lookup[anti_heat_long$Genus_full],
                               levels = rev(anti_display))
anti_heat_long$Sample_ID <- factor(anti_heat_long$Sample_ID,
                                   levels = anti_sample_order)
anti_heat_annotation <- data.frame(
  Sample_ID = factor(anti_sample_order, levels = anti_sample_order),
  Site = factor(metadata[anti_sample_order, "Site"], levels = site_levels),
  Annotation = "Site"
)
write.csv(transform(anti_heat_long, Genus = as.character(Genus),
                    Sample_ID = as.character(Sample_ID)),
          file.path(result_dir, "anti_bd_ancombc2_top13_rclr_heatmap_data.csv"),
          row.names = FALSE)
write.csv(transform(anti_heat_annotation, Sample_ID = as.character(Sample_ID)),
          file.path(result_dir, "anti_bd_ancombc2_top13_sample_annotation.csv"),
          row.names = FALSE)

anti_test_annotation <- data.frame(
  Genus_full = anti_genus_order,
  Genus = factor(anti_display, levels = rev(anti_display)),
  Robust_FDR05 = anti_global$diff_robust_abn[
    match(anti_genus_order, anti_global$taxon)
  ] %in% TRUE,
  ANCOMBC2_q = anti_global$q_val[match(anti_genus_order, anti_global$taxon)],
  stringsAsFactors = FALSE
)
write.csv(transform(anti_test_annotation, Genus = as.character(Genus)),
          file.path(result_dir, "anti_bd_ancombc2_top13_test_annotation.csv"),
          row.names = FALSE)

p_anti_ann <- ggplot(anti_heat_annotation,
                     aes(Sample_ID, Annotation, fill = Site)) +
  geom_tile() + scale_fill_manual(values = site_colors) +
  labs(x = NULL, y = NULL, fill = "Locality") +
  theme_void(base_size = 7) +
  theme(legend.position = "top", legend.title = element_text(face = "bold"))
p_anti_heat <- ggplot(anti_heat_long,
                      aes(Sample_ID, Genus, fill = Row_scaled_rCLR)) +
  geom_tile(colour = "white", linewidth = 0.16) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 0, limits = c(-2.5, 2.5), oob = scales::squish) +
  scale_y_discrete(drop = FALSE, expand = expansion(add = 0.5)) +
  labs(x = "Individual samples grouped by locality", y = NULL,
       fill = "Row-scaled\nrCLR") +
  theme_minimal(base_size = 7, base_family = "sans") +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.2),
        axis.text.y = element_text(face = "italic", colour = "black", size = 6.3),
        axis.title.x = element_text(face = "bold"),
        plot.margin = margin(0, 2, 0, 0))
p_anti_test <- ggplot(anti_test_annotation,
                      aes(1, Genus, colour = Robust_FDR05)) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#BDBDBD"),
                      guide = "none") +
  scale_y_discrete(drop = FALSE, expand = expansion(add = 0.5)) +
  scale_x_continuous(limits = c(0.7, 1.3), breaks = NULL) +
  labs(x = NULL, y = NULL, title = "Robust\nFDR < 0.05") +
  theme_minimal(base_size = 6, base_family = "sans") +
  theme(panel.grid = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 5.8, face = "bold", hjust = 0.5),
        plot.margin = margin(0, 0, 0, 1))
anti_top_row <- (p_anti_ann | plot_spacer()) + plot_layout(widths = c(24, 1))
anti_body <- (p_anti_heat | p_anti_test) + plot_layout(widths = c(24, 1))
figure_anti_heat <- (anti_top_row / anti_body) +
  plot_layout(heights = c(0.55, 8.5), guides = "collect") +
  plot_annotation(
    title = "Anti-Bd genus profiles from the ANCOM-BC2 locality analysis",
    subtitle = paste0("All 13 estimable genera are shown; rows cluster anti-Bd rCLR ",
                      "profiles and orange markers denote robust global FDR < 0.05."),
    caption = paste0("All 47 frogs are shown. This anti-Bd-only ANCOM-BC2 analysis is ",
                     "exploratory because fewer than 50 genera were available for sampling-fraction estimation."),
    theme = theme(plot.title = element_text(size = 10, face = "bold"),
                  plot.subtitle = element_text(size = 7),
                  plot.caption = element_text(size = 6, hjust = 0))
  ) & theme(legend.position = "top")

# Final manuscript heat map: Top 30 whole-community genera, with a separate
# row annotation showing whether each genus contains at least one high-
# confidence Woodhams inhibitory ASV. This avoids using a significance-like
# asterisk for functional annotation.
top30_whole <- head(global$taxon, min(30L, nrow(global)))
whole_sample_order <- unlist(lapply(site_levels, function(s) {
  ids <- sample_ids[metadata$Site == s]
  if (length(ids) < 3L) return(ids)
  ids[hclust(dist(rclr[ids, top30_whole, drop = FALSE]), method = "average")$order]
}), use.names = FALSE)
whole_heat <- t(rclr[whole_sample_order, top30_whole, drop = FALSE])
# Selection and ordering answer different questions: ANCOM-BC2 selects the 30
# locality-associated genera, then average-linkage clustering orders those rows
# by their rCLR patterns across all 47 frogs.
whole_genus_order <- rownames(whole_heat)[
  hclust(dist(whole_heat), method = "average")$order
]
whole_heat <- whole_heat[whole_genus_order, , drop = FALSE]
whole_heat_z <- t(scale(t(whole_heat)))
whole_heat_z[!is.finite(whole_heat_z)] <- 0

display_top30 <- gsub("^Unclassified_", "Uncl. ", whole_genus_order)
display_top30[nchar(display_top30) > 34] <- paste0(
  substr(display_top30[nchar(display_top30) > 34], 1, 31), "..."
)
display_lookup <- setNames(display_top30, whole_genus_order)
whole_heat_long <- as.data.frame(as.table(whole_heat_z), stringsAsFactors = FALSE)
names(whole_heat_long) <- c("Genus_full", "Sample_ID", "Row_scaled_rCLR")
whole_heat_long$Genus <- factor(display_lookup[whole_heat_long$Genus_full],
                                levels = rev(display_top30))
whole_heat_long$Sample_ID <- factor(whole_heat_long$Sample_ID,
                                    levels = whole_sample_order)
whole_heat_annotation <- data.frame(
  Sample_ID = factor(whole_sample_order, levels = whole_sample_order),
  Site = factor(metadata[whole_sample_order, "Site"], levels = site_levels),
  Annotation = "Site"
)
inhibitory_genera <- unique(ann_inhib$Genus_label)
whole_row_annotation <- data.frame(
  Genus_full = whole_genus_order,
  Genus = factor(display_top30, levels = rev(display_top30)),
  Woodhams_inhibitory_match = whole_genus_order %in% inhibitory_genera,
  stringsAsFactors = FALSE
)

write.csv(transform(whole_heat_long, Genus = as.character(Genus),
                    Sample_ID = as.character(Sample_ID)),
          file.path(result_dir, "top30_whole_community_rclr_heatmap_data.csv"),
          row.names = FALSE)
write.csv(transform(whole_heat_annotation, Sample_ID = as.character(Sample_ID)),
          file.path(result_dir, "top30_whole_community_heatmap_sample_annotation.csv"),
          row.names = FALSE)
write.csv(transform(whole_row_annotation, Genus = as.character(Genus)),
          file.path(result_dir, "top30_whole_community_antibd_row_annotation.csv"),
          row.names = FALSE)

p_whole_ann <- ggplot(whole_heat_annotation,
                      aes(Sample_ID, Annotation, fill = Site)) +
  geom_tile() + scale_fill_manual(values = site_colors) +
  labs(x = NULL, y = NULL, fill = "Locality") +
  theme_void(base_size = 7) +
  theme(legend.position = "top", legend.title = element_text(face = "bold"))
p_whole_heat <- ggplot(whole_heat_long,
                       aes(Sample_ID, Genus, fill = Row_scaled_rCLR)) +
  geom_tile(colour = "white", linewidth = 0.14) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 0, limits = c(-2.5, 2.5), oob = scales::squish) +
  scale_y_discrete(drop = FALSE, expand = c(0, 0)) +
  labs(x = "Individual samples grouped by locality", y = NULL,
       fill = "Row-scaled\nrCLR") +
  theme_minimal(base_size = 7, base_family = "sans") +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
        axis.text.y = element_text(face = "italic", colour = "black", size = 5.8),
        axis.title.x = element_text(face = "bold"),
        plot.margin = margin(0, 2, 0, 0))
p_antibd_marker <- ggplot(whole_row_annotation,
                          aes(1, Genus, colour = Woodhams_inhibitory_match)) +
  geom_point(size = 1.7) +
  scale_colour_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#D0D0D0"),
                      guide = "none") +
  scale_y_discrete(drop = FALSE, expand = c(0, 0)) +
  scale_x_continuous(limits = c(0.7, 1.3), breaks = NULL) +
  labs(x = NULL, y = NULL, title = "Woodhams\ninhibitory") +
  theme_minimal(base_size = 6, base_family = "sans") +
  theme(panel.grid = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), plot.title = element_text(size = 5.8,
        face = "bold", hjust = 0.5), plot.margin = margin(0, 0, 0, 1))

whole_top_row <- (p_whole_ann | plot_spacer()) + plot_layout(widths = c(24, 1))
whole_body <- (p_whole_heat | p_antibd_marker) + plot_layout(widths = c(24, 1))
figure_whole_heat <- (whole_top_row / whole_body) +
  plot_layout(heights = c(0.55, 12), guides = "collect") +
  plot_annotation(
    title = "Whole-community genus profiles across individual frogs",
    subtitle = paste0("Top 30 genera ranked by the ANCOM-BC2 global locality test; ",
                      "rows cluster rCLR profiles; orange markers denote Woodhams inhibitory matches."),
    caption = paste0("All 47 frogs are shown. Values are rCLR coordinates standardized ",
                     "within genus; grey markers indicate no inhibitory match."),
    theme = theme(plot.title = element_text(size = 10, face = "bold"),
                  plot.subtitle = element_text(size = 7),
                  plot.caption = element_text(size = 6, hjust = 0))
  ) & theme(legend.position = "top")

save_pub <- function(plot, stem, width_mm, height_mm) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  ggsave(paste0(stem, ".png"), plot, width = w, height = h, dpi = 600, bg = "white")
  svglite::svglite(paste0(stem, ".svg"), width = w, height = h); print(plot); dev.off()
  tryCatch({
    grDevices::pdf(paste0(stem, ".pdf"), width = w, height = h,
                   family = "Helvetica", useDingbats = FALSE)
    print(plot); dev.off()
  }, error = function(e) {
    if (grDevices::dev.cur() > 1L) try(grDevices::dev.off(), silent = TRUE)
    warning("PDF export failed after PNG/SVG succeeded: ", conditionMessage(e))
  })
  tryCatch({
    ragg::agg_tiff(paste0(stem, ".tiff"), width = w, height = h,
                   units = "in", res = 600, background = "white",
                   compression = "lzw")
    print(plot); dev.off()
  }, error = function(e) {
    if (grDevices::dev.cur() > 1L) try(grDevices::dev.off(), silent = TRUE)
    warning("TIFF export failed after PNG/PDF/SVG succeeded: ", conditionMessage(e))
  })
}
old_top20 <- c(
  file.path(result_dir, "top20_rclr_heatmap_data.csv"),
  file.path(result_dir, "top20_rclr_heatmap_annotation.csv"),
  file.path(figure_dir, paste0("figure_top20_genus_rclr_heatmap.",
                              c("png", "svg", "pdf", "tiff")))
)
unlink(old_top20[file.exists(old_top20)])
old_anti_heat <- c(
  file.path(result_dir, "top15_anti_bd_rclr_heatmap_data.csv"),
  file.path(result_dir, "top15_anti_bd_rclr_heatmap_annotation.csv"),
  file.path(figure_dir, paste0("figure_top15_anti_bd_genus_rclr_heatmap.",
                              c("png", "svg", "pdf", "tiff")))
)
unlink(old_anti_heat[file.exists(old_anti_heat)])
save_pub(figure_whole_heat,
         file.path(figure_dir, "figure_top30_whole_community_rclr_heatmap_with_antibd"),
         183, 155)
save_pub(figure_anti_heat,
         file.path(figure_dir, "figure_anti_bd_ancombc2_top13_rclr_heatmap"),
         183, 112)
save_pub(p_violin, file.path(figure_dir, "figure_top15_anti_bd_genus_distributions"), 183, 150)

summary <- data.frame(
  Metric = c("Frog samples", "Genera before filtering", "Genera after filtering",
             "Genera modeled by ANCOM-BC2",
             "Prevalence threshold", "ANCOM-BC2 robust global FDR < 0.05",
             "ALDEx2 global FDR < 0.05", "Significant by both methods"),
  Value = c(ncol(whole_genus), nrow(whole_genus), nrow(filtered), nrow(global),
            prevalence_cutoff,
            sum(global$ANCOMBC2_robust_FDR05, na.rm = TRUE),
            sum(global$ALDEx2_FDR05, na.rm = TRUE),
            sum(global$Method_agreement == "Both", na.rm = TRUE)),
  stringsAsFactors = FALSE
)
write.csv(summary, file.path(result_dir, "analysis_summary.csv"), row.names = FALSE)
write.csv(data.frame(Package = c("ANCOMBC", "ALDEx2", "vegan"),
                     Version = c(as.character(packageVersion("ANCOMBC")),
                                 as.character(packageVersion("ALDEx2")),
                                 as.character(packageVersion("vegan")))),
          file.path(result_dir, "software_versions.csv"), row.names = FALSE)
print(summary)
