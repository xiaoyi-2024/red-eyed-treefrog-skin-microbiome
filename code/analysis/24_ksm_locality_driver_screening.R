#!/usr/bin/env Rscript

# Kolmogorov-Smirnov Measure (KSM) screening of locality-associated genera.
#
# Purpose
#   KSM is used as a complementary, non-parametric feature-ranking method.
#   It is not treated as a differential-abundance test and no KSM p-values are
#   reported. The analysis is run separately for the whole community and for
#   Woodhams-inhibitory (anti-Bd) reads.
#
# Design
#   - 47 independent frog samples; four locality groups.
#   - genus-level features; prevalence >= 20% within the relevant matrix.
#   - whole-community abundance = genus reads / complete-library reads.
#   - anti-Bd abundance = inhibitory genus reads / complete-library reads.
#   - 1,000 locality-stratified bootstrap resamples quantify uncertainty and
#     stability of entry into the top 10 and top 20 KSM ranks.
#   - ANCOM-BC2 and ALDEx2 results are joined for triangulation only.
#
# Reference
#   Loftus SC et al. (2015), Dimension Reduction for Multinomial Models Via a
#   Kolmogorov-Smirnov Measure (KSM), Virginia Tech Technical Report 15-1.

suppressPackageStartupMessages({
  library(phyloseq)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
})

set.seed(20260816)

result_dir <- "results/ksm_locality_driver_screening"
figure_dir <- file.path(result_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

whole_file <- "data/processed/ps_clean_47_rebuilt_from_01_filtered.rds"
anti_file <- "data/processed/blast_97id_90cov/ps_anti_bd_47_rebuilt_97id_90cov.rds"
method_dir <- "results/genus_ancombc2_aldex2"
whole_method_file <- file.path(method_dir, "ancombc2_aldex2_global_comparison.csv")
anti_method_file <- file.path(method_dir, "anti_bd_ancombc2_aldex2_global_comparison.csv")
whole_aldex_file <- file.path(method_dir, "aldex2_global_kw_sensitivity.csv")
anti_aldex_file <- file.path(method_dir, "anti_bd_aldex2_global_kw_sensitivity.csv")
stopifnot(all(file.exists(c(whole_file, anti_file, whole_method_file,
                            anti_method_file, whole_aldex_file,
                            anti_aldex_file))))

site_levels <- c("Lost Iguana", "Veragua", "Altos de Campana", "Soberanía")
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

# Two-sample KS distance without significance testing. Ties/zeros are retained
# because differences in zero frequency are part of the observed distribution.
ks_distance <- function(x, y) {
  nx <- length(x)
  ny <- length(y)
  z <- c(x, y)
  ord <- order(z)
  z_sorted <- z[ord]
  fx <- cumsum(ord <= nx) / nx
  fy <- cumsum(ord > nx) / ny
  # Evaluate after the last member of every tied block. This is equivalent to
  # comparing the two empirical CDFs at every unique observed abundance.
  block_end <- !duplicated(z_sorted, fromLast = TRUE)
  max(abs(fx[block_end] - fy[block_end]))
}

# Equation 6 of Loftus et al. The six unique locality pairs are weighted by
# (n_k + n_k') / (N * (K - 1)); these weights sum to one.
ksm_one <- function(x, group) {
  lev <- levels(group)
  n <- length(x)
  k <- length(lev)
  cmb <- combn(lev, 2, simplify = FALSE)
  sum(vapply(cmb, function(pair) {
    a <- x[group == pair[1]]
    b <- x[group == pair[2]]
    ((length(a) + length(b)) / (n * (k - 1))) * ks_distance(a, b)
  }, numeric(1)))
}

ksm_vector <- function(mat, group) {
  # mat: features x samples
  vapply(seq_len(nrow(mat)), function(i) ksm_one(mat[i, ], group), numeric(1))
}

pairwise_ks <- function(mat, group, community) {
  pairs <- combn(levels(group), 2, simplify = FALSE)
  out <- do.call(rbind, lapply(pairs, function(pair) {
    data.frame(
      Community = community,
      Genus = rownames(mat),
      Site_1 = pair[1], Site_2 = pair[2],
      KS_distance = vapply(seq_len(nrow(mat)), function(i) {
        ks_distance(mat[i, group == pair[1]], mat[i, group == pair[2]])
      }, numeric(1)),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

bootstrap_ksm <- function(mat, group, b = 1000L) {
  site_ids <- split(seq_along(group), group, drop = TRUE)
  features <- rownames(mat)
  one_bootstrap <- function(iter) {
    idx <- unlist(lapply(site_ids, function(z) sample(z, length(z), replace = TRUE)),
                  use.names = FALSE)
    g_boot <- factor(as.character(group[idx]), levels = levels(group))
    ksm_vector(mat[, idx, drop = FALSE], g_boot)
  }
  cores <- max(1L, min(4L, parallel::detectCores(logical = FALSE)))
  values <- simplify2array(parallel::mclapply(
    seq_len(b), one_bootstrap, mc.cores = cores, mc.preschedule = TRUE,
    mc.set.seed = TRUE
  ))
  rownames(values) <- features
  ranks <- apply(values, 2, function(v) rank(-v, ties.method = "first"))
  top10 <- rowSums(ranks <= min(10L, nrow(values)))
  top20 <- rowSums(ranks <= min(20L, nrow(values)))
  data.frame(
    Genus = features,
    Bootstrap_median = apply(values, 1, median),
    Bootstrap_CI_low = apply(values, 1, quantile, probs = 0.025, names = FALSE),
    Bootstrap_CI_high = apply(values, 1, quantile, probs = 0.975, names = FALSE),
    Top10_frequency = top10 / b,
    Top20_frequency = top20 / b,
    Bootstrap_iterations = b,
    stringsAsFactors = FALSE
  )
}

method_support <- function(path, aldex_path) {
  z <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  ancom <- data.frame(
    Genus = z$taxon,
    ANCOMBC2_robust_FDR05 = z$ANCOMBC2_robust_FDR05 %in% TRUE,
    ANCOMBC2_q = z$q_val,
    stringsAsFactors = FALSE
  )
  a <- read.csv(aldex_path, stringsAsFactors = FALSE, check.names = FALSE)
  aldex <- data.frame(Genus = a$taxon, ALDEx2_q = a$kw.eBH,
                      ALDEx2_FDR05 = !is.na(a$kw.eBH) & a$kw.eBH < 0.05,
                      stringsAsFactors = FALSE)
  merge(ancom, aldex, by = "Genus", all = TRUE, sort = FALSE)
}

refresh_support <- function(ans, support) {
  support_columns <- c("ANCOMBC2_robust_FDR05", "ALDEx2_FDR05",
                       "ANCOMBC2_q", "ALDEx2_q", "Method_agreement",
                       "Cross_method_support")
  ans <- ans[, setdiff(names(ans), support_columns), drop = FALSE]
  ans <- merge(ans, support, by = "Genus", all.x = TRUE, sort = FALSE)
  ans$ANCOMBC2_robust_FDR05[is.na(ans$ANCOMBC2_robust_FDR05)] <- FALSE
  ans$ALDEx2_FDR05[is.na(ans$ALDEx2_FDR05)] <- FALSE
  ans$Method_agreement <- ifelse(
    ans$ANCOMBC2_robust_FDR05 & ans$ALDEx2_FDR05, "Both",
    ifelse(ans$ANCOMBC2_robust_FDR05, "ANCOM-BC2 only",
           ifelse(ans$ALDEx2_FDR05, "ALDEx2 only", "Neither"))
  )
  ans$Cross_method_support <- ifelse(
    ans$ANCOMBC2_robust_FDR05 & ans$ALDEx2_FDR05, "ANCOM-BC2 + ALDEx2",
    ifelse(ans$ANCOMBC2_robust_FDR05, "ANCOM-BC2 only",
           ifelse(ans$ALDEx2_FDR05, "ALDEx2 only", "KSM only"))
  )
  ans[order(ans$KSM_rank), ]
}

run_community <- function(mat, group, community, support, anti_genera,
                          prevalence_cutoff = 0.20, b = 1000L) {
  prevalence <- rowMeans(mat > 0)
  keep <- prevalence >= prevalence_cutoff
  filtered <- mat[keep, , drop = FALSE]
  observed <- ksm_vector(filtered, group)
  boot <- bootstrap_ksm(filtered, group, b = b)
  ans <- merge(
    data.frame(Genus = rownames(filtered), Community = community,
               KSM = observed, Prevalence = prevalence[rownames(filtered)],
               Detected_samples = rowSums(filtered > 0),
               Mean_relative_abundance_percent = 100 * rowMeans(filtered),
               Is_anti_Bd = rownames(filtered) %in% anti_genera,
               stringsAsFactors = FALSE),
    boot, by = "Genus", all.x = TRUE, sort = FALSE
  )
  ans <- merge(ans, support, by = "Genus", all.x = TRUE, sort = FALSE)
  ans$ANCOMBC2_robust_FDR05[is.na(ans$ANCOMBC2_robust_FDR05)] <- FALSE
  ans$ALDEx2_FDR05[is.na(ans$ALDEx2_FDR05)] <- FALSE
  ans$Method_agreement[is.na(ans$Method_agreement)] <- "Not retained/tested"
  ans$Cross_method_support <- ifelse(
    ans$ANCOMBC2_robust_FDR05 & ans$ALDEx2_FDR05, "ANCOM-BC2 + ALDEx2",
    ifelse(ans$ANCOMBC2_robust_FDR05, "ANCOM-BC2 only",
           ifelse(ans$ALDEx2_FDR05, "ALDEx2 only", "KSM only"))
  )
  ans <- ans[order(-ans$KSM, -ans$Top20_frequency, ans$Genus), ]
  ans$KSM_rank <- seq_len(nrow(ans))
  list(table = ans, filtered = filtered,
       pairwise = pairwise_ks(filtered, group, community))
}

ps_whole <- readRDS(whole_file)
ps_anti <- readRDS(anti_file)
stopifnot(nsamples(ps_whole) == 47L,
          identical(sort(sample_names(ps_whole)), sort(sample_names(ps_anti))))

whole_counts <- aggregate_genus(ps_whole)
anti_counts <- aggregate_genus(ps_anti)
sample_ids <- colnames(whole_counts)
anti_counts <- anti_counts[, sample_ids, drop = FALSE]
metadata <- data.frame(sample_data(ps_whole), check.names = FALSE,
                       stringsAsFactors = FALSE)[sample_ids, , drop = FALSE]
metadata$Site <- factor(unname(site_key[metadata$Locality]), levels = site_levels)
stopifnot(!anyNA(metadata$Site), all(colSums(whole_counts) > 0))

library_size <- colSums(whole_counts)
whole_relative <- sweep(whole_counts, 2, library_size, "/")
# The denominator remains all bacterial reads so KSM captures the ecological
# abundance of each inhibitory genus, not merely its share of the anti-Bd subset.
anti_relative <- sweep(anti_counts, 2, library_size, "/")
anti_genera <- rownames(anti_counts)[rowSums(anti_counts) > 0]

whole_rank_file <- file.path(result_dir, "whole_genus_ksm_ranked.csv")
anti_rank_file <- file.path(result_dir, "anti_bd_genus_ksm_ranked.csv")
pairwise_file <- file.path(result_dir, "genus_pairwise_site_ks_distances.csv")
if (all(file.exists(c(whole_rank_file, anti_rank_file, pairwise_file)))) {
  pairwise_cached <- read.csv(pairwise_file, stringsAsFactors = FALSE,
                              check.names = FALSE)
  whole <- list(
    table = refresh_support(
      read.csv(whole_rank_file, stringsAsFactors = FALSE, check.names = FALSE),
      method_support(whole_method_file, whole_aldex_file)),
    filtered = whole_relative[rowMeans(whole_relative > 0) >= 0.20, , drop = FALSE],
    pairwise = pairwise_cached[pairwise_cached$Community == "Whole community", ]
  )
  anti <- list(
    table = refresh_support(
      read.csv(anti_rank_file, stringsAsFactors = FALSE, check.names = FALSE),
      method_support(anti_method_file, anti_aldex_file)),
    filtered = anti_relative[rowMeans(anti_relative > 0) >= 0.20, , drop = FALSE],
    pairwise = pairwise_cached[pairwise_cached$Community == "Anti-Bd community", ]
  )
} else {
  whole <- run_community(
    whole_relative, metadata$Site, "Whole community",
    method_support(whole_method_file, whole_aldex_file), anti_genera, b = 1000L
  )
  anti <- run_community(
    anti_relative, metadata$Site, "Anti-Bd community",
    method_support(anti_method_file, anti_aldex_file), anti_genera, b = 1000L
  )
}

all_results <- rbind(whole$table, anti$table)
write.csv(whole$table, file.path(result_dir, "whole_genus_ksm_ranked.csv"),
          row.names = FALSE)
write.csv(anti$table, file.path(result_dir, "anti_bd_genus_ksm_ranked.csv"),
          row.names = FALSE)
write.csv(all_results, file.path(result_dir, "whole_and_anti_bd_ksm_ranked.csv"),
          row.names = FALSE)
write.csv(rbind(whole$pairwise, anti$pairwise),
          file.path(result_dir, "genus_pairwise_site_ks_distances.csv"),
          row.names = FALSE)

top_n <- 20L
plot_data <- rbind(head(whole$table, top_n), head(anti$table, top_n))
plot_data$Community <- factor(plot_data$Community,
                              levels = c("Whole community", "Anti-Bd community"))
plot_data$Display_genus <- gsub("^Unclassified_", "Uncl. ", plot_data$Genus)
plot_data$Display_genus[plot_data$Display_genus ==
  "Burkholderia-Caballeronia-Paraburkholderia"] <- "Burkholderia complex"
plot_data$Plot_genus <- paste(plot_data$Community, plot_data$Display_genus,
                              sep = "___")
plot_data$Plot_genus <- factor(
  plot_data$Plot_genus,
  levels = plot_data$Plot_genus[order(plot_data$Community, plot_data$KSM)]
)
plot_data$Cross_method_support <- factor(
  plot_data$Cross_method_support,
  levels = c("ANCOM-BC2 + ALDEx2", "ANCOM-BC2 only", "ALDEx2 only", "KSM only")
)
write.csv(plot_data, file.path(result_dir, "figure_top20_ksm_source_data.csv"),
          row.names = FALSE)

# Compact manuscript-facing table for the genus heat-map / locality-driver
# section. KSM is a distributional ranking score; q-values retain their exact
# inferential meaning and are not converted into significance stars.
driver_table <- plot_data[, c(
  "Community", "KSM_rank", "Genus", "KSM", "Bootstrap_CI_low",
  "Bootstrap_CI_high", "Top10_frequency", "Prevalence",
  "Mean_relative_abundance_percent", "Is_anti_Bd",
  "ANCOMBC2_q", "ANCOMBC2_robust_FDR05", "ALDEx2_q", "ALDEx2_FDR05",
  "Cross_method_support"
)]
names(driver_table) <- c(
  "Community", "KSM_rank", "Genus", "KSM", "KSM_bootstrap_CI_low",
  "KSM_bootstrap_CI_high", "Bootstrap_top10_frequency", "Prevalence",
  "Mean_relative_abundance_percent", "Woodhams_inhibitory",
  "ANCOMBC2_BH_q", "ANCOMBC2_robust_FDR05", "ALDEx2_BH_q",
  "ALDEx2_FDR05", "Evidence_class"
)
driver_table <- driver_table[order(driver_table$Community,
                                   driver_table$KSM_rank), ]
write.csv(driver_table,
          file.path(result_dir, "table_top20_locality_driver_genera.csv"),
          row.names = FALSE)
write.csv(driver_table[driver_table$Community == "Whole community", ],
          file.path(result_dir, "table_top20_whole_locality_driver_genera.csv"),
          row.names = FALSE)
write.csv(driver_table[driver_table$Community == "Anti-Bd community", ],
          file.path(result_dir, "table_top20_anti_bd_locality_driver_genera.csv"),
          row.names = FALSE)

support_colors <- c("ANCOM-BC2 + ALDEx2" = "#0072B2",
                    "ANCOM-BC2 only" = "#D55E00",
                    "ALDEx2 only" = "#009E73",
                    "KSM only" = "#8A8A8A")

p <- ggplot(plot_data,
            aes(x = KSM, y = Plot_genus)) +
  geom_errorbar(aes(xmin = Bootstrap_CI_low, xmax = Bootstrap_CI_high),
                orientation = "y", width = 0, colour = "grey68", linewidth = 0.5) +
  geom_point(aes(fill = Cross_method_support, size = Top10_frequency,
                 shape = Is_anti_Bd), colour = "black", stroke = 0.45) +
  facet_wrap(~Community, scales = "free_y", ncol = 2) +
  scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
  scale_fill_manual(values = support_colors, drop = FALSE) +
  scale_shape_manual(values = c(`TRUE` = 24, `FALSE` = 21),
                     labels = c(`TRUE` = "Woodhams inhibitory",
                                `FALSE` = "Not annotated inhibitory")) +
  scale_size_continuous(range = c(2.0, 5.2), limits = c(0, 1),
                        breaks = c(0.25, 0.50, 0.75, 1.00),
                        labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "KSM screening of locality-associated bacterial genera",
    subtitle = paste0("Top 20 per community after a 20% prevalence filter.\n",
                      "Points are observed KSM; bars are stratified-bootstrap 95% intervals."),
    x = "Kolmogorov-Smirnov Measure (KSM)",
    y = NULL, fill = "Cross-method support",
    shape = "Anti-Bd annotation", size = "Bootstrap top-10\nfrequency",
    caption = paste0(
      "n = 47 independent frogs; 1,000 bootstrap resamples within locality. ",
      "KSM is a feature-ranking measure, not a significance test."
    )
  ) +
  theme_classic(base_size = 8.2, base_family = "sans") +
  theme(
    strip.background = element_rect(fill = "#F1F1F1", colour = NA),
    strip.text = element_text(face = "bold", size = 9),
    axis.text.y = element_text(face = "italic", colour = "black", size = 6.8),
    axis.title.x = element_text(face = "bold"),
    legend.position = "bottom", legend.box = "vertical",
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 8.2),
    plot.caption = element_text(size = 6.5, hjust = 0),
    panel.spacing.x = grid::unit(1.5, "lines")
  ) +
  guides(fill = guide_legend(order = 1, nrow = 1,
                             override.aes = list(shape = 21, size = 3.2)),
         shape = guide_legend(order = 2),
         size = guide_legend(order = 3, nrow = 1))

save_pub <- function(plot, stem, width_mm = 183, height_mm = 180) {
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

save_pub(p, file.path(figure_dir, "figure_ksm_locality_driver_screening"))

summary_table <- data.frame(
  Community = c("Whole community", "Anti-Bd community"),
  Samples = 47L,
  Localities = 4L,
  Genera_before_filter = c(nrow(whole_relative), nrow(anti_relative)),
  Genera_after_20pct_filter = c(nrow(whole$filtered), nrow(anti$filtered)),
  ANCOMBC2_robust_among_top20 = c(
    sum(head(whole$table, 20)$ANCOMBC2_robust_FDR05),
    sum(head(anti$table, 20)$ANCOMBC2_robust_FDR05)
  ),
  ALDEx2_FDR05_among_top20 = c(
    sum(head(whole$table, 20)$ALDEx2_FDR05),
    sum(head(anti$table, 20)$ALDEx2_FDR05)
  ),
  Both_methods_among_top20 = c(
    sum(head(whole$table, 20)$Cross_method_support == "ANCOM-BC2 + ALDEx2"),
    sum(head(anti$table, 20)$Cross_method_support == "ANCOM-BC2 + ALDEx2")
  ),
  stringsAsFactors = FALSE
)
write.csv(summary_table, file.path(result_dir, "analysis_summary.csv"), row.names = FALSE)
write.csv(data.frame(Package = c("R", "phyloseq", "ggplot2"),
                     Version = c(R.version.string,
                                 as.character(packageVersion("phyloseq")),
                                 as.character(packageVersion("ggplot2")))),
          file.path(result_dir, "software_versions.csv"), row.names = FALSE)

cat("KSM analysis complete.\n")
print(summary_table)
cat("\nWhole-community top 10:\n")
print(whole$table[1:10, c("KSM_rank", "Genus", "KSM", "Top20_frequency",
                          "Is_anti_Bd", "Cross_method_support")], row.names = FALSE)
cat("\nAnti-Bd top 10:\n")
print(anti$table[1:10, c("KSM_rank", "Genus", "KSM", "Top20_frequency",
                         "Cross_method_support")], row.names = FALSE)
