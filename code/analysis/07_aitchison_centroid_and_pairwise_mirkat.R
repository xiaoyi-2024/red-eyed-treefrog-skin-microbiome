suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
})

if (!requireNamespace("MiRKAT", quietly = TRUE)) {
  stop("The CRAN package MiRKAT is required: install.packages('MiRKAT')")
}

set.seed(20260814)
input_file <- "data/processed/ps_clean_47_rebuilt_from_01_filtered.rds"
result_dir <- "results/aitchison_centroid_mirkat"
figure_dir <- file.path(result_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

ps <- readRDS(input_file)
otu <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) otu <- t(otu)
tax <- as(tax_table(ps), "matrix")[colnames(otu), , drop = FALSE]
meta <- data.frame(sample_data(ps), check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(nrow(otu) == 47L, identical(rownames(otu), rownames(meta)))

genus <- as.character(tax[, "Genus"])
missing <- is.na(genus) | trimws(genus) == ""
family <- as.character(tax[, "Family"])
family[is.na(family) | trimws(family) == ""] <- "Bacteria"
genus[missing] <- paste0("Unclassified_", family[missing])
genus_counts <- t(rowsum(t(otu), group = genus, reorder = FALSE))
prevalence_min <- ceiling(0.20 * nrow(genus_counts))
keep <- colSums(genus_counts > 0) >= prevalence_min
filtered <- genus_counts[, keep, drop = FALSE]

site_key <- c(
  "Lost Iguana Resort, La Fortuna, San Carlos" = "Lost Iguana",
  "Veragua Rainforest, Las Brisas de Veragua" = "Veragua",
  "PN Altos de Campana" = "Altos de Campana",
  "PN Soberanía, Las Cumbres" = "Soberanía"
)
site_levels <- c("Lost Iguana", "Veragua", "Altos de Campana", "Soberanía")
site <- factor(unname(site_key[meta$Locality]), levels = site_levels)
stopifnot(!anyNA(site))

# Ordinary Aitchison coordinates use the same pseudocount (1) and CLR
# definition as the primary community-composition analysis.
clr <- log(filtered + 1)
clr <- clr - rowMeans(clr)

# vegan's rclr performs robust CLR with low-rank matrix completion; Euclidean
# distances between these coordinates equal robust Aitchison distances.
rclr <- vegan::decostand(filtered, method = "rclr")
stopifnot(!anyNA(rclr), all(is.finite(rclr)))

site_pairs <- combn(site_levels, 2, simplify = FALSE)
B <- 9999L
bootstrap_pair <- function(coords, pair, metric, seed) {
  i1 <- which(site == pair[1]); i2 <- which(site == pair[2])
  point <- sqrt(sum((colMeans(coords[i1, , drop = FALSE]) -
                     colMeans(coords[i2, , drop = FALSE]))^2))
  set.seed(seed)
  boot <- replicate(B, {
    c1 <- colMeans(coords[sample(i1, length(i1), replace = TRUE), , drop = FALSE])
    c2 <- colMeans(coords[sample(i2, length(i2), replace = TRUE), , drop = FALSE])
    sqrt(sum((c1 - c2)^2))
  })
  percentile_ci <- unname(quantile(boot, c(0.025, 0.975), type = 7))
  # Norm distances are upward biased when noisy high-dimensional centroids are
  # re-estimated. Use the bootstrap distribution to estimate the SE and centre
  # the displayed normal interval on the observed distance. Basic and raw
  # percentile intervals are retained in the CSV for auditing that bias.
  basic_ci <- c(max(0, 2 * point - percentile_ci[2]),
                2 * point - percentile_ci[1])
  boot_se <- sd(boot)
  normal_ci <- c(max(0, point - qnorm(0.975) * boot_se),
                 point + qnorm(0.975) * boot_se)
  data.frame(
    Metric = metric, Site_1 = pair[1], Site_2 = pair[2],
    Pair = paste(pair, collapse = " vs "),
    N_site_1 = length(i1), N_site_2 = length(i2),
    Centroid_distance = point, Bootstrap_SE = boot_se,
    Bootstrap_bias = mean(boot) - point,
    CI95_low = normal_ci[1], CI95_high = normal_ci[2],
    Basic_CI95_low = basic_ci[1], Basic_CI95_high = basic_ci[2],
    Percentile_CI95_low = percentile_ci[1],
    Percentile_CI95_high = percentile_ci[2],
    CI_method = "normal bootstrap using within-locality bootstrap SE",
    Bootstrap_replicates = B,
    Bootstrap_unit = "frog samples resampled independently within locality"
  )
}

centroid_results <- do.call(rbind, c(
  lapply(seq_along(site_pairs), function(i)
    bootstrap_pair(clr, site_pairs[[i]], "CLR-Aitchison", 20260814 + i)),
  lapply(seq_along(site_pairs), function(i)
    bootstrap_pair(rclr, site_pairs[[i]], "Robust Aitchison", 20260914 + i))
))
write.csv(centroid_results,
          file.path(result_dir, "pairwise_aitchison_centroid_bootstrap.csv"),
          row.names = FALSE)

# Convert Euclidean distance matrices into centered positive-semidefinite
# kernels. Both CLR and rCLR distances are Euclidean by construction.
distance_to_kernel <- function(distance) {
  D <- as.matrix(distance); n <- nrow(D)
  H <- diag(n) - matrix(1 / n, n, n)
  K <- -0.5 * H %*% (D^2) %*% H
  K <- (K + t(K)) / 2
  eig <- eigen(K, symmetric = TRUE, only.values = TRUE)$values
  if (min(eig) < -1e-7 * max(abs(eig))) stop("Kernel is not positive semidefinite")
  K
}
K_standard <- distance_to_kernel(dist(clr))
K_robust <- distance_to_kernel(dist(rclr))
# Standard MiRKAT supports a binary outcome. Each of the six locality pairs is
# therefore analysed separately. Permutation inference is used because every
# pair has fewer than 50 independent frogs. BH correction is applied across the
# six locality comparisons separately for each distance kernel.
nperm_mirkat <- 9999L
pairwise_mirkat <- do.call(rbind, lapply(seq_along(site_pairs), function(i) {
  pair <- site_pairs[[i]]
  take <- which(site %in% pair)
  y <- as.integer(site[take] == pair[2])
  pair_kernels <- list(
    `CLR-Aitchison` = K_standard[take, take, drop = FALSE],
    `Robust Aitchison` = K_robust[take, take, drop = FALSE]
  )
  set.seed(20261014 + i)
  fit <- MiRKAT::MiRKAT(
    y = y, X = NULL, Ks = pair_kernels, out_type = "D",
    method = "permutation", omnibus = "permutation", nperm = nperm_mirkat,
    returnKRV = TRUE, returnR2 = TRUE
  )
  data.frame(
    Site_1 = pair[1], Site_2 = pair[2], Pair = paste(pair, collapse = " vs "),
    N_site_1 = sum(y == 0), N_site_2 = sum(y == 1),
    Kernel = names(pair_kernels),
    p = as.numeric(fit$p_values[names(pair_kernels)]),
    KRV = as.numeric(fit$KRV[names(pair_kernels)]),
    R2 = as.numeric(fit$R2[names(pair_kernels)]),
    Permutations = nperm_mirkat,
    Method = "Binary MiRKAT; label permutation; intercept-only model",
    stringsAsFactors = FALSE
  )
}))
pairwise_mirkat$p_FDR_BH <- ave(
  pairwise_mirkat$p, pairwise_mirkat$Kernel,
  FUN = function(x) p.adjust(x, method = "BH")
)
pairwise_mirkat$Significant_FDR_0.05 <- pairwise_mirkat$p_FDR_BH < 0.05
write.csv(pairwise_mirkat, file.path(result_dir, "pairwise_mirkat_results.csv"),
          row.names = FALSE)
write.csv(data.frame(
  Parameter = c("Samples", "Original_ASVs", "Aggregated_genera", "Retained_genera",
                "Prevalence_min_samples", "CLR_pseudocount", "Bootstrap_replicates",
                "MiRKAT_permutations", "vegan_version", "MiRKAT_version"),
  Value = c(nrow(otu), ncol(otu), ncol(genus_counts), ncol(filtered), prevalence_min,
            1, B, nperm_mirkat, as.character(packageVersion("vegan")),
            as.character(packageVersion("MiRKAT")))
), file.path(result_dir, "aitchison_mirkat_parameters.csv"), row.names = FALSE)

short_site <- c("Lost Iguana" = "Lost", "Veragua" = "Veragua",
                "Altos de Campana" = "Campana", "Soberanía" = "Soberanía")
short_pair <- function(x) {
  z <- strsplit(as.character(x), " vs ", fixed = TRUE)
  vapply(z, function(v) paste(short_site[v], collapse = " vs "), character(1))
}
centroid_results$Pair_short <- short_pair(centroid_results$Pair)
theme_pub <- theme_classic(base_size = 8, base_family = "sans") +
  theme(axis.line = element_line(linewidth = 0.4), axis.ticks = element_line(linewidth = 0.4),
        plot.title = element_text(size = 9, face = "bold"),
        plot.subtitle = element_text(size = 7), legend.position = "none")

forest_panel <- function(data, title, colour) {
  data <- data[order(data$Centroid_distance), ]
  data$Pair_short <- factor(data$Pair_short, levels = rev(data$Pair_short))
  ggplot(data, aes(Pair_short, Centroid_distance)) +
    geom_errorbar(aes(ymin = CI95_low, ymax = CI95_high), width = 0.16,
                  colour = colour, linewidth = 0.65) +
    geom_point(shape = 21, size = 2.8, stroke = 0.6, fill = "white", colour = colour) +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.08))) +
    labs(title = title, subtitle = "95% CI from within-locality bootstrap SE",
         x = NULL, y = "Centroid distance") + theme_pub
}
p_standard <- forest_panel(
  centroid_results[centroid_results$Metric == "CLR-Aitchison", ],
  "a  CLR-Aitchison centroid distances", "#7A1F5C"
)
p_robust <- forest_panel(
  centroid_results[centroid_results$Metric == "Robust Aitchison", ],
  "b  Robust Aitchison centroid distances", "#0072B2"
)
forest <- (p_standard | p_robust) +
  plot_annotation(
    title = "Pairwise locality differences in microbiome composition",
    subtitle = "Genus level; 20% prevalence filter; 9,999 bootstrap resamples within each locality.",
    caption = "Distances are metric-specific and should not be compared numerically between panels.",
    theme = theme(plot.title = element_text(size = 11, face = "bold"),
                  plot.subtitle = element_text(size = 7.2),
                  plot.caption = element_text(size = 6.2, hjust = 0))
  )

pairwise_mirkat$Pair_short <- short_pair(pairwise_mirkat$Pair)
pair_order <- unique(pairwise_mirkat$Pair_short[order(pairwise_mirkat$p_FDR_BH, decreasing = TRUE)])
pairwise_mirkat$Pair_short <- factor(pairwise_mirkat$Pair_short, levels = pair_order)
p_mirkat <- ggplot(pairwise_mirkat,
                   aes(Pair_short, -log10(p_FDR_BH), colour = Kernel, shape = Kernel)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey35") +
  geom_point(size = 2.7, position = position_dodge(width = 0.42)) +
  coord_flip() +
  scale_colour_manual(values = c("CLR-Aitchison" = "#7A1F5C",
                                 "Robust Aitchison" = "#0072B2")) +
  scale_shape_manual(values = c("CLR-Aitchison" = 16, "Robust Aitchison" = 17)) +
  labs(title = "Pairwise MiRKAT tests of locality-associated composition",
       subtitle = "Binary tests; 9,999 label permutations; BH correction within each kernel",
       x = NULL, y = expression(-log[10](italic(p)[FDR])), colour = NULL, shape = NULL) +
  theme_pub + theme(legend.position = "top")

save_plot <- function(plot, stem, width_mm, height_mm) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  ggsave(paste0(stem, ".png"), plot, width = w, height = h, dpi = 600, bg = "white")
  ragg::agg_tiff(paste0(stem, ".tiff"), width = w, height = h, units = "in",
                 res = 600, background = "white"); print(plot); dev.off()
  svglite::svglite(paste0(stem, ".svg"), width = w, height = h); print(plot); dev.off()
  grDevices::pdf(paste0(stem, ".pdf"), width = w, height = h,
                 family = "Helvetica", useDingbats = FALSE); print(plot); dev.off()
}
save_plot(forest, file.path(figure_dir, "figure_aitchison_centroid_forest"), 183, 100)
save_plot(p_mirkat, file.path(figure_dir, "figure_pairwise_mirkat"), 150, 92)
print(centroid_results)
print(pairwise_mirkat)
