suppressPackageStartupMessages({
  library(phyloseq)
  library(iNEXT)
  library(ggplot2)
  library(svglite)
  library(ragg)
})

set.seed(20260816)
result_dir <- "results/bd_arrival_community"
figure_dir <- file.path(result_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

full_file <- "data/processed/ps_clean_47_rebuilt_from_01_filtered.rds"
wave_file <- "results/wave_exploratory/site_metrics_single_wave.csv"
stopifnot(file.exists(full_file), file.exists(wave_file))

ps <- readRDS(full_file)
counts <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) counts <- t(counts)
md <- data.frame(sample_data(ps), check.names = FALSE,
                 stringsAsFactors = FALSE)[rownames(counts), , drop = FALSE]

site_key <- c(
  "Lost Iguana Resort, La Fortuna, San Carlos" = "Lost Iguana",
  "Veragua Rainforest, Las Brisas de Veragua" = "Veragua",
  "PN Altos de Campana" = "Altos de Campana",
  "PN Soberanía, Las Cumbres" = "Soberanía"
)
site_levels <- c("Lost Iguana", "Veragua", "Altos de Campana", "Soberanía")
site <- factor(unname(site_key[md$Locality]), levels = site_levels)
stopifnot(!anyNA(site), nrow(counts) == 47L)

wave <- read.csv(wave_file, stringsAsFactors = FALSE, check.names = FALSE)
site_year <- setNames(wave$pred_wave_year, wave$Site)
year <- unname(site_year[as.character(site)])
stopifnot(!anyNA(year))

# Genus-level assemblages. Unclassified genera are retained as family-specific
# groups so that unclassified reads are not silently discarded.
tax <- as(tax_table(ps), "matrix")[colnames(counts), , drop = FALSE]
genus <- as.character(tax[, "Genus"])
family <- as.character(tax[, "Family"])
family[is.na(family) | trimws(family) == ""] <- "Bacteria"
missing_genus <- is.na(genus) | trimws(genus) == ""
genus[missing_genus] <- paste0("Unclassified_", family[missing_genus])
genus_counts <- t(rowsum(t(counts), group = genus, reorder = FALSE))

assemblages <- lapply(seq_len(nrow(genus_counts)), function(i) {
  z <- genus_counts[i, ]
  unname(z[z > 0])
})
names(assemblages) <- rownames(genus_counts)
coverage_info <- iNEXT::DataInfo(assemblages, datatype = "abundance")
common_coverage <- min(coverage_info$SC)
d1_est <- iNEXT::estimateD(assemblages, datatype = "abundance", base = "coverage",
                           level = common_coverage, q = 1, nboot = 0)
d1 <- setNames(d1_est$qD, d1_est$Assemblage)[rownames(counts)]
stopifnot(all(is.finite(d1)))

sample_results <- data.frame(
  Sample_ID = rownames(counts), Bd_arrival_year = year,
  Whole_community_genus_Hill_D1 = unname(d1),
  Locality_internal = as.character(site), stringsAsFactors = FALSE
)
write.csv(sample_results,
          file.path(result_dir, "whole_community_hill_d1_sample_values.csv"),
          row.names = FALSE)
write.csv(coverage_info,
          file.path(result_dir, "whole_community_genus_sample_coverage.csv"),
          row.names = FALSE)

bootstrap_summary <- function(idx, B = 4999L) {
  x <- d1[idx]
  boots <- replicate(B, mean(sample(x, length(x), replace = TRUE)))
  data.frame(
    Bd_arrival_year = unique(year[idx]), Sample_n = length(idx),
    Mean_D1 = mean(x), SE_D1 = sd(x) / sqrt(length(x)),
    CI_low = unname(quantile(boots, 0.025)),
    CI_high = unname(quantile(boots, 0.975))
  )
}
site_summary <- do.call(rbind, lapply(site_levels, function(s) {
  bootstrap_summary(which(site == s))
}))
site_summary <- site_summary[order(site_summary$Bd_arrival_year), ]
write.csv(site_summary,
          file.path(result_dir, "whole_community_hill_d1_year_summary.csv"),
          row.names = FALSE)

permutations <- function(x) {
  if (length(x) == 1L) return(matrix(x, nrow = 1L))
  do.call(rbind, lapply(seq_along(x), function(i) cbind(x[i], permutations(x[-i]))))
}
perm_years <- permutations(unname(site_year[site_levels]))

# Prespecified multiplicative trend: a straight line for log(Hill D1), weighted
# by the precision of each year-specific mean. The year labels are permuted at
# locality level, yielding all 4! = 24 exact assignments.
x10 <- (site_summary$Bd_arrival_year - mean(site_summary$Bd_arrival_year)) / 10
weights <- 1 / pmax(site_summary$SE_D1, .Machine$double.eps)^2
null_fit <- lm(log(Mean_D1) ~ 1, data = site_summary, weights = weights)
alt_fit <- lm(log(Mean_D1) ~ x10, data = site_summary, weights = weights)
obs_stat <- max(0, deviance(null_fit) - deviance(alt_fit))
perm_stat <- numeric(nrow(perm_years))
perm_beta <- numeric(nrow(perm_years))
for (i in seq_len(nrow(perm_years))) {
  # site_summary is ordered by year, so map the permuted locality years to the
  # corresponding locality means before re-ordering is unnecessary for lm().
  mean_by_site <- vapply(site_levels, function(s) mean(d1[site == s]), numeric(1))
  se_by_site <- vapply(site_levels, function(s) sd(d1[site == s]) / sqrt(sum(site == s)), numeric(1))
  xp <- (perm_years[i, ] - mean(site_year)) / 10
  fp <- lm(log(mean_by_site) ~ xp, weights = 1 / pmax(se_by_site, .Machine$double.eps)^2)
  f0 <- lm(log(mean_by_site) ~ 1, weights = 1 / pmax(se_by_site, .Machine$double.eps)^2)
  perm_stat[i] <- max(0, deviance(f0) - deviance(fp))
  perm_beta[i] <- coef(fp)["xp"]
}
exact_p <- mean(perm_stat >= obs_stat - 1e-12)
beta <- unname(coef(alt_fit)["x10"])
ratio <- exp(beta)

test_result <- data.frame(
  Endpoint = "Whole-community genus Hill D1 diversity",
  Model = "Inverse-variance weighted log-linear trend",
  Effect_log_ratio_per_10_years = beta,
  Ratio_per_10_years = ratio,
  Percent_change_per_10_years = 100 * (ratio - 1),
  Statistic = obs_stat, Exact_p = exact_p,
  Locality_n = 4L, Frog_n = 47L, Exact_permutations = 24L,
  Common_sample_coverage = common_coverage,
  Conclusion_at_0_05 = ifelse(exact_p < 0.05,
                               "Evidence of association",
                               "Insufficient evidence to reject no association"),
  stringsAsFactors = FALSE
)
write.csv(test_result,
          file.path(result_dir, "whole_community_hill_d1_bd_year_test.csv"),
          row.names = FALSE)
write.csv(data.frame(Permutation = seq_len(nrow(perm_years)),
                     Statistic = perm_stat,
                     Coefficient_per_10_years = perm_beta),
          file.path(result_dir, "whole_community_hill_d1_exact_permutations.csv"),
          row.names = FALSE)

grid <- data.frame(Bd_arrival_year = seq(min(year), max(year), length.out = 300))
grid$x10 <- (grid$Bd_arrival_year - mean(site_summary$Bd_arrival_year)) / 10
grid$fit <- exp(predict(alt_fit, newdata = grid))

p <- ggplot(sample_results,
            aes(Bd_arrival_year, Whole_community_genus_Hill_D1)) +
  geom_point(position = position_jitter(width = 0.28, height = 0),
             colour = "#8F8F8F", alpha = 0.55, size = 1.5) +
  geom_errorbar(data = site_summary,
                aes(Bd_arrival_year, Mean_D1, ymin = CI_low, ymax = CI_high),
                inherit.aes = FALSE, width = 0.3, linewidth = 0.55,
                colour = "#166AA5") +
  geom_point(data = site_summary, aes(Bd_arrival_year, Mean_D1),
             inherit.aes = FALSE, colour = "#166AA5", size = 3) +
  geom_line(data = grid, aes(Bd_arrival_year, fit), inherit.aes = FALSE,
            colour = "#343434", linewidth = 0.9) +
  scale_x_continuous(breaks = seq(1985, 2010, 5)) +
  labs(
    title = expression("Whole-community Hill "^1*D*" and Bd-arrival year"),
    subtitle = sprintf("Weighted log-linear fit: %.1f%% per 10 years; exact p = %.3f",
                       100 * (ratio - 1), exact_p),
    x = "Estimated Bd arrival year",
    y = expression("Hill "^1*D*" (coverage-standardized)"),
    caption = paste0("Grey points: individual frogs (n = 47). Blue points: year-specific means and bootstrap 95% CIs.\n",
                     "The fitted line uses four independent locality-level means.\n",
                     "Exact p uses all 24 assignments of years among localities.")
  ) +
  theme_classic(base_size = 8, base_family = "sans") +
  theme(axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4),
        plot.title = element_text(size = 9.5, face = "bold"),
        plot.subtitle = element_text(size = 7.2),
        plot.caption = element_text(size = 6.1, hjust = 0))

stem <- file.path(figure_dir, "figure_whole_community_hill_d1_bd_year_fit")
w <- 120 / 25.4; h <- 100 / 25.4
ggsave(paste0(stem, ".png"), p, width = w, height = h, dpi = 600, bg = "white")
svglite::svglite(paste0(stem, ".svg"), width = w, height = h); print(p); dev.off()
grDevices::pdf(paste0(stem, ".pdf"), width = w, height = h,
               family = "Helvetica", useDingbats = FALSE); print(p); dev.off()
ragg::agg_tiff(paste0(stem, ".tiff"), width = w, height = h,
               units = "in", res = 600, background = "white"); print(p); dev.off()

cat("Whole-community Hill D1 trend\n")
print(test_result)
