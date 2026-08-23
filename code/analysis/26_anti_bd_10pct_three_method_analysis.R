#!/usr/bin/env Rscript

# Reanalysis of the anti-Bd genus subset using the ANCOM-BC2 default 10%
# prevalence filter. KSM and ALDEx2 use the identical 10% feature universe so
# that the three methods remain directly comparable.

suppressPackageStartupMessages({
  library(phyloseq)
  library(ANCOMBC)
  library(ALDEx2)
})

set.seed(20260816)
out <- "results/anti_bd_10pct_three_method"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
whole_file <- "data/processed/ps_clean_47_rebuilt_from_01_filtered.rds"
anti_file <- "data/processed/blast_97id_90cov/ps_anti_bd_47_rebuilt_97id_90cov.rds"
stopifnot(all(file.exists(c(whole_file, anti_file))))

site_levels <- c("Lost Iguana", "Veragua", "Altos de Campana", "Soberanía")
site_key <- c(
  "Lost Iguana Resort, La Fortuna, San Carlos" = "Lost Iguana",
  "Veragua Rainforest, Las Brisas de Veragua" = "Veragua",
  "PN Altos de Campana" = "Altos de Campana",
  "PN Soberanía, Las Cumbres" = "Soberanía"
)

otu_rows <- function(ps) {
  x <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) x <- t(x)
  x
}
genus_counts <- function(ps) {
  x <- otu_rows(ps)
  tx <- as(tax_table(ps), "matrix")[rownames(x), , drop = FALSE]
  g <- as.character(tx[, "Genus"])
  f <- as.character(tx[, "Family"])
  f[is.na(f) | trimws(f) == ""] <- "Bacteria"
  miss <- is.na(g) | trimws(g) == ""
  g[miss] <- paste0("Unclassified_", f[miss])
  rowsum(x, g, reorder = FALSE)
}

ps_whole <- readRDS(whole_file)
ps_anti <- readRDS(anti_file)
whole <- genus_counts(ps_whole)
anti <- genus_counts(ps_anti)
ids <- colnames(whole)
anti <- anti[, ids, drop = FALSE]
meta <- data.frame(sample_data(ps_whole), check.names = FALSE,
                   stringsAsFactors = FALSE)[ids, , drop = FALSE]
meta$Site <- factor(unname(site_key[meta$Locality]), levels = site_levels)
stopifnot(ncol(anti) == 47L, !anyNA(meta$Site))

prevalence <- rowMeans(anti > 0)
keep10 <- prevalence >= 0.10
filtered <- anti[keep10, , drop = FALSE]
stopifnot(nrow(filtered) == 47L)
filter_summary <- data.frame(
  Threshold = c(0.10, 0.20),
  Minimum_detected_samples = c(ceiling(47 * 0.10), ceiling(47 * 0.20)),
  Retained_genera = c(sum(prevalence >= 0.10), sum(prevalence >= 0.20)),
  Total_genera = nrow(anti),
  Retained_anti_Bd_reads_percent = c(
    100 * sum(anti[prevalence >= 0.10, ]) / sum(anti),
    100 * sum(anti[prevalence >= 0.20, ]) / sum(anti)
  )
)
write.csv(filter_summary, file.path(out, "filter_10pct_vs_20pct.csv"), row.names = FALSE)
write.csv(data.frame(Genus = rownames(filtered), filtered, check.names = FALSE),
          file.path(out, "anti_bd_genus_counts_10pct.csv"), row.names = FALSE)

# ANCOM-BC2 receives the unfiltered anti-Bd genus matrix and applies its
# documented default prevalence threshold directly (prv_cut = 0.10).
checkpoint <- file.path(out, "ancombc2_10pct_complete_result.rds")
if (file.exists(checkpoint)) {
  fit <- readRDS(checkpoint)
} else {
  fit <- ANCOMBC::ancombc2(
    data = anti, meta_data = meta, taxa_are_rows = TRUE,
    fix_formula = "Site", group = "Site",
    prv_cut = 0.10, lib_cut = 0, p_adj_method = "BH",
    struc_zero = TRUE, neg_lb = TRUE, global = TRUE,
    pairwise = FALSE, pseudo_sens = TRUE, alpha = 0.05,
    n_cl = 1, verbose = FALSE
  )
  saveRDS(fit, checkpoint)
}
write.csv(fit$res_global, file.path(out, "ancombc2_10pct_global_test.csv"),
          row.names = FALSE)
write.csv(fit$res, file.path(out, "ancombc2_10pct_coefficients.csv"),
          row.names = FALSE)
write.csv(fit$ss_tab, file.path(out, "ancombc2_10pct_pseudocount_sensitivity.csv"),
          row.names = FALSE)
write.csv(fit$zero_ind, file.path(out, "ancombc2_10pct_structural_zeros.csv"),
          row.names = FALSE)

# ALDEx2 sensitivity analysis on the same 47 retained genera.
aldex_file <- file.path(out, "aldex2_10pct_global_kw.csv")
if (file.exists(aldex_file)) {
  aldex <- read.csv(aldex_file, stringsAsFactors = FALSE, check.names = FALSE)
} else {
  clr <- ALDEx2::aldex.clr(filtered, as.character(meta$Site), mc.samples = 128,
                          denom = "all", verbose = FALSE)
  aldex <- ALDEx2::aldex.kw(clr, verbose = FALSE)
  aldex$taxon <- rownames(aldex)
  rownames(aldex) <- NULL
  write.csv(aldex, aldex_file, row.names = FALSE)
}

ks_dist <- function(x, y) {
  nx <- length(x); ny <- length(y)
  z <- c(x, y); ord <- order(z); zs <- z[ord]
  fx <- cumsum(ord <= nx) / nx; fy <- cumsum(ord > nx) / ny
  last <- !duplicated(zs, fromLast = TRUE)
  max(abs(fx[last] - fy[last]))
}
ksm_one <- function(x, group) {
  lev <- levels(group); n <- length(x); k <- length(lev)
  pairs <- combn(lev, 2, simplify = FALSE)
  sum(vapply(pairs, function(p) {
    a <- x[group == p[1]]; b <- x[group == p[2]]
    ((length(a) + length(b)) / (n * (k - 1))) * ks_dist(a, b)
  }, numeric(1)))
}
ksm_vec <- function(mat, group) vapply(seq_len(nrow(mat)), function(i)
  ksm_one(mat[i, ], group), numeric(1))

relative <- sweep(filtered, 2, colSums(whole), "/")
observed <- ksm_vec(relative, meta$Site)

# 1,000 locality-stratified bootstrap samples.
site_ids <- split(seq_along(meta$Site), meta$Site, drop = TRUE)
one_boot <- function(iter) {
  idx <- unlist(lapply(site_ids, function(z) sample(z, length(z), replace = TRUE)),
                use.names = FALSE)
  g <- factor(as.character(meta$Site[idx]), levels = site_levels)
  ksm_vec(relative[, idx, drop = FALSE], g)
}
vals <- simplify2array(parallel::mclapply(
  1:1000, one_boot, mc.cores = max(1L, min(4L, parallel::detectCores(FALSE))),
  mc.set.seed = TRUE
))
rownames(vals) <- rownames(relative)
ranks <- apply(vals, 2, function(v) rank(-v, ties.method = "first"))
ksm <- data.frame(
  Genus = rownames(relative), KSM = observed,
  Bootstrap_CI_low = apply(vals, 1, quantile, 0.025, names = FALSE),
  Bootstrap_CI_high = apply(vals, 1, quantile, 0.975, names = FALSE),
  Top10_frequency = rowMeans(ranks <= 10),
  Top20_frequency = rowMeans(ranks <= 20),
  Prevalence = prevalence[rownames(relative)],
  Mean_relative_abundance_percent = 100 * rowMeans(relative),
  stringsAsFactors = FALSE
)

an <- fit$res_global
an2 <- data.frame(Genus = an$taxon, ANCOMBC2_W = an$W,
                  ANCOMBC2_p = an$p_val, ANCOMBC2_q = an$q_val,
                  ANCOMBC2_robust_FDR05 = an$diff_robust_abn %in% TRUE,
                  stringsAsFactors = FALSE)
al2 <- data.frame(Genus = aldex$taxon, ALDEx2_p = aldex$kw.ep,
                  ALDEx2_q = aldex$kw.eBH,
                  ALDEx2_FDR05 = !is.na(aldex$kw.eBH) & aldex$kw.eBH < 0.05,
                  stringsAsFactors = FALSE)
res <- merge(merge(ksm, an2, by = "Genus", all.x = TRUE, sort = FALSE),
             al2, by = "Genus", all.x = TRUE, sort = FALSE)
res$ANCOMBC2_robust_FDR05[is.na(res$ANCOMBC2_robust_FDR05)] <- FALSE
res$ALDEx2_FDR05[is.na(res$ALDEx2_FDR05)] <- FALSE
res$Evidence_class <- ifelse(
  res$ANCOMBC2_robust_FDR05 & res$ALDEx2_FDR05, "ANCOM-BC2 + ALDEx2",
  ifelse(res$ANCOMBC2_robust_FDR05, "ANCOM-BC2 only",
         ifelse(res$ALDEx2_FDR05, "ALDEx2 only", "KSM only"))
)
res <- res[order(-res$KSM, -res$Top10_frequency, res$Genus), ]
res$KSM_rank <- seq_len(nrow(res))
write.csv(res, file.path(out, "anti_bd_10pct_three_method_ranked.csv"),
          row.names = FALSE)

summary <- data.frame(
  Samples = 47L, Localities = 4L, Prevalence_threshold = 0.10,
  Minimum_detected_samples = 5L, Retained_genera = nrow(filtered),
  Retained_anti_Bd_reads_percent = filter_summary$Retained_anti_Bd_reads_percent[1],
  ANCOMBC2_robust_FDR05 = sum(res$ANCOMBC2_robust_FDR05),
  ALDEx2_FDR05 = sum(res$ALDEx2_FDR05),
  Both_methods = sum(res$ANCOMBC2_robust_FDR05 & res$ALDEx2_FDR05)
)
write.csv(summary, file.path(out, "analysis_summary.csv"), row.names = FALSE)
cat("Anti-Bd 10% reanalysis complete.\n")
print(summary)
print(head(res[, c("KSM_rank", "Genus", "KSM", "ANCOMBC2_q",
                   "ANCOMBC2_robust_FDR05", "ALDEx2_q", "ALDEx2_FDR05",
                   "Evidence_class")], 15), row.names = FALSE)
