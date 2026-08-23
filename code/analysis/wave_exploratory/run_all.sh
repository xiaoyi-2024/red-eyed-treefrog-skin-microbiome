#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$ROOT"
mkdir -p results/wave_exploratory
: > results/wave_exploratory/full_run.log
echo "Running the single directed Bd-wave model" | tee -a results/wave_exploratory/full_run.log
Rscript code/analysis/wave_exploratory/make_wave_exploratory_figure.R 2>&1 |
  tee -a results/wave_exploratory/full_run.log
echo "Running the original v3_dow first-panel reproduction" | tee -a results/wave_exploratory/full_run.log
Rscript code/analysis/wave_exploratory/make_wave_exploratory_figure_v3_dow_first_panel.R 2>&1 |
  tee -a results/wave_exploratory/full_run.log
echo "Preparing beta-binomial and Hill-D1 locality metrics" | tee -a results/wave_exploratory/full_run.log
Rscript code/analysis/wave_exploratory/weighted_site_regressions.R 2>&1 |
  tee -a results/wave_exploratory/full_run.log
echo "Testing the Bd-arrival hypothesis with locality-level exact inference" | tee -a results/wave_exploratory/full_run.log
Rscript code/analysis/wave_exploratory/test_bd_arrival_hypothesis.R 2>&1 |
  tee -a results/wave_exploratory/full_run.log
echo "Analysing anti-Bd abundance and whole/anti-Bd Hill D0 against Bd-arrival year" | tee -a results/wave_exploratory/full_run.log
Rscript code/analysis/wave_exploratory/analyse_bd_arrival_abundance_and_d0.R 2>&1 |
  tee -a results/wave_exploratory/full_run.log
echo "Running PAM clustering of all samples" | tee -a results/wave_exploratory/full_run.log
Rscript code/analysis/04_pam_clustering_all_samples.R 2>&1 |
  tee -a results/wave_exploratory/full_run.log
echo "Running Bray-Curtis hierarchical clustering without a fixed k" | tee -a results/wave_exploratory/full_run.log
Rscript code/analysis/05_bray_curtis_hierarchical_clustering.R 2>&1 |
  tee -a results/wave_exploratory/full_run.log
echo "Running pairwise locality Bray-Curtis distance analyses" | tee -a results/wave_exploratory/full_run.log
Rscript code/analysis/06_pairwise_locality_bray_distance.R 2>&1 |
  tee -a results/wave_exploratory/full_run.log
echo "Running Aitchison centroid bootstrap and pairwise MiRKAT" | tee -a results/wave_exploratory/full_run.log
Rscript code/analysis/07_aitchison_centroid_and_pairwise_mirkat.R 2>&1 |
  tee -a results/wave_exploratory/full_run.log
echo "Running pairwise CLR and robust Aitchison comparisons" | tee -a results/wave_exploratory/full_run.log
Rscript code/analysis/08_pairwise_aitchison_comparison.R 2>&1 |
  tee -a results/wave_exploratory/full_run.log
echo "Running anti-Bd PAM, Bray-Curtis clustering, and locality distance analyses" | tee -a results/wave_exploratory/full_run.log
bash code/analysis/09_run_antibd_clustering_distances.sh 2>&1 |
  tee -a results/wave_exploratory/full_run.log
echo "Comparing whole-community and anti-Bd Bray-Curtis distances" | tee -a results/wave_exploratory/full_run.log
Rscript code/analysis/10_compare_whole_antibd_bray_distance.R 2>&1 |
  tee -a results/wave_exploratory/full_run.log
echo "Comparing whole-community and anti-Bd Aitchison distances" | tee -a results/wave_exploratory/full_run.log
Rscript code/analysis/11_compare_whole_antibd_aitchison_distances.R 2>&1 |
  tee -a results/wave_exploratory/full_run.log
echo "Single Bd-wave model completed." | tee -a results/wave_exploratory/full_run.log
