#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"

ANTI_BD_INPUT="data/processed/blast_97id_90cov/ps_anti_bd_47_rebuilt_97id_90cov.rds"

PHYLOSEQ_INPUT="$ANTI_BD_INPUT" \
ANALYSIS_RESULT_DIR="results/anti_bd_pam_clustering" \
COMMUNITY_LABEL="anti-Bd community" \
Rscript code/analysis/04_pam_clustering_all_samples.R

PHYLOSEQ_INPUT="$ANTI_BD_INPUT" \
ANALYSIS_RESULT_DIR="results/anti_bd_bray_curtis_clustering" \
COMMUNITY_LABEL="anti-Bd bacterial community (47 frogs)" \
Rscript code/analysis/05_bray_curtis_hierarchical_clustering.R

PHYLOSEQ_INPUT="$ANTI_BD_INPUT" \
ANALYSIS_RESULT_DIR="results/anti_bd_pairwise_bray_distance" \
COMMUNITY_LABEL="anti-Bd bacterial community (47 frogs)" \
Rscript code/analysis/06_pairwise_locality_bray_distance.R

PHYLOSEQ_INPUT="$ANTI_BD_INPUT" \
ANALYSIS_RESULT_DIR="results/anti_bd_pairwise_aitchison_comparison" \
COMMUNITY_LABEL="anti-Bd bacterial community (47 frogs)" \
Rscript code/analysis/08_pairwise_aitchison_comparison.R
