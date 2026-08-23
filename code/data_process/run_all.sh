#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"
Rscript code/data_process/01_filter_samples_depth1000.R
Rscript code/data_process/02_build_phyloseq_from_filtered.R
Rscript code/data_process/03_blast_woodhams_97id_90cov.R
