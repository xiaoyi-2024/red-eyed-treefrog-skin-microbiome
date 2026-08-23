# Red-eyed tree frog skin microbiome analysis

This repository contains the reproducible R workflow, processed data objects, source tables, statistical outputs, and PNG figures used to analyse skin bacterial communities from 47 *Agalychnis callidryas* individuals sampled at four localities in Costa Rica and Panama.

## Repository structure

- `code/data_process/`: sample-depth filtering, phyloseq reconstruction, SILVA taxonomy integration, and Woodhams anti-Bd annotation.
- `code/analysis/`: community composition, diversity, clustering, locality comparison, differential-abundance, spatial, and Bd-arrival analyses.
- `data/processed/`: compact phyloseq objects and derived analysis inputs.
- `data/reference/`: small reference tables and the Woodhams amphibian-skin isolate sequences.
- `results/`: numerical outputs, source-data tables, and PNG figures.

## Data not stored on GitHub

Raw FASTQ archives, the SILVA 138.1 training database, and large DADA2 intermediate objects are excluded because of their size. They can be regenerated or obtained from the original data sources. Their expected paths are documented in `.gitignore` and in the processing scripts.

## Reproduction

Run commands from the repository root. Data-processing scripts are ordered numerically in `code/data_process/`. Analysis scripts in `code/analysis/` use repository-relative paths and write outputs under `results/`.

The principal compact inputs are:

- `data/processed/ps_clean_47_rebuilt_from_01_filtered.rds`
- `data/processed/blast_97id_90cov/ps_anti_bd_47_rebuilt_97id_90cov.rds`

Package requirements are declared near the beginning of each R script. Random seeds and permutation counts are specified in the corresponding scripts.

## Functional annotation caveat

Candidate anti-Bd ASVs were defined by BLASTn similarity to inhibitory isolates in the Woodhams amphibian-skin antifungal isolate database, using at least 97% sequence identity and 90% query coverage. This annotation indicates sequence similarity to cultured inhibitory isolates; it does not demonstrate antifungal activity of the local strains in vivo.

## Citation

Citation information will be added when the associated thesis or manuscript is publicly available.
