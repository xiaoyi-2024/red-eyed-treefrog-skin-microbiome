# Red-eyed tree frog skin microbiome analysis

This repository contains the R workflow, statistical outputs, source-data tables for figures, and PNG figures used to analyse skin bacterial communities from 47 *Agalychnis callidryas* individuals sampled at four localities in Costa Rica and Panama.

## Repository structure

- `code/data_process/`: sample-depth filtering, phyloseq reconstruction, SILVA taxonomy integration, and Woodhams anti-Bd annotation.
- `code/analysis/`: community composition, diversity, clustering, locality comparison, differential-abundance, spatial, and Bd-arrival analyses.
- `results/`: numerical outputs, source-data tables, and PNG figures.

## Data availability

Research data are not distributed in this GitHub repository. The processing and analysis scripts retain the expected repository-relative input paths. Users must obtain the corresponding raw sequences, reference databases, metadata, and processed phyloseq objects separately before rerunning the workflow.

## Reproduction

Run commands from the repository root. Data-processing scripts are ordered numerically in `code/data_process/`. Analysis scripts in `code/analysis/` use repository-relative paths and write outputs under `results/`.

Package requirements are declared near the beginning of each R script. Random seeds and permutation counts are specified in the corresponding scripts.

## Functional annotation caveat

Candidate anti-Bd ASVs were defined by BLASTn similarity to inhibitory isolates in the Woodhams amphibian-skin antifungal isolate database, using at least 97% sequence identity and 90% query coverage. This annotation indicates sequence similarity to cultured inhibitory isolates; it does not demonstrate antifungal activity of the local strains in vivo.
