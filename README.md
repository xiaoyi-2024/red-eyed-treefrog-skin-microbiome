# Red-eyed tree frog skin microbiome analysis

This repository contains the R workflow, processed ASV-level data, statistical outputs, source-data tables for figures, and PNG figures used to analyse skin bacterial communities from 47 *Agalychnis callidryas* individuals sampled at four localities in Costa Rica and Panama.

## Repository structure

- `code/data_process/`: sample-depth filtering, phyloseq reconstruction, SILVA taxonomy integration, and Woodhams anti-Bd annotation.
- `code/analysis/`: community composition, diversity, clustering, locality comparison, differential-abundance, spatial, and Bd-arrival analyses.
- `data/processed/`: compact whole-community and anti-Bd phyloseq objects, ASV sequences, taxonomy, and functional annotations.
- `results/`: numerical outputs, source-data tables, and PNG figures.

## Data availability

The repository includes processed ASV-level data required for the principal community analyses:

- `data/processed/ps_clean_47_rebuilt_from_01_filtered.rds`: quality-controlled whole-community phyloseq object for 47 frog samples.
- `data/processed/blast_97id_90cov/ps_anti_bd_47_rebuilt_97id_90cov.rds`: candidate anti-Bd phyloseq subset.
- `data/processed/blast_97id_90cov/integrated_asv_taxonomy_function_rebuilt.csv`: ASV sequences, SILVA taxonomy, BLAST statistics, and functional annotations.
- `data/processed/blast_97id_90cov/rebuilt_47_asv_sequences.fasta`: ASV representative sequences.

Raw FASTQ files, filtered FASTQ files, the original sequencing archive, the SILVA training database, BLAST database files, and large DADA2 intermediate objects are not distributed through GitHub.

## Reproduction

Run commands from the repository root. Data-processing scripts are ordered numerically in `code/data_process/`. Analysis scripts in `code/analysis/` use repository-relative paths and write outputs under `results/`.

Package requirements are declared near the beginning of each R script. Random seeds and permutation counts are specified in the corresponding scripts.

## Functional annotation caveat

Candidate anti-Bd ASVs were defined by BLASTn similarity to inhibitory isolates in the Woodhams amphibian-skin antifungal isolate database, using at least 97% sequence identity and 90% query coverage. This annotation indicates sequence similarity to cultured inhibitory isolates; it does not demonstrate antifungal activity of the local strains in vivo.
