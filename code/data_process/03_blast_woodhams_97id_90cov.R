#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(phyloseq); library(Biostrings)})

submission <- normalizePath(".", mustWork = TRUE)
out_dir <- file.path(submission, "data/processed/blast_97id_90cov")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
identity_cutoff <- 97
coverage_cutoff <- 90

ps <- readRDS(file.path(submission, "data/processed/ps_clean_47_rebuilt_from_01_filtered.rds"))
reference <- file.path(submission, "data/reference/Amphibian-skin_bacteria_16S_sequences.fna")
stopifnot(nsamples(ps) == 47L, nzchar(Sys.which("makeblastdb")), nzchar(Sys.which("blastn")))

seqs <- taxa_names(ps); ids <- sprintf("ASV%05d", seq_along(seqs)); names(seqs) <- ids
query <- file.path(out_dir, "rebuilt_47_asv_sequences.fasta")
db <- file.path(out_dir, "woodhams2015")
tsv <- file.path(out_dir, "woodhams_blast_hits_97id_90cov.tsv")
writeXStringSet(DNAStringSet(seqs), query)
stopifnot(system2("makeblastdb", c("-in", shQuote(reference), "-dbtype", "nucl",
  "-parse_seqids", "-out", shQuote(db)), stdout = file.path(out_dir, "makeblastdb.log"),
  stderr = file.path(out_dir, "makeblastdb.log")) == 0)
fmt <- "6 qseqid sseqid pident length mismatch gapopen qlen slen evalue bitscore qcovs"
stopifnot(system2("blastn", c("-query", shQuote(query), "-db", shQuote(db), "-task", "blastn",
  "-strand", "both", "-dust", "no", "-perc_identity", "97", "-qcov_hsp_perc", "90",
  "-max_target_seqs", "100", "-outfmt", shQuote(fmt), "-out", shQuote(tsv)),
  stdout = file.path(out_dir, "blastn.log"), stderr = file.path(out_dir, "blastn.log")) == 0)

cn <- c("ASV_ID", "Reference_ID", "Percent_identity", "Alignment_length", "Mismatches",
  "Gap_opens", "Query_length", "Reference_length", "E_value", "Bit_score", "Query_coverage")
h <- read.delim(tsv, header = FALSE, col.names = cn, stringsAsFactors = FALSE)
y <- tolower(h$Reference_ID)
h$Functional_category <- ifelse(grepl("-inhibitory(_|$)", y), "inhibitory",
  ifelse(grepl("-enhancing(_|$)", y), "enhancing",
  ifelse(grepl("-(ns|non[-_]?significant)(_|$)", y), "non-significant", NA_character_)))
h$High_confidence <- h$Percent_identity >= identity_cutoff & h$Query_coverage >= coverage_cutoff
write.csv(h, file.path(out_dir, "woodhams_blast_all_hits.csv"), row.names = FALSE)
hc <- h[h$High_confidence & !is.na(h$Functional_category), ]
hc <- hc[order(hc$ASV_ID, -hc$Bit_score, -hc$Percent_identity, -hc$Query_coverage, hc$E_value), ]
best <- hc[!duplicated(hc$ASV_ID), ]
write.csv(best, file.path(out_dir, "woodhams_best_high_confidence_hits.csv"), row.names = FALSE)

tax <- as(tax_table(ps), "matrix")
ann <- merge(data.frame(ASV_ID = ids, Sequence = seqs, tax, check.names = FALSE), best,
  by = "ASV_ID", all.x = TRUE, sort = FALSE)
ann$Functional_category[is.na(ann$Functional_category)] <- "unmatched"
ann$Anti_Bd <- ann$Functional_category == "inhibitory"
write.csv(ann, file.path(out_dir, "integrated_asv_taxonomy_function_rebuilt.csv"), row.names = FALSE)
id_to_seq <- setNames(seqs, ids)
anti_seq <- unname(id_to_seq[ann$ASV_ID[ann$Anti_Bd]])
rebuilt_anti <- prune_taxa(anti_seq, ps)
saveRDS(rebuilt_anti, file.path(out_dir, "ps_anti_bd_47_rebuilt_97id_90cov.rds"))
cat("Created anti-Bd phyloseq object with", ntaxa(rebuilt_anti), "ASVs and",
    nsamples(rebuilt_anti), "samples.\n")
