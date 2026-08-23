#!/usr/bin/env Rscript

# Rebuild a phyloseq object from the completed quality-filtered FASTQ checkpoint.
# This does not overwrite the deposited my_phyloseq_data.rds.
suppressPackageStartupMessages({library(dada2); library(phyloseq)})

filtered_dir <- Sys.getenv("FILTERED_DIR", "data/raw/01_filtered_fastq")
metadata_file <- Sys.getenv("METADATA_CSV", "data/raw/Acallidryas_DF.csv")
silva_file <- Sys.getenv("SILVA_FASTA", "data/reference/silva_nr99_v138.1_train_set.fa.gz")
rebuild_dir <- Sys.getenv("REBUILD_DIR", "data/processed/rebuilt_from_01_filtered")
final_file <- Sys.getenv(
  "REBUILT_PHYLOSEQ",
  "data/processed/my_phyloseq_data_rebuilt_from_01_filtered.rds"
)
dir.create(rebuild_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(rebuild_dir, "rebuild_log.txt")
log_msg <- function(...) {
  z <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste(..., collapse = " "))
  cat(z, "\n"); cat(z, "\n", file = log_file, append = TRUE)
}

stopifnot(dir.exists(filtered_dir), file.exists(metadata_file), file.exists(silva_file))
fnF <- sort(list.files(filtered_dir, pattern = "_F_filt[.]fastq[.]gz$", full.names = TRUE))
fnR <- sort(list.files(filtered_dir, pattern = "_R_filt[.]fastq[.]gz$", full.names = TRUE))
idF <- sub("_F_filt[.]fastq[.]gz$", "", basename(fnF))
idR <- sub("_R_filt[.]fastq[.]gz$", "", basename(fnR))
stopifnot(length(fnF) > 0L, length(fnF) == length(fnR), identical(idF, idR))
names(fnF) <- names(fnR) <- idF
set.seed(20260813L)
log_msg("RECOVERY START: rebuilding from completed quality-filtered FASTQ checkpoint;")
log_msg("input pairs =", length(fnF), "; existing processed object will not be overwritten.")

derep_file <- file.path(rebuild_dir, "dereplicated_reads.rds")
if (file.exists(derep_file)) {
  derep <- readRDS(derep_file)
  derepF <- derep$derepF; derepR <- derep$derepR
  log_msg("Loaded dereplication checkpoint.")
} else {
  derepF <- derepFastq(fnF, verbose = TRUE)
  derepR <- derepFastq(fnR, verbose = TRUE)
  input_counts <- data.frame(
    SampleID = idF,
    Filtered_R1_reads = vapply(derepF, function(x) sum(getUniques(x)), numeric(1)),
    Filtered_R2_reads = vapply(derepR, function(x) sum(getUniques(x)), numeric(1)),
    stringsAsFactors = FALSE
  )
  write.csv(input_counts, file.path(rebuild_dir, "filtered_fastq_input_counts.csv"), row.names = FALSE)
  keep <- input_counts$Filtered_R1_reads > 0 & input_counts$Filtered_R2_reads > 0
  if (any(!keep)) {
    write.csv(input_counts[!keep, ], file.path(rebuild_dir, "zero_read_pairs_excluded.csv"), row.names = FALSE)
    derepF <- derepF[keep]; derepR <- derepR[keep]; fnF <- fnF[keep]; fnR <- fnR[keep]; idF <- idF[keep]
  }
  saveRDS(list(derepF = derepF, derepR = derepR), derep_file)
  log_msg("Dereplication complete; eligible pairs =", length(derepF))
}

error_file <- file.path(rebuild_dir, "error_models.rds")
if (file.exists(error_file)) {
  errors <- readRDS(error_file); errF <- errors$errF; errR <- errors$errR
  log_msg("Loaded error-model checkpoint.")
} else {
  set.seed(20260813L)
  errF <- learnErrors(fnF, multithread = 4, randomize = TRUE)
  set.seed(20260813L)
  errR <- learnErrors(fnR, multithread = 4, randomize = TRUE)
  saveRDS(list(errF = errF, errR = errR), error_file)
  log_msg("Error learning complete.")
}

dada_file <- file.path(rebuild_dir, "dada_denoised_reads.rds")
if (file.exists(dada_file)) {
  denoised <- readRDS(dada_file); dadaF <- denoised$dadaF; dadaR <- denoised$dadaR
  log_msg("Loaded denoising checkpoint.")
} else {
  dadaF <- dada(derepF, err = errF, pool = FALSE, multithread = 4)
  dadaR <- dada(derepR, err = errR, pool = FALSE, multithread = 4)
  saveRDS(list(dadaF = dadaF, dadaR = dadaR), dada_file)
  log_msg("DADA2 denoising complete.")
}

merge_file <- file.path(rebuild_dir, "merged_pairs.rds")
if (file.exists(merge_file)) {
  mergers <- readRDS(merge_file)
  log_msg("Loaded merged-pair checkpoint.")
} else {
  mergers <- mergePairs(dadaF, derepF, dadaR, derepR, verbose = TRUE)
  saveRDS(mergers, merge_file)
  log_msg("Paired-read merging complete.")
}

seqtab_file <- file.path(rebuild_dir, "sequence_tables.rds")
if (file.exists(seqtab_file)) {
  tabs <- readRDS(seqtab_file); seqtab <- tabs$seqtab; seqtab_nochim <- tabs$seqtab_nochim
  log_msg("Loaded sequence-table checkpoint.")
} else {
  seqtab <- makeSequenceTable(mergers)
  seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = 4, verbose = TRUE)
  saveRDS(list(seqtab = seqtab, seqtab_nochim = seqtab_nochim), seqtab_file)
  log_msg("Sequence table and chimera removal complete; ASVs =", ncol(seqtab_nochim))
}

taxonomy_file <- file.path(rebuild_dir, "silva_taxonomy.rds")
if (file.exists(taxonomy_file)) {
  tax <- readRDS(taxonomy_file)
  log_msg("Loaded SILVA taxonomy checkpoint.")
} else {
  tax <- assignTaxonomy(seqtab_nochim, silva_file, tryRC = TRUE, multithread = 4)
  saveRDS(tax, taxonomy_file)
  log_msg("SILVA 138.1 taxonomy assignment complete.")
}

meta <- read.csv(metadata_file, check.names = FALSE, fileEncoding = "latin1",
                 stringsAsFactors = FALSE)
meta <- meta[!is.na(meta$Sample.ID) & nzchar(meta$Sample.ID) & !duplicated(meta$Sample.ID), ]
missing_meta <- setdiff(rownames(seqtab_nochim), meta$Sample.ID)
stopifnot(length(missing_meta) == 0L)
meta <- meta[match(rownames(seqtab_nochim), meta$Sample.ID), ]
rownames(meta) <- meta$Sample.ID

ps_rebuilt <- phyloseq(
  otu_table(seqtab_nochim, taxa_are_rows = FALSE),
  tax_table(tax),
  sample_data(meta)
)
saveRDS(ps_rebuilt, final_file)

track <- data.frame(
  SampleID = names(derepF),
  Filtered = vapply(derepF, function(x) sum(getUniques(x)), numeric(1)),
  Denoised_F = vapply(dadaF, function(x) sum(getUniques(x)), numeric(1)),
  Denoised_R = vapply(dadaR, function(x) sum(getUniques(x)), numeric(1)),
  Merged = vapply(mergers, function(x) sum(x$abundance), numeric(1)),
  Nonchim = as.numeric(rowSums(seqtab_nochim)[names(derepF)]),
  stringsAsFactors = FALSE
)
write.csv(track, file.path(rebuild_dir, "read_tracking.csv"), row.names = FALSE)
write.csv(data.frame(SampleID = sample_names(ps_rebuilt),
                     Final_depth = as.numeric(sample_sums(ps_rebuilt))),
          file.path(rebuild_dir, "rebuilt_sample_depths.csv"), row.names = FALSE)

old_file <- "data/processed/my_phyloseq_data.rds"
if (file.exists(old_file)) {
  old <- readRDS(old_file)
  comparison <- data.frame(
    Object = c("Existing deposited object", "Rebuilt from 01_filtered checkpoint"),
    Samples = c(nsamples(old), nsamples(ps_rebuilt)),
    ASVs = c(ntaxa(old), ntaxa(ps_rebuilt)),
    Total_reads = c(sum(sample_sums(old)), sum(sample_sums(ps_rebuilt))),
    stringsAsFactors = FALSE
  )
  write.csv(comparison, file.path(rebuild_dir, "comparison_with_existing_object.csv"), row.names = FALSE)
  print(comparison)
}
log_msg("COMPLETE: saved", final_file, "; samples =", nsamples(ps_rebuilt),
        "; ASVs =", ntaxa(ps_rebuilt), "; reads =", sum(sample_sums(ps_rebuilt)))
