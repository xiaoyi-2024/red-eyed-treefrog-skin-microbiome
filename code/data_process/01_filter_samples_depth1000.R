#!/usr/bin/env Rscript

# Construct the manuscript analysis set from the DADA2-processed phyloseq object.
suppressPackageStartupMessages(library(phyloseq))
out_dir <- Sys.getenv("PROCESSED_DIR", "data/processed")
infile <- file.path(out_dir, "my_phyloseq_data.rds")
stopifnot(file.exists(infile))
ps <- readRDS(infile)

md <- as(sample_data(ps), "data.frame")
depth <- sample_sums(ps)
is_skin <- as.character(md$Sample.type) == "skin swab"
qc <- data.frame(Sample.ID=sample_names(ps), Sample.type=as.character(md$Sample.type),
  Final_depth=as.numeric(depth), Is_skin_sample=is_skin,
  Pass_depth_1000=depth>=1000, stringsAsFactors=FALSE)
qc$Included_final <- qc$Is_skin_sample & qc$Pass_depth_1000
qc$Exclusion_reason <- ifelse(!qc$Is_skin_sample,"control",
  ifelse(!qc$Pass_depth_1000,"depth <1000","included"))
write.csv(qc,file.path(out_dir,"sample_depth_qc.csv"),row.names=FALSE)

ps47 <- prune_samples(qc$Sample.ID[qc$Included_final],ps)
ps47 <- prune_taxa(taxa_sums(ps47)>0,ps47)
stopifnot(nsamples(ps47)==47L)
saveRDS(ps47,file.path(out_dir,"ps_clean_47.rds"))

otu <- as(otu_table(ps47),"matrix"); if(taxa_are_rows(ps47)) otu <- t(otu)
write.csv(data.frame(Sample.ID=rownames(otu),otu,check.names=FALSE),
  file.path(out_dir,"asv_feature_table_47.csv"),row.names=FALSE)
write.csv(data.frame(Sample.ID=sample_names(ps47),as(sample_data(ps47),"data.frame"),
  check.names=FALSE),file.path(out_dir,"sample_metadata_47.csv"),row.names=FALSE)
