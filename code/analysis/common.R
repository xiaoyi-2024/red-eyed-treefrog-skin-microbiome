# Shared configuration for all downstream analyses.
# Run scripts from the submission directory; all paths are package-relative.
suppressPackageStartupMessages(library(phyloseq))

processed_dir <- Sys.getenv("PROCESSED_DIR", "data/processed")
results_root <- "results/analysis"
figures_root <- "figures"
dir.create(results_root, recursive=TRUE, showWarnings=FALSE)
dir.create(figures_root, recursive=TRUE, showWarnings=FALSE)

ps_clean_file <- file.path(processed_dir,"ps_clean_47_rebuilt_from_01_filtered.rds")
ps_anti_file <- file.path(processed_dir,"blast_97id_90cov","ps_anti_bd_47_rebuilt_97id_90cov.rds")
stopifnot(file.exists(ps_clean_file),file.exists(ps_anti_file))
ps_clean <- readRDS(ps_clean_file)
ps_anti_bd <- readRDS(ps_anti_file)
stopifnot(nsamples(ps_clean)==47L,identical(sort(sample_names(ps_clean)),sort(sample_names(ps_anti_bd))))

analysis_seed <- 20260720L
site_levels <- c("Veragua Rainforest","PN Altos de Campana","PN Soberanía","Lost Iguana Resort")
site_colors <- c("Veragua Rainforest"="#3B82A0","PN Altos de Campana"="#5AAE61",
                 "PN Soberanía"="#E6A23C","Lost Iguana Resort"="#C95D63")

short_site <- function(x) {
  z <- iconv(x,from="",to="ASCII//TRANSLIT",sub="")
  ifelse(grepl("^Veragua",z),"Veragua Rainforest",
  ifelse(grepl("Altos de Campana",z),"PN Altos de Campana",
  ifelse(startsWith(z,"PN Sober"),"PN Soberanía",
  ifelse(grepl("Lost Iguana",z),"Lost Iguana Resort",z))))
}
sample_matrix <- function(ps) {
  z<-as(otu_table(ps),"matrix");if(taxa_are_rows(ps))z<-t(z);z
}
site_factor <- function(ps,ids=sample_names(ps)) {
  md<-as(sample_data(ps),"data.frame");factor(short_site(as.character(md[ids,"Locality"])),levels=site_levels)
}
save_plot <- function(p,stem,width,height) {
  ggplot2::ggsave(file.path(figures_root,paste0(stem,".png")),p,width=width,height=height,dpi=400,bg="white")
  # Base PDF avoids an unnecessary XQuartz/Cairo dependency on macOS.
  ggplot2::ggsave(file.path(figures_root,paste0(stem,".pdf")),p,width=width,height=height,device="pdf")
  ggplot2::ggsave(file.path(figures_root,paste0(stem,".svg")),p,width=width,height=height,device=svglite::svglite)
}
