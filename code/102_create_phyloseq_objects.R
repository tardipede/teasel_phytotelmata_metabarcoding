#!/usr/bin/env Rscript

library(phyloseq, quietly = TRUE)
library(Biostrings, quietly = TRUE)
library(ape, quietly = TRUE)
library(magrittr, quietly = TRUE)
options(tidyverse.quiet = TRUE)
library(tidyverse, quietly = TRUE)
source("./code/999_utils.R")


### Eukaryotes dataset

PREFIX = "EUK"
OUT_NAME = "EUK.phyloseq"
MARKER = "18S"
SAMPLES_DATA = "./samples_data/DFP_samples_data.txt"

if(MARKER  == "18S"){taxonomy_ranks = c("Root","Domain",paste0("Rank_",1:17))}
if(MARKER  == "16S"){taxonomy_ranks = c("Root","Domain","Phylum","Class","Order","Family","Genus")}

paths = c(list.files("./intermediates", pattern = PREFIX, full.names = T), list.files("./samples_data", full.names = T))

euk_phyloseq_dataset = make_phyloseq(paths = paths,taxonomy_ranks = taxonomy_ranks) %>%
                 saveRDS_pipe(filename = paste0("./intermediates/",OUT_NAME,"_01_all_reads.rds")) %>%
                 subtract_blank(blank_name = "DFP.000") %>%
                 prune_samples(sample_sums(.)!=0, .) %>%
                 phyloseq::subset_taxa(Domain != "unclassified_Root") %>% # Remove unclassified seqs
                 phyloseq::subset_taxa(Rank_6 != "Embryophyta") %>% # Remove embriophytes
                 saveRDS_pipe(filename = paste0("./intermediates/",OUT_NAME,"_02_noblank_phyloseq.rds")) %>%
                 prune_samples(!sample_names(.) %in% c("DFP.025","DFP.026"), .) %>%  # samples removed as they have too few reads
                 rarefy_even_depth(rngseed = 123456789) %>%
                 saveRDS_pipe(filename = paste0("./intermediates/",OUT_NAME,"_03_rarefied_phyloseq.rds"))




### Bacteria dataset

PREFIX = "BAC"
OUT_NAME = "BAC.phyloseq"
MARKER = "16S"
SAMPLES_DATA = "./samples_data/DFP_samples_data.txt"

if(MARKER  == "18S"){taxonomy_ranks = c("Root","Domain",paste0("Rank_",1:17))}
if(MARKER  == "16S"){taxonomy_ranks = c("Root","Domain","Phylum","Class","Order","Family","Genus")}

paths = c(list.files("./intermediates", pattern = PREFIX, full.names = T), list.files("./samples_data", full.names = T))

bac_phyloseq_dataset = make_phyloseq(paths = paths,taxonomy_ranks = taxonomy_ranks) %>%
                 saveRDS_pipe(filename = paste0("./intermediates/",OUT_NAME,"_01_all_reads.rds")) %>%
                 subtract_blank(blank_name = "DFP.000") %>%
                 prune_samples(sample_sums(.)!=0, .) %>%
                 phyloseq::subset_taxa(Domain == "Bacteria") %>% # Keep only Bacteria
                 phyloseq::subset_taxa(Order != "Chloroplast") %>% # Remove chloroplasts
                 phyloseq::subset_taxa(Family != "Mitochondria") %>% # Remove mitochondria
                 saveRDS_pipe(filename = paste0("./intermediates/",OUT_NAME,"_02_noblank_phyloseq.rds")) %>%
                 prune_samples(!sample_names(.) %in% c("DFP.025","DFP.026"), .) %>%   # Removed as those samples 18S dataset have not enought reads
                 rarefy_even_depth(rngseed = 123456789) %>%
                 saveRDS_pipe(filename = paste0("./intermediates/",OUT_NAME,"_03_rarefied_phyloseq.rds"))
