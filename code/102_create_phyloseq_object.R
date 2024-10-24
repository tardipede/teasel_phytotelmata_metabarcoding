#!/usr/bin/env Rscript

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
# Check if arguments are provided
if (length(args) != 4) {
  stop("Wrong number of arguments provided, provide FILES PREFIX  [arg1], OUT_NAME [arg2], MARKER [arg3](18S or 16S) and SAMPLES_DATA path [arg4]")
}

library(phyloseq, quietly = TRUE)
library(Biostrings, quietly = TRUE)
library(ape, quietly = TRUE)
library(magrittr, quietly = TRUE)
options(tidyverse.quiet = TRUE)
library(tidyverse, quietly = TRUE)

#PREFIX = args[1]
#OUT_NAME = args[2]
#MARKER = args[3]
#SAMPLES_DATA = args[3]

PREFIX = "EUK"
OUT_NAME = "EUK.phyloseq.rds"
MARKER = "18S"
SAMPLES_DATA = "./samples_data/DFP_samples_data.txt"

if(MARKER  == "18S"){taxonomy_ranks = c("Root","Domain",paste0("Rank_",1:17))}
if(MARKER  == "16S"){taxonomy_ranks = c("Root","Domain","Phylum","Class","Order","Family","Genus")}


paths = c(list.files("./intermediates/", pattern = "EUK", full.names = T), list.files("./samples_data/", full.names = T))
path2 = list.files("./samples_data/", full.names = T)

make_phyloseq = function(paths,taxonomy_ranks){
  require(ape)
  require(Biostrings)
  require(phyloseq)
  require(magrittr)
  require(plyr)

  # Get the different files
  otu_file = paths[grepl(".csv",paths)]
  tree_file = paths[grepl(".tre",paths)]
  seqs_file = paths[grepl(".fas",paths)]
  metadata_file = paths[grepl(".txt",paths)]


  # Load files
  otu_tab = read.table(otu_file, sep = "\t", header = T);  rownames(otu_tab) = otu_tab$ASV
  sample_dat = data.frame(readxl::read_xlsx(metadata_file)); rownames(sample_dat) = sample_dat$sample
  tree = ape::read.tree(tree_file)
  seqs = Biostrings::readDNAStringSet(seqs_file)
