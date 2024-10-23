#!/usr/bin/env Rscript

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
# Check if arguments are provided
if (length(args) == 0) {
  stop("No arguments provided, provide run name [arg1], data folder [arg2] and marker (18S or 16S) [arg3]")
}


## Libraryes
library(tidyverse)
library(magrittr)
library(dada2)
library(ShortRead)  
library(Biostrings)  
library(DECIPHER)

# Assign provided arguments
RUN_NAME <- args[1]
DATA_FOLDER <- args[2]
MARKER <- args[3]

#RUN_NAME <- "SSU_test"
#DATA_FOLDER <- "18S"
#MARKER <- "18S"


# Create variables used in the script
DATABASE_PATH <- "SILVA_SSU_r138.rds"

FWD <- case_when(
  MARKER  == "18S" ~ "GCTTGTCTCAAAGATTAAGCC",
  MARKER  == "16S" ~ "GTGCCAGCMGCCGCGGTAA",
  TRUE ~ NA
)

REV <- case_when(
  MARKER  == "18S" ~ "GCCTGCTGCCTTCCTTGGA",
  MARKER  == "16S" ~ "GACTACNVGGGTWTCTAAT",
  TRUE ~ NA
)

min_amplicon_size <- case_when(
  MARKER  == "18S" ~ 350,
  MARKER  == "16S" ~ 250,
  TRUE ~ NA
)

max_amplicon_size <- case_when(
  MARKER  == "18S" ~ 450,
  MARKER  == "16S" ~ 330,
  TRUE ~ NA
)


#### Get the samples organized ####

#define paths
path <- "."
seq_path <- file.path(path, "data", DATA_FOLDER)
raw_path <- file.path(seq_path, "01_raw")
trim_path <- file.path(seq_path, "02_filter")
filt_path <- file.path(seq_path, "03_trim")
backup_path <- file.path(path, "backups")

# create paths if missing
if (!dir.exists(backup_path)) dir.create(backup_path, recursive = TRUE)
if (!dir.exists(filt_path)) dir.create(filt_path, recursive = TRUE)
if (!dir.exists(trim_path)) dir.create(trim_path, recursive = TRUE)

# find files
sample_table <- tibble::tibble(
  fastq_R1 = sort(list.files(raw_path, ".*R1(_001)?.fastq.gz",full.names = T)),
  fastq_R2 = sort(list.files(raw_path, ".*R2(_001)?.fastq.gz",full.names = T))
)  %>%  dplyr::mutate(runcode = stringr::str_extract(fastq_R1,"(ID[0-9]{4})"),
                      sample = stringr::str_extract(fastq_R1,"((Blank|DFP)-.{3})"))


# generate filenames for trimmed and filtered reads
sample_table <- sample_table %>%  dplyr::mutate(
  trim_R1 = file.path(trim_path, paste(runcode, sample, "R1_trim.fastq.gz", sep = "_")),
  trim_R2 = file.path(trim_path, paste(runcode, sample, "R2_trim.fastq.gz", sep = "_")),
  filt_R1 = file.path(filt_path, paste(runcode, sample, "R1_filt.fastq.gz", sep = "_")),
  filt_R2 = file.path(filt_path, paste(runcode, sample, "R2_filt.fastq.gz", sep = "_"))
)

# Check that there are no duplicates
assertthat::assert_that(
  !any(duplicated(sample_table$fastq_R1)),
  !any(duplicated(sample_table$fastq_R2)),
  !any(duplicated(sample_table$trim_R1)),
  !any(duplicated(sample_table$trim_R2)),
  !any(duplicated(sample_table$filt_R1)),
  !any(duplicated(sample_table$filt_R2))
)

cat("Found", nrow(sample_table), "samples.\n")


# Setting flags for Cutadapt and trimming primers, filtering reads without primers
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

CUTADAPT_FLAGS<-paste0("cutadapt -a ^",FWD,"...",REV.RC," -A ^",REV,"...",FWD.RC," --discard-untrimmed -o")

for(i in seq_along(sample_table$fastq_R1)) {
  system(paste(CUTADAPT_FLAGS, sample_table$trim_R1[i],"-p",sample_table$trim_R2[i],sample_table$fastq_R1[i], sample_table$fastq_R2[i]))
}


# Quality-filtering of reads
out <- filterAndTrim(sample_table$trim_R1, sample_table$filt_R1,
                     sample_table$trim_R2, sample_table$filt_R2, 
                     maxN=0, truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=T,minLen=190 ) 


# Error-rate estimation, merging and chimera filtering done separately for each run
error_results = list()
backups = list()
for(i in 1:length(unique(sample_table$runcode))){
  
  sample_table_temp = subset(sample_table, runcode == unique(sample_table$runcode)[i])
  
  ## calculating errors - on a subset of data to save time (https://benjjneb.github.io/dada2/bigdata.html)
  set.seed(100)
  errF <- learnErrors(sample_table_temp$filt_R1, multithread=TRUE, nbases = 1e+8, randomize=TRUE)
  errR <- learnErrors(sample_table_temp$filt_R2, multithread=TRUE, nbases = 1e+8, randomize=TRUE)
  # Sample inference
  dadaFs <- dada(sample_table_temp$filt_R1, err=errF, multithread=T)
  dadaRs <- dada(sample_table_temp$filt_R2, err=errR, multithread=T)
  # Merging paired reads
  mergers <- mergePairs(dadaFs, sample_table_temp$filt_R1, dadaRs, sample_table_temp$filt_R2,   minOverlap = 8)
  #Creating sequence table
  seqtab <- makeSequenceTable(mergers)
  # Removing chimeras
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
  error_results[[i]] = seqtab.nochim
  backups[[i]] = list(errF,errR,dadaFs,dadaRs, mergers, seqtab)
  
  rm(sample_table_temp)
}

# Merge the results from different runs together
if(length(error_results) == 1){seqtab.nochim = error_results[[1]]}
if(length(error_results) != 1){
  seqtab.nochim = dada2::mergeSequenceTables(tables = error_results)
}

# Save the sequences table
saveRDS(seqtab.nochim, file.path(backup_path,"savepoint1_seqtable.rds"))
seqtab.nochim = readRDS(file.path(backup_path,"savepoint1_seqtable.rds"))

# Abundance filtering, length filtering
# This step removes lots of erronneous sequences and saves a lot of RAM to the next steps
seqtab.nochim.filtered_length <-  seqtab.nochim[ ,nchar(colnames(seqtab.nochim)) >= min_amplicon_size & nchar(colnames(seqtab.nochim)) <= max_amplicon_size] 
seqtab.nochim.filtered_length.filtered_min_abundance<-seqtab.nochim.filtered_length[ , apply(seqtab.nochim.filtered_length,2,sum) >= 10]   #filtering reads with counts < 10 
# sorting by abundance
seqtab.nochim.filtered_length.filtered_min_abundance<-seqtab.nochim.filtered_length.filtered_min_abundance[, order(-colSums(seqtab.nochim.filtered_length.filtered_min_abundance))]  

