library(magrittr)
library(dada2)  
library(ShortRead)  
library(Biostrings)  
library(future)
library(dplyr)
library(data.table)
library(seqinr)
library(tibble)
library(DECIPHER)
library(lulu)
source("util.R", encoding = 'UTF-8')

# General options

RUN_NAME <- "DPT_SSU_2024"

FWD <- "GCTTGTCTCAAAGATTAAGCC"    # Primers sequences
REV <- "GCCTGCTGCCTTCCTTGGA"

#DATABASE_PATH <- "pr2_version_5.0.0_SSU.decipher.trained.rds"
DATABASE_PATH <- "SILVA_SSU_r138.rds"

#### Get the samples organized ####

#define paths
path <- "."
seq_path <- file.path(path, "sequences")
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
  trim_R1 = file.path(trim_path,
                      paste(runcode, sample, "R1_trim.fastq.gz", sep = "_")),
  trim_R2 = file.path(trim_path,
                      paste(runcode, sample, "R2_trim.fastq.gz", sep = "_")),
  filt_R1 = file.path(filt_path,
                      paste(runcode, sample, "R1_filt.fastq.gz", sep = "_")),
  filt_R2 = file.path(filt_path,
                      paste(runcode, sample, "R2_filt.fastq.gz", sep = "_"))
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

CUTADAPT_FLAGS<-paste0("conda run -n trim_galore cutadapt -a ^",FWD,"...",REV.RC," -A ^",REV,"...",FWD.RC," --discard-untrimmed -o")

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
seqtab.nochim.filtered_length <-  seqtab.nochim[ ,nchar(colnames(seqtab.nochim)) >= 350 & nchar(colnames(seqtab.nochim)) <= 450] 
seqtab.nochim.filtered_length.filtered_min_abundance<-seqtab.nochim.filtered_length[ , apply(seqtab.nochim.filtered_length,2,sum) >= 10]   #filtering reads with counts < 10 
# sorting by abundance
seqtab.nochim.filtered_length.filtered_min_abundance<-seqtab.nochim.filtered_length.filtered_min_abundance[, order(-colSums(seqtab.nochim.filtered_length.filtered_min_abundance))]  



# ASV curation with lulu
seqtable = seqtab.nochim.filtered_length.filtered_min_abundance
asv_sequences <- colnames(seqtable)
sample_names <- rownames(seqtable)
dna = Biostrings::DNAStringSet(asv_sequences)
names(dna) = paste0("OTU",1:length(asv_sequences))
dir.create("blast")
writeXStringSet(dna, './blast/seqs.fasta')
#First produce a blastdatabase with the OTUs
system("makeblastdb -in ./blast/seqs.fasta -parse_seqids -dbtype nucl -out ./blast/seqs")
# Then blast the OTUs against the database
system("blastn -db ./blast/seqs -outfmt 6  -out ./blast/match_list.txt -qcov_hsp_perc 80 -perc_identity 90 -query ./blast/seqs.fasta")

# Load blast results
blast_res = read.table("./blast/match_list.txt", header = F, sep = "\t")
# Remove blast folder
unlink("./blast", recursive = TRUE)
# Format blast results
blast_res = blast_res[,1:3]
# Formatu otutable
otutable_lulu = t(seqtable)
rownames(otutable_lulu) = names(dna)
#otutable_lulu = cbind(data.frame(OTUid = names(aln)), otutable_lulu)
otutable_lulu = as.data.frame(otutable_lulu)
# Run lulu
curated_result <- lulu(otutable_lulu, blast_res, minimum_match = 99, minimum_relative_cooccurence = 1, minimum_ratio = 10)
# Extract curated ASV table
curated_result = curated_result$curated_table
dna = dna[names(dna) %in% rownames(curated_result)]
dna = dna[rownames(curated_result)]


# Clustering
#seqtab_clust = cluster_seqs(seqtab.nochim.filtered_length.filtered_min_abundance, nproc = 12, cutoff = 0.005, align = T)
#otu_seqtab = seqtab_clust[[1]]

# Saving OTU sequences to fasta file
#seqs_otu<-otu_seqtab$seq
#names_otu<-otu_seqtab$otu
#write.fasta(as.list(seqs_otu),names_otu,paste0(RUN_NAME,"_OTUS.fasta"))

#dna <- Biostrings::DNAStringSet(seqs_otu)


#zOTU classification using SILVA 18s
trainingSet =   readRDS(DATABASE_PATH) #18s database loading
ids <- DECIPHER::IdTaxa(dna,   trainingSet,   strand="both", processors = 12)


#exporting classification results
output <- sapply(ids,
                 function (id) {
                   paste(id$taxon,
                         " (",
                         round(id$confidence, digits=1),
                         "%)",
                         sep="",
                         collapse="; ")})

otu_seqtab = curated_result

otutable_classified<-cbind(data.frame(output),data.frame(OTU = paste0("OTU",1:length(output))),otu_seqtab)
otutable_classified<-otutable_classified %>%select(output, everything()) 
write.table(otutable_classified,file=paste0(RUN_NAME,"_classified_otus.csv"), quote = F)
writexl::write_xlsx(otutable_classified,path=paste0(RUN_NAME,"_classified_otus.xlsx"))


# Align OTUs and save alignment
aln = DECIPHER::AlignSeqs(dna, processors = 12)
names(aln) = otutable_classified$OTU
Biostrings::writeXStringSet(aln, paste0(RUN_NAME,"_aligned_otus.fas"))


# Make OTUs phylogenetic tree with Fastree
system(paste0("FastTree -gtr -nt ",paste0(RUN_NAME,"_aligned_otus.fas")," > ",RUN_NAME,"_otus_phylotree.tre"))