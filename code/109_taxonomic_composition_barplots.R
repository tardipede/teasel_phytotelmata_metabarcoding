#!/usr/bin/env Rscript

library(phyloseq, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(tidyr, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(magrittr, quietly = TRUE)
library(patchwork, quietly = TRUE)
library(Biostrings, quietly = TRUE)
library(taxonomizr, quietly = TRUE)
source("./code/999_utils.R")

# Load phyloseq objects
euk_dataset <- readRDS("intermediates/EUK.phyloseq_03_rarefied_phyloseq.rds")
bac_dataset <- readRDS("intermediates/BAC.phyloseq_03_rarefied_phyloseq.rds")


##### Plots EUK groups
metazoa <- colSums(otu_table(subset_taxa(euk_dataset, 
                                         Rank_6 == "Metazoa")))
fungi <- colSums(otu_table(subset_taxa(euk_dataset, Rank_5 == "Fungi")))
algae <- colSums(otu_table(subset_taxa(euk_dataset, 
                                       (Rank_2 == "Chloroplastida" | 
                                          Rank_3 == "Ochrophyta") & 
                                         Rank_5 != "Chromulinales")))

prot_data <- subset_taxa(euk_dataset, Domain == "Eukaryota" & 
                           (((Rank_6 != "Metazoa") & (Rank_5 != "Fungi") & 
                               (Rank_2 != "Chloroplastida") & (Rank_3 != "Ochrophyta")) | 
                              Rank_5 == "Chromulinales") & (Rank_1 != "unclassified_Eukaryota"))

protists <- colSums(otu_table(prot_data))

plot_grups_euk <- data.frame(
  sample = names(metazoa),
  metazoa = metazoa,
  fungi = fungi,
  algae = algae,
  protists = protists
) %>%
  pivot_longer(2:5) %>%
  ggplot() +
  theme_bw() +
  geom_bar(aes(x = sample, 
               y = value, 
               fill = name), 
           stat = "identity", 
           position = "fill",
           col = "black", 
           alpha = 0.9) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


##### Plot BAC groups
bac_tax <- data.frame(tax_table(bac_dataset))
bac_tax <- bac_tax$Phylum

bac_otus <- data.frame(otu_table(bac_dataset))
bac_otus <- cbind(bac_tax, bac_otus)

bac_otus <- pivot_longer(bac_otus, 2:ncol(bac_otus))
bac_otus <- aggregate(bac_otus$value, 
                      by = list(bac_otus$bac_tax, bac_otus$name), 
                      FUN = sum)
colnames(bac_otus) <- c("taxa", "sample", "reads")

plot_grups_bac <- ggplot(bac_otus) +
  theme_bw() +
  geom_bar(aes(x = sample, 
               y = reads, 
               fill = taxa), 
           stat = "identity", 
           position = "fill", 
           col = "black", 
           alpha = 0.9) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

##### Plot Metazoa with updated classification
temp_dataset <- subset_taxa(euk_dataset, Rank_6 == "Metazoa")

# Save ASV sequences to update the taxonomy with Blast
dir.create("./intermediates/metazoa_updated_taxonomy")
Biostrings::writeXStringSet(DECIPHER::RemoveGaps(refseq(temp_dataset), 
                                                 removeGaps = "all"), 
                            filepath = "./intermediates/metazoa_updated_taxonomy/otus.fasta")

# Run blastn with a remote database, get the subject taxids
system("blastn -db nt -remote  -query ./intermediates/metazoa_updated_taxonomy/otus.fasta -out ./intermediates/metazoa_updated_taxonomy/blast_results.txt -max_target_seqs 50 -perc_identity 80  -outfmt \"6 qseqid sseqid pident qlen length evalue staxids \" ")

# Download data from NCBI and set up SQLite database
taxonomizr::prepareDatabase(sqlFile = "./databases/nameNode.sqlite", 
                            tmpDir = "./databases/", 
                            getAccessions = FALSE)

# Load blastn results
blast_results <- read.table("./intermediates/metazoa_updated_taxonomy/blast_results.txt", 
                            header = FALSE, 
                            sep = "\t")
colnames(blast_results) <- c("qseqid", "sseqid", "pident", "qlen", "length", "evalue", "staxids")

# Calculate query coverage
blast_results$qseq_coverage <- (blast_results$qlen / blast_results$length) * 100
blast_results <- subset(blast_results, 
                        qseq_coverage >= 95) # Keep only matches with >95% query coverage

# Split blast results by query ID
blast_results_list <- split(blast_results, 
                            blast_results$qseqid)


# Keep only best matches and (best matches - 2% similarity)
blast_results_list = lapply(blast_results_list, 
  FUN = function(x){
     x = x[order(x$pident, decreasing = T),]  
    tresh = max(x$pident)-2                # Keep only matches not lower than best match -2%
     x = subset(x, pident > tresh)
     return(x)}
  )

# Functions to classify

# Classify based on consensus taxonomy of retained matched (keep only a  classification when >90% of matches has the same clasification at one level)
classify_otu <- function(x) {
  taxaId <- x$staxids
  tax <- getTaxonomy(taxaId, sqlFile = "./databases/nameNode.sqlite")
  return(cbind(x[, c("qseqid", "pident")], tax))
}

blast_results_tax <- lapply(blast_results_list, FUN = classify_otu)

classify_base <- function(x) {
  temp <- table(x)
  x <- x[!is.na(x)]
  if (length(temp) == 0) {
    tax_out <- NA
    tax_score <- 0
  }
  if (length(temp) == 1) {
    tax_out <- names(temp)
    tax_score <- 100
  }
  if (length(temp) > 1) {
    tax_out <- names(which.max(temp))
    tax_score <- round((max(temp) / sum(temp)) * 100, 0)
  }
  return(c(tax_out, tax_score))
}

classify_table <- function(x) {
  otu_num <- x[1, 1]
  x <- x[, 3:ncol(x)]
  tax <- apply(x, MARGIN = 2, FUN = classify_base)
  tax <- data.frame(tax)
  tax[1, ][as.numeric(tax[2, ]) < 90] <- ""
  tax[2, ][as.numeric(tax[2, ]) < 90] <- ""
  tax <- c(otu_num, unlist((tax[1, ])), unlist((tax[2, ])))
  return(tax)
}

tax_results = lapply(blast_results_tax, FUN = classify_table)
tax_table = do.call(rbind,tax_results)
tax_table = tax_table[,c("phylum","class","order","family")]
class(tax_table) = "matrix"

new_metazoa_tax = tax_table[match(rownames(tax_table), rownames(tax_table(temp_dataset))),]
new_metazoa_tax = tax_table(new_metazoa_tax)
tax_table(temp_dataset) = new_metazoa_tax
temp_dataset = subset_taxa(temp_dataset, order !="")

temp_dataset =transform_sample_counts(temp_dataset, function(x) x/sum(x)*100)

temp_tax = cbind(rownames(tax_table(temp_dataset)),data.frame(tax_table(temp_dataset)))
colnames(temp_tax)[1] = "ASV"

temp_tax = apply(temp_tax, MARGIN = 1, FUN = function(x){paste0(x[3],";",x[4])})

data_reads = data.frame(otu_table(temp_dataset))
data_reads = cbind(data_reads, temp_tax) %>% 
  pivot_longer(1:(ncol(.)-1))%>%
  group_by(temp_tax, name)


data_reads= aggregate(data_reads$value, 
  by = list(data_reads$temp_tax, data_reads$name), 
  FUN = sum)

colnames(data_reads) = c("temp_tax","name","value")

plot_meta = ggplot(data_reads)+
  theme_bw()+
  geom_bar(aes(x = name, 
  y = value, 
  fill = temp_tax), 
  stat = "identity", 
  position = "stack", 
  col = "black", 
  alpha = 0.9)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

### Plot composition
plot_all = plot_grups_bac/plot_grups_euk/plot_meta
ggsave("./results/figures/figure_barplot_all.pdf", plot = plot_all, height = 15, width = 15)

## Insects prevalence table
temp_tax = cbind(rownames(tax_table(temp_dataset)),data.frame(tax_table(temp_dataset)))
temp_tax = apply(temp_tax, MARGIN = 1, FUN = function(x){paste0(x[3],";",x[4],";",x[5])})

data_reads = data.frame(otu_table(temp_dataset))
data_reads = cbind(data_reads, temp_tax) %>% 
  pivot_longer(1:(ncol(.)-1))%>%
  group_by(temp_tax, name)

data_reads= aggregate(data_reads$value, 
  by = list(data_reads$temp_tax, 
    data_reads$name), 
    FUN = sum)
colnames(data_reads) = c("temp_tax","name","value")

data_to_subset = aggregate(data_reads$value, 
  by = list(data_reads$temp_tax), 
  FUN = max)
data_to_subset = subset(data_to_subset, x > 0)

data_reads = subset(data_reads, temp_tax %in% data_to_subset$Group.1)

data_reads = data_reads[grepl("Insecta", data_reads$temp_tax, fixed = T),]

data_reads_summary = data_reads
data_reads_summary$value[data_reads_summary$value != 0] = 1
prevalence = aggregate(data_reads_summary$value, by = list(data_reads_summary$temp_tax), FUN = mean)
prevalence$x  = round(prevalence$x*100,2)
prevalence = prevalence[order(prevalence$x, decreasing = T),]

writexl::write_xlsx(prevalence, "./results/tables/insect_prevalence.xlsx")