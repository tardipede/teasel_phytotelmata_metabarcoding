#!/usr/bin/env Rscript

library(phyloseq, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(tidyr, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(magrittr, quietly = TRUE)
library(patchwork, quietly = TRUE)
library(Biostrings, quietly = TRUE)
source("./code/999_utils.R")

# Load phyloseq objects
euk_dataset = readRDS("intermediates/EUK.phyloseq_03_rarefied_phyloseq.rds")
bac_dataset = readRDS("intermediates/BAC.phyloseq_03_rarefied_phyloseq.rds")


##### Plots EUK groups
metazoa = colSums(otu_table(subset_taxa(euk_dataset, Rank_6 == "Metazoa")))
fungi = colSums(otu_table(subset_taxa(euk_dataset, Rank_5 == "Fungi")))
algae = colSums(otu_table(subset_taxa(euk_dataset, (Rank_2 == "Chloroplastida" | Rank_3 == "Ochrophyta") & Rank_5 != "Chromulinales")))

prot_data = subset_taxa(euk_dataset, Domain == "Eukaryota" & (((Rank_6 != "Metazoa") &
                                                                 (Rank_5 != "Fungi") & (Rank_2 != "Chloroplastida") & ( Rank_3 != "Ochrophyta"))| Rank_5 == "Chromulinales") & (Rank_1 != "unclassified_Eukaryota"))
protists = colSums(otu_table(prot_data))

plot_grups_euk = data.frame(sample = names(metazoa),
                            metazoa = metazoa,
                            fungi = fungi,
                            algae = algae,
                            protists = protists ) %>%
      pivot_longer(2:5) %>%
  ggplot()+
  theme_bw()+
  geom_bar(aes(x = sample, y = value, fill = name), stat = "identity", position = "fill", col = "black", alpha = 0.9)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


##### Plot BAC groups
bac_tax = data.frame(tax_table(bac_dataset))
bac_tax = bac_tax$Phylum

bac_otus = data.frame(otu_table(bac_dataset))
bac_otus = cbind(bac_tax, bac_otus)

bac_otus = pivot_longer(bac_otus, 2:ncol(bac_otus))
bac_otus = aggregate(bac_otus$value, by = list(bac_otus$bac_tax, bac_otus$name), FUN = sum)
colnames(bac_otus) = c("taxa","sample","reads")

plot_grups_bac = ggplot(bac_otus)+
  theme_bw()+
  geom_bar(aes(x = sample, y = reads, fill = taxa), stat = "identity", position = "fill", col = "black", alpha = 0.9)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

##### Plot Metazoa with updated classification
temp_dataset = subset_taxa(euk_dataset, Rank_6 == "Metazoa")

# Save ASV sequences to update the taxonomy with Blast
dir.create("./intermediates/metazoa_updated_taxonomy")
Biostrings::writeXStringSet(DECIPHER::RemoveGaps(refseq(temp_dataset),removeGaps = "all"), filepath = "./intermediates/metazoa_updated_taxonomy/otus.fasta")

system(
"blastn -db nt -remote  -query ./intermediates/metazoa_updated_taxonomy/otus.fasta -out ./intermediates/metazoa_updated_taxonomy/blast_results.txt -max_target_seqs 50 -perc_identity 80 -outfmt 6"
)
