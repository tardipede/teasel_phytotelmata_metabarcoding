#!/usr/bin/env Rscript
library(phyloseq, , quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(patchwork, quietly = TRUE)
library(writexl, quietly = TRUE)
library(magrittr, quietly = TRUE)
options(tidyverse.quiet = TRUE)
library(tidyverse, quietly = TRUE)
source("./code/999_utils.R")

# Load phyloseq objects
euk_dataset = readRDS("intermediates/EUK.phyloseq_03_rarefied_phyloseq.rds")
bac_dataset = readRDS("intermediates/BAC.phyloseq_03_rarefied_phyloseq.rds")


# Calculate richness
meta_richness = estimate_richness(subset_taxa(euk_dataset, Rank_6 == "Metazoa"), measures = c("Observed", "Shannon", "FaithPD")) %>% 
  mutate(sample = rownames(.)) %>%
  merge(.,data.frame(sample_data(subset_taxa(euk_dataset, Rank_6 == "Metazoa"))), by = "sample")

fungi_richness = estimate_richness(subset_taxa(euk_dataset, Rank_5 == "Fungi"), measures = c("Observed", "Shannon", "FaithPD")) %>% 
  mutate(sample = rownames(.)) %>%
  merge(.,data.frame(sample_data(subset_taxa(euk_dataset, Rank_5 == "Fungi"))), by = "sample")

algae_richness = estimate_richness(subset_taxa(euk_dataset, (Rank_2 == "Chloroplastida" | Rank_3 == "Ochrophyta") & Rank_5 != "Chromulinales"), measures = c("Observed", "Shannon", "FaithPD")) %>% 
  mutate(sample = rownames(.)) %>%
  merge(.,data.frame(sample_data(subset_taxa(euk_dataset, (Rank_2 == "Chloroplastida" | Rank_3 == "Ochrophyta") & Rank_5 != "Chromulinales"))), by = "sample")

prot_data = subset_taxa(euk_dataset, Domain == "Eukaryota" & (((Rank_6 != "Metazoa") &
                          (Rank_5 != "Fungi") & (Rank_2 != "Chloroplastida") & ( Rank_3 != "Ochrophyta"))| Rank_5 == "Chromulinales") & (Rank_1 != "unclassified_Eukaryota"))
prot_richness = estimate_richness(prot_data, measures = c("Observed", "Shannon", "FaithPD")) %>% 
  mutate(sample = rownames(.)) %>%
  merge(.,data.frame(sample_data(prot_data)), by = "sample")

bac_richness = estimate_richness(bac_dataset, measures = c("Observed", "Shannon", "FaithPD")) %>% 
  mutate(sample = rownames(.)) %>%
  merge(.,data.frame(sample_data(bac_dataset)), by = "sample")

meta_richness$taxa = rep("Metazoa", nrow(meta_richness))
fungi_richness$taxa = rep("Fungi", nrow(fungi_richness))
algae_richness$taxa = rep("Algae", nrow(algae_richness))
prot_richness$taxa = rep("Protists", nrow(prot_richness))
bac_richness$taxa = rep("Bacteria", nrow(prot_richness))
richness_all = rbind(meta_richness, fungi_richness, algae_richness, prot_richness,bac_richness)
richness_all$Evenness = richness_all$Shannon/log(richness_all$Observed)

write.table(richness_all, "./intermediates/richness_tab.txt", sep = "\t",  quote = F, row.names = F)

### Save richness summary table

richness_summary =  aggregate(richness_all[,c("Observed","Evenness","FaithPD")], by = list(richness_all$site, richness_all$level, richness_all$taxa), FUN = function(x){round(c(mean(x), min(x), max(x)),2)})
richness_summary = do.call(data.frame, richness_summary) 

colnames(richness_summary) = c("Site","Level","taxa","Obs-mean","Obs-min","Obs-max","Eve-mean","Eve-min","Eve-max","PD-mean","PD-min","PD-mad")

writexl::write_xlsx(richness_summary, "./results/tables/richness_summary.xlsx")