#!/usr/bin/env Rscript

library(betapart, quietly = TRUE)
library(phyloseq, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(magrittr, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(phytools, quietly = TRUE)
library(patchwork, quietly = TRUE)
source("./code/999_utils.R")

# Load phyloseq objects
euk_dataset = readRDS("intermediates/EUK.phyloseq_03_rarefied_phyloseq.rds")
bac_dataset = readRDS("intermediates/BAC.phyloseq_03_rarefied_phyloseq.rds")

# Make plots with custom function in code/999_utils.R
p_nest_meta = plot_nestedness(subset_taxa(euk_dataset, Rank_6 == "Metazoa"), "nonphylo", "jaccard","plot") + ggtitle("Metazoa")+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), plot.title = element_text(hjust = 0.5))

p_nest_fungi = plot_nestedness(subset_taxa(euk_dataset, Rank_5 == "Fungi"), "nonphylo", "jaccard", "plot") + ggtitle("Fungi")+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), plot.title = element_text(hjust = 0.5))


prot_data = subset_taxa(euk_dataset, Domain == "Eukaryota" & (((Rank_6 != "Metazoa") &
                                                                 (Rank_5 != "Fungi") & (Rank_2 != "Chloroplastida") & ( Rank_3 != "Ochrophyta"))| Rank_5 == "Chromulinales") & (Rank_1 != "unclassified_Eukaryota"))

p_nest_prot = plot_nestedness(prot_data, "nonphylo", "jaccard","plot") + ggtitle("Protists")+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), plot.title = element_text(hjust = 0.5))

p_nest_bac = plot_nestedness(bac_dataset, "nonphylo", "jaccard","plot") + ggtitle("Bacteria")+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), plot.title = element_text(hjust = 0.5))


plot_nestedness = p_nest_meta|p_nest_fungi|p_nest_prot|p_nest_bac

ggsave("./results/figures/plot_nestedness.pdf", plot = plot_nestedness, height = 7.5, width = 20)
