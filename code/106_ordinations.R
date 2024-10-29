#!/usr/bin/env Rscript

library(phyloseq, quietly = TRUE)
library(vegan, quietly = TRUE)
library(betapart, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(patchwork, quietly = TRUE)

# Load phyloseq objects
euk_dataset <- readRDS("intermediates/EUK.phyloseq_03_rarefied_phyloseq.rds")
bac_dataset <- readRDS("intermediates/BAC.phyloseq_03_rarefied_phyloseq.rds")

# Estract and combine ASV tables
euk_tab <- data.frame(otu_table(euk_dataset))
bac_tab <- data.frame(otu_table(bac_dataset))

all_tab <- rbind(euk_tab, bac_tab)
all_tab <- t(all_tab)
all_tab[all_tab != 0] <- 1

# Calculate distances
distances_betapart <- beta.pair(all_tab, index.family = "jaccard")
jac_dist <- as.dist(distances_betapart$beta.jac)
tur_dist <- as.dist(distances_betapart$beta.jtu)

# NMDS ordination
set.seed(123678)
ord_jac <- metaMDS(jac_dist, k = 2, trymax = 100)
ord_tur <- metaMDS(tur_dist, k = 2, trymax = 100)

# Plot Jaccard
plot_jac <- data.frame(ord_jac$points)
plot_jac$sample <- rownames(plot_jac)
plot_jac <- merge(plot_jac, data.frame(sample_data(euk_dataset)))
plot_jac <- plot_jac[order(plot_jac$level), ]
plot_all_jac <- ggplot(plot_jac) +
  theme_bw() +
  geom_path(aes(x = MDS1, 
                y = MDS2, 
                col = as.factor(plant), 
                group = as.factor(plant))) +
  geom_point(aes(x = MDS1, 
                 y = MDS2, 
                 shape = as.factor(level), 
                 fill = as.factor(plant)), 
             size = 4, 
             color = "black", 
             alpha = 0.75) +
  scale_shape_manual(values = c(21, 22, 24)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  coord_fixed(ratio = diff(range(plot_jac$MDS1)) / diff(range(plot_jac$MDS2))) +
  ggtitle(paste0("Jaccard - Stress: ", round(ord_jac$stress, 4)))

# Plot Turnover
plot_tur <- data.frame(ord_tur$points)
plot_tur$sample <- rownames(plot_tur)
plot_tur <- merge(plot_tur, data.frame(sample_data(euk_dataset)))
plot_tur <- plot_tur[order(plot_tur$level), ]
plot_all_tur <- ggplot(plot_tur) +
  theme_bw() +
  geom_path(aes(x = MDS1, 
                y = MDS2, 
                col = as.factor(plant), 
                group = as.factor(plant))) +
  geom_point(aes(x = MDS1, 
                 y = MDS2, 
                 shape = as.factor(level), 
                 fill = as.factor(plant)), 
             size = 4, 
             color = "black", 
             alpha = 0.75) +
  scale_shape_manual(values = c(21, 22, 24)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  coord_fixed(ratio = diff(range(plot_tur$MDS1)) / diff(range(plot_tur$MDS2))) +
  ggtitle(paste0("Turnover - Stress: ", round(ord_tur$stress, 4)))

plot_nmds <- plot_all_jac | plot_all_tur
ggsave("./results/figures/plot_NMDS.pdf", 
       plot = plot_nmds, 
       height = 10, 
       width = 14)
