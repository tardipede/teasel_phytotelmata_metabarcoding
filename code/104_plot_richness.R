#!/usr/bin/env Rscript
library(ggplot2, quietly = TRUE)
library(patchwork, quietly = TRUE)
library(magrittr, quietly = TRUE)
options(tidyverse.quiet = TRUE)
library(tidyverse, quietly = TRUE)
source("./code/999_utils.R")

# Load richness table
richness_all = read.table("intermediates/richness_tab.txt", sep = "\t", header = TRUE)

# Remove Algae
richness_all = subset(richness_all , taxa != "Algae")

# Transform into long format
richness_all =  tidyr::pivot_longer(richness_all, c("Observed","Evenness","FaithPD")) %>%
  mutate(name = factor(name, levels = c("Observed","FaithPD","Evenness")),
         taxa = factor(taxa, levels = c("Bacteria","Fungi","Protists","Algae","Metazoa")))

plot_richness = ggplot(richness_all)+
  theme_bw()+
  geom_line(aes(x = level, y = value, group = as.factor(plant)), alpha = 0.3)+
  geom_point(aes(x = level, y = value, fill = site, shape = as.factor(level)), col = "black", size = 3, alpha = 0.75)+
  scale_x_continuous(breaks = c(0,1,2), labels = c(0,1,2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")+
  facet_wrap(name~taxa, scales = "free", nrow = 3) + xlab("Phytotelm level") + ylab("")+
  scale_shape_manual(values = c(21,22,24))

ggsave("./results/figures/plot_richness.pdf", plot = plot_richness, height = 10, width = 12)