#!/usr/bin/env Rscript

library(betapart, quietly = TRUE)
library(phyloseq, quietly = TRUE)
library(magrittr, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(plyr, quietly = TRUE)
library(brms, quietly = TRUE)
library(writexl, quietly = TRUE)
source("./code/999_utils.R")

# Load phyloseq objects
euk_dataset = readRDS("intermediates/EUK.phyloseq_03_rarefied_phyloseq.rds")
bac_dataset = readRDS("intermediates/BAC.phyloseq_03_rarefied_phyloseq.rds")

# Calculate prop of nestedness
d_nest_meta = plot_nestedness(subset_taxa(euk_dataset, Rank_6 == "Metazoa"), "nonphylo", "jaccard","data")
d_nest_fungi = plot_nestedness(subset_taxa(euk_dataset, Rank_5 == "Fungi"), "nonphylo", "jaccard", "data")

prot_data = subset_taxa(euk_dataset, Domain == "Eukaryota" & (((Rank_6 != "Metazoa") &
                                                                 (Rank_5 != "Fungi") & (Rank_2 != "Chloroplastida") & ( Rank_3 != "Ochrophyta"))| Rank_5 == "Chromulinales") & (Rank_1 != "unclassified_Eukaryota"))

d_nest_prot = plot_nestedness(prot_data, "nonphylo", "jaccard","data")
d_nest_bac = plot_nestedness(bac_dataset, "nonphylo", "jaccard","data")

d_nest_meta$type = mapvalues(d_nest_meta$type, from = c(110,100,101,1), to = c("A","B","C","D"))
d_nest_fungi$type = mapvalues(d_nest_fungi$type, from = c(110,100,101,1), to = c("A","B","C","D"))
d_nest_prot$type = mapvalues(d_nest_prot$type, from = c(110,100,101,1), to = c("A","B","C","D"))
d_nest_bac$type = mapvalues(d_nest_bac$type, from = c(110,100,101,1), to = c("A","B","C","D"))

d_nest_meta$type = factor(d_nest_meta$type, levels = c("A","B","C","D"), ordered = TRUE)
d_nest_fungi$type = factor(d_nest_fungi$type, levels = c("A","B","C","D"), ordered = TRUE)
d_nest_prot$type = factor(d_nest_prot$type, levels = c("A","B","C","D"), ordered = TRUE)
d_nest_bac$type = factor(d_nest_bac$type, levels = c("A","B","C","D"), ordered = TRUE)

d_nest_meta$value[d_nest_meta$value == 0] = exp(-20) ; d_nest_meta$value[d_nest_meta$value == 1] = 1- exp(-20) # Beta fmily does not accept 0 and 1, so they are replaced with ver close values
d_nest_fungi$value[d_nest_fungi$value == 0] = exp(-20) ; d_nest_fungi$value[d_nest_fungi$value == 1] = 1- exp(-20) # Beta fmily does not accept 0 and 1, so they are replaced with ver close values
d_nest_prot$value[d_nest_prot$value == 0] = exp(-20) ; d_nest_prot$value[d_nest_prot$value == 1] = 1- exp(-20) # Beta fmily does not accept 0 and 1, so they are replaced with ver close values
d_nest_bac$value[d_nest_bac$value == 0] = exp(-20) ; d_nest_bac$value[d_nest_bac$value == 1] = 1- exp(-20) # Beta fmily does not accept 0 and 1, so they are replaced with ver close values


## RUN Models
#If the monotonic effect is used in a linear model, 
#b can be interpreted as the expected average difference between two adjacent categories of the ordinal predictor.
# https://cran.r-project.org/web/packages/brms/vignettes/brms_monotonic.html#:~:text=If%20the%20monotonic%20effect%20is,categories%20of%20the%20ordinal%20predictor.


mod_nest_meta = brm(value ~ mo(type), data = d_nest_meta, family = Beta(),
               seed = 4321, iter = 10000, warmup = 5000, thin = 10)

mod_nest_fungi = brm(value ~ mo(type), data = d_nest_fungi, family = Beta(),
                seed = 4321, iter = 10000, warmup = 5000, thin = 10)

mod_nest_prot = brm(value ~ mo(type), data = d_nest_prot, family = Beta(),
               seed = 4321, iter = 10000, warmup = 5000, thin = 10)

mod_nest_bac = brm(value ~ mo(type), data = d_nest_bac, family = Beta(),
              seed = 4321, iter = 10000, warmup = 5000, thin = 10)



results_mods = list(
mod_nest_meta = clean_brm_output2(mod_nest_meta),
mod_nest_fungi = clean_brm_output2(mod_nest_fungi),
mod_nest_prot = clean_brm_output2(mod_nest_prot),
mod_nest_bac = clean_brm_output2(mod_nest_bac))

writexl::write_xlsx(results_mods, "./results/tables/models_nestedness.xlsx")


### Save all models in rds
all_models = list(
    mod_nest_meta,
    mod_nest_fungi,
    mod_nest_prot,
    mod_nest_bac
)

saveRDS(all_models, "./intermediates/brms_models_nestedness.rds")
