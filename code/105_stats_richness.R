#!/usr/bin/env Rscript

library(brms, quietly = TRUE)
library(performance, quietly = TRUE)
library(bayestestR, quietly = TRUE)
source("./code/999_utils.R")

# Load richness table
richness_all = read.table("intermediates/richness_tab.txt", sep = "\t", header = TRUE)
richness_all$level = factor(richness_all$level, levels = c(0,1,2), ordered = T)

########## RUN MODELS


### OBSERVED DIVERSITY

model_obs_meta = brm(Observed ~ mo(level) + (1|site) + (1|plant), data = subset(richness_all, taxa == "Metazoa"), family = poisson,
                     seed = 4321,
                     control = list(adapt_delta = 0.99, stepsize = 0.001, max_treedepth = 20),
                     chains = 4, iter = 40000, warmup = 5000, thin = 30)

model_obs_fungi = brm(Observed ~ mo(level) + (1|site) + (1|plant), data = subset(richness_all, taxa == "Fungi"), family = poisson,
                     seed = 4321,
                     control = list(adapt_delta = 0.99, stepsize = 0.001, max_treedepth = 20),
                     chains = 4, iter = 40000, warmup = 5000, thin = 30)

model_obs_algae = brm(Observed ~ mo(level) + (1|site) + (1|plant), data = subset(richness_all, taxa == "Algae"), family = poisson,
                     seed = 4321,
                     control = list(adapt_delta = 0.99, stepsize = 0.001, max_treedepth = 20),
                     chains = 4, iter = 40000, warmup = 5000, thin = 30)

model_obs_prot = brm(Observed ~ mo(level) + (1|site) + (1|plant), data = subset(richness_all, taxa == "Protists"), family = poisson,
                     seed = 4321,
                     control = list(adapt_delta = 0.99, stepsize = 0.001, max_treedepth = 20),
                     chains = 4, iter = 40000, warmup = 5000, thin = 30)

model_obs_bac = brm(Observed ~ mo(level) + (1|site) + (1|plant), data = subset(richness_all, taxa == "Bacteria"), family = poisson,
                     seed = 4321,
                     control = list(adapt_delta = 0.99, stepsize = 0.001, max_treedepth = 20),
                     chains = 4, iter = 40000, warmup = 5000, thin = 30)

