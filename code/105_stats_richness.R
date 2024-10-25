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

model_obs_prot = brm(Observed ~ mo(level) + (1|site) + (1|plant), data = subset(richness_all, taxa == "Protists"), family = poisson,
                     seed = 4321,
                     control = list(adapt_delta = 0.99, stepsize = 0.001, max_treedepth = 20),
                     chains = 4, iter = 40000, warmup = 5000, thin = 30)

model_obs_bac = brm(Observed ~ mo(level) + (1|site) + (1|plant), data = subset(richness_all, taxa == "Bacteria"), family = poisson,
                     seed = 4321,
                     control = list(adapt_delta = 0.99, stepsize = 0.001, max_treedepth = 20),
                     chains = 4, iter = 40000, warmup = 5000, thin = 30)

models_observed = list(Metazoa = clean_brm_output(model_obs_meta),
                       Fungi = clean_brm_output(model_obs_fungi),
                       Protists = clean_brm_output(model_obs_prot),
                       Bacteria = clean_brm_output(model_obs_bac))

writexl::write_xlsx(models_observed, "./results/tables/models_observed.xlsx")



### PHYLOGENETIC DIVERSITY


model_pd_meta = brm(FaithPD ~ mo(level) + (1|site) + (1|plant), data = subset(richness_all, taxa == "Metazoa"), family = gaussian,
                    seed = 4321,
                    control = list(adapt_delta = 0.99, stepsize = 0.001, max_treedepth = 20),
                    chains = 4, iter = 40000, warmup = 5000, thin = 30)

model_pd_fungi = brm(FaithPD ~ mo(level) + (1|site) + (1|plant), data = subset(richness_all, taxa == "Fungi"), family = gaussian,
                     seed = 4321,
                     control = list(adapt_delta = 0.99, stepsize = 0.001, max_treedepth = 20),
                     chains = 4, iter = 40000, warmup = 5000, thin = 30)


model_pd_prot = brm(FaithPD ~ mo(level) + (1|site) + (1|plant), data = subset(richness_all, taxa == "Protists"), family = gaussian,
                    seed = 4321,
                    control = list(adapt_delta = 0.99, stepsize = 0.001, max_treedepth = 20),
                    chains = 4, iter = 40000, warmup = 5000, thin = 30)

model_pd_bac = brm(FaithPD ~ mo(level) + (1|site) + (1|plant), data = subset(richness_all, taxa == "Bacteria"), family = gaussian,
                   seed = 4321,
                   control = list(adapt_delta = 0.99, stepsize = 0.001, max_treedepth = 20),
                   chains = 4, iter = 40000, warmup = 5000, thin = 30)




models_pd = list(Metazoa = clean_brm_output(model_pd_meta),
                 Fungi = clean_brm_output(model_pd_fungi),
                 Protists = clean_brm_output(model_pd_prot),
                 Bacteria = clean_brm_output(model_pd_bac))


writexl::write_xlsx(models_pd, "./results/tables/models_pd.xlsx")

### EVENNESS


model_eve_meta = brm(Evenness ~ mo(level) + (1|site) + (1|plant), data = subset(richness_all, taxa == "Metazoa"), family = Beta(link = "logit", link_phi = "log"),
                     seed = 4321,
                     control = list(adapt_delta = 0.99, stepsize = 0.001, max_treedepth = 20),
                     chains = 4, iter = 40000, warmup = 5000, thin = 30)

model_eve_fungi = brm(Evenness ~ mo(level) + (1|site) + (1|plant), data = subset(richness_all, taxa == "Fungi"), family = Beta(link = "logit", link_phi = "log"),
                      seed = 4321,
                      control = list(adapt_delta = 0.99, stepsize = 0.001, max_treedepth = 20),
                      chains = 4, iter = 40000, warmup = 5000, thin = 30)


model_eve_prot = brm(Evenness ~ mo(level) + (1|site) + (1|plant), data = subset(richness_all, taxa == "Protists"), family = c_bernoulli, stanvars = stanvars,
                     seed = 4321,
                     control = list(adapt_delta = 0.99, stepsize = 0.001, max_treedepth = 20),
                     chains = 4, iter = 40000, warmup = 5000, thin = 30)

model_eve_bac = brm(Evenness ~ mo(level) + (1|site) + (1|plant), data = subset(richness_all, taxa == "Bacteria"), family = Beta(link = "logit", link_phi = "log"),
                    seed = 4321,
                    control = list(adapt_delta = 0.99, stepsize = 0.001, max_treedepth = 20),
                    chains = 4, iter = 40000, warmup = 5000, thin = 30)




models_eve = list(Metazoa = clean_brm_output(model_eve_meta),
                  Fungi = clean_brm_output(model_eve_fungi),
                  Protists = clean_brm_output(model_eve_prot),
                  Bacteria = clean_brm_output(model_eve_bac))


writexl::write_xlsx(models_eve, "./results/tables/models_evenness.xlsx")


### Save all models in rds
all_models = list(
    model_obs_meta,
    model_obs_fungi,
    model_obs_prot,
    model_obs_bac,
    model_pd_meta,
    model_pd_fungi,
    model_pd_prot,
    model_pd_bac,
    model_pd_meta,
    model_pd_fungi,
    model_pd_prot,
    model_pd_bac)

saveRDS(all_models, "./intermediates/brms_models_alpha_diversity.rds")