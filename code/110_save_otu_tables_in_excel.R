#!/usr/bin/env Rscript

library(phyloseq, , quietly = TRUE)
library(writexl, , quietly = TRUE)

# Load phyloseq objects
euk_dataset = readRDS("intermediates/EUK.phyloseq_03_rarefied_phyloseq.rds")
bac_dataset = readRDS("intermediates/BAC.phyloseq_03_rarefied_phyloseq.rds")

# EUK
asv_table = data.frame(otu_table(euk_dataset)); asv_table = cbind(data.frame(ASV = rownames(asv_table)), asv_table)
tax_table = data.frame(tax_table(euk_dataset)); tax_table = cbind(data.frame(ASV = rownames(tax_table)), tax_table)
euk_table = merge(asv_table, tax_table, by = "ASV")
writexl::write_xlsx(list(EUK_otu_table = euk_table), "./results/otu_tables_with_taxonomy/EUK_otu_table.xlsx")

# BAC
asv_table = data.frame(otu_table(bac_dataset)); asv_table = cbind(data.frame(ASV = rownames(asv_table)), asv_table)
tax_table = data.frame(tax_table(bac_dataset)); tax_table = cbind(data.frame(ASV = rownames(tax_table)), tax_table)
bac_table = merge(asv_table, tax_table, by = "ASV")
writexl::write_xlsx(list(BAC_otu_table = bac_table), "./results/otu_tables_with_taxonomy/BAC_otu_table.xlsx")
