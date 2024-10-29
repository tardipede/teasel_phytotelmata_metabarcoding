#!/usr/bin/bash

conda activate DP_metabar # activate conda env (created with "mamba env create -f environment.yml")

# Process reads for 18S and 16S datasets
Rscript code/101_setup_reads_process.R "EUK" "18S" "18S" "(ID[0-9]{4})" "((DFP)-.{3})"
Rscript code/101_setup_reads_process.R "BAC" "16S" "16S" "(ID[0-9]{4})" "((DFP)-.{3})"

# Create R-phyloseq objects
Rscript code/102_create_phyloseq_objects.R

# Calculate alpha diversity
Rscript code/103_calculate_richness.R

# Plot alpha richness
Rscript code/104_plot_richness.R

# Do stats on richness
Rscript code/105_stats_richness.R

# NMDS ordinations
Rscript code/106_ordinations.R

# Plot nestedness proportion
Rscript code/107_plot_nestedness.R

# Do stats on nestedness
Rscript code/108_stats_nestedness.R

# Taxonomic composition plots
Rscript code/109_taxonomic_compositon_barplots.R

# Save otu tables with taxonomy as excel files
Rscript code/110_save_otu_tables_in_excel.R