#!/usr/bin/bash

conda activate DP_metabar # activate conda env (created with mamba env create -f environment.yml)




# Process reads for 18S and 16S datasets
Rscript code/101_setup_reads_process.R "EUK" "18S" "18S"
Rscript code/101_setup_reads_process.R "BAC" "16S" "16S"

# Create R-phyloseq objects
Rscript code/102_create_phyloseq_objects.R