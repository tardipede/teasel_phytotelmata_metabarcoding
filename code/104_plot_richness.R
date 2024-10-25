#!/usr/bin/env Rscript
library(ggplot2, quietly = TRUE)
library(patchwork, quietly = TRUE)
library(magrittr, quietly = TRUE)
options(tidyverse.quiet = TRUE)
library(tidyverse, quietly = TRUE)
source("./code/999_utils.R")

# Load richness table
richness_all = read.table("intermediates/richness_tab.txt", sep = "\t", header = T)