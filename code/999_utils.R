#!/usr/bin/env Rscript

## LULU function from: https://github.com/tobiasgf/lulu/blob/master/R/Functions.R
## Modified to not print anything on screen and to save the log file in a specified path

lulu <- function(otutable, matchlist, minimum_ratio_type = "min", minimum_ratio = 1, minimum_match = 84, minimum_relative_cooccurence = 0.95, log_directory) {
  require(dplyr)
  start.time <- Sys.time()
  colnames(matchlist) <- c("OTUid", "hit", "match")
  # remove no hits (vsearch)
  matchlist = matchlist[which(matchlist$hit != "*"), ]
  # remove self-hits
  matchlist = matchlist[which(matchlist$hit != matchlist$OTUid), ]

  # Making a separate table with stats (total readcount and spread).
  statistics_table <- otutable[, 0]
  statistics_table$total <- rowSums(otutable)

  # calculating spread (number of presences (samples with 1+ read) pr OTU)
  statistics_table$spread <- rowSums(otutable > 0)
  statistics_table <- statistics_table[with(statistics_table,
                                            order(spread,
                                                  total,
                                                  decreasing = TRUE)), ]
  otutable <- otutable[match(row.names(statistics_table),
                             row.names(otutable)), ]

  statistics_table$parent_id <- "NA"
  log_con <- file(paste0(log_directory,"lulu.log_", format(start.time, "%Y%m%d_%H%M%S")),
                  open = "a")
  for (line in seq(1:nrow(statistics_table))) {
    # make a progressline
    #print(paste0("progress: ",
    #             round(((line/nrow(statistics_table)) * 100), 0), "%"))
    potential_parent_id <- row.names(otutable)[line]
    cat(paste0("\n", "####processing: ", potential_parent_id, " #####"),
        file = log_con)
    daughter_samples <- otutable[line, ]
    hits <- matchlist[which(matchlist$OTUid == potential_parent_id &
                              matchlist$match > minimum_match), "hit"]
    cat(paste0("\n", "---hits: ", hits), file = log_con)
    last_relevant_entry <- sum(statistics_table$spread >=
                                 statistics_table$spread[line])
    potential_parents <- which(row.names(otutable)[1:last_relevant_entry]
                               %in% hits)
    cat(paste0("\n", "---potential parent: ",
               row.names(statistics_table)[potential_parents]), file = log_con)
    success <- FALSE
    if (length(potential_parents) > 0) {
      for (line2 in potential_parents) {
        cat(paste0("\n", "------checking: ", row.names(statistics_table)[line2]),
            file = log_con)
        if (!success) {
          relative_cooccurence <-
            sum((daughter_samples[otutable[line2, ] > 0]) > 0)/
            sum(daughter_samples > 0)
          cat(paste0("\n", "------relative cooccurence: ",
                     relative_cooccurence), file = log_con)
          if (relative_cooccurence >= minimum_relative_cooccurence) {
            cat(paste0(" which is sufficient!"), file = log_con)
            if (minimum_ratio_type == "avg") {
              relative_abundance <-
                mean(otutable[line2, ][daughter_samples > 0]/
                       daughter_samples[daughter_samples > 0])
              cat(paste0("\n", "------mean avg abundance: ",
                         relative_abundance), file = log_con)
            } else {
              relative_abundance <-
                min(otutable[line2, ][daughter_samples > 0]/
                      daughter_samples[daughter_samples > 0])
              cat(paste0("\n", "------min avg abundance: ",
                         relative_abundance), file = log_con)
            }
            if (relative_abundance > minimum_ratio) {
              cat(paste0(" which is OK!"), file = log_con)
              if (line2 < line) {
                statistics_table$parent_id[line] <-
                  statistics_table[row.names(otutable)[line2],"parent_id"]
                cat(paste0("\n", "SETTING ",
                           potential_parent_id, " to be an ERROR of ",
                           (statistics_table[row.names(otutable)[line2],
                                             "parent_id"]), "\n"),
                    file = log_con)
              } else {
                statistics_table$parent_id[line] <- row.names(otutable)[line2]
                cat(paste0("\n", "SETTING ", potential_parent_id,
                           " to be an ERROR of ", (row.names(otutable)[line2]),
                           "\n"), file = log_con)
              }
              success <- TRUE
            }
          }
        }
      }
    }
    if (!success) {
      statistics_table$parent_id[line] <- row.names(statistics_table)[line]
      cat(paste0("\n", "No parent found!", "\n"), file = log_con)
    }
  }

  close(log_con)
  total_abundances <- rowSums(otutable)
  curation_table <- cbind(nOTUid = statistics_table$parent_id, otutable)
  statistics_table$curated <- "merged"
  curate_index <- row.names(statistics_table) == statistics_table$parent_id
  statistics_table$curated[curate_index] <- "parent"
  statistics_table <- transform(statistics_table,
                                rank = ave(total,FUN = function(x)
                                  rank(-x, ties.method = "first")))
  curation_table <- as.data.frame(curation_table %>%
                                    group_by(nOTUid) %>%
                                    summarise_all(list(sum)))
  row.names(curation_table) <- as.character(curation_table$nOTUid)
  curation_table <- curation_table[, -1]
  curated_otus <- names(table(statistics_table$parent_id))
  curated_count <- length(curated_otus)
  discarded_otus <- setdiff(row.names(statistics_table), curated_otus)
  discarded_count <- length(discarded_otus)
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  result <- list(curated_table = curation_table,
                 curated_count = curated_count,
                 curated_otus = curated_otus,
                 discarded_count = discarded_count,
                 discarded_otus = discarded_otus,
                 runtime = time.taken,
                 minimum_match = minimum_match,
                 minimum_relative_cooccurence = minimum_relative_cooccurence,
                 otu_map = statistics_table,
                 original_table = otutable)

  return(result)
}


##### Import al files in a folder to create a phyloseq object
### In the folder you need 4 files
## OTU table with taxonomy (file name ends with .csv) and must contain:
# Column named "output" with taxonomy 
# Column names "OTU" with Otu names
# Columns with samples (columns names are sample names)
## Sample data (file name ends with.txt) and must contain:
# Column named "sample" with samples names
# Any other column with additional samples data
## Phylogenetic tree (file name ends with .tre)
## Aligned OTUs (file name ends with .fas)
###### path indicates the folder path where all these files are present
make_phyloseq = function(paths,taxonomy_ranks){
  require(ape)
  require(Biostrings)
  require(phyloseq)
  require(magrittr)
  require(tidyverse)

  # Get the different files
  otu_file = paths[grepl(".csv",paths)]
  tree_file = paths[grepl(".tre",paths)]
  seqs_file = paths[grepl(".fas",paths)]
  metadata_file = paths[grepl(".txt",paths)]


  # Load files
  otu_tab = read.table(otu_file, sep = "\t", header = T);  rownames(otu_tab) = otu_tab$ASV
  sample_dat = read.table(metadata_file, sep = "\t", header = T); rownames(sample_dat) = sample_dat$sample
  tree = ape::read.tree(tree_file)
  seqs = Biostrings::readDNAStringSet(seqs_file)
  
  # Prepare the components to build a phyloseq object
  
  ## OTU table
  phyloseq_otu_tab = phyloseq::otu_table(otu_tab[,colnames(otu_tab) %in% sample_dat$sample], taxa_are_rows = T)
  
  ## Taxonomy table
  ### Extract taxonomy from otutable
  taxonomy_table  = otu_tab$output %>%
    gsub("\\s*\\([^\\)]+\\)", "", .) %>%
    gsub(" ","",.) %>%
    strsplit(.,";") %>%
    lapply(., FUN = function(x){data.frame(t(x))}) %>%
    do.call(plyr::rbind.fill, .)
  
  taxonomy_table = as.matrix(taxonomy_table)
  
  colnames(taxonomy_table) = taxonomy_ranks
  rownames(taxonomy_table) = otu_tab$ASV
  
  phyloseq_tax_tab = phyloseq::tax_table(taxonomy_table)
  
  ## Sample data
  phyloseq_sample_data = phyloseq::sample_data(sample_dat)
  
  ## Phylogenetic tree
  phyloseq_tree = phyloseq::phy_tree(tree)
  
  ## Sequences
  phyloseq_seqs = phyloseq::refseq(seqs)
  
  # Assembly phyoseq object
  dataset = phyloseq(phyloseq_otu_tab, phyloseq_sample_data, phyloseq_tax_tab, phyloseq_tree, phyloseq_seqs)
    
  return(dataset)
}

##### Subtract reads in blank from all the samples
# Following https://doi.org/10.1002/edn3.372
subtract_blank = function(phyloseq_object, blank_name){
  require(phyloseq)

  
  otu_tab_temp = data.frame(otu_table(phyloseq_object))
  otu_tab_temp = sweep(otu_tab_temp, MARGIN = 1, STATS = otu_tab_temp[,blank_name], FUN = "-")
  otu_tab_temp[otu_tab_temp<0] = 0
  new_otu_table = phyloseq::otu_table(otu_tab_temp, taxa_are_rows = T)
  
  phyloseq_object_new = phyloseq(new_otu_table, sample_data(phyloseq_object), phy_tree(phyloseq_object), tax_table(phyloseq_object), refseq(phyloseq_object))
  
 # samples_to_keep = data.frame(sample_data(phyloseq_object_new))$sample != blank_name
  
  #otu_table(phyloseq_object) = new_otu_table
  #phyloseq_object_new = phyloseq::subset_samples(phyloseq_object_new, samples_to_keep)
  
  return(phyloseq_object_new)
}

##### Save RDS and return object (to use with magrittr %>%)
saveRDS_pipe = function(object_to_save, filename){
  saveRDS(object_to_save, filename)
  return(object_to_save)
}
