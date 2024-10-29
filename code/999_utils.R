#!/usr/bin/env Rscript

## LULU function from: https://github.com/tobiasgf/lulu/blob/master/R/Functions.R
## Modified to not print anything on screen and to save the log file in a specified path

lulu <- function(otutable, 
                 matchlist, 
                 minimum_ratio_type = "min", 
                 minimum_ratio = 1, 
                 minimum_match = 84, 
                 minimum_relative_cooccurence = 0.95, 
                 log_directory) {
  require(dplyr)
  start.time <- Sys.time()
  colnames(matchlist) <- c("OTUid", "hit", "match")
  # remove no hits (vsearch)
  matchlist <- matchlist[which(matchlist$hit != "*"), ]
  # remove self-hits
  matchlist <- matchlist[which(matchlist$hit != matchlist$OTUid), ]

  # Making a separate table with stats (total readcount and spread).
  statistics_table <- otutable[, 0]
  statistics_table$total <- rowSums(otutable)

  # calculating spread (number of presences (samples with 1+ read) pr OTU)
  statistics_table$spread <- rowSums(otutable > 0)
  statistics_table <- statistics_table[with(
    statistics_table,
    order(spread,
      total,
      decreasing = TRUE
    )
  ), ]
  otutable <- otutable[match(
    row.names(statistics_table),
    row.names(otutable)
  ), ]

  statistics_table$parent_id <- "NA"
  log_con <- file(paste0(log_directory, 
                         "lulu.log_", 
                         format(start.time, "%Y%m%d_%H%M%S")),
    open = "a"
  )
  for (line in seq(1:nrow(statistics_table))) {
    # make a progressline
    # print(paste0("progress: ",
    #             round(((line/nrow(statistics_table)) * 100), 0), "%"))
    potential_parent_id <- row.names(otutable)[line]
    cat(paste0("\n", "####processing: ", potential_parent_id, " #####"),
      file = log_con
    )
    daughter_samples <- otutable[line, ]
    hits <- matchlist[which(matchlist$OTUid == potential_parent_id &
      matchlist$match > minimum_match), "hit"]
    cat(paste0("\n", "---hits: ", hits), file = log_con)
    last_relevant_entry <- sum(statistics_table$spread >=
      statistics_table$spread[line])
    potential_parents <- which(row.names(otutable)[1:last_relevant_entry]
    %in% hits)
    cat(paste0(
      "\n", "---potential parent: ",
      row.names(statistics_table)[potential_parents]
    ), file = log_con)
    success <- FALSE
    if (length(potential_parents) > 0) {
      for (line2 in potential_parents) {
        cat(paste0("\n", "------checking: ", 
                   row.names(statistics_table)[line2]),
          file = log_con
        )
        if (!success) {
          relative_cooccurence <-
            sum((daughter_samples[otutable[line2, ] > 0]) > 0) /
              sum(daughter_samples > 0)
          cat(paste0(
            "\n", "------relative cooccurence: ",
            relative_cooccurence
          ), file = log_con)
          if (relative_cooccurence >= minimum_relative_cooccurence) {
            cat(paste0(" which is sufficient!"), file = log_con)
            if (minimum_ratio_type == "avg") {
              relative_abundance <-
                mean(otutable[line2, ][daughter_samples > 0] /
                  daughter_samples[daughter_samples > 0])
              cat(paste0(
                "\n", "------mean avg abundance: ",
                relative_abundance
              ), file = log_con)
            } else {
              relative_abundance <-
                min(otutable[line2, ][daughter_samples > 0] /
                  daughter_samples[daughter_samples > 0])
              cat(paste0(
                "\n", "------min avg abundance: ",
                relative_abundance
              ), file = log_con)
            }
            if (relative_abundance > minimum_ratio) {
              cat(paste0(" which is OK!"), file = log_con)
              if (line2 < line) {
                statistics_table$parent_id[line] <-
                  statistics_table[row.names(otutable)[line2], "parent_id"]
                cat(
                  paste0(
                    "\n", "SETTING ",
                    potential_parent_id, " to be an ERROR of ",
                    (statistics_table[
                      row.names(otutable)[line2],
                      "parent_id"
                    ]), "\n"
                  ),
                  file = log_con
                )
              } else {
                statistics_table$parent_id[line] <- row.names(otutable)[line2]
                cat(paste0(
                  "\n", "SETTING ", potential_parent_id,
                  " to be an ERROR of ", (row.names(otutable)[line2]),
                  "\n"
                ), file = log_con)
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
    rank = ave(total, FUN = function(x) {
      rank(-x, ties.method = "first")
    })
  )
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
  result <- list(
    curated_table = curation_table,
    curated_count = curated_count,
    curated_otus = curated_otus,
    discarded_count = discarded_count,
    discarded_otus = discarded_otus,
    runtime = time.taken,
    minimum_match = minimum_match,
    minimum_relative_cooccurence = minimum_relative_cooccurence,
    otu_map = statistics_table,
    original_table = otutable
  )

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
make_phyloseq <- function(paths, taxonomy_ranks) {
  require(ape)
  require(Biostrings)
  require(phyloseq)
  require(magrittr)
  require(tidyverse)

  # Get the different files
  otu_file <- paths[grepl(".csv", paths)]
  tree_file <- paths[grepl(".tre", paths)]
  seqs_file <- paths[grepl(".fas", paths)]
  metadata_file <- paths[grepl(".txt", paths)]


  # Load files
  otu_tab <- read.table(otu_file, 
                        sep = "\t", 
                        header = TRUE)
  rownames(otu_tab) <- otu_tab$ASV
  sample_dat <- read.table(metadata_file, 
                           sep = "\t", 
                           header = TRUE)
  rownames(sample_dat) <- sample_dat$sample
  tree <- ape::read.tree(tree_file)
  seqs <- Biostrings::readDNAStringSet(seqs_file)

  # Prepare the components to build a phyloseq object

  ## OTU table
  phyloseq_otu_tab <- phyloseq::otu_table(otu_tab[, colnames(otu_tab) %in% sample_dat$sample], 
                                          taxa_are_rows = TRUE)

  ## Taxonomy table
  ### Extract taxonomy from otutable
  taxonomy_table <- otu_tab$output %>%
    gsub("\\s*\\([^\\)]+\\)", "", .) %>%
    gsub(" ", "", .) %>%
    strsplit(., ";") %>%
    lapply(., FUN = function(x) {
      data.frame(t(x))
    }) %>%
    do.call(plyr::rbind.fill, .)

  taxonomy_table <- as.matrix(taxonomy_table)

  colnames(taxonomy_table) <- taxonomy_ranks
  rownames(taxonomy_table) <- otu_tab$ASV

  phyloseq_tax_tab <- phyloseq::tax_table(taxonomy_table)

  ## Sample data
  phyloseq_sample_data <- phyloseq::sample_data(sample_dat)

  ## Phylogenetic tree
  phyloseq_tree <- phyloseq::phy_tree(tree)

  ## Sequences
  phyloseq_seqs <- phyloseq::refseq(seqs)

  # Assembly phyoseq object
  dataset <- phyloseq(phyloseq_otu_tab, 
                      phyloseq_sample_data, 
                      phyloseq_tax_tab, 
                      phyloseq_tree, 
                      phyloseq_seqs)

  return(dataset)
}

##### Subtract reads in blank from all the samples
# Following https://doi.org/10.1002/edn3.372
subtract_blank <- function(phyloseq_object, blank_name) {
  require(phyloseq)


  otu_tab_temp <- data.frame(otu_table(phyloseq_object))
  otu_tab_temp <- sweep(otu_tab_temp, 
                        MARGIN = 1, 
                        STATS = otu_tab_temp[, blank_name], 
                        FUN = "-")
  otu_tab_temp[otu_tab_temp < 0] <- 0
  new_otu_table <- phyloseq::otu_table(otu_tab_temp, 
                                       taxa_are_rows = TRUE)

  phyloseq_object_new <- phyloseq(new_otu_table, 
                                  sample_data(phyloseq_object), 
                                  phy_tree(phyloseq_object), 
                                  tax_table(phyloseq_object), 
                                  refseq(phyloseq_object))

  return(phyloseq_object_new)
}

##### Save RDS and return object (to use with magrittr %>%)
saveRDS_pipe <- function(object_to_save, filename) {
  saveRDS(object_to_save, filename)
  return(object_to_save)
}

##### Phyloseq estimate_richness modified to calculate also Faith's Phylogenetic Distance 
## (see https://github.com/joey711/phyloseq/issues/661)
estimate_richness <- function(physeq, 
                              split = TRUE, 
                              measures = NULL) {
  require(picante)
  if (!any(otu_table(physeq) == 1)) {
    warning(
      "The data you have provided does not have\n",
      "any singletons. This is highly suspicious. Results of richness\n",
      "estimates (for example) are probably unreliable, or wrong, if you have already\n",
      "trimmed low-abundance taxa from the data.\n",
      "\n", "We recommended that you find the un-trimmed data and retry."
    )
  }
  if (!split) {
    OTU <- taxa_sums(physeq)
  } else if (split) {
    OTU <- as(otu_table(physeq), "matrix")
    if (taxa_are_rows(physeq)) {
      OTU <- t(OTU)
    }
  }
  renamevec <- c(
    "Observed", "Chao1", "ACE",
    "Shannon", "Simpson", "InvSimpson",
    "Fisher"
  )
  names(renamevec) <- c(
    "S.obs", "S.chao1", "S.ACE",
    "shannon", "simpson", "invsimpson",
    "fisher"
  )
  if (is.null(measures)) {
    measures <- as.character(renamevec)
  }
  if (any(measures %in% names(renamevec))) {
    measures[measures %in% names(renamevec)] <- renamevec[names(renamevec) %in%
      measures]
  }
  if (!any(measures %in% renamevec)) {
    stop("None of the `measures` you provided are supported. Try default `NULL` instead.")
  }
  outlist <- vector("list")
  estimRmeas <- c("Chao1", "Observed", "ACE")
  if (any(estimRmeas %in% measures)) {
    outlist <- c(outlist, list(t(data.frame(estimateR(OTU)))))
  }
  if ("Shannon" %in% measures) {
    outlist <- c(outlist, list(shannon = diversity(OTU, index = "shannon")))
  }
  if ("Simpson" %in% measures) {
    outlist <- c(outlist, list(simpson = diversity(OTU, index = "simpson")))
  }
  if ("InvSimpson" %in% measures) {
    outlist <- c(outlist, list(invsimpson = diversity(OTU,
      index = "invsimpson"
    )))
  }
  if ("FaithPD" %in% measures) {
    outlist <- c(outlist, list(FaithPD = t(picante::ses.pd(samp = OTU, 
                                                           tree = phytools::midpoint_root(phy_tree(physeq)), 
                                                           include.root = T, 
                                                           null.model = "richness"))["pd.obs.z", ]))
  }
  if ("Fisher" %in% measures) {
    fisher <- tryCatch(fisher.alpha(OTU, se = TRUE), warning = function(w) {
      warning("phyloseq::estimate_richness: Warning in fisher.alpha(). See `?fisher.fit` or ?`fisher.alpha`. Treat fisher results with caution")
      suppressWarnings(fisher.alpha(OTU, se = TRUE)[, c(
        "alpha",
        "se"
      )])
    })
    if (!is.null(dim(fisher))) {
      colnames(fisher)[1:2] <- c("Fisher", "se.fisher")
      outlist <- c(outlist, list(fisher))
    } else {
      outlist <- c(outlist, Fisher = list(fisher))
    }
  }
  out <- do.call("cbind", outlist)
  namechange <- intersect(colnames(out), names(renamevec))
  colnames(out)[colnames(out) %in% namechange] <- renamevec[namechange]
  colkeep <- sapply(paste0("(se\\.){0,}", measures), grep,
    colnames(out),
    ignore.case = TRUE
  )
  out <- out[, sort(unique(unlist(colkeep))), drop = FALSE]
  out <- as.data.frame(out)
  return(out)
}


#### Clean brmsfit output for stat richness
clean_brm_output <- function(brm_fit) {
  require(performance)
  require(bayestestR)

  chains <- as.data.frame(brm_fit)
  chains <- chains[, c("b_Intercept", 
                       "bsp_molevel", 
                       "sd_plant__Intercept", 
                       "sd_site__Intercept")]

  summary_chains <- as.data.frame(t(apply(chains, 
                                          MARGIN = 2, 
                                          FUN = function(x) {
    c(mean(x), as.numeric(bayestestR::hdi(x, ci = 0.95, verbose = TRUE))[2:3])})))
  colnames(summary_chains) <- c("Average", "CI95low", "CI95high")
  colnames(summary_chains) <- c("Average", "CI95low", "CI95high")

  p_dirs <- p_direction(chains[, 1:2])
  summary_chains$p <- c(pd_to_p(p_dirs$pd), NA, NA)

  r2vals <- r2(brm_fit)

  r2vals <- data.frame(rbind(
    c(r2vals$R2_Bayes_marginal, as.numeric(attributes(r2vals)$CI$R2_Bayes_marginal)[2:3]),
    c(r2vals$R2_Bayes, as.numeric(attributes(r2vals)$CI$R2_Bayes)[2:3]),
    rep(NA, 3), rep(NA, 3)
  ))
  colnames(r2vals) <- c("Average", "CI95low", "CI95high")
  r2vals <- cbind(data.frame(name = c("R2_marginal", 
                                      "R2_conditional", 
                                      "", 
                                      "")), 
                  r2vals)


  summary_chains <- cbind(data.frame(predictor = rownames(summary_chains)), 
                          summary_chains, r2vals)

  return(summary_chains)
}


### Function from https://github.com/vmikk/metagMisc/blob/master/R/dist2list.R
#' @title Convert distance matrix to data frame
#' @description This function takes a distance matrix (of class 'dist') and transforms it to a data.frame, where each row represents a single pairwise comparison.
#' @param dist Distance matrix (object of class 'dist')
#' @param tri Logical, if TRUE - only lower triangular part of dist will be returned
dist2list <- function(dist, tri = TRUE) {
  if (!class(dist) == "dist") {
    stop("Error: The input data must be a dist object.\n")
  }

  dat <- as.data.frame(as.matrix(dist))
  if (is.null(names(dat))) {
    rownames(dat) <- paste(1:nrow(dat))
  }
  value <- stack(dat)$values
  rnames <- rownames(dat)
  namecol <- expand.grid(rnames, rnames)
  colnames(namecol) <- c("col", "row")
  res <- data.frame(namecol, value)

  if (tri == TRUE) { # return only lower triangular part of dist
    res <- res[-which(upper.tri(as.matrix(dist), diag = TRUE)), ]
  }

  return(res)
}


### Function to calculate, extract and plot the % of nestedness in beta diversity
# type_distance options: nonphylo, phylo
# fam options = jaccard, sorensen
# return_what options: plot, data
plot_nestedness <- function(phyloseq_object, 
                            type_distance = "nonphylo", 
                            fam = "sorensen", 
                            return_what = "plot") {
  require(betapart)
  require(phyloseq)
  require(ggplot2)
  require(magrittr)
  require(dplyr)
  require(phytools)



  # Get sample data
  data_samples <- data.frame(sample_data(phyloseq_object))

  # Extract and binarize otu table
  table_binary <- t(data.frame(otu_table(phyloseq_object)))
  table_binary[table_binary != 0] <- 1

  # If distance = Jaccard
  if (type_distance == "nonphylo") {
    beta_res <- beta.pair(table_binary, index.family = fam)
    beta_res_nest_prop <- beta_res[[2]] / (beta_res[[2]] + beta_res[[1]])
  }

  # If distance = Unifrac
  if (type_distance == "phylo") {
    beta_res <- phylo.beta.pair(table_binary, 
                                tree = midpoint_root(phy_tree(phyloseq_object)), 
                                index.family = fam)
    beta_res_nest_prop <- beta_res[[2]] / (beta_res[[2]] + beta_res[[1]])
  }


  # Turn distance matrix into a list and prepare the dataframe with samples information
  beta_res_nest_prop <- dist2list(beta_res_nest_prop)
  beta_res_nest_prop <- merge(beta_res_nest_prop, 
                              data_samples[, c("sample", "site", "plant", "level")], 
                              by.x = "col", 
                              by.y = "sample")
  beta_res_nest_prop <- merge(beta_res_nest_prop, 
                              data_samples[, c("sample", "site", "plant", "level")], 
                              by.x = "row", 
                              by.y = "sample")

  beta_res_nest_prop <- beta_res_nest_prop %>%
    mutate(
      same_site = ifelse(site.x == site.y, 1, 0),
      same_plant = ifelse(plant.x == plant.y, 1, 0),
      same_level = ifelse(level.x == level.y, 1, 0),
      type = same_site * 100 + same_plant * 10 + same_level
    ) %>%
    mutate(type = ifelse(type %in% c(0, 1), 1, type)) %>%
    mutate(type = factor(type, 
                         levels = c("110", "100", "101", "1")))

  # Make plot

  res_plot <- ggplot(beta_res_nest_prop) +
    theme_bw() +
    geom_point(aes(x = as.factor(type), y = value, col = type),
      position = position_jitter(width = 0.25, seed = 987654321),
      alpha = 0.5, show.legend = F, size = 2
    ) +
    geom_boxplot(aes(x = as.factor(type), 
                     y = value), 
                 fill = NA, 
                 outlier.alpha = 0) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    ylim(c(-0.01, 1.01)) +
    ylab("Proportion nestedness") +
    xlab("") +
    scale_x_discrete(breaks = c("110", "101", "100", "1"), 
                     labels = c("Same plant", "Same site,\nsame level", "Same site,\ndifferent level", "Different site"))


  if (return_what == "plot") {
    return(res_plot)
  }
  if (return_what == "data") {
    return(beta_res_nest_prop)
  }
}

#### Clean brmsfit output for nestedness stats
clean_brm_output2 <- function(brm_fit) {
  require(performance)
  require(bayestestR)

  chains <- as.data.frame(brm_fit)
  chains <- chains[, c("b_Intercept", 
                       "bsp_motype", 
                       "simo_motype1[1]", 
                       "simo_motype1[2]", 
                       "simo_motype1[3]")]

  summary_chains <- as.data.frame(t(apply(chains, 
                                          MARGIN = 2, 
                                          FUN = function(x) {
    c(mean(x), as.numeric(bayestestR::hdi(x, ci = 0.95, verbose = TRUE))[2:3])})))
  colnames(summary_chains) <- c("Average", 
                                "CI95low", 
                                "CI95high")
  rownames(summary_chains) <- c("Intercept", 
                                "type", 
                                "simplexAB", 
                                "simplexBC", 
                                "simplexCD")

  p_dirs <- p_direction(chains)
  summary_chains$p <- pd_to_p(p_dirs$pd)
  summary_chains$p[3:5] <- rep(NA, 3)

  r2vals <- performance::r2(brm_fit)

  r2vals <- data.frame(rbind(
    c(r2vals$R2_Bayes, as.numeric(attributes(r2vals)$CI$R2_Bayes)[2:3]),
    rep(NA, 3), rep(NA, 3), rep(NA, 3), rep(NA, 3)
  ))
  colnames(r2vals) <- c("Average", 
                        "CI95low", 
                        "CI95high")
  r2vals <- cbind(data.frame(name = c("R2", 
                                      "", 
                                      "", 
                                      "", 
                                      "")), 
                  r2vals)


  summary_chains <- cbind(data.frame(predictor = rownames(summary_chains)), 
                          summary_chains, r2vals)

  return(summary_chains)
}
