####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 2-prepare-phylogenies.R
# Created August 2022
# Last Updated March 2024

####### Import Libraries and External Files #######

library(ape)
library(magrittr)

source("src/functions/generate-phylo-corr.R")

####### Read Data #################################

phylo_tree <- ape::read.nexus(file = "data/raw/all_species.nex")

####### Main Code #################################

# Generate the correlation matrix for the cross validation (takes a long time)
corr_matrix <- generate_phylo_corr(phylo_tree = phylo_tree)

save(corr_matrix, file = "data/generated/corr_matrix_all.rda")
