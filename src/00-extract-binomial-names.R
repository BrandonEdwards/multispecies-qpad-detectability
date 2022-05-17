####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 00-extract-binomial-names.R
# Created March 2022
# Last Updated April 2022

####### Import Libraries and External Files #######


####### Read Data #################################

#' This file comes directly from NA-POPS: analysis/01-combine-data.R
#' which is not publicly available. Will be loaded in as "project_counts"
load("data/raw/counts.rda")

codes <- read.csv("data/raw/IBP-Alpha-Codes20.csv")
families <- read.csv("data/raw/NACC_list_species.csv")
birdtree <- read.csv("data/raw/birdtree_taxonomy.csv")

####### Main Code #################################

# Extract all binomial names for species that NA-POPS has data for
codes_red <- codes[which(codes$SPEC %in% unique(project_counts$Species)), ]
codes_family <- merge(codes_red, families[, c("species", "family")], 
                      by.x = "SCINAME", by.y = "species")

# Merge with BirdTree master taxonomy so we can use their scientific names
birdtree_red <- birdtree[which(birdtree$English %in% codes_family$COMMONNAME), ]
codes_family <- merge(codes_family, birdtree_red[, c("English", "Scientific")],
                      by.x = "COMMONNAME", by.y = "English", all = TRUE)

binomial <- codes_family[, c("SPEC", "COMMONNAME", "SCINAME", "family")]
names(binomial) <- c("Code", "English", "Scientific", "Family")

binomial[which(binomial$Scientific == "Ixoreus naevius"), "Scientific"] <- "Zoothera naevia"


####### Output ####################################

write.table(binomial,
            file = "data/generated/binomial_names.csv",
            sep = ",",
            row.names = FALSE)
