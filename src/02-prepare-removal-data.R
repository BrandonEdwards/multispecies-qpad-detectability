####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 02-prepare-removal-data.R
# Created May 2022
# Last Updated May 2022

####### Import Libraries and External Files #######

library(plyr)
set.seed(2022)

####### Read Data #################################

load("data/raw/time_count_matrix.rda")
load("data/raw/time_design.rda")
load("data/generated/turdidae_corr_matrix.rda")
binomial <- read.csv("data/generated/binomial_names.csv")

####### Create counts & design lists per species ##

# Most of this code adopted from Edwards et al. 2022

# Drop the 99 column
time_count_matrix <- time_count_matrix[, c(1:ncol(time_count_matrix) - 1)]

max_bands <- ncol(time_design) - 2
count_names <- c("Sample_ID", "Species", "Time_Method",
                 paste0(rep("Int", times = max_bands), 1:max_bands))
names(time_count_matrix) <- count_names

design_names <- c("Time_Method", "Max_Duration",
                  paste0(rep("Interval", times = max_bands), 1:max_bands))
names(time_design) <- design_names

# Join data
count_design <- plyr::join(time_count_matrix, time_design,
                           by = "Time_Method", type = "left")

# Filter out methods with only 1 interval
to_remove <- which(is.na(count_design$Interval2))
count_design <- count_design[-c(to_remove), ]

# Create separate data frames
counts <- count_design[,c("Sample_ID", "Species", "Time_Method", count_names[4:length(count_names)])]
design <- count_design[,c("Sample_ID", "Species", "Time_Method", design_names[3:length(design_names)])]
names(design) <- names(counts)
col_names <- names(counts)[4:length(names(counts))]

# Change 0s to NA in counts table where appropriate based on design table
for (i in col_names)
{
  indices <- which(counts[, i] == 0)
  
  counts[indices, i] <- ifelse(is.na(design[indices, i]), NA, 0)
}

# Re-order the species so that they match the correlation matrix
species <- sort(as.character(unique(counts$Species)))
bin_turdidae <- binomial[which(binomial$Code %in% species), ]
bin_turdidae$Scientific <- gsub(" ", "_", bin_turdidae$Scientific)
bin_turdidae <- bin_turdidae[match(colnames(corr_matrix), bin_turdidae$Scientific),]
species <- bin_turdidae$Code

# List of input counts
Y <- vector(mode = "list", length = length(species))
names(Y) <- species

# List of input design
D <- vector(mode = "list", length = length(species))
names(D) <- species

# List of design methods
methods <- vector(mode = "list", length = length(species))
names(methods) <- species

# Species indicator list
sp_list <- vector(mode = "list", length = length(species))
names(sp_list) <- species

for (s in species)
{
  n_counts <- nrow(counts[counts$Species == s, ])
  
  Y[[s]] <- as.matrix(counts[counts$Species==s, col_names])
  D[[s]] <- as.matrix(design[design$Species==s, col_names])
  methods[[s]] <- design[design$Species == s, "Time_Method"]
  sp_list[[s]] <- data.frame(Species = rep(s, n_counts))
}

####### Create Cross-validation Folds #############

CV_fold <- vector(mode = "list", length = length(species))
names(CV_fold) <- species

for (s in species)
{
  #' For each species, check which methods have less than 10 samples, and remove those.
  #' That way, when we do a stratified 10-fold CV, we can have at least one data point from
  #' each survey method contributing to the CV
  counts_per_method <- as.data.frame(table(methods[[s]]))
  methods_to_remove <- counts_per_method[which(counts_per_method$Freq < 10), "Var1"]
  indices_to_remove <- which(methods[[s]] %in% methods_to_remove)
  
  if (length(indices_to_remove) > 0)
  {
    Y[[s]] <- Y[[s]][-indices_to_remove, ]
    D[[s]] <- D[[s]][-indices_to_remove, ]
    methods[[s]] <- methods[[s]][-indices_to_remove]    
  }
  
  #' Now, for each method per species, assign a number between 1 and 10,
  #' inclusive, for each data point, and then shuffle WITHIN methods. Then
  #' put those vectors together

  counts_per_method <- as.data.frame(table(methods[[s]]))
  cv_fold_df <- data.frame(Method = methods[[s]],
                           Fold = NA)

  for (m in counts_per_method$Var1)
  {
    n_m <- counts_per_method[which(counts_per_method$Var1 == m), "Freq"]
    n_per_fold <- rep(trunc(n_m / 10), 10)
    remainder <- n_m %% 10 # modulus
    # Distribute the remainder
    i <- 1
    while (remainder != 0)
    {
      n_per_fold[i] <- n_per_fold[i] + 1
      remainder <- remainder - 1
      i <- i + 1
    }

    cv_method <- rep(seq(1,10),
                     times = n_per_fold)
    cv_method_shuffled <- sample(x = cv_method,
                                 size = length(cv_method),
                                 replace = FALSE)
    
    cv_fold_df[which(cv_fold_df$Method == m), "Fold"] <- cv_method_shuffled
    
  }
  
  CV_fold[[s]] <- cv_fold_df$Fold
}

####### Output ####################################

save(Y, D, CV_fold, sp_list, species, file = "removal_data.rda")
