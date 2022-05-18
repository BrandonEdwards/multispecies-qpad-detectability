####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 03-removal-model-full.R
# Created April 2022
# Last Updated May 2022

####### Import Libraries and External Files #######

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

####### Read Data #################################

load("data/generated/removal_data.rda")
load("data/generated/turdidae_corr_matrix.rda")
binomial <- read.csv("data/generated/binomial_names.csv")

####### Wrangle Data for Modelling ################

Y <- Y[lengths(Y) != 0]; Y <- do.call(rbind, Y)
D <- D[lengths(D) != 0]; D <- do.call(rbind, D)
sp_list <- sp_list[lengths(sp_list) != 0]; sp_list <- do.call(rbind, sp_list)

#' Corresponds with "bands_per_sample" in removal.stan
time_bands_per_sample <- unname(apply(Y, 1, function(x) sum(!is.na(x))))

#' Total species abundance per sampling event.
#' I.e., this is the sum of Y_sij over j
#' Corresponds with "abund_per_sample" in removal.stan
total_abund_per_sample <- unname(apply(Y, 1, function(x) sum(x, na.rm = TRUE)))

#' Factored list of species
#' Corresponds with "species" in removal.stan
sp_list_numeric <- as.numeric(factor(sp_list[,1],
                                     levels = species))

#' Corresponds with "abund_per_band" in removal.stan
abundance_per_band <- Y
abundance_per_band[is.na(abundance_per_band)] <- 0

#' Corresponds with "max_time" in removal.stan
max_time <- D
max_time[is.na(max_time)] <- 0

n_samples <- nrow(Y)
n_species <- length(unique(sp_list[,1]))
max_intervals <- ncol(Y)

stan_data <- list(n_samples = n_samples,
                  n_species = n_species,
                  max_intervals = max_intervals,
                  species = sp_list_numeric,
                  abund_per_band = abundance_per_band,
                  abund_per_sample = total_abund_per_sample,
                  bands_per_sample = time_bands_per_sample,
                  max_time = max_time,
                  phylo_corr = corr_matrix)

# Output this here in case the stan model goes silly
save(stan_data, file = "data/generated/removal_turdidae_data.rda")

####### Stan Modelling ############################

model <- stan_model(file = "models/removal.stan")

stan_job <- sampling(model,
                          data = stan_data,
                          verbose = TRUE,
                          chains = 4,
                          iter = 1000,
                          warmup = 500,
                          cores = 4,
                          pars = c("log_phi", "tau", "mu"),
                          control = list(adapt_delta = 0.8,
                                         max_treedepth = 15))

####### Output ####################################

save(stan_job, file = "data/generated/removal_turdidae_model_full.rda")
