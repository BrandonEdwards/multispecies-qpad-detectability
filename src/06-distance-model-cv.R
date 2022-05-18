####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 06-distance-model-cv.R
# Created May 2022
# Last Updated May 2022

####### Import Libraries and External Files #######

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

####### Read Data #################################

load("data/generated/distance_data.rda")
load("data/generated/turdidae_corr_matrix.rda")
binomial <- read.csv("data/generated/binomial_names.csv")
model <- stan_model(file = "models/distance.stan")

####### Wrangle Data for Modelling ################

Y <- Y[lengths(Y) != 0]; Y <- do.call(rbind, Y)
D <- D[lengths(D) != 0]; D <- do.call(rbind, D)
sp_list <- sp_list[lengths(sp_list) != 0]; sp_list <- do.call(rbind, sp_list)
CV_fold <- CV_fold[lengths(CV_fold) != 0]; CV_fold <- unname(unlist(CV_fold))

for (fold in 1:10)
{
  fold_indices <- which(CV_fold == fold)
  Y_train <- Y[-indices_to_remove, ]
  D_train <- D[-indices_to_remove, ]
  sp_list_train <- sp_list[-indices_to_remove, ]
  
  #' Corresponds with "bands_per_sample" in distance.stan
  dist_bands_per_sample <- unname(apply(Y_train, 1, function(x) sum(!is.na(x))))
  
  #' Total species abundance per sampling event.
  #' I.e., this is the sum of Y_sik over k
  #' Corresponds with "abund_per_sample" in distance.stan
  total_abund_per_sample <- unname(apply(Y_train, 1, function(x) sum(x, na.rm = TRUE)))
  
  #' Factored list of species
  #' Corresponds with "species" in distance.stan
  sp_list_numeric <- as.numeric(factor(sp_list_train[,1],
                                       levels = species))
  
  #' Corresponds with "abund_per_band" in distance.stan
  abundance_per_band <- Y_train
  abundance_per_band[is.na(abundance_per_band)] <- 0
  
  #' Corresponds with "max_dist" in distance.stan
  max_dist <- D_train
  max_dist[is.na(max_dist)] <- 0
  
  n_samples <- nrow(Y_train)
  n_species <- length(unique(sp_list_train[,1]))
  max_intervals <- ncol(Y_train)
  
  stan_data <- list(n_samples = n_samples,
                    n_species = n_species,
                    max_intervals = max_intervals,
                    species = sp_list_numeric,
                    abund_per_band = abundance_per_band,
                    abund_per_sample = total_abund_per_sample,
                    bands_per_sample = dist_bands_per_sample,
                    max_dist = max_dist,
                    phylo_corr = corr_matrix)
  
  ####### Output ####################################
  
  model <- stan_model(file = "models/distance.stan")
  
  stan_job <- sampling(model,
                       data = stan_data,
                       verbose = TRUE,
                       chains = 4,
                       iter = 1000,
                       warmup = 500,
                       cores = 4,
                       pars = c("log_tau", "sigma", "mu"),
                       #control = list(adapt_delta = 0.8,
                       #              max_treedepth = 15),
                       init = list(list(log_tau = rnorm(11, 5.5, 0.2),
                                        mu = rnorm(1, 5.5, 0.2),
                                        sigma = rep(1,11)),
                                   list(log_tau = rnorm(11, 5.5, 0.2),
                                        mu = rnorm(1, 5.5, 0.2),
                                        sigma = rep(1,11)),
                                   list(log_tau = rnorm(11, 5.5, 0.2),
                                        mu = rnorm(1, 5.5, 0.2),
                                        sigma = rep(1,11)),
                                   list(log_tau = rnorm(11, 5.5, 0.2),
                                        mu = rnorm(1, 5, 0.2),
                                        sigma = rep(1,11))))
  
  ####### Output ####################################
  
  save(stan_job, file = paste0("data/generated/distance_turdidae_model_cv",
                               fold,
                               ".rda"))

}
