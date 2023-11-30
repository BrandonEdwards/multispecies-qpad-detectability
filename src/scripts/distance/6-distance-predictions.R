####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 6-distance-predictions.R
# Created October 2022
# Last Updated October 2023

####### Import Libraries and External Files #######

library(cmdstanr)
library(dplyr)

source("src/functions/generate-distance-inits.R")
source("src/functions/subset-distance-data.R")

####### Load Data #################################

load("data/generated/distance_stan_data.rda")

####### Set Constants #############################

# Stan settings
n_iter <- 2000
n_warmup <- 1000
n_chains <- 4
refresh <- 10
threads_per_chain <- 3

ruhu <- subset_distance_data(distance_stan_data = distance_stan_data,
                                           sps = "RUHU")

alhu <- subset_distance_data(distance_stan_data = distance_stan_data,
                             sps = "ALHU")

rthu <- subset_distance_data(distance_stan_data = distance_stan_data,
                             sps = "RTHU")

bchu <- subset_distance_data(distance_stan_data = distance_stan_data,
                             sps = "BCHU")

distance_stan_data$grainsize <- 1

distance_stan_data$sp_list <- NULL
distance_stan_data$sp_all <- NULL

# Scale the maximum distances for computational ease
distance_stan_data$max_dist <- distance_stan_data$max_dist/10

####### Run Model #################################

ss_model <- cmdstan_model(stan_file = "models/distance_single.stan",
                          cpp_options = list(stan_threads = TRUE))

stan_run <- ss_model$sample(
  data = distance_stan_data,
  iter_warmup = n_warmup,
  iter_sampling = n_iter,
  chains = n_chains,
  parallel_chains = n_chains,
  refresh = refresh,
  threads_per_chain = threads_per_chain
)

draws <- stan_run$sampler_diagnostics(format = "df")
draws$log_tau <- stan_run$draws(variables = "log_tau[1]", format = "df")$`log_tau[1]`
