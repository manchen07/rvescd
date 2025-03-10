# set up simulations
args <- commandArgs(trailingOnly = TRUE)

library(dplyr, quietly = TRUE, warn.conflicts = FALSE)
library(tidyr)
library(purrr)
library(metafor)
library(clubSandwich)
library(simhelpers)
library(mvtnorm)
library(puniform)
library(stringr)
library(plyr)


#-------------------------------------
# Parse command line arguments
#-------------------------------------
extract_int_arg <- function(arg_string, arg_name, default) {
  res <- 
    arg_string |>
    stringr::str_extract(paste(arg_name, "[0-9]+")) |>
    stringr::str_sub(stringr::str_length(arg_name) + 2, -1) |>
    as.integer()
  
  if (is.na(res)) res <- default
  res
}

arg_string <- paste(args, collapse = " ")

cores <- extract_int_arg(arg_string, "-cores", 2L)
batch <- extract_int_arg(arg_string, "-batch", 1L)
reps <- extract_int_arg(arg_string, "-reps", 2L)

#-------------------------------------
# Source the function
#-------------------------------------

source("SR-sim-functions.R") # condor
# source("simulations/SR-sim-functions.R") # local

source_obj <- ls()

#-------------------------------------
# Experimental Design
#-------------------------------------

set.seed(20240518) 

# design factors
design_factors <- list(
  mu = c(0, 0.2, 0.5, 0.8),
  tau = c(0, 0.1, 0.2, 0.4),
  k = c(10, 30, 60, 100),
  cor_mu = c(0.2, 0.4, 0.8),
  cor_sd = c(0.05),
  cut_vals = .025,
  weights = c(1, 0.8, 0.5, 0.25, 0.125, .05),
  id_start = 0L,
  paste_ids = TRUE
)

design_factors_3PSM <- 
  list(
    mu = c(0, 0.2, 0.5, 0.8),
    tau = c(0, 0.1, 0.2, 0.4),
    k = c(10, 30, 60, 100),
    cor_mu = c(0.2, 0.4, 0.8),
    cor_sd = c(0.05),
    cut_vals = .025,
    weights = c(1, 0.8, 0.5, 0.25, 0.125, .05),
    id_start = 0L,
    paste_ids = TRUE
  )

suppressWarnings(
params_3PSM <- 
  purrr::cross_df(design_factors) |> 
  mutate(
    omega = tau / sqrt(2),
    k_multiplier = 1 / weights
  )
)

design_factors_4PSM <- list(
  mu = c(0, 0.2, 0.5, 0.8),
  tau = c(0, 0.1, 0.2, 0.4),
  k = c(10, 30, 60, 100),
  cor_mu = c(0.2, 0.4, 0.8),
  cor_sd = c(0.05),
  weight_labs = c("one", "two", "three", "four", "five", "six"),
  id_start = 0L,
  paste_ids = TRUE
)

weights_tab <- tibble(
  weight_labs = c("one", "two", "three", "four", "five", "six"),
  weights = list(c(1, .5), c(.8, .4), c(.5, .25),
                 c(.25, .125), c(.125, .0625), c(.05, .025)),
  cut_vals = list(c(.025, .50))
)

suppressWarnings(
params_4PSM <- 
  purrr::cross_df(design_factors_4PSM) |> 
  left_join(weights_tab, by = "weight_labs") |> 
  select(-weight_labs) |> 
  mutate(
    omega = tau / sqrt(2),
    k_multiplier = 1 / purrr::map_dbl(weights, function(x) x[2])
  )
)

params <- rbind(params_3PSM, params_4PSM)
nrow(params) # 2304

# combine into a design set
params <- 
  params |>
  dplyr::mutate(
    iterations = reps,
    seed = 20240518 + 1:n()
  )

params
range(params$seed)

n_ES_data <- readRDS("empirical-sample-sizes.rds") # condor
# n_ES_data <- readRDS("simulations/empirical-sample-sizes.rds") # local
n_ES_sim <- n_ES_empirical(n_ES_data)

#--------------------------------------------------------
# run simulations in parallel
#--------------------------------------------------------

tm <- system.time(
  results <- plyr::mdply(params[batch, ], .fun = runSim,
                         n_ES_sim = n_ES_empirical(n_ES_data),
                         k_stop = 5)
)

tm

#--------------------------------------------------------
# Save results and details
#--------------------------------------------------------

session_info <- sessionInfo()
run_date <- date()

results_name <- paste0("Simulation-results-batch",batch,".Rdata")
save(results, tm, session_info, run_date, file = results_name)
