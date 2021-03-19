args <- commandArgs(trailingOnly = TRUE)

library(dplyr, quietly = TRUE, warn.conflicts = FALSE)
library(tidyr)
library(purrr)
library(skellam)
library(metafor)
library(clubSandwich)
library(simhelpers)


#-------------------------------------
# Parse command line arguments
#-------------------------------------
extract_int_arg <- function(arg_string, arg_name, default) {
  res <- 
    arg_string %>%
    stringr::str_extract(paste(arg_name, "[0-9]+")) %>%
    stringr::str_sub(stringr::str_length(arg_name) + 2, -1) %>%
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

source("RVE-sim-functions.R")

source_obj <- ls()

#-------------------------------------
# Experimental Design
#-------------------------------------

set.seed(20210316) 

# design factors
design_factors <- list(
  K = c(10, 30),
  J_max = 5,
  theta = seq(-0.8, 0.8, 0.2),
  tau = c(0.1, 0.3),
  omega = c(0.1, 0.2),
  mu_A_shape = 2,
  mu_A_scale = 7,
  m_mean = 7,
  n_mean = 7,
  phi = c(0, 0.2, 0.4),
  metric = c("LRR", "WC-SMD", "TAU")
)

# combine into a design set
params <- 
  purrr::cross_df(design_factors) %>%
  dplyr::mutate(
    iterations = reps,
    seed = round(runif(1) * 2^30) + 1:n() + (batch - 1) * nrow(.)
  )

range(params$seed)
params

#--------------------------------------------------------
# run simulations in parallel
#--------------------------------------------------------

tm <- system.time(
  results <- plyr::mdply(params[batch, ], .fun = run_sim)
)

tm

#--------------------------------------------------------
# Save results and details
#--------------------------------------------------------

session_info <- sessionInfo()
run_date <- date()

results_name <- paste0("Simulation-results-",batch,".Rdata")
save(results, tm, session_info, run_date, file = results_name)
