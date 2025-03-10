rm(list=ls())

#----------------------------
# data-generating function 
#----------------------------

r_mvt_items <- function(n, p, icc, df) {
  V_mat <- icc + diag(1 - icc, nrow = p)
  X <- mvtnorm::rmvt(n = n, sigma = V_mat, df = df)
  colnames(X) <- LETTERS[1:p]
  X
}

small_sample <- r_mvt_items(n = 8, p = 3, icc = 0.7, df = 5)
# small_sample

#----------------------------
# Estimation function
#----------------------------

alpha_CI <- function(dat, coverage = .95) {
  V <- cov(dat)
  p <- ncol(dat)
  n <- nrow(dat)
  A <- p / (p - 1) * (1 - sum(diag(V)) / sum(V))
  B <- log(1 - A) / 2
  SE_B <- sqrt(p / (2 * n * (p - 1)))
  z <- qnorm((1 - coverage) / 2)
  CI_B <- B + data.frame(L = -1, U = 1) * SE_B * z
  1 - exp(2 * CI_B)
}

# alpha_CI(small_sample)

#----------------------------
# simulation driver
#----------------------------

simulate_alpha <- function(reps, n, p, alpha, df, ..., seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  icc <- alpha / (p - alpha * (p - 1))
  purrr::rerun(reps, {
    dat <- r_mvt_items(n = n, p = p, icc = icc, df = df)
    alpha_CI(dat)
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(cover = L < alpha & alpha < U) %>%
    dplyr::summarise(coverage = mean(cover))
}

source_obj <- ls()

#----------------------------
# experimental design
#----------------------------
library(tidyr)
library(dplyr)

alpha_design <- 
  crossing(
    alpha = seq(0.5, 0.9, 0.1),
    n = seq(20, 100, 20), 
    p = c(3, 9, 12)
  ) %>%
  mutate(
    reps = 1000,
    df = 10,
    seed = row_number()
  )

#----------------------------
# simulate in parallel
#----------------------------
# library(parallel)
# 
# system.time(
#   results <- 
#     alpha_design %>%
#     mutate(
#       res = do.call(mcmapply, c(FUN = simulate_alpha, alpha_design, mc.cores = 18L))
#     ) %>%
#     unnest(res)
# )
# 
# results

# library(future)
# library(furrr)
# options(mc.cores = 18)
# plan(multicore)
# options(future.rng.onMisuse = "ignore")
# 
# system.time(
#   results <-
#     alpha_design %>%
#     mutate(
#       res = future_pmap(alpha_design, .f = simulate_alpha)
#     ) %>%
#     unnest(res)
# )
# 
# results

doParallel::registerDoParallel(cl = 18)

system.time(
  results <- plyr::mdply(alpha_design, .fun = simulate_alpha, .parallel = TRUE)
)

results