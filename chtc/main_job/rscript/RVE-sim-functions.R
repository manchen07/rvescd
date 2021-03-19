#------------------------------------------------------
# Data Generating Model
#------------------------------------------------------

# simulate AR(1) poisson series
r_pois_AR1 <- function(mu, phi) {
  
  len <- length(mu)
  mu <- rep(mu, length.out = len)
  ar_vec <- ifelse(mu[-1] - phi * mu[-len] >= 0, phi, mu[-1] / mu[-len])
  lambda <- pmax(0, mu[-1] - ar_vec * mu[-len])
  
  # binomial thinning
  y <- vector(mode = "numeric", length = len)
  y[1] <- rpois(1, lambda = mu[1])
  for (i in 1:(len-1)) 
    y[i+1] <- rbinom(1, size = y[i], prob = ar_vec[i]) + rpois(1, lambda = lambda[i])
  
  y
}


# calculate mu_b
f_func <- function(mu_A, nap) function(x) (1 - skellam::pskellam(0, mu_A, x) + 0.5 * skellam::dskellam(0, mu_A, x) - nap)

mu_b_nap <- function(mu_A, nap, tol = 10^-5) {
  
  nap <- min(max(nap, tol), 1 - tol)
  
  if (nap == 0.5) {
    mu_A
  } else {
    interval <- if (nap > 1/ 2L) c(0, mu_A) else c(mu_A, (sqrt(mu_A) + sqrt(-log(nap)))^2)
    
    f <- f_func(mu_A, nap)
    mu_b_interval <- f(interval)
    
    if (all(mu_b_interval > 0)) {
      interval[1]
    } else if (all(mu_b_interval < 0)) {
      interval[2]
    } else {
      uniroot(f, interval = interval)$root    
    }
  }
}


# simulate LRRs
r_LRR <- function(lambda, x, mu_A, phi) {
  
  mu_B <- exp(lambda) * mu_A
  mu_vec <- mu_A * (1 - x) + mu_B * x
  y <- r_pois_AR1(mu_vec, phi = phi)
  
  n <- table(x)
  M <- pmax(tapply(y, x, mean), 1 / (2 * n))
  V <- pmax(tapply(y, x, var), 1 / n^3)
  
  LRR1 <- as.numeric(diff(log(M)))
  BC <- log(M) + V / (2 * n * M^2)
  LRR2 <- BC[[2]] - BC[[1]]
  V_LRR <- sum(V / (n * M^2))

  data.frame(LRR1 = LRR1, LRR2 = LRR2, V_LRR = V_LRR)
}


# simulate WC-SMD

r_SMD <- function(delta, x, mu_A, phi) {
  
  mu_B <- mu_A + delta * sqrt(mu_A)
  mu_vec <- mu_A * (1 - x) + mu_B * x
  y <- r_pois_AR1(mu_vec, phi = phi)
  
  n <- table(x)
  M <- tapply(y, x, mean)
  V <- pmax(tapply(y, x, var), 1 / n^3)
  
  if (V[[1]] > 0) {
    d <- (M[[2]] - M[[1]]) / sqrt(V[[1]])
    J <- (1 - 3 / (4 * n[[1]] - 5))
    g <- J * d
    V_g <- J^2 * (1 / n[[1]] + V[[2]] / (n[[2]] * V[[1]]) + d^2 / (2 * n[[1]] -2))
  } else {
    d <- g <- V_g <- NA
  }
  
  data.frame(d = d, g = g, V_g = V_g)
}


# simulate TAU

r_NAP <- function(nap, x, mu_A, phi) {
  
  mu_B <- mu_b_nap(mu_A, nap)
  mu_vec <- mu_A * (1 - x) + mu_B * x
  y <- r_pois_AR1(mu_vec, phi = phi)
  dat <- split(y, x)
  m <- length(dat[[1]])
  n <- length(dat[[2]])
  Q_mat <- sapply(dat[[2]], function(j) (j < dat[[1]]) + 0.5 * (j == dat[[1]]))
  NAP <- mean(Q_mat)
  
  Q1 <- sum(rowSums(Q_mat)^2) / (m * n^2)
  Q2 <- sum(colSums(Q_mat)^2) / (m^2 * n)
  
  trunc <- 0.5 / (m * n)
  NAP_trunc <- min(max(NAP, trunc), 1 - trunc)
  
  X <- sum(Q_mat^2)/(m * n)
  V_un <- (NAP_trunc - (m + n - 1) * NAP_trunc^2 + n * Q1 + m * Q2 - 2 * X) / ((m - 1) * (n - 1))
  
  V_HM <- (NAP_trunc * (1 - NAP_trunc) + (n-1) * (Q1 - NAP_trunc^2) + (m-1) * (Q2 - NAP_trunc^2)) / (m * n)
  h <- (m + n) / 2 - 1
  V_Newc <- NAP_trunc * (1 - NAP_trunc) * (1 + h * (1 - NAP_trunc) / (2 - NAP_trunc) + h * NAP_trunc / (1 + NAP_trunc)) / (m * n)
  
  data.frame(NAP = NAP, V_un = V_un, V_HM = V_HM, V_Newc = V_Newc)
}


r_TAU <- function(TAU, x, mu_A, phi) {
  
  nap <- (TAU + 1) / 2
  res_nap <- r_NAP(nap = nap, x = x, mu_A = mu_A, phi = phi)
  TAU_nap <- as.numeric(2 * res_nap$NAP - 1)
  V_TAU <- 4 * res_nap[-1]

  data.frame(TAU = TAU_nap, V_TAU_un = V_TAU$V_un, V_TAU_HM = V_TAU$V_HM, V_TAU_Newc = V_TAU$V_Newc)
  
}


r_es <- function(theta_k, J, omega, mu_A, m_mean, n_mean, phi, metric, return_theta_jk = FALSE) {
  
  case <- c(1:J)
  m <- 3 + rpois(J, m_mean - 3) # at least three series
  n <- 3 + rpois(J, n_mean - 3)
  x <- purrr::map2(m, n, ~ c(rep(0, .x), rep(1, .y)))
  
  # generate case-level effect sizes
  if (metric == "Tau") {
    
    if (omega^2 <= (1 + theta_k) * (1 - theta_k)) {
      
      alpha_k <- (theta_k + 1)^2 * (1 - theta_k) / (2 * omega^2) - (theta_k + 1) / 2
      beta_k <- alpha_k * (2 / (theta_k + 1) - 1)
      theta_jk <- 2 * rbeta(n = J, shape1 = alpha_k, shape2 = beta_k) - 1
      
    } else {
      theta_jk <- rep(theta_k, times = J)
      
    }
    
  } else {
    theta_jk <- rnorm(n = J, mean = theta_k, sd = omega)
    
  }
  
  if (return_theta_jk) return(theta_jk)
  
  # generate data for each effect size metric
  if (metric == "LRR") {
    dat <- purrr::map2_dfr(theta_jk, x, ~ r_LRR(lambda = .x, x = .y, mu_A = mu_A, phi = phi))
  } else if (metric == "WC-SMD") {
    dat <- purrr::map2_dfr(theta_jk, x, ~ r_SMD(delta = .x, x = .y, mu_A = mu_A, phi = phi))
  } else {
    dat <- purrr::map2_dfr(theta_jk, x, ~ r_TAU(TAU = .x, x = .y, mu_A = mu_A, phi = phi))
  } 
  
  return(data.frame(case = case, mu_A = mu_A, m = m, n = n, dat))
}


# simulate meta-analytic data
r_meta <- function(K, J_max, theta, tau, omega, mu_A_shape, mu_A_scale,
                   m_mean, n_mean, phi, metric, return_study_data = FALSE) {
  
  if (metric == "Tau") {
    alpha <- (theta + 1)^2 * (1 - theta) / (2 * tau^2) - (theta + 1) / 2
    beta <- alpha * (2 / (theta + 1) - 1)
    theta_k <- 2 * rbeta(n = K, shape1 = alpha, shape2 = beta) - 1
  } else {
    theta_k <- rnorm(n = K, mean = theta, sd = tau)
  }
  
  study_data <- 
    data.frame(
      theta_k = theta_k, 
      J = sample(1:J_max, size = K, replace = TRUE),
      mu_A = rgamma(K, shape = mu_A_shape, scale = mu_A_scale)
    )
  
  if (return_study_data) return(study_data)
  
  purrr::pmap_df(.l = study_data, .f = r_es, omega = omega, m = m_mean, n = n_mean, 
                 phi = phi, metric = metric, .id = "study")
} 



#------------------------------------------------------
# Model-fitting/estimation functions
#------------------------------------------------------

estimate_meta <- function(dat, y, v) {
  
  require(metafor, quietly = TRUE, warn.conflicts = FALSE)
  require(clubSandwich, quietly = TRUE, warn.conflicts = FALSE)
  
  names(dat)[which(names(dat) == y)] <- "y"
  names(dat)[which(names(dat) == v)] <- "v"
  
  rma_fit <- rma.mv(yi = y, V = v, 
                    random = ~ 1 | study / case, data = dat,
                    test = "t", method = "REML", sparse = TRUE)
  rma_rve <- coef_test(rma_fit, vcov = "CR2", cluster= dat$study)
  rma_rve_ci <- conf_int(rma_fit, vcov = "CR2", cluster = dat$study)
  ci_rve_l <- rma_rve_ci$CI_L
  ci_rve_u <- rma_rve_ci$CI_U
  
  ols_fit <- lm(y ~ 1, data = dat)
  ols_pval <- summary(ols_fit)$coefficients[, 4]
  ols_rve <- coef_test(ols_fit, vcov = "CR2", cluster = dat$study)
  ols_rve_ci <- conf_int(ols_fit, vcov = "CR2", cluster = dat$study)
  ci_ols_rve_l <- ols_rve_ci$CI_L
  ci_ols_rve_u <- ols_rve_ci$CI_U
  
  data.frame(y = y,
             v = v,
             tausq_est = rma_fit$sigma2[1],
             omegasq_est = rma_fit$sigma2[2],
             rma_est = rma_rve$beta,
             rma_mod_l = rma_fit$ci.lb,
             rma_mod_u = rma_fit$ci.ub,
             rma_rve_l = ci_rve_l,
             rma_rve_u = ci_rve_u,
             ols_est = ols_rve$beta,
             ols_rve_l = ci_ols_rve_l,
             ols_rve_u = ci_ols_rve_u,
             rma_pval = rma_fit$pval,
             rve_pval = rma_rve$p_Satt,
             ols_pval = ols_pval,
             ols_rve_pval = ols_rve$p_Satt
  )
}

#------------------------------------------------------
# Calculate performance measures
#------------------------------------------------------

calc_perf <- function(res_dat) {
  
  require(dplyr, quietly = TRUE, warn.conflicts = FALSE)
  require(tidyr, quietly = TRUE, warn.conflicts = FALSE)
  require(simhelpers, quietly = TRUE, warn.conflicts = FALSE)
  
  perf_abs <- 
    res_dat %>% 
    select(y, v, theta, rma_est, ols_est) %>%
    gather(key = "est_method", value = "es_est", rma_est, ols_est) %>%
    group_by(y, v, est_method) %>%
    do(calc_absolute(., estimates = es_est, true_param = theta,
                                 perfm_criteria = c("bias", "variance"))) %>%
    gather(key = "k", value = "val", bias:var_mcse) %>%
    unite(est_method, k, est_method) %>%
    spread(key = est_method, value = val)
  
  perf_rbs <- 
    res_dat %>% 
    select(y, v, tausq, omegasq, tausq_est, omegasq_est) %>%
    gather(key = "het_type", value = "value", tausq_est, omegasq_est) %>%
    mutate(true_het_var = if_else(het_type == "tausq_est", tausq, omegasq)) %>% 
    group_by(y, v, het_type) %>%
    do(calc_relative(., estimates = value, true_param = true_het_var,
                     perfm_criteria = c("relative bias"))) %>% 
    gather(key = "k", value = "val", rel_bias:rel_bias_mcse) %>%
    unite(het_type, k, het_type) %>%
    spread(key = het_type, value = val)
  
  perf_ci <-
    res_dat %>%
    select(y, v, theta, ends_with("_l"), ends_with("_u")) %>%
    group_by(y,v) %>%
    mutate(rep = row_number()) %>%
    gather(key = "est_method", value = "val", rma_mod_l:ols_rve_u) %>%
    separate(est_method, c("est_method", "b"), sep = -2) %>%
    mutate(b = substr(b, 2, 2)) %>%
    spread(key = b, value = val) %>%
    group_by(y, v, est_method) %>%
    do(calc_coverage(., lower_bound = l, upper_bound = u, true_param = theta)) %>%
    gather(key = "k", value = "val", coverage:width_mcse) %>%
    unite(est_method, k, est_method) %>%
    spread(key = est_method, value = val)

  perf_rej_rate <- 
    res_dat %>% 
    select(y, v, ends_with("_pval")) %>%
    gather(key = "pval_method", value = "pval", rma_pval:ols_rve_pval) %>%
    group_by(y, v, pval_method) %>%
    do(calc_rejection(., p_values = pval)) %>% 
    gather(key = "k", value = "val", rej_rate:rej_rate_mcse) %>%
    unite(pval_method, k, pval_method) %>%
    spread(key = pval_method, value = val)

  performance <- 
    perf_abs %>% 
    left_join(perf_rbs, by = c("y", "v", "K")) %>% 
    left_join(perf_ci, by = c("y", "v", "K")) %>% 
    left_join(perf_rej_rate, by = c("y", "v", "K")) %>% 
    rename(n_converged = K)
  
  return(performance)
}



#------------------------------------------------------
# Simulation Driver
#------------------------------------------------------

run_sim <- function(iterations, K, J_max, theta, tau, omega, 
                    mu_A_shape, mu_A_scale, m_mean, n_mean, phi, metric, 
                    seed = NULL, summarize_results = TRUE) {
  
  if (!is.null(seed)) set.seed(seed)
  
  require(purrr, quietly = TRUE, warn.conflicts = FALSE)
  require(dplyr, quietly = TRUE, warn.conflicts = FALSE)
  
  possibly_meta <- possibly(.f = estimate_meta, otherwise = NULL)
  
  if (metric == "LRR") {
    
    results <- replicate(iterations, {
      dat <- r_meta(K, J_max, theta, tau, omega, mu_A_shape, mu_A_scale, 
                    m_mean, n_mean, phi, metric = "LRR")
      est_res <- rbind(possibly_meta(dat, y = "LRR1", v = "V_LRR"),
                       possibly_meta(dat, y = "LRR2", v = "V_LRR"))
      return(est_res)
      
    }, simplify = FALSE)
    
  } else if (metric == "WC-SMD") {
    
    results <- replicate(iterations, {
      dat <- r_meta(K, J_max, theta, tau, omega, mu_A_shape, mu_A_scale, 
                    m_mean, n_mean, phi, metric = "WC-SMD")
      est_res <- possibly_meta(dat, y = "g", v = "V_g")
      return(est_res)
      
    }, simplify = FALSE)
    
  } else {
    
    results <- replicate(iterations, {
      dat <- r_meta(K, J_max, theta, tau, omega, mu_A_shape, mu_A_scale, 
                    m_mean, n_mean, phi, metric = "TAU")
      est_res <- possibly_meta(dat, y = "TAU", v = "V_TAU_un")
      return(est_res)
      
    }, simplify = FALSE)
    
  }
  
  res_dat <- 
    do.call(rbind, args = results) %>%
    mutate(theta = theta, tausq = tau^2, omegasq = omega^2)
  
  if (summarize_results) {
    performance <- calc_perf(res_dat = res_dat)  
    return(performance)
  } else {
    return(res_dat)
  }
  
}

