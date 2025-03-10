#-------------------------------------
#        Data Generating Model  
#-------------------------------------

# Simulate empirical distribution of sample size and number of effect size (Functional)

## James: n_ES_empirical() is a functional (a function that returns a function). 
## You give it a dataset with primary study sample sizes and numbers of effect sizes per study. 
## It returns a function that generates random samples from the dataset — random study characteristics.

n_ES_empirical <- function(dat) {
  d <- dat
  function(k) d[sample(NROW(dat), size = k, replace = TRUE), ]
  
}

n_ES_param <- function(N_mean, J_mean, min_N = 20L) {
  function(k) {
    data.frame(
      N = min_N + rpois(k, lambda = N_mean - min_N), # N = 20 + 2 * rpois(k, (N_mean - 20) / 2)
      J = 1L + rpois(k, lambda = J_mean - 1L) # J = 2 + rpois(k, J_mean - 2)
    )
  }
}

# Reparameterization
## Man: r_corr gives you one correlation among outcomes for each study
## return k correlations for k studies (assume equal-correlation within study)
## use cor_mu = 0.4, 0.8 and sd_mu = 0.05 following Rodgers & Pustejovsky

r_corr <- function(k, cor_mu, cor_sd) {
  if (cor_sd > 0 & cor_mu > 0) {
    alpha <- cor_mu * ((cor_mu * (1 - cor_mu) / cor_sd^2 ) - 1)
    beta <-  (1 - cor_mu) * ((cor_mu * (1 - cor_mu) / cor_sd^2 ) - 1)
    rbeta(n = k, shape1 = alpha, shape2 = beta)
  } else {
    rep(cor_mu, k)
  }
  
}


# Generate primary data for study k (raw data approach)
## n is the sample size per primary study; n_ES is number of outcomes per individual or number of effects per study
## r_k is correlation among outcomes (equal correlation within study), used for the sigma_matrix

r_study <- function(delta_k, # study-level effect size(s)
                    r_k, # correlation between outcomes
                    n, # sample size
                    n_ES) { # number of es per individual 
  
  # create treatment and control groups
  # divide by 2 and round for uneven sample sizes
  n_tx <- round((n / 2), 0)
  
  # Create Sigma matrix for outcomes, assuming equal correlation
  Sigma_mat <- r_k + diag(rep(1 - r_k, n_ES), nrow = n_ES)
  Sigma_mat <- as.matrix(Sigma_mat)
  
  # Generate treatment/control summary statistics
  y_tx <- mvtnorm::rmvnorm(n_tx, mean = rep(delta_k, length.out = n_ES), sigma = Sigma_mat) # outcomes for each individual in trt grp
  ybar_tx <- colMeans(y_tx)
  var_tx <- diag(cov(y_tx))
  
  y_ctl <- mvtnorm::rmvnorm(n - n_tx, mean = rep(0, n_ES), sigma = Sigma_mat)
  ybar_ctl <- colMeans(y_ctl)
  var_ctl <- diag(cov(y_ctl))
  
  # Calculate deltas and other necessary statistics
  sd_pooled <- ((n_tx - 1) * var_tx + (n - n_tx - 1) * var_ctl) / (n - 2)
  J <- 1 - (3 / (4 * (n - 2) - 1))
  d <- (ybar_tx - ybar_ctl) / sqrt(sd_pooled)
  g <- d * J
  var_g <- (4 / n + d^2 / (2 * (n - 2))) * J^2
  sd_g <- sqrt(var_g)
  Va <- 4 / n # Man: do we need Va and sda? might matter for PET-PEESE?
  sda <- sqrt(Va)
  df <- n - 2
  t_i <- d / sqrt(Va) 
  p_onesided <- pt(t_i, df = df, lower.tail = FALSE) 
  
  # Keep only necessary values
  data.frame(n, n_ES, g, var_g, sd_g, Va, sda, t_i, p_onesided, esid = 1:n_ES)
}


# Censoring functions ---------------------------------------------
# Man: returns a function of p-vals that give you the weight that an effect with a particular p-value was selected
# cut_vals and weights? weights: .3 or .6 in Banks et al., 2016
step_fun <- function(cut_vals = .025, weights = 1) {
  
  if (length(cut_vals) != length(weights)) stop("cut_vals and weights must be the same length, doofus!")
  
  wt_vec <- c(1, weights) # the weight for effects with p-values in each interval being selected
  cut_vals_full <- c(0, cut_vals, 1) # p-value intervals where cut_vals are the thresholds or cutpoints
  
  function(pvals) {
    
    pval_buckets <- cut(pvals, cut_vals_full, labels = FALSE, include.lowest = TRUE)
    wt_vec[pval_buckets] 
    
  }
}


# Generate meta analytic data 

r_meta <- function(mu, # effect size parameter (true mean of effect sizes/overall average effect size)
                   tau, # between-study heterogeneity
                   omega, # within-study heterogeneity
                   k, # number of primary studies
                   cor_mu, # average correlation between outcomes
                   cor_sd, # sd correlation between outcomes
                   n_ES_sim, # distribution of sample sizes and number of outcomes
                   cut_vals,
                   weights,
                   k_multiplier = 2, # related to the severity of the selection, if the selection is 50%, then double k_multiplier to get the target (number of studies)
                   id_start = 0L,
                   paste_ids = TRUE) {   

  # Calculate buffer for number of studies to generate (number of studies before selection)
  k_star <- round(k * k_multiplier)
  
  # Sample correlation values for matrix
  r_k <- r_corr(k = k_star, cor_mu, cor_sd) 

  # Simulate delta_k for each individual within each study. 
  delta_k <- rnorm(k_star, mean = mu, sd = tau)
  
  # Sample a sample-size and n_ES combo 
  n_select <- as.data.frame(n_ES_sim(k_star))
  
  # generate study design parameters
  study_parms <- data.frame(delta_k = delta_k, r_k = r_k, n = n_select$n, n_ES = n_select$n_ES)
  
  # add within-study heterogeneity 
  delta_jk <- purrr::pmap(study_parms, .f = \(delta_k, n_ES, ...) rnorm(n_ES, mean = delta_k, sd = omega))
  study_parms$delta_k <- delta_jk # a vector if multiple outcomes/effects per study
  
  # generate the studies
  studies <- purrr::pmap_dfr(study_parms, .f = r_study, .id = "studyid")
  studies$studyid <- id_start + as.integer(studies$studyid)
  
  # apply censoring process
   K_star <- nrow(studies) # number of studies generated, why not using length(p_onesided)?
    censor_fun <- step_fun(cut_vals = cut_vals, weights = weights)
    p_sel <- censor_fun(studies$p_onesided) # probability of being selected, 1 for sig effects, weight for non-sig effects
    observed <- runif(K_star) <= p_sel # supposed to be less than right? Man: I think so
  
  studies <- studies[observed,,drop=FALSE]

  # Count how many unique studies are kept
  n_studies <- length(unique(studies$studyid))

  # recurse to get desired number of studies per meta-analytic dataset
  if (n_studies < k) {

    more_studies <- r_meta(
      mu = mu,
      tau = tau, omega = omega,
      k = k - n_studies,
      cor_mu = cor_mu, cor_sd = cor_sd,
      n_ES_sim = n_ES_sim,
      cut_vals = cut_vals, weights = weights,
      k_multiplier = k_multiplier,
      id_start = id_start + k_star,
      paste_ids = FALSE
    )

    if (nrow(studies) > 0L) {
      studies <- rbind(studies, more_studies)
    } else {
      studies <- more_studies
    }
  }

  # keep only desired number of studies
  if (paste_ids) studies$esid <- paste(studies$studyid, studies$esid, sep = "-")
  studies$studyid <- factor(studies$studyid)
  studies <- droplevels(studies[studies$studyid %in% (levels(studies$studyid)[1:k]),])
  
  return(studies)
  
  # var_avg?
  # studies <- 
  #   studies %>% 
  #   group_by(studyid) %>% 
  #   mutate(var_avg = mean(var_g)) %>% 
  #   ungroup() %>% 
  #   dplyr::select(studyid, esid, everything())
  # 
  # return(as.data.frame(studies))
  
}


#-------------------------------------
#        Estimation functions
#-------------------------------------
estimate_meta_uni <- function(data,
                              moderator = ~ 1,  # change to sd or var for pet peese
                              W = FALSE, 
                              method, # this is the adjustment method
                              meta_method, # this is the estimation method for meta-analysis
                              returnMod = FALSE) {
  require(metafor)
  
  if (W) W <- 1 / data$Va else W <- NULL
  
  optimizers <- c("nlminb","nloptr","Rvmmin","BFGS")
  mod <- "Non-converged"
  i <- 1L
  
  while (!inherits(mod, "rma.uni") & i <= 4L) {
    if (!is.null(W)) {
      mod <- tryCatch(
        rma.uni(
          yi = g,
          vi = Va,
          weights = W,
          mods = moderator,
          method = meta_method,
          data = data,
          control = list(optimizer=optimizers[i])
        ),
        error = function(e) "Non-converged"
      )
    } else {
      mod <- tryCatch(
        rma.uni(
          yi = g,
          vi = Va,
          mods = moderator,
          method = meta_method,
          data = data,
          control = list(optimizer=optimizers[i])
        ),
        error = function(e) "Non-converged"
      )
    }
    i <- i + 1L
  }
  
  if (inherits(mod, "rma.uni")) { # rma.uni for EK, PET, PEESE, WAAP, WILS
    mod <- robust(mod, cluster = studyid, clubSandwich = TRUE) # always apply RVE
    
    if (returnMod == TRUE) {
      res <- mod
    } else {
      res <- data.frame(
        param = "beta.intrcpt",
        est = as.numeric(mod$beta[1]),
        se = mod$se[1],
        lo = mod$ci.lb[1],
        hi = mod$ci.ub[1],
        p_val = mod$pval[1],
        R_conv = NA,
        method = method
      )
    }
  } else {
    res <- NULL
  }
  
  return(res)
  
}


estimate_meta <- function(data,
                          subset = NULL,
                          rho = .8,
                          rand = NULL,
                          moderator = ~ 1,  # change to sd or var for pet peese
                          V_mat = NULL,
                          W = FALSE, # true for correlated fixed effects model
                          method,
                          modified = TRUE,
                          returnMod = FALSE) {
  require(metafor)
  
  if (!is.null(subset)) {data <- subset(data, subset)}
  
  if (modified == TRUE & is.null(V_mat)) {
    V_mat <- vcalc(Va, cluster = studyid, obs = esid, data = data, rho = rho)
  } else if (modified == FALSE & is.null(V_mat)) {
    V_mat <- vcalc(var_g, cluster = studyid, obs = esid, data = data, rho = rho)
  }
  
  if(W == TRUE) W <- solve(V_mat) else {W <- NULL}
  
  optimizers <- c("nlminb","nloptr","Rvmmin","BFGS")
  mod <- "Non-converged"
  i <- 1L
  
  while (!inherits(mod, "rma.mv") & i <= 4L) {
    if (!is.null(W)) {
      mod <- tryCatch(
        rma.mv(
          yi = g,
          V = V_mat,
          W = W, # did not know how to remove W argument
          random = rand,
          mods = moderator,
          data = data,
          sparse = TRUE,
          control = list(optimizer=optimizers[i])
        ),
        error = function(e) "Non-converged"
      )
    } else {
      mod <- tryCatch(
        rma.mv(
          yi = g,
          V = V_mat,
          random = rand,
          mods = moderator,
          data = data,
          sparse = TRUE,
          control = list(optimizer=optimizers[i])
        ),
        error = function(e) "Non-converged"
      )
    }
    i <- i + 1L
  }
  
  if (inherits(mod, "rma.mv")) { # rma.mv for CHE, PET, PEESE, WAAP, WILS
    mod <- robust(mod, cluster = studyid, clubSandwich = TRUE) # always apply RVE
    
    if (returnMod == TRUE) {
      res <- mod
    } else {
      res <- data.frame(
        param = "beta.intrcpt",
        est = as.numeric(mod$beta[1]),
        se = mod$se[1],
        lo = mod$ci.lb[1],
        hi = mod$ci.ub[1],
        p_val = mod$pval[1],
        R_conv = NA,
        method = method
      )
    }
  } else {
    res <- NULL
  }
  
  return(res)
  
}


estimate_mv_wils <- function(data, modified = TRUE, k_stop = 5) {
  
  require(metafor)
  require(clubSandwich)
  
  if (modified) {
    method <- "mv_wils"
  } else  {
    data$sda <- data$sd_g
    data$Va <- data$var_g
    method <- "mv_wils_um"
  }
  
  k <- length(unique(data$studyid)) # k is the number of studies
  N_ES <- nrow(data)
  drop <- 1L
  i <- 1L
  
  while (i <= 5L && drop > 0 && drop < N_ES && k > k_stop) { 
    # correlated random effects with fixed weights
    mv_uwls <- estimate_meta(data = data, W = TRUE, 
                             rand = ~ 1 | studyid / esid,
                             method = "CHE_ISCW", modified = modified,
                             returnMod = TRUE)
    
    if (!is.null(mv_uwls)) {
      data$resid <- (data$g - as.numeric(mv_uwls$beta)) / data$sda
      data <- data[order(data$resid), ]
      omegasq <- mv_uwls$sigma2[2]
      tausq <- mv_uwls$sigma2[1]
    } else { # remove the rand part
      mv_uwls <- estimate_meta(data = data, W = TRUE, 
                               method = "FE", modified = modified,
                               returnMod = TRUE)
      data$resid <- (data$g - as.numeric(mv_uwls$beta)) / data$sda
      data <- data[order(data$resid), ]
      if (sum(table(data$studyid) > 1L) >= 3) {
        res_CHE <- estimate_meta(data = data,
                                 rand = ~ 1 | studyid / esid,
                                 modified = modified,
                                 method = "CHE",
                                 returnMod = TRUE)
        tausq <- res_CHE$sigma2[1]
        omegasq <- res_CHE$sigma2[2]
      } else {
        res_CHE <- estimate_meta(data = data,
                                 rand = ~ 1 | studyid,
                                 modified = modified,
                                 method = "CHE",
                                 returnMod = TRUE)
        tausq <- res_CHE$sigma2[1]
        omegasq <- 0
      }
    }

    # calculate ESS
    zz <- (1.96 * data$sda - as.numeric(mv_uwls$beta)) / sqrt(data$Va + omegasq + tausq) # eq.4
    data$Esig <- 1 - pnorm(zz) # power: the probability an effect is sig given its sampling dist
    data$SS <- 0
    data$SS[data$t_i > 1.96] <- 1
    data$ESS <- data$SS - data$Esig
    Nsize <- nrow(data)
    drop <- round(sum(data$ESS))
    
    if (drop > 0 & drop < Nsize) {
      keep_index <- Nsize - drop
      keep_dat <- data[1:keep_index, ]
      # make sure there remain at least k_stop studies
      if (length(unique(keep_dat$studyid)) >= k_stop) {
        data <- keep_dat
      } else {
        keep_studyids <- unique(data$studyid)[1:k_stop] # select the first k_stop
        data <- subset(data, studyid %in% keep_studyids)
      }
    } else {
      data <- data
    }
    
    k <- length(unique(data$studyid))
    N_ES <- nrow(data)
    i <- i + 1L
  }
  
  mv_uwls <- estimate_meta(data = data,
                           W = TRUE,
                           rand = ~ 1 | studyid / esid,
                           method = method,
                           modified = modified,
                           returnMod = FALSE)
  
  return(mv_uwls)
  
}


estimate_wils <- function(data, modified = TRUE, k_stop = 5) {
  
  require(metafor)
  require(clubSandwich)
  
  if (modified) {
    method <- "uni_wils"
  } else  {
    data$sda <- data$sd_g
    data$Va <- data$var_g
    method <- "uni_wils_um"
  }
  
  k <- length(unique(data$studyid)) # number of studies
  N_ES <- nrow(data) # number of effect sizes in total
  drop <- 1L
  i <- 1L
  
  while (i <= 5L && drop > 0 && drop < N_ES && k > k_stop) {
    # need to correct the way resid is calculated 
    UWLS_fit <- rma.uni(yi = g, vi = Va, weights = 1 / Va,
                        data = data, method =  c("REML","DL"))
    data$resid <- (data$g - as.numeric(UWLS_fit$beta)) / data$sda
    data <- data[order(data$resid), ]
    tausq <- UWLS_fit$tau2
    
    # calculate ESS: Eq.4 & 5 in Stanley 2022
    zz <- (1.96 * data$sda - as.numeric(UWLS_fit$beta)) / sqrt(data$Va + tausq) 
    data$Esig <- 1 - pnorm(zz) # Esig_i probability the effect is sig
    data$SS <- 0
    data$SS[data$t_i > 1.96] <- 1 # observed significance
    data$ESS <- data$SS - data$Esig # excessive significance
    Nsize <- nrow(data)
    drop <- round(sum(data$ESS))
    
    if (drop > 0 & drop < Nsize) {
      keep_index <- Nsize - drop
      keep_dat <- data[1:keep_index, ]
      # keep_dat$studyid <- droplevels(keep_dat$studyid)
      
      # make sure there are at least k_stop studies in the dataset
      if (length(unique(keep_dat$studyid)) >= k_stop) {
        data <- keep_dat
      } else {
        keep_studyids <- unique(data$studyid)[1:k_stop] # select the first k_stop
        data <- subset(data, studyid %in% keep_studyids)
        # data$studyid <- droplevels(data$studyid)
      }
    } else {
      data <- data
    }
    
    k <- length(unique(data$studyid))
    N_ES <- nrow(data)
    
    i <- i + 1L
    
  }
  
  # Point est same as WLS est, SE differs, does not matter for this paper
  res_wils <- estimate_meta_uni(data = data, W = TRUE, 
                                method = method, meta_method = c("REML","DL"))

  return(res_wils)
  
}


# regression-based univariate
fit_reg_uni <- function(data, modified = TRUE, k_stop = 5) {
   
  require(metafor)
  require(dplyr) # why do i need this? for filter, could just use subset
  
  if (modified == FALSE) {
    data$Va <- data$var_g
    data$sda <- data$sd_g
  }
  
  # FE 
  FE <- estimate_meta_uni(data, method = "uni_FE", meta_method = "FE")
  
  # RE
  RE <- estimate_meta_uni(data, method = "uni_RE", meta_method = "REML")
  
  # RE-ISW
  RE_ISW <- estimate_meta_uni(data, W = TRUE, method = "uni_RE_ISW", meta_method = "REML")
  
  # DL estimated tau^2 for EK method (differ slightly from the formula in the paper)
  fit_RE_DL <- rma.uni(yi = g, vi = Va, weights = 1 / Va, data = data, method = "DL")
  tausq_hat <- fit_RE_DL$tau2
  tau_hat <- sqrt(tausq_hat)
  
  # PET-PEESE
  PET <- estimate_meta_uni(data, W = TRUE, moderator = ~ I(sda), 
                           method = "uni_pet", meta_method = c("REML","DL"))
  PEESE <- estimate_meta_uni(data, W = TRUE, moderator = ~ I(Va), 
                             method = "uni_peese", meta_method = c("REML","DL"))
  PET_PEESE <- if (PET$p_val < 0.1 & PET$est > 0) PEESE else PET
  PET_PEESE$method <- "uni_pet_peese"
  
  # EK with RVE
  beta0_hat <- PET_PEESE$est
  a_est <- (beta0_hat^2 - 1.96^2 * tausq_hat) / ((1.96 + 1.96) * beta0_hat) # dc if it is 1.96?
  kink <- if (beta0_hat <= 1.96 * tau_hat) 0 else a_est
  data$kr_mod <- (data$sda - kink) * (data$sda > kink)
  
  EK <- estimate_meta_uni(data = data, W = TRUE, moderator = ~ I(kr_mod),
                          method = "uni_ek", meta_method = c("REML","DL"))
  
  # WAAP
  # use the FE (or WLS) as the proxy for true effect
  fit_FE <- rma.uni(yi = g, vi = Va, data = data, method = "FE")
  waap_cutoff <- as.numeric(abs(fit_FE$beta / 2.8)) # power for a fixed true effect
  data$waapInd <- data$sda <= waap_cutoff
  waap_data <- subset(data, waapInd == TRUE) 
  if (length(unique(waap_data$studyid)) >= 2) {
    # Note FE SE is different from WLS by RMSE
    WAAP <- estimate_meta_uni(data = waap_data, W = TRUE, 
                              method = "uni_waap", meta_method = c("REML","DL")) # FE 
  } else {
    WAAP <- estimate_meta_uni(data = data, W = TRUE, 
                              method = "uni_waap", meta_method = c("REML","DL"))
  }
  
  # WILS
  # WILS <- estimate_wils(data = data, modified = modified, k_stop = k_stop)
  
  return(rbind(FE, RE, RE_ISW, PET, PEESE, PET_PEESE, EK, WAAP))
  
}


# CHE, PET, PEESE, WAAP, WILS
fit_reg_mv <- function(data, rho = 0.8, modified = TRUE, k_stop = 5) {
  
  require(metafor)
  require(dplyr)
  
  if (modified == TRUE) {
    V_mat <- vcalc(Va, cluster = studyid, obs = esid, data = data, rho = rho)
  } else {
    V_mat <- vcalc(var_g, cluster = studyid, obs = esid, data = data, rho = rho)
  }
  
  # W <- solve(V_mat)
  
  # CHE
  if (sum(table(data$studyid) > 1L) >= 3) {
    res_CHE <- estimate_meta(data = data,
                             rand = ~ 1 | studyid / esid,
                             V_mat = V_mat,
                             modified = modified,
                             method = "CHE")
  } else {
    res_CHE <- estimate_meta(data = data,
                             rand = ~ 1 | studyid,
                             V_mat = V_mat,
                             modified = modified,
                             method = "CHE")
  }
  
  # CHE-ISCW
  res_CHE_ISCW <- estimate_meta(data = data, 
                           rand = ~ 1 | studyid / esid,
                           V_mat = V_mat, 
                           W = TRUE,
                           method = "CHE-ISCW", modified = modified)
  if (is.null(res_CHE_ISCW)) {
    res_CHE_ISCW <- estimate_meta(data = data, 
                             rand = ~ 1 | studyid,
                             V_mat = V_mat, 
                             W = TRUE,
                             method = "CHE-ISCW", modified = modified)
  }
  
  # multivariate PET-PEESE
  pet <- estimate_meta(data = data, V_mat = V_mat, W = TRUE, moderator = ~ sda, 
                       rand = ~ 1 | studyid / esid, 
                       modified = modified, method = "mv_pet")
  if (is.null(pet)) { # make sure pet and peese converge
    pet <- estimate_meta(data = data, V_mat = V_mat, W = TRUE, moderator = ~ sda, 
                         rand = ~ 1 | studyid, 
                         modified = modified, method = "mv_pet")
  }
  
  peese <- estimate_meta(data = data, V_mat = V_mat, W = TRUE, moderator = ~ Va, 
                         rand = ~ 1 | studyid / esid, 
                         modified = modified, method = "mv_peese")
  if (is.null(peese)) {
    peese <- estimate_meta(data = data, V_mat = V_mat, W = TRUE, moderator = ~ Va, 
                           rand = ~ 1 | studyid, 
                           modified = modified, method = "mv_peese")
  }
  
  pet_peese <- if (pet$p_val < 0.1 & pet$est > 0) peese else pet # .1 for one sided, .05 for two sided
  pet_peese$method <- "mv_pet_peese"
  res_pet_peese <- bind_rows(pet, peese, pet_peese)
  
  # multivariate EK
  mod_CHE_ISCW <- estimate_meta(data = data, 
                           rand = ~ 1 | studyid / esid,
                           V_mat = V_mat, 
                           W = TRUE,
                           method = "CHE_ISCW", modified = modified,
                           returnMod = TRUE)
  if (!is.null(mod_CHE_ISCW)) {
    tausq_hat <- mod_CHE_ISCW$sigma2[1] # updated on 6/13
    omegasq_hat <- mod_CHE_ISCW$sigma2[2]
    tau_hat <- sqrt(tausq_hat)
    omega_hat <- sqrt(omegasq_hat)
  } else {
    mod_CHE_ISCW <- estimate_meta(data = data, 
                             rand = ~ 1 | studyid,
                             V_mat = V_mat, 
                             W = TRUE,
                             method = "CHE_ISCW", modified = modified,
                             returnMod = TRUE)
    tausq_hat <- mod_CHE_ISCW$sigma2
    omegasq_hat <- 0
    tau_hat <- sqrt(tausq_hat)
    omega_hat <- sqrt(omegasq_hat)
  }
  beta0_hat <- pet_peese$est
  a_est <- (beta0_hat^2 - 1.96^2 * (tausq_hat + omegasq_hat)) / ((1.96 + 1.96) * beta0_hat)
  kink <- if (beta0_hat <= 1.96 * sqrt(tau_hat^2 + omega_hat^2)) 0 else a_est
  data$kr_mod <- (data$sda - kink) * (data$sda > kink)
  
  EK <- estimate_meta(data = data, V_mat = V_mat, 
                      W = TRUE, moderator = ~ I(kr_mod),
                      rand = ~ 1 | studyid / esid,
                      modified = modified, method = "mv_ek")
  
  # multivariate waap is multivariate version of WLS for ES with adequate power
  ## use CHE_ISCW as initial, determine power based on that
  waap_cutoff <- as.numeric(abs(res_CHE_ISCW$est / 2.8))
  data$waapInd <- data$sda <= waap_cutoff
  waap_data <- subset(data, waapInd == TRUE)
  waap_j <- length(unique(waap_data$studyid))
  if (waap_j >= 2) { # if it does not converge, it does not converge.
    # should not specify V_mat bc it should be a subset of the original V_mat
    res_waap <- estimate_meta(data = waap_data, 
                              W = TRUE,
                              rand = ~ 1 | studyid / esid,
                              method = "mv_waap",
                              modified = modified)
  } else {
    res_waap <- res_CHE_ISCW
    res_waap$method <- "mv_waap"
  }
  
  # multivariate WILS
  # res_wils <- estimate_mv_wils(data = data, modified = modified, k_stop = k_stop)
  
  return(rbind(res_CHE, res_CHE_ISCW, res_pet_peese, EK, res_waap))
  
}


# trim-and-fill test 
fit_TF <- function(data) {
  
  model <- rma.uni(yi = g, vi = Va, data = data, method = "FE")
  TF_se_R0 <- trimfill(model, side = "left", estimator = "R0")
  
  with(TF_se_R0, data.frame(param = "beta.intrcpt",
                            est = as.numeric(b),
                            se = se,
                            lo = ci.lb,
                            hi = ci.ub,
                            p_val = pval, 
                            R_conv = NA,
                            method = "TF_R0_left"))  
}

# p-uniform and p-uniform star
fit_puniform <- function(data, otherwise = NULL) {
  
  safe_puniform <- purrr::possibly(puniform, otherwise = otherwise)
  # Man: not sure about side?
  fit <- safe_puniform(tobs = data$t_i, n1i = data$n/2, n2i = data$n/2, side = "right", plot = FALSE)
  
  if (is.null(fit)) {
    NULL
  } else {
    data.frame(param = "beta.intrcpt",
               est = fit$est,
               se = NA,
               lo = fit$ci.lb,
               hi = fit$ci.ub,
               p_val = fit$pval.0, 
               R_conv = NA,
               method = "p-uniform") # pval for testing publication bias
  }
  
}


fit_punistar <- function(data, otherwise = NULL) {
  
  # possibly is not working.
  # safe_punistar <- purrr::possibly(puni_star, otherwise = otherwise)
  # fit <- safe_punistar(ni = data$n, tobs = data$t_i, side = "right", plot = FALSE)
  
  fit <- tryCatch(puniform::puni_star(tobs = data$t_i, n1i = data$n/2, n2i = data$n/2, side = "right"), 
                   error = function(e) NULL)
  
  if (is.null(fit)) {
    NULL
  } else {
    data.frame(param = "beta.intrcpt",
               est = fit$est,
               se = NA,
               lo = fit$ci.lb,
               hi = fit$ci.ub,
               p_val = fit$pval.0, 
               R_conv = NA,
               method = "p-unistar") # pval for testing publication bias
  }
  
}


# 3PSM: Ignoring dependence
fit_sel <- function(data, type = "stepfun", 
                     alternative = "greater", 
                     steps = .025, otherwise = NULL) {

  model <- tryCatch(metafor::rma.uni(yi = g, vi = Va, data = data, method = c("REML","DL")), 
                  error = function(e) NULL)
  
  if (is.null(model)) {
    model <- rma.uni(yi = g, vi = Va, weights = 1/Va, data = data, method = c("REML","DL"))
  }
  
  safe_selmodel <- purrr::possibly(selmodel, otherwise = otherwise)
  
  fit <- safe_selmodel(model, type = type, 
                       alternative = alternative, 
                       steps = steps)
  
  if (length(steps) == 1) method <- "3PSM" else method <- "4PSM"
  
  if (is.null(fit)) {
    NULL
  } else {
    data.frame(param = "beta.intrcpt",
               est = as.numeric(fit$beta),
               se = fit$se, 
               lo = fit$ci.lb,
               hi = fit$ci.ub,
               p_val = fit$pval,
               R_conv = NA,
               method = method)
  }
}


estimate_models <- function(data, 
                            rho = 0.8, 
                            k_stop = 5, # for WILS
                            smooth_vi = TRUE, 
                            ignoring = TRUE, # ignore dependence, set this way if need to run uni and mv separately
                            dependence = TRUE,
                            modified = TRUE) {
  
  suppressPackageStartupMessages(require(metafor, quietly = TRUE, warn.conflicts = FALSE))
  require(clubSandwich, quietly = TRUE, warn.conflicts = FALSE)
  require(dplyr, quietly = TRUE, warn.conflicts = FALSE)
  require(purrr, quietly = TRUE, warn.conflicts = FALSE)
  
  res <- data.frame()
  
  if (ignoring) {
    if (modified == FALSE) {data$Va <- data$var_g}
    # trim and fill
    est_TF <- fit_TF(data)
    # p-uniform
    est_puniform <- fit_puniform(data)
    est_punistar <- fit_punistar(data)
    # 3PSM and 4PSM
    est_3PSM <- fit_sel(data, steps = .025) # set steps as the cut_vals??
    est_4PSM <- fit_sel(data, steps = c(.025, .50))
    # reg uni
    est_reg_uni <- fit_reg_uni(data, modified = modified, k_stop = k_stop)
    
    res_ignoring <- 
      rbind(
        est_TF,
        est_puniform,
        est_punistar,
        est_3PSM,
        est_4PSM,
        est_reg_uni
      ) 
    
    res <- rbind(res, res_ignoring)
  }
  
  if (dependence) {
    res_dep <- fit_reg_mv(data, rho = rho, modified = modified, k_stop = k_stop)
    res <- rbind(res, res_dep)
  }
  
  return(res)
  
}


#-------------------------------------
#        Simulation drive
#-------------------------------------


calc_performance <- function(res_dat) {
  
  require(dplyr, quietly = TRUE, warn.conflicts = FALSE)
  require(tidyr, quietly = TRUE, warn.conflicts = FALSE)
  require(simhelpers, quietly = TRUE, warn.conflicts = FALSE)
  
  abs <- 
    res_dat %>%
    group_by(method) %>% 
    group_modify(~ calc_absolute(.x, 
                                 estimates = est, 
                                 true_param = mu,
                                 criteria = c("bias", "variance", "mse", "rmse"))) %>%
    ungroup() %>% 
    dplyr::rename(K = K_absolute)
  
  ci_cov <- 
    res_dat %>%
    group_by(method) %>% 
    group_modify(~ calc_coverage(.x, 
                                 lower_bound = lo, 
                                 upper_bound = hi, 
                                 true_param = mu,
                                 criteria = c("coverage","width"))) %>%
    ungroup() %>% 
    dplyr::rename(K = K_coverage)
  
  rej <- 
    res_dat %>% 
    group_by(method) %>% 
    group_modify(~ calc_rejection(.x, p_values = p_val)) %>% 
    dplyr::rename(K = K_rejection)
  
  left_join(abs, ci_cov, by = c("method", "K")) %>% 
    left_join(rej, by = c("method", "K")) %>% 
    dplyr::rename(n_converged = K)
  
}


#-------------------------------------
#        Simulation drive
#-------------------------------------

runSim <- function(iterations,
                   k, mu, tau, omega, cor_mu, cor_sd, 
                   n_ES_sim, 
                   cut_vals, weights,
                   k_stop = 5,
                   ignoring = TRUE, dependence = TRUE,
                   modified = TRUE, returnData = FALSE,
                   seed = NULL, ...) {
  
  suppressPackageStartupMessages(require(purrr))
  suppressPackageStartupMessages(require(dplyr))
  
  if (!is.null(seed)) set.seed(seed)
  
  if (returnData) {
    results <- replicate(iterations, {
      dat <- r_meta(k = k, mu = mu, tau = tau, omega = omega, 
                    cor_mu = cor_mu, cor_sd = cor_sd,
                    n_ES_sim = n_ES_sim, 
                    cut_vals = unlist(cut_vals), weights = unlist(weights))
    }, simplify = FALSE)
    
    res_dat <-
      results %>%
      bind_rows(.id = "rep")
    
    return(res_dat)
    
  } else {
    
    results <- replicate(iterations, {
      dat <- r_meta(k = k, mu = mu, tau = tau, omega = omega, 
                    cor_mu = cor_mu, cor_sd = cor_sd,
                    n_ES_sim = n_ES_sim, 
                    cut_vals = unlist(cut_vals), weights = unlist(weights))
      res <- estimate_models(data = dat, k_stop = k_stop, ignoring = ignoring,
                             dependence = dependence, modified = modified)
      
    }, simplify = FALSE)
    
    res_dat <-
      results %>%
      bind_rows(.id = "rep") %>% 
      mutate(
        mu = mu,
        tau = tau,
        tausq = tau^2,
        omega = omega,
        omegasq = omega^2
      )
    
    calc_performance(res_dat)
    
  }
  
}
