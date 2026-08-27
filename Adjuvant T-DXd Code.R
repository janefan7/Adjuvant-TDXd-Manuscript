# Load packages
library(msm)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ggridges)
library(GGally)
# install IMIS from CRAN archived packages using devtools
# devtools::install_version("IMIS", version = "0.1", 
#                           repos = "http://cran.us.r-project.org")
library(IMIS)

# Global model parameters
tx_window_months <- 10
age_init <- 50
delta <- 1/12

# Strategy names
strategy_names <- c("TDM1", "TDXd")

# Death states
death_states <- c("OC_Death", "BC_Death", "ILD_Death")

# Names of recurrence-free (RF) states after ILD discontinuation by month
post_ild_state_names <- paste0("RF_postILD_m", 1:tx_window_months)

# All state names
state_names <- c("RF_on", post_ild_state_names, "RF_offtx_complete", "LRR", "DR", 
                 "OC_Death", "BC_Death", "ILD_Death")
n_states <- length(state_names)

# 2023 CDC life table for women
lt <- read.csv("data/life_table_2023.csv") %>%
  transmute(age = as.numeric(stringr::str_extract(age, "^\\d+")),
            qx  = as.numeric(qx)) %>%
  arrange(age) %>%
  mutate(qx = pmin(qx, 1 - 1e-8), haz_mo = -log(1 - qx) / 12)

# Function to obtain other-cause mortality (background mortality) from life table
haz_OC <- function(age_years) {
  a <- floor(age_years)
  lt$haz_mo[match(a, lt$age)]
}

#haz_OC <- function(age_years) {
#  a <- floor(age_years)
#  idx <- which(lt$age == a)
#  if (length(idx) == 0) {
#    return(NA)
#  }
#  lt$haz_mo[idx]
#}
  
# ILD inputs from DB-05
p_any_ILD_TDXd_base <- 0.096
p_any_ILD_TDM1_base <- 0.016
any_ild_TDXd <- 77
any_ild_TDM1 <- 13
dis_TDXd_nonfatal <- 59
fatal_TDXd <- 2
dis_TDM1 <- 5
n_TDM1 <- 801

## ILD inputs in deterministic analyses, base case ##
# Overall fatal ILD probability among all T-DXd-assigned patients
p_fatal_ILD_all_TDXd_base <- 0.01

# Fatal ILD among any-grade ILD cases, T-DXd
prop_fatal_ILD_TDXd <- p_fatal_ILD_all_TDXd_base / p_any_ILD_TDXd_base

# DB-05 discontinuation probability among nonfatal ILD cases
prop_dis_given_nonfatal_ILD_TDXd <- dis_TDXd_nonfatal / (any_ild_TDXd - fatal_TDXd)

# Nonfatal ILD discontinuation among any-grade ILD cases, T-DXd
prop_dis_ILD_TDXd <- (1 - prop_fatal_ILD_TDXd) * prop_dis_given_nonfatal_ILD_TDXd

# T-DM1 fatal ILD assumed 0 in base case
prop_fatal_ILD_TDM1 <- 0

# Discontinuation among any-grade ILD cases, T-DM1
prop_dis_ILD_TDM1 <- dis_TDM1 / any_ild_TDM1

# No fatal ILD for T-DM1 in DB-05
prop_dis_given_nonfatal_ILD_TDM1 <- prop_dis_ILD_TDM1

## Helpers
# Clip probabilities between 0 and 1
clip01 <- function(p, eps = 1e-6) {
  pmin(pmax(p, eps), 1 - eps)
}

# Compute log-likelihood contribution for one calibration target
ll_logitnorm <- function(p_hat, p_obs, se_logit_obs) {
  p_hat <- clip01(p_hat)
  p_obs <- clip01(p_obs)
  dnorm(qlogis(p_obs),
        mean = qlogis(p_hat),
        sd   = se_logit_obs,
        log  = TRUE
  )
}

# Obtain cycle-specific HR for RF to DR transition
get_hr_RF_DR <- function(t, arm, params) {
  
  # HR for T-DM1 is always 1 (reference)
  if (arm == "TDM1") {
    return(1)
  }
  
  # Full treatment effect
  hr_initial <- params$hr_RF_DR_TDXd_vs_TDM1
  
  # Evaluate time-varying HR at the midpoint of each monthly cycle
  time_years <- (t - 0.5) * params$delta
  
  # Treatment effect persists through year given by hr_full_effect_years
  if (time_years <= params$hr_full_effect_years) {
    hr_t <- hr_initial
  }
  # Treatment effect has fully waned by hr_wane_end_year
  else if (time_years >= params$hr_wane_end_year) {
    hr_t <- 1
  }
  else {
    # Linear waning from hr_full_effect_years to hr_wane_end_year
    time_years_since_waning <- time_years - params$hr_full_effect_years
    linear_slope <- (1 - hr_initial) / 
      (params$hr_wane_end_year - params$hr_full_effect_years)
    hr_t <- hr_initial + time_years_since_waning * linear_slope
  }
  hr_t
}

# Build cycle-specific transition intensity matrix Q
build_Q_matrix <- function(t, arm, params) {
  
  # Other-cause mortality (from life table) at current age
  age <- params$age_init + (t - 0.5) * params$delta
  lam_oc_mo <- params$haz_OC(age)
  
  # Indicator for 10-month treatment window 
  in_tx_window <- t <= params$tx_window_months
  
  # Cycle-specific HR for RF to DR
  hr_t <- get_hr_RF_DR(t = t, arm = arm, params = params)
  
  # Apply arm-specific ILD discontinuation hazard if in treatment window
  lam_ild_dis <- if (in_tx_window && arm == "TDXd") {
    params$lam_ild_dis_TDXd_mo
  } 
  else if (in_tx_window && arm == "TDM1") {
    params$lam_ild_dis_TDM1_mo
  } 
  else {
    0
  }
  
  # Apply arm-specific fatal ILD hazard if in treatment window
  lam_ild_fatal <- if (in_tx_window && arm == "TDXd") {
    params$lam_ild_fatal_TDXd_mo
  } 
  else if (in_tx_window && arm == "TDM1") {
    params$lam_ild_fatal_TDM1_mo
  } 
  else {
    0
  }
  
  # Specify transition intensity matrix
  Q <- matrix(0, nrow = n_states, ncol = n_states, 
              dimnames = list(state_names, state_names))
  
  # From RF on-treatment states
  Q["RF_on", "LRR"]       <- params$haz_RF_LRR
  Q["RF_on", "DR"]        <- hr_t * params$haz_RF_DR_TDM1
  Q["RF_on", "OC_Death"]  <- lam_oc_mo
  Q["RF_on", "ILD_Death"] <- lam_ild_fatal
  if (in_tx_window){
    Q["RF_on", paste0("RF_postILD_m", t)] <- lam_ild_dis 
  }
  
  # Fraction of benefit retained in the case of ILD discontinuation
  for (st in post_ild_state_names) {
    month_discontinued <- as.integer(gsub("RF_postILD_m", "", st))
    if (arm == "TDXd") {
      hr_postILD <- hr_t ^ params$benefit_frac_postILD[month_discontinued]
    } 
    else {
      hr_postILD <- 1
    }
    # From RF post-ILD discontinuation states
    Q[st, "LRR"]      <- params$haz_RF_LRR
    Q[st, "DR"]       <- hr_postILD * params$haz_RF_DR_TDM1
    Q[st, "OC_Death"] <- lam_oc_mo
  }
  
  # From RF post-treatment completion states
  Q["RF_offtx_complete", "LRR"]      <- params$haz_RF_LRR
  Q["RF_offtx_complete", "DR"]       <- hr_t * params$haz_RF_DR_TDM1
  Q["RF_offtx_complete", "OC_Death"] <- lam_oc_mo
  
  # From LRR
  Q["LRR", "DR"]       <- params$haz_LRR_DR
  Q["LRR", "OC_Death"] <- lam_oc_mo
  
  # From DR
  Q["DR", "BC_Death"] <- params$haz_DR_BCd
  Q["DR", "OC_Death"] <- lam_oc_mo
  
  diag(Q) <- -rowSums(Q)
  Q
}

# Run Markov cohort model for 10-year OS
run_markov_model <- function(params, compute_drfi = FALSE) {
  
  # Any-grade ILD hazard
  lam_ild_any_TDM1_mo <- -log(1 - params$p_any_ILD_TDM1) / params$tx_window_months
  lam_ild_any_TDXd_mo <- -log(1 - params$p_any_ILD_TDXd) / params$tx_window_months
  
  # T-DM1 fatal and nonfatal ILD hazards
  lam_ild_fatal_TDM1_mo <- 0 # no fatal ILD for T-DM1 reported in DB-05
  lam_ild_nonfatal_TDM1_mo <- lam_ild_any_TDM1_mo
  
  # T-DM1 ILD discontinuation hazard (discontinuation is a subset of nonfatal ILD)
  lam_ild_dis_TDM1_mo <- lam_ild_nonfatal_TDM1_mo * params$prop_dis_given_nonfatal_ILD_TDM1
  
  # T-DXd fatal and nonfatal ILD hazards
  lam_ild_fatal_TDXd_mo <- lam_ild_any_TDXd_mo * params$prop_fatal_ILD_TDXd
  lam_ild_nonfatal_TDXd_mo <- lam_ild_any_TDXd_mo * (1 - params$prop_fatal_ILD_TDXd)
  
  # T-DXd ILD discontinuation hazard
  lam_ild_dis_TDXd_mo <- lam_ild_nonfatal_TDXd_mo * params$prop_dis_given_nonfatal_ILD_TDXd
  
  # Store monthly ILD hazards
  params$lam_ild_dis_TDM1_mo   <- lam_ild_dis_TDM1_mo
  params$lam_ild_dis_TDXd_mo   <- lam_ild_dis_TDXd_mo
  params$lam_ild_fatal_TDM1_mo <- lam_ild_fatal_TDM1_mo
  params$lam_ild_fatal_TDXd_mo <- lam_ild_fatal_TDXd_mo
  
  # Initial state vector
  v_m_init <- rep(0, n_states)
  names(v_m_init) <- state_names
  v_m_init["RF_on"] <- 1
  
  # Initialize output
  out_os <- rep(0, length(strategy_names))
  names(out_os) <- strategy_names
  out_idfs <- rep(0, length(strategy_names))
  names(out_idfs) <- strategy_names
  out_drfi <- rep(0, length(strategy_names))
  names(out_drfi) <- strategy_names
  out_bc_death <- rep(0, length(strategy_names))
  names(out_bc_death) <- strategy_names
  out_ild_death <- rep(0, length(strategy_names))
  names(out_ild_death) <- strategy_names
  out_oc_death <- rep(0, length(strategy_names))
  names(out_oc_death) <- strategy_names
  
  # All RF states
  rf_states <- c("RF_on", post_ild_state_names, "RF_offtx_complete")
  
  # At-risk states for first DR
  dr_risk_states <- c(rf_states, "LRR")
  
  # Number of cycles
  n_cycles <- as.integer(round(params$horizon_years / params$delta))
  
  # Obtain cohort trace and clinical outcomes by arm
  for (arm in strategy_names) {
    
    # Initialize cohort trace
    m_trace <- matrix(0, nrow = n_cycles + 1, ncol = n_states, 
                      dimnames = list(0:n_cycles, state_names))
    m_trace[1, ] <- v_m_init
    
    # Initialize cumulative cause-specific hazard for first DR
    # to calculate DRFI among patients at risk of first DR if compute_DRFI = T
    cum_haz_DR <- 0
    
    for (t in seq_len(n_cycles)) {
      
      # Obtain transition intensity matrix
      Q_t <- build_Q_matrix(t = t, arm = arm, params = params)
      
      # Distribution halfway through the cycle for the DR risk set for a better
      # approximation
      if (compute_drfi) {
        P_mid_t <- msm::MatrixExp(Q_t / 2)
        m_mid_t <- as.numeric(m_trace[t, ] %*% P_mid_t)
        names(m_mid_t) <- state_names
        
        # Cause-specific DR hazard midway through cycle among patients
        # in DR at-risk states
        at_risk_mid <- sum(m_mid_t[dr_risk_states])
        
        # Compute cause-specific DR hazard midway through cycle
        if (at_risk_mid > 0) {
          haz_DR_mid <- sum(m_mid_t[dr_risk_states] * 
                              Q_t[dr_risk_states, "DR"]) / at_risk_mid
          # Update cumulative DR hazard
          cum_haz_DR <- cum_haz_DR + haz_DR_mid
        }
      }
      
      # Advance cohort trace
      P_t <- msm::MatrixExp(Q_t)
      m_trace[t + 1, ] <- as.numeric(m_trace[t, ] %*% P_t)
      
      # At the end of the treatment window remaining RF_on patients move to 
      # RF_offtx_complete
      if (t == params$tx_window_months) {
        m_trace[t + 1, "RF_offtx_complete"] <- m_trace[t + 1, "RF_offtx_complete"] + 
          m_trace[t + 1, "RF_on"]
        m_trace[t + 1, "RF_on"] <- 0
      }
    }
    
    # Obtain cohort trace at 10 years
    final_trace <- m_trace[n_cycles + 1, ]
    
    # Obtain arm-specific clinical outcomes at 10 years
    out_os[arm]   <- 1 - sum(final_trace[death_states])
    out_idfs[arm] <- sum(final_trace[rf_states])
    out_drfi[arm] <- if (compute_drfi) {
      exp(-cum_haz_DR)
    } else {
      NA_real_
    }
    out_bc_death[arm] <- sum(final_trace["BC_Death"])
    out_ild_death[arm] <- sum(final_trace["ILD_Death"])
    out_oc_death[arm] <- sum(final_trace["OC_Death"])
  }
  
  # Return outputs
  list(OS = out_os, iDFS = out_idfs, DRFI = out_drfi, 
       BC_death = out_bc_death, ILD_death = out_ild_death, OC_death = out_oc_death)
}

## Setup for IMIS calibration
# Specify DB-05 3-year OS, iDFS, and DRFI endpoints as the main targets
targets_endpoints <- data.frame(
  arm = c("TDM1", "TDM1", "TDM1", "TDXd", "TDXd", "TDXd"),
  outcome = c("iDFS", "DRFI", "OS", "iDFS", "DRFI", "OS"),
  time = c(3, 3, 3, 3, 3, 3),
  value = c(0.837, 0.861, 0.957, 0.924, 0.939, 0.974),
  lower = c(0.802, 0.825, 0.935, 0.897, 0.914, 0.958),
  upper = c(0.867, 0.891, 0.972, 0.944, 0.957, 0.984)
) %>%
  mutate(se_logit = (qlogis(upper) - qlogis(lower)) / (2 * 1.96),
         arm = factor(arm, levels = c("TDM1", "TDXd")),
         outcome = factor(outcome, levels = c("iDFS", "DRFI", "OS")))

# Specify LRR share informed by DB-05 as a separate target
target_LRR_share <- data.frame(value = 0.10, lower = 0.08, upper = 0.12) %>%
  mutate(se_logit = (qlogis(upper) - qlogis(lower)) / (2 * 1.96))

# Predict calibration endpoints
predict_all_endpoints_calib <- function(par_vec,
                                        base_params,
                                        arm,
                                        horizon_years_calib = 3) {
  # Base parameters
  params <- base_params
  params$horizon_years <- horizon_years_calib
  
  # Hazards
  params$haz_RF_DR_TDM1 <- par_vec["haz_RF_DR_base"]
  params$haz_RF_LRR <- (par_vec["LRR_share"] / (1 - par_vec["LRR_share"])) * 
    params$haz_RF_DR_TDM1 
  params$haz_LRR_DR <- par_vec["LRR_DR_multiplier"] * params$haz_RF_DR_TDM1 
  params$haz_DR_BCd <- par_vec["haz_DR_BCd_base"]
  params$hr_RF_DR_TDXd_vs_TDM1 <- par_vec["hr_RF_DR_TDXd_vs_TDM1"]
  
  # Obtain outputs
  out <- run_markov_model(params, compute_drfi = TRUE)
  list(OS = out$OS[arm], iDFS = out$iDFS[arm], DRFI = out$DRFI[arm])
}

# Parameters to be calibrated
v_param_names <- c("haz_RF_DR_base",
                   "LRR_share",
                   "LRR_DR_multiplier",
                   "haz_DR_BCd_base",
                   "hr_RF_DR_TDXd_vs_TDM1")
n_param <- length(v_param_names)

# Specify lower and upper bounds of independent uniform prior distributions
v_lb <- c(haz_RF_DR_base  = 1e-6,
          LRR_share = 0.01,
          LRR_DR_multiplier = 1,
          haz_DR_BCd_base = 1e-4,
          hr_RF_DR_TDXd_vs_TDM1 = 0.20)

v_ub <- c(haz_RF_DR_base  = 0.01,
          LRR_share = 0.30,
          LRR_DR_multiplier = 2,
          haz_DR_BCd_base = 0.15,
          hr_RF_DR_TDXd_vs_TDM1 = 1.00)

# ILD assumptions for calibration
prop_fatal_ILD_TDXd_calib <- fatal_TDXd / any_ild_TDXd # fatal ILD among any-grade ILD
p_fatal_ILD_all_TDXd_calib <-
  p_any_ILD_TDXd_base * prop_fatal_ILD_TDXd_calib # fatal ILD among all T-DXd-assigned patients
prop_dis_given_nonfatal_ILD_TDXd_calib <-
  dis_TDXd_nonfatal / (any_ild_TDXd - fatal_TDXd) # discontinuation among nonfatal ILD
prop_dis_ILD_TDXd_calib <-
  (1 - prop_fatal_ILD_TDXd_calib) * 
  prop_dis_given_nonfatal_ILD_TDXd_calib # discontinuation among any-grade ILD

# Specify all parameters needed for calibration
l_params_calib <- list(
  # General model settings
  delta            = delta,
  age_init         = age_init,
  tx_window_months = tx_window_months,
  horizon_years    = 3,
  
  # Treatment-effect persistence assumptions
  hr_full_effect_years   = 3,
  hr_wane_end_year       = 10,
  
  # Placeholders overwritten during calibration
  haz_RF_LRR            = NA_real_,
  haz_RF_DR_TDM1        = NA_real_,
  haz_LRR_DR            = NA_real_,
  hr_RF_DR_TDXd_vs_TDM1 = NA_real_,
  haz_DR_BCd            = NA_real_,
  
  # Other-cause mortality
  haz_OC = haz_OC,
  
  # ILD parameters
  p_any_ILD_TDM1 = p_any_ILD_TDM1_base,
  p_any_ILD_TDXd = p_any_ILD_TDXd_base,
  
  p_fatal_ILD_all_TDXd = p_fatal_ILD_all_TDXd_calib,
  prop_fatal_ILD_TDM1  = 0,
  prop_fatal_ILD_TDXd  = prop_fatal_ILD_TDXd_calib,
  
  #prop_dis_ILD_TDM1 = prop_dis_ILD_TDM1,
  #prop_dis_ILD_TDXd = prop_dis_ILD_TDXd_calib,
  
  prop_dis_given_nonfatal_ILD_TDM1 = prop_dis_given_nonfatal_ILD_TDM1,
  prop_dis_given_nonfatal_ILD_TDXd = prop_dis_given_nonfatal_ILD_TDXd_calib,
  
  # A-priori-specified fractional benefit after ILD discontinuation
  benefit_frac_postILD = c(
    0.28, 0.49, 0.64, 0.75, 0.83,
    0.89, 0.94, 0.97, 0.99, 1.00
  )
)

# Latin hypercube sampling
sample_prior <- function(n_samp) {
  u <- lhs::randomLHS(n = n_samp, k = n_param)
  m <- matrix(NA_real_, nrow = n_samp, ncol = n_param,
              dimnames = list(NULL, v_param_names))
  for (i in seq_len(n_param)) {
    m[, i] <- qunif(u[, i], min = v_lb[i], max = v_ub[i])
  }
  m
}

# Joint log prior density
calc_log_prior <- function(v_params) {
  if (is.null(dim(v_params))) v_params <- t(v_params)
  colnames(v_params) <- v_param_names
  lp <- rep(0, nrow(v_params))
  for (i in seq_len(n_param)) {
    lp <- lp + dunif(v_params[, i], min = v_lb[i], max = v_ub[i], log = TRUE)
  }
  lp
}

# Prior density
calc_prior <- function(v_params) {
  exp(calc_log_prior(v_params))
}

# Joint log-likelihood
calc_log_lik <- function(v_params) {
  
  if (is.null(dim(v_params))) v_params <- t(v_params)
  colnames(v_params) <- v_param_names
  
  n_samp <- nrow(v_params)
  
  # Proposals outside the bounded prior support have zero posterior density
  llik <- rep(-Inf, n_samp)
  
  # IMIS Gaussian proposals can fall outside the bounded uniform priors;
  # prevent the Markov model from running for those proposals
  valid <- apply(v_params, 1, function(x) {
    all(is.finite(x)) &&
      all(x >= v_lb[v_param_names]) &&
      all(x <= v_ub[v_param_names])
  })
  
  # Valid draws start with log-likelihood = 0
  llik[valid] <- 0
  
  for (j in which(valid)) {
    
    pred_TDM1 <- predict_all_endpoints_calib(
      par_vec = v_params[j, ],
      base_params = l_params_calib,
      arm = "TDM1"
    )
    
    pred_TDXd <- predict_all_endpoints_calib(
      par_vec = v_params[j, ],
      base_params = l_params_calib,
      arm = "TDXd"
    )
    
    pred <- list(TDM1 = pred_TDM1, TDXd = pred_TDXd)
    
    for (k in seq_len(nrow(targets_endpoints))) {
      
      arm_k <- as.character(targets_endpoints$arm[k])
      outcome_k <- as.character(targets_endpoints$outcome[k])
      
      # Likelihood contribution from main targets (DB-05 3-year endpoints)
      llik[j] <- llik[j] +
        ll_logitnorm(
          p_hat  = pred[[arm_k]][[outcome_k]],
          p_obs  = targets_endpoints$value[k],
          se_logit_obs = targets_endpoints$se_logit[k]
        )
    }
    
    # Likelihood contribution from LRR share target
    llik[j] <- llik[j] +
      ll_logitnorm(
        p_hat  = v_params[j, "LRR_share"],
        p_obs  = target_LRR_share$value,
        se_logit_obs = target_LRR_share$se_logit
      )
  }
  
  llik
}

# Likelihood
calc_likelihood <- function(v_params) {
  exp(calc_log_lik(v_params))
}

# Record start time of calibration
t_init <- Sys.time()

# Bayesian calibration using IMIS
prior        <- calc_prior
likelihood   <- calc_likelihood
sample.prior <- sample_prior

# Run IMIS
set.seed(456)
n_resamp <- 1000
fit_imis <- IMIS::IMIS(
  B = 1000,
  B.re = n_resamp,
  number_k = 10,
  D = 0
)

# Posterior summaries
m_calib_res <- as.data.frame(fit_imis$resample)
colnames(m_calib_res) <- v_param_names
post_draws <- m_calib_res[, v_param_names, drop = FALSE]
m_calib_res_95cr <- matrixStats::colQuantiles(
  as.matrix(post_draws),
  probs = c(0.025, 0.5, 0.975)
)

# Posterior median and mean
v_calib_post_median <- m_calib_res_95cr[, "50%"]
v_calib_post_mean <- colMeans(post_draws)

# Save results
saveRDS(m_calib_res, "m_calib_res.rds")
saveRDS(post_draws, "post_draws.rds")
saveRDS(m_calib_res_95cr, "m_calib_res_95cr.rds")
saveRDS(v_calib_post_median, "v_calib_post_median.rds")
saveRDS(v_calib_post_mean, "v_calib_post_mean.rds")

# Model output at posterior median
pred_median_TDM1 <- predict_all_endpoints_calib(
  par_vec = v_calib_post_median,
  base_params = l_params_calib,
  arm = "TDM1"
)
pred_median_TDXd <- predict_all_endpoints_calib(
  par_vec = v_calib_post_median,
  base_params = l_params_calib,
  arm = "TDXd"
)
pred_median_TDM1
pred_median_TDXd

# Model output at posterior mean
pred_mean_TDM1 <- predict_all_endpoints_calib(
  par_vec = v_calib_post_mean,
  base_params = l_params_calib,
  arm = "TDM1"
)

pred_mean_TDXd <- predict_all_endpoints_calib(
  par_vec = v_calib_post_mean,
  base_params = l_params_calib,
  arm = "TDXd"
)
pred_mean_TDM1
pred_mean_TDXd

## Posterior diagnostics
# Calibrated parameters
param_order <- c(
  "haz_RF_DR_base",
  "LRR_share",
  "LRR_DR_multiplier",
  "haz_DR_BCd_base",
  "hr_RF_DR_TDXd_vs_TDM1"
)

param_labels <- c(
  haz_RF_DR_base = "RF to DR hazard\n(T-DM1)",
  LRR_share = "LRR share\nof first recurrences",
  LRR_DR_multiplier = "LRR to DR multiplier",
  haz_DR_BCd_base = "DR to BC death\nhazard",
  hr_RF_DR_TDXd_vs_TDM1 = "HR for RF to DR\n(T-DXd vs T-DM1)"
)

# Marginal posterior distributions and pairwise correlations plot
gg_post_pairs_corr <- GGally::ggpairs(
  data.frame(m_calib_res[, v_param_names]),
  upper = list(
    continuous = wrap("cor", color = "black", size = 5)
  ),
  diag = list(
    continuous = wrap("barDiag", alpha = 0.8)
  ),
  lower = list(
    continuous = wrap("points", alpha = 0.3, size = 0.7)
  ),
  columnLabels = unname(param_labels[param_order])
) +
  theme_bw(base_size = 14.5) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_text(size = 9.5),
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text = element_text(hjust = 0, lineheight = 0.9)
  )

print(gg_post_pairs_corr)

# Prior vs posterior densities
m_samp_prior <- sample.prior(n_resamp)

df_samp_prior <- reshape2::melt(
  cbind(PDF = "Prior", as.data.frame(m_samp_prior)),
  variable.name = "Parameter"
)

df_samp_post_imis <- reshape2::melt(
  cbind(PDF = "Posterior IMIS", as.data.frame(m_calib_res[, v_param_names])),
  variable.name = "Parameter"
)

df_samp_prior_post <- dplyr::bind_rows(df_samp_prior, df_samp_post_imis)
df_samp_prior_post <-
  df_samp_prior_post %>%
  mutate(
    Parameter = factor(
      Parameter,
      levels = param_order,
      labels = unname(param_labels[param_order])
    )
  )

gg_prior_post_imis <- ggplot(
  df_samp_prior_post,
  aes(x = value, y = after_stat(density), fill = PDF)
) +
  facet_wrap(
    ~ Parameter,
    scales = "free",
    ncol = 2
  ) +
  scale_x_continuous(n.breaks = 6) +
  geom_density(alpha = 0.5) +
  theme_bw(base_size = 16) +
  theme(
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

print(gg_prior_post_imis)

## Internal validation: calibration plot
df_targets <- targets_endpoints %>%
  mutate(
    Type    = "Target",
    Arm     = as.character(arm),
    Outcome = outcome,
    time    = time,
    value   = value,
    lb      = lower,
    ub      = upper
  ) %>%
  select(Type, Arm, Outcome, time, value, lb, ub)
df_targets$Outcome <- factor(df_targets$Outcome, levels = c("iDFS", "DRFI", "OS"))

# Generate endpoint predictions across all posterior draws for each arm
pred_rows <- list()
for (arm in strategy_names) {
  pred_mat <- matrix(0, nrow = nrow(post_draws), ncol = 3)
  colnames(pred_mat) <- c("iDFS", "DRFI", "OS")
  
  for (i in 1:nrow(post_draws)) {
    par_row <- as.numeric(post_draws[i, ])
    names(par_row) <- v_param_names
    out <- predict_all_endpoints_calib(
      par_vec     = par_row,
      base_params = l_params_calib,
      arm         = arm
    )
    
    pred_mat[i, ] <- c(as.numeric(out$iDFS), 
                       as.numeric(out$DRFI),
                       as.numeric(out$OS))
  }
  
  pred_rows[[arm]] <- data.frame(Arm = arm, pred_mat, check.names = FALSE)
}

# Bind rows and reshape to long format
df_pred_wide <- bind_rows(pred_rows)
df_pred_long <- df_pred_wide %>%
  pivot_longer(
    cols = c("iDFS", "DRFI", "OS"),
    names_to = "Outcome",
    values_to = "value"
  ) %>%
  mutate(
    Type = "Model",
    time = 3,
    Outcome = factor(Outcome, levels = c("iDFS", "DRFI", "OS"))
  )

# Summarize predictions to get mean and 95% credible intervals
df_model_sum <- df_pred_long %>%
  group_by(Type, Arm, Outcome, time) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    lb = as.numeric(
      quantile(value, probs = 0.025, na.rm = TRUE)
    ),
    ub = as.numeric(
      quantile(value, probs = 0.975, na.rm = TRUE)
    ),
    n_draws = n(),
    .groups = "drop"
  ) %>%
  rename(value = mean_value)

# Bind target and model summaries for plotting
df_plot <- bind_rows(df_targets, df_model_sum) %>%
  mutate(
    Outcome = factor(Outcome, levels = c("iDFS", "DRFI", "OS")),
    Type    = factor(Type, levels = c("Model", "Target"))
  ) %>%
  select(Type, Arm, Outcome, time, value, lb, ub)

# Plot internal validation calibration plot
ggplot(df_plot, aes(x = Outcome, y = value, ymin = lb, ymax = ub,
                        color = Type, group = Type)) +
  geom_errorbar(position = position_dodge(width = 0.55), width = 0.18,
                linewidth = 0.9) +
  geom_point(position = position_dodge(width = 0.55), size = 2.5) +
  facet_wrap(~ Arm) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(y = "3-year probability") +
  theme_bw(base_size = 18) +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c(Model = "steelblue3", Target = "black"))

# Posterior-mean calibrated or derived parameters
haz_RF_DR_TDM1_post <- as.numeric(v_calib_post_mean["haz_RF_DR_base"])
LRR_share_post <- as.numeric(v_calib_post_mean["LRR_share"])
haz_RF_LRR_post <- (LRR_share_post / (1 - LRR_share_post)) * haz_RF_DR_TDM1_post
LRR_DR_multiplier_post <- v_calib_post_mean["LRR_DR_multiplier"]
haz_LRR_DR_post <-  LRR_DR_multiplier_post * haz_RF_DR_TDM1_post
hr_RF_DR_TDXd_vs_TDM1_post <- as.numeric(v_calib_post_mean["hr_RF_DR_TDXd_vs_TDM1"])
haz_DR_BCd_post <- as.numeric(v_calib_post_mean["haz_DR_BCd_base"])

# Full base-case parameter list
l_params_all <- list(
  # General model settings
  delta            = delta,
  age_init         = age_init,
  tx_window_months = tx_window_months,
  horizon_years    = 10,
  
  # Treatment-effect persistence assumptions
  hr_full_effect_years   = 3,
  hr_wane_end_year       = 10,
  
  # Calibrated parameters
  haz_RF_LRR            = haz_RF_LRR_post,
  haz_RF_DR_TDM1        = haz_RF_DR_TDM1_post,
  LRR_share             = LRR_share_post,
  haz_LRR_DR            = haz_LRR_DR_post,
  LRR_DR_multiplier = LRR_DR_multiplier_post,
  hr_RF_DR_TDXd_vs_TDM1 = hr_RF_DR_TDXd_vs_TDM1_post,
  haz_DR_BCd            = haz_DR_BCd_post,
  
  # Other-cause mortality
  haz_OC = haz_OC,
  
  # ILD parameters
  p_any_ILD_TDM1 = p_any_ILD_TDM1_base,
  p_any_ILD_TDXd = p_any_ILD_TDXd_base,
  
  p_fatal_ILD_all_TDXd = p_fatal_ILD_all_TDXd_base,
  
  prop_fatal_ILD_TDM1 = prop_fatal_ILD_TDM1,
  prop_fatal_ILD_TDXd = prop_fatal_ILD_TDXd,
  
  #prop_dis_ILD_TDM1 = prop_dis_ILD_TDM1,
  #prop_dis_ILD_TDXd = prop_dis_ILD_TDXd,
  
  prop_dis_given_nonfatal_ILD_TDM1 = prop_dis_given_nonfatal_ILD_TDM1,
  prop_dis_given_nonfatal_ILD_TDXd = prop_dis_given_nonfatal_ILD_TDXd,
  
  # A-priori-specified fractional benefit following ILD discontinuation
  benefit_frac_postILD = c(
    0.28, 0.49, 0.64, 0.75, 0.83,
    0.89, 0.94, 0.97, 0.99, 1.00
  )
)
saveRDS(l_params_all, "l_params_all_basecase.rds")

## Map DR hazard to RCB categories and HR subgroups
# 5-year scaled T-DM1 (baseline) EFS from Yau et al
RCB_EFS_scaled <- data.frame(
  RCB = c("RCB-I", "RCB-II", "RCB-III", "RCB-I", "RCB-II", "RCB-III"),
  HR_status = c("HR-/HER2+", "HR-/HER2+", "HR-/HER2+", "HR+/HER2+", "HR+/HER2+", "HR+/HER2+"),
  EFS_5y_TDM1 = c(1 - 0.084, 1 - 0.221, 1 - 0.241, 1 - 0.050, 1 - 0.138, 1 - 0.283)
) %>%
  mutate(
    RCB = factor(RCB, levels = c("RCB-I", "RCB-II", "RCB-III")),
    HR_status = factor(HR_status, levels = c("HR-/HER2+", "HR+/HER2+"))
  )

# Grid search to find multiplier applied to DR hazard that most closely reproduce
# scaled EFS estimates
# recurrence_multiplier_grid <- seq(0.01, 10, length.out = 1000)
# Given that the DR hazard is likely to be small, we can use a log-spaced grid 
# to better capture the lower end of the multiplier range
recurrence_multiplier_grid <- exp(seq(log(1e-3), log(10), length.out = 1000)) # log-spaced

# Store base parameters
l_params_rcb_base <- l_params_all
l_params_rcb_base$horizon_years <- 5
haz_RF_DR_TDM1_base <- l_params_all$haz_RF_DR_TDM1
haz_RF_LRR_base     <- l_params_all$haz_RF_LRR

# Get modeled 5-year T-DM1 EFS (corresponding to iDFS in our model) for one DR multiplier
get_EFS_5y_TDM1 <- function(recurrence_multiplier, params_base) {
  
  params <- params_base
  
  # Apply multiplier to RF to DR hazard
  params$haz_RF_DR_TDM1 <- haz_RF_DR_TDM1_base * recurrence_multiplier
  
  # Apply multiplier to RF to LRR hazard
  params$haz_RF_LRR <- haz_RF_LRR_base * recurrence_multiplier
  
  out <- run_markov_model(params)
  
  as.numeric(out$iDFS["TDM1"])
}

# Root-finding approach
fit_multiplier <- function(target_EFS, params_base) {
  f <- function(m) get_EFS_5y_TDM1(m, params_base) - target_EFS
  uniroot(f, interval = c(1e-4, 20), tol = 1e-8)$root
}

# Do it for all scenarios
RCB_EFS_scaled_multiplier <- RCB_EFS_scaled %>%
  mutate(
    recurrence_multiplier = map_dbl(
      EFS_5y_TDM1,
      ~ fit_multiplier(
        target_EFS = .x,
        params_base = l_params_rcb_base
      )
    ),
    
    haz_RF_DR_TDM1 =
      haz_RF_DR_TDM1_base * recurrence_multiplier,
    
    haz_RF_LRR =
      haz_RF_LRR_base * recurrence_multiplier
  )
RCB_EFS_scaled_multiplier

# Save results
saveRDS(RCB_EFS_scaled_multiplier, "RCB_EFS_scaled_multiplier.rds")

# Highest-risk T-DM1 RCB/HR subgroup
haz_RF_DR_max <- max(RCB_EFS_scaled_multiplier$haz_RF_DR_TDM1)

# Apply calibrated LRR to DR hazard ratio in the highest-risk subgroup
haz_LRR_DR_common <- as.numeric(l_params_all$LRR_DR_multiplier * haz_RF_DR_max)

# Treatment persistence scenarios
treatment_effect_scenarios <- data.frame(persistence_scenario = c("Effect ends at year 3",
                                                                  "Waning completed by year 5",
                                                                  "Waning completed by year 7",
                                                                  "Waning completed by year 10",
                                                                  "Persistent through year 10"),
                                         hr_full_effect_years = c(3,3,3,3,10),
                                         hr_wane_end_year = c(3,5,7,10,10))

# Model input scenarios (base case + 3 clinically relevant edge cases) 
base_case <- function(params) {
  params
}

weaker_hr <- function(params) {
  params$hr_RF_DR_TDXd_vs_TDM1 <- 0.80
  params
}

higher_ild_incidence <- function(params) {
  params$p_any_ILD_TDXd <- min(1 - 1e-12, params$p_any_ILD_TDXd * 2)
  params$p_fatal_ILD_all_TDXd <-
    params$p_any_ILD_TDXd * params$prop_fatal_ILD_TDXd
  
  params
}

combined <- function(params) {
  params <- weaker_hr(params)
  params <- higher_ild_incidence(params)
  params
}

# Named scenarios for optional deterministic analyses
det_scenarios <- list(
  "Base case" = base_case,
  "HR 0.80" = weaker_hr,
  "Any-grade ILD incidence x2" = higher_ild_incidence,
  "Combined" = combined
)
det_scenario_order <- names(det_scenarios)

## Optional deterministic analyses (not described in manuscript)
# Evaluate each scenario across RCB categories and HR subgroups
eval_os_RCB <- function(params,
                        det_scenario_name,
                        det_scenario_modify = function(params) params) {
  # Initialize results list
  results <- vector("list", nrow(treatment_effect_scenarios)*
                      nrow(RCB_EFS_scaled_multiplier))
  result_index <- 1
  
  # Loop over all treatment persistence scenarios
  for (j in seq_len(nrow(treatment_effect_scenarios))) {
    persistence_row <- treatment_effect_scenarios[j, ]
    
    # Loop over all DR multipliers (corresponding to RCB/HR subgroups)
    for (i in seq_len(nrow(RCB_EFS_scaled_multiplier))) {
      
      multiplier_row <- RCB_EFS_scaled_multiplier[i, ]
      
      current_params <- params
      
      # Set treatment persistence scenario
      current_params$hr_full_effect_years <- persistence_row$hr_full_effect_years
      current_params$hr_wane_end_year <- persistence_row$hr_wane_end_year
      
      # Set RCB/HR-specific baseline T-DM1 DR hazard
      current_params$haz_RF_DR_TDM1 <- multiplier_row$haz_RF_DR_TDM1
      
      # Set RCB/HR-specific LRR hazard (same for T-DM1 and T-DXd)
      current_params$haz_RF_LRR <- multiplier_row$haz_RF_LRR
      
      # Set common LRR to DR hazard (same across subgroups and treatment arms)
      current_params$haz_LRR_DR <- haz_LRR_DR_common
        
      # Apply scenario after setting subgroup-specific hazards
      current_params <- det_scenario_modify(current_params)
      
      # Run Markov model
      outcomes_10y <- run_markov_model(params = current_params)
      
      # Output results for this scenario
      results[[result_index]] <- data.frame(
        det_scenario = det_scenario_name,
        persistence_scenario = persistence_row$persistence_scenario,
        HR_status = multiplier_row$HR_status,
        RCB = multiplier_row$RCB,
        
        # Overall survival difference, percentage points
        OS_10y_TDM1 = as.numeric(outcomes_10y$OS["TDM1"]),
        OS_10y_TDXd = as.numeric(outcomes_10y$OS["TDXd"]),
        OS_diff_pp = 100 * (as.numeric(outcomes_10y$OS["TDXd"]) - 
                              as.numeric(outcomes_10y$OS["TDM1"])),
        
        # Breast cancer death difference, percentage points
        BC_death_10y_TDM1 = as.numeric(outcomes_10y$BC_death["TDM1"]),
        BC_death_10y_TDXd = as.numeric(outcomes_10y$BC_death["TDXd"]),
        BC_death_diff_pp = 100 * (as.numeric(outcomes_10y$BC_death["TDM1"]) -
                                    as.numeric(outcomes_10y$BC_death["TDXd"])),
        
        # ILD death difference, percentage points
        ILD_death_10y_TDM1 = as.numeric(outcomes_10y$ILD_death["TDM1"]),
        ILD_death_10y_TDXd = as.numeric(outcomes_10y$ILD_death["TDXd"]),
        ILD_death_diff_pp = 100 * (as.numeric(outcomes_10y$ILD_death["TDM1"]) -
                                     as.numeric(outcomes_10y$ILD_death["TDXd"])),
        
        # Other-cause death difference, percentage points
        OC_death_10y_TDM1 = as.numeric(outcomes_10y$OC_death["TDM1"]),
        OC_death_10y_TDXd = as.numeric(outcomes_10y$OC_death["TDXd"]),
        OC_death_diff_pp = 100 * (as.numeric(outcomes_10y$OC_death["TDM1"]) -
                                    as.numeric(outcomes_10y$OC_death["TDXd"]))
      )
      result_index <- result_index + 1
    }
  }
  dplyr::bind_rows(results)
}

# Run deterministic analysis for each scenario
df_os_RCB_list <- list()

for (s in seq_along(det_scenarios)) {
  det_scenario_name <- names(det_scenarios)[s]
  det_scenario_modify <- det_scenarios[[s]]
  
  cat("Running deterministic analysis:", det_scenario_name, "\n")
  
  df_os_RCB_list[[s]] <- eval_os_RCB(
    params = l_params_all,
    det_scenario_name = det_scenario_name,
    det_scenario_modify = det_scenario_modify
  )
}

# Summarize deterministic analysis
df_os_RCB <- dplyr::bind_rows(df_os_RCB_list) %>%
  mutate(det_scenario = factor(det_scenario, levels = det_scenario_order),
         persistence_scenario = factor(persistence_scenario, levels = 
                                         treatment_effect_scenarios$persistence_scenario)) %>%
  arrange(persistence_scenario, HR_status, RCB, det_scenario)

table_os_RCB <- df_os_RCB %>%
  mutate(OS_diff_pp = round(OS_diff_pp, 2)) %>%
  select(persistence_scenario, det_scenario, HR_status, RCB, 
         OS_diff_pp, BC_death_diff_pp, ILD_death_diff_pp, OC_death_diff_pp) %>%
  arrange(persistence_scenario, det_scenario, HR_status, RCB)
table_os_RCB

# Save results
saveRDS(df_os_RCB, "df_os_RCB.rds")
saveRDS(table_os_RCB, "table_os_RCB.rds")

# Named scenarios for probabilistic sensitivity analysis (PSA)
psa_scenarios <- list(
  "Base case PSA" = base_case,
  "HR 0.80 PSA" = weaker_hr,
  "Any-grade ILD incidence x2 PSA" = higher_ild_incidence,
  "Combined PSA" = combined
)
psa_scenario_order <- names(psa_scenarios)

## PSA
n_sim <- 1000

# Function to sample probabilities from a logit-normal distribution
draw_prob_logitnormal <- function(mu, lower_95, upper_95, n_sim) {
  mean_logit <- qlogis(mu)
  se_logit <- (qlogis(upper_95) - qlogis(lower_95)) / (2 * 1.96)
  plogis(rnorm(n = n_sim, mean = mean_logit, sd = se_logit))
}

# Draw calibrated model parameters from IMIS posterior
set.seed(1)
calib_param_draws <- as.data.frame(post_draws)

# Any-grade ILD probability over treatment window: T-DXd
p_any_ILD_TDXd_draw <- draw_prob_logitnormal(
  mu = 0.096,
  lower_95 = 0.05,
  upper_95 = 0.21,
  n_sim = n_sim
)

# Any-grade ILD probability over treatment window: T-DM1; uniform Beta(1,1) prior
p_any_ILD_TDM1_draw <- rbeta(
  n = n_sim,
  shape1 = any_ild_TDM1 + 1,
  shape2 = n_TDM1 - any_ild_TDM1 + 1
)

# Overall fatal ILD probability among all T-DXd-assigned patients
p_fatal_ILD_all_TDXd_draw <- draw_prob_logitnormal(
  mu = 0.01,
  lower_95 = 0.002,
  upper_95 = 0.030,
  n_sim = n_sim
)

# Histograms of PSA draws
hist(calib_param_draws$haz_RF_DR_base)
hist(calib_param_draws$LRR_share)
hist(calib_param_draws$LRR_DR_multiplier)
hist(calib_param_draws$haz_DR_BCd_base)
hist(calib_param_draws$hr_RF_DR_TDXd_vs_TDM1)
hist(p_any_ILD_TDXd_draw)
hist(p_fatal_ILD_all_TDXd_draw)

# Highest-risk T-DM1 RCB/HR recurrence multiplier
max_recurrence_multiplier <- max(RCB_EFS_scaled_multiplier$recurrence_multiplier)

# Add derived quantities to posterior draws
calib_and_derived_draws <- calib_param_draws %>%
    mutate(
        # RF to LRR hazard
        haz_RF_LRR = (LRR_share / (1 - LRR_share)) * haz_RF_DR_base,
        
        # Highest-risk T-DM1 RF to DR hazard
        haz_RF_DR_max = max_recurrence_multiplier * haz_RF_DR_base,
        
        # Common LRR -> DR hazard used in final RCB/HR model
        haz_LRR_DR = LRR_DR_multiplier * haz_RF_DR_max
    )

# Summarize distributions used in PSA for calibrated and derived parameters
calib_and_derived_summary <- calib_and_derived_draws %>%
    select(haz_RF_DR_base, LRR_share, LRR_DR_multiplier, haz_DR_BCd_base, 
           hr_RF_DR_TDXd_vs_TDM1, haz_RF_LRR, haz_LRR_DR) %>%
    pivot_longer(cols = everything(),
                 names_to = "parameter",
                 values_to = "value") %>%
    group_by(parameter) %>%
    summarise(mean = mean(value, na.rm = TRUE),
              lower_95CrI = quantile(value, 0.025, na.rm = TRUE),
              upper_95CrI = quantile(value, 0.975, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(type = ifelse(parameter %in% c("haz_RF_LRR", "haz_LRR_DR"), "Derived", 
                         "Calibrated"))
  
# Run PSA
eval_os_RCB_PSA <- function(params, psa_scenario_name, psa_scenario_modify,
                            n_sim = 1000, calib_and_derived_draws,
                            p_any_ILD_TDM1_draw, p_any_ILD_TDXd_draw,
                            p_fatal_ILD_all_TDXd_draw) {
    results <- vector("list", nrow(RCB_EFS_scaled_multiplier))
    
    # Loop over RCB/HR subgroups
    for (i in seq_len(nrow(RCB_EFS_scaled_multiplier))) {
        multiplier_row <- RCB_EFS_scaled_multiplier[i, ]
        
        sim_results <- vector("list", n_sim * nrow(treatment_effect_scenarios))
        
        result_index <- 1
        
        for (s in seq_len(n_sim)) {
            current_params <- params
            
            # PSA draw: calibrated parameters
            haz_RF_DR_TDM1_base_s <-
                as.numeric(calib_and_derived_draws$haz_RF_DR_base[s])
            LRR_share_s <-
                as.numeric(calib_and_derived_draws$LRR_share[s])
            haz_DR_BCd_s <-
                as.numeric(calib_and_derived_draws$haz_DR_BCd_base[s])
            hr_RF_DR_s <-
                as.numeric(calib_and_derived_draws$hr_RF_DR_TDXd_vs_TDM1[s])
            
            # RCB/HR-specific RF to DR hazard
            current_params$haz_RF_DR_TDM1 <- haz_RF_DR_TDM1_base_s *
                as.numeric(multiplier_row$recurrence_multiplier)
            
            # RCB/HR-specific RF to LRR hazard
            current_params$haz_RF_LRR <- (LRR_share_s / (1 - LRR_share_s)) *
                (haz_RF_DR_TDM1_base_s *
                     as.numeric(multiplier_row$recurrence_multiplier))
            
            # Common LRR to DR hazard
            current_params$haz_LRR_DR <-
                as.numeric(calib_and_derived_draws$haz_LRR_DR[s])
      
            # PSA draw: DR to BC death hazard
            current_params$haz_DR_BCd <- haz_DR_BCd_s
            
            # PSA draw: HR for RF to DR, T-DXd vs T-DM1
            current_params$hr_RF_DR_TDXd_vs_TDM1 <- hr_RF_DR_s
            
            # PSA draw: any-grade ILD probability for T-DXd
            current_params$p_any_ILD_TDXd <- p_any_ILD_TDXd_draw[s]
            
            # PSA draw: any-grade ILD probability for T-DM1
            current_params$p_any_ILD_TDM1 <- p_any_ILD_TDM1_draw[s]
            
            # PSA draw: overall fatal ILD probability among all T-DXd-assigned patients
            current_params$p_fatal_ILD_all_TDXd <- p_fatal_ILD_all_TDXd_draw[s]
            
            # Ensure fatal ILD remains a subset of any-grade ILD before applying scenarios
            current_params$p_fatal_ILD_all_TDXd <- min(current_params$p_fatal_ILD_all_TDXd,
                                                       0.99 * current_params$p_any_ILD_TDXd)
            
            # Convert sampled overall fatal ILD probability to fatal proportion among ILD
            current_params$prop_fatal_ILD_TDXd <- current_params$p_fatal_ILD_all_TDXd /
              current_params$p_any_ILD_TDXd
            
            # Apply PSA scenarios
            current_params <- psa_scenario_modify(current_params)
            
            # Recalculate ILD quantities after scenario modification
            p_any_ILD_TDXd <- current_params$p_any_ILD_TDXd
            p_fatal_ILD_all_TDXd <- min(current_params$p_fatal_ILD_all_TDXd,
                                        0.99 * p_any_ILD_TDXd)
            p_nonfatal_ILD_all_TDXd <- p_any_ILD_TDXd - p_fatal_ILD_all_TDXd
            
            # Fixed probability of discontinuation among nonfatal ILD cases
            prop_dis_given_nonfatal_ILD_TDXd <-
              current_params$prop_dis_given_nonfatal_ILD_TDXd
            
            # Overall probability of ILD discontinuation among T-DXd-assigned patients
            p_dis_ILD_all_TDXd <- p_nonfatal_ILD_all_TDXd * 
                prop_dis_given_nonfatal_ILD_TDXd
            
            # Convert to proportions among any-grade ILD cases
            prop_fatal_ILD_TDXd <- p_fatal_ILD_all_TDXd / p_any_ILD_TDXd
            prop_dis_ILD_TDXd <- p_dis_ILD_all_TDXd / p_any_ILD_TDXd
            
            # Assign ILD parameters used by the model
            current_params$p_any_ILD_TDXd <- p_any_ILD_TDXd
            current_params$p_fatal_ILD_all_TDXd <- p_fatal_ILD_all_TDXd
            current_params$prop_fatal_ILD_TDXd <- prop_fatal_ILD_TDXd
            current_params$prop_dis_ILD_TDXd <- prop_dis_ILD_TDXd
            current_params$prop_dis_given_nonfatal_ILD_TDXd <-
              prop_dis_given_nonfatal_ILD_TDXd
            
            # Loop over treatment persistence assumptions
            for (j in seq_len(nrow(treatment_effect_scenarios))){
              persistence_row <- treatment_effect_scenarios[j,]
              persistence_params <- current_params
              persistence_params$hr_full_effect_years <- persistence_row$hr_full_effect_years
              persistence_params$hr_wane_end_year <- persistence_row$hr_wane_end_year
              
              # Run Markov model
              outcomes_10y <- run_markov_model(params = persistence_params)
              
              # Output results
              sim_results[[result_index]] <- data.frame(
                sim = s,
                psa_scenario = psa_scenario_name,
                persistence_scenario = persistence_row$persistence_scenario,
                HR_status = multiplier_row$HR_status,
                RCB = multiplier_row$RCB,
                # Overall survival difference, percentage points
                OS_10y_TDM1 = as.numeric(outcomes_10y$OS["TDM1"]),
                OS_10y_TDXd = as.numeric(outcomes_10y$OS["TDXd"]),
                OS_diff_pp = 100 * (as.numeric(outcomes_10y$OS["TDXd"]) - 
                                      as.numeric(outcomes_10y$OS["TDM1"])),
              
                # Breast cancer death difference, percentage points
                BC_death_10y_TDM1 = as.numeric(outcomes_10y$BC_death["TDM1"]),
                BC_death_10y_TDXd = as.numeric(outcomes_10y$BC_death["TDXd"]),
                BC_death_diff_pp = 100 * (as.numeric(outcomes_10y$BC_death["TDM1"]) -
                                            as.numeric(outcomes_10y$BC_death["TDXd"])),
                
                # ILD death difference, percentage points
                ILD_death_10y_TDM1 = as.numeric(outcomes_10y$ILD_death["TDM1"]),
                ILD_death_10y_TDXd = as.numeric(outcomes_10y$ILD_death["TDXd"]),
                ILD_death_diff_pp = 100 * (as.numeric(outcomes_10y$ILD_death["TDM1"]) -
                                             as.numeric(outcomes_10y$ILD_death["TDXd"])),
                
                # Other-cause death difference, percentage points
                OC_death_10y_TDM1 = as.numeric(outcomes_10y$OC_death["TDM1"]),
                OC_death_10y_TDXd = as.numeric(outcomes_10y$OC_death["TDXd"]),
                OC_death_diff_pp = 100 * (as.numeric(outcomes_10y$OC_death["TDM1"]) -
                                            as.numeric(outcomes_10y$OC_death["TDXd"]))
                )
              
              result_index <- result_index + 1
      }
    }
    results[[i]] <- dplyr::bind_rows(sim_results)
  }
  
  dplyr::bind_rows(results)
}

# Conduct probabilistic analyses across scenarios
df_os_RCB_psa_list <- list()

for (s in seq_along(psa_scenarios)) {
  psa_scenario_name <- names(psa_scenarios)[s]
  psa_scenario_modify <- psa_scenarios[[s]]
  
  cat("Running PSA:", psa_scenario_name, "\n")
  
  df_os_RCB_psa_list[[s]] <- eval_os_RCB_PSA(
    params = l_params_all,
    psa_scenario_name = psa_scenario_name,
    psa_scenario_modify = psa_scenario_modify,
    n_sim = n_sim,
    calib_and_derived_draws = calib_and_derived_draws,
    p_any_ILD_TDM1_draw = p_any_ILD_TDM1_draw,
    p_any_ILD_TDXd_draw = p_any_ILD_TDXd_draw,
    p_fatal_ILD_all_TDXd_draw = p_fatal_ILD_all_TDXd_draw
  )
}

# Organize probabilistic analysis results into dataframe
df_os_RCB_psa <- dplyr::bind_rows(df_os_RCB_psa_list) %>%
  mutate(psa_scenario = factor(psa_scenario, levels = psa_scenario_order),
         persistence_scenario = factor(persistence_scenario,
                                       levels = treatment_effect_scenarios$persistence_scenario)) %>%
  arrange(psa_scenario, persistence_scenario, HR_status, RCB, sim)

# Summarize results of probabilistic analysis
table_os_RCB_psa <- df_os_RCB_psa %>%
  group_by(psa_scenario, persistence_scenario, HR_status, RCB) %>%
  summarise(
    # Overall survival difference, percentage points
    mean_OS_diff_pp = mean(OS_diff_pp, na.rm = TRUE),
    lower_OS_diff_pp = quantile(OS_diff_pp, 0.025, na.rm = TRUE),
    upper_OS_diff_pp = quantile(OS_diff_pp, 0.975, na.rm = TRUE),
    
    # Breast cancer death difference, percentage points
    mean_BC_death_diff_pp = mean(BC_death_diff_pp, na.rm = TRUE),
    lower_BC_death_diff_pp = quantile(BC_death_diff_pp, 0.025, na.rm = TRUE),
    upper_BC_death_diff_pp = quantile(BC_death_diff_pp, 0.975, na.rm = TRUE),
    
    # ILD death difference, percentage points
    mean_ILD_death_diff_pp = mean(ILD_death_diff_pp, na.rm = TRUE),
    lower_ILD_death_diff_pp = quantile(ILD_death_diff_pp, 0.025, na.rm = TRUE),
    upper_ILD_death_diff_pp = quantile(ILD_death_diff_pp, 0.975, na.rm = TRUE),
    
    # Other-cause death difference, percentage points
    mean_OC_death_diff_pp = mean(OC_death_diff_pp, na.rm = TRUE),
    lower_OC_death_diff_pp = quantile(OC_death_diff_pp, 0.025, na.rm = TRUE),
    upper_OC_death_diff_pp = quantile(OC_death_diff_pp, 0.975, na.rm = TRUE),
    
    # Percent favoring T-DXd
    pct_favors_TDXd_OS = 100 * mean(OS_diff_pp > 0, na.rm = TRUE),
    pct_BC_death_favors_TDXd = 100 * mean(BC_death_diff_pp > 0, na.rm = TRUE),
    pct_ILD_death_favors_TDXd = 100 * mean(ILD_death_diff_pp > 0, na.rm = TRUE),
    pct_OC_death_favors_TDXd = 100 * mean(OC_death_diff_pp > 0, na.rm = TRUE),

    .groups = "drop"
  ) %>%
  arrange(psa_scenario, persistence_scenario, HR_status, RCB)

table_os_RCB_psa

# Save results
saveRDS(df_os_RCB_psa, "df_os_RCB_PSA.rds")
saveRDS(table_os_RCB_psa, "table_os_RCB_PSA.rds")

## One-way ILD threshold analysis
threshold_RCB <- function(params,
                             scenario_modify = function(params) params,
                             multiplier_grid = RCB_EFS_scaled_multiplier) {
  # Initialize results list
  results <- vector("list", nrow(treatment_effect_scenarios)*
                      nrow(multiplier_grid))
  result_index <- 1
  
  # Loop over treatment persistence assumptions
  for (i in seq_len(nrow(treatment_effect_scenarios))){
    persistence_row <- treatment_effect_scenarios[i, ]
    
    # Loop over RCB/HR subgroups
    for (j in seq_len(nrow(multiplier_grid))){
      multiplier_row <- multiplier_grid[j, ]
      current_params <- params
      
      # Set treatment persistence scenarios
      current_params$hr_full_effect_years <- persistence_row$hr_full_effect_years
      current_params$hr_wane_end_year <- persistence_row$hr_wane_end_year
      
      # Set subgroup-specific hazards
      current_params$haz_RF_DR_TDM1 <- multiplier_row$haz_RF_DR_TDM1
      current_params$haz_RF_LRR <- multiplier_row$haz_RF_LRR
      
      # Set common LRR to DR hazard (same across subgroups and treatment arms)
      current_params$haz_LRR_DR <- haz_LRR_DR_common
      
      # Apply threshold after subgroup-specific hazards are set
      current_params <- scenario_modify(current_params)
      
      # Run Markov model
      outcomes_10y <- run_markov_model(params = current_params)
      
      # Output results
      results[[result_index]] <- data.frame(
        persistence_scenario = persistence_row$persistence_scenario,
        HR_status = multiplier_row$HR_status,
        RCB = multiplier_row$RCB,
        recurrence_multiplier = multiplier_row$recurrence_multiplier,
        
        # Overall survival difference, percentage points
        OS_diff_pp = 100 * (as.numeric(outcomes_10y$OS["TDXd"]) -
                              as.numeric(outcomes_10y$OS["TDM1"]))
      )
      result_index <- result_index + 1
    }
  }
  dplyr::bind_rows(results)
  
}

# Apply multiplier to T-DXd any-grade ILD incidence
modify_param_ild <- function(params, x) {
  params$p_any_ILD_TDXd <- min(1 - 1e-12, params$p_any_ILD_TDXd * x)
  params$p_fatal_ILD_all_TDXd <- params$p_any_ILD_TDXd * params$prop_fatal_ILD_TDXd
  params
}

# Maximum possible multiplier before any-grade ILD reaches 1
max_ild_multiplier <- (1 - 1e-12) / l_params_all$p_any_ILD_TDXd
ild_grid <- seq(from = 1, to = max_ild_multiplier, length.out = 500)

# Evaluate the ILD multiplier grid 
diff_ild_list <- vector("list", length(ild_grid))
for (i in seq_along(ild_grid)){
  diff_ild_list[[i]] <- threshold_RCB(params = l_params_all, 
                                            scenario_modify = function(params) {
                                              modify_param_ild(params = params,
                                                               x = ild_grid[i])}) %>%
    mutate(threshold = ild_grid[i])
}
diff_ild <- dplyr::bind_rows(diff_ild_list)

# Find first ILD multiplier where 10-year OS difference is negative (T-DM1 favored)
first_neg_diff_ild <- diff_ild %>%
  arrange(persistence_scenario, HR_status, RCB, threshold) %>%
  group_by(persistence_scenario, HR_status, RCB) %>%
  summarise(
    # OS difference < 0
    OS_threshold = if (any(OS_diff_pp < 0)) {
        threshold[which(OS_diff_pp < 0)[1]]
      } else {
        NA_real_
      },
    .groups = "drop"
  )

first_neg_diff_ild <- first_neg_diff_ild %>%
  # Threshold near 1 indicates T-DM1 favored in base case
  mutate(
    OS_threshold = case_when(
      is.na(OS_threshold) ~ "No threshold found",
      TRUE ~ paste0(round(OS_threshold, 2), "x")
    )
  )

# Save results
saveRDS(diff_ild, "diff_ild.rds")
saveRDS(first_neg_diff_ild, "first_neg_diff_ild.rds")

## Plots
# Desired top-to-bottom subgroup ordering
subgroup_order_top_to_bottom <- c(
  "RCB-I HR+",
  "RCB-I HR-",
  "RCB-II HR+",
  "RCB-II HR-",
  "RCB-III HR+",
  "RCB-III HR-"
)
subgroup_levels_for_plot <- rev(subgroup_order_top_to_bottom)

# Colors used for the deterministic plot
subgroup_colors <- c(
  `FALSE` = "#1D9E75",  # Other subgroups
  `TRUE`  = "#EF9F27"   # RCB-I HR+
)

# Colors used for the PSA distributions
distribution_side_colors <- c(
  `TRUE`  = "#F2C94C",  # OS difference < 0: favors T-DM1
  `FALSE` = "#1D9E75"   # OS difference >= 0: favors T-DXd
)

# Shared axis labels
deterministic_x_label <- paste0(
  "Difference in 10-year OS, percentage points\n",
  "Positive values favor T-DXd"
)
psa_x_label <- paste0(
  "Difference in 10-year OS, percentage points ",
  "(T-DXd minus T-DM1)"
)

# Shared formatting
format_subgroups <- function(data) {
  data %>%
    mutate(
      HR_short = recode(
        as.character(HR_status),
        "HR-/HER2+" = "HR-",
        "HR+/HER2+" = "HR+"
      ),
      subgroup = factor(
        paste(RCB, HR_short),
        levels = subgroup_levels_for_plot
      )
    )
}

# Treatment-persistence assumptions
persistence_levels <- c(
    "Effect ends at year 3",
    "Waning completed by year 5",
    "Waning completed by year 7",
    "Waning completed by year 10",
    "Persistent through year 10"
)

det_plots <- lapply(
    persistence_levels,
    function(persistence_name) {
        # Deterministic dumbbell plot
        det_plot_df <- df_os_RCB %>%
            filter(
                persistence_scenario == persistence_name,
                det_scenario %in% c("Base case", "Combined")
            ) %>%
            format_subgroups() %>%
            mutate(
                equivocal = RCB == "RCB-I" & HR_status == "HR+/HER2+"
            ) %>%
            select(
                RCB, HR_status, HR_short, subgroup, equivocal, det_scenario,
                OS_diff_pp
            ) %>%
            pivot_wider(
                id_cols = c(RCB, HR_status, HR_short, subgroup, equivocal),
                names_from = det_scenario,
                values_from = OS_diff_pp
            )
        
        gg_det_dumbbell <- ggplot(det_plot_df) +
            geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.5) +
            
            geom_segment(
                aes(
                    x = Combined, xend = `Base case`,
                    y = subgroup, yend = subgroup,
                    colour = equivocal
                ),
                linewidth = 1.1
            ) +
            
            geom_point(
                aes(x = Combined, y = subgroup),
                shape = 21, size = 3, fill = "white",
                stroke = 1, colour = "grey45"
            ) +
            
            geom_point(
                aes(x = `Base case`, y = subgroup, colour = equivocal),
                size = 3.4
            ) +
            
            geom_text(
                aes(
                    x = `Base case`,
                    y = subgroup,
                    label = sprintf("%+.1f", `Base case`)
                ),
                vjust = -1.1,
                size = 3
            ) +
            
            geom_text(
                aes(
                    x = Combined,
                    y = subgroup,
                    label = sprintf("%+.1f", Combined)
                ),
                vjust = -1.1,
                size = 3,
                colour = "grey45"
            ) +
            
            scale_colour_manual(values = subgroup_colors, guide = "none") +
            
            labs(
                x = deterministic_x_label,
                y = NULL,
                title = persistence_name
            ) +
            
            theme_minimal(base_size = 12) +
            
            theme(
                panel.grid.major.y = element_blank(),
                panel.grid.minor = element_blank()
            )
        
        return(gg_det_dumbbell)
    }
)
        
# Name each plot by its persistence assumption
names(det_plots) <- persistence_levels

# Print deterministic plots
# lapply(det_plots, print)
det_plots[["Effect ends at year 3"]]
det_plots[["Waning completed by year 5"]]
det_plots[["Waning completed by year 7"]]
det_plots[["Waning completed by year 10"]]
det_plots[["Persistent through year 10"]]

#ggsave(
#    filename = "det_persistent_year10.jpeg",
#    plot = det_plots[["Persistent through year 10"]],
#    width = 8,
#    height = 5.5,
#    units = "in",
#    dpi = 600,
#    quality = 100
#)

# Create a separate base-case PSA ridge plot
# for each treatment-persistence assumption
psa_plots <- lapply(
    persistence_levels,
    function(persistence_name) {
        
        # Base-case PSA for selected persistence assumption
        psa_plot_df <- df_os_RCB_psa %>%
            filter(
                psa_scenario == "Base case PSA",
                persistence_scenario == persistence_name
            ) %>%
            format_subgroups() %>%
            mutate(
                os_diff = OS_diff_pp
            )
        
        # Probability that OS difference favors T-DM1
        pct_tdm1 <- psa_plot_df %>%
            group_by(subgroup) %>%
            summarise(
                p_tdm1 = mean(os_diff < 0, na.rm = TRUE),
                .groups = "drop"
            )
        
        # X-axis limits specific to this persistence assumption
        x_limits <- range(
            psa_plot_df$os_diff,
            na.rm = TRUE
        )
        
        x_limits <- x_limits +
            c(-1, 1) * 0.08 * diff(x_limits)
        
        
        # Ridge plot
        gg_psa_ridges <- ggplot(
            psa_plot_df,
            aes(x = os_diff, y = subgroup)
        ) +
            ggridges::geom_density_ridges_gradient(
                aes(fill = after_stat(x < 0)),
                scale = 1.6,
                rel_min_height = 0,
                colour = "grey30",
                linewidth = 0.4
            ) +
            
            geom_vline(
                xintercept = 0,
                linetype = 2,
                colour = "grey35",
                linewidth = 0.6
            ) +
            
            geom_text(
                data = pct_tdm1,
                aes(
                    x = x_limits[2],
                    y = subgroup,
                    label = sprintf(
                        "%.0f%% <0",
                        100 * p_tdm1
                    )
                ),
                hjust = 1,
                vjust = -0.4,
                size = 3,
                colour = "grey20",
                inherit.aes = FALSE
            ) +
            
            scale_fill_manual(
                values = distribution_side_colors,
                breaks = c(TRUE, FALSE),
                labels = c(
                    "Favors T-DM1 (<0)",
                    "Favors T-DXd (\u22650)"
                ),
                name = NULL
            ) +
            
            coord_cartesian(
                xlim = x_limits
            ) +
            
            labs(
                title = persistence_name,
                x = psa_x_label,
                y = NULL
            ) +
            
            ggridges::theme_ridges(
                center_axis_labels = TRUE
            ) +
            
            theme(
                legend.position = "bottom",
                plot.title = element_text(
                    face = "bold"
                ),
                panel.grid.major.x = element_line(
                    colour = "grey90",
                    linewidth = 0.3
                ),
                plot.margin = margin(
                    5.5, 15, 5.5, 5.5
                )
            )
        
        return(gg_psa_ridges)
    }
)


# Name each plot based on treatment persistence assumption
names(psa_plots) <- persistence_levels

# Print PSA plots
# lapply(psa_plots, print)
psa_plots[["Effect ends at year 3"]]
psa_plots[["Waning completed by year 5"]]
psa_plots[["Waning completed by year 7"]]
psa_plots[["Waning completed by year 10"]]
psa_plots[["Persistent through year 10"]]

#ggsave(
#    filename = "PSA_persistent_year10.jpeg",
#    plot = psa_plots[["Persistent through year 10"]],
#    width = 8,
#    height = 5.5,
#    units = "in",
#    dpi = 600,
#    quality = 100
#)
                        
