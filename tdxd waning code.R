# Load packages
library(tidyverse)  
library(msm)
library(patchwork)
library(ggridges)
library(GGally)

# Global model parameters
tx_window_months <- 10
age_init <- 50
delta <- 1/12

# Strategy names
strategy_names <- c("TDM1", "TDXd")

# Death states
death_states <- c("OC_Death", "BC_Death", "ILD_Death")

# Names of recurrence-free (RF) states after ILD discontinuation by month
post_ild_state_names <- paste0("RF_postILD_m", seq_len(tx_window_months))

# All state names
state_names <- c("RF_on", post_ild_state_names, "RF_offtx_complete", "LRR", "DR", 
                 "OC_Death", "BC_Death", "ILD_Death")
n_states <- length(state_names)

# 2023 CDC life table for women
lt <- readr::read_csv("life_table_2023.csv") %>%
    transmute(age = as.numeric(stringr::str_extract(age, "^\\d+")),
              qx  = as.numeric(qx)) %>%
    arrange(age) %>%
    mutate(qx = pmin(qx, 1 - 1e-8), haz_mo = -log(1 - qx) / 12)

# Function to obtain other-cause mortality (background mortality) from life table
haz_OC <- function(age_years) {
    a <- floor(age_years)
    lt$haz_mo[match(a, lt$age)]
}

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

# 10-year DFS after LRR from CALOR trial (used to derive LRR to DR transition)
S10_LRR_DFS_base <- 0.70
S10_LRR_DFS_se   <- 0.09

## Helpers
# Clip probabilities between 0 and 1
clip01 <- function(p, eps = 1e-6) {
  pmin(pmax(p, eps), 1 - eps)
}

# Transform SE to logit scale
se_logit <- function(p, se_p, eps = 1e-6) {
  p <- clip01(p, eps)
  se_p / (p * (1 - p))
}

# Method of moments estimator for beta parameters from mean and SE
beta_params_from_mean_se <- function(mean, se) {
  k <- mean * (1 - mean) / se^2 - 1
  c(alpha = mean * k, beta = (1 - mean) * k)
}

# Compute log-likelihood contribution for one calibration target
ll_logitnorm <- function(p_hat, p_obs, se_obs) {
  p_hat <- clip01(p_hat)
  p_obs <- clip01(p_obs)
  dnorm(qlogis(p_obs),
        mean = qlogis(p_hat),
        sd   = se_logit(p_obs, se_obs),
        log  = TRUE
  )
}

# Obtain cycle-specific HR for RF to DR
get_hr_RF_DR <- function(t, arm, params) {
    
    # HR for T-DM1 is always 1 (reference)
    if (arm == "TDM1") {
        return(1)
    }
    
    # Full treatment effect
    hr_initial <- params$hr_RF_DR_TDXd_vs_TDM1

    # Evaluate time-varying HR at the midpoint of each monthly cycle
    time_years <- (t - 0.5) * params$delta

    # Treatment effect persists through year 3
    if (time_years <= params$hr_full_effect_years) {
        hr_t <- hr_initial
    }
    # Treatment effect has fully waned by year 10
    else if (time_years >= params$hr_wane_end_year) {
        hr_t <- 1
    }
    else {
        # Linear waning from year 3 to year 10
        waning_fraction <- (time_years - params$hr_full_effect_years) /
        (params$hr_wane_end_year - params$hr_full_effect_years)
        hr_t <- hr_initial + waning_fraction * (1 - hr_initial)
    }
    hr_t
}

# Build cycle-specific transition intensity matrix Q
build_Q_matrix <- function(t, arm, params) {
    
    # Other-cause mortality (from life table) at current age
    age <- params$age_init + (t - 1) * params$delta
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
    if (in_tx_window) {
    Q["RF_on", paste0("RF_postILD_m", t)] <- lam_ild_dis
    }
    
    # Fraction of benefit retained in the case of ILD discontinuation
    for (st in post_ild_state_names) {
        month_discontinued <- as.integer(sub("RF_postILD_m", "", st))
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
run_markov_model <- function(params) {
    
    # LRR to DR hazard
    params$haz_LRR_DR <- -log(params$S10_LRR_DFS) / 120
    
    # Any-grade ILD hazard
    lam_ild_any_TDM1_mo <- -log(1 - params$p_any_ILD_TDM1) / params$tx_window_months
    lam_ild_any_TDXd_mo <- -log(1 - params$p_any_ILD_TDXd) / params$tx_window_months
    
    # T-DM1 fatal and nonfatal ILD hazards
    lam_ild_fatal_TDM1_mo <- 0 # no fatal ILD for T-DM1 reported in DB-05
    lam_ild_nonfatal_TDM1_mo <- lam_ild_any_TDM1_mo
  
    # T-DM1 ILD discontinuation hazard
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
        # to eventually calculate DRFI
        cum_haz_DR <- 0
        
        for (t in seq_len(n_cycles)) {
            
            # Obtain transition intensity matrix
            Q_t <- build_Q_matrix(t = t, arm = arm, params = params)
            
            # Distribution halfway through the cycle for the DR risk set for a better
            # approximation
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
        out_drfi[arm] <- exp(-cum_haz_DR)
    }
    
    # Return outputs
    list(OS = out_os, iDFS = out_idfs, DRFI = out_drfi)
}

## Setup for IMIS calibration
# Specify DB-05 3-year OS, iDFS, and DRFI endpoints as the main targets
targets_endpoints <- tibble::tribble(~arm,   ~outcome, ~time, ~value, ~lower, ~upper,
                           "TDM1", "iDFS",   3,     0.837,  0.802,  0.867,
                           "TDM1", "DRFI",   3,     0.861,  0.825,  0.891,
                           "TDM1", "OS",     3,     0.957,  0.935,  0.972,
                           "TDXd", "iDFS",   3,     0.924,  0.897,  0.944,
                           "TDXd", "DRFI",   3,     0.939,  0.914,  0.957,
                           "TDXd", "OS",     3,     0.974,  0.958,  0.984) %>%
  mutate(se = (upper - lower) / (2 * 1.96),
         arm = factor(arm, levels = c("TDM1", "TDXd")),
         outcome = factor(outcome, levels = c("iDFS", "DRFI", "OS")))

# Specify LRR share informed by DB-05 as a separate target
target_lrr_share <- tibble::tibble(value = 0.10, lower = 0.08, upper = 0.12) %>%
    mutate(se = (upper - lower) / (2 * 1.96))

# Predict calibration endpoints
predict_all_endpoints_calib <- function(par_vec,
                                        base_params,
                                        arm,
                                        horizon_years_calib = 3) {
  # Base parameters
  params <- base_params
  params$horizon_years <- horizon_years_calib
  
  # Hazards
  params$haz_RF_DR_TDM1 <- par_vec["haz_rf_dr_base"]
  params$haz_RF_LRR <- (par_vec["lrr_share"] / (1 - par_vec["lrr_share"])) * 
    params$haz_RF_DR_TDM1 
  params$haz_DR_BCd <- par_vec["haz_dr_bcd_base"]
  params$hr_RF_DR_TDXd_vs_TDM1 <- par_vec["hr_rf_dr_TDXd_vs_TDM1"]
  
  # Obtain outputs
  out <- run_markov_model(params)
  list(OS = out$OS[arm], iDFS = out$iDFS[arm], DRFI = out$DRFI[arm])
}

# Parameters to be calibrated
v_param_names <- c("haz_rf_dr_base",
                   "lrr_share",
                   "haz_dr_bcd_base",
                   "hr_rf_dr_TDXd_vs_TDM1")
n_param <- length(v_param_names)

# Specify lower and upper bounds of independent uniform prior distributions
v_lb <- c(haz_rf_dr_base  = 1e-6,
          lrr_share = 0.01,
          haz_dr_bcd_base = 1e-4,
          hr_rf_dr_TDXd_vs_TDM1 = 0.20)

v_ub <- c(haz_rf_dr_base  = 0.01,
          lrr_share = 0.30,
          haz_dr_bcd_base = 0.15,
          hr_rf_dr_TDXd_vs_TDM1 = 1.00)

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
  hr_RF_DR_TDXd_vs_TDM1 = NA_real_,
  haz_DR_BCd            = NA_real_,
  
  # 10-year DFS after LRR
  S10_LRR_DFS = S10_LRR_DFS_base,
  
  # Other-cause mortality
  haz_OC = haz_OC,
  
  # ILD parameters
  p_any_ILD_TDM1 = p_any_ILD_TDM1_base,
  p_any_ILD_TDXd = p_any_ILD_TDXd_base,
  
  p_fatal_ILD_all_TDXd = p_fatal_ILD_all_TDXd_calib,
  prop_fatal_ILD_TDM1  = 0,
  prop_fatal_ILD_TDXd  = prop_fatal_ILD_TDXd_calib,
  
  prop_dis_ILD_TDM1 = prop_dis_ILD_TDM1,
  prop_dis_ILD_TDXd = prop_dis_ILD_TDXd_calib,
  
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
  llik <- numeric(n_samp)
  
  for (j in seq_len(n_samp)) {
    
    pred_tdm <- predict_all_endpoints_calib(
      par_vec = v_params[j, ],
      base_params = l_params_calib,
      arm = "TDM1"
    )
    
    pred_tdx <- predict_all_endpoints_calib(
      par_vec = v_params[j, ],
      base_params = l_params_calib,
      arm = "TDXd"
    )
    
    pred <- list(TDM1 = pred_tdm, TDXd = pred_tdx)
    
    for (k in seq_len(nrow(targets_endpoints))) {
      arm_k <- as.character(targets_endpoints$arm[k])
      outcome_k <- as.character(targets_endpoints$outcome[k])
      
      # Likelihood contribution from main targets (DB-05 3-year endpoints)
      llik[j] <- llik[j] +
        ll_logitnorm(
          p_hat  = pred[[arm_k]][[outcome_k]],
          p_obs  = targets_endpoints$value[k],
          se_obs = targets_endpoints$se[k]
        )
    }
    # Likelihood contribution from LRR share target
    llik[j] <- llik[j] +
      ll_logitnorm(
        p_hat  = v_params[j, "lrr_share"],
        p_obs  = target_lrr_share$value,
        se_obs = target_lrr_share$se
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

# # Model output at posterior mean
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
  "haz_rf_dr_base",
  "lrr_share",
  "haz_dr_bcd_base",
  "hr_rf_dr_TDXd_vs_TDM1"
)

param_labels <- c(
  haz_rf_dr_base = "RF→DR hazard\n(T-DM1)",
  lrr_share = "LRR share\nof first recurrences",
  haz_dr_bcd_base = "DR→BC death\nhazard",
  hr_rf_dr_TDXd_vs_TDM1 = "HR for RF→DR\n(T-DXd vs T-DM1)"
)

post_draws_plot <- post_draws %>%
  dplyr::select(dplyr::all_of(param_order))

# Pairs plot
gg_post_pairs_corr <- GGally::ggpairs(
  post_draws_plot,
  upper = list(
    continuous = GGally::wrap("cor", color = "black", size = 4)
  ),
  diag = list(
    continuous = GGally::wrap("densityDiag", alpha = 0.6)
  ),
  lower = list(
    continuous = GGally::wrap("points", alpha = 0.25, size = 0.6)
  ),
  columnLabels = unname(param_labels[param_order])
) +
  theme_bw(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x  = element_text(size = 8.5),
    axis.text.y  = element_text(size = 8.5),
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text = element_text(hjust = 0, lineheight = 0.9)
  )

print(gg_post_pairs_corr)

# Sample from prior for comparison
m_samp_prior <- sample_prior(nrow(post_draws_plot))
df_samp_prior <- as.data.frame(m_samp_prior)
colnames(df_samp_prior) <- param_order

# Prior vs. posterior plot
df_prior_post <- dplyr::bind_rows(
  df_samp_prior %>%
    mutate(Distribution = "Prior"),
  post_draws_plot %>%
    mutate(Distribution = "Posterior")
) %>%
  tidyr::pivot_longer(
    cols = all_of(param_order),
    names_to = "Parameter",
    values_to = "value"
  ) %>%
  mutate(
    Parameter = factor(Parameter, levels = param_order, labels = unname(param_labels)),
    Distribution = factor(Distribution, levels = c("Prior", "Posterior"))
  )

gg_prior_post <- ggplot(
  df_prior_post,
  aes(x = value, y = after_stat(density), fill = Distribution)
) +
  geom_density(alpha = 0.45) +
  facet_wrap(~ Parameter, scales = "free", ncol = 2) +
  theme_bw(base_size = 15) +
  theme(
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  labs(title = "Prior versus posterior distributions")

print(gg_prior_post)

## Internal validation: calibration plot
df_targets <- targets_endpoints %>%
  transmute(
    Type    = "Target",
    Arm     = as.character(arm),
    Outcome = outcome,
    time    = time,
    value   = value,
    lb      = lower,
    ub      = upper
  )
df_targets$Outcome <- factor(df_targets$Outcome, levels = c("iDFS", "DRFI", "OS"))

# Predict endpoints for one posterior draw and arm
pred_one <- function(par_row, arm) {
  
  out <- predict_all_endpoints_calib(
    par_vec     = par_row,
    base_params = l_params_calib,
    arm         = arm
  )
  
  c(
    iDFS = as.numeric(out$iDFS),
    DRFI = as.numeric(out$DRFI),
    OS   = as.numeric(out$OS)
  )
}

post_param_draws <- post_draws

# Generate endpoint predictions across all posterior draws for each arm
pred_rows <- lapply(levels(targets_endpoints$arm), function(arm) {
  
  pred_mat <- t(apply(post_param_draws, 1,
                      function(par_row) {
                        pred_one(par_row = par_row, arm = arm)
                      }))
  
  colnames(pred_mat) <- c("iDFS", "DRFI", "OS")
  
  data.frame(Arm  = arm,
             pred_mat,
             check.names = FALSE)
})

df_pred_wide <- bind_rows(pred_rows)
df_pred_long <- pivot_longer(df_pred_wide,
                             cols = c("iDFS", "DRFI", "OS"),
                             names_to = "Outcome",
                             values_to = "value")

df_pred_long$Type <- "Model"
df_pred_long$time <- 3
df_pred_long$Outcome <- factor(df_pred_long$Outcome, levels = c("iDFS", "DRFI", "OS"))

df_model_sum <- df_pred_long %>%
  group_by(Type, Arm, Outcome, time) %>%
  summarise(mean_value = mean(value, na.rm = TRUE),
            lb = as.numeric(quantile(value, probs = 0.025, na.rm = TRUE, names = FALSE)),
            ub = as.numeric(quantile(value, probs = 0.975, na.rm = TRUE, names = FALSE)),
            n_draws = n(),
            .groups = "drop") %>%
  rename(value = mean_value)

df_plot <- bind_rows(df_targets, df_model_sum)
df_plot$Outcome <- factor(df_plot$Outcome, levels = c("iDFS", "DRFI", "OS"))
df_plot$Type    <- factor(df_plot$Type, levels = c("Model", "Target"))

df_plot_int <- df_plot %>%
  select(Type, Arm, Outcome, time, lb, ub)

ggplot(df_plot_int, aes(x = Outcome, ymin = lb, ymax = ub, group = Type)) +
  geom_errorbar(aes(color = Type),
                position = position_dodge(width = 0.55),
                width = 0.18,
                linewidth = 0.9) +
  facet_wrap(~ Arm) +
  scale_y_continuous("3-year probability", limits = c(0, 1)) +
  theme_bw(base_size = 18) +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c(Model = "steelblue3", Target = "black"))

# Posterior-mean calibrated parameters
haz_RF_DR_TDM1_post <- as.numeric(v_calib_post_mean["haz_rf_dr_base"])
lrr_share_post      <- as.numeric(v_calib_post_mean["lrr_share"])
haz_RF_LRR_post     <- 
  (lrr_share_post / (1 - lrr_share_post)) * haz_RF_DR_TDM1_post
hr_RF_DR_TDXd_vs_TDM1_post <- as.numeric(v_calib_post_mean["hr_rf_dr_TDXd_vs_TDM1"])
haz_DR_BCd_post        <- as.numeric(v_calib_post_mean["haz_dr_bcd_base"])

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
  lrr_share             = lrr_share_post,
  hr_RF_DR_TDXd_vs_TDM1 = hr_RF_DR_TDXd_vs_TDM1_post,
  haz_DR_BCd            = haz_DR_BCd_post,
  
  # 10-year DFS after LRR
  S10_LRR_DFS = S10_LRR_DFS_base,
  
  # Other-cause mortality
  haz_OC = haz_OC,
  
  # ILD parameters
  p_any_ILD_TDM1 = p_any_ILD_TDM1_base,
  p_any_ILD_TDXd = p_any_ILD_TDXd_base,
  
  p_fatal_ILD_all_TDXd = p_fatal_ILD_all_TDXd_base,
  
  prop_fatal_ILD_TDM1 = prop_fatal_ILD_TDM1,
  prop_fatal_ILD_TDXd = prop_fatal_ILD_TDXd,
  
  prop_dis_ILD_TDM1 = prop_dis_ILD_TDM1,
  prop_dis_ILD_TDXd = prop_dis_ILD_TDXd,
  
  prop_dis_given_nonfatal_ILD_TDM1 = prop_dis_given_nonfatal_ILD_TDM1,
  prop_dis_given_nonfatal_ILD_TDXd = prop_dis_given_nonfatal_ILD_TDXd,
  
  # A-priori-specified fractional benefit following ILD discontinuation
  benefit_frac_postILD = c(
    0.28, 0.49, 0.64, 0.75, 0.83,
    0.89, 0.94, 0.97, 0.99, 1.00
  )
)
saveRDS(l_params_all, "l_params_all_basecase.rds")

## Map DR hazard to RCB categories and ER subgroups
# 5-year scaled T-DM1 (baseline) EFS from Yau et al
RCB_EFS_scaled <- tibble::tribble(
  ~RCB,      ~ER_status,    ~EFS_5y_TDM1,
  "RCB-I",   "ER-/HER2+",   1 - 0.084,
  "RCB-II",  "ER-/HER2+",   1 - 0.221,
  "RCB-III", "ER-/HER2+",   1 - 0.241,
  "RCB-I",   "ER+/HER2+",   1 - 0.050,
  "RCB-II",  "ER+/HER2+",   1 - 0.138,
  "RCB-III", "ER+/HER2+",   1 - 0.283
) %>%
  mutate(
    RCB = factor(RCB, levels = c("RCB-I", "RCB-II", "RCB-III")),
    ER_status = factor(ER_status, levels = c("ER-/HER2+", "ER+/HER2+"))
  )

# Grid search to find multiplier applied to DR hazard that most closely reproduce
# scaled EFS estimates
dr_multiplier_grid <- seq(0.01, 10, length.out = 1000)

# Store base parameters
l_params_rcb_base <- l_params_all
l_params_rcb_base$horizon_years <- 5
haz_RF_DR_TDM1_base <- l_params_all$haz_RF_DR_TDM1
haz_RF_LRR_base     <- l_params_all$haz_RF_LRR

# Get modeled 5-year T-DM1 EFS (corresponding to iDFS in our model) for one DR multiplier
get_EFS_5y_TDM1 <- function(dr_multiplier, params_base) {
  
  params <- params_base
  
  # Apply multiplier only to RF to DR hazard
  params$haz_RF_DR_TDM1 <- haz_RF_DR_TDM1_base * dr_multiplier
  
  # Keep RF to LRR fixed
  params$haz_RF_LRR <- haz_RF_LRR_base
  
  out <- run_markov_model(params)
  
  as.numeric(out$iDFS["TDM1"])
}

# Evaluate multiplier grid
grid_results <- tibble::tibble(dr_multiplier = dr_multiplier_grid) %>%
    mutate(EFS_5y_pred = purrr::map_dbl(dr_multiplier, get_EFS_5y_TDM1,
                                        params_base = l_params_rcb_base))

RCB_DR_multipliers_grid <- RCB_EFS_scaled %>%
  rowwise() %>%
  mutate(
    best_idx = which.min(abs(grid_results$EFS_5y_pred - EFS_5y_TDM1)),
    dr_multiplier = grid_results$dr_multiplier[best_idx],
    EFS_5y_pred = grid_results$EFS_5y_pred[best_idx],
    abs_error = abs(EFS_5y_pred - EFS_5y_TDM1),
    haz_RF_DR_TDM1 = haz_RF_DR_TDM1_base * dr_multiplier,
    haz_RF_LRR = haz_RF_LRR_base
  ) %>%
  ungroup() %>%
  select(-best_idx)
RCB_DR_multipliers_grid

# Scenarios 
base_case <- function(params) {
  params
}

weaker_hr <- function(params) {
  params$hr_RF_DR_TDXd_vs_TDM1 <- 0.80
  params
}

higher_ild_incidence <- function(params) {
  params$p_any_ILD_TDXd <- min(1 - 1e-12, params$p_any_ILD_TDXd * 2)
  
  # For reporting purposes only
  params$p_fatal_ILD_all_TDXd <-
    params$p_any_ILD_TDXd * params$prop_fatal_ILD_TDXd
  
  params
}

combined <- function(params) {
  params <- weaker_hr(params)
  params <- higher_ild_incidence(params)
  params
}

# Named scenarios for deterministic analyses
scenarios <- list(
  "Base case" = base_case,
  "HR 0.80" = weaker_hr,
  "Any-grade ILD incidence x2" = higher_ild_incidence,
  "Combined" = combined
)
scenario_order <- names(scenarios)

## Deterministic analyses
# Evaluate each scenario across RCB categories and ER subgroups
eval_os_RCB <- function(params,
                        scenario_name,
                        scenario_modify = function(params) params) {
  
  results <- vector("list", nrow(RCB_DR_multipliers_grid))
  
  # Loop over all DR multipliers (corresponding to RCB categories/ER subgroups)
  for (i in seq_len(nrow(RCB_DR_multipliers_grid))) {
    
    multiplier_row <- RCB_DR_multipliers_grid[i, ]
    
    current_params <- params
    
    # Set RCB/ER-specific baseline T-DM1 DR hazard
    current_params$haz_RF_DR_TDM1 <- multiplier_row$haz_RF_DR_TDM1
    
    # Keep RF to LRR fixed at calibrated base value
    current_params$haz_RF_LRR <- multiplier_row$haz_RF_LRR
    
    # Apply scenario after setting subgroup-specific hazards
    current_params <- scenario_modify(current_params)
    
    # Run Markov model
    outcomes_10y <- run_markov_model(params = current_params)
    
    results[[i]] <- tibble::tibble(
        scenario = scenario_name,
        ER_status = multiplier_row$ER_status,
        RCB = multiplier_row$RCB,
        OS_10y_TDM1 = as.numeric(outcomes_10y$OS["TDM1"]),
        OS_10y_TDXd = as.numeric(outcomes_10y$OS["TDXd"]),
        OS_diff_percent_TDXd_minus_TDM1 =
            100 * (as.numeric(outcomes_10y$OS["TDXd"]) - 
                       as.numeric(outcomes_10y$OS["TDM1"])
            )
    )
  }
  
  dplyr::bind_rows(results)
}

# Run deterministic analysis for each scenario
df_os_RCB_list <- list()
result_index <- 1

for (s in seq_along(scenarios)) {
    scenario_name <- names(scenarios)[s]
    scenario_modify <- scenarios[[s]]

    cat(
      "Running deterministic analysis:",
      scenario_name,
      "\n"
    )

    combined_modify <- function(params) {
      params <- scenario_modify(params)
      params
    }

    df_os_RCB_list[[result_index]] <- eval_os_RCB(
      params = l_params_all,
      scenario_name = scenario_name,
      scenario_modify = combined_modify
    )

    result_index <- result_index + 1
  }

# Summarize deterministic analysis
df_os_RCB <- dplyr::bind_rows(df_os_RCB_list) %>%
    mutate(
        scenario = factor(
            scenario,
            levels = scenario_order
        )
    ) %>%
    arrange(ER_status, RCB, scenario)

table_os_RCB <- df_os_RCB %>%
    mutate(
        OS_diff_percent_TDXd_minus_TDM1 =
            round(OS_diff_percent_TDXd_minus_TDM1, 2)
    ) %>%
    select(
        scenario,
        ER_status,
        RCB,
        OS_diff_percent_TDXd_minus_TDM1
    ) %>%
    arrange(scenario, ER_status, RCB)

table_os_RCB

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

# Function to sample calibrated parameters from posterior draws
sample_calibrated_params <- function(n_sim, post_draws) {
  idx <- sample(seq_len(nrow(post_draws)), size = n_sim, replace = TRUE)
  as.data.frame(post_draws[idx, , drop = FALSE])
}

# Draw calibrated model parameters from IMIS posterior
set.seed(1)
calib_param_draws <- sample_calibrated_params(
  n_sim = n_sim,
  post_draws = post_draws
)

# Any-grade ILD probability over treatment window: T-DXd
p_any_ILD_TDXd_draw <- draw_prob_logitnormal(
  mu = 0.096,
  lower_95 = 0.05,
  upper_95 = 0.21,
  n_sim = n_sim
)

# Any-grade ILD probability over treatment window: T-DM1
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

# Obtain parameters for beta distribution for DFS after LRR
beta_params <- beta_params_from_mean_se(
    mean = S10_LRR_DFS_base,
    se   = S10_LRR_DFS_se
)
alpha_LRR_DFS <- beta_params["alpha"]
beta_LRR_DFS  <- beta_params["beta"]

# 10-year DFS after LRR
S10_LRR_DFS_draw <- rbeta(
  n = n_sim,
  shape1 = alpha_LRR_DR,
  shape2 = beta_LRR_DR
)

# Histograms of PSA draws
hist(calib_param_draws$haz_rf_dr_base)
hist(calib_param_draws$lrr_share)
hist(calib_param_draws$haz_dr_bcd_base)
hist(calib_param_draws$hr_rf_dr_TDXd_vs_TDM1)
hist(p_any_ILD_TDXd_draw)
hist(p_fatal_ILD_all_TDXd_draw)
hist(S10_LRR_DFS_draw)

# Run probabilistic analysis
eval_os_RCB_PSA <- function(
    params,
    scenario_name,
    scenario_modify,
    n_sim = 1000,
    calib_param_draws,
    p_any_ILD_TDM1_draw,
    p_any_ILD_TDXd_draw,
    p_fatal_ILD_all_TDXd_draw,
    S10_LRR_DFS_draw
) {
  
  results <- vector("list", nrow(RCB_DR_multipliers_grid))
  
  for (i in seq_len(nrow(RCB_DR_multipliers_grid))) {
    
    multiplier_row <- RCB_DR_multipliers_grid[i, ]
    
    sim_results <- vector("list", n_sim)
    
    for (s in seq_len(n_sim)) {
      
      current_params <- params
      
      # PSA draw: calibrated base parameters
      haz_RF_DR_TDM1_base_s <- as.numeric(calib_param_draws$haz_rf_dr_base[s])
      lrr_share_s           <- as.numeric(calib_param_draws$lrr_share[s])
      haz_DR_BCd_s          <- as.numeric(calib_param_draws$haz_dr_bcd_base[s])
      hr_RF_DR_s            <- as.numeric(calib_param_draws$hr_rf_dr_TDXd_vs_TDM1[s])
      
      # Apply RCB/ER-specific multiplier only to RF to DR hazard
      current_params$haz_RF_DR_TDM1 <-
        haz_RF_DR_TDM1_base_s * as.numeric(multiplier_row$dr_multiplier)
      
      # LRR hazard is derived from calibrated posterior draw
      current_params$lrr_share <- lrr_share_s
      
      current_params$haz_RF_LRR <-
        (lrr_share_s / (1 - lrr_share_s)) * haz_RF_DR_TDM1_base_s
      
      # PSA draw: DR to BC death hazard
      current_params$haz_DR_BCd <- haz_DR_BCd_s
      
      # PSA draw: 10-year DFS after LRR (used to derive LRR to DR transition)
      current_params$S10_LRR_DFS <- S10_LRR_DFS_draw[s]
      
      # PSA draw: HR for RF to DR, T-DXd vs T-DM1
      current_params$hr_RF_DR_TDXd_vs_TDM1 <- hr_RF_DR_s
      
      # PSA draw: any-grade ILD probability for T-DXd
      current_params$p_any_ILD_TDXd <- p_any_ILD_TDXd_draw[s]
      
      # PSA draw: any-grade ILD probability for T-DM1
      current_params$p_any_ILD_TDM1 <- p_any_ILD_TDM1_draw[s]
      
      # PSA draw: overall fatal ILD probability among all T-DXd-assigned patients
      current_params$p_fatal_ILD_all_TDXd <- p_fatal_ILD_all_TDXd_draw[s]
      
      # Ensure fatal ILD remains a subset of any-grade ILD before applying scenarios
      current_params$p_fatal_ILD_all_TDXd <- min(
        current_params$p_fatal_ILD_all_TDXd,
        0.99 * current_params$p_any_ILD_TDXd
      )
      
      # Convert sampled overall fatal ILD probability to fatal proportion among ILD
      current_params$prop_fatal_ILD_TDXd <-
        current_params$p_fatal_ILD_all_TDXd /
        current_params$p_any_ILD_TDXd
      
      # Apply PSA scenarios
      current_params <- scenario_modify(current_params)
      
      # Recalculate ILD quantities after scenario modification
      p_any_ILD_TDXd <- current_params$p_any_ILD_TDXd
      
      p_fatal_ILD_all_TDXd <- min(
        current_params$p_fatal_ILD_all_TDXd,
        0.99 * p_any_ILD_TDXd
      )
      
      p_nonfatal_ILD_all_TDXd <-
        p_any_ILD_TDXd - p_fatal_ILD_all_TDXd
      
      # Fixed probability of discontinuation among nonfatal ILD cases
      prop_dis_given_nonfatal_ILD_TDXd <-
        current_params$prop_dis_given_nonfatal_ILD_TDXd
      
      # Overall probability of ILD discontinuation among T-DXd-assigned patients
      p_dis_ILD_all_TDXd <-
        p_nonfatal_ILD_all_TDXd * prop_dis_given_nonfatal_ILD_TDXd
      
      # Convert to proportions among any-grade ILD cases
      prop_fatal_ILD_TDXd <-
        p_fatal_ILD_all_TDXd / p_any_ILD_TDXd
      prop_dis_ILD_TDXd <-
        p_dis_ILD_all_TDXd / p_any_ILD_TDXd
      
      # Assign ILD parameters used by the model
      current_params$p_any_ILD_TDXd <- p_any_ILD_TDXd
      current_params$p_fatal_ILD_all_TDXd <- p_fatal_ILD_all_TDXd
      current_params$prop_fatal_ILD_TDXd <- prop_fatal_ILD_TDXd
      current_params$prop_dis_ILD_TDXd <- prop_dis_ILD_TDXd
      current_params$prop_dis_given_nonfatal_ILD_TDXd <-
        prop_dis_given_nonfatal_ILD_TDXd
      
      # Run Markov model
      outcomes_10y <- run_markov_model(params = current_params)
      
      sim_results[[s]] <- tibble::tibble(
          sim = s,
          scenario = scenario_name,
          ER_status = multiplier_row$ER_status,
          RCB = multiplier_row$RCB,
          OS_10y_TDM1 = as.numeric(outcomes_10y$OS["TDM1"]),
          OS_10y_TDXd = as.numeric(outcomes_10y$OS["TDXd"]),
          OS_diff_percent_TDXd_minus_TDM1 =
              100 * (
                  as.numeric(outcomes_10y$OS["TDXd"]) -
                      as.numeric(outcomes_10y$OS["TDM1"])
              )
      )
    }
    
    results[[i]] <- dplyr::bind_rows(sim_results)
  }
  
  dplyr::bind_rows(results)
}

# Conduct probabilistic analyses across scenarios
df_os_RCB_psa_list <- list()
result_index <- 1

for (s in seq_along(psa_scenarios)) {
    scenario_name <- names(psa_scenarios)[s]
    scenario_modify <- psa_scenarios[[s]]
    
    cat(
        "Running PSA:",
        scenario_name,
        "\n"
    )
    
    combined_modify <- function(params) {
        params <- scenario_modify(params)
        params
    }
    
    df_os_RCB_psa_list[[result_index]] <- eval_os_RCB_PSA(
        params = l_params_all,
        scenario_name = scenario_name,
        scenario_modify = combined_modify,
        n_sim = n_sim,
        calib_param_draws = calib_param_draws,
        p_any_ILD_TDM1_draw = p_any_ILD_TDM1_draw,
        p_any_ILD_TDXd_draw = p_any_ILD_TDXd_draw,
        p_fatal_ILD_all_TDXd_draw = p_fatal_ILD_all_TDXd_draw,
        S10_LRR_DFS_draw = S10_LRR_DFS_draw
    ) 
    
    result_index <- result_index + 1
}

# Organize probabilistic analysis results into dataframe
df_os_RCB_psa <- dplyr::bind_rows(df_os_RCB_psa_list) %>%
    mutate(
        scenario = factor(
            scenario,
            levels = psa_scenario_order
        )
    ) %>%
    arrange(ER_status, RCB, scenario, sim)

# Summarize results of probabilistic analysis
table_os_RCB_psa <- df_os_RCB_psa %>%
    group_by(scenario, ER_status, RCB) %>%
    summarise(
        mean_OS_diff_percent_TDXd_minus_TDM1 =
            round(mean(OS_diff_percent_TDXd_minus_TDM1), 2),
        
        median_OS_diff_percent_TDXd_minus_TDM1 =
            round(median(OS_diff_percent_TDXd_minus_TDM1), 2),
        
        lower_OS_diff_percent_TDXd_minus_TDM1 =
            round(quantile(OS_diff_percent_TDXd_minus_TDM1, 0.025), 2),
        
        upper_OS_diff_percent_TDXd_minus_TDM1 =
            round(quantile(OS_diff_percent_TDXd_minus_TDM1, 0.975), 2),
        
        pct_favors_TDXd =
            round(100 * mean(OS_diff_percent_TDXd_minus_TDM1 > 0), 1),
        
        pct_favors_TDM1 =
            round(100 * mean(OS_diff_percent_TDXd_minus_TDM1 < 0), 1),
        
        .groups = "drop"
    ) %>%
    arrange(scenario, ER_status, RCB)

table_os_RCB_psa

# Save results
saveRDS(df_os_RCB_psa, "df_os_RCB_PSA.rds")
saveRDS(table_os_RCB_psa, "table_os_RCB_PSA.rds")

## One-way ILD threshold analysis
threshold_os_RCB <- function(params,
                             scenario_modify = function(params) params,
                             multiplier_grid = RCB_DR_multipliers_grid) {
  
  purrr::map_dfr(seq_len(nrow(multiplier_grid)), function(i) {
    
    multiplier_row <- multiplier_grid[i, ]
    
    current_params <- params
    current_params$haz_RF_DR_TDM1 <- multiplier_row$haz_RF_DR_TDM1
    current_params$haz_RF_LRR <- multiplier_row$haz_RF_LRR
    
    # Apply threshold after subgroup-specific hazards are set
    current_params <- scenario_modify(current_params)
    
    outcomes_10y <- run_markov_model(params = current_params)
    
    tibble::tibble(
      ER_status = multiplier_row$ER_status,
      RCB = multiplier_row$RCB,
      dr_multiplier = multiplier_row$dr_multiplier,
      OS_diff_percent_TDXd_minus_TDM1 =
        100 * (as.numeric(outcomes_10y$OS["TDXd"]) - as.numeric(outcomes_10y$OS["TDM1"]))
    )
  })
}

# Apply multiplier to T-DXd any-grade ILD incidence
modify_param_ild <- function(params, x) {
    params$p_any_ILD_TDXd <- min(1 - 1e-12, params$p_any_ILD_TDXd * x)
  
  # For reporting purposes only
  params$p_fatal_ILD_all_TDXd <-
    params$p_any_ILD_TDXd *
    params$prop_fatal_ILD_TDXd
  
  params
}

# Maximum possible multiplier before any-grade ILD reaches 1
max_ild_multiplier <- (1 - 1e-12) / l_params_all$p_any_ILD_TDXd
ild_grid <- seq(from = 1, to = max_ild_multiplier, length.out = 500)


# Evaluate the ILD multiplier grid 
os_diff_ild <- purrr::map_dfr(ild_grid, function(x) {
    threshold_os_RCB(params = l_params_all, scenario_modify = function(params) {
        
        # Apply ILD multiplier
        params <- modify_param_ild(params = params, x = x)
        }) %>%
      mutate(threshold = x)
    }
)


# Find first ILD multiplier such that 10-year OS difference is negative (T-DM1 favored)
first_neg_diff_ild <- os_diff_ild %>%
    arrange(ER_status, RCB, threshold) %>%
    group_by(ER_status, RCB) %>%
    filter(OS_diff_percent_TDXd_minus_TDM1 < 0) %>%
    slice(1) %>%
    ungroup()
first_neg_diff_ild

## Plots
# Desired top-to-bottom subgroup ordering
subgroup_order_top_to_bottom <- c(
    "RCB-I ER+",
    "RCB-I ER-",
    "RCB-II ER+",
    "RCB-II ER-",
    "RCB-III ER+",
    "RCB-III ER-"
)
subgroup_levels_for_plot <- rev(subgroup_order_top_to_bottom)

# Colors used for the deterministic plot
subgroup_colors <- c(
    `FALSE` = "#1D9E75",  # Other subgroups
    `TRUE`  = "#EF9F27"   # RCB-I ER+
)

# Colors used for the PSA distributions
distribution_side_colors <- c(
    `TRUE`  = "#F2C94C",  # OS difference < 0: favors T-DM1
    `FALSE` = "#1D9E75"   # OS difference >= 0: favors T-DXd
)

# Shared axis labels
deterministic_x_label <- paste0(
    "Absolute difference in 10-year OS, percentage points\n",
    "Positive values favor T-DXd"
)
psa_x_label <- paste0(
    "Difference in 10-year OS, percentage points ",
    "(T-DXd \u2212 T-DM1)"
)

# Shared formatting
format_subgroups <- function(data) {
    data %>%
        mutate(
            er_short = recode(
                as.character(ER_status),
                "ER-/HER2+" = "ER-",
                "ER+/HER2+" = "ER+"
            ),
            subgroup = factor(
                paste(RCB, er_short),
                levels = subgroup_levels_for_plot
            )
        )
}

# Deterministic dumbbell plot
det_plot_df <- df_os_RCB %>%
    filter(scenario %in% c("Base case", "Combined")) %>%
    format_subgroups() %>%
    mutate(equivocal = RCB == "RCB-I" & ER_status == "ER+/HER2+") %>%
    select(
        RCB, ER_status, er_short, subgroup, equivocal, scenario,
        OS_diff_percent_TDXd_minus_TDM1
    ) %>%
    pivot_wider(
        id_cols = c(RCB, ER_status, er_short, subgroup, equivocal),
        names_from = scenario,
        values_from = OS_diff_percent_TDXd_minus_TDM1
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
        y = NULL
    ) +
    
    theme_minimal(base_size = 12) +
    
    theme(
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank()
    )

gg_det_dumbbell

# Base-case PSA ridge plot
psa_plot_df <- df_os_RCB_psa %>%
    filter(scenario == "Base case PSA") %>%
    format_subgroups() %>%
    mutate(os_diff = OS_diff_percent_TDXd_minus_TDM1)

pct_tdm1 <- psa_plot_df %>%
    group_by(subgroup) %>%
    summarise(
        p_tdm1 = mean(os_diff < 0, na.rm = TRUE),
        .groups = "drop"
    )

x_limits <- range(psa_plot_df$os_diff, na.rm = TRUE)
x_limits <- x_limits + c(-1, 1) * 0.08 * diff(x_limits)

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
            label = sprintf("%.0f%% <0", 100 * p_tdm1)
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
    
    coord_cartesian(xlim = x_limits) +
    
    labs(
        title = "Base-case probabilistic sensitivity analysis",
        x = psa_x_label,
        y = NULL
    ) +
    
    ggridges::theme_ridges(center_axis_labels = TRUE) +
    
    theme(
        legend.position = "bottom",
        plot.title = element_text(face = "bold"),
        panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3),
        plot.margin = margin(5.5, 15, 5.5, 5.5)
    )
gg_psa_ridges
