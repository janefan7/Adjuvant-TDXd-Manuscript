# Load packages
library(msm)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(readr)
library(stringr)
library(patchwork)
library(grid)

# Global model parameters
tx_window_months <- 10
age_init <- 50
delta <- 1/12
horizon_years <- 10

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

# Build cycle-specific transition intensity matrix Q
build_Q_matrix <- function(t, arm, params) {
  
  age <- params$age_init + (t - 1) * params$delta
  in_tx_window <- t <= params$tx_window_months
  
  lam_oc_mo <- params$haz_OC(age)
  hr_full <- if (arm == "TDXd") params$hr_RF_DR_TDXd_vs_TDM1 else 1
  
  lam_ild_dis <- if (in_tx_window && arm == "TDXd") {
    params$lam_ild_dis_TDXd_mo
  } else if (in_tx_window && arm == "TDM1") {
    params$lam_ild_dis_TDM1_mo
  } else {
    0
  }
  
  lam_ild_fatal <- if (in_tx_window && arm == "TDXd") {
    params$lam_ild_fatal_TDXd_mo
  } else if (in_tx_window && arm == "TDM1") {
    params$lam_ild_fatal_TDM1_mo
  } else {
    0
  }
  
  Q <- matrix(
    0,
    nrow = n_states,
    ncol = n_states,
    dimnames = list(state_names, state_names)
  )
  
  Q["RF_on", "LRR"]       <- params$haz_RF_LRR
  Q["RF_on", "DR"]        <- hr_full * params$haz_RF_DR_TDM1
  Q["RF_on", "OC_Death"]  <- lam_oc_mo
  Q["RF_on", "ILD_Death"] <- lam_ild_fatal
  
  if (in_tx_window) {
    Q["RF_on", paste0("RF_postILD_m", t)] <- lam_ild_dis
  }
  
  for (st in post_ild_state_names) {
    
    month_discontinued <- as.integer(sub("RF_postILD_m", "", st))
    
    if (arm == "TDXd") {
      hr_postILD <- hr_full ^ params$benefit_frac_postILD[month_discontinued]
    } else {
      hr_postILD <- 1
    }
    
    Q[st, "LRR"]      <- params$haz_RF_LRR
    Q[st, "DR"]       <- hr_postILD * params$haz_RF_DR_TDM1
    Q[st, "OC_Death"] <- lam_oc_mo
  }
  
  Q["RF_offtx_complete", "LRR"]      <- params$haz_RF_LRR
  Q["RF_offtx_complete", "DR"]       <- hr_full * params$haz_RF_DR_TDM1
  Q["RF_offtx_complete", "OC_Death"] <- lam_oc_mo
  
  Q["LRR", "DR"]       <- params$haz_LRR_DR
  Q["LRR", "OC_Death"] <- lam_oc_mo
  
  Q["DR", "BC_Death"] <- params$haz_DR_BCd
  Q["DR", "OC_Death"] <- lam_oc_mo
  
  diag(Q) <- -rowSums(Q)
  Q
}

# Markov model for 10-year OS
run_markov_model <- function(params) {
  n_cycles <- as.integer(round(params$horizon_years / params$delta))

  # LRR to DR hazard
  params$haz_LRR_DR <- -log(params$S10_LRR_DFS) / 120

  # Any-grade ILD hazard
  lam_ild_any_TDM1_mo <- -log(1 - params$p_any_ILD_TDM1) / params$tx_window_months
  lam_ild_any_TDXd_mo <- -log(1 - params$p_any_ILD_TDXd) / params$tx_window_months
  
  # T-DM1 fatal and nonfatal ILD hazards
  lam_ild_fatal_TDM1_mo <- 0 # no fatal ILD reported in DB-05
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

  # Obtain cohort trace by arm
  for (arm in strategy_names) {
    
    m_trace <- matrix(0, nrow = n_cycles + 1, ncol = n_states, 
                      dimnames = list(0:n_cycles, state_names))
    
    m_trace[1, ] <- v_m_init
    
    for (t in seq_len(n_cycles)) {
      Q_t <- build_Q_matrix(t = t, arm = arm, params = params)
      P_t <- msm::MatrixExp(Q_t)
      m_trace[t + 1, ] <- as.numeric(m_trace[t, ] %*% P_t)
      
      # At the end of the treatment window remaining RF_on patients move to RF_offtx_complete
      if (t == params$tx_window_months) {
        m_trace[t + 1, "RF_offtx_complete"] <- m_trace[t + 1, "RF_offtx_complete"] + 
          m_trace[t + 1, "RF_on"]
        m_trace[t + 1, "RF_on"] <- 0
      }
    }
    
    # Obtain arm-specific output
    out_os[arm] <- 1 - sum(m_trace[n_cycles + 1, death_states])
  }
  out_os
}


# 2023 CDC life table for women
lt <- readr::read_csv("life_table_2023.csv") %>%
  transmute(age = as.numeric(stringr::str_extract(age, "^\\d+")),
            qx  = as.numeric(qx)) %>%
  arrange(age) %>%
  mutate(qx = pmin(qx, 1 - 1e-8),
         haz_mo = -log(1 - qx) / 12)

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

## ILD inputs in deterministic analyses, base case ##
# Overall fatal ILD probability among all T-DXd-assigned patients
p_fatal_ILD_all_TDXd_base <- 0.01

# Fatal ILD among any-grade ILD cases, T-DXd
prop_fatal_ILD_TDXd <- p_fatal_ILD_all_TDXd_base / p_any_ILD_TDXd_base

# DB-05 discontinuation probability among nonfatal ILD cases
prop_dis_given_nonfatal_ILD_TDXd <- dis_TDXd_nonfatal / (any_ild_TDXd - fatal_TDXd)

# Nonfatal ILD discontinuation among any-grade ILD cases, T-DXd
prop_dis_ILD_TDXd <- 
  (1 - prop_fatal_ILD_TDXd) * prop_dis_given_nonfatal_ILD_TDXd

# T-DM1 fatal ILD assumed 0 in base case
prop_fatal_ILD_TDM1 <- 0

# Discontinuation among any-grade ILD cases, T-DM1
prop_dis_ILD_TDM1 <- dis_TDM1 / any_ild_TDM1

# No fatal ILD for T-DM1 in DB-05
prop_dis_given_nonfatal_ILD_TDM1 <- prop_dis_ILD_TDM1

# 10-year DFS estimate from CALOR trial (used to derive LRR to DR hazard)
S10_LRR_DFS_base <- 0.70

# DR to BC death hazard
median_OS_DR_months <- 57.1
haz_DR_BCd_base <- log(2) / median_OS_DR_months

l_params_all <- list(
  delta            = delta,
  age_init         = age_init,
  tx_window_months = tx_window_months,
  horizon_years = horizon_years,
  
  hr_RF_DR_TDXd_vs_TDM1 = 0.49,
  
  S10_LRR_DFS = S10_LRR_DFS_base,
  
  median_OS_DR_months = median_OS_DR_months,
  haz_DR_BCd = haz_DR_BCd_base,
  
  haz_OC = haz_OC,
  
  p_any_ILD_TDM1 = p_any_ILD_TDM1_base,
  p_any_ILD_TDXd = p_any_ILD_TDXd_base,
  
  # Overall fatal ILD probability among all T-DXd-assigned patients
  p_fatal_ILD_all_TDXd = p_fatal_ILD_all_TDXd_base,
  
  prop_fatal_ILD_TDM1 = prop_fatal_ILD_TDM1,
  prop_fatal_ILD_TDXd = prop_fatal_ILD_TDXd,
  
  prop_dis_ILD_TDM1 = prop_dis_ILD_TDM1,
  prop_dis_ILD_TDXd = prop_dis_ILD_TDXd,
  
  prop_dis_given_nonfatal_ILD_TDM1 = prop_dis_given_nonfatal_ILD_TDM1,
  prop_dis_given_nonfatal_ILD_TDXd = prop_dis_given_nonfatal_ILD_TDXd,
  
  benefit_frac_postILD = c(
    0.28, 0.49, 0.64, 0.75, 0.83,
    0.89, 0.94, 0.97, 0.99, 1.00
  )
)
saveRDS(l_params_all, "l_params_all_basecase.rds")

# Deterministic scenarios
base_case <- function(params) {
  params
}
weaker_hr <- function(params) {
  params$hr_RF_DR_TDXd_vs_TDM1 <- 0.80
  params
}
higher_ild_incidence <- function(params) {
  params$p_any_ILD_TDXd <- min(1 - 1e-12, params$p_any_ILD_TDXd * 2)
  params
}
combined <- function(params) {
  params <- weaker_hr(params)
  params <- higher_ild_incidence(params)
  params
}
scenarios <- list(
  "Base case" = base_case,
  "HR 0.80" = weaker_hr,
  "Any-grade ILD incidence x2" = higher_ild_incidence,
  "Combined" = combined
)
scenario_order <- names(scenarios)

## Deterministic analyses
# RCB- and ER_status-specific 5-year scaled total recurrence risks
RCB_risk_values_scaled <- tibble::tribble(
  ~RCB, ~ER_status,    ~risk_5y_total_rec_TDM1,
  "RCB-I",    "ER-/HER2+", 0.084,
  "RCB-II",   "ER-/HER2+", 0.221,
  "RCB-III",  "ER-/HER2+", 0.241,
  "RCB-I",    "ER+/HER2+", 0.050,
  "RCB-II",   "ER+/HER2+", 0.138,
  "RCB-III",  "ER+/HER2+", 0.283
) %>%
  mutate(
    # Convert 5-year total recurrence risk to monthly total recurrence hazard
    haz_RF_rec_TDM1 = -log(1 - risk_5y_total_rec_TDM1) / 60,
    
    # Split total recurrence hazard into 90% DR and 10% LRR
    haz_RF_DR_TDM1 = 0.90 * haz_RF_rec_TDM1,
    haz_RF_LRR     = 0.10 * haz_RF_rec_TDM1,
    
    RCB = factor(RCB, levels = c("RCB-I", "RCB-II", "RCB-III")),
    ER_status = factor(ER_status, levels = c("ER-/HER2+", "ER+/HER2+"))
  )

# Function to evaluate each scenario across RCB categories and ER subgroups
eval_os_RCB <- function(params,
                        scenario_name,
                        scenario_modify = function(params) params) {
  
  results <- vector("list", nrow(RCB_risk_values_scaled))
  
  for (i in seq_len(nrow(RCB_risk_values_scaled))) {
    
    risk_row <- RCB_risk_values_scaled[i, ]
    
    # Apply scenario
    current_params <- scenario_modify(params)
    
    # Set baseline (T-DM1) LRR and DR hazards
    current_params$haz_RF_DR_TDM1 <- risk_row$haz_RF_DR_TDM1
    current_params$haz_RF_LRR <- risk_row$haz_RF_LRR
    
    # Run Markov model to obtain 10-year OS
    os_10y <- run_markov_model(params = current_params)
  
    results[[i]] <- tibble::tibble(
      scenario = scenario_name,
      ER_status = risk_row$ER_status,
      RCB = risk_row$RCB,
      
      risk_5y_total_rec_TDM1 = risk_row$risk_5y_total_rec_TDM1,
      haz_RF_rec_TDM1 = risk_row$haz_RF_rec_TDM1,
      haz_RF_DR_TDM1 = current_params$haz_RF_DR_TDM1,
      haz_RF_LRR = current_params$haz_RF_LRR,
      
      OS_10y_TDM1 = os_10y["TDM1"],
      OS_10y_TDXd = os_10y["TDXd"],
      OS_diff_percent_TDXd_minus_TDM1 =
        100 * (os_10y["TDXd"] - os_10y["TDM1"])
    )
  }
  
  dplyr::bind_rows(results)
}

# Run deterministic analysis for each scenario
df_os_RCB_list <- vector("list", length(scenarios))
for (s in seq_along(scenarios)) {
  scenario_name <- names(scenarios)[s]
  cat("Running deterministic scenario:", scenario_name, "\n")
  df_os_RCB_list[[s]] <- eval_os_RCB(
    params = l_params_all,
    scenario_name = scenario_name,
    scenario_modify = scenarios[[s]]
  )
}

# Summarize results of deterministic analysis
df_os_RCB <- dplyr::bind_rows(df_os_RCB_list) %>%
  mutate(scenario = factor(scenario, levels = scenario_order)) %>%
  arrange(ER_status, RCB, scenario)

table_os_RCB <- df_os_RCB %>%
  mutate(risk_5y_total_rec_TDM1_percent = 100 * risk_5y_total_rec_TDM1,
         OS_diff_percent_TDXd_minus_TDM1 = round(OS_diff_percent_TDXd_minus_TDM1, 2)) %>%
  select(scenario, ER_status, RCB, risk_5y_total_rec_TDM1_percent, OS_diff_percent_TDXd_minus_TDM1) %>%
  arrange(ER_status, RCB, scenario)

# Define scenarios for probabilistic analysis
psa_scenarios <- list(
  "Base case PSA" = base_case,
  "HR 0.80 PSA" = weaker_hr,
  "Any-grade ILD incidence x2 PSA" = higher_ild_incidence,
  "Combined PSA" = combined
)
psa_scenario_order <- names(psa_scenarios)

# Probabilistic analysis
set.seed(123)
n_sim <- 1000

# Function to sample probabilities from a logit-normal distribution
draw_prob_logitnormal <- function(mu, lower_95, upper_95, n_sim) {
  mean_logit <- qlogis(mu)
  se_logit <- (qlogis(upper_95) - qlogis(lower_95)) / (2 * 1.96)
  plogis(rnorm(n = n_sim, mean = mean_logit, sd = se_logit))
}

# Obtained using Yau et al EFS estimates and 95% CI (unscaled)
RCB_risk_values_raw <- tibble::tribble(
  ~RCB, ~ER_status,    ~risk_5y, ~lower_95, ~upper_95,
  "RCB-I",    "ER-/HER2+", 0.15,     0.04,      0.24,
  "RCB-II",   "ER-/HER2+", 0.37,     0.25,      0.48,
  "RCB-III",  "ER-/HER2+", 0.40,     0.14,      0.58,
  "RCB-I",    "ER+/HER2+", 0.09,     0.04,      0.15,
  "RCB-II",   "ER+/HER2+", 0.24,     0.18,      0.30,
  "RCB-III",  "ER+/HER2+", 0.46,     0.29,      0.60
) %>%
  mutate(RCB = factor(RCB, levels = c("RCB-I", "RCB-II", "RCB-III")),
         ER_status = factor(ER_status, levels = c("ER-/HER2+", "ER+/HER2+")))

# Any-grade ILD probability over treatment window: T-DXd
p_any_ILD_TDXd_draw <- draw_prob_logitnormal(mu = 0.096, lower_95 = 0.05,
                                                upper_95 = 0.21, n_sim = n_sim)

# Overall fatal ILD probability among all T-DXd-treated patients
p_fatal_ILD_all_TDXd_draw <- draw_prob_logitnormal(mu = 0.01, lower_95 = 0.002,
                                                   upper_95 = 0.030, n_sim = n_sim)

# HR PSA draw: DB-05 DRFI HR and 95% CI for T-DXd vs. T-DM1
hr_mean  <- 0.49
hr_lower <- 0.34
hr_upper <- 0.71
se_log_hr <- (log(hr_upper) - log(hr_lower)) / (2 * 1.96)
hr_RF_DR_TDXd_vs_TDM1_draw <- exp(rnorm(n = n_sim, mean = log(hr_mean), sd = se_log_hr))

# Histograms of PSA draws
hist(p_any_ILD_TDXd_draw)
hist(p_fatal_ILD_all_TDXd_draw)
hist(hr_RF_DR_TDXd_vs_TDM1_draw)

# Function to run probabilistic analysis
eval_os_RCB_uncertainty <- function(params,
                                    scenario_name,
                                    scenario_modify,
                                    n_sim = 1000,
                                    hr_scale = 0.54,
                                    hr_RF_DR_TDXd_vs_TDM1_draw,
                                    p_any_ILD_TDXd_draw,
                                    p_fatal_ILD_all_TDXd_draw) {
  
  results <- vector("list", nrow(RCB_risk_values_raw))
  
  for (i in seq_len(nrow(RCB_risk_values_raw))) {
    
    # Sample Yau et al 5-year recurrence risks
    risk_sim_yau <- draw_prob_logitnormal(
      mu  = RCB_risk_values_raw$risk_5y[i],
      lower_95 = RCB_risk_values_raw$lower_95[i],
      upper_95 = RCB_risk_values_raw$upper_95[i],
      n_sim = n_sim
    )
    
    # Scale each sampled risk using KATHERINE HR
    # Equivalent: risk_sim_tdm1 <- 1 - (1 - risk_sim_yau)^hr_scale
    risk_sim_tdm1 <- 1 - exp(-hr_scale * (-log(1 - risk_sim_yau)))
    
    sim_results <- vector("list", n_sim)
    
    for (s in seq_len(n_sim)) {
      
      current_params <- params
      
      # PSA draw: convert baseline total recurrence risk into monthly recurrence hazard 
      # and split into 90% DR and 10% LRR
      haz_RF_rec_TDM1 <- -log(1 - risk_sim_tdm1[s]) / 60
      current_params$haz_RF_DR_TDM1 <- 0.90 * haz_RF_rec_TDM1
      current_params$haz_RF_LRR <- 0.10 * haz_RF_rec_TDM1
      
      # PSA draw: HR for DR, T-DXd vs T-DM1
      current_params$hr_RF_DR_TDXd_vs_TDM1 <- hr_RF_DR_TDXd_vs_TDM1_draw[s]
      
      # PSA draw: any-grade ILD probability for T-DXd over treatment window
      current_params$p_any_ILD_TDXd <- p_any_ILD_TDXd_draw[s]
      
      # PSA draw: overall fatal ILD probability among all T-DXd-assigned patients
      current_params$p_fatal_ILD_all_TDXd <- p_fatal_ILD_all_TDXd_draw[s]
      
      # Apply PSA scenarios
      current_params <- scenario_modify(current_params)
      
      # Ensure fatal ILD remains a subset of any-grade ILD
      p_any_ILD_TDXd <- current_params$p_any_ILD_TDXd
      p_fatal_ILD_all_TDXd <- current_params$p_fatal_ILD_all_TDXd
      p_fatal_ILD_all_TDXd <- min(p_fatal_ILD_all_TDXd, 0.99 * p_any_ILD_TDXd)
      
      # Nonfatal ILD among all T-DXd-assigned patients
      p_nonfatal_ILD_all_TDXd <- p_any_ILD_TDXd - p_fatal_ILD_all_TDXd
      
      # Fixed probability of discontinuation among nonfatal ILD cases
      prop_dis_given_nonfatal_ILD_TDXd <- current_params$prop_dis_given_nonfatal_ILD_TDXd
      
      # Overall probability of ILD discontinuation among T-DXd-assigned patients
      p_dis_ILD_all_TDXd <- p_nonfatal_ILD_all_TDXd * prop_dis_given_nonfatal_ILD_TDXd
      
      # Convert to proportions among any-grade ILD cases
      prop_fatal_ILD_TDXd <- p_fatal_ILD_all_TDXd / p_any_ILD_TDXd
      prop_dis_ILD_TDXd <- p_dis_ILD_all_TDXd / p_any_ILD_TDXd
      
      # Assign ILD parameters used by the model
      current_params$p_any_ILD_TDXd <- p_any_ILD_TDXd
      current_params$prop_fatal_ILD_TDXd <- prop_fatal_ILD_TDXd
      current_params$prop_dis_ILD_TDXd <- prop_dis_ILD_TDXd
      current_params$prop_dis_given_nonfatal_ILD_TDXd <- prop_dis_given_nonfatal_ILD_TDXd
      
      os_10y <- run_markov_model(params = current_params)
      
      sim_results[[s]] <- tibble::tibble(
        sim = s,
        scenario = scenario_name,
        ER_status = RCB_risk_values_raw$ER_status[i],
        RCB = RCB_risk_values_raw$RCB[i],
        
        risk_5y_yau_sim = risk_sim_yau[s],
        risk_5y_tdm1_sim = risk_sim_tdm1[s],
        
        haz_RF_rec_TDM1 = haz_RF_rec_TDM1,
        haz_RF_DR_TDM1 = current_params$haz_RF_DR_TDM1,
        haz_RF_LRR = current_params$haz_RF_LRR,
        
        HR_RF_DR_TDXd_vs_TDM1 = current_params$hr_RF_DR_TDXd_vs_TDM1,
        
        p_any_ILD_TDXd = p_any_ILD_TDXd,
        p_fatal_ILD_all_TDXd = p_fatal_ILD_all_TDXd,
        p_nonfatal_ILD_all_TDXd = p_nonfatal_ILD_all_TDXd,
        p_dis_ILD_all_TDXd = p_dis_ILD_all_TDXd,
        prop_fatal_ILD_TDXd = prop_fatal_ILD_TDXd,
        prop_dis_ILD_TDXd = prop_dis_ILD_TDXd,
        prop_dis_given_nonfatal_ILD_TDXd = prop_dis_given_nonfatal_ILD_TDXd,
        
        OS_10y_TDM1 = os_10y["TDM1"],
        OS_10y_TDXd = os_10y["TDXd"],
        OS_diff_percent_TDXd_minus_TDM1 =
          100 * (os_10y["TDXd"] - os_10y["TDM1"]),
        favors_TDXd =
          100 * (os_10y["TDXd"] - os_10y["TDM1"]) > 0,
        favors_TDM1 =
          100 * (os_10y["TDXd"] - os_10y["TDM1"]) < 0
      )
    }
    
    results[[i]] <- dplyr::bind_rows(sim_results)
  }
  
  dplyr::bind_rows(results)
}

# Conduct probabilistic analyses
df_os_RCB_psa_list <- vector("list", length(psa_scenarios))
for (s in seq_along(psa_scenarios)) {
  scenario_name <- names(psa_scenarios)[s]
  cat("Running PSA scenario:", scenario_name, "\n")
  df_os_RCB_psa_list[[s]] <- eval_os_RCB_uncertainty(
    params = l_params_all,
    scenario_name = scenario_name,
    scenario_modify = psa_scenarios[[s]],
    n_sim = n_sim,
    hr_scale = 0.54,
    hr_RF_DR_TDXd_vs_TDM1_draw = hr_RF_DR_TDXd_vs_TDM1_draw,
    p_any_ILD_TDXd_draw = p_any_ILD_TDXd_draw,
    p_fatal_ILD_all_TDXd_draw = p_fatal_ILD_all_TDXd_draw
  )
}

# Organize probabilistic analysis results into dataframe
df_os_RCB_psa <- dplyr::bind_rows(df_os_RCB_psa_list) %>%
  mutate(scenario = factor(scenario, levels = psa_scenario_order)) %>%
  arrange(ER_status, RCB, scenario)

# Determine minimum 5-year total recurrence risk in base case PSA among all simulations favoring T-DXd
df_os_RCB_psa %>% 
  filter(OS_diff_percent_TDXd_minus_TDM1 > 0 & scenario == "Base case PSA") %>%
  group_by(RCB, ER_status) %>%
  summarise(min(risk_5y_tdm1_sim)*100)

# Summarize results of probabilistic analysis 
table_os_RCB_psa <- df_os_RCB_psa %>%
  group_by(scenario, ER_status, RCB) %>%
  summarise(
    mean_risk_5y_tdm1 = mean(risk_5y_tdm1_sim),
    mean_haz_RF_rec_TDM1 = mean(haz_RF_rec_TDM1),
    mean_haz_RF_DR_TDM1 = mean(haz_RF_DR_TDM1),
    mean_haz_RF_LRR = mean(haz_RF_LRR),
    mean_OS_diff_percent_TDXd_minus_TDM1 =
      round(mean(OS_diff_percent_TDXd_minus_TDM1), 2),
    median_OS_diff_percent_TDXd_minus_TDM1 =
      round(median(OS_diff_percent_TDXd_minus_TDM1), 2),
    lower_OS_diff_percent_TDXd_minus_TDM1 =
      round(quantile(OS_diff_percent_TDXd_minus_TDM1, 0.025), 2),
    upper_OS_diff_percent_TDXd_minus_TDM1 =
      round(quantile(OS_diff_percent_TDXd_minus_TDM1, 0.975), 2),
    pct_favors_TDXd =
      100 * mean(OS_diff_percent_TDXd_minus_TDM1 > 0),
    pct_favors_TDM1 =
      100 * mean(OS_diff_percent_TDXd_minus_TDM1 < 0),
    .groups = "drop"
  ) %>%
  arrange(scenario, ER_status, RCB)
table_os_RCB_psa
  
# Save results
saveRDS(df_os_RCB_psa, "df_os_RCB_PSA.rds")
saveRDS(table_os_RCB_psa, "table_os_RCB_PSA.rds")

# Visualization of probabilistic analysis results
# Data frame for forest plot
df_forest_psa <- table_os_RCB_psa %>%
  filter(scenario == "Base case PSA") %>%
  mutate(
    y_pos = case_when(
      RCB == "RCB-I"   ~ 3,
      RCB == "RCB-II"  ~ 2,
      RCB == "RCB-III" ~ 1
    ),
    mean_ui_label = sprintf(
      "%.1f (%.1f, %.1f)",
      mean_OS_diff_percent_TDXd_minus_TDM1,
      lower_OS_diff_percent_TDXd_minus_TDM1,
      upper_OS_diff_percent_TDXd_minus_TDM1
    ),
    pct_favors_label = sprintf("%.1f%%", pct_favors_TDXd)
  )

# Ensure rows align exactly across panels
y_limits <- c(0.5, 4.15)
y_header <- 3.78

# Center positions for text columns
x_mean <- 0.75
x_favor <- 2.25

df_header <- tibble::tibble(ER_status = factor(c("ER-/HER2+", "ER+/HER2+"),
                                               levels = c("ER-/HER2+", "ER+/HER2+")),
                            y_pos = y_header)

# Add arrows
arrow_df <- tibble::tibble(ER_status = factor(c("ER+/HER2+", "ER+/HER2+"),
                                              levels = c("ER-/HER2+", "ER+/HER2+")),
  x = c(-0.10,  0.10),
  xend    = c(-2.20, 15.30),
  y       = c(-0.08, -0.08),
  yend    = c(-0.08, -0.08),
  label   = c("Favors\nT-DM1", "Favors T-DXd"),
  label_x = c(-1.25, 7.60),
  label_y = c(-0.28, -0.28)
)

p_forest <- ggplot(
  df_forest_psa,
  aes(
    x = mean_OS_diff_percent_TDXd_minus_TDM1,
    y = y_pos
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.6,
    color = "gray40"
  ) +
  geom_errorbarh(
    aes(
      xmin = lower_OS_diff_percent_TDXd_minus_TDM1,
      xmax = upper_OS_diff_percent_TDXd_minus_TDM1
    ),
    height = 0.12,
    linewidth = 0.8
  ) +
  geom_point(size = 3.2) +
  facet_grid(ER_status ~ ., switch = "y") +
  
  # Directional arrows
  geom_segment(
    data = arrow_df,
    aes(x = x, xend = xend, y = y, yend = yend),
    inherit.aes = FALSE,
    linewidth = 0.8,
    arrow = arrow(length = unit(0.18, "cm"), type = "closed")
  ) +
  geom_text(
    data = arrow_df,
    aes(x = label_x, y = label_y, label = label),
    inherit.aes = FALSE,
    hjust = 0.5,
    vjust = 1,
    size = 4,
    fontface = "bold",
    lineheight = 0.9
  ) +
  
  scale_y_continuous(
    breaks = c(3, 2, 1),
    labels = c("RCB-I", "RCB-II", "RCB-III"),
    expand = c(0, 0)
  ) +
  scale_x_continuous(
    limits = c(-2.5, 15.5),
    breaks = c(0, 5, 10, 15)
  ) +
  coord_cartesian(
    ylim = y_limits,
    clip = "off"
  ) +
  labs(
    x = "Absolute difference in 10-year OS\n(T-DXd minus T-DM1, percentage points)",
    y = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "gray95"),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, face = "bold"),
    legend.position = "none",
    panel.spacing = unit(0.8, "lines"),
    
    # Extra space under the plot
    axis.title.x = element_text(margin = margin(t = 42)),
    plot.margin = margin(5.5, 5.5, 38, 5.5)
  )

p_table <- ggplot(
  df_forest_psa,
  aes(y = y_pos)
) +
  geom_text(
    aes(x = x_mean, label = mean_ui_label),
    hjust = 0.5,
    size = 4
  ) +
  geom_text(
    aes(x = x_favor, label = pct_favors_label),
    hjust = 0.5,
    size = 4
  ) +
  geom_text(
    data = df_header,
    aes(
      x = x_mean,
      y = y_pos,
      label = "Mean (95%\nuncertainty interval)"
    ),
    hjust = 0.5,
    vjust = 0.5,
    fontface = "bold",
    size = 4,
    lineheight = 0.9
  ) +
  geom_text(
    data = df_header,
    aes(
      x = x_favor,
      y = y_pos,
      label = "% of simulations\nfavoring T-DXd"
    ),
    hjust = 0.5,
    vjust = 0.5,
    fontface = "bold",
    size = 4,
    lineheight = 0.9
  ) +
  facet_grid(ER_status ~ .) +
  scale_y_continuous(
    limits = y_limits,
    breaks = NULL,
    expand = c(0, 0)
  ) +
  scale_x_continuous(
    limits = c(0, 3.05),
    expand = c(0, 0)
  ) +
  coord_cartesian(clip = "off") +
  theme_void(base_size = 13) +
  theme(
    strip.text = element_blank(),
    panel.spacing = unit(0.8, "lines"),
    plot.margin = margin(15, 35, 5.5, 10)
  )

final_plot <- patchwork::wrap_plots(
  p_forest,
  p_table,
  ncol = 2,
  widths = c(3.7, 2.6)
)

final_plot

# One-way threshold analysis
threshold_os_RCB <- function(params,
                             scenario_modify = function(params) params) {
  
  results <- vector("list", nrow(RCB_risk_values_scaled))
  
  for (i in seq_len(nrow(RCB_risk_values_scaled))) {
    
    risk_row <- RCB_risk_values_scaled[i, ]
    
    current_params <- params
    
    # Set baseline RCB/subtype-specific recurrence hazards
    current_params$haz_RF_DR_TDM1 <- risk_row$haz_RF_DR_TDM1
    current_params$haz_RF_LRR     <- risk_row$haz_RF_LRR
    
    # Apply threshold/scenario modification AFTER baseline hazards are set
    current_params <- scenario_modify(current_params)
    
    os_10y <- run_markov_model(params = current_params)
    
    results[[i]] <- tibble::tibble(
      ER_status = risk_row$ER_status,
      RCB = risk_row$RCB,
      OS_diff_percent_TDXd_minus_TDM1 =
        100 * (os_10y["TDXd"] - os_10y["TDM1"])
    )
  }
  
  dplyr::bind_rows(results)
}

# Grid search for ILD incidence multiplier threshold
modify_param_ild <- function(params, x) {
  params$p_any_ILD_TDXd <- params$p_any_ILD_TDXd * x
  params
}

max_ild_multiplier <- (1- 1e-12) / l_params_all$p_any_ILD_TDXd
ild_grid <- seq(1, max_ild_multiplier, length.out = 500)
os_diff_ild <- vector("list", length(ild_grid))

for (g in seq_along(ild_grid)) {
  tmp <- threshold_os_RCB(
    params = l_params_all,
    scenario_modify = function(params) {
      modify_param_ild(params, ild_grid[g])
    }
  )
  tmp$threshold <- ild_grid[g]
  os_diff_ild[[g]] <- tmp
}

os_diff_new_ild <- dplyr::bind_rows(os_diff_ild)
first_neg_diff_ild <- os_diff_new_ild %>%
  arrange(ER_status, RCB, threshold) %>%
  group_by(ER_status, RCB) %>%
  filter(OS_diff_percent_TDXd_minus_TDM1 < 0) %>%
  slice(1) %>%
  ungroup()
first_neg_diff_ild

# Grid search for recurrence risk threshold
modify_param_rec <- function(params, p5y_total_rec_TDM1) {
  haz_RF_rec_TDM1 <- -log(1 - p5y_total_rec_TDM1) / 60
  params$haz_RF_DR_TDM1 <- 0.90 * haz_RF_rec_TDM1
  params$haz_RF_LRR     <- 0.10 * haz_RF_rec_TDM1
  params
}

rec_grid <- seq(0, 0.10, length.out = 500)
os_diff_rec <- vector("list", length(rec_grid))

for (g in seq_along(rec_grid)) {
  tmp <- threshold_os_RCB(
    params = l_params_all,
    scenario_modify = function(params) {
      modify_param_rec(params, rec_grid[g])
    }
  )
  tmp$threshold <- rec_grid[g]
  tmp$threshold_percent <- 100 * rec_grid[g]
  os_diff_rec[[g]] <- tmp
}

os_diff_new_rec <- dplyr::bind_rows(os_diff_rec)
rec_threshold <- os_diff_new_rec %>%
  filter(OS_diff_percent_TDXd_minus_TDM1 < 0) %>%
  summarise(
    recurrence_risk_threshold_prob = max(threshold, na.rm = TRUE),
    recurrence_risk_threshold_percent = 100 * recurrence_risk_threshold_prob
  )
rec_threshold
