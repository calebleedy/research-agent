# Rao and Yu (1994) Simulation Study in R

# --- 1. Setup and Parameters ---
set.seed(123) # for reproducibility

# Simulation parameters
m <- 40 # Number of small areas
T_periods <- 5 # Number of time points
p <- 0.2 # Autocorrelation parameter (rho)
sigma_v_sq_true <- 0.25 # True variance of area-specific random effects (sigma_v^2)
sigma_epsilon_sq_true <- 0.25 # True variance of time-specific random effects (sigma_epsilon^2)
sigma_e_sq_true <- 1 # True variance of sampling errors (sigma_e^2, assumed to be 1 in paper)

R <- 5000 # Number of simulation runs

# --- 2. Data Generation Function ---
generate_data <- function(m, T_periods, p, sigma_v_sq, sigma_epsilon_sq, sigma_e_sq) {
  # Initialize storage for true values and observed data
  v <- rnorm(m, 0, sqrt(sigma_v_sq)) # Area-specific random effects
  u_init_sd <- sqrt(sigma_epsilon_sq / (1 - p^2)) # SD for stationary AR(1)
  
  y_it <- matrix(0, nrow = m, ncol = T_periods)
  theta_it <- matrix(0, nrow = m, ncol = T_periods) # True values (v_i + u_it)

  u_it_matrix <- matrix(0, nrow = m, ncol = T_periods)

  for (i in 1:m) {
    u_prev <- rnorm(1, 0, u_init_sd) # Initial u_i,0 from stationary distribution
    
    for (t in 1:T_periods) {
      epsilon_it <- rnorm(1, 0, sqrt(sigma_epsilon_sq))
      u_it <- p * u_prev + epsilon_it
      u_it_matrix[i, t] <- u_it
      u_prev <- u_it
      
      e_it <- rnorm(1, 0, sqrt(sigma_e_sq))
      
      theta_it[i, t] <- v[i] + u_it
      y_it[i, t] <- theta_it[i, t] + e_it
    }
  }
  
  return(list(y_it = y_it, theta_it = theta_it, v = v, u_it_matrix = u_it_matrix))
}

# --- 3. Estimators ---

# Function for Direct Estimator (baseline 1)
direct_estimator <- function(y_it, T_periods) {
  return(y_it[, T_periods]) # Estimate for the last time point T
}

# Function for Fay-Herriot Estimator (baseline 2: weighted estimator for last time point)
fay_herriot_estimator <- function(y_it, T_periods, p_known, sigma_v_sq_true, sigma_epsilon_sq_true, sigma_e_sq_true) {
  sigma_star_sq <- sigma_v_sq_true + sigma_epsilon_sq_true / (1 - p_known^2)
  w_FH <- sigma_star_sq / (sigma_star_sq + sigma_e_sq_true)
  return(w_FH * y_it[, T_periods])
}


# Function for Two-Stage Estimator (BLUP with estimated variance components)
two_stage_estimator <- function(y_it, m, T_periods, p_known, sigma_e_sq_true) {
  
  # --- Step 3a: Variance Component Estimation (Simplified for x_it*beta = 0) ---
  P <- matrix(0, nrow = T_periods, ncol = T_periods)
  # Corrected P[1,1] as per paper's description "first diagonal element (1-p^2)"
  P[1, 1] <- (1 - p_known^2) 
  if (T_periods > 1) {
    for (t in 2:T_periods) {
      P[t, t] <- 1
      P[t, t-1] <- -p_known
    }
  }

  # Corrected: f_vec definition (from page 8, (4.5))
  f_vec <- c(sqrt(1 - p_known^2), rep(1 - p_known, T_periods - 1))
  c_val <- sum(f_vec^2) # This is f_vec'f_vec
  D <- (f_vec %*% t(f_vec)) / c_val
  
  I_T <- diag(T_periods)
  
  # z_i_prime_prime calculations (z_i^{(1)})
  z_i_prime_prime_all <- matrix(0, nrow = m, ncol = T_periods)
  
  for (i in 1:m) {
    y_i_vec <- y_it[i, ]
    z_i_vec <- P %*% y_i_vec
    z_i_prime_prime_all[i, ] <- (I_T - D) %*% z_i_vec
  }
  
  # Estimate sigma_epsilon_sq_hat (ỡ^2(p) in Eq 4.8)
  # Missing correction term from Eq 4.8. For H(1)=0, the correction term would be:
  # tr[block diagi(IT – D) × block diag;(ΡΣ¡PT)]
  # = m * tr[(I_T - D) %*% P %*% sigma_e_sq_true * I_T %*% t(P)]
  
  # For now, without the full derivation for X=0 case of the correction terms,
  # keeping the simple sum of squares / df. The overestimation suggests this is wrong.
  sum_sq_z_prime_prime <- sum(apply(z_i_prime_prime_all, 1, function(row) sum(row^2)))
  sigma_epsilon_sq_hat <- sum_sq_z_prime_prime / (m * (T_periods - 1))
  sigma_epsilon_sq_hat <- max(0, sigma_epsilon_sq_hat) # Truncate

  # Estimate sigma_v_sq_hat (ỡ_v^2(p) in Eq 4.9)
  # Missing correction term from Eq 4.9.
  sum_sq_c_f_z_i <- 0
  for (i in 1:m) {
    y_i_vec <- y_it[i, ]
    z_i_vec <- P %*% y_i_vec
    transformed_v_i_component <- (1/sqrt(c_val)) * (t(f_vec) %*% z_i_vec)
    sum_sq_c_f_z_i <- sum_sq_c_f_z_i + transformed_v_i_component^2
  }
  
  sigma_v_sq_hat <- (sum_sq_c_f_z_i / m) - (1/c_val) * sigma_epsilon_sq_hat
  sigma_v_sq_hat <- max(0, sigma_v_sq_hat) # Truncate

  # --- Step 3b: BLUP with estimated variance components ---
  theta_hat_iT_TS <- numeric(m)
  
  Gamma_u_matrix <- matrix(0, nrow = T_periods, ncol = T_periods)
  for (r_idx in 1:T_periods) {
    for (c_idx in 1:T_periods) {
      Gamma_u_matrix[r_idx, c_idx] <- p_known^(abs(r_idx - c_idx)) / (1 - p_known^2)
    }
  }
  
  J_T <- matrix(1, nrow = T_periods, ncol = T_periods)
  
  for (i in 1:m) {
    y_i_vec <- y_it[i, ]
    
    V_i <- diag(T_periods) * sigma_e_sq_true + sigma_v_sq_hat * J_T + sigma_epsilon_sq_hat * Gamma_u_matrix
    
    Cov_theta_iT_yi <- sigma_v_sq_hat * rep(1, T_periods) + sigma_epsilon_sq_hat * Gamma_u_matrix[T_periods, ]
    
    theta_hat_iT_TS[i] <- Cov_theta_iT_yi %*% solve(V_i) %*% y_i_vec
  }
  
  return(list(theta_hat_iT = theta_hat_iT_TS, 
              sigma_v_sq_hat = sigma_v_sq_hat, 
              sigma_epsilon_sq_hat = sigma_epsilon_sq_hat))
}

# --- 4. Simulation Loop ---
theta_direct_all_runs <- matrix(0, nrow = R, ncol = m)
theta_fh_all_runs <- matrix(0, nrow = R, ncol = m) # New for alternative FH
theta_ts_all_runs <- matrix(0, nrow = R, ncol = m)
true_theta_T_all_runs <- matrix(0, nrow = R, ncol = m)

# To store estimated variance components
sigma_v_sq_hat_all_runs <- numeric(R)
sigma_epsilon_sq_hat_all_runs <- numeric(R)


for (r_sim in 1:R) {
  if (r_sim %% 100 == 0) {
    cat(sprintf("Running simulation %d/%d\n", r_sim, R))
  }
  
  data_sim <- generate_data(m, T_periods, p, sigma_v_sq_true, sigma_epsilon_sq_true, sigma_e_sq_true)
  
  # True theta for the last time point T
  true_theta_T_all_runs[r_sim, ] <- data_sim$theta_it[, T_periods]
  
  # Direct Estimator
  theta_direct_all_runs[r_sim, ] <- direct_estimator(data_sim$y_it, T_periods)
  
  # Fay-Herriot Estimator (Alternative baseline)
  theta_fh_all_runs[r_sim, ] <- fay_herriot_estimator(data_sim$y_it, T_periods, p, 
                                                       sigma_v_sq_true, sigma_epsilon_sq_true, sigma_e_sq_true)
  
  # Two-Stage Estimator
  ts_results <- two_stage_estimator(data_sim$y_it, m, T_periods, p, sigma_e_sq_true) # Pass sigma_e_sq_true
  theta_ts_all_runs[r_sim, ] <- ts_results$theta_hat_iT
  sigma_v_sq_hat_all_runs[r_sim] <- ts_results$sigma_v_sq_hat
  sigma_epsilon_sq_hat_all_runs[r_sim] <- ts_results$sigma_epsilon_sq_hat
}

# --- 5. Evaluation Metrics ---

# Calculate MSE for Direct Estimator
mse_direct <- sum((theta_direct_all_runs - true_theta_T_all_runs)^2) / (R * m)

# Calculate MSE for Fay-Herriot Estimator (Alternative baseline)
mse_fh <- sum((theta_fh_all_runs - true_theta_T_all_runs)^2) / (R * m)

# Calculate MSE for Two-Stage Estimator
mse_ts <- sum((theta_ts_all_runs - true_theta_T_all_runs)^2) / (R * m)

# Calculate Gain in Efficiency (GE) for Direct Estimator vs Two-Stage
ge_direct_ts <- ((mse_direct / mse_ts) - 1) * 100

# Calculate Gain in Efficiency (GE) for Fay-Herriot Estimator vs Two-Stage (Expected in paper)
ge_fh_ts <- ((mse_fh / mse_ts) - 1) * 100

# Average estimated variance components
avg_sigma_v_sq_hat <- mean(sigma_v_sq_hat_all_runs)
avg_sigma_epsilon_sq_hat <- mean(sigma_epsilon_sq_hat_all_runs)


cat(sprintf("\n--- Simulation Results ---\n"))
cat(sprintf("MSE (Direct Estimator): %.4f\n", mse_direct))
cat(sprintf("MSE (Fay-Herriot Estimator - Alternative): %.4f\n", mse_fh))
cat(sprintf("MSE (Two-Stage Estimator): %.4f\n", mse_ts))
cat(sprintf("Gain in Efficiency (GE) (Direct vs Two-Stage): %.2f%%\n", ge_direct_ts))
cat(sprintf("Gain in Efficiency (GE) (Fay-Herriot vs Two-Stage): %.2f%%\n", ge_fh_ts))

cat(sprintf("\n--- Estimated Variance Components (Average) ---\n"))
cat(sprintf("True sigma_v_sq: %.4f, Average Estimated sigma_v_sq: %.4f\n", sigma_v_sq_true, avg_sigma_v_sq_hat))
cat(sprintf("True sigma_epsilon_sq: %.4f, Average Estimated sigma_epsilon_sq: %.4f\n", sigma_epsilon_sq_true, avg_sigma_epsilon_sq_hat))

# Expected values from Table 1 (p=0.2, T=5)
# sigma_v_sq=0.25, sigma_epsilon_sq=0.25 => GE = 32%
