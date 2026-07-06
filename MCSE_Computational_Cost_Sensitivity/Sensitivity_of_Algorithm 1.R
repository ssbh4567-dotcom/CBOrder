# ==============================================================================
# SENSITIVITY TO INITIALIZATION SIMULATION (Algorithm 1) --- Table S3.2
# Reproducibility Script 
# ==============================================================================
# Description: 
# This script evaluates the sensitivity of the proposed algorithm (Algorithm 1) to its
# starting values. For each generated data set, it compares 
# the final estimates of the principled unconstrained MLE initialization 
# against a randomly perturbed initialization.
#
# ==============================================================================

# Check each package one by one. If it's missing, install it.
if (!require("Iso")) install.packages("Iso")
# Load them
library(Iso)
# ------------------------------------------------------------------------------
# 1. CORE ALGORITHMS 
# ------------------------------------------------------------------------------
CB_MLE <- function(X, w, k1) {
  X1 <- X[1:k1]
  X2 <- X[(k1+1):(length(X))]
  w1 <- w[1:k1]
  w2 <- w[(k1+1):(length(X))]
  sort_ind1 <- order(X1)
  sort_ind2 <- order(X2)
  w1_sorted <- w1[sort_ind1]
  w2_sorted <- w2[sort_ind2]
  sortX1 <- sort(X1)
  sortX2 <- sort(X2)
  X_new <- c(sortX1,sortX2)
  w_new <- c(w1_sorted,w2_sorted)
  X_star <- pava(X_new,w_new)
  X_star1 <- X_star[1:k1]
  X_star2 <- X_star[(k1+1):(length(X))]
  CB1 <- X_star1[sort_ind1]
  CB2 <- X_star2[sort_ind2]
  CB_mle <- c(CB1,CB2)
  return(CB_mle)
}
CB_H1_custom_init <- function(sample_data_list, k2, init_mu, init_var) {
  n <- sapply(sample_data_list, length)
  
  mu0 <- init_mu
  var0 <- init_var
  w0 <- n / var0
  
  repeat {
    new_mu0 <- CB_MLE(sapply(sample_data_list, mean), w0, k2)
    new_var0 <- sapply(1:length(sample_data_list), function(i) {
      sum((sample_data_list[[i]] - new_mu0[i])^2) / n[[i]]
    })
    new_w0 <- n / new_var0
    
    if (max(abs(new_mu0 - mu0)) <= 0.0001) {
      break 
    }
    
    w0 <- new_w0
    mu0 <- new_mu0
    var0 <- new_var0
  }
  
  return(list(mu_hat = mu0, var_hat = var0))
}

# ------------------------------------------------------------------------------
# 2. SENSITIVITY SIMULATION
# ------------------------------------------------------------------------------
check_sensitivity_rate <- function(n, true_mu, true_sds, k2, Q = 10000, seed = 456) {
  p <- length(n)
  if (length(true_mu) != p || length(true_sds) != p) {
    stop("Input error: 'n', 'true_mu', and 'true_sds' must have the exact same length.")
  }
  
  set.seed(seed)
  
  simulations <- replicate(Q, {
    # Generate Data
    sample_data <- lapply(1:p, function(i) rnorm(n[i], mean = true_mu[i], sd = true_sds[i]))
    
    # A. Baseline Initialization (Unconstrained MLE)
    base_mu <- sapply(sample_data, mean)
    base_var <- sapply(1:p, function(i) sum((sample_data[[i]] - base_mu[i])^2) / n[i])
    
    # B. Baseline Run
    base_res <- CB_H1_custom_init(sample_data, k2, base_mu, base_var)
    
    # C. Perturbed Initialization (Random normal noise added to means, variance scaled)
    pert_mu <- base_mu + rnorm(p, mean = 0, sd = 3.0)
    pert_var <- base_var * runif(p, min = 0.1, max = 5)
    
    # D. Perturbed Run
    pert_res <- CB_H1_custom_init(sample_data, k2, pert_mu, pert_var)
    
    # E. Compare final answers
    diff_mu <- max(abs(pert_res$mu_hat - base_res$mu_hat))
    diff_var <- max(abs(pert_res$var_hat - base_res$var_hat))
    
    return(diff_mu < 1e-4 && diff_var < 1e-4)
  })
  
  success_count <- sum(simulations)
  
  return(data.frame(
    Groups = p,
    Total_Simulations = format(Q, scientific = FALSE),
    Times_Matched = format(success_count, scientific = FALSE),
    Identical_sol = paste0(round((success_count / Q) * 100, 2), "%")
  ))
}

# ------------------------------------------------------------------------------
# 3. DEFINE TRUE MEAN VECTORS
# ------------------------------------------------------------------------------
mu_p6 <- c(1.0, 0.5, 1.4, 1.6, 1.1, 1.8)
mu_p8 <- c(1.0, 0.5, 0.2, 1.0, 1.4, 1.8, 2.5, 2.0)

# ==============================================================================
# 4. EXECUTE SIMULATIONS ROW-BY-ROW (Q = 10000)
# ==============================================================================
cat("Starting Sensitivity Simulations (Q = 10000 replications per row)...\n\n")

# ---------------------------------------------------------
# SCENARIOS 1-4: p = 6 (k2 = 3), Homoscedastic (Var = 9)
# ---------------------------------------------------------
cat("Row 1: p=6, Homoscedastic, Balanced Small\n")
print(check_sensitivity_rate(n = c(15, 15, 15, 15, 15, 15), k2 = 3, 
                             true_mu = mu_p6, true_sds = sqrt(c(9, 9, 9, 9, 9, 9))))

cat("\nRow 2: p=6, Homoscedastic, Unbalanced Small\n")
print(check_sensitivity_rate(n = c(12, 10, 15, 18, 16, 8), k2 = 3,  
                             true_mu = mu_p6, true_sds = sqrt(c(9, 9, 9, 9, 9, 9))))

cat("\nRow 3: p=6, Homoscedastic, Balanced Moderate\n")
print(check_sensitivity_rate(n = c(25, 25, 25, 25, 25, 25), k2 = 3, 
                             true_mu = mu_p6, true_sds = sqrt(c(9, 9, 9, 9, 9, 9))))

cat("\nRow 4: p=6, Homoscedastic, Unbalanced Moderate\n")
print(check_sensitivity_rate(n = c(32, 35, 28, 40, 38, 30), k2 = 3, 
                             true_mu = mu_p6, true_sds = sqrt(c(9, 9, 9, 9, 9, 9))))


# ---------------------------------------------------------
# SCENARIOS 5-8: p = 6 (k2 = 3), Heteroscedastic 
# ---------------------------------------------------------
cat("\nRow 5: p=6, Heteroscedastic, Balanced Small\n")
print(check_sensitivity_rate(n = c(15, 15, 15, 15, 15, 15), k2 = 3, 
                             true_mu = mu_p6, true_sds = sqrt(c(12, 15, 10, 8, 14, 9))))

cat("\nRow 6: p=6, Heteroscedastic, Unbalanced Small\n")
print(check_sensitivity_rate(n = c(12, 10, 15, 18, 16, 8), k2 = 3,  
                             true_mu = mu_p6, true_sds = sqrt(c(12, 15, 10, 8, 14, 9))))

cat("\nRow 7: p=6, Heteroscedastic, Balanced Moderate\n")
print(check_sensitivity_rate(n = c(25, 25, 25, 25, 25, 25), k2 = 3, 
                             true_mu = mu_p6, true_sds = sqrt(c(12, 15, 10, 8, 14, 9))))

cat("\nRow 8: p=6, Heteroscedastic, Unbalanced Moderate\n")
print(check_sensitivity_rate(n = c(32, 35, 28, 40, 38, 30), k2 = 3, 
                             true_mu = mu_p6, true_sds = sqrt(c(12, 15, 10, 8, 14, 9))))


# ---------------------------------------------------------
# SCENARIOS 9-12: p = 8 (k2 = 4), Homoscedastic (Var = 8)
# ---------------------------------------------------------
cat("\nRow 9: p=8, Homoscedastic, Balanced Small\n")
print(check_sensitivity_rate(n = c(20, 20, 20, 20, 20, 20, 20, 20), k2 = 4, 
                             true_mu = mu_p8, true_sds = sqrt(c(8, 8, 8, 8, 8, 8, 8, 8))))

cat("\nRow 10: p=8, Homoscedastic, Unbalanced Small\n")
print(check_sensitivity_rate(n = c(15, 12, 14, 18, 16, 17, 15, 19), k2 = 4, 
                             true_mu = mu_p8, true_sds = sqrt(c(8, 8, 8, 8, 8, 8, 8, 8))))

cat("\nRow 11: p=8, Homoscedastic, Balanced Moderate\n")
print(check_sensitivity_rate(n = c(40, 40, 40, 40, 40, 40, 40, 40), k2 = 4, 
                             true_mu = mu_p8, true_sds = sqrt(c(8, 8, 8, 8, 8, 8, 8, 8))))

cat("\nRow 12: p=8, Homoscedastic, Unbalanced Moderate\n")
print(check_sensitivity_rate(n = c(35, 32, 30, 37, 38, 39, 36, 29), k2 = 4, 
                             true_mu = mu_p8, true_sds = sqrt(c(8, 8, 8, 8, 8, 8, 8, 8))))


# ---------------------------------------------------------
# SCENARIOS 13-16: p = 8 (k2 = 4), Heteroscedastic
# ---------------------------------------------------------
cat("\nRow 13: p=8, Heteroscedastic, Balanced Small\n")
print(check_sensitivity_rate(n = c(20, 20, 20, 20, 20, 20, 20, 20), k2 = 4, 
                             true_mu = mu_p8, true_sds = sqrt(c(10, 12, 11, 9, 15, 9, 16, 14))))

cat("\nRow 14: p=8, Heteroscedastic, Unbalanced Small\n")
print(check_sensitivity_rate(n = c(15, 12, 14, 18, 16, 17, 15, 19), k2 = 4, 
                             true_mu = mu_p8, true_sds = sqrt(c(10, 12, 11, 9, 15, 9, 16, 14))))

cat("\nRow 15: p=8, Heteroscedastic, Balanced Moderate\n")
print(check_sensitivity_rate(n = c(40, 40, 40, 40, 40, 40, 40, 40), k2 = 4, 
                             true_mu = mu_p8, true_sds = sqrt(c(10, 12, 11, 9, 15, 9, 16, 14))))

cat("\nRow 16: p=8, Heteroscedastic, Unbalanced Moderate\n")
print(check_sensitivity_rate(n = c(35, 32, 30, 37, 38, 39, 36, 29), k2 = 4, 
                             true_mu = mu_p8, true_sds = sqrt(c(10, 12, 11, 9, 15, 9, 16, 14))))

cat("\n=== All Simulations Completed Successfully ===\n")




