# ==============================================================================
# EMPIRICAL FREQUENCY OF THE MLE UNIQUENESS CONDITION --- Table S2.1 of supplementary
# Reproducibility Script for Table S2.1
# ==============================================================================
# Description: 
# This script contains the Monte Carlo simulation function used to verify the 
# data-dependent uniqueness condition (Equation 2.5 of the manuscript) of the Maximum Likelihood 
# Estimator under the Complete Bipartite (CB) order restriction.
#
# Running this script will output the exact results row-by-row for Table S2.1.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. CORE SIMULATION FUNCTION
# ------------------------------------------------------------------------------
check_uniqueness_rate <- function(n, true_mu, true_sds, M = 10000, seed = 456) {
  
  # Safety check to ensure dimension consistency
  p <- length(n)
  if (length(true_mu) != p || length(true_sds) != p) {
    stop("Input error: 'n', 'true_mu', and 'true_sds' must have the exact same length.")
  }
  
  # Set seed for exact peer-review reproducibility
  set.seed(seed)
  
  # Run the Monte Carlo simulation M times
  simulations <- replicate(M, {
    
    y_bar <- numeric(p)
    gamma <- numeric(p)
    
    for(i in 1:p) {
      # Generate normally distributed data for group i
      y <- rnorm(n[i], mean = true_mu[i], sd = true_sds[i])
      
      # Calculate sample mean and Maximum Likelihood variance (divided by n)
      y_bar[i] <- mean(y)
      gamma[i] <- mean((y - y_bar[i])^2) 
    }
    
    # Test Theorem 1 Condition (Eq 2.5): Gamma_i > max(...)
    rhs <- pmax((y_bar - min(y_bar))^2, (y_bar - max(y_bar))^2)
    
    # Return TRUE if the condition holds universally across all p groups
    return(all(gamma > rhs))
  })
  
  success_count <- sum(simulations)
  
  # Output a clean 1-row data frame
  return(data.frame(
    Groups = p,
    Total_Simulations = format(M, scientific = FALSE),
    Times_Condition_Held = format(success_count, scientific = FALSE),
    Equation_2.5_satisfied = paste0((success_count / M) * 100, "%")
  ))
}

# ------------------------------------------------------------------------------
# 2. DEFINE TRUE MEAN VECTORS
# ------------------------------------------------------------------------------
# Fixed underlying mean vectors used for p=6 and p=8 configurations
mu_p6 <- c(1.0, 0.5, 1.4, 1.6, 1.1, 1.8)
mu_p8 <- c(1.0, 0.5, 0.2, 1.0, 1.4, 1.8, 2.5, 2.0)

# ==============================================================================
# 3. EXECUTE SIMULATIONS ROW-BY-ROW (Matches Table S2.1)
# ==============================================================================
cat("Starting Monte Carlo Simulations (M = 10,000 replications per row)...\n\n")

# ---------------------------------------------------------
# SCENARIOS 1-4: p = 6, Homoscedastic
# ---------------------------------------------------------
cat("Row 1: p=6, Homoscedastic, Balanced Small\n")
print(check_uniqueness_rate(n = c(15, 15, 15, 15, 15, 15), 
                            true_mu = mu_p6, true_sds = c(9, 9, 9, 9, 9, 9)))

cat("\nRow 2: p=6, Homoscedastic, Unbalanced Small\n")
print(check_uniqueness_rate(n = c(12, 10, 15, 18, 16, 8),  
                            true_mu = mu_p6, true_sds = c(9, 9, 9, 9, 9, 9)))

cat("\nRow 3: p=6, Homoscedastic, Balanced Moderate\n")
print(check_uniqueness_rate(n = c(25, 25, 25, 25, 25, 25), 
                            true_mu = mu_p6, true_sds = c(9, 9, 9, 9, 9, 9)))

cat("\nRow 4: p=6, Homoscedastic, Unbalanced Moderate\n")
print(check_uniqueness_rate(n = c(32, 35, 28, 40, 38, 30), 
                            true_mu = mu_p6, true_sds = c(9, 9, 9, 9, 9, 9)))


# ---------------------------------------------------------
# SCENARIOS 5-8: p = 6, Heteroscedastic
# ---------------------------------------------------------
cat("\nRow 5: p=6, Heteroscedastic, Balanced Small\n")
print(check_uniqueness_rate(n = c(15, 15, 15, 15, 15, 15), 
                            true_mu = mu_p6, true_sds = c(12, 15, 10, 8, 14, 9)))

cat("\nRow 6: p=6, Heteroscedastic, Unbalanced Small\n")
print(check_uniqueness_rate(n = c(12, 10, 15, 18, 16, 8),  
                            true_mu = mu_p6, true_sds = c(12, 15, 10, 8, 14, 9)))

cat("\nRow 7: p=6, Heteroscedastic, Balanced Moderate\n")
print(check_uniqueness_rate(n = c(25, 25, 25, 25, 25, 25), 
                            true_mu = mu_p6, true_sds = c(12, 15, 10, 8, 14, 9)))

cat("\nRow 8: p=6, Heteroscedastic, Unbalanced Moderate\n")
print(check_uniqueness_rate(n = c(32, 35, 28, 40, 38, 30), 
                            true_mu = mu_p6, true_sds = c(12, 15, 10, 8, 14, 9)))


# ---------------------------------------------------------
# SCENARIOS 9-12: p = 8, Homoscedastic
# ---------------------------------------------------------
cat("\nRow 9: p=8, Homoscedastic, Balanced Small\n")
print(check_uniqueness_rate(n = c(20, 20, 20, 20, 20, 20, 20, 20), 
                            true_mu = mu_p8, true_sds = c(8, 8, 8, 8, 8, 8, 8, 8)))

cat("\nRow 10: p=8, Homoscedastic, Unbalanced Small\n")
print(check_uniqueness_rate(n = c(15, 12, 14, 18, 16, 17, 15, 19), 
                            true_mu = mu_p8, true_sds = c(8, 8, 8, 8, 8, 8, 8, 8)))

cat("\nRow 11: p=8, Homoscedastic, Balanced Moderate\n")
print(check_uniqueness_rate(n = c(40, 40, 40, 40, 40, 40, 40, 40), 
                            true_mu = mu_p8, true_sds = c(8, 8, 8, 8, 8, 8, 8, 8)))

cat("\nRow 12: p=8, Homoscedastic, Unbalanced Moderate\n")
print(check_uniqueness_rate(n = c(35, 32, 30, 37, 38, 39, 36, 29), 
                            true_mu = mu_p8, true_sds = c(8, 8, 8, 8, 8, 8, 8, 8)))


# ---------------------------------------------------------
# SCENARIOS 13-16: p = 8, Heteroscedastic
# ---------------------------------------------------------
cat("\nRow 13: p=8, Heteroscedastic, Balanced Small\n")
print(check_uniqueness_rate(n = c(20, 20, 20, 20, 20, 20, 20, 20), 
                            true_mu = mu_p8, true_sds = c(10, 12, 11, 9, 15, 9, 16, 14)))

cat("\nRow 14: p=8, Heteroscedastic, Unbalanced Small\n")
print(check_uniqueness_rate(n = c(15, 12, 14, 18, 16, 17, 15, 19), 
                            true_mu = mu_p8, true_sds = c(10, 12, 11, 9, 15, 9, 16, 14)))

cat("\nRow 15: p=8, Heteroscedastic, Balanced Moderate\n")
print(check_uniqueness_rate(n = c(40, 40, 40, 40, 40, 40, 40, 40), 
                            true_mu = mu_p8, true_sds = c(10, 12, 11, 9, 15, 9, 16, 14)))

cat("\nRow 16: p=8, Heteroscedastic, Unbalanced Moderate\n")
print(check_uniqueness_rate(n = c(35, 32, 30, 37, 38, 39, 36, 29), 
                            true_mu = mu_p8, true_sds = c(10, 12, 11, 9, 15, 9, 16, 14)))

cat("\n=== All Simulations Completed Successfully ===\n")