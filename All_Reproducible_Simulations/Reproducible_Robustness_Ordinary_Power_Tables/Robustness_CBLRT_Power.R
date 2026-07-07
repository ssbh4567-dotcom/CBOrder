# ==============================================================================================================
# REPLICATION SCRIPT: Robustness Analysis (Table S7.20 of Supplementary Material --- CBLRT Ordinary Power Table)
#
# The functions called in this script (CB_LRT, CB_LRT_mix1, CB_LRT_mix2, CB_LRT_mix3, CB_LRT_mix4, CB_LRT_Lap)
# are designed to test the robustness of the Complete Bipartite Likelihood
# Ratio Test (CBLRT) under deviations from normality.
#
#
# RETURN VALUE:
# Each function returns the estimated probability of Type I Error or power.
#
# DISTRIBUTION MAPPINGS (As defined in Section: Robustness Analysis):
# - CB_LRT      : Normal distribution (N).
# - CB_LRT_mix1 : Mixture of N(mu_i + 4, (sigma_i + 4)^2) and N(mu_i, sigma_i^2)
#                 with mixture weights (0.15, 0.85).
# - CB_LRT_mix2 : Mixture of N(mu_i + 4, (sigma_i + 7)^2) and N(mu_i, sigma_i^2)
#                 with mixture weights (0.15, 0.85).
# - CB_LRT_mix3 : Mixture of N(mu_i + 4, (sigma_i + 4)^2) and N(mu_i, sigma_i^2)
#                 with mixture weights (0.30, 0.70).
# - CB_LRT_mix4 : Mixture of N(mu_i + 4, (sigma_i + 7)^2) and N(mu_i, sigma_i^2)
#                 with mixture weights (0.30, 0.70).
# - CB_LRT_Lap  : Laplace distribution with location parameter -2 and scale
#                 parameter 3.

# NOTE ON COMPUTATIONAL TIME:
# The results reported in the manuscript were obtained on a 64-core
# high-performance computing server using M = 5000 and Q = 5000.
# Under these settings, computing a single empirical size value requires
# approximately 30 minutes for the CBLRT, CB_LRT_mix1, CB_LRT_mix2, CB_LRT_mix3, CB_LRT_mix4, and CB_LRT_Lap.
#
# *** WARNING ***
# Reproducing the manuscript results with M = 5000 and Q = 5000 on
# a standard desktop or laptop (e.g., 4–16 CPU cores) can result in extremely
# long execution times, potentially requiring several days to complete.
#
# To facilitate quick reproducibility, this script automatically detects the
# number of available CPU cores and, by default, performs a rapid test run with
# 100 bootstrap repetitions and 100 simulation repetitions. To exactly
# reproduce the empirical size results reported in the manuscript/supplementary, set both
# M and Q to 5000.

# ==============================================================================

# Check each package one by one. If it's missing, install it.
if (!require("Iso")) install.packages("Iso")
if (!require("parallel")) install.packages("parallel")
if (!require("doParallel")) install.packages("doParallel")
if (!require("foreach")) install.packages("foreach")
if (!require("LaplacesDemon")) install.packages("LaplacesDemon")


# Load them
library(Iso)
library(parallel)
library(doParallel)
library(foreach)
library(LaplacesDemon)

# Number of bootstrap repetitions
M <- 100  # Set to 5000 to exactly reproduce the empirical size results reported in the supplementary

# Number of simulation repetitions
Q <- 100   # Set to 5000 to exactly reproduce the empirical size results reported in the supplementary

# Number of CPU cores to use
CORE <- max(1, parallel::detectCores() - 1)  # The manuscript used CORE = 64.




#Main Functions

#MLE Under H1 for known variance
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
#MLE under H0 under heteroscedasticity
CB_H0_new <- function(sample_data_list) {
  means <- sapply(sample_data_list, mean)
  sample_sizes <- sapply(sample_data_list, length)
  S <- unlist(sample_data_list)
  mu1 <- mean(S)
  var1 <- sapply(1:length(sample_data_list), function(i) sum((sample_data_list[[i]] - means[i])^2) / sample_sizes[i])
  u1 <- sample_sizes / var1
  
  repeat {
    new_mu1 <- (sum(u1 * means)) / sum(u1)
    new_var1 <- sapply(1:length(sample_data_list), function(i) sum((sample_data_list[[i]] - new_mu1)^2) / sample_sizes[i])
    new_u1 <- sample_sizes / new_var1
    
    if (max(abs(new_mu1 - mu1)) <= 0.0001) {
      break  # Exit the loop if the difference is less than epsilon
    }
    
    u1 <- new_u1
    mu1 <- new_mu1
    var1 <- new_var1
  }
  
  return(var1)
}
#MLE under H0 under heteroscedasticity
CB_H1_new <- function(sample_data_list, k2) {
  n <- sapply(sample_data_list, length)
  mu0 <- sapply(sample_data_list, mean)
  var0 <- sapply(1:length(sample_data_list), function(i) sum((sample_data_list[[i]] - mu0[i])^2) / n[[i]])
  w0 <- n / var0
  repeat {
    new_mu0 <- CB_MLE(sapply(sample_data_list, mean), w0, k2)
    new_var0 <- sapply(1:length(sample_data_list), function(i) sum((sample_data_list[[i]] - new_mu0[i])^2) / n[[i]])
    new_w0 <- n/new_var0
    
    if (max(abs(new_mu0 - mu0)) <= 0.0001) {
      break  # Exit the loop if the difference is less than epsilon
    }
    
    w0 <- new_w0
    mu0 <- new_mu0
    var0 <- new_var0
  }
  
  return(var0)
}

#Normal
CB_LRT <- function(num_datasets, n, mean, sd, k0, num_repetitions = M, outnum_repetitions = Q, num_cores = CORE) {
  cl <- makeCluster(num_cores)
  on.exit({
    stopCluster(cl)
    registerDoSEQ()
    closeAllConnections()
  }, add = TRUE)
  registerDoParallel(cl)
  clusterExport(cl, c("pava", "CB_H0_new", "CB_H1_new", "CB_MLE"))
  results <- numeric(1)
  for (a in 1:1) {
    set.seed(123 + a)
    results[a] <- mean(foreach(rep1 = 1:outnum_repetitions, .combine = "c", .packages = c("parallel", "foreach")) %dopar% {
      set.seed(123 + a*1000 + rep1)
      sample_data_list1 <- lapply(1:num_datasets, function(i) rnorm(n[i], mean[i], sd[i]))
      lambda_values_star <- numeric(num_repetitions)
      for (rep in 1:num_repetitions) {
        sample_data_list2 <- lapply(1:num_datasets, function(i) rnorm(n[i], 0, sqrt(var(sample_data_list1[[i]]))))
        V_R_star <- CB_H1_new(sample_data_list2, k0) / CB_H0_new(sample_data_list2)
        weights <- sapply(1:num_datasets, function(i) V_R_star[i]^(length(sample_data_list2[[i]]) / 2))
        lambda_values_star[rep] <- prod(weights)
      }
      sort_lambda_star <- sort(lambda_values_star)
      quantile_value <- quantile(sort_lambda_star, probs = 0.05)
      V_R <- CB_H1_new(sample_data_list1, k0) / CB_H0_new(sample_data_list1)
      weights <- sapply(1:num_datasets, function(i) V_R[i]^(length(sample_data_list1[[i]]) / 2))
      lambda_value <- prod(weights)
      alpha_LRT <- ifelse(quantile_value > lambda_value, 1, 0)
      return(alpha_LRT)
    })
  }
  return(results)
}

#Mix-normal with mixture weights = c(0.15,0.85),mu = c((mean[i]+4),mean[i]),sigma = c((sd[i]+4),sd[i])
CB_LRT_mix1 <- function(num_datasets, n, mean, sd, k0, num_repetitions = M, outnum_repetitions = Q, num_cores = CORE) {
  cl <- makeCluster(num_cores)
  on.exit({
    stopCluster(cl)
    registerDoSEQ()
    closeAllConnections()
  }, add = TRUE)
  registerDoParallel(cl)
  clusterExport(cl, c("pava", "CB_H0_new", "CB_H1_new", "CB_MLE"))
  results <- numeric(1)
  for (a in 1:1) {
    set.seed(123 + a)
    results[a] <- mean(foreach(rep1 = 1:outnum_repetitions, .combine = "c", .packages = c("parallel", "foreach")) %dopar% {
      set.seed(123 + a*1000 + rep1)
      library(LaplacesDemon)
      sample <- lapply(1:num_datasets, function(i) ((rnormm(n = n[i], p = c(0.15,0.85),mu = c((mean[i]+4),mean[i]),sigma = c((sd[i]+4),sd[i]))
                                                     -(0.15*(mean[i]+4)+0.85*mean[i]))/sqrt((((mean[i]+4)^2 + (sd[i]+4)^2)*0.15) +
                                                                                              ((mean[i]^2 + sd[i]^2)*0.85) - (0.15*(mean[i]+4)+0.85*mean[i])^2)))
      sample_data_list1 <- lapply(1:num_datasets, function(i) (sample[[i]]*sd[i])+mean[i])
      lambda_values_star <- numeric(num_repetitions)
      for (rep in 1:num_repetitions) {
        sample_data_list2 <- lapply(1:num_datasets, function(i) rnorm(n[i], 0, sqrt(var(sample_data_list1[[i]]))))
        V_R_star <- CB_H1_new(sample_data_list2, k0) / CB_H0_new(sample_data_list2)
        weights <- sapply(1:num_datasets, function(i) V_R_star[i]^(length(sample_data_list2[[i]]) / 2))
        lambda_values_star[rep] <- prod(weights)
      }
      sort_lambda_star <- sort(lambda_values_star)
      quantile_value <- quantile(sort_lambda_star, probs = 0.05)
      V_R <- CB_H1_new(sample_data_list1, k0) / CB_H0_new(sample_data_list1)
      weights <- sapply(1:num_datasets, function(i) V_R[i]^(length(sample_data_list1[[i]]) / 2))
      lambda_value <- prod(weights)
      alpha_LRT <- ifelse(quantile_value > lambda_value, 1, 0)
      return(alpha_LRT)
    })
  }
  return(results)
}

#Mix-normal with mixture weights = c(0.15,0.85),mu = c((mean[i]+4),mean[i]),sigma = c((sd[i]+7),sd[i])
CB_LRT_mix2 <- function(num_datasets, n, mean, sd, k0, num_repetitions = M, outnum_repetitions = Q, num_cores = CORE) {
  cl <- makeCluster(num_cores)
  on.exit({
    stopCluster(cl)
    registerDoSEQ()
    closeAllConnections()
  }, add = TRUE)
  registerDoParallel(cl)
  clusterExport(cl, c("pava", "CB_H0_new", "CB_H1_new", "CB_MLE"))
  results <- numeric(1)
  for (a in 1:1) {
    set.seed(123 + a)
    results[a] <- mean(foreach(rep1 = 1:outnum_repetitions, .combine = "c", .packages = c("parallel", "foreach")) %dopar% {
      set.seed(123 + a*1000 + rep1)
      library(LaplacesDemon)
      sample <- lapply(1:num_datasets, function(i) ((rnormm(n = n[i], p = c(0.15,0.85),mu = c((mean[i]+4),mean[i]),sigma = c((sd[i]+7),sd[i]))
                                                     -(0.15*(mean[i]+4)+0.85*mean[i]))/sqrt((((mean[i]+4)^2 + (sd[i]+7)^2)*0.15) +
                                                                                              ((mean[i]^2 + sd[i]^2)*0.85) - (0.15*(mean[i]+4)+0.85*mean[i])^2)))
      sample_data_list1 <- lapply(1:num_datasets, function(i) (sample[[i]]*sd[i])+mean[i])
      lambda_values_star <- numeric(num_repetitions)
      for (rep in 1:num_repetitions) {
        sample_data_list2 <- lapply(1:num_datasets, function(i) rnorm(n[i], 0, sqrt(var(sample_data_list1[[i]]))))
        V_R_star <- CB_H1_new(sample_data_list2, k0) / CB_H0_new(sample_data_list2)
        weights <- sapply(1:num_datasets, function(i) V_R_star[i]^(length(sample_data_list2[[i]]) / 2))
        lambda_values_star[rep] <- prod(weights)
      }
      sort_lambda_star <- sort(lambda_values_star)
      quantile_value <- quantile(sort_lambda_star, probs = 0.05)
      V_R <- CB_H1_new(sample_data_list1, k0) / CB_H0_new(sample_data_list1)
      weights <- sapply(1:num_datasets, function(i) V_R[i]^(length(sample_data_list1[[i]]) / 2))
      lambda_value <- prod(weights)
      alpha_LRT <- ifelse(quantile_value > lambda_value, 1, 0)
      return(alpha_LRT)
    })
  }
  return(results)
}

#Mix-normal with mixture weights = c(0.30,0.70),mu = c((mean[i]+4),mean[i]),sigma = c((sd[i]+4),sd[i])
CB_LRT_mix3 <- function(num_datasets, n, mean, sd, k0, num_repetitions = M, outnum_repetitions = Q, num_cores = CORE) {
  cl <- makeCluster(num_cores)
  on.exit({
    stopCluster(cl)
    registerDoSEQ()
    closeAllConnections()
  }, add = TRUE)
  registerDoParallel(cl)
  clusterExport(cl, c("pava", "CB_H0_new", "CB_H1_new", "CB_MLE"))
  results <- numeric(1)
  for (a in 1:1) {
    set.seed(123 + a)
    results[a] <- mean(foreach(rep1 = 1:outnum_repetitions, .combine = "c", .packages = c("parallel", "foreach")) %dopar% {
      set.seed(123 + a*1000 + rep1)
      library(LaplacesDemon)
      sample <- lapply(1:num_datasets, function(i) ((rnormm(n = n[i], p = c(0.30,0.70),mu = c((mean[i]+4),mean[i]),sigma = c((sd[i]+4),sd[i]))
                                                     -(0.30*(mean[i]+4)+0.70*mean[i]))/sqrt((((mean[i]+4)^2 + (sd[i]+4)^2)*0.30) +
                                                                                              ((mean[i]^2 + sd[i]^2)*0.70) - (0.30*(mean[i]+4)+0.70*mean[i])^2)))
      sample_data_list1 <- lapply(1:num_datasets, function(i) (sample[[i]]*sd[i])+mean[i])
      lambda_values_star <- numeric(num_repetitions)
      for (rep in 1:num_repetitions) {
        sample_data_list2 <- lapply(1:num_datasets, function(i) rnorm(n[i], 0, sqrt(var(sample_data_list1[[i]]))))
        V_R_star <- CB_H1_new(sample_data_list2, k0) / CB_H0_new(sample_data_list2)
        weights <- sapply(1:num_datasets, function(i) V_R_star[i]^(length(sample_data_list2[[i]]) / 2))
        lambda_values_star[rep] <- prod(weights)
      }
      sort_lambda_star <- sort(lambda_values_star)
      quantile_value <- quantile(sort_lambda_star, probs = 0.05)
      V_R <- CB_H1_new(sample_data_list1, k0) / CB_H0_new(sample_data_list1)
      weights <- sapply(1:num_datasets, function(i) V_R[i]^(length(sample_data_list1[[i]]) / 2))
      lambda_value <- prod(weights)
      alpha_LRT <- ifelse(quantile_value > lambda_value, 1, 0)
      return(alpha_LRT)
    })
  }
  return(results)
}

#Mix-normal with mixture weights = c(0.30,0.70),mu = c((mean[i]+4),mean[i]),sigma = c((sd[i]+7),sd[i])
CB_LRT_mix4 <- function(num_datasets, n, mean, sd, k0, num_repetitions = M, outnum_repetitions = Q, num_cores = CORE) {
  cl <- makeCluster(num_cores)
  on.exit({
    stopCluster(cl)
    registerDoSEQ()
    closeAllConnections()
  }, add = TRUE)
  registerDoParallel(cl)
  clusterExport(cl, c("pava", "CB_H0_new", "CB_H1_new", "CB_MLE"))
  results <- numeric(1)
  for (a in 1:1) {
    set.seed(123 + a)
    results[a] <- mean(foreach(rep1 = 1:outnum_repetitions, .combine = "c", .packages = c("parallel", "foreach")) %dopar% {
      set.seed(123 + a*1000 + rep1)
      library(LaplacesDemon)
      sample <- lapply(1:num_datasets, function(i) ((rnormm(n = n[i], p = c(0.30,0.70),mu = c((mean[i]+4),mean[i]),sigma = c((sd[i]+7),sd[i]))
                                                     -(0.30*(mean[i]+4)+0.70*mean[i]))/sqrt((((mean[i]+4)^2 + (sd[i]+7)^2)*0.30) +
                                                                                              ((mean[i]^2 + sd[i]^2)*0.70) - (0.30*(mean[i]+4)+0.70*mean[i])^2)))
      sample_data_list1 <- lapply(1:num_datasets, function(i) (sample[[i]]*sd[i])+mean[i])
      lambda_values_star <- numeric(num_repetitions)
      for (rep in 1:num_repetitions) {
        sample_data_list2 <- lapply(1:num_datasets, function(i) rnorm(n[i], 0, sqrt(var(sample_data_list1[[i]]))))
        V_R_star <- CB_H1_new(sample_data_list2, k0) / CB_H0_new(sample_data_list2)
        weights <- sapply(1:num_datasets, function(i) V_R_star[i]^(length(sample_data_list2[[i]]) / 2))
        lambda_values_star[rep] <- prod(weights)
      }
      sort_lambda_star <- sort(lambda_values_star)
      quantile_value <- quantile(sort_lambda_star, probs = 0.05)
      V_R <- CB_H1_new(sample_data_list1, k0) / CB_H0_new(sample_data_list1)
      weights <- sapply(1:num_datasets, function(i) V_R[i]^(length(sample_data_list1[[i]]) / 2))
      lambda_value <- prod(weights)
      alpha_LRT <- ifelse(quantile_value > lambda_value, 1, 0)
      return(alpha_LRT)
    })
  }
  return(results)
}

#Laplace with location = -2, scale = 3
CB_LRT_Lap <- function(num_datasets, n, mean, sd, k0, num_repetitions = M, outnum_repetitions = Q, num_cores = CORE) {
  cl <- makeCluster(num_cores)
  on.exit({
    stopCluster(cl)
    registerDoSEQ()
    closeAllConnections()
  }, add = TRUE)
  registerDoParallel(cl)
  clusterExport(cl, c("pava", "CB_H0_new", "CB_H1_new", "CB_MLE"))
  results <- numeric(1)
  for (a in 1:1) {
    set.seed(123 + a)
    results[a] <- mean(foreach(rep1 = 1:outnum_repetitions, .combine = "c", .packages = c("parallel", "foreach")) %dopar% {
      set.seed(123 + a*1000 + rep1)
      library(LaplacesDemon)
      sample <- lapply(1:num_datasets, function(i) ((rlaplace(n = n[i], location = -2, scale = 3)+2)/sqrt(18)))
      sample_data_list1 <- lapply(1:num_datasets, function(i) (sample[[i]]*sd[i])+mean[i])
      lambda_values_star <- numeric(num_repetitions)
      for (rep in 1:num_repetitions) {
        sample_data_list2 <- lapply(1:num_datasets, function(i) rnorm(n[i], 0, sqrt(var(sample_data_list1[[i]]))))
        V_R_star <- CB_H1_new(sample_data_list2, k0) / CB_H0_new(sample_data_list2)
        weights <- sapply(1:num_datasets, function(i) V_R_star[i]^(length(sample_data_list2[[i]]) / 2))
        lambda_values_star[rep] <- prod(weights)
      }
      sort_lambda_star <- sort(lambda_values_star)
      quantile_value <- quantile(sort_lambda_star, probs = 0.05)
      V_R <- CB_H1_new(sample_data_list1, k0) / CB_H0_new(sample_data_list1)
      weights <- sapply(1:num_datasets, function(i) V_R[i]^(length(sample_data_list1[[i]]) / 2))
      lambda_value <- prod(weights)
      alpha_LRT <- ifelse(quantile_value > lambda_value, 1, 0)
      return(alpha_LRT)
    })
  }
  return(results)
}



cat("\n==========================================================================\n")
cat("STARTING ROBUSTNESS ANALYSIS CBLRT: Table S7.20 of Supplementary Material (Ordinary Power Tabel)\n")
cat("\n==========================================================================\n")


cat("\n===================================================================================================\n")
cat("Parameters: p=6, k=2, sample sizes = (12, 10, 15, 18, 16, 8), mean = c*(1, 0.5, 1.4, 1.6, 1.1, 1.8), 
    variance = (12, 15, 10, 8, 14, 9), alpha=0.05\n") 
cat("=====================================================================================================\n")

num_datasets = 6
n = c(12,10,15,18,16,8)
mean = c(1,0.5,1.4,1.6,1.1,1.8)
sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9))


# ==============================================================================
cat("\n--- Distribution: Normal (N) ---\n")
# ==============================================================================
res <- CB_LRT(num_datasets, n, 0*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 0    | Size/Power: ", res, "\n")
res <- CB_LRT(num_datasets, n, 1*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 1    | Size/Power: ", res, "\n")
res <- CB_LRT(num_datasets, n, 1.25*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 1.25 | Size/Power: ", res, "\n")
res <- CB_LRT(num_datasets, n, 1.5*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 1.5  | Size/Power: ", res, "\n")
res <- CB_LRT(num_datasets, n, 1.75*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 1.75 | Size/Power: ", res, "\n")
res <- CB_LRT(num_datasets, n, 2*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 2    | Size/Power: ", res, "\n")
res <- CB_LRT(num_datasets, n, 2.25*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 2.25 | Size/Power: ", res, "\n")
res <- CB_LRT(num_datasets, n, 2.5*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 2.5  | Size/Power: ", res, "\n")
res <- CB_LRT(num_datasets, n, 2.75*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 2.75 | Size/Power: ", res, "\n")
res <- CB_LRT(num_datasets, n, 3*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 3    | Size/Power: ", res, "\n")
res <- CB_LRT(num_datasets, n, 3.25*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 3.25 | Size/Power: ", res, "\n")
res <- CB_LRT(num_datasets, n, 3.5*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 3.5  | Size/Power: ", res, "\n")
res <- CB_LRT(num_datasets, n, 3.75*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 3.75 | Size/Power: ", res, "\n")
res <- CB_LRT(num_datasets, n, 4*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 4    | Size/Power: ", res, "\n")
res <- CB_LRT(num_datasets, n, 4.25*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 4.25 | Size/Power: ", res, "\n")
res <- CB_LRT(num_datasets, n, 4.5*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 4.5  | Size/Power: ", res, "\n")
res <- CB_LRT(num_datasets, n, 4.75*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 4.75 | Size/Power: ", res, "\n")
res <- CB_LRT(num_datasets, n, 5*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 5    | Size/Power: ", res, "\n")


# ==============================================================================
cat("\n--- Distribution: Mix-1 ---\n")
# ==============================================================================
res <- CB_LRT_mix1(num_datasets, n, 0*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 0    | Size/Power: ", res, "\n")
res <- CB_LRT_mix1(num_datasets, n, 1*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 1    | Size/Power: ", res, "\n")
res <- CB_LRT_mix1(num_datasets, n, 1.25*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 1.25 | Size/Power: ", res, "\n")
res <- CB_LRT_mix1(num_datasets, n, 1.5*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 1.5  | Size/Power: ", res, "\n")
res <- CB_LRT_mix1(num_datasets, n, 1.75*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 1.75 | Size/Power: ", res, "\n")
res <- CB_LRT_mix1(num_datasets, n, 2*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 2    | Size/Power: ", res, "\n")
res <- CB_LRT_mix1(num_datasets, n, 2.25*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 2.25 | Size/Power: ", res, "\n")
res <- CB_LRT_mix1(num_datasets, n, 2.5*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 2.5  | Size/Power: ", res, "\n")
res <- CB_LRT_mix1(num_datasets, n, 2.75*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 2.75 | Size/Power: ", res, "\n")
res <- CB_LRT_mix1(num_datasets, n, 3*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 3    | Size/Power: ", res, "\n")
res <- CB_LRT_mix1(num_datasets, n, 3.25*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 3.25 | Size/Power: ", res, "\n")
res <- CB_LRT_mix1(num_datasets, n, 3.5*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 3.5  | Size/Power: ", res, "\n")
res <- CB_LRT_mix1(num_datasets, n, 3.75*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 3.75 | Size/Power: ", res, "\n")
res <- CB_LRT_mix1(num_datasets, n, 4*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 4    | Size/Power: ", res, "\n")
res <- CB_LRT_mix1(num_datasets, n, 4.25*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 4.25 | Size/Power: ", res, "\n")
res <- CB_LRT_mix1(num_datasets, n, 4.5*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 4.5  | Size/Power: ", res, "\n")
res <- CB_LRT_mix1(num_datasets, n, 4.75*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 4.75 | Size/Power: ", res, "\n")
res <- CB_LRT_mix1(num_datasets, n, 5*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 5    | Size/Power: ", res, "\n")


# ==============================================================================
cat("\n--- Distribution: Mix-2 ---\n")
# ==============================================================================
res <- CB_LRT_mix2(num_datasets, n, 0*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 0    | Size/Power: ", res, "\n")
res <- CB_LRT_mix2(num_datasets, n, 1*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 1    | Size/Power: ", res, "\n")
res <- CB_LRT_mix2(num_datasets, n, 1.25*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 1.25 | Size/Power: ", res, "\n")
res <- CB_LRT_mix2(num_datasets, n, 1.5*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 1.5  | Size/Power: ", res, "\n")
res <- CB_LRT_mix2(num_datasets, n, 1.75*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 1.75 | Size/Power: ", res, "\n")
res <- CB_LRT_mix2(num_datasets, n, 2*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 2    | Size/Power: ", res, "\n")
res <- CB_LRT_mix2(num_datasets, n, 2.25*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 2.25 | Size/Power: ", res, "\n")
res <- CB_LRT_mix2(num_datasets, n, 2.5*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 2.5  | Size/Power: ", res, "\n")
res <- CB_LRT_mix2(num_datasets, n, 2.75*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 2.75 | Size/Power: ", res, "\n")
res <- CB_LRT_mix2(num_datasets, n, 3*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 3    | Size/Power: ", res, "\n")
res <- CB_LRT_mix2(num_datasets, n, 3.25*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 3.25 | Size/Power: ", res, "\n")
res <- CB_LRT_mix2(num_datasets, n, 3.5*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 3.5  | Size/Power: ", res, "\n")
res <- CB_LRT_mix2(num_datasets, n, 3.75*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 3.75 | Size/Power: ", res, "\n")
res <- CB_LRT_mix2(num_datasets, n, 4*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 4    | Size/Power: ", res, "\n")
res <- CB_LRT_mix2(num_datasets, n, 4.25*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 4.25 | Size/Power: ", res, "\n")
res <- CB_LRT_mix2(num_datasets, n, 4.5*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 4.5  | Size/Power: ", res, "\n")
res <- CB_LRT_mix2(num_datasets, n, 4.75*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 4.75 | Size/Power: ", res, "\n")
res <- CB_LRT_mix2(num_datasets, n, 5*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 5    | Size/Power: ", res, "\n")


# ==============================================================================
cat("\n--- Distribution: Mix-3 ---\n")
# ==============================================================================
res <- CB_LRT_mix3(num_datasets, n, 0*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 0    | Size/Power: ", res, "\n")
res <- CB_LRT_mix3(num_datasets, n, 1*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 1    | Size/Power: ", res, "\n")
res <- CB_LRT_mix3(num_datasets, n, 1.25*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 1.25 | Size/Power: ", res, "\n")
res <- CB_LRT_mix3(num_datasets, n, 1.5*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 1.5  | Size/Power: ", res, "\n")
res <- CB_LRT_mix3(num_datasets, n, 1.75*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 1.75 | Size/Power: ", res, "\n")
res <- CB_LRT_mix3(num_datasets, n, 2*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 2    | Size/Power: ", res, "\n")
res <- CB_LRT_mix3(num_datasets, n, 2.25*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 2.25 | Size/Power: ", res, "\n")
res <- CB_LRT_mix3(num_datasets, n, 2.5*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 2.5  | Size/Power: ", res, "\n")
res <- CB_LRT_mix3(num_datasets, n, 2.75*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 2.75 | Size/Power: ", res, "\n")
res <- CB_LRT_mix3(num_datasets, n, 3*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 3    | Size/Power: ", res, "\n")
res <- CB_LRT_mix3(num_datasets, n, 3.25*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 3.25 | Size/Power: ", res, "\n")
res <- CB_LRT_mix3(num_datasets, n, 3.5*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 3.5  | Size/Power: ", res, "\n")
res <- CB_LRT_mix3(num_datasets, n, 3.75*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 3.75 | Size/Power: ", res, "\n")
res <- CB_LRT_mix3(num_datasets, n, 4*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 4    | Size/Power: ", res, "\n")
res <- CB_LRT_mix3(num_datasets, n, 4.25*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 4.25 | Size/Power: ", res, "\n")
res <- CB_LRT_mix3(num_datasets, n, 4.5*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 4.5  | Size/Power: ", res, "\n")
res <- CB_LRT_mix3(num_datasets, n, 4.75*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 4.75 | Size/Power: ", res, "\n")
res <- CB_LRT_mix3(num_datasets, n, 5*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 5    | Size/Power: ", res, "\n")


# ==============================================================================
cat("\n--- Distribution: Mix-4 ---\n")
# ==============================================================================
res <- CB_LRT_mix4(num_datasets, n, 0*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 0    | Size/Power: ", res, "\n")
res <- CB_LRT_mix4(num_datasets, n, 1*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 1    | Size/Power: ", res, "\n")
res <- CB_LRT_mix4(num_datasets, n, 1.25*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 1.25 | Size/Power: ", res, "\n")
res <- CB_LRT_mix4(num_datasets, n, 1.5*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 1.5  | Size/Power: ", res, "\n")
res <- CB_LRT_mix4(num_datasets, n, 1.75*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 1.75 | Size/Power: ", res, "\n")
res <- CB_LRT_mix4(num_datasets, n, 2*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 2    | Size/Power: ", res, "\n")
res <- CB_LRT_mix4(num_datasets, n, 2.25*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 2.25 | Size/Power: ", res, "\n")
res <- CB_LRT_mix4(num_datasets, n, 2.5*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 2.5  | Size/Power: ", res, "\n")
res <- CB_LRT_mix4(num_datasets, n, 2.75*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 2.75 | Size/Power: ", res, "\n")
res <- CB_LRT_mix4(num_datasets, n, 3*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 3    | Size/Power: ", res, "\n")
res <- CB_LRT_mix4(num_datasets, n, 3.25*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 3.25 | Size/Power: ", res, "\n")
res <- CB_LRT_mix4(num_datasets, n, 3.5*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 3.5  | Size/Power: ", res, "\n")
res <- CB_LRT_mix4(num_datasets, n, 3.75*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 3.75 | Size/Power: ", res, "\n")
res <- CB_LRT_mix4(num_datasets, n, 4*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 4    | Size/Power: ", res, "\n")
res <- CB_LRT_mix4(num_datasets, n, 4.25*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 4.25 | Size/Power: ", res, "\n")
res <- CB_LRT_mix4(num_datasets, n, 4.5*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 4.5  | Size/Power: ", res, "\n")
res <- CB_LRT_mix4(num_datasets, n, 4.75*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 4.75 | Size/Power: ", res, "\n")
res <- CB_LRT_mix4(num_datasets, n, 5*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 5    | Size/Power: ", res, "\n")


# ==============================================================================
cat("\n--- Distribution: Laplace ---\n")
# ==============================================================================
res <- CB_LRT_Lap(num_datasets, n, 0*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 0    | Size/Power: ", res, "\n")
res <- CB_LRT_Lap(num_datasets, n, 1*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 1    | Size/Power: ", res, "\n")
res <- CB_LRT_Lap(num_datasets, n, 1.25*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 1.25 | Size/Power: ", res, "\n")
res <- CB_LRT_Lap(num_datasets, n, 1.5*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 1.5  | Size/Power: ", res, "\n")
res <- CB_LRT_Lap(num_datasets, n, 1.75*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 1.75 | Size/Power: ", res, "\n")
res <- CB_LRT_Lap(num_datasets, n, 2*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 2    | Size/Power: ", res, "\n")
res <- CB_LRT_Lap(num_datasets, n, 2.25*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 2.25 | Size/Power: ", res, "\n")
res <- CB_LRT_Lap(num_datasets, n, 2.5*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 2.5  | Size/Power: ", res, "\n")
res <- CB_LRT_Lap(num_datasets, n, 2.75*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 2.75 | Size/Power: ", res, "\n")
res <- CB_LRT_Lap(num_datasets, n, 3*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 3    | Size/Power: ", res, "\n")
res <- CB_LRT_Lap(num_datasets, n, 3.25*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 3.25 | Size/Power: ", res, "\n")
res <- CB_LRT_Lap(num_datasets, n, 3.5*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 3.5  | Size/Power: ", res, "\n")
res <- CB_LRT_Lap(num_datasets, n, 3.75*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 3.75 | Size/Power: ", res, "\n")
res <- CB_LRT_Lap(num_datasets, n, 4*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 4    | Size/Power: ", res, "\n")
res <- CB_LRT_Lap(num_datasets, n, 4.25*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 4.25 | Size/Power: ", res, "\n")
res <- CB_LRT_Lap(num_datasets, n, 4.5*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 4.5  | Size/Power: ", res, "\n")
res <- CB_LRT_Lap(num_datasets, n, 4.75*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 4.75 | Size/Power: ", res, "\n")
res <- CB_LRT_Lap(num_datasets, n, 5*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 5    | Size/Power: ", res, "\n")




cat("\n===================================================================================================\n")
cat("Parameters: p=6, k=2, sample sizes = (32, 35, 28, 40, 38, 30), mean = c*(1, 0.5, 1.4, 1.6, 1.1, 1.8), 
    variance = (12, 15, 10, 8, 14, 9), alpha=0.05\n") 
cat("=====================================================================================================\n")

num_datasets = 6
n = c(32,35,28,40,38,30)
mean = c(1,0.5,1.4,1.6,1.1,1.8)
sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9))



# ==============================================================================
cat("\n--- Distribution: Normal (N) ---\n")
# ==============================================================================
res <- CB_LRT(num_datasets, n, 0*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 0    | Size/Power: ", res, "\n")
res <- CB_LRT(num_datasets, n, 1*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 1    | Size/Power: ", res, "\n")
res <- CB_LRT(num_datasets, n, 1.25*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 1.25 | Size/Power: ", res, "\n")
res <- CB_LRT(num_datasets, n, 1.5*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 1.5  | Size/Power: ", res, "\n")
res <- CB_LRT(num_datasets, n, 1.75*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 1.75 | Size/Power: ", res, "\n")
res <- CB_LRT(num_datasets, n, 2*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 2    | Size/Power: ", res, "\n")
res <- CB_LRT(num_datasets, n, 2.25*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 2.25 | Size/Power: ", res, "\n")
res <- CB_LRT(num_datasets, n, 2.5*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 2.5  | Size/Power: ", res, "\n")
res <- CB_LRT(num_datasets, n, 2.75*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 2.75 | Size/Power: ", res, "\n")
res <- CB_LRT(num_datasets, n, 3*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 3    | Size/Power: ", res, "\n")
res <- CB_LRT(num_datasets, n, 3.25*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 3.25 | Size/Power: ", res, "\n")
res <- CB_LRT(num_datasets, n, 3.5*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 3.5  | Size/Power: ", res, "\n")
res <- CB_LRT(num_datasets, n, 3.75*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 3.75 | Size/Power: ", res, "\n")
res <- CB_LRT(num_datasets, n, 4*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 4    | Size/Power: ", res, "\n")
res <- CB_LRT(num_datasets, n, 4.25*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 4.25 | Size/Power: ", res, "\n")
res <- CB_LRT(num_datasets, n, 4.5*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 4.5  | Size/Power: ", res, "\n")
res <- CB_LRT(num_datasets, n, 4.75*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 4.75 | Size/Power: ", res, "\n")
res <- CB_LRT(num_datasets, n, 5*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 5    | Size/Power: ", res, "\n")


# ==============================================================================
cat("\n--- Distribution: Mix-1 ---\n")
# ==============================================================================
res <- CB_LRT_mix1(num_datasets, n, 0*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 0    | Size/Power: ", res, "\n")
res <- CB_LRT_mix1(num_datasets, n, 1*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 1    | Size/Power: ", res, "\n")
res <- CB_LRT_mix1(num_datasets, n, 1.25*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 1.25 | Size/Power: ", res, "\n")
res <- CB_LRT_mix1(num_datasets, n, 1.5*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 1.5  | Size/Power: ", res, "\n")
res <- CB_LRT_mix1(num_datasets, n, 1.75*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 1.75 | Size/Power: ", res, "\n")
res <- CB_LRT_mix1(num_datasets, n, 2*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 2    | Size/Power: ", res, "\n")
res <- CB_LRT_mix1(num_datasets, n, 2.25*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 2.25 | Size/Power: ", res, "\n")
res <- CB_LRT_mix1(num_datasets, n, 2.5*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 2.5  | Size/Power: ", res, "\n")
res <- CB_LRT_mix1(num_datasets, n, 2.75*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 2.75 | Size/Power: ", res, "\n")
res <- CB_LRT_mix1(num_datasets, n, 3*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 3    | Size/Power: ", res, "\n")
res <- CB_LRT_mix1(num_datasets, n, 3.25*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 3.25 | Size/Power: ", res, "\n")
res <- CB_LRT_mix1(num_datasets, n, 3.5*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 3.5  | Size/Power: ", res, "\n")
res <- CB_LRT_mix1(num_datasets, n, 3.75*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 3.75 | Size/Power: ", res, "\n")
res <- CB_LRT_mix1(num_datasets, n, 4*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 4    | Size/Power: ", res, "\n")
res <- CB_LRT_mix1(num_datasets, n, 4.25*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 4.25 | Size/Power: ", res, "\n")
res <- CB_LRT_mix1(num_datasets, n, 4.5*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 4.5  | Size/Power: ", res, "\n")
res <- CB_LRT_mix1(num_datasets, n, 4.75*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 4.75 | Size/Power: ", res, "\n")
res <- CB_LRT_mix1(num_datasets, n, 5*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 5    | Size/Power: ", res, "\n")


# ==============================================================================
cat("\n--- Distribution: Mix-2 ---\n")
# ==============================================================================
res <- CB_LRT_mix2(num_datasets, n, 0*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 0    | Size/Power: ", res, "\n")
res <- CB_LRT_mix2(num_datasets, n, 1*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 1    | Size/Power: ", res, "\n")
res <- CB_LRT_mix2(num_datasets, n, 1.25*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 1.25 | Size/Power: ", res, "\n")
res <- CB_LRT_mix2(num_datasets, n, 1.5*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 1.5  | Size/Power: ", res, "\n")
res <- CB_LRT_mix2(num_datasets, n, 1.75*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 1.75 | Size/Power: ", res, "\n")
res <- CB_LRT_mix2(num_datasets, n, 2*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 2    | Size/Power: ", res, "\n")
res <- CB_LRT_mix2(num_datasets, n, 2.25*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 2.25 | Size/Power: ", res, "\n")
res <- CB_LRT_mix2(num_datasets, n, 2.5*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 2.5  | Size/Power: ", res, "\n")
res <- CB_LRT_mix2(num_datasets, n, 2.75*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 2.75 | Size/Power: ", res, "\n")
res <- CB_LRT_mix2(num_datasets, n, 3*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 3    | Size/Power: ", res, "\n")
res <- CB_LRT_mix2(num_datasets, n, 3.25*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 3.25 | Size/Power: ", res, "\n")
res <- CB_LRT_mix2(num_datasets, n, 3.5*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 3.5  | Size/Power: ", res, "\n")
res <- CB_LRT_mix2(num_datasets, n, 3.75*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 3.75 | Size/Power: ", res, "\n")
res <- CB_LRT_mix2(num_datasets, n, 4*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 4    | Size/Power: ", res, "\n")
res <- CB_LRT_mix2(num_datasets, n, 4.25*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 4.25 | Size/Power: ", res, "\n")
res <- CB_LRT_mix2(num_datasets, n, 4.5*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 4.5  | Size/Power: ", res, "\n")
res <- CB_LRT_mix2(num_datasets, n, 4.75*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 4.75 | Size/Power: ", res, "\n")
res <- CB_LRT_mix2(num_datasets, n, 5*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 5    | Size/Power: ", res, "\n")


# ==============================================================================
cat("\n--- Distribution: Mix-3 ---\n")
# ==============================================================================
res <- CB_LRT_mix3(num_datasets, n, 0*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 0    | Size/Power: ", res, "\n")
res <- CB_LRT_mix3(num_datasets, n, 1*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 1    | Size/Power: ", res, "\n")
res <- CB_LRT_mix3(num_datasets, n, 1.25*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 1.25 | Size/Power: ", res, "\n")
res <- CB_LRT_mix3(num_datasets, n, 1.5*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 1.5  | Size/Power: ", res, "\n")
res <- CB_LRT_mix3(num_datasets, n, 1.75*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 1.75 | Size/Power: ", res, "\n")
res <- CB_LRT_mix3(num_datasets, n, 2*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 2    | Size/Power: ", res, "\n")
res <- CB_LRT_mix3(num_datasets, n, 2.25*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 2.25 | Size/Power: ", res, "\n")
res <- CB_LRT_mix3(num_datasets, n, 2.5*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 2.5  | Size/Power: ", res, "\n")
res <- CB_LRT_mix3(num_datasets, n, 2.75*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 2.75 | Size/Power: ", res, "\n")
res <- CB_LRT_mix3(num_datasets, n, 3*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 3    | Size/Power: ", res, "\n")
res <- CB_LRT_mix3(num_datasets, n, 3.25*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 3.25 | Size/Power: ", res, "\n")
res <- CB_LRT_mix3(num_datasets, n, 3.5*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 3.5  | Size/Power: ", res, "\n")
res <- CB_LRT_mix3(num_datasets, n, 3.75*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 3.75 | Size/Power: ", res, "\n")
res <- CB_LRT_mix3(num_datasets, n, 4*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 4    | Size/Power: ", res, "\n")
res <- CB_LRT_mix3(num_datasets, n, 4.25*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 4.25 | Size/Power: ", res, "\n")
res <- CB_LRT_mix3(num_datasets, n, 4.5*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 4.5  | Size/Power: ", res, "\n")
res <- CB_LRT_mix3(num_datasets, n, 4.75*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 4.75 | Size/Power: ", res, "\n")
res <- CB_LRT_mix3(num_datasets, n, 5*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 5    | Size/Power: ", res, "\n")


# ==============================================================================
cat("\n--- Distribution: Mix-4 ---\n")
# ==============================================================================
res <- CB_LRT_mix4(num_datasets, n, 0*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 0    | Size/Power: ", res, "\n")
res <- CB_LRT_mix4(num_datasets, n, 1*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 1    | Size/Power: ", res, "\n")
res <- CB_LRT_mix4(num_datasets, n, 1.25*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 1.25 | Size/Power: ", res, "\n")
res <- CB_LRT_mix4(num_datasets, n, 1.5*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 1.5  | Size/Power: ", res, "\n")
res <- CB_LRT_mix4(num_datasets, n, 1.75*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 1.75 | Size/Power: ", res, "\n")
res <- CB_LRT_mix4(num_datasets, n, 2*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 2    | Size/Power: ", res, "\n")
res <- CB_LRT_mix4(num_datasets, n, 2.25*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 2.25 | Size/Power: ", res, "\n")
res <- CB_LRT_mix4(num_datasets, n, 2.5*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 2.5  | Size/Power: ", res, "\n")
res <- CB_LRT_mix4(num_datasets, n, 2.75*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 2.75 | Size/Power: ", res, "\n")
res <- CB_LRT_mix4(num_datasets, n, 3*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 3    | Size/Power: ", res, "\n")
res <- CB_LRT_mix4(num_datasets, n, 3.25*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 3.25 | Size/Power: ", res, "\n")
res <- CB_LRT_mix4(num_datasets, n, 3.5*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 3.5  | Size/Power: ", res, "\n")
res <- CB_LRT_mix4(num_datasets, n, 3.75*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 3.75 | Size/Power: ", res, "\n")
res <- CB_LRT_mix4(num_datasets, n, 4*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 4    | Size/Power: ", res, "\n")
res <- CB_LRT_mix4(num_datasets, n, 4.25*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 4.25 | Size/Power: ", res, "\n")
res <- CB_LRT_mix4(num_datasets, n, 4.5*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 4.5  | Size/Power: ", res, "\n")
res <- CB_LRT_mix4(num_datasets, n, 4.75*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 4.75 | Size/Power: ", res, "\n")
res <- CB_LRT_mix4(num_datasets, n, 5*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 5    | Size/Power: ", res, "\n")


# ==============================================================================
cat("\n--- Distribution: Laplace ---\n")
# ==============================================================================
res <- CB_LRT_Lap(num_datasets, n, 0*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 0    | Size/Power: ", res, "\n")
res <- CB_LRT_Lap(num_datasets, n, 1*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 1    | Size/Power: ", res, "\n")
res <- CB_LRT_Lap(num_datasets, n, 1.25*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 1.25 | Size/Power: ", res, "\n")
res <- CB_LRT_Lap(num_datasets, n, 1.5*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 1.5  | Size/Power: ", res, "\n")
res <- CB_LRT_Lap(num_datasets, n, 1.75*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 1.75 | Size/Power: ", res, "\n")
res <- CB_LRT_Lap(num_datasets, n, 2*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 2    | Size/Power: ", res, "\n")
res <- CB_LRT_Lap(num_datasets, n, 2.25*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 2.25 | Size/Power: ", res, "\n")
res <- CB_LRT_Lap(num_datasets, n, 2.5*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 2.5  | Size/Power: ", res, "\n")
res <- CB_LRT_Lap(num_datasets, n, 2.75*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 2.75 | Size/Power: ", res, "\n")
res <- CB_LRT_Lap(num_datasets, n, 3*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 3    | Size/Power: ", res, "\n")
res <- CB_LRT_Lap(num_datasets, n, 3.25*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 3.25 | Size/Power: ", res, "\n")
res <- CB_LRT_Lap(num_datasets, n, 3.5*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 3.5  | Size/Power: ", res, "\n")
res <- CB_LRT_Lap(num_datasets, n, 3.75*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 3.75 | Size/Power: ", res, "\n")
res <- CB_LRT_Lap(num_datasets, n, 4*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 4    | Size/Power: ", res, "\n")
res <- CB_LRT_Lap(num_datasets, n, 4.25*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 4.25 | Size/Power: ", res, "\n")
res <- CB_LRT_Lap(num_datasets, n, 4.5*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 4.5  | Size/Power: ", res, "\n")
res <- CB_LRT_Lap(num_datasets, n, 4.75*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 4.75 | Size/Power: ", res, "\n")
res <- CB_LRT_Lap(num_datasets, n, 5*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 5    | Size/Power: ", res, "\n")





cat("\n==========================================================================\n")
cat("TABLE S7.20 REPLICATION COMPLETE\n")
cat("==========================================================================\n")





