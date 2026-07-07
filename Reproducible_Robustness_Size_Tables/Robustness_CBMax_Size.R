# ==============================================================================
# REPLICATION SCRIPT: Robustness Analysis (CBMax Size Table S7.16 of Supplementary Material)
#
# The functions called in this script (CB_Max, CB_Max_mix1, CB_Max_mix2, CB_Max_mix3, CB_Max_mix4, CB_Max_Lap)
# are designed to test the robustness of the CBMax test under deviations from normality.
#
#
# RETURN VALUE:
# Each function returns the estimated probability of a Type I Error or power, specifically.
#
# DISTRIBUTION MAPPINGS (As defined in Section: Robustness Analysis):
# - CB_Max      : Normal distribution (N).
# - CB_Max_mix1 : Mixture of N(mu_i + 4, (sigma_i + 4)^2) and N(mu_i, sigma_i^2)
#                 with mixture weights (0.15, 0.85).
# - CB_Max_mix2 : Mixture of N(mu_i + 4, (sigma_i + 7)^2) and N(mu_i, sigma_i^2)
#                 with mixture weights (0.15, 0.85).
# - CB_Max_mix3 : Mixture of N(mu_i + 4, (sigma_i + 4)^2) and N(mu_i, sigma_i^2)
#                 with mixture weights (0.30, 0.70).
# - CB_Max_mix4 : Mixture of N(mu_i + 4, (sigma_i + 7)^2) and N(mu_i, sigma_i^2)
#                 with mixture weights (0.30, 0.70).
# - CB_Max_Lap  : Laplace distribution with location parameter -2 and scale
#                 parameter 3.
#
# NOTE ON COMPUTATIONAL TIME:
# The results reported in the manuscript were obtained on a 64-core
# high-performance computing server using M = 5000 and Q = 5000.
# Under these settings, computing a single empirical size value requires
# approximately 9 minutes for each of the CB_Max, CB_Max_mix1, CB_Max_mix2, CB_Max_mix3, CB_Max_mix4, and CB_Max_Lap.
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
if (!require("parallel")) install.packages("parallel")
if (!require("doParallel")) install.packages("doParallel")
if (!require("foreach")) install.packages("foreach")
if (!require("LaplacesDemon")) install.packages("LaplacesDemon")


# Load them
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

#Normal
CB_Max <- function(num_datasets, n, mean, sd, k0, num_samples = M, outnum_repetitions = Q, num_cores = CORE) {
  cl <- makeCluster(num_cores)
  on.exit({
    stopCluster(cl)
    registerDoSEQ()
    closeAllConnections()
  }, add = TRUE)
  registerDoParallel(cl)
  results <- numeric(1)
  for (a in 1:1) {
    set.seed(123 + a)
    results[a] <- mean(foreach(rep1 = 1:outnum_repetitions, .combine = "c", .packages = c("parallel", "foreach")) %dopar% {
      set.seed(123 + a * 1000 + rep1)
      sample_data <- lapply(1:num_datasets, function(i) rnorm(n = n[i], mean = mean[i], sd = sd[i]))
      T_max <- numeric(num_samples)
      for (u in 1:num_samples) {
        bootstrap_samples <- lapply(1:num_datasets, function(i) vector("list"))
        for (j in 1:num_datasets) {
          bootstrap_samples[[j]] <- rnorm(n = n[j], mean = 0, sd = sqrt(var(sample_data[[j]])))
        }
        T <- numeric(k0)
        CB_values <- numeric(k0)
        for (p in 1:k0) {
          T[p] <- max(sapply((k0+1):num_datasets, function(i) {
            
            (mean(bootstrap_samples[[i]]) - mean(bootstrap_samples[[p]])) /
              
              sqrt(
                
                (var(bootstrap_samples[[i]]) / length(bootstrap_samples[[i]])) +
                  
                  (var(bootstrap_samples[[p]]) / length(bootstrap_samples[[p]]))
              )
          }))
        }
        T_max[u] <- max(T)
      }
      sort_T_max <- sort(T_max)
      for (p in 1:k0) {
        CB_values[p] <- max(sapply((k0+1):num_datasets, function(i) {
          
          (mean(sample_data[[i]]) - mean(sample_data[[p]])) /
            
            sqrt(
              
              (var(sample_data[[i]]) / length(sample_data[[i]])) +
                
                (var(sample_data[[p]]) / length(sample_data[[p]]))
              
            )
          
        }))
      }
      CB <- max(CB_values)
      quantile_values <- quantile(sort_T_max, probs = 0.95)
      alpha_LRT <- ifelse(quantile_values < CB, 1, 0)
      return(alpha_LRT)
    })
  }
  return((results))
}

#Mix-normal with mixture weights = c(0.15,0.85),mu = c((mean[i]+4),mean[i]),sigma = c((sd[i]+4),sd[i])
CB_Max_mix1 <- function(num_datasets, n, mean, sd, k0, num_samples = M, outnum_repetitions = Q, num_cores = CORE) {
  cl <- makeCluster(num_cores)
  on.exit({
    stopCluster(cl)
    registerDoSEQ()
    closeAllConnections()
  }, add = TRUE)
  registerDoParallel(cl)
  results <- numeric(1)
  for (a in 1:1) {
    set.seed(123 + a)
    results[a] <- mean(foreach(rep1 = 1:outnum_repetitions, .combine = "c", .packages = c("parallel", "foreach")) %dopar% {
      set.seed(123 + a * 1000 + rep1)
      library(LaplacesDemon)
      sample <- lapply(1:num_datasets, function(i) ((rnormm(n = n[i], p = c(0.15,0.85),mu = c((mean[i]+4),mean[i]),sigma = c((sd[i]+4),sd[i]))
                                                     -(0.15*(mean[i]+4)+0.85*mean[i]))/sqrt((((mean[i]+4)^2 + (sd[i]+4)^2)*0.15) +
                                                                                              ((mean[i]^2 + sd[i]^2)*0.85) - (0.15*(mean[i]+4)+0.85*mean[i])^2)))
      sample_data <- lapply(1:num_datasets, function(i) (sample[[i]]*sd[i])+mean[i])
      T_max <- numeric(num_samples)
      for (u in 1:num_samples) {
        bootstrap_samples <- lapply(1:num_datasets, function(i) vector("list"))
        for (j in 1:num_datasets) {
          bootstrap_samples[[j]] <- rnorm(n = n[j], mean = 0, sd = sqrt(var(sample_data[[j]])))
        }
        T <- numeric(k0)
        CB_values <- numeric(k0)
        for (p in 1:k0) {
          T[p] <- max(sapply((k0+1):num_datasets, function(i) {
            
            (mean(bootstrap_samples[[i]]) - mean(bootstrap_samples[[p]])) /
              
              sqrt(
                
                (var(bootstrap_samples[[i]]) / length(bootstrap_samples[[i]])) +
                  
                  (var(bootstrap_samples[[p]]) / length(bootstrap_samples[[p]]))
              )
          }))
        }
        T_max[u] <- max(T)
      }
      sort_T_max <- sort(T_max)
      for (p in 1:k0) {
        CB_values[p] <- max(sapply((k0+1):num_datasets, function(i) {
          
          (mean(sample_data[[i]]) - mean(sample_data[[p]])) /
            
            sqrt(
              
              (var(sample_data[[i]]) / length(sample_data[[i]])) +
                
                (var(sample_data[[p]]) / length(sample_data[[p]]))
              
            )
          
        }))
      }
      CB <- max(CB_values)
      quantile_values <- quantile(sort_T_max, probs = 0.95)
      alpha_LRT <- ifelse(quantile_values < CB, 1, 0)
      return(alpha_LRT)
    })
  }
  return((results))
}

#Mix-normal with mixture weights = c(0.15,0.85),mu = c((mean[i]+4),mean[i]),sigma = c((sd[i]+7),sd[i])
CB_Max_mix2 <- function(num_datasets, n, mean, sd, k0, num_samples = M, outnum_repetitions = Q, num_cores = CORE) {
  cl <- makeCluster(num_cores)
  on.exit({
    stopCluster(cl)
    registerDoSEQ()
    closeAllConnections()
  }, add = TRUE)
  registerDoParallel(cl)
  results <- numeric(1)
  for (a in 1:1) {
    set.seed(123 + a)
    results[a] <- mean(foreach(rep1 = 1:outnum_repetitions, .combine = "c", .packages = c("parallel", "foreach")) %dopar% {
      set.seed(123 + a * 1000 + rep1)
      library(LaplacesDemon)
      sample <- lapply(1:num_datasets, function(i) ((rnormm(n = n[i], p = c(0.15,0.85),mu = c((mean[i]+4),mean[i]),sigma = c((sd[i]+7),sd[i]))
                                                     -(0.15*(mean[i]+4)+0.85*mean[i]))/sqrt((((mean[i]+4)^2 + (sd[i]+7)^2)*0.15) +
                                                                                              ((mean[i]^2 + sd[i]^2)*0.85) - (0.15*(mean[i]+4)+0.85*mean[i])^2)))
      sample_data <- lapply(1:num_datasets, function(i) (sample[[i]]*sd[i])+mean[i])
      T_max <- numeric(num_samples)
      for (u in 1:num_samples) {
        bootstrap_samples <- lapply(1:num_datasets, function(i) vector("list"))
        for (j in 1:num_datasets) {
          bootstrap_samples[[j]] <- rnorm(n = n[j], mean = 0, sd = sqrt(var(sample_data[[j]])))
        }
        T <- numeric(k0)
        CB_values <- numeric(k0)
        for (p in 1:k0) {
          T[p] <- max(sapply((k0+1):num_datasets, function(i) {
            
            (mean(bootstrap_samples[[i]]) - mean(bootstrap_samples[[p]])) /
              
              sqrt(
                
                (var(bootstrap_samples[[i]]) / length(bootstrap_samples[[i]])) +
                  
                  (var(bootstrap_samples[[p]]) / length(bootstrap_samples[[p]]))
              )
          }))
        }
        T_max[u] <- max(T)
      }
      sort_T_max <- sort(T_max)
      for (p in 1:k0) {
        CB_values[p] <- max(sapply((k0+1):num_datasets, function(i) {
          
          (mean(sample_data[[i]]) - mean(sample_data[[p]])) /
            
            sqrt(
              
              (var(sample_data[[i]]) / length(sample_data[[i]])) +
                
                (var(sample_data[[p]]) / length(sample_data[[p]]))
              
            )
          
        }))
      }
      CB <- max(CB_values)
      quantile_values <- quantile(sort_T_max, probs = 0.95)
      alpha_LRT <- ifelse(quantile_values < CB, 1, 0)
      return(alpha_LRT)
    })
  }
  return((results))
}

#Mix-normal with mixture weights = c(0.30,0.70),mu = c((mean[i]+4),mean[i]),sigma = c((sd[i]+4),sd[i])
CB_Max_mix3 <- function(num_datasets, n, mean, sd, k0, num_samples = M, outnum_repetitions = Q, num_cores = CORE) {
  cl <- makeCluster(num_cores)
  on.exit({
    stopCluster(cl)
    registerDoSEQ()
    closeAllConnections()
  }, add = TRUE)
  registerDoParallel(cl)
  results <- numeric(1)
  for (a in 1:1) {
    set.seed(123 + a)
    results[a] <- mean(foreach(rep1 = 1:outnum_repetitions, .combine = "c", .packages = c("parallel", "foreach")) %dopar% {
      set.seed(123 + a * 1000 + rep1)
      library(LaplacesDemon)
      sample <- lapply(1:num_datasets, function(i) ((rnormm(n = n[i], p = c(0.30,0.70),mu = c((mean[i]+4),mean[i]),sigma = c((sd[i]+4),sd[i]))
                                                     -(0.30*(mean[i]+4)+0.70*mean[i]))/sqrt((((mean[i]+4)^2 + (sd[i]+4)^2)*0.30) +
                                                                                              ((mean[i]^2 + sd[i]^2)*0.70) - (0.30*(mean[i]+4)+0.70*mean[i])^2)))
      sample_data <- lapply(1:num_datasets, function(i) (sample[[i]]*sd[i])+mean[i])
      T_max <- numeric(num_samples)
      for (u in 1:num_samples) {
        bootstrap_samples <- lapply(1:num_datasets, function(i) vector("list"))
        for (j in 1:num_datasets) {
          bootstrap_samples[[j]] <- rnorm(n = n[j], mean = 0, sd = sqrt(var(sample_data[[j]])))
        }
        T <- numeric(k0)
        CB_values <- numeric(k0)
        for (p in 1:k0) {
          T[p] <- max(sapply((k0+1):num_datasets, function(i) {
            
            (mean(bootstrap_samples[[i]]) - mean(bootstrap_samples[[p]])) /
              
              sqrt(
                
                (var(bootstrap_samples[[i]]) / length(bootstrap_samples[[i]])) +
                  
                  (var(bootstrap_samples[[p]]) / length(bootstrap_samples[[p]]))
              )
          }))
        }
        T_max[u] <- max(T)
      }
      sort_T_max <- sort(T_max)
      for (p in 1:k0) {
        CB_values[p] <- max(sapply((k0+1):num_datasets, function(i) {
          
          (mean(sample_data[[i]]) - mean(sample_data[[p]])) /
            
            sqrt(
              
              (var(sample_data[[i]]) / length(sample_data[[i]])) +
                
                (var(sample_data[[p]]) / length(sample_data[[p]]))
              
            )
          
        }))
      }
      CB <- max(CB_values)
      quantile_values <- quantile(sort_T_max, probs = 0.95)
      alpha_LRT <- ifelse(quantile_values < CB, 1, 0)
      return(alpha_LRT)
    })
  }
  return((results))
}

#Mix-normal with mixture weights = c(0.30,0.70),mu = c((mean[i]+4),mean[i]),sigma = c((sd[i]+7),sd[i])
CB_Max_mix4 <- function(num_datasets, n, mean, sd, k0, num_samples = M, outnum_repetitions = Q, num_cores = CORE) {
  cl <- makeCluster(num_cores)
  on.exit({
    stopCluster(cl)
    registerDoSEQ()
    closeAllConnections()
  }, add = TRUE)
  registerDoParallel(cl)
  results <- numeric(1)
  for (a in 1:1) {
    set.seed(123 + a)
    results[a] <- mean(foreach(rep1 = 1:outnum_repetitions, .combine = "c", .packages = c("parallel", "foreach")) %dopar% {
      set.seed(123 + a * 1000 + rep1)
      library(LaplacesDemon)
      sample <- lapply(1:num_datasets, function(i) ((rnormm(n = n[i], p = c(0.30,0.70),mu = c((mean[i]+4),mean[i]),sigma = c((sd[i]+7),sd[i]))
                                                     -(0.30*(mean[i]+4)+0.70*mean[i]))/sqrt((((mean[i]+4)^2 + (sd[i]+7)^2)*0.30) +
                                                                                              ((mean[i]^2 + sd[i]^2)*0.70) - (0.30*(mean[i]+4)+0.70*mean[i])^2)))
      sample_data <- lapply(1:num_datasets, function(i) (sample[[i]]*sd[i])+mean[i])
      T_max <- numeric(num_samples)
      for (u in 1:num_samples) {
        bootstrap_samples <- lapply(1:num_datasets, function(i) vector("list"))
        for (j in 1:num_datasets) {
          bootstrap_samples[[j]] <- rnorm(n = n[j], mean = 0, sd = sqrt(var(sample_data[[j]])))
        }
        T <- numeric(k0)
        CB_values <- numeric(k0)
        for (p in 1:k0) {
          T[p] <- max(sapply((k0+1):num_datasets, function(i) {
            
            (mean(bootstrap_samples[[i]]) - mean(bootstrap_samples[[p]])) /
              
              sqrt(
                
                (var(bootstrap_samples[[i]]) / length(bootstrap_samples[[i]])) +
                  
                  (var(bootstrap_samples[[p]]) / length(bootstrap_samples[[p]]))
              )
          }))
        }
        T_max[u] <- max(T)
      }
      sort_T_max <- sort(T_max)
      for (p in 1:k0) {
        CB_values[p] <- max(sapply((k0+1):num_datasets, function(i) {
          
          (mean(sample_data[[i]]) - mean(sample_data[[p]])) /
            
            sqrt(
              
              (var(sample_data[[i]]) / length(sample_data[[i]])) +
                
                (var(sample_data[[p]]) / length(sample_data[[p]]))
              
            )
          
        }))
      }
      CB <- max(CB_values)
      quantile_values <- quantile(sort_T_max, probs = 0.95)
      alpha_LRT <- ifelse(quantile_values < CB, 1, 0)
      return(alpha_LRT)
    })
  }
  return((results))
}

#Laplace with location = -2, scale = 3
CB_Max_Lap <- function(num_datasets, n, mean, sd, k0, num_samples = M, outnum_repetitions = Q, num_cores = CORE) {
  cl <- makeCluster(num_cores)
  on.exit({
    stopCluster(cl)
    registerDoSEQ()
    closeAllConnections()
  }, add = TRUE)
  registerDoParallel(cl)
  results <- numeric(1)
  for (a in 1:1) {
    set.seed(123 + a)
    results[a] <- mean(foreach(rep1 = 1:outnum_repetitions, .combine = "c", .packages = c("parallel", "foreach")) %dopar% {
      set.seed(123 + a * 1000 + rep1)
      library(LaplacesDemon)
      sample <- lapply(1:num_datasets, function(i) ((rlaplace(n = n[i], location = -2, scale = 3)+2)/sqrt(18)))
      sample_data <- lapply(1:num_datasets, function(i) (sample[[i]]*sd[i])+mean[i])
      T_max <- numeric(num_samples)
      for (u in 1:num_samples) {
        bootstrap_samples <- lapply(1:num_datasets, function(i) vector("list"))
        for (j in 1:num_datasets) {
          bootstrap_samples[[j]] <- rnorm(n = n[j], mean = 0, sd = sqrt(var(sample_data[[j]])))
        }
        T <- numeric(k0)
        CB_values <- numeric(k0)
        for (p in 1:k0) {
          T[p] <- max(sapply((k0+1):num_datasets, function(i) {
            
            (mean(bootstrap_samples[[i]]) - mean(bootstrap_samples[[p]])) /
              
              sqrt(
                
                (var(bootstrap_samples[[i]]) / length(bootstrap_samples[[i]])) +
                  
                  (var(bootstrap_samples[[p]]) / length(bootstrap_samples[[p]]))
              )
          }))
        }
        T_max[u] <- max(T)
      }
      sort_T_max <- sort(T_max)
      for (p in 1:k0) {
        CB_values[p] <- max(sapply((k0+1):num_datasets, function(i) {
          
          (mean(sample_data[[i]]) - mean(sample_data[[p]])) /
            
            sqrt(
              
              (var(sample_data[[i]]) / length(sample_data[[i]])) +
                
                (var(sample_data[[p]]) / length(sample_data[[p]]))
              
            )
          
        }))
      }
      CB <- max(CB_values)
      quantile_values <- quantile(sort_T_max, probs = 0.95)
      alpha_LRT <- ifelse(quantile_values < CB, 1, 0)
      return(alpha_LRT)
    })
  }
  return((results))
}



cat("\n==========================================================================\n")
cat("STARTING ROBUSTNESS ANALYSIS: S7.16 of Supplementary Material (CBMax under Non-Normal Distributions)\n")
cat("Parameters: p=6, k=2, mean=(1,1,1,1,1,1), variance=(12,15,10,8,14,9), alpha=0.05\n")
cat("==========================================================================\n")


# ==============================================================================
cat("\n--- Row 1: Sample Sizes = (7, 9, 6, 5, 10, 8) ---\n")
# ==============================================================================
res <- CB_Max(num_datasets = 6, n = c(7,9,6,5,10,8), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Normal (N) | Size: ", res, "\n")

res <- CB_Max_mix1(num_datasets = 6, n = c(7,9,6,5,10,8), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-1      | Size: ", res, "\n")

res <- CB_Max_mix2(num_datasets = 6, n = c(7,9,6,5,10,8), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-2      | Size: ", res, "\n")

res <- CB_Max_mix3(num_datasets = 6, n = c(7,9,6,5,10,8), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-3      | Size: ", res, "\n")

res <- CB_Max_mix4(num_datasets = 6, n = c(7,9,6,5,10,8), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-4      | Size: ", res, "\n")

res <- CB_Max_Lap(num_datasets = 6, n = c(7,9,6,5,10,8), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Laplace    | Size: ", res, "\n")


# ==============================================================================
cat("\n--- Row 2: Sample Sizes = (10, 10, 10, 10, 10, 10) ---\n")
# ==============================================================================
res <- CB_Max(num_datasets = 6, n = c(10,10,10,10,10,10), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Normal (N) | Size: ", res, "\n")

res <- CB_Max_mix1(num_datasets = 6, n = c(10,10,10,10,10,10), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-1      | Size: ", res, "\n")

res <- CB_Max_mix2(num_datasets = 6, n = c(10,10,10,10,10,10), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-2      | Size: ", res, "\n")

res <- CB_Max_mix3(num_datasets = 6, n = c(10,10,10,10,10,10), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-3      | Size: ", res, "\n")

res <- CB_Max_mix4(num_datasets = 6, n = c(10,10,10,10,10,10), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-4      | Size: ", res, "\n")

res <- CB_Max_Lap(num_datasets = 6, n = c(10,10,10,10,10,10), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Laplace    | Size: ", res, "\n")


# ==============================================================================
cat("\n--- Row 3: Sample Sizes = (12, 10, 15, 18, 16, 8) ---\n")
# ==============================================================================
res <- CB_Max(num_datasets = 6, n = c(12,10,15,18,16,8), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Normal (N) | Size: ", res, "\n")

res <- CB_Max_mix1(num_datasets = 6, n = c(12,10,15,18,16,8), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-1      | Size: ", res, "\n")

res <- CB_Max_mix2(num_datasets = 6, n = c(12,10,15,18,16,8), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-2      | Size: ", res, "\n")

res <- CB_Max_mix3(num_datasets = 6, n = c(12,10,15,18,16,8), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-3      | Size: ", res, "\n")

res <- CB_Max_mix4(num_datasets = 6, n = c(12,10,15,18,16,8), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-4      | Size: ", res, "\n")

res <- CB_Max_Lap(num_datasets = 6, n = c(12,10,15,18,16,8), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Laplace    | Size: ", res, "\n")


# ==============================================================================
cat("\n--- Row 4: Sample Sizes = (11, 9, 32, 40, 35, 28) ---\n")
# ==============================================================================
res <- CB_Max(num_datasets = 6, n = c(11,9,32,40,35,28), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Normal (N) | Size: ", res, "\n")

res <- CB_Max_mix1(num_datasets = 6, n = c(11,9,32,40,35,28), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-1      | Size: ", res, "\n")

res <- CB_Max_mix2(num_datasets = 6, n = c(11,9,32,40,35,28), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-2      | Size: ", res, "\n")

res <- CB_Max_mix3(num_datasets = 6, n = c(11,9,32,40,35,28), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-3      | Size: ", res, "\n")

res <- CB_Max_mix4(num_datasets = 6, n = c(11,9,32,40,35,28), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-4      | Size: ", res, "\n")

res <- CB_Max_Lap(num_datasets = 6, n = c(11,9,32,40,35,28), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Laplace    | Size: ", res, "\n")


# ==============================================================================
cat("\n--- Row 5: Sample Sizes = (32, 35, 28, 40, 38, 30) ---\n")
# ==============================================================================
res <- CB_Max(num_datasets = 6, n = c(32,35,28,40,38,30), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Normal (N) | Size: ", res, "\n")

res <- CB_Max_mix1(num_datasets = 6, n = c(32,35,28,40,38,30), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-1      | Size: ", res, "\n")

res <- CB_Max_mix2(num_datasets = 6, n = c(32,35,28,40,38,30), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-2      | Size: ", res, "\n")

res <- CB_Max_mix3(num_datasets = 6, n = c(32,35,28,40,38,30), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-3      | Size: ", res, "\n")

res <- CB_Max_mix4(num_datasets = 6, n = c(32,35,28,40,38,30), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-4      | Size: ", res, "\n")

res <- CB_Max_Lap(num_datasets = 6, n = c(32,35,28,40,38,30), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Laplace    | Size: ", res, "\n")


# ==============================================================================
cat("\n--- Row 6: Sample Sizes = (42, 35, 12, 14, 10, 8) ---\n")
# ==============================================================================
res <- CB_Max(num_datasets = 6, n = c(42,35,12,14,10,8), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Normal (N) | Size: ", res, "\n")

res <- CB_Max_mix1(num_datasets = 6, n = c(42,35,12,14,10,8), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-1      | Size: ", res, "\n")

res <- CB_Max_mix2(num_datasets = 6, n = c(42,35,12,14,10,8), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-2      | Size: ", res, "\n")

res <- CB_Max_mix3(num_datasets = 6, n = c(42,35,12,14,10,8), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-3      | Size: ", res, "\n")

res <- CB_Max_mix4(num_datasets = 6, n = c(42,35,12,14,10,8), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-4      | Size: ", res, "\n")

res <- CB_Max_Lap(num_datasets = 6, n = c(42,35,12,14,10,8), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Laplace    | Size: ", res, "\n")


# ==============================================================================
cat("\n--- Row 7: Sample Sizes = (35, 12, 42, 10, 30, 8) ---\n")
# ==============================================================================
res <- CB_Max(num_datasets = 6, n = c(35,12,42,10,30,8), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Normal (N) | Size: ", res, "\n")

res <- CB_Max_mix1(num_datasets = 6, n = c(35,12,42,10,30,8), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-1      | Size: ", res, "\n")

res <- CB_Max_mix2(num_datasets = 6, n = c(35,12,42,10,30,8), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-2      | Size: ", res, "\n")

res <- CB_Max_mix3(num_datasets = 6, n = c(35,12,42,10,30,8), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-3      | Size: ", res, "\n")

res <- CB_Max_mix4(num_datasets = 6, n = c(35,12,42,10,30,8), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-4      | Size: ", res, "\n")

res <- CB_Max_Lap(num_datasets = 6, n = c(35,12,42,10,30,8), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Laplace    | Size: ", res, "\n")


# ==============================================================================
cat("\n--- Row 8: Sample Sizes = (11, 32, 8, 38, 6, 28) ---\n")
# ==============================================================================
res <- CB_Max(num_datasets = 6, n = c(11,32,8,38,6,28), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Normal (N) | Size: ", res, "\n")

res <- CB_Max_mix1(num_datasets = 6, n = c(11,32,8,38,6,28), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-1      | Size: ", res, "\n")

res <- CB_Max_mix2(num_datasets = 6, n = c(11,32,8,38,6,28), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-2      | Size: ", res, "\n")

res <- CB_Max_mix3(num_datasets = 6, n = c(11,32,8,38,6,28), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-3      | Size: ", res, "\n")

res <- CB_Max_mix4(num_datasets = 6, n = c(11,32,8,38,6,28), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-4      | Size: ", res, "\n")

res <- CB_Max_Lap(num_datasets = 6, n = c(11,32,8,38,6,28), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Laplace    | Size: ", res, "\n")


# ==============================================================================
cat("\n--- Row 9: Sample Sizes = (25, 25, 25, 25, 25, 25) ---\n")
# ==============================================================================
res <- CB_Max(num_datasets = 6, n = c(25,25,25,25,25,25), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Normal (N) | Size: ", res, "\n")

res <- CB_Max_mix1(num_datasets = 6, n = c(25,25,25,25,25,25), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-1      | Size: ", res, "\n")

res <- CB_Max_mix2(num_datasets = 6, n = c(25,25,25,25,25,25), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-2      | Size: ", res, "\n")

res <- CB_Max_mix3(num_datasets = 6, n = c(25,25,25,25,25,25), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-3      | Size: ", res, "\n")

res <- CB_Max_mix4(num_datasets = 6, n = c(25,25,25,25,25,25), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-4      | Size: ", res, "\n")

res <- CB_Max_Lap(num_datasets = 6, n = c(25,25,25,25,25,25), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Laplace    | Size: ", res, "\n")


# ==============================================================================
cat("\n--- Row 10: Sample Sizes = (12, 15, 17, 20, 25, 32) ---\n")
# ==============================================================================
res <- CB_Max(num_datasets = 6, n = c(12,15,17,20,25,32), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Normal (N) | Size: ", res, "\n")

res <- CB_Max_mix1(num_datasets = 6, n = c(12,15,17,20,25,32), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-1      | Size: ", res, "\n")

res <- CB_Max_mix2(num_datasets = 6, n = c(12,15,17,20,25,32), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-2      | Size: ", res, "\n")

res <- CB_Max_mix3(num_datasets = 6, n = c(12,15,17,20,25,32), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-3      | Size: ", res, "\n")

res <- CB_Max_mix4(num_datasets = 6, n = c(12,15,17,20,25,32), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-4      | Size: ", res, "\n")

res <- CB_Max_Lap(num_datasets = 6, n = c(12,15,17,20,25,32), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Laplace    | Size: ", res, "\n")


# ==============================================================================
cat("\n--- Row 11: Sample Sizes = (30, 28, 26, 20, 15, 10) ---\n")
# ==============================================================================
res <- CB_Max(num_datasets = 6, n = c(30,28,26,20,15,10), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Normal (N) | Size: ", res, "\n")

res <- CB_Max_mix1(num_datasets = 6, n = c(30,28,26,20,15,10), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-1      | Size: ", res, "\n")

res <- CB_Max_mix2(num_datasets = 6, n = c(30,28,26,20,15,10), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-2      | Size: ", res, "\n")

res <- CB_Max_mix3(num_datasets = 6, n = c(30,28,26,20,15,10), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-3      | Size: ", res, "\n")

res <- CB_Max_mix4(num_datasets = 6, n = c(30,28,26,20,15,10), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-4      | Size: ", res, "\n")

res <- CB_Max_Lap(num_datasets = 6, n = c(30,28,26,20,15,10), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Laplace    | Size: ", res, "\n")


# ==============================================================================
cat("\n--- Row 12: Sample Sizes = (8, 35, 32, 30, 42, 38) ---\n")
# ==============================================================================
res <- CB_Max(num_datasets = 6, n = c(8,35,32,30,42,38), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Normal (N) | Size: ", res, "\n")

res <- CB_Max_mix1(num_datasets = 6, n = c(8,35,32,30,42,38), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-1      | Size: ", res, "\n")

res <- CB_Max_mix2(num_datasets = 6, n = c(8,35,32,30,42,38), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-2      | Size: ", res, "\n")

res <- CB_Max_mix3(num_datasets = 6, n = c(8,35,32,30,42,38), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-3      | Size: ", res, "\n")

res <- CB_Max_mix4(num_datasets = 6, n = c(8,35,32,30,42,38), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Mix-4      | Size: ", res, "\n")

res <- CB_Max_Lap(num_datasets = 6, n = c(8,35,32,30,42,38), mean = rep(1,6), sd = c(sqrt(12), sqrt(15), sqrt(10), sqrt(8), sqrt(14), sqrt(9)), k0 = 2)
cat("Laplace    | Size: ", res, "\n")

cat("\n==========================================================================\n")
cat("TABLE S7.16 REPLICATION COMPLETE\n")
cat("==========================================================================\n")
