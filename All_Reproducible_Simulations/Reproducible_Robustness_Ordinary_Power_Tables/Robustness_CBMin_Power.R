# ==============================================================================================================
# REPLICATION SCRIPT: Robustness Analysis (Table S7.22 of Supplementary Material --- CBMin Ordinary Power Table)
#
# The functions called in this script (CB_Min, CB_Min_mix1, CB_Min_mix2, CB_Min_mix3, CB_Min_mix4, CB_Min_Lap)
# are designed to test the robustness of the CBMin test under deviations from normality.
#
#
# RETURN VALUE:
# Each function returns the estimated probability of Type I Error or power.
#
# DISTRIBUTION MAPPINGS (As defined in Section: Robustness Analysis):
# - CB_Min      : Normal distribution (N).
# - CB_Min_mix1 : Mixture of N(mu_i + 4, (sigma_i + 4)^2) and N(mu_i, sigma_i^2)
#                 with mixture weights (0.15, 0.85).
# - CB_Min_mix2 : Mixture of N(mu_i + 4, (sigma_i + 7)^2) and N(mu_i, sigma_i^2)
#                 with mixture weights (0.15, 0.85).
# - CB_Min_mix3 : Mixture of N(mu_i + 4, (sigma_i + 4)^2) and N(mu_i, sigma_i^2)
#                 with mixture weights (0.30, 0.70).
# - CB_Min_mix4 : Mixture of N(mu_i + 4, (sigma_i + 7)^2) and N(mu_i, sigma_i^2)
#                 with mixture weights (0.30, 0.70).
# - CB_Min_Lap  : Laplace distribution with location parameter -2 and scale
#                 parameter 3.

# NOTE ON COMPUTATIONAL TIME:
# The results reported in the manuscript were obtained on a 64-core
# high-performance computing server using M = 5000 and Q = 5000.
# Under these settings, computing a single empirical size value requires
# approximately 9 minutes for the CBMin, CB_Min_mix1, CB_Min_mix2, CB_Min_mix3, CB_Min_mix4, and CB_Min_Lap.
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
CB_Min <- function(num_datasets, n, mean, sd, k0, num_samples = M, outnum_repetitions = Q, num_cores = CORE) {
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
      T_min <- numeric(num_samples)
      for (u in 1:num_samples) {
        bootstrap_samples <- lapply(1:num_datasets, function(i) vector("list"))
        for (j in 1:num_datasets) {
          bootstrap_samples[[j]] <- rnorm(n = n[j], mean = 0, sd = sqrt(var(sample_data[[j]])))
        }
        T <- numeric(k0)
        CB_values <- numeric(k0)
        for (p in 1:k0) {
          T[p] <- min(sapply((k0+1):num_datasets, function(i) {
            
            (mean(bootstrap_samples[[i]]) - mean(bootstrap_samples[[p]])) /
              
              sqrt(
                
                (var(bootstrap_samples[[i]]) / length(bootstrap_samples[[i]])) +
                  
                  (var(bootstrap_samples[[p]]) / length(bootstrap_samples[[p]]))
              )
          }))
        }
        T_min[u] <- min(T)
      }
      sort_T_min <- sort(T_min)
      for (p in 1:k0) {
        CB_values[p] <- min(sapply((k0+1):num_datasets, function(i) {
          
          (mean(sample_data[[i]]) - mean(sample_data[[p]])) /
            
            sqrt(
              
              (var(sample_data[[i]]) / length(sample_data[[i]])) +
                
                (var(sample_data[[p]]) / length(sample_data[[p]]))
              
            )
          
        }))
      }
      CB <- min(CB_values)
      quantile_values <- quantile(sort_T_min, probs = 0.95)
      alpha_CBmin <- ifelse(quantile_values < CB, 1, 0)
      return(alpha_CBmin)
    })
  }
  return((results))
}

#Mix-normal with mixture weights = c(0.15,0.85),mu = c((mean[i]+4),mean[i]),sigma = c((sd[i]+4),sd[i])
CB_Min_mix1 <- function(num_datasets, n, mean, sd, k0, num_samples = M, outnum_repetitions = Q, num_cores = CORE) {
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
      T_min <- numeric(num_samples)
      for (u in 1:num_samples) {
        bootstrap_samples <- lapply(1:num_datasets, function(i) vector("list"))
        for (j in 1:num_datasets) {
          bootstrap_samples[[j]] <- rnorm(n = n[j], mean = 0, sd = sqrt(var(sample_data[[j]])))
        }
        T <- numeric(k0)
        CB_values <- numeric(k0)
        for (p in 1:k0) {
          T[p] <- min(sapply((k0+1):num_datasets, function(i) {
            
            (mean(bootstrap_samples[[i]]) - mean(bootstrap_samples[[p]])) /
              
              sqrt(
                
                (var(bootstrap_samples[[i]]) / length(bootstrap_samples[[i]])) +
                  
                  (var(bootstrap_samples[[p]]) / length(bootstrap_samples[[p]]))
              )
          }))
        }
        T_min[u] <- min(T)
      }
      sort_T_min <- sort(T_min)
      for (p in 1:k0) {
        CB_values[p] <- min(sapply((k0+1):num_datasets, function(i) {
          
          (mean(sample_data[[i]]) - mean(sample_data[[p]])) /
            
            sqrt(
              
              (var(sample_data[[i]]) / length(sample_data[[i]])) +
                
                (var(sample_data[[p]]) / length(sample_data[[p]]))
              
            )
          
        }))
      }
      CB <- min(CB_values)
      quantile_values <- quantile(sort_T_min, probs = 0.95)
      alpha_CBmin <- ifelse(quantile_values < CB, 1, 0)
      return(alpha_CBmin)
    })
  }
  return((results))
}

#Mix-normal with mixture weights = c(0.15,0.85),mu = c((mean[i]+4),mean[i]),sigma = c((sd[i]+7),sd[i])
CB_Min_mix2 <- function(num_datasets, n, mean, sd, k0, num_samples = M, outnum_repetitions = Q, num_cores = CORE) {
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
      T_min <- numeric(num_samples)
      for (u in 1:num_samples) {
        bootstrap_samples <- lapply(1:num_datasets, function(i) vector("list"))
        for (j in 1:num_datasets) {
          bootstrap_samples[[j]] <- rnorm(n = n[j], mean = 0, sd = sqrt(var(sample_data[[j]])))
        }
        T <- numeric(k0)
        CB_values <- numeric(k0)
        for (p in 1:k0) {
          T[p] <- min(sapply((k0+1):num_datasets, function(i) {
            
            (mean(bootstrap_samples[[i]]) - mean(bootstrap_samples[[p]])) /
              
              sqrt(
                
                (var(bootstrap_samples[[i]]) / length(bootstrap_samples[[i]])) +
                  
                  (var(bootstrap_samples[[p]]) / length(bootstrap_samples[[p]]))
              )
          }))
        }
        T_min[u] <- min(T)
      }
      sort_T_min <- sort(T_min)
      for (p in 1:k0) {
        CB_values[p] <- min(sapply((k0+1):num_datasets, function(i) {
          
          (mean(sample_data[[i]]) - mean(sample_data[[p]])) /
            
            sqrt(
              
              (var(sample_data[[i]]) / length(sample_data[[i]])) +
                
                (var(sample_data[[p]]) / length(sample_data[[p]]))
              
            )
          
        }))
      }
      CB <- min(CB_values)
      quantile_values <- quantile(sort_T_min, probs = 0.95)
      alpha_CBmin <- ifelse(quantile_values < CB, 1, 0)
      return(alpha_CBmin)
    })
  }
  return((results))
}

#Mix-normal with mixture weights = c(0.30,0.70),mu = c((mean[i]+4),mean[i]),sigma = c((sd[i]+4),sd[i])
CB_Min_mix3 <- function(num_datasets, n, mean, sd, k0, num_samples = M, outnum_repetitions = Q, num_cores = CORE) {
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
      T_min <- numeric(num_samples)
      for (u in 1:num_samples) {
        bootstrap_samples <- lapply(1:num_datasets, function(i) vector("list"))
        for (j in 1:num_datasets) {
          bootstrap_samples[[j]] <- rnorm(n = n[j], mean = 0, sd = sqrt(var(sample_data[[j]])))
        }
        T <- numeric(k0)
        CB_values <- numeric(k0)
        for (p in 1:k0) {
          T[p] <- min(sapply((k0+1):num_datasets, function(i) {
            
            (mean(bootstrap_samples[[i]]) - mean(bootstrap_samples[[p]])) /
              
              sqrt(
                
                (var(bootstrap_samples[[i]]) / length(bootstrap_samples[[i]])) +
                  
                  (var(bootstrap_samples[[p]]) / length(bootstrap_samples[[p]]))
              )
          }))
        }
        T_min[u] <- min(T)
      }
      sort_T_min <- sort(T_min)
      for (p in 1:k0) {
        CB_values[p] <- min(sapply((k0+1):num_datasets, function(i) {
          
          (mean(sample_data[[i]]) - mean(sample_data[[p]])) /
            
            sqrt(
              
              (var(sample_data[[i]]) / length(sample_data[[i]])) +
                
                (var(sample_data[[p]]) / length(sample_data[[p]]))
              
            )
          
        }))
      }
      CB <- min(CB_values)
      quantile_values <- quantile(sort_T_min, probs = 0.95)
      alpha_CBmin <- ifelse(quantile_values < CB, 1, 0)
      return(alpha_CBmin)
    })
  }
  return((results))
}

#Mix-normal with mixture weights = c(0.30,0.70),mu = c((mean[i]+4),mean[i]),sigma = c((sd[i]+7),sd[i])
CB_Min_mix4 <- function(num_datasets, n, mean, sd, k0, num_samples = M, outnum_repetitions = Q, num_cores = CORE) {
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
      T_min <- numeric(num_samples)
      for (u in 1:num_samples) {
        bootstrap_samples <- lapply(1:num_datasets, function(i) vector("list"))
        for (j in 1:num_datasets) {
          bootstrap_samples[[j]] <- rnorm(n = n[j], mean = 0, sd = sqrt(var(sample_data[[j]])))
        }
        T <- numeric(k0)
        CB_values <- numeric(k0)
        for (p in 1:k0) {
          T[p] <- min(sapply((k0+1):num_datasets, function(i) {
            
            (mean(bootstrap_samples[[i]]) - mean(bootstrap_samples[[p]])) /
              
              sqrt(
                
                (var(bootstrap_samples[[i]]) / length(bootstrap_samples[[i]])) +
                  
                  (var(bootstrap_samples[[p]]) / length(bootstrap_samples[[p]]))
              )
          }))
        }
        T_min[u] <- min(T)
      }
      sort_T_min <- sort(T_min)
      for (p in 1:k0) {
        CB_values[p] <- min(sapply((k0+1):num_datasets, function(i) {
          
          (mean(sample_data[[i]]) - mean(sample_data[[p]])) /
            
            sqrt(
              
              (var(sample_data[[i]]) / length(sample_data[[i]])) +
                
                (var(sample_data[[p]]) / length(sample_data[[p]]))
              
            )
          
        }))
      }
      CB <- min(CB_values)
      quantile_values <- quantile(sort_T_min, probs = 0.95)
      alpha_CBmin <- ifelse(quantile_values < CB, 1, 0)
      return(alpha_CBmin)
    })
  }
  return((results))
}

#Laplace with location = -2, scale = 3
CB_Min_Lap <- function(num_datasets, n, mean, sd, k0, num_samples = M, outnum_repetitions = Q, num_cores = CORE) {
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
      T_min <- numeric(num_samples)
      for (u in 1:num_samples) {
        bootstrap_samples <- lapply(1:num_datasets, function(i) vector("list"))
        for (j in 1:num_datasets) {
          bootstrap_samples[[j]] <- rnorm(n = n[j], mean = 0, sd = sqrt(var(sample_data[[j]])))
        }
        T <- numeric(k0)
        CB_values <- numeric(k0)
        for (p in 1:k0) {
          T[p] <- min(sapply((k0+1):num_datasets, function(i) {
            
            (mean(bootstrap_samples[[i]]) - mean(bootstrap_samples[[p]])) /
              
              sqrt(
                
                (var(bootstrap_samples[[i]]) / length(bootstrap_samples[[i]])) +
                  
                  (var(bootstrap_samples[[p]]) / length(bootstrap_samples[[p]]))
              )
          }))
        }
        T_min[u] <- min(T)
      }
      sort_T_min <- sort(T_min)
      for (p in 1:k0) {
        CB_values[p] <- min(sapply((k0+1):num_datasets, function(i) {
          
          (mean(sample_data[[i]]) - mean(sample_data[[p]])) /
            
            sqrt(
              
              (var(sample_data[[i]]) / length(sample_data[[i]])) +
                
                (var(sample_data[[p]]) / length(sample_data[[p]]))
              
            )
          
        }))
      }
      CB <- min(CB_values)
      quantile_values <- quantile(sort_T_min, probs = 0.95)
      alpha_CBmin <- ifelse(quantile_values < CB, 1, 0)
      return(alpha_CBmin)
    })
  }
  return((results))
}



cat("\n==========================================================================\n")
cat("STARTING ROBUSTNESS ANALYSIS CBMin: Table S7.22 of Supplementary Material (Ordinary Power Tabel)\n")
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
res <- CB_Min(num_datasets, n, 0*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 0    | Size/Power: ", res, "\n")
res <- CB_Min(num_datasets, n, 1*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 1    | Size/Power: ", res, "\n")
res <- CB_Min(num_datasets, n, 1.25*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 1.25 | Size/Power: ", res, "\n")
res <- CB_Min(num_datasets, n, 1.5*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 1.5  | Size/Power: ", res, "\n")
res <- CB_Min(num_datasets, n, 1.75*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 1.75 | Size/Power: ", res, "\n")
res <- CB_Min(num_datasets, n, 2*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 2    | Size/Power: ", res, "\n")
res <- CB_Min(num_datasets, n, 2.25*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 2.25 | Size/Power: ", res, "\n")
res <- CB_Min(num_datasets, n, 2.5*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 2.5  | Size/Power: ", res, "\n")
res <- CB_Min(num_datasets, n, 2.75*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 2.75 | Size/Power: ", res, "\n")
res <- CB_Min(num_datasets, n, 3*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 3    | Size/Power: ", res, "\n")
res <- CB_Min(num_datasets, n, 3.25*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 3.25 | Size/Power: ", res, "\n")
res <- CB_Min(num_datasets, n, 3.5*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 3.5  | Size/Power: ", res, "\n")
res <- CB_Min(num_datasets, n, 3.75*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 3.75 | Size/Power: ", res, "\n")
res <- CB_Min(num_datasets, n, 4*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 4    | Size/Power: ", res, "\n")
res <- CB_Min(num_datasets, n, 4.25*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 4.25 | Size/Power: ", res, "\n")
res <- CB_Min(num_datasets, n, 4.5*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 4.5  | Size/Power: ", res, "\n")
res <- CB_Min(num_datasets, n, 4.75*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 4.75 | Size/Power: ", res, "\n")
res <- CB_Min(num_datasets, n, 5*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 5    | Size/Power: ", res, "\n")


# ==============================================================================
cat("\n--- Distribution: Mix-1 ---\n")
# ==============================================================================
res <- CB_Min_mix1(num_datasets, n, 0*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 0    | Size/Power: ", res, "\n")
res <- CB_Min_mix1(num_datasets, n, 1*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 1    | Size/Power: ", res, "\n")
res <- CB_Min_mix1(num_datasets, n, 1.25*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 1.25 | Size/Power: ", res, "\n")
res <- CB_Min_mix1(num_datasets, n, 1.5*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 1.5  | Size/Power: ", res, "\n")
res <- CB_Min_mix1(num_datasets, n, 1.75*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 1.75 | Size/Power: ", res, "\n")
res <- CB_Min_mix1(num_datasets, n, 2*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 2    | Size/Power: ", res, "\n")
res <- CB_Min_mix1(num_datasets, n, 2.25*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 2.25 | Size/Power: ", res, "\n")
res <- CB_Min_mix1(num_datasets, n, 2.5*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 2.5  | Size/Power: ", res, "\n")
res <- CB_Min_mix1(num_datasets, n, 2.75*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 2.75 | Size/Power: ", res, "\n")
res <- CB_Min_mix1(num_datasets, n, 3*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 3    | Size/Power: ", res, "\n")
res <- CB_Min_mix1(num_datasets, n, 3.25*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 3.25 | Size/Power: ", res, "\n")
res <- CB_Min_mix1(num_datasets, n, 3.5*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 3.5  | Size/Power: ", res, "\n")
res <- CB_Min_mix1(num_datasets, n, 3.75*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 3.75 | Size/Power: ", res, "\n")
res <- CB_Min_mix1(num_datasets, n, 4*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 4    | Size/Power: ", res, "\n")
res <- CB_Min_mix1(num_datasets, n, 4.25*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 4.25 | Size/Power: ", res, "\n")
res <- CB_Min_mix1(num_datasets, n, 4.5*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 4.5  | Size/Power: ", res, "\n")
res <- CB_Min_mix1(num_datasets, n, 4.75*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 4.75 | Size/Power: ", res, "\n")
res <- CB_Min_mix1(num_datasets, n, 5*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 5    | Size/Power: ", res, "\n")


# ==============================================================================
cat("\n--- Distribution: Mix-2 ---\n")
# ==============================================================================
res <- CB_Min_mix2(num_datasets, n, 0*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 0    | Size/Power: ", res, "\n")
res <- CB_Min_mix2(num_datasets, n, 1*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 1    | Size/Power: ", res, "\n")
res <- CB_Min_mix2(num_datasets, n, 1.25*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 1.25 | Size/Power: ", res, "\n")
res <- CB_Min_mix2(num_datasets, n, 1.5*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 1.5  | Size/Power: ", res, "\n")
res <- CB_Min_mix2(num_datasets, n, 1.75*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 1.75 | Size/Power: ", res, "\n")
res <- CB_Min_mix2(num_datasets, n, 2*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 2    | Size/Power: ", res, "\n")
res <- CB_Min_mix2(num_datasets, n, 2.25*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 2.25 | Size/Power: ", res, "\n")
res <- CB_Min_mix2(num_datasets, n, 2.5*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 2.5  | Size/Power: ", res, "\n")
res <- CB_Min_mix2(num_datasets, n, 2.75*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 2.75 | Size/Power: ", res, "\n")
res <- CB_Min_mix2(num_datasets, n, 3*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 3    | Size/Power: ", res, "\n")
res <- CB_Min_mix2(num_datasets, n, 3.25*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 3.25 | Size/Power: ", res, "\n")
res <- CB_Min_mix2(num_datasets, n, 3.5*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 3.5  | Size/Power: ", res, "\n")
res <- CB_Min_mix2(num_datasets, n, 3.75*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 3.75 | Size/Power: ", res, "\n")
res <- CB_Min_mix2(num_datasets, n, 4*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 4    | Size/Power: ", res, "\n")
res <- CB_Min_mix2(num_datasets, n, 4.25*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 4.25 | Size/Power: ", res, "\n")
res <- CB_Min_mix2(num_datasets, n, 4.5*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 4.5  | Size/Power: ", res, "\n")
res <- CB_Min_mix2(num_datasets, n, 4.75*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 4.75 | Size/Power: ", res, "\n")
res <- CB_Min_mix2(num_datasets, n, 5*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 5    | Size/Power: ", res, "\n")


# ==============================================================================
cat("\n--- Distribution: Mix-3 ---\n")
# ==============================================================================
res <- CB_Min_mix3(num_datasets, n, 0*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 0    | Size/Power: ", res, "\n")
res <- CB_Min_mix3(num_datasets, n, 1*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 1    | Size/Power: ", res, "\n")
res <- CB_Min_mix3(num_datasets, n, 1.25*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 1.25 | Size/Power: ", res, "\n")
res <- CB_Min_mix3(num_datasets, n, 1.5*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 1.5  | Size/Power: ", res, "\n")
res <- CB_Min_mix3(num_datasets, n, 1.75*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 1.75 | Size/Power: ", res, "\n")
res <- CB_Min_mix3(num_datasets, n, 2*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 2    | Size/Power: ", res, "\n")
res <- CB_Min_mix3(num_datasets, n, 2.25*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 2.25 | Size/Power: ", res, "\n")
res <- CB_Min_mix3(num_datasets, n, 2.5*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 2.5  | Size/Power: ", res, "\n")
res <- CB_Min_mix3(num_datasets, n, 2.75*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 2.75 | Size/Power: ", res, "\n")
res <- CB_Min_mix3(num_datasets, n, 3*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 3    | Size/Power: ", res, "\n")
res <- CB_Min_mix3(num_datasets, n, 3.25*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 3.25 | Size/Power: ", res, "\n")
res <- CB_Min_mix3(num_datasets, n, 3.5*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 3.5  | Size/Power: ", res, "\n")
res <- CB_Min_mix3(num_datasets, n, 3.75*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 3.75 | Size/Power: ", res, "\n")
res <- CB_Min_mix3(num_datasets, n, 4*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 4    | Size/Power: ", res, "\n")
res <- CB_Min_mix3(num_datasets, n, 4.25*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 4.25 | Size/Power: ", res, "\n")
res <- CB_Min_mix3(num_datasets, n, 4.5*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 4.5  | Size/Power: ", res, "\n")
res <- CB_Min_mix3(num_datasets, n, 4.75*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 4.75 | Size/Power: ", res, "\n")
res <- CB_Min_mix3(num_datasets, n, 5*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 5    | Size/Power: ", res, "\n")


# ==============================================================================
cat("\n--- Distribution: Mix-4 ---\n")
# ==============================================================================
res <- CB_Min_mix4(num_datasets, n, 0*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 0    | Size/Power: ", res, "\n")
res <- CB_Min_mix4(num_datasets, n, 1*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 1    | Size/Power: ", res, "\n")
res <- CB_Min_mix4(num_datasets, n, 1.25*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 1.25 | Size/Power: ", res, "\n")
res <- CB_Min_mix4(num_datasets, n, 1.5*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 1.5  | Size/Power: ", res, "\n")
res <- CB_Min_mix4(num_datasets, n, 1.75*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 1.75 | Size/Power: ", res, "\n")
res <- CB_Min_mix4(num_datasets, n, 2*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 2    | Size/Power: ", res, "\n")
res <- CB_Min_mix4(num_datasets, n, 2.25*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 2.25 | Size/Power: ", res, "\n")
res <- CB_Min_mix4(num_datasets, n, 2.5*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 2.5  | Size/Power: ", res, "\n")
res <- CB_Min_mix4(num_datasets, n, 2.75*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 2.75 | Size/Power: ", res, "\n")
res <- CB_Min_mix4(num_datasets, n, 3*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 3    | Size/Power: ", res, "\n")
res <- CB_Min_mix4(num_datasets, n, 3.25*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 3.25 | Size/Power: ", res, "\n")
res <- CB_Min_mix4(num_datasets, n, 3.5*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 3.5  | Size/Power: ", res, "\n")
res <- CB_Min_mix4(num_datasets, n, 3.75*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 3.75 | Size/Power: ", res, "\n")
res <- CB_Min_mix4(num_datasets, n, 4*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 4    | Size/Power: ", res, "\n")
res <- CB_Min_mix4(num_datasets, n, 4.25*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 4.25 | Size/Power: ", res, "\n")
res <- CB_Min_mix4(num_datasets, n, 4.5*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 4.5  | Size/Power: ", res, "\n")
res <- CB_Min_mix4(num_datasets, n, 4.75*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 4.75 | Size/Power: ", res, "\n")
res <- CB_Min_mix4(num_datasets, n, 5*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 5    | Size/Power: ", res, "\n")


# ==============================================================================
cat("\n--- Distribution: Laplace ---\n")
# ==============================================================================
res <- CB_Min_Lap(num_datasets, n, 0*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 0    | Size/Power: ", res, "\n")
res <- CB_Min_Lap(num_datasets, n, 1*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 1    | Size/Power: ", res, "\n")
res <- CB_Min_Lap(num_datasets, n, 1.25*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 1.25 | Size/Power: ", res, "\n")
res <- CB_Min_Lap(num_datasets, n, 1.5*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 1.5  | Size/Power: ", res, "\n")
res <- CB_Min_Lap(num_datasets, n, 1.75*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 1.75 | Size/Power: ", res, "\n")
res <- CB_Min_Lap(num_datasets, n, 2*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 2    | Size/Power: ", res, "\n")
res <- CB_Min_Lap(num_datasets, n, 2.25*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 2.25 | Size/Power: ", res, "\n")
res <- CB_Min_Lap(num_datasets, n, 2.5*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 2.5  | Size/Power: ", res, "\n")
res <- CB_Min_Lap(num_datasets, n, 2.75*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 2.75 | Size/Power: ", res, "\n")
res <- CB_Min_Lap(num_datasets, n, 3*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 3    | Size/Power: ", res, "\n")
res <- CB_Min_Lap(num_datasets, n, 3.25*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 3.25 | Size/Power: ", res, "\n")
res <- CB_Min_Lap(num_datasets, n, 3.5*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 3.5  | Size/Power: ", res, "\n")
res <- CB_Min_Lap(num_datasets, n, 3.75*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 3.75 | Size/Power: ", res, "\n")
res <- CB_Min_Lap(num_datasets, n, 4*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 4    | Size/Power: ", res, "\n")
res <- CB_Min_Lap(num_datasets, n, 4.25*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 4.25 | Size/Power: ", res, "\n")
res <- CB_Min_Lap(num_datasets, n, 4.5*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 4.5  | Size/Power: ", res, "\n")
res <- CB_Min_Lap(num_datasets, n, 4.75*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 4.75 | Size/Power: ", res, "\n")
res <- CB_Min_Lap(num_datasets, n, 5*mean, sd, k0=2)
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
res <- CB_Min(num_datasets, n, 0*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 0    | Size/Power: ", res, "\n")
res <- CB_Min(num_datasets, n, 1*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 1    | Size/Power: ", res, "\n")
res <- CB_Min(num_datasets, n, 1.25*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 1.25 | Size/Power: ", res, "\n")
res <- CB_Min(num_datasets, n, 1.5*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 1.5  | Size/Power: ", res, "\n")
res <- CB_Min(num_datasets, n, 1.75*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 1.75 | Size/Power: ", res, "\n")
res <- CB_Min(num_datasets, n, 2*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 2    | Size/Power: ", res, "\n")
res <- CB_Min(num_datasets, n, 2.25*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 2.25 | Size/Power: ", res, "\n")
res <- CB_Min(num_datasets, n, 2.5*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 2.5  | Size/Power: ", res, "\n")
res <- CB_Min(num_datasets, n, 2.75*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 2.75 | Size/Power: ", res, "\n")
res <- CB_Min(num_datasets, n, 3*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 3    | Size/Power: ", res, "\n")
res <- CB_Min(num_datasets, n, 3.25*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 3.25 | Size/Power: ", res, "\n")
res <- CB_Min(num_datasets, n, 3.5*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 3.5  | Size/Power: ", res, "\n")
res <- CB_Min(num_datasets, n, 3.75*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 3.75 | Size/Power: ", res, "\n")
res <- CB_Min(num_datasets, n, 4*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 4    | Size/Power: ", res, "\n")
res <- CB_Min(num_datasets, n, 4.25*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 4.25 | Size/Power: ", res, "\n")
res <- CB_Min(num_datasets, n, 4.5*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 4.5  | Size/Power: ", res, "\n")
res <- CB_Min(num_datasets, n, 4.75*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 4.75 | Size/Power: ", res, "\n")
res <- CB_Min(num_datasets, n, 5*mean, sd, k0=2)
cat("Normal (N) | Multiplier (c): 5    | Size/Power: ", res, "\n")


# ==============================================================================
cat("\n--- Distribution: Mix-1 ---\n")
# ==============================================================================
res <- CB_Min_mix1(num_datasets, n, 0*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 0    | Size/Power: ", res, "\n")
res <- CB_Min_mix1(num_datasets, n, 1*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 1    | Size/Power: ", res, "\n")
res <- CB_Min_mix1(num_datasets, n, 1.25*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 1.25 | Size/Power: ", res, "\n")
res <- CB_Min_mix1(num_datasets, n, 1.5*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 1.5  | Size/Power: ", res, "\n")
res <- CB_Min_mix1(num_datasets, n, 1.75*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 1.75 | Size/Power: ", res, "\n")
res <- CB_Min_mix1(num_datasets, n, 2*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 2    | Size/Power: ", res, "\n")
res <- CB_Min_mix1(num_datasets, n, 2.25*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 2.25 | Size/Power: ", res, "\n")
res <- CB_Min_mix1(num_datasets, n, 2.5*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 2.5  | Size/Power: ", res, "\n")
res <- CB_Min_mix1(num_datasets, n, 2.75*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 2.75 | Size/Power: ", res, "\n")
res <- CB_Min_mix1(num_datasets, n, 3*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 3    | Size/Power: ", res, "\n")
res <- CB_Min_mix1(num_datasets, n, 3.25*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 3.25 | Size/Power: ", res, "\n")
res <- CB_Min_mix1(num_datasets, n, 3.5*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 3.5  | Size/Power: ", res, "\n")
res <- CB_Min_mix1(num_datasets, n, 3.75*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 3.75 | Size/Power: ", res, "\n")
res <- CB_Min_mix1(num_datasets, n, 4*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 4    | Size/Power: ", res, "\n")
res <- CB_Min_mix1(num_datasets, n, 4.25*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 4.25 | Size/Power: ", res, "\n")
res <- CB_Min_mix1(num_datasets, n, 4.5*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 4.5  | Size/Power: ", res, "\n")
res <- CB_Min_mix1(num_datasets, n, 4.75*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 4.75 | Size/Power: ", res, "\n")
res <- CB_Min_mix1(num_datasets, n, 5*mean, sd, k0=2)
cat("Mix-1      | Multiplier (c): 5    | Size/Power: ", res, "\n")


# ==============================================================================
cat("\n--- Distribution: Mix-2 ---\n")
# ==============================================================================
res <- CB_Min_mix2(num_datasets, n, 0*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 0    | Size/Power: ", res, "\n")
res <- CB_Min_mix2(num_datasets, n, 1*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 1    | Size/Power: ", res, "\n")
res <- CB_Min_mix2(num_datasets, n, 1.25*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 1.25 | Size/Power: ", res, "\n")
res <- CB_Min_mix2(num_datasets, n, 1.5*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 1.5  | Size/Power: ", res, "\n")
res <- CB_Min_mix2(num_datasets, n, 1.75*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 1.75 | Size/Power: ", res, "\n")
res <- CB_Min_mix2(num_datasets, n, 2*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 2    | Size/Power: ", res, "\n")
res <- CB_Min_mix2(num_datasets, n, 2.25*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 2.25 | Size/Power: ", res, "\n")
res <- CB_Min_mix2(num_datasets, n, 2.5*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 2.5  | Size/Power: ", res, "\n")
res <- CB_Min_mix2(num_datasets, n, 2.75*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 2.75 | Size/Power: ", res, "\n")
res <- CB_Min_mix2(num_datasets, n, 3*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 3    | Size/Power: ", res, "\n")
res <- CB_Min_mix2(num_datasets, n, 3.25*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 3.25 | Size/Power: ", res, "\n")
res <- CB_Min_mix2(num_datasets, n, 3.5*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 3.5  | Size/Power: ", res, "\n")
res <- CB_Min_mix2(num_datasets, n, 3.75*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 3.75 | Size/Power: ", res, "\n")
res <- CB_Min_mix2(num_datasets, n, 4*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 4    | Size/Power: ", res, "\n")
res <- CB_Min_mix2(num_datasets, n, 4.25*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 4.25 | Size/Power: ", res, "\n")
res <- CB_Min_mix2(num_datasets, n, 4.5*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 4.5  | Size/Power: ", res, "\n")
res <- CB_Min_mix2(num_datasets, n, 4.75*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 4.75 | Size/Power: ", res, "\n")
res <- CB_Min_mix2(num_datasets, n, 5*mean, sd, k0=2)
cat("Mix-2      | Multiplier (c): 5    | Size/Power: ", res, "\n")


# ==============================================================================
cat("\n--- Distribution: Mix-3 ---\n")
# ==============================================================================
res <- CB_Min_mix3(num_datasets, n, 0*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 0    | Size/Power: ", res, "\n")
res <- CB_Min_mix3(num_datasets, n, 1*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 1    | Size/Power: ", res, "\n")
res <- CB_Min_mix3(num_datasets, n, 1.25*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 1.25 | Size/Power: ", res, "\n")
res <- CB_Min_mix3(num_datasets, n, 1.5*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 1.5  | Size/Power: ", res, "\n")
res <- CB_Min_mix3(num_datasets, n, 1.75*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 1.75 | Size/Power: ", res, "\n")
res <- CB_Min_mix3(num_datasets, n, 2*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 2    | Size/Power: ", res, "\n")
res <- CB_Min_mix3(num_datasets, n, 2.25*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 2.25 | Size/Power: ", res, "\n")
res <- CB_Min_mix3(num_datasets, n, 2.5*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 2.5  | Size/Power: ", res, "\n")
res <- CB_Min_mix3(num_datasets, n, 2.75*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 2.75 | Size/Power: ", res, "\n")
res <- CB_Min_mix3(num_datasets, n, 3*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 3    | Size/Power: ", res, "\n")
res <- CB_Min_mix3(num_datasets, n, 3.25*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 3.25 | Size/Power: ", res, "\n")
res <- CB_Min_mix3(num_datasets, n, 3.5*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 3.5  | Size/Power: ", res, "\n")
res <- CB_Min_mix3(num_datasets, n, 3.75*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 3.75 | Size/Power: ", res, "\n")
res <- CB_Min_mix3(num_datasets, n, 4*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 4    | Size/Power: ", res, "\n")
res <- CB_Min_mix3(num_datasets, n, 4.25*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 4.25 | Size/Power: ", res, "\n")
res <- CB_Min_mix3(num_datasets, n, 4.5*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 4.5  | Size/Power: ", res, "\n")
res <- CB_Min_mix3(num_datasets, n, 4.75*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 4.75 | Size/Power: ", res, "\n")
res <- CB_Min_mix3(num_datasets, n, 5*mean, sd, k0=2)
cat("Mix-3      | Multiplier (c): 5    | Size/Power: ", res, "\n")


# ==============================================================================
cat("\n--- Distribution: Mix-4 ---\n")
# ==============================================================================
res <- CB_Min_mix4(num_datasets, n, 0*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 0    | Size/Power: ", res, "\n")
res <- CB_Min_mix4(num_datasets, n, 1*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 1    | Size/Power: ", res, "\n")
res <- CB_Min_mix4(num_datasets, n, 1.25*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 1.25 | Size/Power: ", res, "\n")
res <- CB_Min_mix4(num_datasets, n, 1.5*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 1.5  | Size/Power: ", res, "\n")
res <- CB_Min_mix4(num_datasets, n, 1.75*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 1.75 | Size/Power: ", res, "\n")
res <- CB_Min_mix4(num_datasets, n, 2*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 2    | Size/Power: ", res, "\n")
res <- CB_Min_mix4(num_datasets, n, 2.25*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 2.25 | Size/Power: ", res, "\n")
res <- CB_Min_mix4(num_datasets, n, 2.5*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 2.5  | Size/Power: ", res, "\n")
res <- CB_Min_mix4(num_datasets, n, 2.75*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 2.75 | Size/Power: ", res, "\n")
res <- CB_Min_mix4(num_datasets, n, 3*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 3    | Size/Power: ", res, "\n")
res <- CB_Min_mix4(num_datasets, n, 3.25*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 3.25 | Size/Power: ", res, "\n")
res <- CB_Min_mix4(num_datasets, n, 3.5*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 3.5  | Size/Power: ", res, "\n")
res <- CB_Min_mix4(num_datasets, n, 3.75*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 3.75 | Size/Power: ", res, "\n")
res <- CB_Min_mix4(num_datasets, n, 4*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 4    | Size/Power: ", res, "\n")
res <- CB_Min_mix4(num_datasets, n, 4.25*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 4.25 | Size/Power: ", res, "\n")
res <- CB_Min_mix4(num_datasets, n, 4.5*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 4.5  | Size/Power: ", res, "\n")
res <- CB_Min_mix4(num_datasets, n, 4.75*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 4.75 | Size/Power: ", res, "\n")
res <- CB_Min_mix4(num_datasets, n, 5*mean, sd, k0=2)
cat("Mix-4      | Multiplier (c): 5    | Size/Power: ", res, "\n")


# ==============================================================================
cat("\n--- Distribution: Laplace ---\n")
# ==============================================================================
res <- CB_Min_Lap(num_datasets, n, 0*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 0    | Size/Power: ", res, "\n")
res <- CB_Min_Lap(num_datasets, n, 1*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 1    | Size/Power: ", res, "\n")
res <- CB_Min_Lap(num_datasets, n, 1.25*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 1.25 | Size/Power: ", res, "\n")
res <- CB_Min_Lap(num_datasets, n, 1.5*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 1.5  | Size/Power: ", res, "\n")
res <- CB_Min_Lap(num_datasets, n, 1.75*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 1.75 | Size/Power: ", res, "\n")
res <- CB_Min_Lap(num_datasets, n, 2*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 2    | Size/Power: ", res, "\n")
res <- CB_Min_Lap(num_datasets, n, 2.25*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 2.25 | Size/Power: ", res, "\n")
res <- CB_Min_Lap(num_datasets, n, 2.5*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 2.5  | Size/Power: ", res, "\n")
res <- CB_Min_Lap(num_datasets, n, 2.75*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 2.75 | Size/Power: ", res, "\n")
res <- CB_Min_Lap(num_datasets, n, 3*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 3    | Size/Power: ", res, "\n")
res <- CB_Min_Lap(num_datasets, n, 3.25*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 3.25 | Size/Power: ", res, "\n")
res <- CB_Min_Lap(num_datasets, n, 3.5*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 3.5  | Size/Power: ", res, "\n")
res <- CB_Min_Lap(num_datasets, n, 3.75*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 3.75 | Size/Power: ", res, "\n")
res <- CB_Min_Lap(num_datasets, n, 4*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 4    | Size/Power: ", res, "\n")
res <- CB_Min_Lap(num_datasets, n, 4.25*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 4.25 | Size/Power: ", res, "\n")
res <- CB_Min_Lap(num_datasets, n, 4.5*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 4.5  | Size/Power: ", res, "\n")
res <- CB_Min_Lap(num_datasets, n, 4.75*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 4.75 | Size/Power: ", res, "\n")
res <- CB_Min_Lap(num_datasets, n, 5*mean, sd, k0=2)
cat("Laplace    | Multiplier (c): 5    | Size/Power: ", res, "\n")





cat("\n==========================================================================\n")
cat("TABLE S7.22 REPLICATION COMPLETE\n")
cat("==========================================================================\n")





