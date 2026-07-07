# ==============================================================================
# REPLICATION SCRIPT: Empirical Size and Power – Comparison with Benchmark Test Tables 
#(Table 8 of the Manuscript and Table S6.13 of Supplementary Material)
#
# NOTE ON COMPUTATIONAL TIME:
# The results reported in the manuscript were obtained on a 64-core
# high-performance computing server using M = 5000 and Q = 5000.
# Under these settings, computing a single empirical size value requires
# approximately 30 minutes for the CBLRT and approximately 9 minutes for each
# of the other proposed tests.
#
# *** WARNING ***
# Reproducing the manuscript results with M = 5000 and Q = 5000 on
# a standard desktop or laptop (e.g., 4–16 CPU cores) can result in extremely
# long execution times, potentially requiring several days to complete.
#
# To facilitate quick reproducibility, this script automatically detects the
# number of available CPU cores and, by default, performs a rapid test run with
# 100 bootstrap repetitions and 100 simulation repetitions. To exactly
# reproduce the empirical size results reported in the manuscript, set both
# M and Q to 5000.
# ==============================================================================

# Check each package one by one. If it's missing, install it.
if (!require("Iso")) install.packages("Iso")
if (!require("parallel")) install.packages("parallel")
if (!require("doParallel")) install.packages("doParallel")
if (!require("foreach")) install.packages("foreach")

# Load them
library(Iso)
library(parallel)
library(doParallel)
library(foreach)


# Number of bootstrap repetitions
M <- 100  # Set to 5000 to exactly reproduce the empirical size results reported in the manuscript.

# Number of simulation repetitions
Q <- 100   # Set to 5000 to exactly reproduce the empirical size results reported in the manuscript.

# Number of CPU cores to use
CORE <- max(1, parallel::detectCores() - 1)  # The manuscript used CORE = 64.


##MAIN CODE
CB_LRT <- function(num_datasets, n, mean, sd, k0, num_repetitions = M, outnum_repetitions = Q, num_cores = CORE) {
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
  cl <- makeCluster(num_cores)
  on.exit({
    stopCluster(cl)
    registerDoSEQ()
    closeAllConnections()
  }, add = TRUE)
  registerDoParallel(cl)
  clusterExport(cl, c("pava", "CB_H0_new", "CB_H1_new", "CB_MLE"), envir = environment())
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
CB_MaxMin <- function(num_datasets, n, mean, sd, k0, num_samples = M, outnum_repetitions = Q, num_cores = CORE) {
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
          T[p] <- min(sapply((k0+1):num_datasets, function(i) {
            
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
        CB_values[p] <- min(sapply((k0+1):num_datasets, function(i) {
          
          (mean(sample_data[[i]]) - mean(sample_data[[p]])) /
            
            sqrt(
              
              (var(sample_data[[i]]) / length(sample_data[[i]])) +
                
                (var(sample_data[[p]]) / length(sample_data[[p]]))
              
            )
          
        }))
      }
      CB <- max(CB_values)
      quantile_values <- quantile(sort_T_max, probs = 0.95)
      alpha_CBmax <- ifelse(quantile_values < CB, 1, 0)
      return(alpha_CBmax)
    })
  }
  return((results))
}
CB_MinMax <- function(num_datasets, n, mean, sd, k0, num_samples = M, outnum_repetitions = Q, num_cores = CORE) {
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
          T[p] <- max(sapply((k0+1):num_datasets, function(i) {
            
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
        CB_values[p] <- max(sapply((k0+1):num_datasets, function(i) {
          
          (mean(sample_data[[i]]) - mean(sample_data[[p]])) /
            
            sqrt(
              
              (var(sample_data[[i]]) / length(sample_data[[i]])) +
                
                (var(sample_data[[p]]) / length(sample_data[[p]]))
              
            )
          
        }))
      }
      CB <- min(CB_values)
      quantile_values <- quantile(sort_T_min, probs = 0.95)
      alpha_CBmax <- ifelse(quantile_values < CB, 1, 0)
      return(alpha_CBmax)
    })
  }
  return((results))
}


# SIMULATION RESULTS

cat("\n=======================================================\n")
cat("STARTING SIMULATION of Empirical Size Comparison with Benchmark Test Table (Table 8 of the Manuscript)\n")
cat("=======================================================\n\n")

NUM <- 6
MEAN <- c(0, 0, 0, 0, 0, 0)
SD1 <- c(sqrt(1), sqrt(1), sqrt(1), sqrt(1), sqrt(1), sqrt(1))
SD2 <- c(sqrt(0.1), sqrt(0.1), sqrt(0.1), sqrt(0.5), sqrt(0.5), sqrt(0.5))
SD3 <- c(sqrt(1), sqrt(1), sqrt(1), sqrt(0.5), sqrt(0.5), sqrt(0.5))
SD4 <- c(sqrt(0.1), sqrt(0.2), sqrt(0.3), sqrt(0.4), sqrt(0.5), sqrt(1))
SD5 <- c(sqrt(0.3), sqrt(0.9), sqrt(0.4), sqrt(0.7), sqrt(0.5), sqrt(1))
SD6 <- c(sqrt(0.01), sqrt(0.1), sqrt(0.1), sqrt(0.1), sqrt(0.1), sqrt(1))

# ---------------------------------------------------------
cat("\n--- Generating LRT column of comparison with existing test table ---\n")
# ---------------------------------------------------------

res <- CB_LRT(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical size of the CBLRT for parameter configuration N1_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical size of the CBLRT for parameter configuration N2_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical size of the CBLRT for parameter configuration N3_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical size of the CBLRT for parameter configuration N4_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical size of the CBLRT for parameter configuration N5_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")


res <- CB_LRT(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical size of the CBLRT for parameter configuration N1_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical size of the CBLRT for parameter configuration N2_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical size of the CBLRT for parameter configuration N3_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical size of the CBLRT for parameter configuration N4_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical size of the CBLRT for parameter configuration N5_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")


res <- CB_LRT(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical size of the CBLRT for parameter configuration N1_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical size of the CBLRT for parameter configuration N2_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical size of the CBLRT for parameter configuration N3_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical size of the CBLRT for parameter configuration N4_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical size of the CBLRT for parameter configuration N5_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")


res <- CB_LRT(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical size of the CBLRT for parameter configuration N1_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical size of the CBLRT for parameter configuration N2_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical size of the CBLRT for parameter configuration N3_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical size of the CBLRT for parameter configuration N4_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical size of the CBLRT for parameter configuration N5_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")


res <- CB_LRT(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical size of the CBLRT for parameter configuration N1_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical size of the CBLRT for parameter configuration N2_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical size of the CBLRT for parameter configuration N3_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical size of the CBLRT for parameter configuration N4_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical size of the CBLRT for parameter configuration N5_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")


res <- CB_LRT(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical size of the CBLRT for parameter configuration N1_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical size of the CBLRT for parameter configuration N2_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical size of the CBLRT for parameter configuration N3_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical size of the CBLRT for parameter configuration N4_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical size of the CBLRT for parameter configuration N5_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")



# ---------------------------------------------------------
cat("\n--- Generating Max column of comparison with existing test table ---\n")
# ---------------------------------------------------------

res <- CB_Max(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical size of the CBMax for parameter configuration N1_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical size of the CBMax for parameter configuration N2_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical size of the CBMax for parameter configuration N3_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical size of the CBMax for parameter configuration N4_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical size of the CBMax for parameter configuration N5_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")


res <- CB_Max(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical size of the CBMax for parameter configuration N1_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical size of the CBMax for parameter configuration N2_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical size of the CBMax for parameter configuration N3_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical size of the CBMax for parameter configuration N4_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical size of the CBMax for parameter configuration N5_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")


res <- CB_Max(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical size of the CBMax for parameter configuration N1_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical size of the CBMax for parameter configuration N2_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical size of the CBMax for parameter configuration N3_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical size of the CBMax for parameter configuration N4_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical size of the CBMax for parameter configuration N5_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")


res <- CB_Max(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical size of the CBMax for parameter configuration N1_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical size of the CBMax for parameter configuration N2_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical size of the CBMax for parameter configuration N3_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical size of the CBMax for parameter configuration N4_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical size of the CBMax for parameter configuration N5_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")


res <- CB_Max(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical size of the CBMax for parameter configuration N1_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical size of the CBMax for parameter configuration N2_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical size of the CBMax for parameter configuration N3_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical size of the CBMax for parameter configuration N4_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical size of the CBMax for parameter configuration N5_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")


res <- CB_Max(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical size of the CBMax for parameter configuration N1_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical size of the CBMax for parameter configuration N2_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical size of the CBMax for parameter configuration N3_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical size of the CBMax for parameter configuration N4_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical size of the CBMax for parameter configuration N5_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")


# ---------------------------------------------------------
cat("\n--- Generating Min column of comparison with existing test table ---\n")
# ---------------------------------------------------------

res <- CB_Min(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical size of the CBMin for parameter configuration N1_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical size of the CBMin for parameter configuration N2_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical size of the CBMin for parameter configuration N3_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical size of the CBMin for parameter configuration N4_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical size of the CBMin for parameter configuration N5_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")


res <- CB_Min(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical size of the CBMin for parameter configuration N1_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical size of the CBMin for parameter configuration N2_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical size of the CBMin for parameter configuration N3_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical size of the CBMin for parameter configuration N4_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical size of the CBMin for parameter configuration N5_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")


res <- CB_Min(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical size of the CBMin for parameter configuration N1_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical size of the CBMin for parameter configuration N2_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical size of the CBMin for parameter configuration N3_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical size of the CBMin for parameter configuration N4_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical size of the CBMin for parameter configuration N5_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")


res <- CB_Min(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical size of the CBMin for parameter configuration N1_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical size of the CBMin for parameter configuration N2_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical size of the CBMin for parameter configuration N3_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical size of the CBMin for parameter configuration N4_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical size of the CBMin for parameter configuration N5_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")


res <- CB_Min(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical size of the CBMin for parameter configuration N1_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical size of the CBMin for parameter configuration N2_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical size of the CBMin for parameter configuration N3_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical size of the CBMin for parameter configuration N4_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical size of the CBMin for parameter configuration N5_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")


res <- CB_Min(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical size of the CBMin for parameter configuration N1_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical size of the CBMin for parameter configuration N2_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical size of the CBMin for parameter configuration N3_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical size of the CBMin for parameter configuration N4_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical size of the CBMin for parameter configuration N5_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")


# ---------------------------------------------------------
cat("\n--- Generating MaxMin column of comparison with existing test table ---\n")
# ---------------------------------------------------------


res <- CB_MaxMin(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical size of the CBMaxMin for parameter configuration N1_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical size of the CBMaxMin for parameter configuration N2_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical size of the CBMaxMin for parameter configuration N3_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical size of the CBMaxMin for parameter configuration N4_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical size of the CBMaxMin for parameter configuration N5_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")


res <- CB_MaxMin(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical size of the CBMaxMin for parameter configuration N1_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical size of the CBMaxMin for parameter configuration N2_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical size of the CBMaxMin for parameter configuration N3_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical size of the CBMaxMin for parameter configuration N4_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical size of the CBMaxMin for parameter configuration N5_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")


res <- CB_MaxMin(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical size of the CBMaxMin for parameter configuration N1_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical size of the CBMaxMin for parameter configuration N2_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical size of the CBMaxMin for parameter configuration N3_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical size of the CBMaxMin for parameter configuration N4_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical size of the CBMaxMin for parameter configuration N5_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")


res <- CB_MaxMin(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical size of the CBMaxMin for parameter configuration N1_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical size of the CBMaxMin for parameter configuration N2_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical size of the CBMaxMin for parameter configuration N3_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical size of the CBMaxMin for parameter configuration N4_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical size of the CBMaxMin for parameter configuration N5_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")


res <- CB_MaxMin(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical size of the CBMaxMin for parameter configuration N1_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical size of the CBMaxMin for parameter configuration N2_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical size of the CBMaxMin for parameter configuration N3_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical size of the CBMaxMin for parameter configuration N4_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical size of the CBMaxMin for parameter configuration N5_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")


res <- CB_MaxMin(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical size of the CBMaxMin for parameter configuration N1_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical size of the CBMaxMin for parameter configuration N2_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical size of the CBMaxMin for parameter configuration N3_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical size of the CBMaxMin for parameter configuration N4_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical size of the CBMaxMin for parameter configuration N5_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")



# ---------------------------------------------------------
cat("\n--- Generating MinMax column of comparison with existing test table ---\n")
# ---------------------------------------------------------


res <- CB_MinMax(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical size of the CBMinMax for parameter configuration N1_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical size of the CBMinMax for parameter configuration N2_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical size of the CBMinMax for parameter configuration N3_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical size of the CBMinMax for parameter configuration N4_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical size of the CBMinMax for parameter configuration N5_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")


res <- CB_MinMax(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical size of the CBMinMax for parameter configuration N1_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical size of the CBMinMax for parameter configuration N2_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical size of the CBMinMax for parameter configuration N3_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical size of the CBMinMax for parameter configuration N4_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical size of the CBMinMax for parameter configuration N5_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")


res <- CB_MinMax(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical size of the CBMinMax for parameter configuration N1_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical size of the CBMinMax for parameter configuration N2_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical size of the CBMinMax for parameter configuration N3_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical size of the CBMinMax for parameter configuration N4_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical size of the CBMinMax for parameter configuration N5_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")


res <- CB_MinMax(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical size of the CBMinMax for parameter configuration N1_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical size of the CBMinMax for parameter configuration N2_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical size of the CBMinMax for parameter configuration N3_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical size of the CBMinMax for parameter configuration N4_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical size of the CBMinMax for parameter configuration N5_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")


res <- CB_MinMax(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical size of the CBMinMax for parameter configuration N1_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical size of the CBMinMax for parameter configuration N2_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical size of the CBMinMax for parameter configuration N3_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical size of the CBMinMax for parameter configuration N4_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical size of the CBMinMax for parameter configuration N5_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")


res <- CB_MinMax(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical size of the CBMinMax for parameter configuration N1_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical size of the CBMinMax for parameter configuration N2_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical size of the CBMinMax for parameter configuration N3_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical size of the CBMinMax for parameter configuration N4_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical size of the CBMinMax for parameter configuration N5_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")



#################################################################################################################################################################



cat("\n=======================================================\n")
cat("STARTING SIMULATION of Empirical Ordinary Power Comparison with Benchmark Test Table (Table S6.13 of Supplementary Material)\n")
cat("=======================================================\n\n")

NUM <- 6
MEAN <- c(0, 0, -0.2, 0.2, 0.4, 0.8)
SD1 <- c(sqrt(1), sqrt(1), sqrt(1), sqrt(1), sqrt(1), sqrt(1))
SD2 <- c(sqrt(0.1), sqrt(0.1), sqrt(0.1), sqrt(0.5), sqrt(0.5), sqrt(0.5))
SD3 <- c(sqrt(1), sqrt(1), sqrt(1), sqrt(0.5), sqrt(0.5), sqrt(0.5))
SD4 <- c(sqrt(0.1), sqrt(0.2), sqrt(0.3), sqrt(0.4), sqrt(0.5), sqrt(1))
SD5 <- c(sqrt(0.3), sqrt(0.9), sqrt(0.4), sqrt(0.7), sqrt(0.5), sqrt(1))
SD6 <- c(sqrt(0.01), sqrt(0.1), sqrt(0.1), sqrt(0.1), sqrt(0.1), sqrt(1))

# ---------------------------------------------------------
cat("\n--- Generating LRT column of comparison with existing test table ---\n")
# ---------------------------------------------------------

res <- CB_LRT(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical power of the CBLRT for parameter configuration N1_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical power of the CBLRT for parameter configuration N2_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical power of the CBLRT for parameter configuration N3_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical power of the CBLRT for parameter configuration N4_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical power of the CBLRT for parameter configuration N5_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")


res <- CB_LRT(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical power of the CBLRT for parameter configuration N1_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical power of the CBLRT for parameter configuration N2_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical power of the CBLRT for parameter configuration N3_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical power of the CBLRT for parameter configuration N4_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical power of the CBLRT for parameter configuration N5_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")


res <- CB_LRT(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical power of the CBLRT for parameter configuration N1_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical power of the CBLRT for parameter configuration N2_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical power of the CBLRT for parameter configuration N3_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical power of the CBLRT for parameter configuration N4_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical power of the CBLRT for parameter configuration N5_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")


res <- CB_LRT(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical power of the CBLRT for parameter configuration N1_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical power of the CBLRT for parameter configuration N2_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical power of the CBLRT for parameter configuration N3_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical power of the CBLRT for parameter configuration N4_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical power of the CBLRT for parameter configuration N5_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")


res <- CB_LRT(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical power of the CBLRT for parameter configuration N1_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical power of the CBLRT for parameter configuration N2_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical power of the CBLRT for parameter configuration N3_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical power of the CBLRT for parameter configuration N4_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical power of the CBLRT for parameter configuration N5_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")


res <- CB_LRT(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical power of the CBLRT for parameter configuration N1_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical power of the CBLRT for parameter configuration N2_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical power of the CBLRT for parameter configuration N3_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical power of the CBLRT for parameter configuration N4_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")
res <- CB_LRT(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical power of the CBLRT for parameter configuration N5_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")



# ---------------------------------------------------------
cat("\n--- Generating Max column of comparison with existing test table ---\n")
# ---------------------------------------------------------

res <- CB_Max(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical power of the CBMax for parameter configuration N1_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical power of the CBMax for parameter configuration N2_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical power of the CBMax for parameter configuration N3_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical power of the CBMax for parameter configuration N4_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical power of the CBMax for parameter configuration N5_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")


res <- CB_Max(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical power of the CBMax for parameter configuration N1_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical power of the CBMax for parameter configuration N2_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical power of the CBMax for parameter configuration N3_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical power of the CBMax for parameter configuration N4_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical power of the CBMax for parameter configuration N5_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")


res <- CB_Max(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical power of the CBMax for parameter configuration N1_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical power of the CBMax for parameter configuration N2_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical power of the CBMax for parameter configuration N3_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical power of the CBMax for parameter configuration N4_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical power of the CBMax for parameter configuration N5_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")


res <- CB_Max(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical power of the CBMax for parameter configuration N1_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical power of the CBMax for parameter configuration N2_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical power of the CBMax for parameter configuration N3_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical power of the CBMax for parameter configuration N4_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical power of the CBMax for parameter configuration N5_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")


res <- CB_Max(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical power of the CBMax for parameter configuration N1_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical power of the CBMax for parameter configuration N2_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical power of the CBMax for parameter configuration N3_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical power of the CBMax for parameter configuration N4_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical power of the CBMax for parameter configuration N5_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")


res <- CB_Max(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical power of the CBMax for parameter configuration N1_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical power of the CBMax for parameter configuration N2_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical power of the CBMax for parameter configuration N3_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical power of the CBMax for parameter configuration N4_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")
res <- CB_Max(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical power of the CBMax for parameter configuration N5_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")


# ---------------------------------------------------------
cat("\n--- Generating Min column of comparison with existing test table ---\n")
# ---------------------------------------------------------

res <- CB_Min(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical power of the CBMin for parameter configuration N1_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical power of the CBMin for parameter configuration N2_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical power of the CBMin for parameter configuration N3_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical power of the CBMin for parameter configuration N4_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical power of the CBMin for parameter configuration N5_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")


res <- CB_Min(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical power of the CBMin for parameter configuration N1_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical power of the CBMin for parameter configuration N2_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical power of the CBMin for parameter configuration N3_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical power of the CBMin for parameter configuration N4_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical power of the CBMin for parameter configuration N5_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")


res <- CB_Min(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical power of the CBMin for parameter configuration N1_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical power of the CBMin for parameter configuration N2_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical power of the CBMin for parameter configuration N3_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical power of the CBMin for parameter configuration N4_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical power of the CBMin for parameter configuration N5_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")


res <- CB_Min(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical power of the CBMin for parameter configuration N1_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical power of the CBMin for parameter configuration N2_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical power of the CBMin for parameter configuration N3_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical power of the CBMin for parameter configuration N4_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical power of the CBMin for parameter configuration N5_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")


res <- CB_Min(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical power of the CBMin for parameter configuration N1_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical power of the CBMin for parameter configuration N2_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical power of the CBMin for parameter configuration N3_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical power of the CBMin for parameter configuration N4_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical power of the CBMin for parameter configuration N5_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")


res <- CB_Min(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical power of the CBMin for parameter configuration N1_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical power of the CBMin for parameter configuration N2_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical power of the CBMin for parameter configuration N3_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical power of the CBMin for parameter configuration N4_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")
res <- CB_Min(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical power of the CBMin for parameter configuration N5_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")


# ---------------------------------------------------------
cat("\n--- Generating MaxMin column of comparison with existing test table ---\n")
# ---------------------------------------------------------


res <- CB_MaxMin(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical power of the CBMaxMin for parameter configuration N1_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical power of the CBMaxMin for parameter configuration N2_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical power of the CBMaxMin for parameter configuration N3_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical power of the CBMaxMin for parameter configuration N4_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical power of the CBMaxMin for parameter configuration N5_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")


res <- CB_MaxMin(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical power of the CBMaxMin for parameter configuration N1_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical power of the CBMaxMin for parameter configuration N2_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical power of the CBMaxMin for parameter configuration N3_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical power of the CBMaxMin for parameter configuration N4_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical power of the CBMaxMin for parameter configuration N5_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")


res <- CB_MaxMin(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical power of the CBMaxMin for parameter configuration N1_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical power of the CBMaxMin for parameter configuration N2_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical power of the CBMaxMin for parameter configuration N3_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical power of the CBMaxMin for parameter configuration N4_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical power of the CBMaxMin for parameter configuration N5_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")


res <- CB_MaxMin(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical power of the CBMaxMin for parameter configuration N1_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical power of the CBMaxMin for parameter configuration N2_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical power of the CBMaxMin for parameter configuration N3_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical power of the CBMaxMin for parameter configuration N4_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical power of the CBMaxMin for parameter configuration N5_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")


res <- CB_MaxMin(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical power of the CBMaxMin for parameter configuration N1_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical power of the CBMaxMin for parameter configuration N2_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical power of the CBMaxMin for parameter configuration N3_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical power of the CBMaxMin for parameter configuration N4_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical power of the CBMaxMin for parameter configuration N5_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")


res <- CB_MaxMin(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical power of the CBMaxMin for parameter configuration N1_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical power of the CBMaxMin for parameter configuration N2_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical power of the CBMaxMin for parameter configuration N3_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical power of the CBMaxMin for parameter configuration N4_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")
res <- CB_MaxMin(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical power of the CBMaxMin for parameter configuration N5_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")



# ---------------------------------------------------------
cat("\n--- Generating MinMax column of comparison with existing test table ---\n")
# ---------------------------------------------------------


res <- CB_MinMax(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical power of the CBMinMax for parameter configuration N1_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical power of the CBMinMax for parameter configuration N2_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical power of the CBMinMax for parameter configuration N3_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical power of the CBMinMax for parameter configuration N4_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD1, k0 = 3)
cat("Empirical power of the CBMinMax for parameter configuration N5_star with variance configuration (1, 1, 1, 1, 1, 1):        ", res, "\n")


res <- CB_MinMax(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical power of the CBMinMax for parameter configuration N1_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical power of the CBMinMax for parameter configuration N2_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical power of the CBMinMax for parameter configuration N3_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical power of the CBMinMax for parameter configuration N4_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD2, k0 = 3)
cat("Empirical power of the CBMinMax for parameter configuration N5_star with variance configuration (0.1,0.1,0.1,0.5,0.5,0.5):        ", res, "\n")


res <- CB_MinMax(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical power of the CBMinMax for parameter configuration N1_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical power of the CBMinMax for parameter configuration N2_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical power of the CBMinMax for parameter configuration N3_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical power of the CBMinMax for parameter configuration N4_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD3, k0 = 3)
cat("Empirical power of the CBMinMax for parameter configuration N5_star with variance configuration (1,1,1,0.5,0.5,0.5):        ", res, "\n")


res <- CB_MinMax(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical power of the CBMinMax for parameter configuration N1_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical power of the CBMinMax for parameter configuration N2_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical power of the CBMinMax for parameter configuration N3_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical power of the CBMinMax for parameter configuration N4_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD4, k0 = 3)
cat("Empirical power of the CBMinMax for parameter configuration N5_star with variance configuration (0.1,0.2,0.3,0.4,0.5,1):        ", res, "\n")


res <- CB_MinMax(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical power of the CBMinMax for parameter configuration N1_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical power of the CBMinMax for parameter configuration N2_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical power of the CBMinMax for parameter configuration N3_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical power of the CBMinMax for parameter configuration N4_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD5, k0 = 3)
cat("Empirical power of the CBMinMax for parameter configuration N5_star with variance configuration (0.3,0.9,0.4,0.7,0.5,1):        ", res, "\n")


res <- CB_MinMax(num_datasets = NUM, n = c(5,5,5,5,5,5), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical power of the CBMinMax for parameter configuration N1_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(10,10,10,10,10,10), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical power of the CBMinMax for parameter configuration N2_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(3,3,4,5,6,6), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical power of the CBMinMax for parameter configuration N3_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(4,6,8,12,16,20), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical power of the CBMinMax for parameter configuration N4_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")
res <- CB_MinMax(num_datasets = NUM, n = c(35,35,35,35,35,35), mean = MEAN, sd = SD6, k0 = 3)
cat("Empirical power of the CBMinMax for parameter configuration N5_star with variance configuration (0.01,0.1,0.1,0.1,0.1,1):        ", res, "\n")









