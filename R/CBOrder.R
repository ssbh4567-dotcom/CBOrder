#' Likelihood Ratio Test for Complete Bipartite (CB) Ordered Alternatives
#'
#' Performs a likelihood ratio test (LRT) to evaluate the equality of means
#' across multiple groups under the Complete Bipartite (CB) ordered alternatives.
#' The test compares the null hypothesis of equal means for all groups
#' to the alternative that the control group means are less than or equal
#' to the treatment group means under the CB order restriction.
#'
#' @param sample_data A list of numeric vectors, where the first \code{k0} elements
#' represent control groups and the remaining elements represent treatment groups.
#' @param significance_level A numeric value between 0 and 1 specifying the
#' significance level for the test (e.g., \code{0.05}).
#' @param k0 An integer specifying the number of control groups.
#' @param n.boot Number of bootstrap replications to estimate the critical value
#' (default is \code{100000}).
#' @param seed Optional random seed for reproducibility.
#'
#' @return A character string summarizing the results, including the critical value,
#' CBLRT test statistic, estimated p-value, and the test decision.
#'
#' @details The likelihood ratio statistic is computed using constrained
#' maximum likelihood estimates under the null and CB ordered alternative hypotheses.
#' The critical value and p-value are estimated by a bootstrap procedure.
#' The test evaluates the hypotheses H_0: μ_i = μ_j versus H_1: μ_i ≤ μ_j (at least one strict inequality)
#' for i = 1, ..., k0 (control groups) and j = k0 + 1, ..., p (total treatment groups).
#'
#' @seealso Halder, Mondal, and Kumar (2025)
#' "Testing Against Complete Bipartite Ordered Alternatives in One-way ANOVA"
#' (forthcoming manuscript).
#'
#' @importFrom stats quantile rnorm var
#' @importFrom Iso pava
#' @export
#'
#' @author Subha Halder
#'
#' @examples
#' # Generate sample data
#' set.seed(456)
#' control1 <- rnorm(20, mean = 5)
#' control2 <- rnorm(20, mean = 5)
#' treatment1 <- rnorm(20, mean = 6)
#' treatment2 <- rnorm(20, mean = 7)
#'
#' data_list <- list(control1, control2, treatment1, treatment2)
#'
#' # Run CB LRT with fewer bootstrap replications for demonstration
#' CBLRT(data_list, 0.05, k0 = 2, n.boot = 1000)
#'
#' \donttest{
#' # Recommended: Use 100000 for research-level accuracy
#' CBLRT(data_list, 0.05, k0 = 2, n.boot = 100000)
#' }
CBLRT <- function(sample_data, significance_level, k0, n.boot = 100000, seed = NULL) {
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

      if (max(abs(new_mu1 - mu1)) <= 0.000001) {
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

      if (max(abs(new_mu0 - mu0)) <= 0.000001) {
        break  # Exit the loop if the difference is less than epsilon
      }

      w0 <- new_w0
      mu0 <- new_mu0
      var0 <- new_var0
    }

    return(var0)
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  sample_data <- lapply(sample_data, function(x) x[!is.na(x)])
  num_samples <- n.boot
  num_datasets <- length(sample_data)
  n <- sapply(sample_data, length)
  lambda_values_star <- numeric(num_samples)

  # Bootstrap under null hypothesis
  for (i in 1:num_samples) {
    bootstrap_samples <- lapply(sample_data, function(x)
      rnorm(n = length(x), mean = 0, sd = sqrt(var(x))))

    V_R_star <- CB_H1_new(bootstrap_samples, k0) / CB_H0_new(bootstrap_samples)
    weights <- sapply(1:num_datasets, function(j)
      V_R_star[j]^(length(bootstrap_samples[[j]]) / 2))
    lambda_values_star[i] <- prod(weights)
  }

  # Critical value and test statistic
  quantile_value <- quantile(lambda_values_star, probs = significance_level)
  V_R <- CB_H1_new(sample_data, k0) / CB_H0_new(sample_data)
  weights <- sapply(1:num_datasets, function(i)
    V_R[i]^(length(sample_data[[i]]) / 2))
  lambda <- prod(weights)

  # Bootstrap-based p-value
  p_value <- mean(lambda_values_star <= lambda)

  # Decision
  if (lambda < quantile_value) {
    result <- "Reject null hypothesis"
  } else {
    result <- "Do not reject null hypothesis"
  }

  # Output
  return(paste(
    "Critical value:", quantile_value,
    "; CBLRT Test statistic:", lambda,
    "; p-value:", p_value,
    "; Result:", result
  ))
}

#' Maximum Difference-based Test for Complete Bipartite (CB) Ordered Alternatives
#'
#' Performs a Maximum Difference-based (CBMax) test to evaluate the equality of group means
#' under the Complete Bipartite (CB) order restriction. This test compares the null
#' hypothesis of equal means across all groups against the alternative that the means
#' of the control groups are less than or equal to those of the treatment groups.
#'
#' @param sample_data A list of numeric vectors, where the first \code{k0} elements
#' correspond to the control groups and the remaining elements correspond to the
#' treatment groups.
#' @param significance_level A numeric value between 0 and 1 specifying the
#' significance level of the test (e.g., \code{0.05}).
#' @param k0 An integer specifying the number of control groups.
#' @param n.boot Number of bootstrap replications used to estimate the critical value
#' (default is \code{100000}).
#' @param seed Optional random seed for reproducibility.
#'
#' @return A character string summarizing the results, including the critical value,
#' the CBMax test statistic, the estimated p-value, and the test decision.
#'
#' @details The test statistic is computed as the maximum standardized difference
#' between each treatment group mean and each control group mean under the CB order restriction.
#' The critical value and p-value are estimated via bootstrap resampling from the null distribution.
#' The test evaluates the hypotheses H_0: μ_i = μ_j versus H_1: μ_i ≤ μ_j (at least one strict inequality)
#' for i = 1, ..., k0 (control groups) and j = k0 + 1, ..., p (total treatment groups).
#'
#' @importFrom stats quantile rnorm var
#'
#' @author Subha Halder
#' @export
#'
#' @examples
#' # Generate sample data
#' set.seed(456)
#' control1 <- rnorm(20, mean = 5)
#' control2 <- rnorm(20, mean = 5)
#' treatment1 <- rnorm(20, mean = 6)
#' treatment2 <- rnorm(20, mean = 7)
#' data_list <- list(control1, control2, treatment1, treatment2)
#'
#' # Run CBMax test with fewer bootstrap samples for demonstration
#' CBMax(sample_data = data_list, significance_level = 0.05, k0 = 2, n.boot = 1000)
#'
#' \donttest{
#' # Recommended: use 100000 bootstrap samples for higher precision
#' CBMax(sample_data = data_list, significance_level = 0.05, k0 = 2, n.boot = 100000)
#' }
CBMax <- function(sample_data, significance_level, k0, n.boot = 100000, seed = NULL) {
  if (!is.null(seed)) { set.seed(seed) }

  sample_data <- lapply(sample_data, function(x) x[!is.na(x)])
  num_samples <- n.boot
  num_datasets <- length(sample_data)
  n <- sapply(sample_data, length)
  group_sds <- sapply(sample_data, function(x) sqrt(var(x)))
  T_max <- numeric(num_samples)

  for (u in 1:num_samples) {
    bootstrap_samples <- vector("list", num_datasets)   # Corrected initialization
    for (j in 1:num_datasets) {
      bootstrap_samples[[j]] <- rnorm(n = n[j], mean = 0, sd = group_sds[j])
    }
    T <- numeric(k0)
    for (p in 1:k0) {
      T[p] <- max(sapply((k0+1):num_datasets, function(i) {
        (mean(bootstrap_samples[[i]]) - mean(bootstrap_samples[[p]])) /
          sqrt((var(bootstrap_samples[[i]]) / length(bootstrap_samples[[i]])) +
                 (var(bootstrap_samples[[p]]) / length(bootstrap_samples[[p]])))
      }))
    }
    T_max[u] <- max(T)
  }

  CB_values <- numeric(k0)   # Corrected initialization
  for (p in 1:k0) {
    CB_values[p] <- max(sapply((k0+1):num_datasets, function(i) {
      (mean(sample_data[[i]]) - mean(sample_data[[p]])) /
        sqrt((var(sample_data[[i]]) / length(sample_data[[i]])) +
               (var(sample_data[[p]]) / length(sample_data[[p]])))
    }))
  }
  CB <- max(CB_values)

  quantile_value <- quantile(T_max, probs = 1 - significance_level, na.rm = TRUE)
  p_value <- mean(T_max >= CB, na.rm = TRUE)

  # Decision
  if (CB > quantile_value) { result <- "Reject null hypothesis" }
  else { result <- "Do not reject null hypothesis" }

  # Output summary
  return(paste(
    "Critical value:", quantile_value,
    "; CBMax Test statistic:", CB,
    "; p-value:", p_value,
    "; Result:", result
  ))
}



#' Minimum Difference-based Test for Complete Bipartite (CB) Ordered Alternatives
#'
#' Performs a Minimum Difference-based (CBMin) test to evaluate the equality of group means
#' under the Complete Bipartite (CB) order restriction. This test compares the null
#' hypothesis of equal means across all groups against the alternative that the means
#' of the control groups are less than or equal to those of the treatment groups.
#'
#' @param sample_data A list of numeric vectors, where the first \code{k0} elements
#' correspond to the control groups and the remaining elements correspond to the
#' treatment groups.
#' @param significance_level A numeric value between 0 and 1 specifying the
#' significance level of the test (e.g., \code{0.05}).
#' @param k0 An integer specifying the number of control groups.
#' @param n.boot Number of bootstrap replications used to estimate the critical value
#' (default is \code{100000}).
#' @param seed Optional random seed for reproducibility.
#'
#' @return A character string summarizing the results, including the critical value,
#' the CBMin test statistic, the estimated p-value, and the test decision.
#'
#' @details The test statistic is computed as the minimum standardized difference
#' between each treatment group mean and each control group mean under the CB order restriction.
#' The critical value and p-value are estimated via bootstrap resampling from the null distribution.
#' The test evaluates the hypotheses H_0: μ_i = μ_j versus H_1: μ_i ≤ μ_j (at least one strict inequality)
#' for i = 1, ..., k0 (control groups) and j = k0 + 1, ..., p (total treatment groups).
#'
#' @importFrom stats quantile rnorm var
#'
#' @author Subha Halder
#' @export
#'
#' @examples
#' # Generate sample data
#' set.seed(456)
#' control1 <- rnorm(20, mean = 5)
#' control2 <- rnorm(20, mean = 5)
#' treatment1 <- rnorm(20, mean = 6)
#' treatment2 <- rnorm(20, mean = 7)
#' data_list <- list(control1, control2, treatment1, treatment2)
#'
#' # Run CBMin test with fewer bootstrap samples for demonstration
#' CBMin(sample_data = data_list, significance_level = 0.05, k0 = 2, n.boot = 1000)
#'
#' \donttest{
#' # Recommended: use 100000 bootstrap samples for higher precision
#' CBMin(sample_data = data_list, significance_level = 0.05, k0 = 2, n.boot = 100000)
#' }
CBMin <- function(sample_data, significance_level, k0, n.boot = 100000, seed = NULL) {
  if (!is.null(seed)) { set.seed(seed) }

  sample_data <- lapply(sample_data, function(x) x[!is.na(x)])
  num_samples <- n.boot
  num_datasets <- length(sample_data)
  n <- sapply(sample_data, length)
  group_sds <- sapply(sample_data, function(x) sqrt(var(x)))
  T_min <- numeric(num_samples)

  for (u in 1:num_samples) {
    bootstrap_samples <- vector("list", num_datasets)   # Corrected initialization
    for (j in 1:num_datasets) {
      bootstrap_samples[[j]] <- rnorm(n = n[j], mean = 0, sd = group_sds[j])
    }
    T <- numeric(k0)
    for (p in 1:k0) {
      T[p] <- min(sapply((k0+1):num_datasets, function(i) {
        (mean(bootstrap_samples[[i]]) - mean(bootstrap_samples[[p]])) /
          sqrt((var(bootstrap_samples[[i]]) / length(bootstrap_samples[[i]])) +
                 (var(bootstrap_samples[[p]]) / length(bootstrap_samples[[p]])))
      }))
    }
    T_min[u] <- min(T)
  }

  CB_values <- numeric(k0)   # Corrected initialization
  for (p in 1:k0) {
    CB_values[p] <- min(sapply((k0+1):num_datasets, function(i) {
      (mean(sample_data[[i]]) - mean(sample_data[[p]])) /
        sqrt((var(sample_data[[i]]) / length(sample_data[[i]])) +
               (var(sample_data[[p]]) / length(sample_data[[p]])))
    }))
  }
  CB <- min(CB_values)

  quantile_value <- quantile(T_min, probs = 1 - significance_level, na.rm = TRUE)
  p_value <- mean(T_min >= CB, na.rm = TRUE)

  # Decision
  if (CB > quantile_value) { result <- "Reject null hypothesis" }
  else { result <- "Do not reject null hypothesis" }

  # Output summary
  return(paste(
    "Critical value:", quantile_value,
    "; CBMin Test statistic:", CB,
    "; p-value:", p_value,
    "; Result:", result
  ))
}



#' Max-Min Difference-based Test for Complete Bipartite (CB) Ordered Alternatives
#'
#' Performs a Max-Min Difference-based (CBMaxMin) test to evaluate the equality of group means
#' under the Complete Bipartite (CB) order restriction. This test compares the null
#' hypothesis of equal means across all groups against the alternative that the means
#' of the control groups are less than or equal to those of the treatment groups.
#'
#' @param sample_data A list of numeric vectors, where the first \code{k0} elements
#' correspond to the control groups and the remaining elements correspond to the
#' treatment groups.
#' @param significance_level A numeric value between 0 and 1 specifying the
#' significance level of the test (e.g., \code{0.05}).
#' @param k0 An integer specifying the number of control groups.
#' @param n.boot Number of bootstrap replications used to estimate the critical value
#' (default is \code{100000}).
#' @param seed Optional random seed for reproducibility.
#'
#' @return A character string summarizing the results, including the critical value,
#' the CBMaxMin test statistic, the estimated p-value, and the test decision.
#'
#' @details The CBMaxMin statistic is the maximum of the minimum standardized differences
#' between each control group and all treatment groups. Critical value and p-value
#' are estimated via bootstrap under the null of equal means.
#' The test evaluates the hypotheses H_0: μ_i = μ_j versus H_1: μ_i ≤ μ_j (at least one strict inequality)
#' for i = 1, ..., k0 (control groups) and j = k0 + 1, ..., p (total treatment groups).
#'
#' @importFrom stats quantile rnorm var
#'
#' @author Subha Halder
#' @export
#'
#' @examples
#' # Generate sample data
#' set.seed(456)
#' control1 <- rnorm(20, mean = 5)
#' control2 <- rnorm(20, mean = 5)
#' treatment1 <- rnorm(20, mean = 6)
#' treatment2 <- rnorm(20, mean = 7)
#' data_list <- list(control1, control2, treatment1, treatment2)
#'
#' # Run CBMaxMin test with fewer bootstrap samples for demonstration
#' CBMaxMin(sample_data = data_list, significance_level = 0.05, k0 = 2, n.boot = 1000)
#'
#' \donttest{
#' # Recommended: use 100000 bootstrap samples for higher precision
#' CBMaxMin(sample_data = data_list, significance_level = 0.05, k0 = 2, n.boot = 100000)
#' }
CBMaxMin <- function(sample_data, significance_level, k0, n.boot = 100000, seed = NULL) {
  if (!is.null(seed)) { set.seed(seed) }

  sample_data <- lapply(sample_data, function(x) x[!is.na(x)])
  num_datasets <- length(sample_data)
  n <- sapply(sample_data, length)
  group_sds <- sapply(sample_data, function(x) sqrt(var(x)))
  T_maxmin <- numeric(n.boot)

  for (u in 1:n.boot) {
    bootstrap_samples <- vector("list", num_datasets)
    for (j in 1:num_datasets) {
      bootstrap_samples[[j]] <- rnorm(n = n[j], mean = 0, sd = group_sds[j])
    }
    T <- numeric(k0)
    for (p in 1:k0) {
      T[p] <- min(sapply((k0+1):num_datasets, function(i) {
        (mean(bootstrap_samples[[i]]) - mean(bootstrap_samples[[p]])) /
          sqrt((var(bootstrap_samples[[i]]) / length(bootstrap_samples[[i]])) +
                 (var(bootstrap_samples[[p]]) / length(bootstrap_samples[[p]])))
      }))
    }
    T_maxmin[u] <- max(T)
  }

  CB_values <- numeric(k0)
  for (p in 1:k0) {
    CB_values[p] <- min(sapply((k0+1):num_datasets, function(i) {
      (mean(sample_data[[i]]) - mean(sample_data[[p]])) /
        sqrt((var(sample_data[[i]]) / length(sample_data[[i]])) +
               (var(sample_data[[p]]) / length(sample_data[[p]])))
    }))
  }
  CB <- max(CB_values)

  quantile_value <- quantile(T_maxmin, probs = 1 - significance_level, na.rm = TRUE)
  p_value <- mean(T_maxmin >= CB, na.rm = TRUE)

  # Decision
  result <- if (CB > quantile_value) "Reject null hypothesis" else "Do not reject null hypothesis"

  # Output summary
  return(paste(
    "Critical value:", quantile_value,
    "; CBMaxMin Test statistic:", CB,
    "; p-value:", p_value,
    "; Result:", result
  ))
}


#' Min-Max Difference-based Test for Complete Bipartite (CB) Ordered Alternatives
#'
#' Performs a Min-Max Difference-based (CBMinMax) test to evaluate the equality of group means
#' under the Complete Bipartite (CB) order restriction. This test compares the null
#' hypothesis of equal means across all groups against the alternative that all treatment
#' groups are larger than the control groups.
#'
#' @param sample_data A list of numeric vectors, where the first \code{k0} elements
#' correspond to the control groups and the remaining elements correspond to the
#' treatment groups.
#' @param significance_level A numeric value between 0 and 1 specifying the
#' significance level of the test (e.g., \code{0.05}).
#' @param k0 An integer specifying the number of control groups.
#' @param n.boot Number of bootstrap replications used to estimate the critical value
#' (default is \code{100000}).
#' @param seed Optional random seed for reproducibility.
#'
#' @return A character string summarizing the results, including the critical value,
#' the CBMinMax test statistic, the estimated p-value, and the test decision.
#'
#' @details The CBMinMax statistic is the minimum of the maximum standardized differences
#' between each control group and all treatment groups. Critical value and p-value
#' are estimated via bootstrap under the null of equal means.
#' The test evaluates the hypotheses H_0: μ_i = μ_j versus H_1: μ_i ≤ μ_j (at least one strict inequality)
#' for i = 1, ..., k0 (control groups) and j = k0 + 1, ..., p (total treatment groups).
#'
#' @importFrom stats quantile rnorm var
#'
#' @author Subha Halder
#' @export
#'
#' @examples
#' # Generate sample data
#' set.seed(456)
#' control1 <- rnorm(20, mean = 5)
#' control2 <- rnorm(20, mean = 5)
#' treatment1 <- rnorm(20, mean = 6)
#' treatment2 <- rnorm(20, mean = 7)
#' data_list <- list(control1, control2, treatment1, treatment2)
#'
#' # Run CBMinMax test with fewer bootstrap samples for demonstration
#' CBMinMax(sample_data = data_list, significance_level = 0.05, k0 = 2, n.boot = 1000)
#'
#' \donttest{
#' # Recommended: use 100000 bootstrap samples for higher precision
#' CBMinMax(sample_data = data_list, significance_level = 0.05, k0 = 2, n.boot = 100000)
#' }
CBMinMax <- function(sample_data, significance_level, k0, n.boot = 100000, seed = NULL) {
  if (!is.null(seed)) { set.seed(seed) }

  sample_data <- lapply(sample_data, function(x) x[!is.na(x)])
  num_datasets <- length(sample_data)
  n <- sapply(sample_data, length)
  group_sds <- sapply(sample_data, function(x) sqrt(var(x)))
  T_minmax <- numeric(n.boot)

  for (u in 1:n.boot) {
    bootstrap_samples <- vector("list", num_datasets)
    for (j in 1:num_datasets) {
      bootstrap_samples[[j]] <- rnorm(n = n[j], mean = 0, sd = group_sds[j])
    }
    T <- numeric(k0)
    for (p in 1:k0) {
      T[p] <- max(sapply((k0+1):num_datasets, function(i) {
        (mean(bootstrap_samples[[i]]) - mean(bootstrap_samples[[p]])) /
          sqrt((var(bootstrap_samples[[i]]) / length(bootstrap_samples[[i]])) +
                 (var(bootstrap_samples[[p]]) / length(bootstrap_samples[[p]])))
      }))
    }
    T_minmax[u] <- min(T)
  }

  CB_values <- numeric(k0)
  for (p in 1:k0) {
    CB_values[p] <- max(sapply((k0+1):num_datasets, function(i) {
      (mean(sample_data[[i]]) - mean(sample_data[[p]])) /
        sqrt((var(sample_data[[i]]) / length(sample_data[[i]])) +
               (var(sample_data[[p]]) / length(sample_data[[p]])))
    }))
  }
  CB <- min(CB_values)

  quantile_value <- quantile(T_minmax, probs = 1 - significance_level, na.rm = TRUE)
  p_value <- mean(T_minmax >= CB, na.rm = TRUE)

  # Decision
  result <- if (CB > quantile_value) "Reject null hypothesis" else "Do not reject null hypothesis"

  # Output summary
  return(paste(
    "Critical value:", quantile_value,
    "; CBMinMax Test statistic:", CB, 
    "; p-value:", p_value,
    "; Result:", result
  ))
}








