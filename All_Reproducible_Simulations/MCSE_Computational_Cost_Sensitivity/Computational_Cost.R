# ==============================================================================
# REPLICATION SCRIPT: Computational Cost Benchmark (Time Only) --- Table 6 of manuscript
# ==============================================================================
# ==============================================================================
# NOTE ON COMPUTATIONAL TIME (SCALING WITH 'p')
# ==============================================================================
# The execution time of these tests scales directly with the number of
# treatment levels (p). Based on our benchmarks using an intensive 100,000
# parametric bootstrap resamples:
#
# 1. Low-dimensional setting ($p=4$): The pairwise tests (CBMax, CBMin, CBMaxMin,
# and CBMinMax) are computationally efficient, each requiring approximately 
# 0.40 minutes, whereas the CBLRT is the most computationally intensive,
# requiring approximately 2 minutes.

# 2. High-dimensional setting ($p=20$): The execution time increases to 
# approximately 2.6 minutes for the pairwise tests and 5.0 minutes for the CBLRT.

# Note: The reported execution times are based on our computing system and may 
# vary depending on the system.
#
# Overall, the four pairwise procedures exhibit nearly identical
# computational costs at any given level of 'p' and consistently remain
# substantially faster than the CBLRT.


# ==============================================================================
# Increase timeout limit to 10 minutes to prevent GitHub download failures 
options(timeout = 600)

# Ensure the custom package is loaded using 'remotes' (lighter and faster than devtools)
if (!requireNamespace("CBOrder", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  message("CBOrder not found. Installing from GitHub...")
  remotes::install_github("ssbh4567-dotcom/CBOrder", force = TRUE)
}

library(CBOrder)
# ==============================================================================

alpha <- 0.05
k0 <- 2
n_obs <- 50

# Helper function to run the tests and ONLY print the time
benchmark_time_only <- function(p) {
  cat("\n=======================================================\n")
  cat(" TIME (IN MINUTES) FOR p =", p, "\n")
  cat("=======================================================\n")

  set.seed(456)
  data_list <- lapply(1:p, function(x) rnorm(n_obs, mean = 0, sd = 4))

  time_cblrt <- system.time({
    invisible(CBLRT(data_list, alpha, k0, seed = 456))
  })["elapsed"] / 60
  cat(sprintf("CBLRT    : %.4f\n", time_cblrt))

  time_cbmax <- system.time({
    invisible(CBMax(data_list, alpha, k0, seed = 456))
  })["elapsed"] / 60
  cat(sprintf("CBMax    : %.4f\n", time_cbmax))

  time_cbmin <- system.time({
    invisible(CBMin(data_list, alpha, k0, seed = 456))
  })["elapsed"] / 60
  cat(sprintf("CBMin    : %.4f\n", time_cbmin))

  time_cbmaxmin <- system.time({
    invisible(CBMaxMin(data_list, alpha, k0, seed = 456))
  })["elapsed"] / 60
  cat(sprintf("CBMaxMin : %.4f\n", time_cbmaxmin))

  time_cbminmax <- system.time({
    invisible(CBMinMax(data_list, alpha, k0, seed = 456))
  })["elapsed"] / 60
  cat(sprintf("CBMinMax : %.4f\n", time_cbminmax))
}

# Run the benchmark for all treatment levels matching the manuscript (Table 6)
benchmark_time_only(4)  #Row 1 of Table 6
benchmark_time_only(8)  #Row 2 of Table 6
benchmark_time_only(12) #Row 3 of Table 6
benchmark_time_only(16) #Row 4 of Table 6
benchmark_time_only(20) #Row 5 of Table 6

