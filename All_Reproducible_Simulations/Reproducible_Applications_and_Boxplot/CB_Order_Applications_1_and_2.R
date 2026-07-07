# ==============================================================================
# R Script for Applications
# ==============================================================================


# ==============================================================================
# Application 1: Vitamin C Tooth Growth Analysis Complete Bipartite 
# Order-Restricted Inference
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. DATA PREPARATION (P1: Controls, P2: Treatments)
# ------------------------------------------------------------------------------
# P1: CONTROL GROUPS (0.5 mg Baseline Doses)
vc_05 <- c(4.2, 11.5, 7.3, 5.8, 6.4, 10.0, 11.2, 11.2, 5.2, 7.0)  # mu_1
oj_05 <- c(15.2, 21.5, 17.6, 9.7, 14.5, 10.0, 8.2, 9.4, 16.1, 9.7) # mu_2

# P2: TREATMENT GROUPS (1.0 mg and 2.0 mg Higher Doses)
vc_10 <- c(16.1, 16.1, 15.2, 17.3, 22.1, 17.3, 13.6, 14.5, 18.8, 15.5) # mu_3
vc_20 <- c(23.6, 18.5, 33.9, 25.5, 26.4, 32.1, 26.7, 21.5, 23.3, 29.1) # mu_4
oj_10 <- c(19.7, 23.3, 23.6, 26.4, 20.0, 25.2, 25.8, 21.2, 14.5, 27.3) # mu_5
oj_20 <- c(25.5, 26.4, 22.4, 24.5, 24.8, 30.9, 26.4, 27.3, 29.4, 23.0) # mu_6

# Combine into a single list corresponding to mu_1 through mu_6
data_list <- list(vc_05, oj_05, vc_10, vc_20, oj_10, oj_20)
group_names <- c("VC (0.5 mg)", "OJ (0.5 mg)", "VC (1.0 mg)",
                 "VC (2.0 mg)", "OJ (1.0 mg)", "OJ (2.0 mg)")

# ------------------------------------------------------------------------------
# 2. DESCRIPTIVE STATISTICS (Matches Table 10)
# ------------------------------------------------------------------------------
cat("\n--- Table 10: Summary Statistics ---\n")
summary_table <- data.frame(
  Treatment = group_names,
  SampleSize = sapply(data_list, length),
  SampleMean = round(sapply(data_list, mean), 2),
  SampleVariance = round(sapply(data_list, var), 5)
)
print(summary_table, row.names = FALSE)

# ------------------------------------------------------------------------------
# 3. NORMALITY ASSUMPTION CHECKS
# ------------------------------------------------------------------------------
cat("\n--- Shapiro-Wilk Test p-values ---\n")
# Extracting exact p-values to match the manuscript text
p_values <- sapply(data_list, function(x) shapiro.test(x)$p.value)
names(p_values) <- group_names
print(round(p_values, 4))

# ------------------------------------------------------------------------------
# 4. EXPLORATORY VISUALIZATION (Boxplot)
# ------------------------------------------------------------------------------
# Consolidated data frame for plotting
plot_data <- data.frame(
  length = unlist(data_list),
  group = factor(rep(group_names, times = sapply(data_list, length)), levels = group_names)
)

# ------------------------------------------------------------------------------
# 5. COMPLETE BIPARTITE TESTS (Matches Table 11)
# ------------------------------------------------------------------------------
# Ensure the custom package is installed
if (!require("CBOrder")) {
  if (!require("devtools")) install.packages("devtools")
  devtools::install_github("ssbh4567-dotcom/CBOrder", force = TRUE)
  library(CBOrder)
}

cat("\n--- Table 11: Complete Bipartite Order-Restricted Tests ---\n")
# Note: Ensure the function documentation specifies that 100,000 bootstrap
# repetitions are utilized (either internally or via an argument like B=100000).

alpha <- 0.05
k0 <- 2  # First two elements in data_list are the P1 controls

cat("\n--- # NOTE: Computational time is approx 10 min (total) for Table 11 ---\n")

print(CBLRT(data_list, alpha, k0), seed = 456)
print(CBMax(data_list, alpha, k0), seed = 456)
print(CBMin(data_list, alpha, k0), seed = 456)
print(CBMaxMin(data_list, alpha, k0), seed = 456)
print(CBMinMax(data_list, alpha, k0), seed = 456)

# ------------------------------------------------------------------------------
# 6. SIMULTANEOUS CONFIDENCE INTERVALS
# ------------------------------------------------------------------------------
cat("\n--- Simultaneous Lower Confidence Bounds ---\n")

# Critical value corresponding to the CBMax test
c_val <- 2.68250

# Helper function to compute the lower bound of the contrast: mu_t - mu_s
get_lower_bound <- function(trt, ctrl) {
  point_estimate <- mean(trt) - mean(ctrl)
  se <- sqrt((var(trt)/length(trt)) + (var(ctrl)/length(ctrl)))
  return(point_estimate - (c_val * se))
}

# Calculate bounds for all t in P2 (Treatments) vs s in P1 (Controls)
bounds <- data.frame(
  Contrast = c("mu_3 - mu_1 (VC 1.0 vs VC 0.5)", "mu_4 - mu_1 (VC 2.0 vs VC 0.5)",
               "mu_5 - mu_1 (OJ 1.0 vs VC 0.5)", "mu_6 - mu_1 (OJ 2.0 vs VC 0.5)",
               "mu_3 - mu_2 (VC 1.0 vs OJ 0.5)", "mu_4 - mu_2 (VC 2.0 vs OJ 0.5)",
               "mu_5 - mu_2 (OJ 1.0 vs OJ 0.5)", "mu_6 - mu_2 (OJ 2.0 vs OJ 0.5)"),
  LowerBound = round(c(
    get_lower_bound(vc_10, vc_05), get_lower_bound(vc_20, vc_05),
    get_lower_bound(oj_10, vc_05), get_lower_bound(oj_20, vc_05),
    get_lower_bound(vc_10, oj_05), get_lower_bound(vc_20, oj_05),
    get_lower_bound(oj_10, oj_05), get_lower_bound(oj_20, oj_05)
  ), 5),
  UpperBound = rep("Inf", 8)
)

print(bounds, row.names = FALSE)







##################################################################################################







# ==============================================================================
# Application 2: Cholesterol Reduction Analysis
# Complete Bipartite Order-Restricted Inference
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. DATA PREPARATION (P1: Test Drugs, P2: Active Comparators)
# ------------------------------------------------------------------------------
# P1: "CONTROL" GROUPS (Test Drug Regimens) - k0 = 3
one_time   <- c(3.8612, 10.3868,  5.9059,  3.0609,  7.7204,
                2.7139,  4.9243,  2.3039,  7.5301,  9.4123)  # mu_1 (20 mg)
two_times  <- c(10.3993,  8.6027, 13.6320,  3.5054,  7.7703,
                8.6266,  9.2274,  6.3159, 15.8258,  8.3443)  # mu_2 (10 mg)
four_times <- c(13.9621, 13.9606, 13.9176,  8.0534, 11.0432,
                12.3692, 10.3921,  9.0286, 12.8416, 18.1794) # mu_3 (5 mg)

# P2: "TREATMENT" GROUPS (Active Comparators)
drug_D     <- c(16.9819, 15.4576, 19.9793, 14.7389, 13.5850,
                10.8648, 17.5897,  8.8194, 17.9635, 17.6316) # mu_D
drug_E     <- c(21.5119, 27.2445, 20.5199, 15.7707, 22.8850,
                23.9527, 21.5925, 18.3058, 20.3851, 17.3071) # mu_E

# Combine into a single list (Must be in order: P1 groups, then P2 groups)
data_list <- list(one_time, two_times, four_times, drug_D, drug_E)
group_names <- c("Test Drug (20 mg)", "Test Drug (10 mg)", "Test Drug (5 mg)",
                 "Drug D", "Drug E")

# ------------------------------------------------------------------------------
# 2. DESCRIPTIVE STATISTICS (Matches Table 12
# ------------------------------------------------------------------------------
cat("\n--- Table 12: Summary Statistics ---\n")
summary_table <- data.frame(
  Treatment = group_names,
  SampleSize = sapply(data_list, length),
  SampleMean = round(sapply(data_list, mean), 5),
  SampleVariance = round(sapply(data_list, var), 5)
)
# Reordering visually for console to match the exact row order of Table 12 in text
print(summary_table[c(4, 5, 1, 2, 3), ], row.names = FALSE)

# ------------------------------------------------------------------------------
# 3. STATISTICAL ASSUMPTION CHECKS
# ------------------------------------------------------------------------------
cat("\n--- Shapiro-Wilk Test p-values (Normality) ---\n")
p_values <- sapply(data_list, function(x) shapiro.test(x)$p.value)
names(p_values) <- group_names
print(round(p_values, 4))



# ------------------------------------------------------------------------------
# 5. COMPLETE BIPARTITE TESTS (Matches Table 13)
# ------------------------------------------------------------------------------

# Ensure the custom package is installed
if (!require("CBOrder")) {
  if (!require("devtools")) install.packages("devtools")
  devtools::install_github("ssbh4567-dotcom/CBOrder", force = TRUE)
  library(CBOrder)
}

cat("\n--- Table 13: Complete Bipartite Order-Restricted Tests ---\n")
alpha <- 0.05
k0 <- 3  # First three elements in data_list are the P1 test drug regimens

cat("\n--- # NOTE: Computational time is approx 10 min (total) for Table 13 ---\n")
print(CBLRT(data_list, alpha, k0), seed = 456)
print(CBMax(data_list, alpha, k0), seed = 456)
print(CBMin(data_list, alpha, k0), seed = 456)
print(CBMaxMin(data_list, alpha, k0), seed = 456)
print(CBMinMax(data_list, alpha, k0), seed = 456)


# ------------------------------------------------------------------------------
# 6. SIMULTANEOUS CONFIDENCE INTERVALS
# ------------------------------------------------------------------------------
cat("\n--- Simultaneous Lower Confidence Bounds ---\n")

# Critical value corresponding to the CBMax test
c_val <- 2.54549

# Helper function to compute the lower bound of the contrast: mu_t - mu_s
get_lower_bound <- function(trt, ctrl) {
  point_estimate <- mean(trt) - mean(ctrl)
  se <- sqrt((var(trt)/length(trt)) + (var(ctrl)/length(ctrl)))
  return(point_estimate - (c_val * se))
}

# Calculate bounds for all t in P2 (Drugs D, E) vs s in P1 (Test Drugs 1, 2, 3)
bounds <- data.frame(
  Contrast = c("mu_D - mu_1 (Drug D vs 20 mg)", "mu_D - mu_2 (Drug D vs 10 mg)", "mu_D - mu_3 (Drug D vs 5 mg)",
               "mu_E - mu_1 (Drug E vs 20 mg)", "mu_E - mu_2 (Drug E vs 10 mg)", "mu_E - mu_3 (Drug E vs 5 mg)"),
  LowerBound = round(c(
    get_lower_bound(drug_D, one_time), get_lower_bound(drug_D, two_times), get_lower_bound(drug_D, four_times),
    get_lower_bound(drug_E, one_time), get_lower_bound(drug_E, two_times), get_lower_bound(drug_E, four_times)
  ), 5),
  UpperBound = rep("Inf", 6)
)

print(bounds, row.names = FALSE)







##################################################################################################






# ------------------------------------------------------------------------------
# 2. CONFIGURE PLOT LAYOUT AND COLORS
# ------------------------------------------------------------------------------
# Set up a 1x2 grid.
# The 'mar' argument adds extra space at the bottom (margin 1) so the
# perpendicular labels (las=2) do not get cut off.
par(mfrow = c(1, 2), mar = c(10, 5, 4, 2) + 0.1)

cols1 <- c("violet", "blue", "skyblue", "green", "yellow")
cols2 <- c("violet", "blue", "skyblue", "green", "yellow", "orange")

# ------------------------------------------------------------------------------
# 3. GENERATE PLOT 1: TOOTH GROWTH
# ------------------------------------------------------------------------------
boxplot(vc_05, oj_05, vc_10, vc_20, oj_10, oj_20,
        col = cols2,
        main = "Boxplot of odontoblast lengths",
        xlab = "Different Vitamin C delivery methods and dosages",
        ylab = "Odontoblast length",
        notch = FALSE,
        border = "black",
        las = 2,           # Perpendicular labels
        lwd = 2,           # Thicker box borders and whiskers
        cex.main = 1.8,    # Enlarge main title
        cex.lab = 1.6,     # Enlarge axis labels
        cex.axis = 1.3)    # Enlarge axis tick labels


# Use mtext to place the X-axis label lower down, avoiding text overlap
mtext("Different Vitamin C delivery methods and dosages", side = 1, line = 8, cex = 1.2)

legend("topleft",
       legend = c("VC (0.5 mg)",
                  "OJ (0.5 mg)",
                  "VC (1.0 mg)",
                  "VC (2.0 mg)",
                  "OJ (1.0 mg)",
                  "OJ (2.0 mg)"),
       fill = cols2,
       cex = 1.4,
       border = "black",
       bty = "n")

# ------------------------------------------------------------------------------
# 4. GENERATE PLOT 2: CHOLESTEROL REDUCTION
# ------------------------------------------------------------------------------
boxplot(one_time, two_times, four_times, drug_E, drug_D,
        col    = cols1,
        main = "Boxplot of cholesterol reduction",
        xlab = "Different test drug dosages and two competing drugs",
        ylab = "Reduction of cholesterol",
        notch = FALSE,
        border = "black",
        las = 2,           # Perpendicular labels
        lwd = 2,           # Thicker box borders and whiskers
        cex.main = 1.8,    # Enlarge main title
        cex.lab = 1.6,     # Enlarge Y-axis label
        cex.axis = 1.3)    # Enlarge X and Y axis tick numbers/names


# Use mtext to place the X-axis label lower down, avoiding text overlap
mtext("Different test drug dosages and two competing drugs", side = 1, line = 8, cex = 1.2)

legend("topleft",
       legend = c("One time daily (test drug)","Two time daily (test drug)",
                  "Four time daily (test drug)","Drug E", "Drug D"),
       fill = cols1,
       cex = 1.4,
       border = "black",
       bty = "n")
























