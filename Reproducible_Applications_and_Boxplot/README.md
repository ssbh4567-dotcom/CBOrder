# `CB_Order_Applications_1_and_2.R`

This script reproduces both real-data applications in the paper using the `CBOrder` package: Application 1 (Vitamin C tooth growth) and Application 2 (cholesterol reduction).

## Requirements
- R (version used for the paper: [4.3.3])
- `CBOrder` package — installed automatically from GitHub if not already available:
  ```r
  devtools::install_github("ssbh4567-dotcom/CBOrder")
  ```
- `devtools` (installed automatically if missing)

## Data
Input data are hardcoded directly in the script as numeric vectors (see `Applcations_Data/` for the same values as standalone CSV files with source citations).

## What the script does

**Application 1 — Tooth Growth (lines ~6–115)**
1. Defines the 6 treatment groups (μ₁–μ₆) split into control partition P1 (`k0 = 2`) and treatment partition P2.
2. Computes summary statistics → **Table 10**.
3. Runs Shapiro–Wilk normality tests for each group.
4. Runs the complete bipartite order-restricted tests (`CBLRT`, `CBMax`, `CBMin`, `CBMaxMin`, `CBMinMax`) → **Table 11**. Uses `seed = 456`; ~10 minutes total runtime (100,000 bootstrap repetitions).
5. Computes simultaneous lower confidence bounds for all P2-vs-P1 contrasts.

**Application 2 — Cholesterol Reduction (lines ~130–229)**
1. Defines the 5 groups (μ₁–μ₃, μ_D, μ_E) split into test-drug partition P1 (`k0 = 3`) and comparator partition P2.
2. Computes summary statistics → **Table 12**.
3. Runs Shapiro–Wilk normality tests for each group.
4. Runs the same five complete bipartite tests → **Table 13**. Uses `seed = 456`; ~10 minutes total runtime.
5. Computes simultaneous lower confidence bounds for all P2-vs-P1 contrasts.

**Plots (lines ~244–end)**
Generates the two boxplots used in the paper (odontoblast lengths; cholesterol reduction), arranged side by side in a 1×2 layout.

## How to run
One can download the R script and then open in RStudio or R, and run the full code.

Runtime: approximately 20 minutes total (10 min per application), due to bootstrap-based critical value estimation.

## Output
All summary tables and test results print directly to the R console in the order they appear in the paper (Tables 10–13); the two boxplots render in the graphics device.

## Notes
- Critical values for the confidence bounds (`c_val`) are hardcoded per application — these correspond to the `CBMax` test's critical value at α = 0.05 and must match the value produced by `CBMax()` above them.
- Both applications must be run within the same R session/order if reusing variable names, since `data_list` and `group_names` are reassigned between Application 1 and Application 2.
