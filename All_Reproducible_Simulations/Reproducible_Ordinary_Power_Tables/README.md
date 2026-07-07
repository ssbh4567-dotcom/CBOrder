# Reproducible_Ordinary_Power_Tables

Reproduces the empirical ordinary power results for the five complete bipartite tests (`CBLRT`, `CBMax`, `CBMin`, `CBMaxMin`, `CBMinMax`), and their comparison against a benchmark test.

## Contents (in paper order)

| Script | Reproduces |
|---|---|
| `run_table1_power.R` … `run_table10_power.R` | Tables S6.3 – S6.12 (Supplementary Material) — empirical ordinary power under a range of sample-size, variance, and mean-shift configurations |
| `run_comparison_with_benchmark_size_power.R` | Table 8 (Manuscript) and Table S6.13 (Supplementary Material) — empirical size and power compared with a benchmark test, across 5 sample-size and 6 variance configurations |

Each script is fully self-contained: the five test functions are defined within the script itself, so no other files need to be sourced.

## Requirements

* R
* Packages (installed automatically if missing): `Iso`, `parallel`, `doParallel`, `foreach`
* Multiple CPU cores recommended (uses `doParallel`; manuscript used 64 cores)

## Data

No external data — all samples are simulated internally from normal distributions with hardcoded means, standard deviations, and sample sizes.

## What the scripts do

* **`run_table1_power.R` – `run_table10_power.R`**: each fixes one or more sample-size/variance configurations and a base mean-shift vector, then computes empirical power under each of the five tests as the mean shift is scaled from 0 up to 5× (17 steps), one column at a time (LRT, Max, Min, MaxMin, MinMax).
* **`run_comparison_with_benchmark_size_power.R`**: computes empirical size (mean = 0) and power (non-zero mean shift) for the five tests across 5 sample-size configurations (`N1_star`–`N5_star`) and 6 variance configurations (`SD1`–`SD6`), for direct comparison with benchmark tests.

## Key settings

* `M` (bootstrap reps) and `Q` (simulation reps) default to 100 for a quick test run.
* Set both to 5000 to exactly reproduce the manuscript/supplementary results.
* `CORE <- max(1, parallel::detectCores() - 1)` by default (manuscript used `CORE = 64`).

## How to run

Download and run each script independently in R/RStudio; they do not depend on one another. With defaults (`M = Q = 100`), runtime is short. Reproducing exactly (`M = Q = 5000`) is expensive: a single empirical power value takes ~30 minutes for `CB_LRT` and ~9 minutes for each of the other tests on a 64-core server, so a full script can take **several days** on a standard desktop/laptop.

## Output

Console-only: each script prints a labeled line per test/configuration/mean-shift combination with the empirical size or power value, in the same order the corresponding table's columns appear in the paper.

## Notes

* `M`, `Q`, and `CORE` can be edited at the top of each script.
* Sample sizes, means, and variance configurations are hardcoded to match each script's corresponding table exactly.
