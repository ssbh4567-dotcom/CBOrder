# Reproducible_Robustness_Ordinary_Power_Tables

Reproduces the robustness-to-non-normality analysis for the empirical ordinary power of the five complete bipartite tests, evaluating each test under deviations from the assumed normal distribution as the mean shift is scaled up.

## Contents (in paper order)

| Script | Reproduces |
|---|---|
| `Robustness_CBLRT_Power.R` | Table S7.20 (Supplementary Material) — CBLRT |
| `Robustness_CBMax_Power.R` | Table S7.21 (Supplementary Material) — CBMax |
| `Robustness_CBMin_Power.R` | Table S7.22 (Supplementary Material) — CBMin |
| `Robustness_CBMaxMin_Power.R` | Table S7.23 (Supplementary Material) — CBMaxMin |
| `Robustness_CBMinMax_Power.R` | Table S7.24 (Supplementary Material) — CBMinMax |

Each script is fully self-contained: all required functions are defined within the script itself, so no other files need to be sourced.

## Requirements

* R
* Packages (installed automatically if missing): `parallel`, `doParallel`, `foreach`, `LaplacesDemon` (`Robustness_CBLRT_Power.R` also uses `Iso`)
* Multiple CPU cores recommended (uses `doParallel`; manuscript used 64 cores)

## Data

No external data — all samples are simulated internally, under two sample-size configurations (`p = 6`, `k0 = 2`) and a fixed base mean and variance vector.

## What the scripts do

For each of two sample-size configurations, the script fixes a base mean vector and variance vector, then computes the empirical power of the corresponding test as the mean shift is scaled by a multiplier `c` from 0 up to 5 (17 steps), under six data-generating scenarios:

* **Normal (N)** — the assumed distribution, $X_i \sim N(\mu_i, \sigma_i^2)$
* **Mix-1** — $0.15 \cdot N(\mu_i + 4, (\sigma_i + 4)^2) + 0.85 \cdot N(\mu_i, \sigma_i^2)$
* **Mix-2** — $0.15 \cdot N(\mu_i + 4, (\sigma_i + 7)^2) + 0.85 \cdot N(\mu_i, \sigma_i^2)$
* **Mix-3** — $0.30 \cdot N(\mu_i + 4, (\sigma_i + 4)^2) + 0.70 \cdot N(\mu_i, \sigma_i^2)$
* **Mix-4** — $0.30 \cdot N(\mu_i + 4, (\sigma_i + 7)^2) + 0.70 \cdot N(\mu_i, \sigma_i^2)$
* **Laplace** — $X_i \sim \text{Laplace}(\text{location} = -2, \text{scale} = 3)$

## Key settings

* `M` (bootstrap reps) and `Q` (simulation reps) default to 100 for a quick test run.
* Set both to 5000 to exactly reproduce the manuscript/supplementary results.
* `CORE <- max(1, parallel::detectCores() - 1)` by default (manuscript used `CORE = 64`).

## How to run

Download and run each script independently in R/RStudio; they do not depend on one another. With defaults (`M = Q = 100`), runtime is short. Reproducing exactly (`M = Q = 5000`) is expensive: a single empirical power value takes ~30 minutes (CBLRT) or ~9 minutes (CBMax/CBMin/CBMaxMin/CBMinMax) per distribution on a 64-core server, so a full script can take **several days** on a standard desktop/laptop.

## Output

Console-only: for each sample-size configuration and distribution, prints the empirical power at each mean multiplier `c`, in the same order the corresponding table's rows and columns appear in the Supplementary Material.

## Notes

* `M`, `Q`, and `CORE` can be edited at the top of each script.
* Sample sizes, means, variances, and mixture parameters are hardcoded to match each script's corresponding table exactly.
