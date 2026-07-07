# Reproducible_Size_Tables

Reproduces the empirical size results for the five complete bipartite tests (`CBLRT`, `CBMax`, `CBMin`, `CBMaxMin`, `CBMinMax`), reported in Tables 2 and 3 of the manuscript.

## Contents (in paper order)

| Script | Reproduces |
|---|---|
| `run_table1_size.R` | Table 2 — empirical size, `p = 6` groups, `k0 = 2` |
| `run_table2_size.R` | Table 3 — empirical size, `p = 7` groups, `k0 = 4` |

Each script is fully self-contained: the five test functions (`CB_LRT`, `CB_Max`, `CB_Min`, `CB_MaxMin`, `CB_MinMax`) are defined within the script itself, so no other files need to be sourced.

## Requirements

* R
* Packages (installed automatically if missing): `Iso`, `parallel`, `doParallel`, `foreach`
* Multiple CPU cores recommended (uses `doParallel`; manuscript used 64 cores)

## Data

No external data — all samples are simulated internally from normal distributions with hardcoded means, standard deviations, and sample sizes for each configuration.

## What the scripts do

For a fixed sample-size split (`k0` control groups vs. the rest) and equal means (null case), each script runs the five tests across several variance configurations and prints the empirical size for each, one column at a time (LRT, Max, Min, MaxMin, MinMax).

## Key settings

* `M` (bootstrap reps) and `Q` (simulation reps) default to 100 for a quick test run.
* Set both `M` and `Q` to 5000 to exactly reproduce the manuscript.
* `CORE <- max(1, parallel::detectCores() - 1)` by default (manuscript used `CORE = 64`).

## How to run

Download and run the full script in R/RStudio. With defaults (`M = Q = 100`), runtime is short. Reproducing the manuscript exactly (`M = Q = 5000`) is expensive: per the script's own notes, a single empirical size value takes ~30 minutes for `CB_LRT` and ~9 minutes for each of the other tests on a 64-core server — so a full table can take **several days** on a standard desktop/laptop.

## Output

Console-only: for each variance configuration and each test, prints a labeled line with the empirical size value, in the same order the columns appear in the manuscript tables.

## Notes

* `M`, `Q`, and `CORE` can be edited at the top of each script.
* Sample sizes, means, and variance configurations are hardcoded to match the manuscript's table rows exactly.
