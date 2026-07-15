# CBOrder
 
**Testing Against Complete Bipartite (CB) Order Alternatives in One-Way ANOVA under Heteroscedasticity**
 
R package implementing the complete bipartite order-restricted tests (`CBLRT`, `CBMax`, `CBMin`, `CBMaxMin`, `CBMinMax`) and simultaneous lower confidence bounds, plus all code to reproduce the paper's simulations, applications, and figures.
  
## Installation
 
```r
remotes::install_github("ssbh4567-dotcom/CBOrder", force = TRUE)
```
 
## Usage of R package for implementation of the five tests
 
The following steps outline the procedure to implement the proposed tests using the `CBOrder` package, which includes the functions for the Likelihood Ratio Test (LRT) and the four pairwise standardized mean difference-based tests (`CBMax`, `CBMin`, `CBMaxMin`, and `CBMinMax`).
 
```r
# Increase timeout limit to prevent download failures 
options(timeout = 600)

# Safely install and load the CBOrder package from GitHub using 'remotes'
if (!requireNamespace("CBOrder", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("ssbh4567-dotcom/CBOrder", force = TRUE)
}
library(CBOrder)
 
# Perform the proposed tests
CBLRT(data, alpha, k0)
CBMax(data, alpha, k0)
CBMin(data, alpha, k0)
CBMaxMin(data, alpha, k0)
CBMinMax(data, alpha, k0)
```
 
Here, `data` should be supplied as a list of vectors $(\underline{X}_1, \underline{X}_2, \ldots, \underline{X}_k)$, where $\underline{X}_i$ denotes the vector of observations corresponding to the $i$-th treatment group. The argument `alpha` represents the significance level, and `k0` specifies the number of control groups. Each function performs hypothesis testing of homogeneity of mean effects against the Complete Bipartite (CB) order restriction within a heteroscedastic one-way ANOVA framework.
 
## Repository structure
 
```
CBOrder/
├── R/, man/, DESCRIPTION, NAMESPACE     # Package source and documentation
└── All_Reproducible_Simulations/
    ├── Applications_Data/                              # Real-data application datasets (CSV)
    ├── Reproducible_Applications_and_Boxplot/          # Real-data applications & boxplots (Tables 10-13)
    ├── Reproducible_Size_Tables/                       # Empirical size (Tables 2-3)
    ├── Reproducible_Ordinary_Power_Tables/             # Empirical power (Table 8, S6.3-S6.13)
    ├── Reproducible_Penalized_Power_Plot/              # Penalized power figures (Figs 3-7, MATLAB)
    ├── MCSE_Computational_Cost_Sensitivity/            # MCSE, runtime, sensitivity (Tables 4-7, S2.1, S3.2)
    ├── Reproducible_Robustness_Size_Tables/            # Robustness of size to non-normality (S7.15-S7.19)
    ├── Reproducible_Robustness_Ordinary_Power_Tables/  # Robustness of power to non-normality (S7.20-S7.24)
    └── Reproducible_Robustness_Penalized_Power_Plot/   # Robustness power figures (Figs S7.1-S7.5, MATLAB)
```
 
Each subfolder has its own `README.md` with script-level details and run instructions.
 
## Requirements
 
* R (version: 4.3.3); MATLAB for the two penalized-power-plot folders
* R packages (installed automatically where needed): `Iso`, `parallel`, `doParallel`, `foreach`, `LaplacesDemon`, `remotes`
* Multiple CPU cores recommended for bootstrap-based procedures (paper used 64 cores)

## Computational cost
 
With `M = Q = 5000` (bootstrap replications and Monte Carlo simulations, as used in the paper), on a 64-core server, computing a single empirical size/power value takes roughly **~30 minutes for `CBLRT`** and **~9–10 minutes for each of `CBMax`, `CBMin`, `CBMaxMin`, and `CBMinMax`**. Runtime also scales with the number of treatment groups `p`: at `p = 4` each pairwise test takes ~0.6 min and `CBLRT` ~3.0 min, rising to ~4.0 min and ~7.0 min respectively at `p = 20` (100,000 bootstrap resamples; see `MCSE_Computational_Cost_Sensitivity/`). Full reproduction of all tables is therefore expensive and can take several days on a standard desktop/laptop — reduce `M`/`Q` for quick testing.
 

## License
 
GNU General Public License v3.0.
