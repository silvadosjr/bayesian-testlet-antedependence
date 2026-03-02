# Bayesian Testlet Antedependence (Probit 2PP)

This repository contains **Stan** + **R** code for Bayesian estimation and evaluation of **probit 2-parameter IRT (2PP)** models with **local item dependence (LID)** in **testlet-based assessments**. The key model extends the standard locally independent 2PP by introducing **within-testlet residual correlation** using structured correlation blocks (e.g., **HU / compound symmetry** and **HT / Toeplitz-by-lag**) and an **antedependence-style conditional representation** for efficient likelihood construction.

In addition to model fitting, the repo includes a full workflow for:
- **Simulation studies** (data generation under Toeplitz / non-Toeplitz blocks, fitting competing models, parameter recovery)
- **Real-data analysis** (fit Testlet-2PP vs standard 2PP, posterior summaries, predictive checks, Q3 envelopes/heatmaps, LOO/WAIC)

---

## Repository structure (high level)

> Folder names follow the layout used in the scripts.

- `Programs/`  
  Stan programs for:
  - Standard probit **2PP** (vectorized versions and diagnostic variants)
  - **Testlet 2PP** with HU/HT correlation blocks (and diagnostic variants)
  - Simulation Stan programs (latent-variable data generation)

- `Simulation Study/`  
  End-to-end simulation pipeline:
  - `simData/` simulated datasets (RDS)
  - `fitData/` fitted model summaries (RDS)
  - `RecoveryResults/` aggregated metrics + plots (bias, |bias|, RMSE, ESS, Rhat; optionally including ρ)

- `Real Data Analysis/`  
  Scripts and helper functions for large-scale educational assessment applications:
  - posterior summaries and plots for `a`, `b`, `theta`, and `rho_global`
  - holdout predictive checks and residual diagnostics (e.g., Q3 envelopes, ΔQ3 heatmaps)
  - model comparison using **LOO** / **WAIC**

---

## Models

### 1) Standard 2PP (probit)
Assumes local independence across items:
$$P(Y_{ni}=1 \mid \theta_n)=\Phi\big(a_i(\theta_n-b_i)\big).$$

### 2) Testlet 2PP (probit) with residual dependence
Items are partitioned into testlets; within each testlet, residuals are correlated via a block correlation matrix:
- **HU**: one correlation parameter per testlet (compound symmetry)
- **HT**: Toeplitz-by-lag (correlation depends on the lag)

Correlation parameters are stored in a single vector `rho_global`, indexed by `rho_len` and `rho_start` (ragged structure across testlets).

---

## Requirements

### R packages
Core:
- `rstan`, `rstansim`, `here`, `tidyverse` (`dplyr`, `tidyr`, `ggplot2`, `stringr`, `purrr`, `tibble`)
Diagnostics / comparison:
- `loo`, `mvtnorm`, `tmvtnorm`, `Matrix`
Optional:
- `shinystan`, `data.table`

Install (example):
```r
install.packages(c(
  "here","tidyverse","stringr","purrr","tibble",
  "Matrix","loo","mvtnorm","tmvtnorm","data.table","shinystan"
))
install.packages("rstan", repos = "https://cloud.r-project.org")
```
---
## Reproducibility

All simulation studies and empirical analyses are fully reproducible using the scripts provided in this repository. Paths are defined relative to the project root using the `here` package whenever possible.

---

## Citation

If you use this repository in academic work, please cite the following manuscript:

Santos, J. R. S. dos, & Andrade, J. A. A. (2026). Bayesian Modeling of Local Item Dependence in IRT Testlet Data Using Antedependence Models. Journal of Educational and Behavioral Statistics, 0(0).

