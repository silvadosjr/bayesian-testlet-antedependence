## =============================================================================
## Real Data Analysis — Small Reading Test (Cross-Validation / Holdout Students)
## Bayesian Testlet-2PP (probit) with local dependence vs. standard 2PP (probit)
##
## This script:
##  1) Loads a TRAIN split of the Small Reading test responses (ld_train.rds),
##  2) Defines a 3-testlet structure (K=3, nk=(2,2,2)) and required indexing,
##  3) Fits two Stan models on TRAIN:
##     - Testlet-2PP with local dependence (Testlet2PP_HT_Diag.stan),
##     - Standard 2PP under local independence (2PPModelVec_Diag.stan),
##  4) Saves fitted objects as .rds,
##  5) Loads a TEST (holdout) split of students (ld_test.rds),
##  6) Computes posterior predictive residual dependence diagnostics on TEST:
##     - Q3 ECDF envelopes by testlet and globally,
##     - Delta-Q3 heatmaps (Testlet vs 2PP),
##     - Q3bar replicate distributions and model superiority probabilities,
##  7) Saves all figures to Fit/Figures.
##
## Required project files:
##  - Real Data Analysis/HelpersRealDataAnal.R
##  - Real Data Analysis/ld_train.rds
##  - Real Data Analysis/ld_test.rds
##  - Programs/Testlet2PP_HT_Diag.stan
##  - Programs/2PPModelVec_Diag.stan
##
## Notes:
##  - Missing responses are set to 0 for Stan input (check if appropriate).
##  - Item columns are reordered to match the declared testlet structure.
##  - This script uses 1 chain + thinning for speed; increase chains for reporting.
## =============================================================================


## -------------------- 0) Paths and helper sources ----------------------------
root_local <- "~/GitHub/bayesian-testlet-antedependence"
path_project <- "~/GitHub/bayesian-testlet-antedependence/Real Data Analysis/Small Reading Test"
stopifnot(dir.exists(root_local), dir.exists(path_project))

setwd(path_project)


## ------------------------- 1) Setup and packages -----------------------------
suppressPackageStartupMessages({
  library(rstan)
  library(rstansim)
  library(shinystan)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tmvtnorm)
  library(mvtnorm)
  library(here)
  library(loo)
  library(ggplot2)
  library(purrr)
  library(tibble)
})

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## Simple timing utilities for logging blocks
tic <- function(msg) { cat(sprintf("\n[START] %s ... %s\n", msg, Sys.time())); Sys.time() }
toc <- function(t0)  { cat(sprintf("[ END ] Elapsed: %s\n", Sys.time() - t0)) }

## Path to Stan programs
pathProgram <- here(root_local, "Programs")

## Source helper functions used for indexing + diagnostics (Q3, envelopes, heatmaps, etc.)
source(here(path_project, "HelpersRealDataAnal.R"))

## Output directories
pathFit <- here(path_project, "Fit")
dir.create(pathFit, showWarnings = FALSE, recursive = TRUE)
stopifnot(dir.exists(pathFit))

saveFigures <- here(pathFit, "Figures")
dir.create(saveFigures, showWarnings = FALSE, recursive = TRUE)
stopifnot(dir.exists(saveFigures))


## ------------------- 2) Load TRAIN data — Small Reading test -----------------
## TRAIN split (students)
file_train <- file.path(path_project, "ld_train.rds")
stopifnot(file.exists(file_train))
base_train <- readRDS(file_train)

## Response matrix (drop first column, typically an ID)
mYc   <- base_train[, -1]
mYcNA <- as.data.frame(mYc)

## Replace missing responses with 0 for Stan input.
## IMPORTANT: confirm this missingness strategy is appropriate for your application.
mYcNA[is.na(mYc)] <- 0L

## Basic checks
n  <- nrow(mYcNA)
vI <- ncol(mYcNA)
stopifnot(n > 0, vI > 0)

## Reorder items to match the intended testlet ordering (must be consistent across train/test)
mYcNA <- mYcNA[, c(1, 6, 4, 5, 2, 3)]


## -------------------- 3) Testlet structure and indices -----------------------
## Declared testlet structure:
##  - K: number of testlets
##  - nk: lengths of each testlet
##  - dk: starting positions (if required by the Stan program)
##  - is_HU: flags controlling HU/HT-style rho indexing used by helper functions
K  <- 3
nk <- c(2, 2, 2)
dk <- c(1, 3, 5)
is_HU <- c(TRUE, TRUE, TRUE)  # kept as provided

## Build rho indexing (rho_len, rho_start, etc.) used by Stan + postprocessing
idx <- make_rho_index(nk, is_HU)

## Items that do NOT belong to any testlet (none in this example)
ind_items <- integer(0L)

## Consistency check: testlet items are complement of ind_items
all_items     <- seq_len(vI)
testlet_items <- setdiff(all_items, ind_items)
stopifnot(length(testlet_items) == sum(nk))

## Split items into testlets according to nk (order matters for diagnostics)
idx_testlets <- split(sort(testlet_items), rep(seq_along(nk), times = nk))
stopifnot(
  sum(lengths(idx_testlets)) == length(testlet_items),
  length(idx_testlets) == K
)


## --------------------- 4) Data lists for Stan --------------------------------
## Testlet-2PP (probit) with local dependence (Diag version, includes log_lik)
data_testlet <- list(
  I = vI, N = n, K = K, dk = dk, nk = nk,
  ind_items = ind_items, n_ind = vI - sum(nk),
  Y = mYcNA,
  sigma_a = .6, sigma_b = 4, sigma_rho = 1,
  rho_len = idx$rho_len, S_mc = 200,
  rho_start = idx$rho_start
)

## Standard 2PP (probit), no local dependence (includes log_lik)
data_2pp <- list(
  N = n, vI = vI, Y = mYcNA,
  sigma_a = .6, sigma_b = 4
)


## ---------------------- 5) Initial values and sampling -----------------------
## Standardized raw score used to initialize theta
scores <- scale(rowSums(mYcNA))[ ,1]

## (i) Testlet-2PP includes rho_global
init_testlet <- function() {
  list(
    theta = as.numeric(scores),
    a     = rep(0.1, vI),
    b     = rep(0.1, vI),
    rho_global = rep(0.1, sum(idx$rho_len))
  )
}

## (ii) Standard 2PP has no rho_global
init_2pp <- function() {
  list(
    theta = as.numeric(scores),
    a     = rep(0.1, vI),
    b     = rep(0.1, vI)
  )
}

## Monitored parameters (log_lik included for LOO / model comparison)
pars_testlet <- c("a", "b", "theta", "rho_global", "log_lik")
pars_2pp     <- c("a", "b", "theta", "log_lik")

## NUTS settings (kept modest for speed)
nChains       <- 1
burnInSteps   <- 1000
thinSteps     <- 20
numSavedSteps <- 1000
nIter         <- ceiling(burnInSteps + numSavedSteps * thinSteps)
ctrl_nuts     <- list(adapt_delta = 0.8, max_treedepth = 10)


## --------------------------- 6) Fit models on TRAIN --------------------------
## (a) Testlet-2PP
t0 <- tic("Fitting Testlet-2PP")
fit_testlet <- stan(
  data   = data_testlet,
  file   = file.path(pathProgram, "Testlet2PP_HT_Diag.stan"),
  init   = init_testlet,
  chains = nChains, pars = pars_testlet,
  iter   = nIter, warmup = burnInSteps, thin = thinSteps,
  control = ctrl_nuts
)
toc(t0)

## (b) Standard 2PP
t0 <- tic("Fitting 2PP")
fit_2pp <- stan(
  data   = data_2pp,
  file   = file.path(pathProgram, "2PPModelVec_Diag.stan"),
  init   = init_2pp,
  chains = nChains, pars = pars_2pp,
  iter   = nIter, warmup = burnInSteps, thin = thinSteps,
  control = list(adapt_delta = 0.8)
)
toc(t0)

## Save fitted objects
saveRDS(fit_testlet, file.path(pathFit, "ResultTestlet2PP_SmallReading_CrossValidation.rds"))
saveRDS(fit_2pp,     file.path(pathFit, "Result2PP_SmallReading_CrossValidation.rds"))

## Reload (useful if running only postprocessing)
fit_testlet <- readRDS(file.path(pathFit, "ResultTestlet2PP_SmallReading_CrossValidation.rds"))
fit_2pp     <- readRDS(file.path(pathFit, "Result2PP_SmallReading_CrossValidation.rds"))


## ------------------------- 7) Holdout (students) -----------------------------
## Load TEST split (students)
file_test <- file.path(path_project, "ld_test.rds")
stopifnot(file.exists(file_test))
base_test <- readRDS(file_test)

## Keep same item columns (drop ID) and apply same missingness handling
mYc_test   <- base_test[, -1]
mYcNA_test <- as.data.frame(mYc_test)
mYcNA_test[is.na(mYc_test)] <- 0L
N_test <- nrow(mYcNA_test)  # number of holdout students

## IMPORTANT: reorder items identically to TRAIN to preserve testlet mapping
mYcNA_test <- mYcNA_test[, c(1, 6, 4, 5, 2, 3)]


## ------------------------- 8) Basic post-processing --------------------------
## Extract posterior draws (arrays) needed by helper functions:
##  - build_params_testlet_from_arrays()
##  - build_params_2pp_from_arrays()
post_testlet <- rstan::extract(
  fit_testlet,
  pars = c("a", "b", "rho_global"),
  permuted = TRUE
)

post_2pp <- rstan::extract(
  fit_2pp,
  pars = c("a", "b"),
  permuted = TRUE
)


## -------------------- Q3 envelopes by TESTLET (facets) -----------------------
## Compute ECDF envelopes of Q3 residual correlations on holdout data,
## by testlet (global = FALSE).
env_test <- q3_envelope_data(
  post_testlet,
  function(draw_id, ext) build_params_testlet_from_arrays(draw_id, ext, idx, nk, K),
  mYc_test,
  idx_testlets,
  model_label = "Testlet 2PP",
  ind_items   = integer(0),
  R_eval      = 200,
  seed        = 123,
  global      = FALSE
)

env_2pp <- q3_envelope_data(
  post_2pp,
  function(draw_id, ext) build_params_2pp_from_arrays(draw_id, ext),
  mYc_test,
  idx_testlets,
  model_label = "2PP",
  ind_items   = integer(0),
  R_eval      = 200,
  seed        = 123,
  global      = FALSE
)

## Combine for plotting
df_plot <- dplyr::bind_rows(env_test, env_2pp)

## Plot in grayscale and save
p_env <- plot_q3_ecdf_envelope_gray(df_plot)
ggsave(
  file.path(saveFigures, "ECDF_SmallReading.pdf"),
  p_env,
  width  = 10,
  height = 5,
  units  = "in"
)


## ------------------------- Q3 envelopes — GLOBAL -----------------------------
## Same ECDF envelopes but pooling correlations globally (global = TRUE).
env_test_g <- q3_envelope_data(
  post_testlet,
  function(draw_id, ext) build_params_testlet_from_arrays(draw_id, ext, idx, nk, K),
  mYc_test,
  idx_testlets,
  model_label = "Testlet 2PP",
  ind_items   = integer(0),
  R_eval      = 200,
  global      = TRUE,
  seed        = 123
)

env_2pp_g <- q3_envelope_data(
  post_2pp,
  function(draw_id, ext) build_params_2pp_from_arrays(draw_id, ext),
  mYc_test,
  idx_testlets,
  model_label = "2PP",
  ind_items   = integer(0),
  R_eval      = 200,
  global      = TRUE,
  seed        = 123
)

p_env_g <- plot_q3_ecdf_envelope_gray(bind_rows(env_test_g, env_2pp_g))

ggsave(
  file.path(saveFigures, "ECDF_Global_SmallReading.pdf"),
  p_env_g,
  width  = 18,
  height = 12,
  units  = "cm"
)


## ---------------------- Delta-Q3 heatmaps (TEST) -----------------------------
## Compute delta-Q3 by testlet:
##  - For each model, estimate Q3 structure on holdout data,
##  - Compare via heatmaps (both full and differences-only).
df_testlet <- delta_q3_testlet(
  post_testlet,
  function(draw_id, ext) build_params_testlet_from_arrays(draw_id, ext, idx, nk, K),
  mYc_test,
  idx_testlets,
  R_eval      = 100,
  model_label = "Testlet 2PP"
)

df_2pp <- delta_q3_testlet(
  post_2pp,
  function(draw_id, ext) build_params_2pp_from_arrays(draw_id, ext),
  mYc_test,
  idx_testlets,
  R_eval      = 100,
  model_label = "2PP"
)

## Save heatmaps: full comparison
save_delta_q3_heatmap_pdf(
  df_testlet,
  df_2pp,
  diff_only = FALSE,
  file      = file.path(saveFigures, "Heatmap_DeltaQ3_Testlet_vs_2PP_SmallReading.pdf")
)

## Save heatmaps: differences only
save_delta_q3_heatmap_pdf(
  df_testlet,
  df_2pp,
  diff_only = TRUE,
  file      = file.path(saveFigures, "Heatmap_DeltaQ3_Diferencas_SmallReading.pdf")
)


## ------------------- Q3bar replicates + superiority --------------------------
## Generate replicate distributions of Q3bar (can increase R_eval for stability).
rep_T <- q3bar_replicates(
  post_testlet,
  function(draw_id, ext) build_params_testlet_from_arrays(draw_id, ext, idx, nk, K),
  mYc_test,
  idx_testlets,
  R_eval = 300,
  seed   = 123
)

rep_2P <- q3bar_replicates(
  post_2pp,
  function(draw_id, ext) build_params_2pp_from_arrays(draw_id, ext),
  mYc_test,
  idx_testlets,
  R_eval = 300,
  seed   = 123
)

## Superiority probability of Testlet-2PP over 2PP:
## returns a tibble with columns: Scope, n, p_sup, se, lo95, hi95
tab_sup <- model_superiority_probability(rep_T, rep_2P, idx_testlets)
tab_sup

## Save a grayscale summary plot/table to PDF
save_superiority_pdf(tab_sup, "Prob_Superioridade_Testlet_vs_2PP.pdf")


## --------------------- 9) Model comparison (LOO/WAIC) -----------------
## LOO with moment matching (to improve stability)
loo_with_mm <- function(fit, cores = 4) {
  ll <- loo::extract_log_lik(fit, merge_chains = FALSE)  # S x C x N
  S  <- dim(ll)[1]; C <- dim(ll)[2]; N <- dim(ll)[3]
  dim(ll) <- c(S*C, N)                                  # draws x N
  r_eff <- loo::relative_eff(exp(ll), chain_id = rep(1:C, each = S))
  loo::loo(ll, r_eff = r_eff, moment_match = TRUE, cores = cores)
}

loo_testlet_mm <- loo_with_mm(fit_testlet)
loo_2pp_mm     <- loo_with_mm(fit_2pp)

## Useful tables
loo::pareto_k_table(loo_testlet_mm)
loo::pareto_k_table(loo_2pp_mm)

## Comparison 
comp_loo  <- loo::loo_compare(loo_testlet_mm, loo_2pp_mm)
print(comp_loo, simplify = FALSE, digits = 3)

## WAIC 
waic_testlet <- loo::waic(loo::extract_log_lik(fit_testlet))
waic_2pp     <- loo::waic(loo::extract_log_lik(fit_2pp))
comp_waic    <- loo::loo_compare(waic_testlet, waic_2pp)
print(comp_waic)

