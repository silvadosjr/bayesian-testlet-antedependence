## -------------------- 0) Paths and helper sources ----------------------------
## This script runs a full train/test workflow for a large-scale educational assessment:
##   1) Fit two Bayesian probit IRT models on a TRAIN sample:
##        - Testlet 2PP (accounts for local dependence inside predefined testlets)
##        - Standard 2PP (local independence)
##   2) Evaluate out-of-sample residual dependence diagnostics on a TEST sample:
##        - Q3 ECDF envelopes by testlet (and global)
##        - Delta-Q3 heatmaps (including difference-only version)
##        - Q3bar superiority probability summaries and plot
##   3) Compare predictive performance on the TRAIN sample using LOO (moment matching)
##      and (optionally) WAIC.
##
## Notes:
##   - This analysis assumes that the item ordering/columns in TRAIN and TEST match.
##   - The testlet structure (dk, nk, ind_items, is_HU) must be consistent with the data.
##   - Stan programs are read from ~/GitHub/bayesian-testlet-antedependence/Programs.
##
## Adjust the paths below if your local repository structure differs.

root_local   <- "~/GitHub/bayesian-testlet-antedependence"
path_project <- "~/GitHub/bayesian-testlet-antedependence/Real Data Analysis/Large-Scale Educational Assessment"
stopifnot(dir.exists(root_local), dir.exists(path_project))

setwd(path_project)  # keeps compatibility with existing scripts


## ------------------------- 1) Setup and packages -----------------------------
## rstan: model fitting
## rstansim: simulation utilities used elsewhere in the repo
## shinystan: interactive diagnostics (optional, but loaded here for convenience)
## tidyverse stack + mvtnorm/tmvtnorm: used by helper functions and post-processing
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
tic <- function(msg) {
  cat(sprintf("\n[START] %s ... %s\n", msg, Sys.time()))
  Sys.time()
}
toc <- function(t0) {
  cat(sprintf("[ END ] Elapsed: %s\n", Sys.time() - t0))
}

## Path to Stan programs
pathProgram <- here(root_local, "Programs")

## Source helper functions used for indexing + diagnostics (Q3, envelopes, heatmaps, etc.)
## The helper file is expected to live in the repository under "Real Data Analysis/"
source(here("Real Data Analysis", "HelpersRealDataAnal.R"))

## Output directories
pathFit <- here(path_project, "Fit")
dir.create(pathFit, showWarnings = FALSE, recursive = TRUE)
stopifnot(dir.exists(pathFit))

saveFigures <- here(pathFit, "Figures")
dir.create(saveFigures, showWarnings = FALSE, recursive = TRUE)
stopifnot(dir.exists(saveFigures))


## ------------------- 2) Load TRAIN data (3rd year HS) ------------------------
## TRAIN response matrix (0/1, same columns/items as in the TEST file)
file_train <- file.path(path_project, "Binary_1EM_train_n2k.RDS")
stopifnot(file.exists(file_train))

mYcNA <- readRDS(file_train)

## Basic checks
n  <- nrow(mYcNA)  # number of students in TRAIN
vI <- ncol(mYcNA)  # number of items
stopifnot(n > 0, vI > 0)


## -------------------- 3) Testlet structure and indexing ----------------------
## User-defined testlet structure:
##   K      : number of testlets
##   nk[k]  : size (number of items) in testlet k
##   dk[k]  : 1-based starting position of testlet k in the item vector
##   is_HU  : indicator of within-testlet correlation structure
##            - TRUE  : HU (compound symmetry; a single rho for all pairs in the testlet)
##            - FALSE : HT (Toeplitz-by-lag; nk[k]-1 rhos, one per lag)
K  <- 5
nk <- c(3, 2, 2, 2, 2)             # testlet sizes (in the declared order)
dk <- c(2, 8, 12, 14, 16)          # starting indices for each testlet (1-based)
is_HU <- c(FALSE, TRUE, TRUE, TRUE, TRUE)  # testlet 1 is HT; remaining testlets are HU

## Build the "ragged" index mapping for rho_global:
##   idx$rho_len[k]   : number of rho parameters used by testlet k
##                      (1 if HU; nk[k]-1 if HT)
##   idx$rho_start[k] : 1-based start position of testlet k's rho block in rho_global
idx <- make_rho_index(nk, is_HU)

## Independent items (not belonging to any testlet)
ind_items <- c(1, 5:7, 10, 11, 18:26)

## Consistency check: testlet items must match complement of independent items
all_items     <- seq_len(vI)
testlet_items <- setdiff(all_items, ind_items)
stopifnot(length(testlet_items) == sum(nk))

## Partition items by testlet (for plotting and Q3 diagnostics)
idx_testlets <- split(sort(testlet_items), rep(seq_along(nk), times = nk))
stopifnot(
  sum(lengths(idx_testlets)) == length(testlet_items),
  length(idx_testlets) == K
)


## --------------------- 4) Stan data lists ------------------------------------
## (a) Testlet-2PP (probit) with diagnostic log-likelihood (log_lik in GQs)
##     S_mc controls Monte Carlo samples used inside generated quantities.
data_testlet <- list(
  I = vI, N = n, K = K, dk = dk, nk = nk,
  ind_items = ind_items, n_ind = vI - sum(nk),
  Y = mYcNA,
  sigma_a = 0.6, sigma_b = 4, sigma_rho = 1,
  rho_len = idx$rho_len, S_mc = 200,
  rho_start = idx$rho_start
)

## (b) Standard 2PP (probit), with log_lik in generated quantities
data_2pp <- list(
  N = n, vI = vI, Y = mYcNA,
  sigma_a = 0.6, sigma_b = 4
)


## ---------------------- 5) Initial values -----------------------------------
## Initialize theta using standardized raw scores from TRAIN.
scores <- scale(rowSums(mYcNA))[, 1]

## Testlet-2PP init (includes rho_global)
init_testlet <- function() {
  list(
    theta      = as.numeric(scores),
    a          = rep(0.1, vI),
    b          = rep(0.1, vI),
    rho_global = rep(0.1, sum(idx$rho_len))
  )
}

## 2PP init (no rho parameters)
init_2pp <- function() {
  list(
    theta = as.numeric(scores),
    a     = rep(0.1, vI),
    b     = rep(0.1, vI)
  )
}

## Parameters to monitor (include log_lik for LOO/WAIC)
pars_testlet <- c("a", "b", "theta", "rho_global", "log_lik")
pars_2pp     <- c("a", "b", "theta", "log_lik")

## NUTS settings (single chain; adjust for production runs)
nChains       <- 1
burnInSteps   <- 1000
thinSteps     <- 1
numSavedSteps <- 1000
nIter         <- ceiling(burnInSteps + numSavedSteps * thinSteps)
ctrl_nuts     <- list(adapt_delta = 0.8, max_treedepth = 10)


## --------------------------- 6) Fit the models -------------------------------
## (a) Testlet-2PP
t0 <- tic("Fitting Testlet-2PP")
fit_testlet <- stan(
  data    = data_testlet,
  file    = file.path(pathProgram, "Testlet2PP_HT_Diag.stan"),
  init    = init_testlet,
  chains  = nChains,
  pars    = pars_testlet,
  iter    = nIter,
  warmup  = burnInSteps,
  thin    = thinSteps,
  control = ctrl_nuts
)
toc(t0)

## (b) Standard 2PP
t0 <- tic("Fitting 2PP")
fit_2pp <- stan(
  data    = data_2pp,
  file    = file.path(pathProgram, "2PPModelVec_Diag.stan"),
  init    = init_2pp,
  chains  = nChains,
  pars    = pars_2pp,
  iter    = nIter,
  warmup  = burnInSteps,
  thin    = thinSteps,
  control = list(adapt_delta = 0.8)  # default treedepth usually OK
)
toc(t0)

## Save fitted objects (avoid refitting when only diagnostics are needed)
saveRDS(fit_testlet, file.path(pathFit, "ResultTestlet2PP_1EM_CrossValidation2k.rds"))
saveRDS(fit_2pp,     file.path(pathFit, "Result2PP_1EM_CrossValidation2k.rds"))

## Reload fits (convenient for modular execution)
fit_testlet <- readRDS(file.path(pathFit, "ResultTestlet2PP_1EM_CrossValidation2k.rds"))
fit_2pp     <- readRDS(file.path(pathFit, "Result2PP_1EM_CrossValidation2k.rds"))


## ------------------------- 7) Holdout students (TEST) ------------------------
## Load TEST response matrix. Must have the same items/columns as TRAIN.
file_test <- file.path(path_project, "Binary_1EM_test_n2k.RDS")
stopifnot(file.exists(file_test))

mYcNA_test <- readRDS(file_test)
N_test <- nrow(mYcNA_test)


## ------------------------- 8) Posterior draw extraction ----------------------
## Extract posterior draws as arrays for downstream helper functions.
## We include theta because some predictive checks may rely on person parameters.
post_testlet <- rstan::extract(
  fit_testlet,
  pars = c("a", "b", "theta", "rho_global"),
  permuted = TRUE
)

post_2pp <- rstan::extract(
  fit_2pp,
  pars = c("a", "b", "theta"),
  permuted = TRUE
)


## ------------------------- 9) Out-of-sample Q3 diagnostics -------------------
## The functions below are defined in HelpersRealDataAnal.R.
## They typically:
##   - Build predicted response dependence implied by posterior draws
##   - Compute Q3-type residual correlation summaries
##   - Create ECDF envelope plots and delta-Q3 heatmaps
##
## IMPORTANT:
##   - idx_testlets encodes which items belong to each testlet
##   - ind_items is set to integer(0) here because we are evaluating *within-testlet*
##     residual dependence (the helper functions handle the relevant pairs)

## ---- Testlet-level ECDF envelopes (facet per testlet) ----
env_test <- q3_envelope_data(
  post_testlet,
  function(draw_id, ext) build_params_testlet_from_arrays(draw_id, ext, idx, nk, K),
  mYcNA_test,
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
  mYcNA_test,
  idx_testlets,
  model_label = "2PP",
  ind_items   = integer(0),
  R_eval      = 200,
  seed        = 123,
  global      = FALSE
)

df_plot <- dplyr::bind_rows(env_test, env_2pp)

## Improve readability for dense facets by reducing x tick count and spacing panels
p_env <- plot_q3_ecdf_envelope_gray(df_plot) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
  theme(
    axis.text.x     = element_text(size = 8),
    panel.spacing.x = unit(1.5, "lines")
  )

ggsave(file.path(saveFigures, "ECDF_1EM.pdf"), p_env, width = 12, height = 5, units = "in")

## ---- Global ECDF envelope (aggregated across testlets) ----
env_test_g <- q3_envelope_data(
  post_testlet,
  function(draw_id, ext) build_params_testlet_from_arrays(draw_id, ext, idx, nk, K),
  mYcNA_test,
  idx_testlets,
  model_label = "Testlet 2PP",
  ind_items   = integer(0),
  R_eval      = 200,
  global      = TRUE,
  seed        = 234
)

env_2pp_g <- q3_envelope_data(
  post_2pp,
  function(draw_id, ext) build_params_2pp_from_arrays(draw_id, ext),
  mYcNA_test,
  idx_testlets,
  model_label = "2PP",
  ind_items   = integer(0),
  R_eval      = 200,
  global      = TRUE
)

p_env_g <- plot_q3_ecdf_envelope_gray(dplyr::bind_rows(env_test_g, env_2pp_g))

## Note: filename mentions "25k" in your original line; kept as-is for compatibility,
## but this run uses the n2k split (train/test). Rename if preferred.
ggsave(file.path(saveFigures, "ECDF_Global_1EM.pdf"), p_env_g, width = 18, height = 12, units = "cm")


## ---- Delta-Q3 heatmaps (model vs model, and differences only) ----
df_testlet <- delta_q3_testlet(
  post_testlet,
  function(draw_id, ext) build_params_testlet_from_arrays(draw_id, ext, idx, nk, K),
  mYcNA_test,
  idx_testlets,
  R_eval      = 100,
  model_label = "Testlet"
)

df_2pp <- delta_q3_testlet(
  post_2pp,
  function(draw_id, ext) build_params_2pp_from_arrays(draw_id, ext),
  mYcNA_test,
  idx_testlets,
  R_eval      = 100,
  model_label = "2PP"
)

## Heatmap with both models
save_delta_q3_heatmap_pdf(
  df_testlet,
  df_2pp,
  diff_only = FALSE,
  file      = file.path(saveFigures, "Heatmap_DeltaQ3_Testlet_vs_2PP_1EM_2k.pdf")
)

## Heatmap showing only differences
save_delta_q3_heatmap_pdf(
  df_testlet,
  df_2pp,
  diff_only = TRUE,
  file      = file.path(saveFigures, "Heatmap_DeltaQ3_Diferencas_1EM.pdf")
)


## ---- Q3bar replicates and superiority probabilities ----
## These replicates can be used to estimate P(Testlet better than 2PP) globally and by testlet.
rep_T <- q3bar_replicates(
  post_testlet,
  function(draw_id, ext) build_params_testlet_from_arrays(draw_id, ext, idx, nk, K),
  mYcNA_test,
  idx_testlets,
  R_eval = 300,
  seed   = 123
)

rep_2P <- q3bar_replicates(
  post_2pp,
  function(draw_id, ext) build_params_2pp_from_arrays(draw_id, ext),
  mYcNA_test,
  idx_testlets,
  R_eval = 300,
  seed   = 123
)

## Superiority probability table:
## Expected columns: Scope, n, p_sup, se, lo95, hi95
tab_sup <- model_superiority_probability(rep_T, rep_2P, idx_testlets)
tab_sup

## Save a grayscale plot summarizing superiority probabilities
save_superiority_pdf(
  tab_sup,
  file.path(saveFigures, "Prob_Superioridade_Testlet_vs_2PP.pdf")
)


## --------------------- 10) Predictive model comparison (LOO / WAIC) ----------
## We compute LOO using moment matching for better stability in problematic Pareto-k cases.
## IMPORTANT:
##   - log_lik must be present in generated quantities of each Stan program.
##   - merge_chains = FALSE keeps chain structure for relative_eff computation.
loo_with_mm <- function(fit, cores = 4) {
  ll <- loo::extract_log_lik(fit, merge_chains = FALSE)  # S x C x N
  S  <- dim(ll)[1]
  C  <- dim(ll)[2]
  N  <- dim(ll)[3]
  
  ## Convert to (draws x N) for loo()
  dim(ll) <- c(S * C, N)
  
  ## Relative efficiency (needed for r_eff in loo)
  r_eff <- loo::relative_eff(exp(ll), chain_id = rep(1:C, each = S))
  
  loo::loo(ll, r_eff = r_eff, moment_match = TRUE, cores = cores)
}

loo_testlet_mm <- loo_with_mm(fit_testlet)
loo_2pp_mm     <- loo_with_mm(fit_2pp)

## Pareto-k diagnostics tables
loo::pareto_k_table(loo_testlet_mm)
loo::pareto_k_table(loo_2pp_mm)

## LOO comparison (first argument is the "best" baseline in loo_compare output)
comp_loo <- loo::loo_compare(loo_testlet_mm, loo_2pp_mm)
print(comp_loo, simplify = FALSE, digits = 3)

## Optional: WAIC (often less stable than LOO, but included for completeness)
waic_testlet <- loo::waic(loo::extract_log_lik(fit_testlet))
waic_2pp     <- loo::waic(loo::extract_log_lik(fit_2pp))
comp_waic    <- loo::loo_compare(waic_testlet, waic_2pp)
print(comp_waic)

cat("\n[OK] Pipeline completed.\n")
