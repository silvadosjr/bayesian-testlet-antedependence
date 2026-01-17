# ==============================================================================
# Simulation: 2PP with residual dependence (Toeplitz within each testlet)
# + option to also fit the standard 2PP model (local independence)
# ==============================================================================
# This script reads pre-generated simulated datasets (saved as .rds files) and
# fits:
#   (i)  Testlet-2PP (probit) model with within-testlet residual dependence, and/or
#   (ii) Standard 2PP (probit) model assuming local independence,
# across multiple scenario configurations and replications.
#
# Key features:
#  - Supports different prior "strength" levels (weak/moderate/strong).
#  - Pre-compiles Stan models once (faster in loops).
#  - Builds Stan data lists automatically from K, nk, dk and a HU/HT rule.
#  - Saves posterior summaries to disk (one file per replication).
#
# Requirements (project structure):
#  - Programs/Testlet2PP_HTGeral.stan
#  - Programs/2PPModelVec.stan
#  - Simulation Study/simData/  (simulated datasets, created elsewhere)
#  - Simulation Study/fitData/  (created by this script for outputs)
#
# Notes:
#  - The simulation filenames must match the pattern used below:
#      SimTestlet2PP_{K}_{nk}_{tag}_{GEN_TYPE}_{r}.rds
#    and each file must contain a response matrix in dat$Y.
#  - This script uses 1 chain and a long warmup (3000) by default; adjust as needed.
# ==============================================================================

## -----------------------------
## 0) Setup
## -----------------------------
setwd('~/GitHub/bayesian-testlet-antedependence')

suppressPackageStartupMessages({
  library(here)
  library(rstan)
  library(rstansim)
})

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## (Optional) Global reproducibility
# set.seed(1234)

## -----------------------------
## Generation config (used only to build filenames/tags here)
## -----------------------------
## GEN_TYPE options:
##   - "UNS"  : unstructured (non-Toeplitz) within testlets
##   - "TOEP" : Toeplitz-by-lag within testlets
##   - "MIX"  : alternating patterns by testlet
GEN_TYPE <- "TOEP"
MIX_INT  <- 0.85
JITTER   <- 0.07

## Which model(s) to fit?
##   "testlet" : fit only the testlet model
##   "2pp"     : fit only the standard 2PP model
##   "both"    : fit both models
model_to_fit <- "testlet"

## -----------------------------
## 1) Paths
## -----------------------------
programas_dir <- here("Programs")
pathSim       <- here("Simulation Study", "simData")
pathFit       <- here("Simulation Study", "fitData")
dir.create(pathFit, showWarnings = FALSE, recursive = TRUE)

## -----------------------------
## 2) Stan files (compiled once)
## -----------------------------
stan_testlet_file <- file.path(programas_dir, "Testlet2PP_HTGeral.stan")
stan_2pp_file     <- file.path(programas_dir, "2PPModelVec.stan")

## Compile Testlet model if needed
if (model_to_fit %in% c("testlet", "both")) {
  stopifnot(file.exists(stan_testlet_file))
  message("Compiling Testlet model: ", basename(stan_testlet_file))
  sm_testlet <- stan_model(stan_testlet_file)
}

## Compile 2PP model if needed
if (model_to_fit %in% c("2pp", "both")) {
  stopifnot(file.exists(stan_2pp_file))
  message("Compiling 2PP model: ", basename(stan_2pp_file))
  sm_2pp <- stan_model(stan_2pp_file)
}

## -----------------------------
## 3) Prior levels and NUTS control
## -----------------------------
## Choose a prior strength level:
##  - "weak"     : broader priors, more conservative sampling controls
##  - "moderate" : default (often a good balance)
##  - "strong"   : tighter priors, slightly looser sampling controls
prior_level <- "moderate"

## Map prior level to hyperparameters
get_hypers <- function(nivel = c("weak", "moderate", "strong")) {
  nivel <- match.arg(nivel)
  switch(
    nivel,
    "weak"     = list(sigma_a = 2.0, sigma_b = 10, sigma_rho = 10),
    "moderate" = list(sigma_a = 1.0, sigma_b = 4,  sigma_rho = 1),
    "strong"   = list(sigma_a = 0.6, sigma_b = 1,  sigma_rho = 0.5)
  )
}

## Map prior level to NUTS tuning
get_nuts_ctrl <- function(nivel = c("weak", "moderate", "strong")) {
  nivel <- match.arg(nivel)
  switch(
    nivel,
    "weak"     = list(adapt_delta = 0.98, max_treedepth = 12),
    "moderate" = list(adapt_delta = 0.95, max_treedepth = 12),
    "strong"   = list(adapt_delta = 0.90, max_treedepth = 12)
  )
}

hyper <- get_hypers(prior_level)
nutsc <- get_nuts_ctrl(prior_level)

## -----------------------------
## 4) Helpers to build Stan data lists
## -----------------------------
## Build data list for Testlet-2PP (HU/HT within testlets)
montar_data_list_testlet <- function(
    Y, K, nk, dk,
    structure = rep("HT", K),
    hypers = list(sigma_a = 1, sigma_b = 1, sigma_rho = 0.35)
) {
  stopifnot(is.matrix(Y))
  storage.mode(Y) <- "integer"
  N <- nrow(Y)
  I <- ncol(Y)
  
  stopifnot(length(nk) == K, length(dk) == K, length(structure) == K)
  stopifnot(all(dk >= 1L), all(nk >= 1L))
  stopifnot(all(dk + nk - 1L <= I))  # ensure blocks are within 1..I
  
  structure <- toupper(structure)
  stopifnot(all(structure %in% c("HU", "HT")))
  
  ## Identify which items are inside testlets vs independent
  idx_testlet <- unlist(
    Map(function(s, L) seq.int(s, s + L - 1L), dk, nk),
    use.names = FALSE
  )
  idx_all   <- seq_len(I)
  ind_items <- setdiff(idx_all, idx_testlet)
  n_ind     <- length(ind_items)
  
  ## Ragged rho indexing:
  ##  - HU: one correlation parameter
  ##  - HT: L-1 lag correlations
  rho_len <- ifelse(structure == "HU", 1L, nk - 1L)
  if (any(rho_len < 0L)) stop("A testlet has nk[k]=1 with structure='HT'. Use 'HU' in this case.")
  
  ## Starting positions (1-based) for each segment inside rho_global
  rho_start <- integer(K)
  if (K > 0) {
    rho_start[1] <- 1L
    if (K > 1) for (k in 2:K) rho_start[k] <- rho_start[k - 1] + rho_len[k - 1]
  }
  
  stopifnot(all(c("sigma_a", "sigma_b", "sigma_rho") %in% names(hypers)))
  
  list(
    I = as.integer(I),
    K = as.integer(K),
    N = as.integer(N),
    dk = as.array(as.integer(dk)),
    nk = as.array(as.integer(nk)),
    n_ind = as.integer(n_ind),
    ind_items = if (n_ind > 0) as.array(as.integer(ind_items)) else array(1L, 0),
    Y = Y,
    rho_len = as.array(as.integer(rho_len)),
    rho_start = as.array(as.integer(rho_start)),
    sigma_a = as.numeric(hypers$sigma_a),
    sigma_b = as.numeric(hypers$sigma_b),
    sigma_rho = as.numeric(hypers$sigma_rho)
  )
}

## Build data list for the standard 2PP model (no testlets)
montar_data_list_2pp <- function(Y, hypers = list(sigma_a = 1, sigma_b = 1)) {
  stopifnot(is.matrix(Y))
  storage.mode(Y) <- "integer"
  N <- nrow(Y)
  I <- ncol(Y)
  
  stopifnot(all(c("sigma_a", "sigma_b") %in% names(hypers)))
  
  list(
    N = as.integer(N),
    vI = as.integer(I),
    Y = Y,
    sigma_a = as.numeric(hypers$sigma_a),
    sigma_b = as.numeric(hypers$sigma_b)
  )
}

## -----------------------------
## 5) Scenarios (K, nk, rho_start, rho_end)
##    [0.10, 0.40] weak; [0.40, 0.70] moderate; [0.70, 0.90] strong
## -----------------------------
Scenario <- rbind(
  c(4, 2, .9, .7),  # strong
  c(5, 2, .7, .4),  # moderate
  c(6, 2, 0,   0),  # null
  c(4, 4, .7, .4),  # moderate
  c(5, 4, 0,   0),  # null
  c(6, 4, .9, .7),  # strong
  c(4, 5, 0,   0),  # null
  c(5, 5, .9, .7),  # strong
  c(6, 5, .7, .4)   # moderate
)

## -----------------------------
## 6) Output helper
## -----------------------------
## Save a tidy posterior summary (mean, sd, quantiles, etc.) as an RDS file.
save_fit_summary <- function(fit, out_rds, pars = c("a", "b", "theta")) {
  sum_list <- summary(fit, pars = pars, probs = c(0.025, 0.5, 0.975))$summary
  sum_df <- data.frame(parameter = rownames(sum_list), sum_list, check.names = FALSE)
  rownames(sum_df) <- NULL
  saveRDS(sum_df, out_rds)
}

## -----------------------------
## 7) Sampler controls
## -----------------------------
nChains       <- 1
burnInSteps   <- 3000
thinSteps     <- 1
numSavedSteps <- 1000
nIter         <- ceiling(burnInSteps + (numSavedSteps * thinSteps))

## -----------------------------
## 8) Scenario and replication loop
## -----------------------------
R <- 20  # number of replications per scenario

for (j in seq_len(nrow(Scenario))) {
  
  ## Scenario-specific testlet layout
  K  <- Scenario[j, 1]
  nk <- rep(Scenario[j, 2], K)
  dk <- seq(1L, by = unique(nk), length.out = K)
  
  ## Scenario "true" rho sequence (used only for tagging filenames here)
  if (Scenario[j, 3] == 0 && Scenario[j, 4] == 0) {
    rho_verd <- rep(0, nk[1] - 1)
  } else {
    rho_verd <- round(seq(Scenario[j, 3], to = Scenario[j, 4], length.out = nk[1] - 1), 2)
  }
  tag <- if (length(rho_verd) == 0 || all(rho_verd == 0)) "rho0" else sprintf("%.2f", rho_verd[1])
  
  message("\n================= Scenario ", j, " =================")
  message("K = ", K, " | nk = ", nk[1], " | tag = ", tag)
  
  for (r in seq_len(R)) {
    
    ## Expected simulation file name (must match the generator script)
    sim_file <- sprintf("SimTestlet2PP_%d_%d_%s_%s_%d.rds", K, nk[1], tag, GEN_TYPE, r)
    sim_path <- file.path(pathSim, sim_file)
    
    ## Skip if simulation file is missing
    if (!file.exists(sim_path)) {
      message("Simulation file not found, skipping: ", sim_file)
      next
    }
    
    ## Load simulated dataset
    dat <- readRDS(sim_path)
    
    ## Expect response matrix in dat$Y
    if (is.null(dat$Y)) {
      message("dat$Y missing in: ", sim_file, " — skipping.")
      next
    }
    
    Y <- dat$Y
    if (!is.matrix(Y)) Y <- as.matrix(Y)
    storage.mode(Y) <- "integer"
    N <- nrow(Y)
    I <- ncol(Y)
    
    ## HU/HT rule: HU if nk==2, otherwise HT (simple heuristic)
    structure_vec <- rep(ifelse(nk[1] == 2, "HU", "HT"), K)
    
    ## Hyperparameters (selected prior level)
    hyper <- get_hypers(prior_level)
    
    ## Shared initial values:
    ##   - theta initialized from standardized raw scores
    ##   - a and b initialized near 0.1
    scores <- rowSums(Y)
    sc <- sd(scores)
    if (is.na(sc) || sc == 0) sc <- 1
    scores <- (scores - mean(scores)) / sc
    
    ## ===================== FIT TESTLET MODEL (optional) =======================
    if (model_to_fit %in% c("testlet", "both")) {
      
      ## Build Stan data list for testlet model
      data_list_testlet <- montar_data_list_testlet(
        Y = Y, K = K, nk = nk, dk = dk,
        structure = structure_vec,
        hypers = hyper
      )
      
      ## Total number of rho parameters across all testlets
      nRho <- sum(data_list_testlet$rho_len)
      
      ## Initial values for testlet model
      ini_testlet <- function() {
        list(
          theta      = as.numeric(scores),
          a          = rep(0.1, I),
          b          = rep(0.1, I),
          rho_global = if (nRho > 0) rep(0.1, nRho) else numeric(0)
        )
      }
      
      ## Parameters to monitor (rho included)
      monitor_pars_testlet <- c("theta", "a", "b", "rho_global")
      
      ## Fit with tryCatch to avoid aborting the whole loop on one failure
      fit_testlet <- tryCatch(
        {
          sampling(
            object  = sm_testlet,
            data    = data_list_testlet,
            chains  = nChains,
            iter    = nIter,
            thin    = thinSteps,
            warmup  = burnInSteps,
            pars    = monitor_pars_testlet,
            control = list(
              adapt_delta    = nutsc$adapt_delta,
              max_treedepth  = nutsc$max_treedepth
            ),
            init    = ini_testlet,
            refresh = 200
          )
        },
        error = function(e) {
          message("Sampling failed (Testlet) for ", sim_file, ": ", conditionMessage(e))
          return(NULL)
        }
      )
      
      ## Save summary if fitting succeeded
      if (!is.null(fit_testlet)) {
        summary_file <- file.path(
          pathFit,
          sprintf("Summary_Testlet_%d_%d_%s_%d.rds", K, nk[1], tag, r)
        )
        save_fit_summary(fit_testlet, summary_file, pars = monitor_pars_testlet)
        message("Saved (Testlet): ", basename(summary_file))
      }
    }
    
    ## ===================== FIT STANDARD 2PP (optional) ========================
    if (model_to_fit %in% c("2pp", "both")) {
      
      ## Build Stan data list for 2PP
      data_list_2pp <- montar_data_list_2pp(
        Y = Y,
        hypers = list(sigma_a = hyper$sigma_a, sigma_b = hyper$sigma_b)
      )
      
      ## Initial values for 2PP
      ini_2pp <- function() {
        list(
          theta = as.numeric(scores),
          a     = rep(0.1, data_list_2pp$vI),
          b     = rep(0.1, data_list_2pp$vI)
        )
      }
      
      ## Parameters to monitor
      monitor_pars_2pp <- c("theta", "a", "b")
      
      ## Fit with tryCatch
      fit_2pp <- tryCatch(
        {
          sampling(
            object  = sm_2pp,
            data    = data_list_2pp,
            chains  = nChains,
            iter    = nIter,
            thin    = thinSteps,
            warmup  = burnInSteps,
            pars    = monitor_pars_2pp,
            control = list(
              adapt_delta   = nutsc$adapt_delta,
              max_treedepth = nutsc$max_treedepth
            ),
            init    = ini_2pp,
            refresh = 200
          )
        },
        error = function(e) {
          message("Sampling failed (2PP) for ", sim_file, ": ", conditionMessage(e))
          return(NULL)
        }
      )
      
      ## Save summary if fitting succeeded
      if (!is.null(fit_2pp)) {
        summary_file_2pp <- file.path(
          pathFit,
          sprintf("Summary_2PP_%d_%d_%s_%d.rds", K, nk[1], tag, r)
        )
        save_fit_summary(fit_2pp, summary_file_2pp, pars = monitor_pars_2pp)
        message("Saved (2PP): ", basename(summary_file_2pp))
      }
    }
    
  } # end replication loop
} # end scenario loop

message("\n===== Finished all runs =====")
