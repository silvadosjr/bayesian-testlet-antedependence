## -------------------- 0) Paths and helper sources ----------------------------
## This script fits two Bayesian probit IRT models to a large-scale assessment dataset:
##   (1) Testlet 2PP (accounts for local item dependence within predefined testlets)
##   (2) Standard 2PP (assumes local independence across all items)
##
## It then produces posterior summaries and comparison plots for:
##   - Discrimination parameters (a)
##   - Difficulty parameters (b)
##   - Testlet correlation parameters (rho_global)
##
## Repository paths (edit if your local clone differs)
root_local    <- "~/GitHub/bayesian-testlet-antedependence"
path_project  <- "~/GitHub/bayesian-testlet-antedependence/Real Data Analysis/Large-Scale Educational Assessment"
stopifnot(dir.exists(root_local), dir.exists(path_project))

setwd(path_project)

## ------------------------- 1) Setup and packages -----------------------------
## Core: rstan (model fitting) + tidyverse stack (data wrangling + ggplot)
## Extra: mvtnorm/tmvtnorm, loo, shinystan are loaded because helper functions
## and/or downstream diagnostics may rely on them.
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

## Stan configuration
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

## Path to Stan programs (shared across analyses)
pathProgram <- here(root_local, "Programs")

## Source helper functions used for indexing + diagnostics (Q3, envelopes, heatmaps, etc.)
## This file is expected to be part of the repository.
source(here("Real Data Analysis", "HelpersRealDataAnal.R"))

## Output directories (model fits and figures)
pathFit <- here(path_project, "Fit")
dir.create(pathFit, showWarnings = FALSE, recursive = TRUE)
stopifnot(dir.exists(pathFit))

saveFigures <- here(pathFit, "Figures")
dir.create(saveFigures, showWarnings = FALSE, recursive = TRUE)
stopifnot(dir.exists(saveFigures))

## ------------------- 2) Load data (3rd year high school) ---------------------
## Binary response matrix with 0/1 entries (and possibly already NA-handled).
file_full <- here(path_project, "Binary_1EM_municipio_n25k.rds")
stopifnot(file.exists(file_full))

mYcNA <- readRDS(file_full)

## Basic checks
n  <- nrow(mYcNA)  # number of students
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

## Consistency check: testlet items should match the complement of ind_items
all_items      <- seq_len(vI)
testlet_items  <- setdiff(all_items, ind_items)
stopifnot(length(testlet_items) == sum(nk))

## Partition items by testlet in the declared nk order
idx_testlets <- split(sort(testlet_items), rep(seq_along(nk), times = nk))
stopifnot(sum(lengths(idx_testlets)) == length(testlet_items),
          length(idx_testlets) == K)

## --------------------- 4) Stan data lists ------------------------------------
## (a) Testlet-2PP (probit), with local dependence inside testlets
##     S_mc is used in Stan generated quantities (Monte Carlo approximation)
data_testlet <- list(
  I = vI, N = n, K = K, dk = dk, nk = nk,
  ind_items = ind_items, n_ind = vI - sum(nk),
  Y = mYcNA,
  sigma_a = .6, sigma_b = 4, sigma_rho = 1,
  rho_len = idx$rho_len, S_mc = 200,
  rho_start = idx$rho_start
)

## (b) Standard 2PP (probit) with local independence
data_2pp <- list(
  N = n, vI = vI, Y = mYcNA,
  sigma_a = .6, sigma_b = 4
)

## ---------------------- 5) Initial values -----------------------------------
## Use standardized raw scores as an informative initialization for theta.
scores <- scale(rowSums(mYcNA))[, 1]

## Initial values for Testlet-2PP (includes rho_global)
init_testlet <- function() {
  list(
    theta      = as.numeric(scores),
    a          = rep(0.1, vI),
    b          = rep(0.1, vI),
    rho_global = rep(0.1, sum(idx$rho_len))
  )
}

## Initial values for 2PP (no rho parameters)
init_2pp <- function() {
  list(
    theta = as.numeric(scores),
    a     = rep(0.1, vI),
    b     = rep(0.1, vI)
  )
}

## Parameters monitored (kept in the fitted object)
pars_testlet <- c("a", "b", "theta", "rho_global")
pars_2pp     <- c("a", "b", "theta")

## NUTS settings (single chain configuration for large N)
nChains       <- 1
burnInSteps   <- 1000
thinSteps     <- 1
numSavedSteps <- 1000
nIter         <- ceiling(burnInSteps + numSavedSteps * thinSteps)

ctrl_nuts <- list(adapt_delta = 0.8, max_treedepth = 10)

## --------------------------- 6) Fit the models -------------------------------
## (a) Testlet-2PP
t0 <- tic("Fitting Testlet-2PP")
fit_testlet <- stan(
  data   = data_testlet,
  file   = file.path(pathProgram, "Testlet2PP_HTGeral.stan"),
  init   = init_testlet,
  chains = nChains, pars = pars_testlet,
  iter   = nIter, warmup = burnInSteps, thin = thinSteps,
  control = ctrl_nuts
)
toc(t0)

## (b) 2PP
t0 <- tic("Fitting 2PP")
fit_2pp <- stan(
  data   = data_2pp,
  file   = file.path(pathProgram, "2PPModelVec.stan"),
  init   = init_2pp,
  chains = nChains, pars = pars_2pp,
  iter   = nIter, warmup = burnInSteps, thin = thinSteps,
  control = list(adapt_delta = 0.8)  # depth default ok
)
toc(t0)

## Save fitted objects for reuse (avoid refitting)
saveRDS(fit_testlet, file.path(pathFit, "ResultTestlet2PP_1EM_municipio25k.rds"))
saveRDS(fit_2pp,     file.path(pathFit, "Result2PP_1EM_municipio25k.rds"))


## Reload fits (useful when running only the plotting section)
fit_testlet<-readRDS(file.path(pathFit,'ResultTestlet2PP_1EM_municipio25k.rds'))
fit_2pp<-readRDS(file.path(pathFit,"Result2PP_1EM_municipio25k.rds"))


## ======================= Posterior summaries and plots =======================

## Helper: extract rstan summary as a data.frame with an explicit Parameter column
post_summary <- function(fit, pars) {
  ss <- rstan::summary(fit, pars = pars)$summary
  out <- as.data.frame(ss)
  out$Parameter <- rownames(ss)
  rownames(out) <- NULL
  out
}

## Posterior summaries: discrimination (a), difficulty (b), latent traits (theta),
## and testlet correlations (rho_global, Testlet model only).
sum_a1     <- post_summary(fit_testlet, "a")
sum_b1     <- post_summary(fit_testlet, "b")
sum_theta1 <- post_summary(fit_testlet, "theta")
sum_rho    <- post_summary(fit_testlet, "rho_global")

sum_a2     <- post_summary(fit_2pp, "a")
sum_b2     <- post_summary(fit_2pp, "b")
sum_theta2 <- post_summary(fit_2pp, "theta")

## ---- Robust ordering of parameter rows by index ----
## rstan returns parameter names like "a[13]" and "theta[200]".
## These utilities parse the integer index and sort accordingly.
get_param_index <- function(param_vec, par) {
  mat <- stringr::str_match(param_vec, sprintf("^%s\\[(\\d+)\\]$", par))
  suppressWarnings(as.integer(mat[, 2]))
}

order_param_df <- function(df, par) {
  idx <- get_param_index(df$Parameter, par)
  keep <- !is.na(idx)
  df2 <- df[keep, , drop = FALSE]
  idx2 <- idx[keep]
  df2[order(idx2), , drop = FALSE]
}

## Apply ordering for a, b, and theta (both models)
sum_a1     <- order_param_df(sum_a1, "a")
sum_a2     <- order_param_df(sum_a2, "a")
sum_b1     <- order_param_df(sum_b1, "b")
sum_b2     <- order_param_df(sum_b2, "b")
sum_theta1 <- order_param_df(sum_theta1, "theta")
sum_theta2 <- order_param_df(sum_theta2, "theta")

## ---- Labels for rho parameters ----
## We label each rho as T_k(rho[r]) where:
##   k = testlet index
##   r = lag/correlation index within that testlet (1..rho_len[k])
build_rho_labels <- function(rho_len) {
  exprs <- list()
  for (k in seq_along(rho_len)) {
    for (r in seq_len(rho_len[k])) {
      exprs[[length(exprs) + 1]] <- bquote(T[.(k)](rho[.(r)]))
    }
  }
  as.expression(exprs)
}

rho_labels <- build_rho_labels(idx$rho_len)

rho_data <- data.frame(
  Rho      = seq_len(sum(idx$rho_len)),
  Estimate = sum_rho$mean,
  CI_Low   = sum_rho$`2.5%`,
  CI_High  = sum_rho$`97.5%`
)

## ---- Item labels with testlet annotation ----
## Map each item i to its testlet k (or NA if independent).
testlet_of_item <- rep(NA_integer_, vI)
for (k in seq_along(idx_testlets)) {
  testlet_of_item[idx_testlets[[k]]] <- k
}

## Build expression labels:
##   - independent items: just "i"
##   - testlet items:     "T_k(i)"
build_item_labels <- function(vI, testlet_of_item) {
  labs <- vector("list", vI)
  for (i in seq_len(vI)) {
    k <- testlet_of_item[i]
    labs[[i]] <- if (is.na(k)) {
      bquote(.(i))
    } else {
      bquote(T[.(k)](.(i)))
    }
  }
  as.expression(labs)
}

item_labels <- build_item_labels(vI, testlet_of_item)

## ---- Plot palette (grayscale) ----
## Here the Testlet 2PP is black and the 2PP is a lighter gray.
scale_cols <- scale_color_manual(values = c("Testlet 2PP" = "black", "2PP" = "grey50"))

## -----------------------------
## Plot 1: rho estimates (Testlet model only)
## -----------------------------
p_rho <- ggplot(rho_data, aes(x = Rho, y = Estimate)) +
  geom_point(size = 2, color = "black") +
  geom_errorbar(aes(ymin = CI_Low, ymax = CI_High), width = 0.15, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal(base_size = 25) +
  scale_x_continuous(breaks = rho_data$Rho, labels = rho_labels) +
  labs(x = NULL, y = "Estimate")

ggsave(file.path(pathFit, "Fig_Rho_Testlet.pdf"), p_rho, width = 9, height = 3.6)

## -----------------------------
## Plot 2: discrimination (a) comparison
## -----------------------------
## Creates a long data.frame with one row per (item, model) containing:
##   estimate + 95% credible interval.
disc_df <- data.frame(
  Item     = rep(seq_len(vI), 2),
  Estimate = c(sum_a1$mean, sum_a2$mean),
  CI_Low   = c(sum_a1$`2.5%`, sum_a2$`2.5%`),
  CI_High  = c(sum_a1$`97.5%`, sum_a2$`97.5%`),
  Model    = rep(c("Testlet 2PP", "2PP"), each = vI)
)

p_a <- ggplot(disc_df, aes(Item, Estimate, color = Model)) +
  geom_point(size = 1.7) +
  geom_errorbar(aes(ymin = CI_Low, ymax = CI_High), width = 0.15) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_cols +
  theme_minimal(base_size = 14) +
  scale_x_continuous(breaks = 1:vI, labels = item_labels) +
  labs(x = "Item", y = "Estimate", color = "Model") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(pathFit, "Fig_Discrimination_gray.pdf"), p_a, width = 11, height = 4.4)

## -----------------------------
## Plot 3: difficulty (b) comparison
## -----------------------------
diff_df <- data.frame(
  Item     = rep(seq_len(vI), 2),
  Estimate = c(sum_b1$mean, sum_b2$mean),
  CI_Low   = c(sum_b1$`2.5%`, sum_b2$`2.5%`),
  CI_High  = c(sum_b1$`97.5%`, sum_b2$`97.5%`),
  Model    = rep(c("Testlet 2PP", "2PP"), each = vI)
)

p_b <- ggplot(diff_df, aes(Item, Estimate, color = Model)) +
  geom_point(size = 1.7) +
  geom_errorbar(aes(ymin = CI_Low, ymax = CI_High), width = 0.15) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_cols +
  theme_minimal(base_size = 14) +
  scale_x_continuous(breaks = 1:vI, labels = item_labels) +
  labs(x = "Item", y = "Estimate", color = "Model") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(pathFit, "Fig_Difficulty_gray.pdf"), p_b, width = 11, height = 4.4)


