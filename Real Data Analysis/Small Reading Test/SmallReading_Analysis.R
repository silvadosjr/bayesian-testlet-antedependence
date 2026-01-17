## =============================================================================
## Real Data Analysis — Small Reading Test
## Bayesian Testlet-2PP (probit) with local dependence vs. standard 2PP (probit)
##
## This script:
##  1) Loads the "Small Reading test" binary responses (0/1 with possible NAs),
##  2) Defines a testlet structure (K=3, each with nk=2 items),
##  3) Fits:
##     - A Testlet-2PP probit model with (HU/HT-style) local dependence parameters,
##     - A standard 2PP probit model (local independence),
##  4) Saves fitted Stan objects,
##  5) Produces posterior summaries and comparison plots for:
##     - rho parameters (testlet dependence),
##     - item discrimination (a),
##     - item difficulty (b).
##
## Requirements:
##  - The project directory must exist and contain:
##      Programs/Testlet2PP_HTGeral.stan
##      Programs/2PPModelVec.stan
##      Real Data Analysis/HelpersRealDataAnal.R
##      Real Data Analysis/ld_Small_reading_test.txt
##  - R packages listed below.
##
## Notes:
##  - This script uses 1 chain and thinning for speed; for final inference, increase
##    chains and consider stronger NUTS settings.
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

## Simple timing utilities for logging
tic <- function(msg) { cat(sprintf("\n[START] %s ... %s\n", msg, Sys.time())); Sys.time() }
toc <- function(t0)  { cat(sprintf("[ END ] Elapsed: %s\n", Sys.time() - t0)) }

## Path to Stan programs
pathProgram <- here(root_local, "Programs")

## Source helper functions used below (e.g., make_rho_index)
source(here(path_project, "HelpersRealDataAnal.R"))

## Output directories
pathFit <- here(path_project, "Fit")
dir.create(pathFit, showWarnings = FALSE, recursive = TRUE)
stopifnot(dir.exists(pathFit))

saveFigures <- here(pathFit, "Figures")
dir.create(saveFigures, showWarnings = FALSE, recursive = TRUE)
stopifnot(dir.exists(saveFigures))


## ------------------- 2) Load data — Small Reading test -----------------------
## (a) Full dataset (used here as "training" set)
file_full <- file.path(path_project, "ld_Small_reading_test.txt")
stopifnot(file.exists(file_full))
base_full <- read.table(file_full)

## Response matrix:
##  - remove first column (typically an ID),
##  - keep 0/1 with NAs for missing.
mYc   <- base_full[, -1]
mYcNA <- as.data.frame(mYc)

## Replace missing with 0 for Stan input (common in these pipelines).
## IMPORTANT: confirm that missingness treatment is appropriate for your study.
mYcNA[is.na(mYc)] <- 0L

## Basic checks
n  <- nrow(mYcNA)  # number of respondents
vI <- ncol(mYcNA)  # number of items
stopifnot(n > 0, vI > 0)

## Reorder items to match the intended testlet ordering
mYcNA <- mYcNA[, c(1, 6, 4, 5, 2, 3)]


## -------------------- 3) Testlet structure and indices -----------------------
## Declared structure:
##  - K = number of testlets
##  - nk = number of items in each testlet
##  - dk = starting positions (if required by the Stan program)
##  - is_HU = logical flags controlling the HU/HT structure used by helper code
K  <- 3
nk <- c(2, 2, 2)
dk <- c(1, 3, 5)
is_HU <- c(TRUE, TRUE, TRUE)  # user-defined pattern; kept as-is

## Build indices used to map rho parameters into testlets
idx <- make_rho_index(nk, is_HU)  # provides $rho_len, $rho_start, etc.

## Independent items (items not belonging to any testlet)
ind_items <- integer(0L)

## Consistency checks: testlet items are the complement of independent items
all_items     <- seq_len(vI)
testlet_items <- setdiff(all_items, ind_items)
stopifnot(length(testlet_items) == sum(nk))

## Partition items by testlet according to nk
idx_testlets <- split(sort(testlet_items), rep(seq_along(nk), times = nk))
stopifnot(
  sum(lengths(idx_testlets)) == length(testlet_items),
  length(idx_testlets) == K
)


## --------------------- 4) Data lists for Stan -------------------------------
## Testlet-2PP (probit) — "HT general" version (per Stan file naming)
data_testlet <- list(
  I = vI, N = n, K = K, dk = dk, nk = nk,
  ind_items = ind_items, n_ind = vI - sum(nk),
  Y = mYcNA,
  sigma_a = .6, sigma_b = 4, sigma_rho = 1,
  rho_len = idx$rho_len, S_mc = 200,
  rho_start = idx$rho_start
)

## Standard 2PP (probit) under local independence
data_2pp <- list(
  N = n, vI = vI, Y = mYcNA,
  sigma_a = .6, sigma_b = 4
)


## ---------------------- 5) Initial values -----------------------------------
## Standardized raw score as an initializer for theta
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

## (ii) Standard 2PP does not include rho_global
init_2pp <- function() {
  list(
    theta = as.numeric(scores),
    a     = rep(0.1, vI),
    b     = rep(0.1, vI)
  )
}

## Parameters to monitor
pars_testlet <- c("a", "b", "theta", "rho_global")
pars_2pp     <- c("a", "b", "theta")

## NUTS settings (kept modest for speed)
nChains       <- 1
burnInSteps   <- 1000
thinSteps     <- 20
numSavedSteps <- 1000
nIter         <- ceiling(burnInSteps + numSavedSteps * thinSteps)
ctrl_nuts     <- list(adapt_delta = 0.8, max_treedepth = 10)


## --------------------------- 6) Fit models ----------------------------------
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

## (b) Standard 2PP
t0 <- tic("Fitting 2PP")
fit_2pp <- stan(
  data   = data_2pp,
  file   = file.path(pathProgram, "2PPModelVec.stan"),
  init   = init_2pp,
  chains = nChains, pars = pars_2pp,
  iter   = nIter, warmup = burnInSteps, thin = thinSteps,
  control = list(adapt_delta = 0.8)
)
toc(t0)

## Save fitted objects for reproducibility
saveRDS(fit_testlet, file.path(pathFit, "ResultTestlet2PP_SmallReading.rds"))
saveRDS(fit_2pp,     file.path(pathFit, "Result2PP_SmallReading.rds"))

## Reload (useful when running only the plotting block)
fit_testlet <- readRDS(file.path(pathFit, "ResultTestlet2PP_SmallReading.rds"))
fit_2pp     <- readRDS(file.path(pathFit, "Result2PP_SmallReading.rds"))


## ======================= Posterior summaries and plots =======================
## Helper to extract posterior summaries from rstan::summary()
post_summary <- function(fit, pars) {
  ss <- rstan::summary(fit, pars = pars)$summary
  out <- as.data.frame(ss)
  out$Parameter <- rownames(ss)
  rownames(out) <- NULL
  out
}

## Summaries: Testlet model
sum_a1     <- post_summary(fit_testlet, "a")
sum_b1     <- post_summary(fit_testlet, "b")
sum_theta1 <- post_summary(fit_testlet, "theta")
sum_rho    <- post_summary(fit_testlet, "rho_global")

## Summaries: Standard 2PP
sum_a2     <- post_summary(fit_2pp, "a")
sum_b2     <- post_summary(fit_2pp, "b")
sum_theta2 <- post_summary(fit_2pp, "theta")

## Extract numeric indices from parameter names like "a[13]" or "theta[200]"
get_param_index <- function(param_vec, par) {
  mat <- stringr::str_match(param_vec, sprintf("^%s\\[(\\d+)\\]$", par))
  idx <- suppressWarnings(as.integer(mat[, 2]))  # NA for non-matching rows
  idx
}

## Keep only valid entries and sort by their numeric index
order_param_df <- function(df, par) {
  idx <- get_param_index(df$Parameter, par)
  keep <- !is.na(idx)
  df2 <- df[keep, , drop = FALSE]
  idx2 <- idx[keep]
  df2[order(idx2), , drop = FALSE]
}

## Apply ordering (important for consistent item alignment across models)
sum_a1     <- order_param_df(sum_a1,     "a")
sum_a2     <- order_param_df(sum_a2,     "a")
sum_b1     <- order_param_df(sum_b1,     "b")
sum_b2     <- order_param_df(sum_b2,     "b")
sum_theta1 <- order_param_df(sum_theta1, "theta")
sum_theta2 <- order_param_df(sum_theta2, "theta")


## ------------------ Build labels for rho parameters --------------------------
## Create expressions like T_k(rho[r]) where r runs over rho_len[k]
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
  Rho = seq_len(sum(idx$rho_len)),
  Estimate = sum_rho$mean,
  CI_Low   = sum_rho$`2.5%`,
  CI_High  = sum_rho$`97.5%`
)


## ------------------ Build labels for item axis (testlets) --------------------
## Map each item to its testlet (or NA if independent)
testlet_of_item <- rep(NA_integer_, vI)
for (k in seq_along(idx_testlets)) {
  testlet_of_item[idx_testlets[[k]]] <- k
}

## Item labels:
##  - independent items: "i"
##  - testlet items: T_k(i)
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


## -------------------------- Plot styling -------------------------------------
## Grayscale palette for model comparison
scale_cols <- scale_color_manual(values = c("Testlet 2PP" = "grey50", "2PP" = "black"))


## -------------------------- Plot: rho ----------------------------------------
p_rho <- ggplot(rho_data, aes(x = Rho, y = Estimate)) +
  geom_point(size = 2, color = "black") +
  geom_errorbar(aes(ymin = CI_Low, ymax = CI_High), width = 0.15, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal(base_size = 14) +
  scale_x_continuous(breaks = rho_data$Rho, labels = rho_labels) +
  labs(x = NULL, y = "Estimate")

ggsave(file.path(pathFit, "Fig_Rho_Testlet.pdf"), p_rho, width = 9, height = 3.6)


## ----------------------- Plot: discrimination (a) ----------------------------
disc_df <- data.frame(
  Item = rep(seq_len(vI), 2),
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


## ------------------------- Plot: difficulty (b) ------------------------------
diff_df <- data.frame(
  Item = rep(seq_len(vI), 2),
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





