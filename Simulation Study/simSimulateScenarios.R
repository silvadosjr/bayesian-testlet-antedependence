## =============================================================================
## Simulation Study: 2PP with residual dependence (NON-Toeplitz within testlets)
## =============================================================================
## This script generates binary item-response data under a multivariate probit 2PP
## model with *residual* dependence within testlets. The dependence can be:
##   - "TOEP": Toeplitz-by-lag within each testlet (correlation depends on lag)
##   - "UNS" : Unstructured (non-Toeplitz) correlation within each testlet
##   - "MIX" : Alternates by testlet (e.g., odd testlets unstructured, even Toeplitz)
##
## Output: RDS simulation files created via rstansim::simulate_data(), containing
## simulated responses and (optionally) latent objects depending on vars_keep.
##
## Requirements (project structure):
##   - Programs/SimMult2PP_Chol.stan
##   - Simulation Study/Itens30.txt
##   - Simulation Study/simData/ (created if missing)
##
## Notes:
##   - Dependence is introduced through a correlation matrix Sigma_Y, then passed
##     to Stan through its Cholesky factor (L_Sigma).
##   - For GEN_TYPE = "UNS", PD is enforced via spectral truncation (no nearPD).
## =============================================================================

setwd('~/GitHub/bayesian-testlet-antedependence')

library(here)
library(rstan)
library(rstansim)
library(Matrix)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## ------------------------------ Utilities ------------------------------------

## Build a Toeplitz correlation matrix from lag correlations rhos:
##   - rhos has length L-1
##   - output is L×L with diag=1 and (i,j) correlation = rhos[|i-j|]
make_toeplitz <- function(rhos) {
  L <- length(rhos) + 1L
  C <- diag(L)
  for (r in seq_len(L - 1L)) {
    C[row(C) == col(C) - r] <- rhos[r]
    C[row(C) == col(C) + r] <- rhos[r]
  }
  C
}

## Construct a "well-behaved" unstructured correlation matrix (m×m) whose
## average off-diagonal magnitude is approximately target_rho.
##
## Strategy:
##   1) Build a random SPD matrix via S = A'A and convert to correlation R0.
##   2) Break structured patterns via a random permutation P and add small jitter.
##   3) Rescale off-diagonal entries to match a target mean correlation.
##   4) Enforce positive definiteness via spectral projection (truncate eigenvalues).
##   5) Re-normalize to correlation scale (diag=1).
##
## Args:
##   m          : block size
##   target_rho : desired average off-diagonal correlation (approximately)
##   mix        : shrink factor for the off-diagonal part after rescaling (0..1)
##   jitter     : small symmetric noise added to off-diagonals (pattern breaking)
##   eps        : minimum eigenvalue after truncation (PD safeguard)
make_corr_unstructured <- function(m, target_rho = 0.4, mix = 0.8, jitter = 0.05, eps = 1e-6) {
  stopifnot(m >= 1, mix >= 0, mix <= 1)
  if (m == 1) return(matrix(1, 1, 1))
  
  ## 1) Random SPD base, then correlation
  A  <- matrix(rnorm(m * m), m, m)
  S  <- crossprod(A)          # SPD
  R0 <- cov2cor(S)            # correlation, PD
  
  ## 2) Break patterns: random permutation + small jitter
  P  <- diag(m)
  if (m >= 3) P <- P[sample.int(m), , drop = FALSE]
  R  <- t(P) %*% R0 %*% P
  
  if (jitter > 0) {
    J  <- matrix(runif(m * m, -1, 1), m)
    diag(J) <- 0
    R  <- R + jitter * J
  }
  
  ## 3) Rescale off-diagonals to match target_rho (on average)
  R  <- (R + t(R)) / 2
  diag(R) <- 1
  off <- R - diag(m)
  
  mu <- mean(off[upper.tri(off)])
  sc <- if (!is.na(mu) && abs(mu) > .Machine$double.eps) target_rho / mu else 1
  R1 <- diag(m) + mix * (off * sc)
  
  ## 4) Spectral projection to PD (truncate small/negative eigenvalues)
  R1  <- (R1 + t(R1)) / 2
  eig <- eigen(R1, symmetric = TRUE)
  lam <- pmax(eig$values, eps)
  Rpd <- eig$vectors %*% diag(lam) %*% t(eig$vectors)
  
  ## 5) Normalize back to correlation (diag=1)
  D  <- sqrt(diag(Rpd))
  R2 <- Rpd / outer(D, D)
  diag(R2) <- 1
  
  ## 6) Final symmetrization for numerical stability
  R2 <- (R2 + t(R2)) / 2
  R2
}

## -------------------------------- Paths --------------------------------------
programas_dir <- here("Programs")
pathSim       <- here("Simulation Study", "simData")
dir.create(pathSim, showWarnings = FALSE, recursive = TRUE)

## --------------------------- Fixed dimensions --------------------------------
vI <- 30   # total number of items
n  <- 350  # sample size (respondents)

## ---------------------- Item parameters (a, b) --------------------------------
## Reads a and b values from file and standardizes b to have mean 0 and sd 1.
stopifnot(file.exists(here("Simulation Study", "Itens30.txt")))
zeta_verd <- read.table(here("Simulation Study", "Itens30.txt"))[1:vI, -3]
zeta_verd[, 2] <- round((zeta_verd[, 2] - mean(zeta_verd[, 2])) / sd(zeta_verd[, 2]), 2)

## ---------------------- Scenarios (K, nk, rho_start, rho_end) -----------------
## rho_start -> rho_end defines a sequence of lag correlations (Toeplitz) or
## a target average correlation (unstructured); "nula" means rho=0.
Scenario <- rbind(
  c(4, 2, .9, .7),   # strong
  c(5, 2, .7, .4),   # moderate
  c(6, 2, 0, 0),     # null
  c(4, 4, .7, .4),   # moderate
  c(5, 4, 0, 0),     # null
  c(6, 4, .9, .7),   # strong
  c(4, 5, 0, 0),     # null
  c(5, 5, .9, .7),   # strong
  c(6, 5, .7, .4)    # moderate
)

## ------------------------------ Replications ---------------------------------
R <- 20  # number of simulated datasets per scenario

## ------------------------- Dependence generation config -----------------------
## GEN_TYPE:
##   "TOEP" -> Toeplitz within each testlet
##   "UNS"  -> unstructured within each testlet
##   "MIX"  -> alternate by testlet (odd: UNS, even: TOEP)
GEN_TYPE <- "TOEP"

## Parameters used by make_corr_unstructured() when GEN_TYPE is UNS or MIX:
MIX_INT  <- 0.85  # shrinkage applied to the rescaled off-diagonals
JITTER   <- 0.07  # additional perturbation (breaks patterns)

## Auxiliary data (fixed) passed to simulate_data()
data_simAux <- list(N = n, vI = vI)

## ---------------------------- Scenario loop ----------------------------------
## set.seed(123)  # optional: global reproducibility

for (j in 1:nrow(Scenario)) {
  
  ## Scenario-specific testlet structure
  K  <- Scenario[j, 1]
  nk <- rep(Scenario[j, 2], K)                 # constant block size per testlet
  dk <- seq(1, K * unique(nk), by = unique(nk))# starting indices (1-based)
  
  ## Define "true" lag correlations for Toeplitz blocks:
  ## if nk==1, rho_verd has length 0 (no within-block dependence).
  rho_verd <- round(seq(Scenario[j, 3], to = Scenario[j, 4], length.out = nk[1] - 1), 2)
  
  ## Tag used in filenames (rho0 or first rho value)
  tag_rho <- if (length(rho_verd) == 0 || all(rho_verd == 0)) "rho0" else sprintf("%.2f", rho_verd[1])
  
  ## Output dataset name used by rstansim (used to build file names)
  data_name <- paste0("SimTestlet2PP_", K, "_", unique(nk), "_", tag_rho, "_", GEN_TYPE)
  
  ## =================== Build Sigma_Y (residual correlation) ===================
  ## Sigma_Y is an item-by-item correlation matrix. It is block diagonal:
  ## each testlet defines a within-block correlation, and between-testlet
  ## correlations remain 0 (identity outside blocks).
  Sigma_Y <- diag(vI)
  
  ## For unstructured blocks we target the average off-diagonal correlation
  target_rho <- mean(rho_verd)
  
  for (i in 1:K) {
    idx <- dk[i]:(dk[i] + nk[i] - 1)  # item indices in testlet i
    
    if (GEN_TYPE == "TOEP") {
      ## Toeplitz-by-lag correlation within testlet
      Cblk <- make_toeplitz(rho_verd)
      
    } else if (GEN_TYPE == "MIX") {
      ## Alternate: odd testlets unstructured, even testlets Toeplitz
      if (i %% 2 == 1) {
        Cblk <- make_corr_unstructured(nk[i], target_rho, MIX_INT, JITTER)
      } else {
        Cblk <- make_toeplitz(rho_verd)
      }
      
    } else {
      ## Default: unstructured correlation within testlet
      Cblk <- make_corr_unstructured(nk[i], target_rho, MIX_INT, JITTER)
    }
    
    ## Insert block into the global correlation matrix
    Sigma_Y[idx, idx] <- Cblk
  }
  
  ## ============================= Simulate data ================================
  ## The Stan generator SimMult2PP_Chol.stan is expected to:
  ##  - take N and vI as input,
  ##  - use sim_a, sim_b, and L_Sigma (Cholesky of Sigma_Y) as parameters,
  ##  - output sim_Y (and optionally sim_theta / L_Sigma depending on vars_keep).
  param_list <- list(
    sim_a   = as.numeric(zeta_verd[, 1]),
    sim_b   = as.numeric(zeta_verd[, 2]),
    L_Sigma = t(chol(Sigma_Y))   # upper-triangular or lower depends on Stan file convention
  )
  
  file_stan <- file.path(programas_dir, "SimMult2PP_Chol.stan")
  
  ## Variables to keep from the simulator output
  vars_keep <- c("sim_Y", "sim_theta", "L_Sigma")
  
  ## set.seed(123 + j)  # optional: scenario-level reproducibility
  
  simulate_data(
    file         = file_stan,
    data_name    = data_name,
    input_data   = list(N = n, vI = vI),
    param_values = param_list,
    vars         = vars_keep,
    nsim         = R,
    path         = pathSim
  )
  
  message(sprintf(
    "Scenario %d concluded: K=%d, nk=%d, rho=[%.2f -> %.2f], GEN=%s",
    j, K, nk[1], Scenario[j, 3], Scenario[j, 4], GEN_TYPE
  ))
}

  
# dat<-readRDS(file.path(pathSim,'SimTestlet2PP_4_4_0.70_UNS_1.rds'))  
#   
# L_Sigma<-dat$L_Sigma  
# 
# 
# View(L_Sigma%*%t(L_Sigma))




  

