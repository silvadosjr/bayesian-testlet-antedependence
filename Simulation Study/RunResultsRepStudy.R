## ======================================================================
## Comparative analysis: 2PP vs Testlet2PP
## Metrics: Bias, |Bias|, RMSE, ESS, Rhat (including rho parameters when available)
## ======================================================================
## This script:
##   1) Reads posterior summaries (saved as .rds) produced in the fitting step
##      for both models (2PP and Testlet2PP).
##   2) Aligns posterior estimates with the "true" parameters used in simulation:
##        - item parameters a_i and b_i (from Itens30.txt)
##        - person parameters theta_j (from the corresponding simulation .rds file)
##        - residual dependence parameters rho (scenario-dependent truth)
##   3) Computes per-draw/per-parameter error quantities (bias, abs bias, squared error),
##      and aggregates them by scenario and by index.
##   4) Saves tidy results (.rds) and produces comparison boxplots by scenario
##      for bias, |bias|, RMSE, ESS, and Rhat, with a grayscale palette.
##
## Expected directory structure:
##   Simulation Study/
##     simData/        -> simulation .rds files (contain Y and sim_theta)
##     fitData/
##       2PP/          -> posterior summary .rds files for 2PP
##       Testlet2PP/   -> posterior summary .rds files for Testlet2PP
##     RecoveryResults/ -> output tables and figures created by this script
##
## Important:
##   - Filenames are parsed to recover (K, nk, tag, rep) and model.
##   - GEN_TYPE must match the simulation file naming convention.
## ======================================================================

setwd('~/GitHub/bayesian-testlet-antedependence')

library(tidyverse)
library(stringr)
library(here)

programas_dir <- here("Programs")
pathSim       <- here("Simulation Study", "simData")
pathFit       <- here("Simulation Study", "fitData")
pathFit_2PP   <- file.path(pathFit, "2PP")
pathFit_Test  <- file.path(pathFit, "Testlet2PP")
pathResults   <- here("Simulation Study", "RecoveryResults")

if (!dir.exists(pathResults)) dir.create(pathResults, recursive = TRUE)

stopifnot(dir.exists(pathSim), dir.exists(pathFit_2PP), dir.exists(pathFit_Test))

## ----------------------------------------------------------------------
## True item parameters (a, b)
## ----------------------------------------------------------------------
## Reads Itens30.txt, keeps columns for discrimination (a) and raw difficulty (b_raw),
## then standardizes b to mean 0 and sd 1 to match the simulator convention.
itens_path <- file.path(here("Simulation Study", "Itens30.txt"))
stopifnot(file.exists(itens_path))

zeta_verd <- read.table(itens_path)[, 1:3]
zeta_verd <- zeta_verd[, -3, drop = FALSE]
colnames(zeta_verd) <- c("a_true", "b_raw")

zeta_verd <- zeta_verd %>%
  mutate(b_true = round((b_raw - mean(b_raw)) / sd(b_raw), 2)) %>%
  select(a_true, b_true)

## ----------------------------------------------------------------------
## Utilities: parse file names / infer keys
## ----------------------------------------------------------------------
## Supports filenames such as:
##   Summary_2PP_K_nk_tag_rep.rds
##   Summary_Testlet_K_nk_tag_rep.rds
## Also allows robustness if the directory indicates the model type.
parse_summary_name <- function(fname) {
  b <- basename(fname)
  m <- str_match(b, "^Summary(?:_(2PP|Testlet))?_(\\d+)_(\\d+)_([^_]+)_(\\d+)\\.rds$")
  if (any(is.na(m))) return(NULL)
  list(
    model = dplyr::case_when(
      !is.na(m[2]) & m[2] == "2PP"     ~ "2PP",
      !is.na(m[2]) & m[2] == "Testlet" ~ "Testlet2PP",
      TRUE ~ NA_character_
    ),
    K   = as.integer(m[3]),
    nk  = as.integer(m[4]),
    tag = m[5],
    rep = as.integer(m[6])
  )
}

## If model is not encoded in the filename, infer from the subdirectory
infer_model_from_path <- function(path) {
  if (grepl("/2PP/|\\\\2PP\\\\", path)) return("2PP")
  if (grepl("/Testlet2PP/|\\\\Testlet2PP\\\\", path)) return("Testlet2PP")
  NA_character_
}

## Build expected simulation path from the scenario keys
sim_path_from_keys <- function(K, nk, tag, GEN_TYPE, rep) {
  file.path(pathSim, sprintf("SimTestlet2PP_%d_%d_%s_%s_%d.rds", K, nk, tag, GEN_TYPE, rep))
}

## ----------------------------------------------------------------------
## Utilities: extract parameter families from a summary table
## ----------------------------------------------------------------------
## The summary tables are expected to have at least:
##   parameter, mean, sd
## Some summaries include:
##   n_eff, Rhat
## This helper ensures these columns exist (as NA if missing).
ensure_cols <- function(df) {
  if (!("n_eff" %in% names(df))) df$n_eff <- NA_real_
  if (!("Rhat"  %in% names(df))) df$Rhat  <- NA_real_
  df
}

## Extract a parameter family (e.g., "a", "b", "theta") and return a tidy table
## ordered by the parameter index.
extract_family <- function(df, family_prefix) {
  df <- ensure_cols(df)
  df %>%
    filter(str_detect(parameter, paste0("^", family_prefix, "\\["))) %>%
    arrange(as.integer(str_match(parameter, "\\[(\\d+)\\]")[, 2])) %>%
    transmute(
      index = as.integer(str_match(parameter, "\\[(\\d+)\\]")[, 2]),
      mean, sd, n_eff, Rhat
    )
}

## Simulation type used in filenames (must match the generator script)
GEN_TYPE <- "TOEP"

## ----------------------------------------------------------------------
## Scenario table (used to define the "true" rho per scenario)
## ----------------------------------------------------------------------
Scenario <- rbind(
  c(4, 2, .9, .7),
  c(5, 2, .7, .4),
  c(6, 2, 0,   0),
  c(4, 4, .7, .4),
  c(5, 4, 0,   0),
  c(6, 4, .9, .7),
  c(4, 5, 0,   0),
  c(5, 5, .9, .7),
  c(6, 5, .7, .4)
) %>% as.data.frame()
colnames(Scenario) <- c("K", "nk", "rho_ini", "rho_fim")

## Build the concatenated "true" rho vector for a given scenario:
##  - If tag == "rho0": all correlations are 0.
##  - Otherwise, interpret tag as the first rho (rho_ini) and locate the scenario row.
##  - For Toeplitz within testlets, the per-testlet rho vector has length nk-1,
##    and the global rho vector is repeated K times (one block per testlet).
true_rho_concat <- function(K, nk, tag, Scenario_df) {
  if (tag == "rho0") return(rep(0, K * (nk - 1)))
  
  rstart <- suppressWarnings(as.numeric(tag))
  cand <- Scenario_df %>% filter(K == !!K, nk == !!nk)
  if (nrow(cand) == 0) stop("Scenario not found for K=", K, " nk=", nk)
  
  ## If multiple candidates exist, try to match rho_ini to the numeric tag
  if (nrow(cand) > 1) {
    idx <- which(abs(cand$rho_ini - rstart) < 1e-6)
    cand <- if (length(idx) == 0) cand[1, , drop = FALSE] else cand[idx[1], , drop = FALSE]
  }
  
  v <- seq(cand$rho_ini, cand$rho_fim, length.out = nk - 1)
  rep(v, times = K)
}

## Extract rho parameters from the summary table.
## Tries multiple naming conventions to be robust across Stan files.
extract_rho <- function(df) {
  df <- ensure_cols(df)
  prefixes <- c("^rho_global\\[", "^rho_lag\\[", "^rho\\[")
  for (pref in prefixes) {
    tmp <- df %>%
      filter(str_detect(parameter, pref)) %>%
      arrange(as.integer(str_match(parameter, "\\[(\\d+)\\]")[, 2])) %>%
      transmute(
        index = as.integer(str_match(parameter, "\\[(\\d+)\\]")[, 2]),
        mean, sd, n_eff, Rhat
      )
    if (nrow(tmp) > 0) return(tmp)
  }
  NULL
}

## ----------------------------------------------------------------------
## List posterior summary files for each model directory
## ----------------------------------------------------------------------
list_model_files <- function(dir_path, model_name) {
  if (!dir.exists(dir_path)) return(character(0))
  
  ## Accepts summary files with flexible naming as long as they end with:
  ##   _K_nk_tag_rep.rds
  list.files(dir_path,
             pattern = "^Summary.*_\\d+_\\d+_.+_\\d+\\.rds$",
             full.names = TRUE) %>%
    as_tibble() %>%
    transmute(file = value, model = model_name)
}

summary_files_tbl <- bind_rows(
  list_model_files(pathFit_2PP,  "2PP"),
  list_model_files(pathFit_Test, "Testlet2PP")
)

stopifnot(nrow(summary_files_tbl) > 0)

## ----------------------------------------------------------------------
## Main loop: read each summary and compute error quantities
## ----------------------------------------------------------------------
results_long_all <- list()
rho_results_all  <- list()

for (k in seq_len(nrow(summary_files_tbl))) {
  
  sf_path <- summary_files_tbl$file[k]
  model0  <- summary_files_tbl$model[k]
  
  ## Parse keys from filename (fallback: infer from path and parse trailing fields)
  keys <- parse_summary_name(sf_path)
  if (is.null(keys)) {
    keys <- list(model = infer_model_from_path(sf_path), K = NA, nk = NA, tag = NA, rep = NA)
    m2 <- str_match(basename(sf_path), "_(\\d+)_(\\d+)_([^_]+)_(\\d+)\\.rds$")
    if (!any(is.na(m2))) {
      keys$K   <- as.integer(m2[2])
      keys$nk  <- as.integer(m2[3])
      keys$tag <- m2[4]
      keys$rep <- as.integer(m2[5])
    }
  }
  if (is.na(keys$model)) keys$model <- model0
  
  ## Skip files with unrecognized naming
  if (any(is.na(c(keys$K, keys$nk, keys$tag, keys$rep)))) {
    message("Unexpected filename format (skipping): ", basename(sf_path))
    next
  }
  
  K     <- keys$K
  nk1   <- keys$nk
  tag   <- keys$tag
  rp    <- keys$rep
  model <- keys$model
  
  ## Read summary table
  sum_df <- readRDS(sf_path)
  if (!all(c("parameter", "mean", "sd") %in% names(sum_df))) {
    message("Summary missing expected columns: ", basename(sf_path))
    next
  }
  
  ## Extract parameter families
  est_a     <- extract_family(sum_df, "a")
  est_b     <- extract_family(sum_df, "b")
  est_theta <- extract_family(sum_df, "theta")
  rho_est   <- extract_rho(sum_df)
  
  ## ------------------ rho truth and rho errors (if rho exists) ----------------
  if (!is.null(rho_est) && nrow(rho_est) > 0) {
    rho_true <- true_rho_concat(K, nk1, tag, Scenario)
    
    ## Safety: if lengths disagree, truncate to the minimum
    if (length(rho_true) != nrow(rho_est)) {
      warning("rho length mismatch: est(", nrow(rho_est), ") != true(", length(rho_true),
              ") in ", basename(sf_path), " — truncating.")
      m <- min(nrow(rho_est), length(rho_true))
      rho_est  <- rho_est[seq_len(m), , drop = FALSE]
      rho_true <- rho_true[seq_len(m)]
    }
    
    rho_tbl <- rho_est %>%
      mutate(
        param = "rho",
        true  = rho_true[index],
        bias  = mean - true,
        Abias = abs(mean - true),
        sqerr = (mean - true)^2,
        K = K, nk = nk1, tag = tag, rep = rp, model = model
      )
    
    rho_results_all[[sf_path]] <- rho_tbl
  }
  
  ## ------------------ a/b truth and errors -----------------------------------
  I_true <- nrow(zeta_verd)
  if (nrow(est_a) != I_true || nrow(est_b) != I_true) {
    warning("Item count mismatch: est_a=", nrow(est_a), " est_b=", nrow(est_b),
            " I_true=", I_true, " (", basename(sf_path), ")")
  }
  
  a_true <- zeta_verd$a_true[seq_len(nrow(est_a))]
  b_true <- zeta_verd$b_true[seq_len(nrow(est_b))]
  
  ## ------------------ theta truth (from simulation file) ----------------------
  ## We attempt to load sim_theta from the corresponding simulation .rds.
  sim_path <- sim_path_from_keys(K, nk1, tag, GEN_TYPE, rp)
  theta_true <- NULL
  
  if (file.exists(sim_path)) {
    dat_sim <- readRDS(sim_path)
    cand_names <- c("sim_theta", "theta", "thetas", "Theta", "sim.Theta")
    found <- cand_names[cand_names %in% names(dat_sim)]
    if (length(found) > 0) theta_true <- as.numeric(dat_sim[[found[1]]])
  } else {
    message("Simulation not found for theta: ", basename(sim_path))
  }
  
  ## Build long tables for a and b
  tb_a <- est_a %>%
    mutate(
      param = "a",
      true  = a_true[index],
      bias  = mean - true,
      Abias = abs(mean - true),
      sqerr = (mean - true)^2,
      K = K, nk = nk1, tag = tag, rep = rp, model = model
    )
  
  tb_b <- est_b %>%
    mutate(
      param = "b",
      true  = b_true[index],
      bias  = mean - true,
      Abias = abs(mean - true),
      sqerr = (mean - true)^2,
      K = K, nk = nk1, tag = tag, rep = rp, model = model
    )
  
  ## Theta table (only if lengths match)
  tb_theta <- NULL
  if (!is.null(theta_true) && nrow(est_theta) == length(theta_true)) {
    tb_theta <- est_theta %>%
      mutate(
        param = "theta",
        true  = theta_true[index],
        bias  = mean - true,
        Abias = abs(mean - true),
        sqerr = (mean - true)^2,
        K = K, nk = nk1, tag = tag, rep = rp, model = model
      )
  } else if (!is.null(theta_true)) {
    warning("Theta length mismatch: est=", nrow(est_theta), " true=", length(theta_true),
            " (", basename(sf_path), ")")
  }
  
  results_long_all[[sf_path]] <- bind_rows(tb_a, tb_b, tb_theta)
}

## Bind everything together (including rho if available)
results_long <- bind_rows(results_long_all) %>%
  arrange(model, param, K, nk, tag, rep, index)

rho_results <- bind_rows(rho_results_all)

results_long_ext <- bind_rows(results_long, rho_results) %>%
  arrange(model, param, K, nk, tag, rep, index)

## ----------------------------------------------------------------------
## Aggregate metrics (including rho)
## ----------------------------------------------------------------------
metrics_scenario_ext <- results_long_ext %>%
  group_by(model, param, K, nk, tag) %>%
  summarise(
    n = n(),
    bias_mean  = mean(bias, na.rm = TRUE),
    Abias_mean = mean(Abias, na.rm = TRUE),
    RMSE       = sqrt(mean(sqerr, na.rm = TRUE)),
    ESS_median = median(n_eff, na.rm = TRUE),
    Rhat_max   = max(Rhat, na.rm = TRUE),
    .groups = "drop"
  )

metrics_item_ext <- results_long_ext %>%
  group_by(model, param, K, nk, tag, index) %>%
  summarise(
    reps = n(),
    bias_mean  = mean(bias, na.rm = TRUE),
    Abias_mean = mean(Abias, na.rm = TRUE),
    RMSE       = sqrt(mean(sqerr, na.rm = TRUE)),
    ESS_median = median(n_eff, na.rm = TRUE),
    Rhat_max   = max(Rhat, na.rm = TRUE),
    .groups = "drop"
  )

diag_eff_ext <- results_long_ext %>%
  group_by(model, param, K, nk, tag) %>%
  summarise(
    ESS_min    = min(n_eff, na.rm = TRUE),
    ESS_q25    = quantile(n_eff, 0.25, na.rm = TRUE),
    ESS_median = median(n_eff, na.rm = TRUE),
    ESS_q75    = quantile(n_eff, 0.75, na.rm = TRUE),
    Rhat_max   = max(Rhat, na.rm = TRUE),
    .groups = "drop"
  )

## ----------------------------------------------------------------------
## Save RDS outputs
## ----------------------------------------------------------------------
saveRDS(results_long_ext,     file.path(pathResults, "results_long_params_with_rho_by_model.rds"))
saveRDS(metrics_scenario_ext, file.path(pathResults, "metrics_by_scenario_with_rho_by_model.rds"))
saveRDS(metrics_item_ext,     file.path(pathResults, "metrics_by_item_or_theta_with_rho_by_model.rds"))
saveRDS(diag_eff_ext,         file.path(pathResults, "diagnostics_efficiency_with_rho_by_model.rds"))

## ----------------------------------------------------------------------
## Plot styling helpers
## ----------------------------------------------------------------------
## Scenario labeling based on the tag value
label_rho <- function(tag) {
  case_when(
    tag == "rho0" ~ "Rho = Null",
    tag == "0.90" ~ "Rho = Strong",
    tag == "0.70" ~ "Rho = Moderate",
    TRUE ~ paste0("Rho = ", tag)
  )
}

## Adds a scenario label and orders scenarios by (K, nk, tag)
base_prepare <- function(df) {
  df %>%
    mutate(
      RhoLabel = label_rho(tag),
      scenario = paste0("K=", K, ", nk=", nk, ", ", RhoLabel)
    ) %>%
    group_by(param) %>%
    mutate(scenario = factor(scenario, levels = unique(scenario[order(K, nk, tag)]))) %>%
    ungroup()
}

## Grayscale palette for model comparison
model_levels <- c("2PP", "Testlet2PP")

add_model_greys <- function() {
  ggplot2::scale_fill_manual(
    name   = "Model",
    breaks = model_levels,
    values = c("2PP" = "grey30", "Testlet2PP" = "grey70")
  )
}

## ----------------------------------------------------------------------
## Boxplots: Bias by scenario (comparison)
## ----------------------------------------------------------------------
bias_long <- base_prepare(results_long_ext)

plot_box_bias_param_cmp <- function(df_param, param_name) {
  ggplot(df_param, aes(x = scenario, y = bias, fill = model)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_boxplot(outlier.alpha = 0.4, position = position_dodge(width = 0.8)) +
    labs(title = "", x = "Scenario", y = "Bias", fill = "Model") +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, size = 10),
      plot.title   = element_text(face = "bold")
    ) +
    add_model_greys()
}

## 1) Separate plot per parameter
for (p in sort(unique(bias_long$param))) {
  dfp <- filter(bias_long, param == p)
  if (nrow(dfp) == 0) next
  gp <- plot_box_bias_param_cmp(dfp, p)
  
  ggsave(file.path(pathResults, paste0("cmp_boxplot_bias_por_cenario_", p, ".pdf")),
         gp, width = 11, height = 6.2)
  
  ggsave(file.path(pathResults, paste0("cmp_boxplot_bias_por_cenario_", p, ".png")),
         gp, width = 1400, height = 800, units = "px", dpi = 120)
}


param_labels <- c(
  a     = "Discrimination",
  b     = "Difficulty",
  rho   = "Correlation",
  theta = "Ability"
)

## 2) Combined panel (facet by parameter)

#setEPS()
#postscript(file = file.path(pathResults, "cmp_boxplot_bias_por_cenario_todos.eps"),width =12.5, height = 9)

gp_all <- ggplot(bias_long, aes(x = scenario, y = bias, fill = model)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(outlier.alpha = 1, position = position_dodge(width = 0.8)) +
  labs(title = "", x = "Scenario", y = "Bias") +
  facet_wrap(~ param, scales = "free_y", ncol = 2,
             labeller = as_labeller(param_labels)) +
  theme_minimal(base_size = 16) +  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title  = element_text(face = "bold", size = 18),  
    strip.text  = element_text(size = 14, face = "bold"),  
    axis.title  = element_text(size = 14),                 
    legend.text = element_text(size = 12),                 
    legend.title = element_text(size = 13)                 
  ) +
  add_model_greys()




ggsave(file.path(pathResults, "cmp_boxplot_bias_por_cenario_todos.eps"),
       gp_all, width = 12.5, height = 9)

ggsave(file.path(pathResults, "cmp_boxplot_bias_por_cenario_todos.png"),
       gp_all, width = 1600, height = 1150, units = "px", dpi = 120)

## ----------------------------------------------------------------------
## Boxplots: Absolute Bias (|Bias|) by scenario
## ----------------------------------------------------------------------
plot_box_Abias_param_cmp <- function(df_param, param_name) {
  ggplot(df_param, aes(x = scenario, y = Abias, fill = model)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_boxplot(outlier.alpha = 0.4, position = position_dodge(width = 0.8)) +
    labs(title = "", x = "Scenario", y = "ABias", fill = "Model") +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      plot.title  = element_text(face = "bold")
    ) +
    add_model_greys()
}

for (p in sort(unique(bias_long$param))) {
  dfp <- filter(bias_long, param == p)
  if (nrow(dfp) == 0) next
  gp <- plot_box_Abias_param_cmp(dfp, p)
  
  ggsave(file.path(pathResults, paste0("cmp_boxplot_Abias_por_cenario_", p, ".pdf")),
         gp, width = 11, height = 6.2)
  
  ggsave(file.path(pathResults, paste0("cmp_boxplot_Abias_por_cenario_", p, ".png")),
         gp, width = 1400, height = 800, units = "px", dpi = 120)
}

## ----------------------------------------------------------------------
## Boxplots: RMSE by scenario (computed per replication)
## ----------------------------------------------------------------------
rmse_long <- base_prepare(results_long_ext)

rmse_rep <- rmse_long %>%
  group_by(model, param, K, nk, tag, rep, scenario) %>%
  summarise(RMSE = sqrt(mean(sqerr, na.rm = TRUE)), .groups = "drop")

plot_box_rmse_param_cmp <- function(df_param, param_name) {
  ggplot(df_param, aes(x = scenario, y = RMSE, fill = model)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_boxplot(outlier.alpha = 0.4, position = position_dodge(width = 0.8)) +
    labs(title = "", x = "Scenario", y = "RMSE", fill = "Model") +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      plot.title  = element_text(face = "bold")
    ) +
    add_model_greys()
}

for (p in sort(unique(rmse_rep$param))) {
  dfp <- filter(rmse_rep, param == p)
  if (nrow(dfp) == 0) next
  gp <- plot_box_rmse_param_cmp(dfp, p)
  
  ggsave(file.path(pathResults, paste0("cmp_boxplot_RMSE_por_cenario_", p, ".pdf")),
         gp, width = 11, height = 6.2)
  
  ggsave(file.path(pathResults, paste0("cmp_boxplot_RMSE_por_cenario_", p, ".png")),
         gp, width = 1400, height = 800, units = "px", dpi = 120)
}


gp_rmse_all <- ggplot(rmse_rep, aes(x = scenario, y = RMSE, fill = model)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(outlier.alpha = 1, position = position_dodge(width = 0.8)) +
  labs(title = "", x = "Scenario", y = "RMSE", fill = "Model") +
  facet_wrap(~ param, scales = "free_y", ncol = 2,
             labeller = as_labeller(param_labels)) +
  theme_minimal(base_size = 16) +  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title  = element_text(face = "bold", size = 18),  
    strip.text  = element_text(size = 14, face = "bold"),  
    axis.title  = element_text(size = 14),                 
    legend.text = element_text(size = 12),                 
    legend.title = element_text(size = 13)                 
  ) +
  add_model_greys()




ggsave(file.path(pathResults, "cmp_boxplot_RMSE_por_cenario_todos.eps"),
       gp_rmse_all, width = 12.5, height = 9)

ggsave(file.path(pathResults, "cmp_boxplot_RMSE_por_cenario_todos.png"),
       gp_rmse_all, width = 1600, height = 1150, units = "px", dpi = 120)

## ----------------------------------------------------------------------
## Boxplots: ESS by scenario (median ESS per replication)
## ----------------------------------------------------------------------
ess_long <- base_prepare(results_long_ext)

ess_rep <- ess_long %>%
  group_by(model, param, K, nk, tag, rep, scenario) %>%
  summarise(ESS = stats::median(n_eff, na.rm = TRUE), .groups = "drop")

plot_box_ess_param_cmp <- function(df_param, param_name) {
  ggplot(df_param, aes(x = scenario, y = ESS, fill = model)) +
    geom_hline(yintercept = 100, linetype = "dashed") +
    geom_boxplot(outlier.alpha = 0.4, position = position_dodge(width = 0.8)) +
    labs(title = "", x = "Scenario", y = "ESS", fill = "Model") +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      plot.title  = element_text(face = "bold")
    ) +
    add_model_greys()
}

for (p in sort(unique(ess_rep$param))) {
  dfp <- filter(ess_rep, param == p)
  if (nrow(dfp) == 0) next
  gp <- plot_box_ess_param_cmp(dfp, p)
  
  ggsave(file.path(pathResults, paste0("cmp_boxplot_ESS_por_cenario_", p, ".pdf")),
         gp, width = 11, height = 6.2)
  
  ggsave(file.path(pathResults, paste0("cmp_boxplot_ESS_por_cenario_", p, ".png")),
         gp, width = 1400, height = 800, units = "px", dpi = 120)
}

gp_ess_all <- ggplot(ess_rep, aes(x = scenario, y = ESS, fill = model)) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  geom_boxplot(outlier.alpha = 0.4, position = position_dodge(width = 0.8)) +
  labs(title = "", x = "Scenario", y = "ESS", fill = "Model") +
  facet_wrap(~ param, scales = "free_y", ncol = 2,
             labeller = as_labeller(param_labels)) +
  theme_minimal(base_size = 16) +  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title  = element_text(face = "bold", size = 18),  
    strip.text  = element_text(size = 14, face = "bold"),  
    axis.title  = element_text(size = 14),                 
    legend.text = element_text(size = 12),                 
    legend.title = element_text(size = 13)                 
  ) +
  add_model_greys()

ggsave(file.path(pathResults, "cmp_boxplot_ESS_por_cenario_todos.pdf"),
       gp_ess_all, width = 12.5, height = 9)

ggsave(file.path(pathResults, "cmp_boxplot_ESS_por_cenario_todos.png"),
       gp_ess_all, width = 1600, height = 1150, units = "px", dpi = 120)

## ----------------------------------------------------------------------
## Boxplots: Rhat by scenario (max Rhat per replication)
## ----------------------------------------------------------------------
rhat_long <- base_prepare(results_long_ext)

rhat_rep <- rhat_long %>%
  group_by(model, param, K, nk, tag, rep, scenario) %>%
  summarise(Rhat = max(Rhat, na.rm = TRUE), .groups = "drop")

plot_box_rhat_param_cmp <- function(df_param, param_name) {
  ggplot(df_param, aes(x = scenario, y = Rhat, fill = model)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_boxplot(outlier.alpha = 0.4, position = position_dodge(width = 0.8)) +
    labs(title = "", x = "Scenario", y = "Rhat", fill = "Model") +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      plot.title  = element_text(face = "bold")
    ) +
    add_model_greys()
}

for (p in sort(unique(rhat_rep$param))) {
  dfp <- filter(rhat_rep, param == p)
  if (nrow(dfp) == 0) next
  gp <- plot_box_rhat_param_cmp(dfp, p)
  
  ggsave(file.path(pathResults, paste0("cmp_boxplot_Rhat_por_cenario_", p, ".pdf")),
         gp, width = 11, height = 6.2)
  
  ggsave(file.path(pathResults, paste0("cmp_boxplot_Rhat_por_cenario_", p, ".png")),
         gp, width = 1400, height = 800, units = "px", dpi = 120)
}

gp_rhat_all <- ggplot(rhat_rep, aes(x = scenario, y = Rhat, fill = model)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_boxplot(outlier.alpha = 0.4, position = position_dodge(width = 0.8)) +
  labs(title = "", x = "Scenario", y = "Rhat", fill = "Model") +
  facet_wrap(~ param, scales = "free_y", ncol = 2,
             labeller = as_labeller(param_labels)) +
  theme_minimal(base_size = 16) +  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title  = element_text(face = "bold", size = 18),  
    strip.text  = element_text(size = 14, face = "bold"),  
    axis.title  = element_text(size = 14),                 
    legend.text = element_text(size = 12),                 
    legend.title = element_text(size = 13)                 
  ) +
  add_model_greys()

ggsave(file.path(pathResults, "cmp_boxplot_Rhat_por_cenario_todos.pdf"),
       gp_rhat_all, width = 12.5, height = 9)

ggsave(file.path(pathResults, "cmp_boxplot_Rhat_por_cenario_todos.png"),
       gp_rhat_all, width = 1600, height = 1150, units = "px", dpi = 120)

message("\n===== 2PP vs Testlet2PP comparison finished =====")
