## =============================================================================
## Testlet 2PP vs 2PP em dados reais (Seduc/CE) com HPC e comparação LOO/WAIC
## José R. S. Santos — organização e comentários: 2025-10-31
## =============================================================================


## -------------------- 0) Caminhos e fontes auxiliares ------------------------
## Ajuste estes caminhos se necessário. Mantive os absolutos que você já usa.
root_local <- "~/GitHub/bayesian-testlet-antedependence"
path_project <- "~/GitHub/bayesian-testlet-antedependence/Real Data Analysis/Large-Scale Educational Assessment"

stopifnot(dir.exists(root_local), dir.exists(path_project))

setwd(path_project)  # mantém compatibilidade com seus scripts


## ------------------------- 1) Setup e pacotes --------------------------------
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

## Função utilitária simples para medir tempo de blocos
tic <- function(msg) { cat(sprintf("\n[START] %s ... %s\n", msg, Sys.time())); Sys.time() }
toc <- function(t0)  { cat(sprintf("[ END ] Elapsed: %s\n", Sys.time() - t0)) }


## Path to Stan programs
pathProgram <- here(root_local, "Programs")

## Source helper functions used for indexing + diagnostics (Q3, envelopes, heatmaps, etc.)
source(here('Real Data Analysis', "HelpersRealDataAnal.R"))


## Output directories
pathFit <- here(path_project, "Fit")
dir.create(pathFit, showWarnings = FALSE, recursive = TRUE)
stopifnot(dir.exists(pathFit))

saveFigures <- here(pathFit, "Figures")
dir.create(saveFigures, showWarnings = FALSE, recursive = TRUE)
stopifnot(dir.exists(saveFigures))



## ------------------- 2) Carregar dados de treino (3ª série) ------------------
## (a) Base de treino
file_train <- file.path(path_project,'Binary_3EM_train_n2k.RDS')
stopifnot(file.exists(file_train))
mYcNA <- readRDS(file_train)


## Checagens básicas
n  <- nrow(mYcNA)
vI <- ncol(mYcNA)
stopifnot(n > 0, vI > 0)

## -------------------- 3) Estrutura de testlets e índices ---------------------
## Definição declarada por você:
K  <- 6
nk <- c(3, 2, 2, 2, 3, 2)              # comprimentos dos testlets (ordem natural)
dk <- c(4, 9, 15, 17, 19, 25)          # posições iniciais dos testlets (se necessário)
is_HU <- c(FALSE, TRUE, TRUE, TRUE, FALSE, TRUE)  # k=1 e 5 em HT; demais HU (seu padrão)

## Índices auxiliares p/ estrutura HU/HT usada nas funções auxiliares
idx <- make_rho_index(nk, is_HU)       # fornece $rho_len, $rho_start etc.

## Itens independentes (fora de testlets)
ind_items <- c(1:3, 7, 8, 11:14, 22:24)

## Consistência: itens de testlet são o complemento
all_items      <- seq_len(vI)
testlet_items  <- setdiff(all_items, ind_items)
stopifnot(length(testlet_items) == sum(nk))

## Particionar itens por testlet, na ordem declarada em nk
idx_testlets <- split(sort(testlet_items), rep(seq_along(nk), times = nk))
stopifnot(sum(lengths(idx_testlets)) == length(testlet_items),
          length(idx_testlets) == K)

## --------------------- 4) Listas de dados para o Stan ------------------------
## Testlet-2PP (probit) — versão “HT_Diag” (como no seu arquivo Stan)
data_testlet <- list(
  I = vI, N = n, K = K, dk = dk, nk = nk,
  ind_items = ind_items, n_ind = vI - sum(nk),
  Y = mYcNA,
  sigma_a = .6, sigma_b = 4, sigma_rho = 1,
  rho_len = idx$rho_len, S_mc = 200,
  rho_start = idx$rho_start
)

## 2PP (probit) sem dependência local
data_2pp <- list(
  N = n, vI = vI, Y = mYcNA,
  sigma_a = .6, sigma_b = 4
)

## ---------------------- 5) Valores iniciais dos parâmetros -------------------
## Escore bruto padronizado para inicializar theta
scores <- scale(rowSums(mYcNA))[ ,1]

## (i) Testlet-2PP: inclui rho_global
init_testlet <- function() {
  list(
    theta = as.numeric(scores),
    a     = rep(0.1, vI),
    b     = rep(0.1, vI),
    rho_global = rep(0.1, sum(idx$rho_len))
  )
}

## (ii) 2PP: não há rho_global
init_2pp <- function() {
  list(
    theta = as.numeric(scores),
    a     = rep(0.1, vI),
    b     = rep(0.1, vI)
  )
}

## Parâmetros monitorados
pars_testlet <- c("a", "b", "theta", "rho_global",'log_lik')
pars_2pp     <- c("a", "b", "theta", "log_lik")

## NUTS: configurações comuns
nChains       <- 1
burnInSteps   <- 1000
thinSteps     <- 1
numSavedSteps <- 1000
nIter         <- ceiling(burnInSteps + numSavedSteps * thinSteps)
ctrl_nuts     <- list(adapt_delta = 0.8, max_treedepth = 10)

## --------------------------- 6) Ajustar os modelos ---------------------------
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

## (b) 2PP
t0 <- tic("Fitting 2PP")
fit_2pp <- stan(
  data   = data_2pp,
  file   = file.path(pathProgram, "2PPModelVec_Diag.stan"),
  init   = init_2pp,
  chains = nChains, pars = pars_2pp,
  iter   = nIter, warmup = burnInSteps, thin = thinSteps,
  control = list(adapt_delta = 0.8)  # depth default ok
)
toc(t0)

## Salvar objetos de ajuste
saveRDS(fit_testlet, file.path(pathFit, "ResultTestlet2PP_3EM_municipio2k.rds"))
saveRDS(fit_2pp,     file.path(pathFit, "Result2PP_3EM_municipio2k.rds"))



fit_testlet<-readRDS(file.path(pathFit,'ResultTestlet2PP_3EM_municipio2k.rds'))
fit_2pp<-readRDS(file.path(pathFit,"Result2PP_3EM_municipio2k.rds"))



## ------------------------- 7) Holdout (estudantes) ---------------------------
## (a) Carregar base de *teste* e selecionar itens nas MESMAS colunas
file_test <- file.path(path_project,'Binary_3EM_test_n2k.RDS')
stopifnot(file.exists(file_test))
mYcNA_test <- readRDS(file_test)

N_test <- nrow(mYcNA_test)

#saveRDS(mYcNA_test,'~/GitHub/bayesian-testlet-antedependence/Real Data Analysis/Large-Scale Educational Assessment/Binary_3EM_test_n2k.RDS')

## ------------------------- 8) Pós-processamento básico -----------------------
## Extrair listas “arrays” para funções *build_params_*_from_arrays
post_testlet <- rstan::extract(fit_testlet, pars = c("a", "b",'theta', 'rho_global'), permuted = TRUE)
post_2pp     <- rstan::extract(fit_2pp,     pars = c("a", "b", 'theta'),             permuted = TRUE)


# ---- por TESTLET (facetas T1, T2, ...) ----

env_test <- q3_envelope_data(post_testlet, function(draw_id, ext)
  build_params_testlet_from_arrays(draw_id, ext, idx, nk, K), mYcNA_test, idx_testlets,
  model_label = "Testlet 2PP", ind_items = integer(0),
  R_eval = 200, seed = 123, global = FALSE)

env_2pp  <- q3_envelope_data(post_2pp,  function(draw_id, ext)
  build_params_2pp_from_arrays(draw_id, ext), mYcNA_test, idx_testlets,
  model_label = "2PP",     ind_items = integer(0),
  R_eval = 200, seed = 123, global = FALSE)

df_plot  <- dplyr::bind_rows(env_test, env_2pp)

p_env<-plot_q3_ecdf_envelope_gray(df_plot) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
  theme(
    axis.text.x = element_text(size = 8),
    panel.spacing.x = unit(1.5, "lines")
  )

ggsave(file.path(saveFigures, "ECDF_3EM.pdf"), p_env, width = 12, height = 5,units = 'in')

# ---- Global ----

env_test_g <-q3_envelope_data(post_testlet, function(draw_id, ext)
  build_params_testlet_from_arrays(draw_id, ext, idx, nk, K),mYcNA_test, idx_testlets,
  model_label = "Testlet 2PP", ind_items = integer(0),
  R_eval = 200, global = T,seed = 234)

env_2pp_g  <- q3_envelope_data(post_2pp,  function(draw_id, ext)
  build_params_2pp_from_arrays(draw_id, ext),  mYcNA_test, idx_testlets,
  model_label = "2PP",     ind_items = integer(0),
  R_eval = 200, global = T)

p_env_g<-plot_q3_ecdf_envelope_gray(bind_rows(env_test_g, env_2pp_g))

ggsave(file.path(saveFigures, "ECDF_Global_3EM_25k.pdf"), p_env_g, width = 18, height = 12,units = 'cm')



# Heatmaps differences

df_testlet <- delta_q3_testlet(post_testlet, function(draw_id, ext)
  build_params_testlet_from_arrays(draw_id, ext, idx, nk, K), mYcNA_test, idx_testlets,
  R_eval = 100, model_label = "Testlet")

df_2pp <- delta_q3_testlet(post_2pp, function(draw_id, ext)
  build_params_2pp_from_arrays(draw_id, ext), mYcNA_test, idx_testlets,
  R_eval = 100, model_label = "2PP")

save_delta_q3_heatmap_pdf(df_testlet, df_2pp,
                          diff_only = FALSE,
                          file =file.path(saveFigures, "Heatmap_DeltaQ3_Testlet_vs_2PP_3EM_25k.pdf"))


save_delta_q3_heatmap_pdf(df_testlet, df_2pp,
                          diff_only = TRUE,
                          file =file.path(saveFigures, "Heatmap_DeltaQ3_Diferencas_3EM.pdf"))



# Réplicas para cada modelo (pode aumentar R_eval)
rep_T  <- q3bar_replicates(post_testlet, function(draw_id, ext)
  build_params_testlet_from_arrays(draw_id, ext, idx, nk, K), mYcNA_test, idx_testlets,
  R_eval = 300, seed = 123)

rep_2P <- q3bar_replicates(post_2pp,  function(draw_id, ext)
  build_params_2pp_from_arrays(draw_id, ext),  mYcNA_test, idx_testlets,
  R_eval = 300, seed = 123)

# Probabilidade de superioridade (global e por testlet)
tab_sup <- model_superiority_probability(rep_T, rep_2P, idx_testlets)
tab_sup
#> tibble com colunas: Scope, n, p_sup, se, lo95, hi95

# Gráfico em tons de cinza + salvar PDF
save_superiority_pdf(tab_sup,file.path(saveFigures, "Prob_Superioridade_Testlet_vs_2PP.pdf"))






## --------------------- 10) Comparação de modelos (LOO/WAIC) -----------------
## LOO com moment matching para ambos (melhor estabilidade)
loo_with_mm <- function(fit, cores = 4) {
  ll <- loo::extract_log_lik(fit, merge_chains = FALSE)  # S x C x N
  S  <- dim(ll)[1]; C <- dim(ll)[2]; N <- dim(ll)[3]
  dim(ll) <- c(S*C, N)                                  # draws x N
  r_eff <- loo::relative_eff(exp(ll), chain_id = rep(1:C, each = S))
  loo::loo(ll, r_eff = r_eff, moment_match = TRUE, cores = cores)
}

loo_testlet_mm <- loo_with_mm(fit_testlet)
loo_2pp_mm     <- loo_with_mm(fit_2pp)

## Tabelas úteis
loo::pareto_k_table(loo_testlet_mm)
loo::pareto_k_table(loo_2pp_mm)

## Comparação (quanto menor o elpd_diff melhor o segundo argumento)
comp_loo  <- loo::loo_compare(loo_testlet_mm, loo_2pp_mm)
print(comp_loo, simplify = FALSE, digits = 3)

## WAIC (opcional, por completude)
waic_testlet <- loo::waic(loo::extract_log_lik(fit_testlet))
waic_2pp     <- loo::waic(loo::extract_log_lik(fit_2pp))
comp_waic    <- loo::loo_compare(waic_testlet, waic_2pp)
print(comp_waic)

cat("\n[OK] Pipeline concluído.\n")
