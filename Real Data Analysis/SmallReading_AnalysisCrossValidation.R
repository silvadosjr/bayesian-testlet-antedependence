

## -------------------- 0) Caminhos e fontes auxiliares ------------------------
## Ajuste estes caminhos se necessário. Mantive os absolutos que você já usa.
root_local <- "C:/Users/Usuário/OneDrive/Documentos"
path_project <- "C:/Users/Usuário/OneDrive/Documentos/Artigos/IRT Residual dependency"

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


pathProgram <- here("Programas")
stopifnot(dir.exists(pathProgram))

## Fontes auxiliares (funções usadas abaixo)
#source(here('Programas', "AuxFunctions.R"))
source(here('Programas', "HelpersRealDataAnal.R"))

## Diretórios para saída
pathFit <- 'D:/IRT Residual dependence/Real data analysis/CrossValidation/Fit'
dir.create(pathFit, showWarnings = FALSE, recursive = TRUE)
stopifnot(dir.exists(pathFit))

saveFigures<-'D:/IRT Residual dependence/Real data analysis/CrossValidation/Fit09112025/Figuras'
dir.create(saveFigures, showWarnings = FALSE, recursive = TRUE)
stopifnot(dir.exists(saveFigures))


## ------------------- 2) Carregar dados de treino (3ª série) ------------------
## (a) Base de treino
file_train <- 'C:/Users/Usuário/OneDrive/Documentos/Artigos/IRT Residual dependency/Dados Reais/ld_train.rds'
stopifnot(file.exists(file_train))
base_train <- readRDS(file_train)



mYc   <- base_train[,-1]             # 0/1 com NAs para ausentes
mYcNA <- as.data.frame(mYc)
mYcNA[is.na(mYc)] <- 0L                       # geralmente todos apresentaram, então no-op

## Checagens básicas
n  <- nrow(mYcNA)
vI <- ncol(mYcNA)
stopifnot(n > 0, vI > 0)

## Organizando por testlet

mYcNA<-mYcNA[,c(1,6,4,5,2,3)]


## -------------------- 3) Estrutura de testlets e índices ---------------------
## Definição declarada por você:
K  <- 3
nk <- c(2,2,2)             # comprimentos dos testlets (ordem natural)
dk <- c(1,3,5)          # posições iniciais dos testlets (se necessário)
is_HU <- c(TRUE, TRUE, TRUE)  # k=1 e 5 em HT; demais HU (seu padrão)

## Índices auxiliares p/ estrutura HU/HT usada nas funções auxiliares
idx <- make_rho_index(nk, is_HU)       # fornece $rho_len, $rho_start etc.

## Itens independentes (fora de testlets)
ind_items <- integer(0L)

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
pars_testlet <- c("a", "b", "theta", "rho_global", "log_lik")
pars_2pp     <- c("a", "b", "theta", "log_lik")

## NUTS: configurações comuns
nChains       <- 1
burnInSteps   <- 1000
thinSteps     <- 20
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
saveRDS(fit_testlet, file.path(pathFit, "ResultTestlet2PP_SmallReading.rds"))
saveRDS(fit_2pp,     file.path(pathFit, "Result2PP_SmallReading.rds"))


fit_testlet<-readRDS(file.path(pathFit,"ResultTestlet2PP_SmallReading.rds"))
fit_2pp<-readRDS(file.path(pathFit,"Result2PP_SmallReading.rds"))



## ------------------------- 7) Holdout (estudantes) ---------------------------
## (a) Carregar base de *teste* e selecionar itens nas MESMAS colunas
file_test <- 'C:/Users/Usuário/OneDrive/Documentos/Artigos/IRT Residual dependency/Dados Reais/ld_test.rds'
stopifnot(file.exists(file_test))
base_test <- readRDS(file_test)


## (b) Corrigir itens (mesmo gabarito)
mYc_test   <- base_test[,-1]
mYcNA_test <- as.data.frame(mYc_test)
mYcNA_test[is.na(mYc_test)] <- 0L
N_test <- nrow(mYcNA_test)

## Organizando por testlet

mYcNA_test<-mYcNA_test[,c(1,6,4,5,2,3)]


## ------------------------- 8) Pós-processamento básico -----------------------
## Extrair listas “arrays” para funções *build_params_*_from_arrays
post_testlet <- rstan::extract(fit_testlet, pars = c("a", "b", "rho_global"), permuted = TRUE)
post_2pp     <- rstan::extract(fit_2pp,     pars = c("a", "b"),             permuted = TRUE)

## ----------------------- 9) HPC: logscore e residQ3 --------------------------
## Número de réplicas usadas na HPC
R_eval      <- 250
R_eval_aug  <- 150  # se quiser aumentar, custo ~ linear


## (a) Testlet-2PP

# hpc_logaug_testlet <- compute_hpc_aug_arrays(
#   ext = post_testlet,
#   build_params_fun = function(draw_id, ext)
#     build_params_testlet_from_arrays(draw_id, ext, idx, nk, K),
#   Y_new        = mYcNA_test,
#   idx_testlets = idx_testlets,
#   ind_items    = ind_items,
#   R_eval       = R_eval,   # pode usar R_eval_aug se quiser expandir
#   seed         = 2025
# )


hpc_q3_testlet <- compute_hpc_arrays(
  ext = post_testlet,
  build_params_fun = function(draw_id, ext)
    build_params_testlet_from_arrays(draw_id, ext, idx, nk, K),
  Y_new        = mYcNA_test,
  idx_testlets = idx_testlets,
  ind_items    = ind_items,
  diag         = "residQ3",
  R_eval       = R_eval,
  by_testlet = T,
  seed         = 70
)


## (b) 2PP (atenção: builder do 2PP em TODAS as chamadas)

# hpc_logaug_2pp <- compute_hpc_aug_arrays(
#   ext = post_2pp,
#   build_params_fun = function(draw_id, ext)
#     build_params_2pp_from_arrays(draw_id, ext),
#   Y_new        = mYcNA_test,
#   idx_testlets = idx_testlets,
#   ind_items    = ind_items,
#   R_eval       = R_eval,   # idem observação acima
#   seed         = 2025
# )


hpc_q3_2pp <- compute_hpc_arrays(
  ext = post_2pp,
  build_params_fun = function(draw_id, ext)
    build_params_2pp_from_arrays(draw_id, ext), 
  Y_new        = mYcNA_test,
  idx_testlets = idx_testlets,
  ind_items    = integer(0),   # todos independentes
  diag         = "residQ3",
  R_eval       = R_eval,
  by_testlet = T,
  seed         = 1002
)

## (c) Resumo HPC
hpc_summary <- list(
  R_eval = R_eval,
  N_test = N_test,
  logscore = list(
    Testlet2PP = hpc_logaug_testlet$p_HPC,   # usamos a versão augmented
    M2PP       = hpc_logaug_2pp$p_HPC
  ),
  residQ3 = list(
    Testlet2PP = hpc_q3_testlet$p_HPC,
    M2PP       = hpc_q3_2pp$p_HPC
  )
)

cat("\n===== HPC (augmented) p-values =====\n")
print(list(
  logscore_aug = list(
    Testlet2PP = hpc_logaug_testlet$p_HPC,
    M2PP       = hpc_logaug_2pp$p_HPC
  )
))

## Salvar artefatos da HPC
saveRDS(hpc_logaug_testlet, file.path(pathFit, "HPC_logscore_aug_Testlet2PP_n2k.RDS"))
saveRDS(hpc_logaug_2pp,     file.path(pathFit, "HPC_logscore_aug_2PP_n2k.RDS"))
saveRDS(hpc_q3_testlet,     file.path(pathFit, "HPC_residQ3_Testlet2PP_n2k.RDS"))
saveRDS(hpc_q3_2pp,         file.path(pathFit, "HPC_residQ3_2PP_n2k.RDS"))





# ---- por TESTLET (facetas T1, T2, ...) ----

env_test <- q3_envelope_data(post_testlet, function(draw_id, ext)
  build_params_testlet_from_arrays(draw_id, ext, idx, nk, K), mYc_test, idx_testlets,
  model_label = "Testlet 2PP", ind_items = integer(0),
  R_eval = 200, seed = 123, global = FALSE)

env_2pp  <- q3_envelope_data(post_2pp,  function(draw_id, ext)
  build_params_2pp_from_arrays(draw_id, ext), mYc_test, idx_testlets,
  model_label = "2PP",     ind_items = integer(0),
  R_eval = 200, seed = 123, global = FALSE)

df_plot  <- dplyr::bind_rows(env_test, env_2pp)

p_env<-plot_q3_ecdf_envelope_gray(df_plot)

ggsave(file.path(saveFigures, "ECDF_SmallReading.pdf"), p_env, width = 10, height = 5,units = 'in')

# ---- Global ----

env_test_g <-q3_envelope_data(post_testlet, function(draw_id, ext)
  build_params_testlet_from_arrays(draw_id, ext, idx, nk, K), mYc_test, idx_testlets,
  model_label = "Testlet 2PP", ind_items = integer(0),
  R_eval = 200, global = T,seed = 123)

env_2pp_g  <- q3_envelope_data(post_2pp,  function(draw_id, ext)
  build_params_2pp_from_arrays(draw_id, ext),  mYc_test, idx_testlets,
  model_label = "2PP",     ind_items = integer(0),
  R_eval = 200, global = T,seed = 123)

p_env_g<-plot_q3_ecdf_envelope_gray(bind_rows(env_test_g, env_2pp_g))

ggsave(file.path(saveFigures, "ECDF_Global_SmallReading.pdf"), p_env_g, width = 18, height = 12,units = 'cm')

# Heatmaps differences

df_testlet <- delta_q3_testlet(post_testlet, function(draw_id, ext)
  build_params_testlet_from_arrays(draw_id, ext, idx, nk, K),mYc_test, idx_testlets,
  R_eval = 100, model_label = "Testlet 2PP")

df_2pp <- delta_q3_testlet(post_2pp, function(draw_id, ext)
  build_params_2pp_from_arrays(draw_id, ext),mYc_test, idx_testlets,
  R_eval = 100, model_label = "2PP")

save_delta_q3_heatmap_pdf(df_testlet, df_2pp,
                          diff_only = FALSE,
                          file =file.path(saveFigures, "Heatmap_DeltaQ3_Testlet_vs_2PP_SmallReading.pdf"))


save_delta_q3_heatmap_pdf(df_testlet, df_2pp,
                          diff_only = TRUE,
                          file =file.path(saveFigures, "Heatmap_DeltaQ3_Diferencas_SmallReading.pdf"))



# Réplicas para cada modelo (pode aumentar R_eval)
rep_T  <- q3bar_replicates(post_testlet, function(draw_id, ext)
  build_params_testlet_from_arrays(draw_id, ext, idx, nk, K),mYc_test, idx_testlets,
  R_eval = 300, seed = 123)

rep_2P <- q3bar_replicates(post_2pp,  function(draw_id, ext)
  build_params_2pp_from_arrays(draw_id, ext),  mYc_test, idx_testlets,
  R_eval = 300, seed = 123)

# Probabilidade de superioridade (global e por testlet)
tab_sup <- model_superiority_probability(rep_T, rep_2P, idx_testlets)
tab_sup
#> tibble com colunas: Scope, n, p_sup, se, lo95, hi95

# Gráfico em tons de cinza + salvar PDF
save_superiority_pdf(tab_sup, "Prob_Superioridade_Testlet_vs_2PP.pdf")





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


## ======================= Posterior estimates ==============================##

post_summary <- function(fit, pars) {
  ss <- rstan::summary(fit, pars = pars)$summary
  out <- as.data.frame(ss)
  out$Parameter <- rownames(ss)
  rownames(out) <- NULL
  out
}

## --- Resumos (a, b, theta, rho) ---
sum_a1     <- post_summary(fit_testlet, "a")
sum_b1     <- post_summary(fit_testlet, "b")
sum_theta1 <- post_summary(fit_testlet, "theta")
sum_rho    <- post_summary(fit_testlet, "rho_global")

sum_a2     <- post_summary(fit_2pp, "a")
sum_b2     <- post_summary(fit_2pp, "b")
sum_theta2 <- post_summary(fit_2pp, "theta")

# Extrai com segurança o índice numérico de nomes tipo "a[13]", "theta[200]" etc.
get_param_index <- function(param_vec, par) {
  # casa apenas strings EXATAS do tipo par[<inteiro>], ignorando outras linhas
  mat <- stringr::str_match(param_vec, sprintf("^%s\\[(\\d+)\\]$", par))
  idx <- suppressWarnings(as.integer(mat[, 2]))  # NA para os que não casarem
  idx
}

# Retorna data.frame filtrado para entradas válidas e ordenado por índice
order_param_df <- function(df, par) {
  idx <- get_param_index(df$Parameter, par)
  keep <- !is.na(idx)
  df2 <- df[keep, , drop = FALSE]
  idx2 <- idx[keep]
  df2[order(idx2), , drop = FALSE]
}

# Aplicar
sum_a1     <- order_param_df(sum_a1,     "a")
sum_a2     <- order_param_df(sum_a2,     "a")
sum_b1     <- order_param_df(sum_b1,     "b")
sum_b2     <- order_param_df(sum_b2,     "b")
sum_theta1 <- order_param_df(sum_theta1, "theta")
sum_theta2 <- order_param_df(sum_theta2, "theta")



## Constrói expressões T_k(rho[r]) com base em idx$rho_len
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



## Mapeia cada item ao testlet k (ou NA se independente)
testlet_of_item <- rep(NA_integer_, vI)
for (k in seq_along(idx_testlets)) {
  testlet_of_item[idx_testlets[[k]]] <- k
}

## Constrói rótulos como expression: independentes -> número; testlets -> T[k](i)
build_item_labels <- function(vI, testlet_of_item) {
  labs <- vector("list", vI)
  for (i in seq_len(vI)) {
    k <- testlet_of_item[i]
    labs[[i]] <- if (is.na(k)) {
      bquote(.(i))                # número simples como expressão
    } else {
      bquote(T[.(k)](.(i)))       # T_k(i)
    }
  }
  as.expression(labs)
}

item_labels <- build_item_labels(vI, testlet_of_item)


## Paleta cinza
scale_cols <- scale_color_manual(values = c("Testlet 2PP" = "black", "2PP" = "grey50"))

## ρ
p_rho <- ggplot(rho_data, aes(x = Rho, y = Estimate)) +
  geom_point(size = 2, color = "black") +
  geom_errorbar(aes(ymin = CI_Low, ymax = CI_High), width = 0.15, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal(base_size = 14) +
  scale_x_continuous(breaks = rho_data$Rho, labels = rho_labels) +
  labs(x = NULL, y = "Estimate")
ggsave(file.path(pathFit, "Fig_Rho_Testlet.pdf"), p_rho, width = 9, height = 3.6)

## a (discriminação)
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

## b (dificuldade)
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













