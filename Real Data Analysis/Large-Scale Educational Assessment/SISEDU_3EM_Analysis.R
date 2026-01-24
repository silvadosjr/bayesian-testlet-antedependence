## -------------------- 0) Paths and helper sources ----------------------------
root_local <- "~/GitHub/bayesian-testlet-antedependence"
path_project <- "~/GitHub/bayesian-testlet-antedependence/Real Data Analysis/Large-Scale Educational Assessment"
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
source(here('Real Data Analysis', "HelpersRealDataAnal.R"))


## Output directories
pathFit <- here(path_project, "Fit")
dir.create(pathFit, showWarnings = FALSE, recursive = TRUE)
stopifnot(dir.exists(pathFit))

saveFigures <- here(pathFit, "Figures")
dir.create(saveFigures, showWarnings = FALSE, recursive = TRUE)
stopifnot(dir.exists(saveFigures))


## ------------------- 2) Carregar dados (3ª série) ------------------

file_full <- here(path_project,'Binary_3EM_municipio_n25k.rds')
stopifnot(file.exists(file_full))
mYcNA <- readRDS(file_full)

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
pars_testlet <- c("a", "b", "theta", "rho_global")
pars_2pp     <- c("a", "b", "theta")

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

## Salvar objetos de ajuste
saveRDS(fit_testlet, file.path(pathFit, "ResultTestlet2PP_3EM_municipio25k.rds"))
saveRDS(fit_2pp,     file.path(pathFit, "Result2PP_3EM_municipio25k.rds"))



fit_testlet<-readRDS(file.path(pathFit,'ResultTestlet2PP_3EM_municipio25k.rds'))
fit_2pp<-readRDS(file.path(pathFit,"Result2PP_3EM_municipio25k.rds"))



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
  theme_minimal(base_size = 25) +
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












