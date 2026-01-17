corrigeItens <-function (respostas, gabarito)
{
  #  gabarito <- validarClave(respostas, gabarito)
  respostas[is.na(respostas)] <- 0
  output <- matrix(NA, nrow(respostas), ncol(respostas))
  colnames(output) <- colnames(respostas)
  rownames(output) <- rownames(respostas)
  for (i in 1:ncol(respostas)) {
    output[,i] <- (respostas[,i] == gabarito[i]) * 1
  }
  return(output)
}



rtrunc_norm_scalar <- function(mu, sd, lower, upper) {
  a <- (lower - mu)/sd
  b <- (upper - mu)/sd
  ua <- pnorm(a); ub <- pnorm(b)
  # proteção numérica
  if (ua >= ub) return(if (is.finite(lower)) lower else mu)
  u <- runif(1, ua, ub)
  mu + sd * qnorm(u)
}



# Um draw GHK para Z ~ N(mu, Sigma) com truncação componente a componente:
# lower[t] (tip. 0 ou -Inf), upper[t] (tip. Inf ou 0)
# Implementação sequencial via fatoração de Cholesky
ghk_draw_trunc_mvn <- function(mu, Sigma, lower, upper, max_try = 50L) {
  L <- tryCatch(chol(Sigma), error = function(e) NULL)
  if (is.null(L)) stop("GHK: Cholesky falhou (Sigma não SPD).")
  p <- length(mu)
  
  for (attempt in seq_len(max_try)) {
    z <- numeric(p)
    # simulamos e ~ N(0,I) com truncações transformadas
    # Z = mu + L %*% e  =>  e = solve(L, Z - mu)
    # No passo t, condicionamos em z[1..t-1]
    ok <- TRUE
    e <- numeric(p)
    for (t in 1:p) {
      # Média e variância condicionais de Z_t dado Z_{1..t-1} = z_{1..t-1}
      if (t == 1L) {
        m_t <- mu[1L]
        s_t <- sqrt(Sigma[1L, 1L])
      } else {
        # particiona Sigma:
        S11 <- Sigma[1:(t-1), 1:(t-1), drop = FALSE]
        S12 <- Sigma[1:(t-1), t, drop = FALSE]
        S21 <- Sigma[t, 1:(t-1), drop = FALSE]
        S22 <- Sigma[t, t]
        # condicional: Z_t | Z_{1..t-1}=z_{1..t-1} ~ N(m_t, s_t^2)
        # m_t = mu_t + S21 S11^{-1} (z1 - mu1)
        # s_t^2 = S22 - S21 S11^{-1} S12
        Sinv_11 <- tryCatch(solve(S11), error = function(e) NULL)
        if (is.null(Sinv_11)) { ok <- FALSE; break }
        delta   <- z[1:(t-1)] - mu[1:(t-1)]
        m_t     <- as.numeric(mu[t] + S21 %*% (Sinv_11 %*% delta))
        s2_t    <- as.numeric(S22 - S21 %*% (Sinv_11 %*% S12))
        s_t     <- sqrt(max(s2_t, 1e-12))
      }
      lo <- lower[t]; up <- upper[t]
      z[t] <- rtrunc_norm_scalar(m_t, s_t, lo, up)
      # continua
    }
    if (ok) return(z)
  }
  stop("GHK: não conseguiu amostrar Z consistente dentro de max_try.")
}



# log-score aumentado: soma de log N(Z_block; mu_block, Ck) + log N(Z_ind; mu_i, 1)
# - para Y_new: Z_new é um draw consistente (via GHK / trunc escalar)
# - para Y_rep: Z_rep é o próprio draw do modelo (sem truncação)
diag_logscore_aug_arrays <- function(params, theta_vec, Y, idx_testlets,
                                     Z_given_Y = TRUE) {
  # params: list(a,b,C_list)
  # theta_vec: comprimento N
  # Y: N x I
  N <- nrow(Y); I <- ncol(Y)
  a <- params$a; b <- params$b
  
  ll <- 0
  
  for (j in 1:N) {
    # --- independentes ---
    ind_items <- setdiff(seq_len(I), unlist(idx_testlets, use.names = FALSE))
    if (length(ind_items)) {
      mu_ind <- a[ind_items] * (theta_vec[j] - b[ind_items])
      if (Z_given_Y) {
        # amostra Z|Y, escalar por escalar
        for (ii in ind_items) {
          mui <- a[ii] * (theta_vec[j] - b[ii])
          if (Y[j, ii] == 1L) {
            zi <- rtrunc_norm_scalar(mui, 1, 0, Inf)
          } else {
            zi <- rtrunc_norm_scalar(mui, 1, -Inf, 0)
          }
          ll <- ll + dnorm(zi, mean = mui, sd = 1, log = TRUE)
        }
      } else {
        # Z rep já vem direto do modelo — aqui, por simetria, quem chama
        # deve passar Z já sorteado (ver função HPC abaixo).
        stop("diag_logscore_aug_arrays: modo Z_given_Y=FALSE deve receber Z externo.")
      }
    }
    
    # --- testlets ---
    if (!is.null(params$C_list)) {
      for (k in seq_along(idx_testlets)) {
        items_k <- idx_testlets[[k]]
        if (!length(items_k)) next
        Ck <- params$C_list[[k]]
        muk <- a[items_k] * (theta_vec[j] - b[items_k])
        
        if (Z_given_Y) {
          # limites por item: (0,Inf) se Y=1; (-Inf,0] se Y=0
          Lk <- length(items_k)
          lower <- ifelse(Y[j, items_k] == 1L, 0, -Inf)
          upper <- ifelse(Y[j, items_k] == 1L, Inf, 0)
          zk <- ghk_draw_trunc_mvn(muk, Ck, lower, upper)
          ll <- ll + mvtnorm::dmvnorm(zk, mean = muk, sigma = Ck, log = TRUE)
        } else {
          stop("diag_logscore_aug_arrays: modo Z_given_Y=FALSE deve receber Z externo.")
        }
      }
    }
  }
  ll
}


# ext: lista vinda de rstan::extract (ext1 para Testlet, ext2 para 2PP)
# build_params_fun: uma das funções 'build_params_*_from_arrays'
# idx_testlets, ind_items: como no restante do seu código
# ATENÇÃO: para d_rep, vamos gerar Z_rep diretamente do modelo,
#          e reutilizar a MESMA métrica aumentada (log N(Z; mu, Sigma)).
compute_hpc_aug_arrays <- function(ext, build_params_fun, Y_new, idx_testlets,
                                   ind_items = NULL, R_eval = 250, seed = 123,
                                   # controle GHK:
                                   max_try = 50L) {
  set.seed(seed)
  
  R_avail <- nrow(ext$a)
  if (is.null(R_avail) || R_avail < 1L) stop("extract: sem draws.")
  draw_ids <- sample.int(R_avail, size = min(R_eval, R_avail), replace = FALSE)
  
  N_new <- nrow(Y_new); I <- ncol(Y_new)
  if (is.null(ind_items)) ind_items <- integer(0)
  
  d_rep <- numeric(length(draw_ids))
  d_new <- numeric(length(draw_ids))
  
  for (t in seq_along(draw_ids)) {
    dr <- draw_ids[t]
    params <- build_params_fun(dr, ext)
    
    # --- thetas novos da priori ---
    theta_new <- rnorm(n = N_new, mean = 0, sd = 1)
    
    # --- d_new: usa Z|Y por GHK (augmented observed) ---
    d_new[t] <- diag_logscore_aug_arrays(params, theta_new, Y_new, idx_testlets,
                                         Z_given_Y = TRUE)
    
    # --- d_rep: gera Z_rep do modelo e Y_rep por threshold ---
    # Independentes
    Z_rep <- matrix(0, nrow = N_new, ncol = I)
    a <- params$a; b <- params$b
    
    ind_out <- setdiff(seq_len(I), unlist(idx_testlets, use.names = FALSE))
    if (length(ind_out)) {
      MU_ind <- outer(theta_new, a[ind_out], "*") -
        matrix(rep(a[ind_out] * b[ind_out], each = N_new), nrow = N_new)
      # Z ~ N(mu, 1)
      Z_rep[, ind_out] <- matrix(stats::rnorm(length(MU_ind), mean = MU_ind, sd = 1),
                                 nrow = N_new)
    }
    
    # Testlets
    if (!is.null(params$C_list)) {
      for (k in seq_along(idx_testlets)) {
        items_k <- idx_testlets[[k]]
        if (!length(items_k)) next
        Ck <- params$C_list[[k]]
        for (j in seq_len(N_new)) {
          muk <- a[items_k] * (theta_new[j] - b[items_k])
          Z_rep[j, items_k] <- as.numeric(mvtnorm::rmvnorm(1, mean = muk, sigma = Ck))
        }
      }
    }
    
    # Y_rep apenas se você quiser guardar; para o logscore aumentado basta Z_rep
    # Y_rep <- (Z_rep > 0L)
    
    # Agora avalie o logscore aumentado em Z_rep:
    ll <- 0
    # independentes
    if (length(ind_out)) {
      mu_ind_all <- outer(theta_new, a[ind_out], "*") -
        matrix(rep(a[ind_out]*b[ind_out], each = N_new), nrow = N_new)
      ll <- ll + sum(dnorm(Z_rep[, ind_out], mean = mu_ind_all, sd = 1, log = TRUE))
    }
    # blocos
    if (!is.null(params$C_list)) {
      for (k in seq_along(idx_testlets)) {
        items_k <- idx_testlets[[k]]
        if (!length(items_k)) next
        Ck <- params$C_list[[k]]
        for (j in seq_len(N_new)) {
          muk <- a[items_k] * (theta_new[j] - b[items_k])
          ll  <- ll + mvtnorm::dmvnorm(Z_rep[j, items_k], mean = muk, sigma = Ck, log = TRUE)
        }
      }
    }
    d_rep[t] <- ll
  }
  
  list(p_HPC = mean(d_rep >= d_new), d_rep = d_rep, d_new = d_new)
}




## 6) Helpers: Toeplitz e log-prob de um testlet
make_toeplitz <- function(rhos) {
  L <- length(rhos) + 1L
  C <- diag(L)
  for (r in seq_len(L-1L)) {
    C[row(C) == col(C) - r] <- rhos[r]
    C[row(C) == col(C) + r] <- rhos[r]
  }
  C
}

make_rho_index <- function(nk, is_HU) {
  stopifnot(length(nk) == length(is_HU))
  K <- length(nk)
  rho_len <- integer(K)
  for (k in seq_len(K)) {
    if (is.na(is_HU[k])) {
      rho_len[k] <- 0L
    } else if (is_HU[k]) {
      rho_len[k] <- 1L
    } else {
      rho_len[k] <- nk[k] - 1L
    }
  }
  rho_start <- integer(K)
  if (K > 0) {
    rho_start[1] <- if (rho_len[1] > 0) 1L else 1L  # começa em 1
    for (k in 2:K) {
      rho_start[k] <- rho_start[k-1] + max(rho_len[k-1], 0L)
    }
  }
  list(rho_len = rho_len, rho_start = rho_start)
}



logprob_testlet <- function(y_vec, mu_vec, C) {
  # y=1 -> Z>0; y=0 -> Z<=0; probit multivariado
  L <- length(y_vec)
  lower <- ifelse(y_vec == 1, 0, -Inf)
  upper <- ifelse(y_vec == 1, Inf, 0)
  # pmvnorm retorna probabilidade; proteger contra underflow
  p <- mvtnorm::pmvnorm(lower=lower, upper=upper, mean=mu_vec, sigma=C)[1]
  if (!is.finite(p) || p <= 0) return(log(1e-300))
  log(p)
}


## 8) Diagnóstico 1: log-score médio no holdout
diag_logscore <- function(params, theta, Y, idx_testlets, ind_items) {
  # params: list(a, b, C_list) no Testlet; no 2PP, C_list=NULL
  # theta: vetor N (amostrado de N(0,1))
  N <- nrow(Y); I <- ncol(Y)
  a <- params$a; b <- params$b
  # por-aluno, acumula log-prob de todos os itens
  ll <- 0
  # Itens independentes (modelo 2PP OU marginais no Testlet)
  # para 2PP: todos ind; para Testlet: apenas 'ind_items'
  if (!is.null(ind_items) && length(ind_items) > 0) {
    MU_ind <- outer(theta, a[ind_items], "*") - matrix(rep(a[ind_items]*b[ind_items], each=N), nrow=N)
    P_ind  <- pnorm(MU_ind)
    P_ind  <- pmin(pmax(P_ind, 1e-12), 1-1e-12)
    Yi     <- as.matrix(Y[, ind_items, drop=FALSE])
    ll <- ll + sum( Yi*log(P_ind) + (1-Yi)*log(1-P_ind) )
  }
  # Testlets
  if (!is.null(idx_testlets)) {
    for (k in seq_along(idx_testlets)) {
      items_k <- idx_testlets[[k]]
      Lk <- length(items_k)
      if (Lk == 0) next
      if (is.null(params$C_list)) {
        # 2PP: itens tratados como independentes
        MU_k <- outer(theta, a[items_k], "*") - matrix(rep(a[items_k]*b[items_k], each=N), nrow=N)
        P_k  <- pnorm(MU_k)
        P_k  <- pmin(pmax(P_k, 1e-12), 1-1e-12)
        Yk   <- as.matrix(Y[, items_k, drop=FALSE])
        ll   <- ll + sum( Yk*log(P_k) + (1-Yk)*log(1-P_k) )
      } else {
        # Testlet-2PP: probit multivariado
        Ck <- params$C_list[[k]]
        for (j in seq_len(N)) {
          mu_jk <- a[items_k]*(theta[j] - b[items_k])
          y_jk  <- as.numeric(Y[j, items_k])
          ll <- ll + logprob_testlet(y_jk, mu_jk, Ck)
        }
      }
    }
  }
  # retorna NEG-log-verossimilhança média
  return( - ll / (N*I) )
}

# === Diagnóstico 2 (max |corr| intra-testlet), versão estável/rápida ===
diag_residQ3 <- function(params, theta, Y, idx_testlets, eps = 1e-9) {
  a <- params$a; b <- params$b
  # mu = a*(theta - b) -> p = Phi(mu)
  MU <- outer(theta, a, "*") - matrix(rep(a*b, each = nrow(Y)), nrow = nrow(Y))
  P  <- pnorm(MU)
  
  # clipping suave evita var ~= 0 quando p ~ 0 ou 1
  P  <- pmin(pmax(P, eps), 1 - eps)
  R  <- Y - P  # resíduos preditivos
  
  max_abs <- 0
  for (k in seq_along(idx_testlets)) {
    items_k <- idx_testlets[[k]]
    if (length(items_k) < 2) next
    
    Rk <- as.matrix(R[, items_k, drop = FALSE])
    
    # centraliza por coluna (item) para remover nível médio do resíduo
    # equivalente a cor(Rk) mas com guardas para sd ~ 0
    Rk <- scale(Rk, center = TRUE, scale = FALSE)
    
    # covariância e conversão p/ correlação com proteção
    S  <- crossprod(Rk) / (nrow(Rk) - 1)        # mais rápido que cov()
    sd <- sqrt(pmax(diag(S), eps))
    Ck <- S / (sd %o% sd)
    
    Ck[!is.finite(Ck)] <- 0
    val <- max(abs(Ck[upper.tri(Ck)]), na.rm = TRUE)
    if (is.finite(val)) max_abs <- max(max_abs, val)
  }
  return(max_abs)
}

diag_residQ3_vec <- function(params, theta, Y, idx_testlets, eps = 1e-9) {
  a <- params$a; b <- params$b
  MU <- outer(theta, a, "*") - matrix(rep(a*b, each = nrow(Y)), nrow = nrow(Y))
  P  <- pnorm(MU)
  P  <- pmin(pmax(P, eps), 1 - eps)
  R  <- Y - P
  
  K <- length(idx_testlets)
  out <- rep(NA_real_, K)
  
  for (k in seq_len(K)) {
    items_k <- idx_testlets[[k]]
    if (length(items_k) < 2) next
    
    Rk <- as.matrix(R[, items_k, drop = FALSE])
    Rk <- scale(Rk, center = TRUE, scale = FALSE)
    
    S  <- crossprod(Rk) / (nrow(Rk) - 1)
    sd <- sqrt(pmax(diag(S), eps))
    Ck <- S / (sd %o% sd)
    Ck[!is.finite(Ck)] <- 0
    
    out[k] <- max(abs(Ck[upper.tri(Ck)]), na.rm = TRUE)
  }
  out
}

## 10) Construtores de parâmetros por draw
build_params_testlet <- function(draw_row) {
  a <- as.numeric(draw_row[a_names]); names(a) <- NULL
  b <- as.numeric(draw_row[b_names]); names(b) <- NULL
  C_list <- vector("list", K)
  for (k in seq_len(K)) {
    Lk <- nk[k]
    if (Lk == 2L) {
      rho_k <- as.numeric(draw_row[rho_names_list[[k]]])
      C_list[[k]] <- make_toeplitz(rhos = rho_k)
    } else {
      rho_k <- as.numeric(draw_row[rho_names_list[[k]]])
      C_list[[k]] <- make_toeplitz(rhos = rho_k)
    }
  }
  list(a=a, b=b, C_list=C_list)
}

build_params_2pp <- function(draw_row) {
  a <- as.numeric(draw_row[a_names]); names(a) <- NULL
  b <- as.numeric(draw_row[b_names]); names(b) <- NULL
  list(a=a, b=b, C_list=NULL)
}

## 11) HPC (p-valor) para um diagnóstico genérico
compute_hpc <- function(post_df, build_params_fun, Y_new, idx_testlets, ind_items=NULL,
                        diag = c("logscore","residQ3"), R_eval = 250, seed=123) {
  set.seed(seed)
  R_avail <- nrow(post_df)
  idx <- sample.int(R_avail, size = min(R_eval, R_avail), replace = FALSE)
  d_rep <- numeric(length(idx))
  d_new <- numeric(length(idx))
  for (t in seq_along(idx)) {
    dr <- idx[t]
    params <- build_params_fun(post_df[dr, , drop=FALSE])
    # θ dos alunos NOVOS (holdout): amostra da priori N(0,1)
    theta_new <- rnorm(n = nrow(Y_new), mean = 0, sd = 1)
    if (diag[1] == "logscore") {
      # d(y) = neg-logverossim média; y_rep simulado do modelo
      d_new[t] <- diag_logscore(params, theta_new, Y_new, idx_testlets, ind_items)
      # simula y_rep:
      # (a) indep items
      I <- ncol(Y_new)
      Y_rep <- matrix(0L, nrow=nrow(Y_new), ncol=I)
      a <- params$a; b <- params$b
      # independentes
      if (!is.null(ind_items) && length(ind_items) > 0) {
        MU_ind <- outer(theta_new, a[ind_items], "*") - matrix(rep(a[ind_items]*b[ind_items], each=nrow(Y_new)), nrow=nrow(Y_new))
        P_ind <- pnorm(MU_ind)
        Y_rep[, ind_items] <- matrix(rbinom(length(P_ind), 1L, P_ind), nrow=nrow(Y_new))
      } else {
        ind_items <- integer(0)
      }
      # testlets
      if (!is.null(params$C_list)) {
        for (k in seq_along(idx_testlets)) {
          items_k <- idx_testlets[[k]]
          if (length(items_k) == 0) next
          Ck <- params$C_list[[k]]
          Lk <- length(items_k)
          for (j in seq_len(nrow(Y_new))) {
            mu_jk <- a[items_k]*(theta_new[j] - b[items_k])
            # simula Z ~ N(mu, C) e aplica threshold
            Z <- as.numeric(mvtnorm::rmvnorm(n=1, mean=mu_jk, sigma=Ck))
            Y_rep[j, items_k] <- as.integer(Z > 0)
          }
        }
      } else {
        # 2PP: tudo independente
        MU_all <- outer(theta_new, a, "*") - matrix(rep(a*b, each=nrow(Y_new)), nrow=nrow(Y_new))
        P_all <- pnorm(MU_all)
        Y_rep <- matrix(rbinom(length(P_all), 1L, P_all), nrow=nrow(Y_new))
      }
      d_rep[t] <- diag_logscore(params, theta_new, Y_rep, idx_testlets, ind_items)
    } else if (diag[1] == "residQ3") {
      # usamos o mesmo theta_new para y_new e y_rep (realized-style)
      d_new[t] <- diag_residQ3(params, theta_new, Y_new, idx_testlets)
      # simula y_rep (como acima)
      I <- ncol(Y_new)
      Y_rep <- matrix(0L, nrow=nrow(Y_new), ncol=I)
      a <- params$a; b <- params$b
      if (!is.null(params$C_list)) {
        # independentes fora dos testlets
        if (!is.null(ind_items) && length(ind_items) > 0) {
          MU_ind <- outer(theta_new, a[ind_items], "*") - matrix(rep(a[ind_items]*b[ind_items], each=nrow(Y_new)), nrow=nrow(Y_new))
          P_ind <- pnorm(MU_ind)
          Y_rep[, ind_items] <- matrix(rbinom(length(P_ind), 1L, P_ind), nrow=nrow(Y_new))
        }
        for (k in seq_along(idx_testlets)) {
          items_k <- idx_testlets[[k]]
          if (length(items_k) == 0) next
          Ck <- params$C_list[[k]]
          for (j in seq_len(nrow(Y_new))) {
            mu_jk <- a[items_k]*(theta_new[j] - b[items_k])
            Z <- as.numeric(mvtnorm::rmvnorm(n=1, mean=mu_jk, sigma=Ck))
            Y_rep[j, items_k] <- as.integer(Z > 0)
          }
        }
      } else {
        MU_all <- outer(theta_new, a, "*") - matrix(rep(a*b, each=nrow(Y_new)), nrow=nrow(Y_new))
        P_all <- pnorm(MU_all)
        Y_rep <- matrix(rbinom(length(P_all), 1L, P_all), nrow=nrow(Y_new))
      }
      d_rep[t] <- diag_residQ3(params, theta_new, Y_rep, idx_testlets)
    } else {
      stop("Diagnóstico não reconhecido.")
    }
  }
  # p_HPC = Pr{ d(y_rep) >= d(y_new) }
  p_val <- mean(d_rep >= d_new)
  list(p_HPC = p_val, d_rep = d_rep, d_new = d_new)
}


## --- Univariada truncada com sd != 1 (por inversão) ---
rtrunc_norm_1d_sd <- function(mu, sd, lower, upper) {
  # Retorna 1 amostra de N(mu, sd^2) truncada em [lower, upper]
  a <- pnorm(lower, mean = mu, sd = sd)
  b <- pnorm(upper, mean = mu, sd = sd)
  # Casos degenerados: se a==b, joga no limite mais próximo
  if (!is.finite(a)) a <- 0
  if (!is.finite(b)) b <- 1
  if (b <= a) {
    # devolve o limite que estiver mais próximo de mu
    return(if (abs(lower - mu) < abs(upper - mu)) lower else upper)
  }
  u <- runif(1, a, b)
  qnorm(u, mean = mu, sd = sd)
}

## --- Gibbs para MVN truncada (dim <= 3), com algumas iterações de burn-in ---
rtrunc_mvn_gibbs <- function(mu, Sigma, lower, upper, n_iter = 300, burn = 100) {
  L <- length(mu)
  # Chute inicial: projeta mu para dentro do retângulo [lower, upper]
  z <- pmin(pmax(mu, lower), upper)
  # Pré-cálculo da inversa/soluções
  for (it in seq_len(n_iter)) {
    for (j in seq_len(L)) {
      idx_other <- setdiff(seq_len(L), j)
      sigma_jj <- Sigma[j, j]
      if (length(idx_other) == 0L) {
        mu_cond <- mu[j]
        sd_cond <- sqrt(sigma_jj)
      } else {
        Sigma_jO <- Sigma[j, idx_other, drop=FALSE]
        Sigma_OO <- Sigma[idx_other, idx_other, drop=FALSE]
        # solve(Sigma_OO) é baratíssimo em L<=3
        invSigma_OO <- solve(Sigma_OO)
        mu_cond <- as.numeric(mu[j] + Sigma_jO %*% invSigma_OO %*% (z[idx_other] - mu[idx_other]))
        var_cond <- sigma_jj - as.numeric(Sigma_jO %*% invSigma_OO %*% t(Sigma_jO))
        sd_cond <- sqrt(max(var_cond, 1e-12))
      }
      z[j] <- rtrunc_norm_1d_sd(mu_cond, sd_cond, lower[j], upper[j])
    }
  }
  z
}

## --- MVN truncada “inteligente”: tenta tmvtnorm -> rejeição -> Gibbs ---
rtrunc_mvn_small <- function(mu, Sigma, y, max_iter = 5000L, eps = 1e-8) {
  L <- length(mu)
  if (L == 0L) return(numeric(0))
  # Intervalos por y (com epsilon)
  lower <- ifelse(y == 1L, 0 + eps, -Inf)
  upper <- ifelse(y == 1L, Inf, 0 - eps)
  
  # 1) Se pacote tmvtnorm estiver disponível, use-o (rápido e estável)
  if (requireNamespace("tmvtnorm", quietly = TRUE)) {
    z <- try(tmvtnorm::rtmvnorm(
      n = 1, mean = mu, sigma = Sigma, lower = lower, upper = upper, algorithm = "gibbs"
    ), silent = TRUE)
    if (!inherits(z, "try-error")) return(as.numeric(z))
  }
  # 2) Rejeição rápida (boa quando probabilidade não é minúscula)
  for (t in 1:max_iter) {
    z <- as.numeric(mvtnorm::rmvnorm(1, mean = mu, sigma = Sigma))
    if (all((y == 1L & z > 0) | (y == 0L & z <= 0))) return(z)
  }
  # 3) Fallback: Gibbs por condicionais univariadas truncadas
  z <- rtrunc_mvn_gibbs(mu, Sigma, lower, upper, n_iter = 400, burn = 150)
  z
}


## ---- Log-densidade de Z sob Normal univariada/multivariada (sem constantes irrelevantes para comparação) ----
log_dnorm1 <- function(z, mu) {
  # log phi(z; mu, 1)
  -0.5 * (z - mu)^2 - 0.5 * log(2*pi)
}

log_dmvnorm_fast <- function(z, mu, Sigma, cholSigma = NULL) {
  # log phi_L(z; mu, Sigma)
  # usa fatoração de Cholesky se fornecida
  L <- length(mu)
  if (is.null(cholSigma)) cholSigma <- chol(Sigma)
  x <- backsolve(cholSigma, z - mu, transpose = TRUE)
  quad <- sum(x^2)
  const <- -0.5 * L * log(2*pi) - sum(log(diag(cholSigma)))
  const - 0.5 * quad
}







## ---- Diagnóstico aumentado: neg-log-densidade média de Z (global) ----
diag_logscore_aug <- function(params, theta, Y, idx_testlets, ind_items) {
  # params: list(a, b, C_list) no Testlet; no 2PP, C_list=NULL
  N <- nrow(Y); I <- ncol(Y)
  a <- params$a; b <- params$b
  ll <- 0
  
  # Independentes
  items_ind <- ind_items
  if (is.null(items_ind)) items_ind <- integer(0)
  if (length(items_ind) > 0) {
    MU_ind <- outer(theta, a[items_ind], "*") - matrix(rep(a[items_ind]*b[items_ind], each=N), nrow=N)
    for (j in seq_len(N)) {
      for (ii in seq_along(items_ind)) {
        i_idx <- items_ind[ii]
        zji <- rtrunc_norm_1d(mu = MU_ind[j, ii], y = Y[j, i_idx])
        ll <- ll + log_dnorm1(zji, MU_ind[j, ii])
      }
    }
  }
  
  # Testlets
  if (!is.null(params$C_list)) {
    for (k in seq_along(idx_testlets)) {
      items_k <- idx_testlets[[k]]
      if (length(items_k) == 0) next
      Ck <- params$C_list[[k]]
      cholCk <- chol(Ck)
      for (j in seq_len(N)) {
        mu_jk <- a[items_k]*(theta[j] - b[items_k])
        y_jk  <- as.integer(Y[j, items_k])
        z_jk  <- rtrunc_mvn_small(mu = mu_jk, Sigma = Ck, y = y_jk)
        ll <- ll + log_dmvnorm_fast(z_jk, mu_jk, Ck, cholCk)
      }
    }
  } else {
    # 2PP: todos independentes -> já tratado acima
    # (se existir item não em ind_items, trate como independente aqui também)
    rest <- setdiff(seq_len(I), items_ind)
    if (length(rest) > 0) {
      MU_rest <- outer(theta, a[rest], "*") - matrix(rep(a[rest]*b[rest], each=N), nrow=N)
      for (j in seq_len(N)) {
        for (ii in seq_along(rest)) {
          i_idx <- rest[ii]
          zji <- rtrunc_norm_1d(mu = MU_rest[j, ii], y = Y[j, i_idx])
          ll <- ll + log_dnorm1(zji, MU_rest[j, ii])
        }
      }
    }
  }
  
  # retorna NEG-log-média por resposta
  return( - ll / (N*I) )
}


## ---- Versão aumentada do compute_hpc (apenas o caso diag="logscore_aug") ----
compute_hpc_aug <- function(post_df, build_params_fun, Y_new, idx_testlets, ind_items=NULL,
                            R_eval = 150, seed=2025) {
  set.seed(seed)
  R_avail <- nrow(post_df)
  idx <- sample.int(R_avail, size = min(R_eval, R_avail), replace = FALSE)
  d_rep <- numeric(length(idx))
  d_new <- numeric(length(idx))
  
  for (t in seq_along(idx)) {
    dr <- idx[t]
    params <- build_params_fun(post_df[dr, , drop=FALSE])
    a <- params$a; b <- params$b
    N <- nrow(Y_new); I <- ncol(Y_new)
    
    # Habilidades no holdout (realized style)
    theta_new <- rnorm(N, 0, 1)
    
    # ---------- d_new: Z|Y (holdout) ----------
    d_new[t] <- diag_logscore_aug(params, theta_new, Y_new, idx_testlets, ind_items)
    
    # ---------- y_rep, z_rep ----------
    # Gera Z_rep do modelo e define Y_rep pelo sinal
    Z_rep <- matrix(0, nrow=N, ncol=I)
    Y_rep <- matrix(0L, nrow=N, ncol=I)
    
    # Independentes
    items_ind <- ind_items
    if (is.null(items_ind)) items_ind <- integer(0)
    
    if (length(items_ind) > 0) {
      MU_ind <- outer(theta_new, a[items_ind], "*") - matrix(rep(a[items_ind]*b[items_ind], each=N), nrow=N)
      Z_rep[, items_ind] <- matrix(rnorm(N*length(items_ind), mean = as.vector(MU_ind), sd = 1),
                                   nrow=N, byrow=FALSE)
      Y_rep[, items_ind] <- (Z_rep[, items_ind] > 0L)
    }
    
    # Testlets
    if (!is.null(params$C_list)) {
      for (k in seq_along(idx_testlets)) {
        items_k <- idx_testlets[[k]]
        if (length(items_k) == 0) next
        Ck <- params$C_list[[k]]
        for (j in seq_len(N)) {
          mu_jk <- a[items_k]*(theta_new[j] - b[items_k])
          z_jk  <- as.numeric(mvtnorm::rmvnorm(1, mean=mu_jk, sigma=Ck))
          Z_rep[j, items_k] <- z_jk
          Y_rep[j, items_k] <- as.integer(z_jk > 0)
        }
      }
    } else {
      # 2PP: todos independentes (garante também para itens fora de ind_items)
      rest <- setdiff(seq_len(I), items_ind)
      if (length(rest) > 0) {
        MU_rest <- outer(theta_new, a[rest], "*") - matrix(rep(a[rest]*b[rest], each=N), nrow=N)
        Z_rep[, rest] <- matrix(rnorm(N*length(rest), mean = as.vector(MU_rest), sd = 1),
                                nrow=N, byrow=FALSE)
        Y_rep[, rest] <- (Z_rep[, rest] > 0L)
      }
    }
    
    # ---------- d_rep: log-densidade de Z_rep ----------
    # Reutiliza a mesma função mas passando Y_rep apenas para interface (não é usado aqui)
    d_rep[t] <- {
      # calculamos a neg-log-média de Z_rep sob o mesmo params
      # versão "rápida": copiar diag_logscore_aug mas sem amostrar Z|Y; usar Z_rep já gerado
      N <- nrow(Y_rep); I <- ncol(Y_rep)
      ll <- 0
      
      # Independentes
      all_ind <- union(items_ind, if (is.null(params$C_list)) setdiff(seq_len(I), items_ind) else integer(0))
      if (length(all_ind) > 0) {
        MU_ind <- outer(theta_new, a[all_ind], "*") - matrix(rep(a[all_ind]*b[all_ind], each=N), nrow=N)
        for (j in seq_len(N)) {
          for (ii in seq_along(all_ind)) {
            i_idx <- all_ind[ii]
            ll <- ll + log_dnorm1(Z_rep[j, i_idx], MU_ind[j, ii])
          }
        }
      }
      
      # Testlets
      if (!is.null(params$C_list)) {
        for (k in seq_along(idx_testlets)) {
          items_k <- idx_testlets[[k]]
          if (length(items_k) == 0) next
          Ck <- params$C_list[[k]]
          cholCk <- chol(Ck)
          for (j in seq_len(N)) {
            mu_jk <- a[items_k]*(theta_new[j] - b[items_k])
            z_jk  <- Z_rep[j, items_k]
            ll <- ll + log_dmvnorm_fast(z_jk, mu_jk, Ck, cholCk)
          }
        }
      }
      
      - ll / (N*I)
    }
  }
  
  p_val <- mean(d_rep >= d_new)
  list(p_HPC = p_val, d_rep = d_rep, d_new = d_new)
}

# Retorna matriz LxL:
# - length(rhos)==0  -> identidade (sem dependência)
# - length(rhos)==1  -> HU (uma única correlação)
# - length(rhos)==L-1-> HT (Toeplitz por defasagem)
make_R_block <- function(L, rhos) {
  R <- diag(L)
  if (L <= 1L || length(rhos) == 0L) return(R)
  if (length(rhos) == 1L) {
    R[row(R) != col(R)] <- rhos[1]
    return(R)
  }
  if (length(rhos) == (L - 1L)) {
    for (d in 1:(L - 1L)) {
      R[row(R) == col(R) + d] <- rhos[d]
      R[row(R) == col(R) - d] <- rhos[d]
    }
    return(R)
  }
  stop("make_R_block: comprimento de 'rhos' incompatível com L.")
}


# draw_id: índice do draw (linha dos arrays extraídos)
# ext     : lista do extract() (ext1 para Testlet; ext2 para 2PP)
# idx     : lista com $rho_start, $rho_len
# nk      : tamanhos dos testlets
# K       : número de testlets
build_params_testlet_from_arrays <- function(draw_id, ext, idx, nk, K) {
  a <- unname(ext$a[draw_id, ])
  b <- unname(ext$b[draw_id, ])
  
  C_list <- vector("list", K)
  for (k in seq_len(K)) {
    Lk  <- nk[k]
    len <- idx$rho_len[k]
    if (len == 0L) {
      C_list[[k]] <- make_R_block(Lk, numeric(0))
    } else {
      j0   <- idx$rho_start[k]
      jend <- j0 + len - 1L
      rhos <- unname(ext$rho_global[draw_id, j0:jend])
      C_list[[k]] <- make_R_block(Lk, rhos)
    }
  }
  list(a = a, b = b, C_list = C_list)
}

build_params_2pp_from_arrays <- function(draw_id, ext) {
  a <- unname(ext$a[draw_id, ])
  b <- unname(ext$b[draw_id, ])
  list(a = a, b = b, C_list = NULL)
}


# === Geração vetorizada de Y_rep dentro do compute_hpc_arrays ===
# Helper: simula bloco de um testlet de uma vez (sem loop em j)
.simulate_testlet_block <- function(theta_vec, a, b, items_k, Ck) {
  n  <- length(theta_vec); Ik <- length(items_k)
  mu_mat <- outer(theta_vec, a[items_k], "*") - matrix(rep(a[items_k]*b[items_k], each = n), nrow = n)
  
  # fator de Cholesky da correlação do testlet (inferior)
  L <- chol(Ck)  # default retorna superior; usamos t(L) para inferior
  # ruído ~ MN(0, I_n ⊗ Ck): G ~ N(0,1)  ->  G %*% t(L) tem cov Ck
  G <- matrix(rnorm(n * Ik), n, Ik)
  Z <- mu_mat + G %*% t(L)
  
  out <- matrix(0L, nrow = n, ncol = Ik)
  out[] <- as.integer(Z > 0)
  out
}

compute_hpc_arrays <- function(ext, build_params_fun, Y_new, idx_testlets, ind_items = NULL,
                               diag = c("logscore","residQ3"),
                               by_testlet = FALSE,     # << NOVO
                               R_eval = 250, seed = 123, ...) {
  set.seed(seed)
  diag <- match.arg(diag)
  
  R_avail <- nrow(ext$a)
  if (is.null(R_avail) || R_avail < 1L) stop("extract: sem draws.")
  
  draw_ids <- sample.int(R_avail, size = min(R_eval, R_avail), replace = FALSE)
  
  N_new <- nrow(Y_new); I <- ncol(Y_new)
  if (is.null(ind_items)) ind_items <- integer(0)
  
  # alocação: se by_testlet=TRUE e diag="residQ3", guardamos matriz R_eval x K
  K <- length(idx_testlets)
  if (diag == "residQ3" && by_testlet) {
    d_rep <- matrix(NA_real_, nrow = length(draw_ids), ncol = K)
    d_new <- rep(NA_real_, K)
  } else {
    d_rep <- numeric(length(draw_ids))
    d_new <- numeric(length(draw_ids))
  }
  
  # helper p/ gerar blocos de testlet sem loop em j
  .simulate_testlet_block <- function(theta_vec, a, b, items_k, Ck) {
    n  <- length(theta_vec); Ik <- length(items_k)
    mu_mat <- outer(theta_vec, a[items_k], "*") -
      matrix(rep(a[items_k]*b[items_k], each = n), nrow = n)
    L <- chol(Ck)
    G <- matrix(rnorm(n * Ik), n, Ik)
    Z <- mu_mat + G %*% t(L)
    out <- matrix(0L, nrow = n, ncol = Ik)
    out[] <- as.integer(Z > 0)
    out
  }
  
  for (t in seq_along(draw_ids)) {
    dr <- draw_ids[t]
    params <- build_params_fun(dr, ext, ...)
    theta_new <-rnorm(n = N_new, mean = 0, sd = 1)
    
    if (diag == "logscore") {
      d_new[t] <- diag_logscore(params, theta_new, Y_new, idx_testlets, ind_items)
      
      # --- gera Y_rep ---
      Y_rep <- matrix(0L, nrow = N_new, ncol = I)
      a <- params$a; b <- params$b
      if (length(ind_items) > 0) {
        MU_ind <- outer(theta_new, a[ind_items], "*") -
          matrix(rep(a[ind_items]*b[ind_items], each = N_new), nrow = N_new)
        P_ind <- pnorm(MU_ind)
        Y_rep[, ind_items] <- matrix(rbinom(length(P_ind), 1L, P_ind), nrow = N_new)
      }
      if (!is.null(params$C_list)) {
        for (k in seq_along(idx_testlets)) {
          items_k <- idx_testlets[[k]]
          if (!length(items_k)) next
          Ck <- params$C_list[[k]]
          Y_rep[, items_k] <- .simulate_testlet_block(theta_new, a, b, items_k, Ck)
        }
      } else {
        MU_all <- outer(theta_new, a, "*") - matrix(rep(a*b, each = N_new), nrow = N_new)
        P_all <- pnorm(MU_all)
        Y_rep <- matrix(rbinom(length(P_all), 1L, P_all), nrow = N_new)
      }
      d_rep[t] <- diag_logscore(params, theta_new, Y_rep, idx_testlets, ind_items)
      
    } else { # residQ3
      if (by_testlet) {
        d_new <- diag_residQ3_vec(params, theta_new, Y_new, idx_testlets)
      } else {
        d_new[t] <- diag_residQ3(params, theta_new, Y_new, idx_testlets)
      }
      
      # --- gera Y_rep ---
      Y_rep <- matrix(0L, nrow = N_new, ncol = I)
      a <- params$a; b <- params$b
      
      if (!is.null(params$C_list)) {
        if (length(ind_items) > 0) {
          MU_ind <- outer(theta_new, a[ind_items], "*") -
            matrix(rep(a[ind_items]*b[ind_items], each = N_new), nrow = N_new)
          P_ind <- pnorm(MU_ind)
          Y_rep[, ind_items] <- matrix(rbinom(length(P_ind), 1L, P_ind), nrow = N_new)
        }
        for (k in seq_along(idx_testlets)) {
          items_k <- idx_testlets[[k]]
          if (!length(items_k)) next
          Ck <- params$C_list[[k]]
          Y_rep[, items_k] <- .simulate_testlet_block(theta_new, a, b, items_k, Ck)
        }
      } else {
        MU_all <- outer(theta_new, a, "*") - matrix(rep(a*b, each = N_new), nrow = N_new)
        P_all <- pnorm(MU_all)
        Y_rep <- matrix(rbinom(length(P_all), 1L, P_all), nrow = N_new)
      }
      
      if (by_testlet) {
        d_rep[t, ] <- diag_residQ3_vec(params, theta_new, Y_rep, idx_testlets)
      } else {
        d_rep[t]   <- diag_residQ3(params, theta_new, Y_rep, idx_testlets)
      }
    }
  }
  
  if (diag == "logscore") {
    return(list(p_HPC = mean(d_rep >= d_new), d_rep = d_rep, d_new = d_new))
  }
  
  # residQ3:
  if (!by_testlet) {
    return(list(p_HPC = mean(d_rep >= d_new), d_rep = d_rep, d_new = d_new))
  } else {
    # p por testlet (ignora NAs de testlets com <2 itens)
    p_by_testlet <- mapply(function(col_rep, val_new) {
      if (is.na(val_new)) return(NA_real_)
      mean(col_rep >= val_new, na.rm = TRUE)
    }, as.data.frame(d_rep), d_new)
    
    # p global opcional (máximo entre testlets)
    p_global <- mean(apply(d_rep, 1, function(v) max(v, na.rm = TRUE)) >= max(d_new, na.rm = TRUE))
    
    names(p_by_testlet) <- paste0("T", seq_along(p_by_testlet))
    return(list(
      p_HPC_global     = p_global,
      p_HPC_by_testlet = p_by_testlet,
      d_rep_mat        = d_rep,  # R_eval x K
      d_new_vec        = d_new   # K
    ))
  }
}



## ============================================================================
## Q3 utilities: observado e preditivo (replicado) por par de itens
## Requer:
##  - Um gerador de Y_rep: via builder_fun(draw_id, ext) -> list(a,b,theta[, rho...])
##    e uma função de simulação P(Y|parâmetros). Aqui definimos um simulador binário
##    probit com dependência local opcional via C_list (testlets).
##  - pairs_list: lista de matrizes 2 x Pk com índices (i,j) por testlet.
##  - mYcNA: matriz 0/1 de respostas (N x vI)
## ----------------------------------------------------------------------------

# (Opcional) aceleração
if (!requireNamespace("matrixStats", quietly = TRUE)) {
  message("[AuxFunctions] Pacote 'matrixStats' não encontrado. Recomenda-se instalar para melhor desempenho.")
}

## ---- A.1) Resíduos de Pearson por par e Q3 observado -----------------------
## Q3 por par = cor(resid_i, resid_j). Aqui resid = (Y_i - p_i)/sqrt(p_i*(1-p_i))
## onde p_i é a proporção marginal observada (modelo nulo simples).
compute_Q3_by_pairs <- function(Y, pairs_list) {
  stopifnot(is.matrix(Y) || is.data.frame(Y))
  Y <- as.matrix(Y)
  N <- nrow(Y); vI <- ncol(Y)
  p_hat <- colMeans(Y)                        # proporções marginais (nulo simples)
  sd_hat <- sqrt(p_hat * (1 - p_hat) + 1e-9)  # evita divisão por zero
  R <- sweep(Y, 2, p_hat, "-")
  R <- sweep(R, 2, sd_hat, "/")               # resíduos "tipo Pearson" simples
  
  # Concatenar Q3 par a par dentro de cada testlet
  out <- numeric(0L)
  for (k in seq_along(pairs_list)) {
    M <- pairs_list[[k]]
    if (length(M) == 0L) next
    Pk <- ncol(M)
    for (p in 1:Pk) {
      i <- M[1, p]; j <- M[2, p]
      # correlação de Pearson entre resíduos (Q3)
      c_ij <- suppressWarnings(cor(R[, i], R[, j], use = "pairwise.complete.obs"))
      if (!is.finite(c_ij)) c_ij <- NA_real_
      out <- c(out, c_ij)
    }
  }
  out
}

## ---- A.2) Simulador Y_rep (probit) -----------------------------------------
## Gera Y_rep dado os parâmetros de um draw. Suporta:
##  - modelo 2PP independente: C_list = NULL  (sem dependência local)
##  - modelo Testlet: C_list = lista de blocos (matrizes de correlação por testlet)
## Argumentos:
##  - par: lista com componentes numéricas: a (vI), b (vI), theta (N),
##         e opcionalmente C_list = list(C_1, ..., C_K), cada C_k de dimensão nk x nk
## Retorna matriz N x vI de 0/1.
simulate_Yrep_probit <- function(par) {
  stopifnot(all(c("a", "b", "theta") %in% names(par)))
  a <- as.numeric(par$a); b <- as.numeric(par$b); theta <- as.numeric(par$theta)
  N <- length(theta); vI <- length(a)
  eta <- outer(theta, a) - matrix(rep(b, each = N), nrow = N, ncol = vI)  # N x vI
  # probit: Y = 1{Z > 0}, Z = eta + epsilon, epsilon ~ N(0,1) (ou estrutura com DL por testlet)
  
  Y_rep <- matrix(0L, nrow = N, ncol = vI)
  
  if (is.null(par$C_list) || length(par$C_list) == 0L) {
    # Independente entre itens (2PP)
    Z <- eta + matrix(rnorm(N * vI), N, vI)
    Y_rep[] <- (Z > 0)
  } else {
    # Dependência local por testlet: gera blocos correlacionados para itens do testlet
    # Precisamos também do mapeamento "itens por testlet" em par$idx_testlets
    stopifnot("idx_testlets" %in% names(par))
    for (n in 1:N) {
      # Começa como ruído independente
      eps <- rnorm(vI)
      # Ajusta blocos com correlação C_k
      for (k in seq_along(par$idx_testlets)) {
        items <- par$idx_testlets[[k]]
        nk <- length(items)
        if (nk < 2L) next
        Ck <- par$C_list[[k]]
        # Gera ruído correlacionado para o bloco e substitui no vetor eps
        L <- tryCatch(chol(Ck), error = function(e) NULL)
        if (is.null(L)) {
          # fallback numérico: regulariza
          eig <- eigen(Ck, symmetric = TRUE, only.values = FALSE)
          eig$values[eig$values < 1e-6] <- 1e-6
          Ck_reg <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
          L <- chol(Ck_reg)
        }
        z_block <- L %*% rnorm(nk)
        eps[items] <- as.numeric(z_block)
      }
      Z_n <- eta[n, ] + eps
      Y_rep[n, ] <- as.integer(Z_n > 0)
    }
  }
  Y_rep
}

## ---- A.3) Construtor de parâmetros a partir de arrays extraídos -------------
## Estes wrappers assumem suas funções existentes build_params_*_from_arrays(),
## apenas adicionando C_list e idx_testlets quando necessário para simulação.
## Ajuste conforme sua implementação atual dessas funções de 'builder'.

# NOTA: Estas funções assumem que 'builder_fun(draw_id, ext)' retorna, no mínimo:
# list(a=..., b=..., theta=..., [e para Testlet] C_list=..., idx_testlets=...)
# Se seu builder não monta C_list/idx_testlets, inclua-os lá ou aqui.

## ---- A.4) Réplicas Q3 por par para um conjunto de draws ---------------------
simulate_Q3_rep_by_pairs <- function(ext, builder_fun, pairs_list, Y_obs, R_eval = 250, seed = 1) {
  set.seed(seed)
  P_total <- sum(sapply(pairs_list, function(M) if (length(M)) ncol(M) else 0L))
  if (P_total == 0L) return(matrix(numeric(0), nrow = 0, ncol = 0))
  N <- nrow(Y_obs); vI <- ncol(Y_obs)
  
  # Prealocação: S x P (S ~ R_eval)
  Q3_rep <- matrix(NA_real_, nrow = R_eval, ncol = P_total)
  
  for (s in 1:R_eval) {
    # constrói parâmetros do draw s
    par_s <- builder_fun(s, ext)
    # para o simulador saber os blocos de testlet (se houver)
    if (!is.null(par_s$idx_testlets)) {
      # ok
    }
    # simula Y_rep e calcula Q3 por par
    Y_rep <- simulate_Yrep_probit(par_s)
    Q3_s  <- compute_Q3_by_pairs(Y_rep, pairs_list)
    Q3_rep[s, ] <- Q3_s
  }
  Q3_rep
}

## ---- A.5) HPC completo: retorna componentes e p-global ----------------------
## Em vez de retornar só p_HPC, devolve também Q3_obs e Q3_rep (S x P)
compute_hpc_q3_components_arrays <- function(
    ext, builder_fun, Y_new, idx_testlets, pairs_list,
    R_eval = 250, seed = 1
) {
  set.seed(seed)
  
  # 1) Q3 observado por par
  Q3_obs <- compute_Q3_by_pairs(Y_new, pairs_list)
  stopifnot(is.numeric(Q3_obs))
  
  # 2) Réplicas preditivas de Q3 por par
  Q3_rep <- simulate_Q3_rep_by_pairs(
    ext         = ext,
    builder_fun = builder_fun,
    pairs_list  = pairs_list,
    Y_obs       = Y_new,
    R_eval      = R_eval,
    seed        = seed
  )
  S <- nrow(Q3_rep)
  if (S < 5L) warning("Poucas réplicas em Q3_rep; aumente R_eval para HPC estável.")
  
  # 3) p-valor global (alvo = média dos Q3 por par)
  Q3_mean_obs <- mean(Q3_obs, na.rm = TRUE)
  Q3_mean_rep <- rowMeans(Q3_rep, na.rm = TRUE)
  mu <- mean(Q3_mean_rep)
  p_global <- mean(abs(Q3_mean_rep - mu) >= abs(Q3_mean_obs - mu))
  
  list(
    Q3_obs   = Q3_obs,      # vetor de tamanho P
    Q3_rep   = Q3_rep,      # matriz S x P
    p_HPC    = p_global,    # p-valor global (média Q3)
    target   = "mean(Q3)",  # estatística-alvo usada no p_global
    R_eval   = R_eval,
    S_draws  = S
  )
}


# -------- util: gera bloco de testlet vetorizado ----------
.simulate_testlet_block <- function(theta_vec, a, b, items_k, Ck) {
  n  <- length(theta_vec); Ik <- length(items_k)
  mu <- outer(theta_vec, a[items_k], "*") -
    matrix(rep(a[items_k]*b[items_k], each = n), nrow = n)
  L  <- chol(Ck)                 # superior
  G  <- matrix(rnorm(n * Ik), n, Ik)
  Z  <- mu + G %*% t(L)          # cov = Ck
  (Z > 0L) * 1L
}

# -------- Q3 pareados dentro de um testlet ----------
# Retorna o vetor de correlações residuais (off-diagonal) para UM testlet
.q3_pairs_one_testlet <- function(Yk, Pk, eps = 1e-9) {
  # resíduos centrados por item
  Pk <- pmin(pmax(Pk, eps), 1 - eps)
  Rk <- Yk - Pk
  Rk <- scale(Rk, center = TRUE, scale = FALSE)
  
  # correlação com proteção de variâncias pequenas
  S  <- crossprod(Rk) / (nrow(Rk) - 1)
  sd <- sqrt(pmax(diag(S), eps))
  Ck <- S / (sd %o% sd)
  Ck[!is.finite(Ck)] <- 0
  Ck[upper.tri(Ck)]
}

# -------- Q3 pareados por testlet OU global ----------
# Retorna lista (por testlet) de vetores Q3, ou um único vetor se global=TRUE
q3_pairs <- function(params, theta, Y, idx_testlets, global = FALSE, eps = 1e-9) {
  a <- params$a; b <- params$b
  MU <- outer(theta, a, "*") - matrix(rep(a*b, each = nrow(Y)), nrow = nrow(Y))
  P  <- pnorm(MU)
  
  keep <- lengths(idx_testlets) >= 2
  if (!any(keep)) return(if (global) numeric(0) else vector("list", length(idx_testlets)))
  
  if (global) {
    out <- c()
    for (k in which(keep)) {
      items_k <- idx_testlets[[k]]
      Yk <- as.matrix(Y[, items_k, drop = FALSE])
      Pk <- as.matrix(P[, items_k, drop = FALSE])
      out <- c(out, .q3_pairs_one_testlet(Yk, Pk, eps))
    }
    return(out)
  } else {
    res <- vector("list", length(idx_testlets))
    for (k in seq_along(idx_testlets)) {
      items_k <- idx_testlets[[k]]
      if (length(items_k) < 2) { res[[k]] <- NA; next }
      Yk <- as.matrix(Y[, items_k, drop = FALSE])
      Pk <- as.matrix(P[, items_k, drop = FALSE])
      res[[k]] <- .q3_pairs_one_testlet(Yk, Pk, eps)
    }
    return(res)
  }
}

# -------- simulador de Y_rep dado um draw --------
simulate_Yrep <- function(params, theta_new, idx_testlets, ind_items = integer(0)) {
  N_new <- length(theta_new); I <- length(params$a)
  Y_rep <- matrix(0L, nrow = N_new, ncol = I)
  a <- params$a; b <- params$b
  
  # itens independentes
  if (length(ind_items) > 0) {
    MU_ind <- outer(theta_new, a[ind_items], "*") -
      matrix(rep(a[ind_items]*b[ind_items], each = N_new), nrow = N_new)
    P_ind <- pnorm(MU_ind)
    Y_rep[, ind_items] <- matrix(rbinom(length(P_ind), 1L, P_ind), nrow = N_new)
  }
  
  # testlets (se houver dependência)
  if (!is.null(params$C_list)) {
    for (k in seq_along(idx_testlets)) {
      items_k <- idx_testlets[[k]]
      if (!length(items_k)) next
      Ck <- params$C_list[[k]]
      Y_rep[, items_k] <- .simulate_testlet_block(theta_new, a, b, items_k, Ck)
    }
  } else {
    # tudo independente (2PP)
    MU_all <- outer(theta_new, a, "*") - matrix(rep(a*b, each = N_new), nrow = N_new)
    P_all  <- pnorm(MU_all)
    Y_rep  <- matrix(rbinom(length(P_all), 1L, P_all), nrow = N_new)
  }
  Y_rep
}

# -------- prepara dados de envelope (um modelo) ----------
# ext: draws extraídos; build_params_fun: (dr, ext, ...) -> list(a,b, C_list?)
# Retorna tibble com (model, facet, x, F_low, F_high, F_obs) p/ plot
q3_envelope_data <- function(ext, build_params_fun, Y_new, idx_testlets,
                             model_label = "Testlet",
                             ind_items = integer(0),
                             R_eval = 200, seed = 123,
                             global = FALSE, grid_len = 200, ...) {
  set.seed(seed)
  
  R_avail <- nrow(ext$a)
  stopifnot(R_avail >= 1L)
  draw_ids <- sample.int(R_avail, size = min(R_eval, R_avail), replace = FALSE)
  
  N_new <- nrow(Y_new)
  theta_new <- rnorm(n = N_new)   # 1 amostra de theta_new fixa por envelope
  
  # --- Q3 observado ---
  # Aqui usamos o mesmo theta_new só para construir P; ECDF só usa Y_new.
  params_1st <- build_params_fun(draw_ids[1], ext, ...)
  q_obs <- q3_pairs(params_1st, theta_new, Y_new, idx_testlets, global = global)
  
  # Normaliza formato: lista por "facet"
  to_list_facets <- function(x) {
    if (global) list(Global = x) else {
      out <- x
      names(out) <- paste0("T", seq_along(x))
      out
    }
  }
  q_obs_list <- to_list_facets(q_obs)
  
  # --- Q3 replicados para cada draw ---
  q_rep_list <- vector("list", length(draw_ids))
  for (t in seq_along(draw_ids)) {
    dr <- draw_ids[t]
    params <- build_params_fun(dr, ext, ...)
    Y_rep  <- simulate_Yrep(params, theta_new, idx_testlets, ind_items)
    q_rep_list[[t]] <- to_list_facets(q3_pairs(params, theta_new, Y_rep, idx_testlets, global))
  }
  
  facets <- names(q_obs_list)
  # grid em x baseado no suporte observado + replicados
  all_vals <- c(unlist(q_obs_list), unlist(q_rep_list, recursive = TRUE, use.names = FALSE))
  x_grid   <- sort(unique(quantile(all_vals, probs = seq(0, 1, length.out = grid_len))))
  
  # função ECDF
  Fn <- function(v, x) mean(v <= x)
  
  # monta tibble com envelopes
  out <- map_dfr(facets, function(fa) {
    # pega lista de vetores replicados para esse facet
    reps_fa <- map(q_rep_list, ~ .x[[fa]])
    # alguns testlets podem ter NA (menos de 2 itens) -> ignorar
    reps_fa <- reps_fa[!map_lgl(reps_fa, ~ all(is.na(.x)) || length(.x) == 0)]
    if (length(reps_fa) == 0) return(tibble())
    
    F_rep_mat <- sapply(x_grid, function(x) vapply(reps_fa, Fn, numeric(1), x = x))
    # F_rep_mat: R' (linhas) x |x_grid| (colunas); transpomos para quantis por x
    F_low  <- apply(F_rep_mat, 2, quantile, probs = 0.025, na.rm = TRUE)
    F_high <- apply(F_rep_mat, 2, quantile, probs = 0.975, na.rm = TRUE)
    
    F_obs <- sapply(x_grid, function(x) Fn(q_obs_list[[fa]], x))
    
    tibble(
      model  = model_label,
      facet  = fa,
      x      = x_grid,
      F_low  = F_low,
      F_high = F_high,
      F_obs  = F_obs
    )
  })
  
  out
}

# -------- plotagem: dois painéis lado a lado por modelo, facet por testlet ----------
plot_q3_ecdf_envelope <- function(df_env, ncol_models = 2) {
  # df_env: bind_rows de saídas de q3_envelope_data para "Testlet" e "2PP"
  df_long <- df_env %>%
    mutate(model = factor(model, levels = unique(df_env$model))) %>%
    pivot_longer(c(F_low, F_high, F_obs),
                 names_to = "what", values_to = "Fval")
  
  # dados para ribbon e curva
  df_ribbon <- df_env
  df_curve  <- df_env %>% select(model, facet, x, F_obs)
  
  ggplot() +
    geom_ribbon(data = df_ribbon,
                aes(x = x, ymin = F_low, ymax = F_high, fill = model),
                alpha = 0.20, inherit.aes = FALSE) +
    geom_line(data = df_curve,
              aes(x = x, y = F_obs, color = model),
              linewidth = 0.9, show.legend = TRUE) +
    facet_grid(~ model + facet, scales = "free_x", switch = "x") +
    labs(x = expression(Q[3]), y = "ECDF",
         title = "ECDF do Q3 com envelope preditivo (95%)",
         subtitle = "Painéis por testlet (ou Global) e comparação entre modelos") +
    theme_minimal(base_size = 13) +
    theme(
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
}

# df_env: bind_rows de q3_envelope_data(...) para "Testlet" e "2PP"
plot_q3_ecdf_envelope_gray <- function(df_env,
                                       line_size = 0.9,
                                       alpha_fill = 0.22,
                                       base_size = 13,
                                       title = '',
                                       subtitle = '') {
  df_env <- df_env %>%
    dplyr::mutate(
      model = factor(model, levels = c("Testlet 2PP", "2PP")),
      facet = factor(facet, levels = unique(facet))
    )
  
  df_ribbon <- df_env
  df_curve  <- df_env %>% dplyr::select(model, facet, x, F_obs)
  
  ggplot() +
    # envelope cinza uniforme
    geom_ribbon(
      data = df_ribbon,
      aes(x = x, ymin = F_low, ymax = F_high,
          group = interaction(model, facet)),
      fill = "gray60", alpha = alpha_fill, inherit.aes = FALSE
    ) +
    # curvas observadas
    geom_line(
      data = df_curve,
      aes(x = x, y = F_obs, linetype = model),
      linewidth = line_size, color = "black"
    ) +
    # mantém disposição lado a lado (Testlet | 2PP),
    # mas oculta os títulos de faixa da variável 'model'
    facet_grid(~ model + facet, scales = "free_x", switch = "x",
               labeller = labeller(model = ~ "", facet = label_value)) +
    scale_linetype_manual(values = c("solid", "dashed"), name = NULL) +
    labs(
      x = expression(Q[3]), y = "ECDF",
      title = title, subtitle = subtitle
    ) +
    theme_minimal(base_size = base_size) +
    theme(
      strip.text = element_text(face = "bold"),
      strip.background.x = element_blank(),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      legend.key.width = unit(1.5, "cm")
    )
}



# Salva em PDF com tamanho fixo (cm) e, se disponível, incorpora fontes (cairo_pdf)
save_q3_ecdf_envelope_pdf <- function(df_env,
                                      file = "Fig_Q3_ECDF_envelope.pdf",
                                      width_cm = 18, height_cm = 10,
                                      embed_fonts = TRUE,
                                      ...) {
  p <- plot_q3_ecdf_envelope_gray(df_env, ...)
  # escolhe dispositivo
  dev_fun <- if (embed_fonts && capabilities("cairo")) cairo_pdf else "pdf"
  ggsave(filename = file, plot = p,
         width = width_cm, height = height_cm, units = "cm",
         device = dev_fun, dpi = 300)
  message(sprintf("Figura salva em '%s' (%gx%g cm).", file, width_cm, height_cm))
  invisible(p)
}





# ---- cálculo de Q3 pareado por testlet ----
.q3_matrix <- function(Y, P, eps = 1e-9) {
  P <- pmin(pmax(P, eps), 1 - eps)
  R <- Y - P
  R <- scale(R, center = TRUE, scale = FALSE)
  S <- crossprod(R) / (nrow(R) - 1)
  sd <- sqrt(pmax(diag(S), eps))
  C <- S / (sd %o% sd)
  diag(C) <- NA
  C
}

# ---- simula replicados e devolve ΔQ3 médio ----
delta_q3_testlet <- function(ext, build_params_fun, Y_new, idx_testlets,
                             ind_items = integer(0), R_eval = 100,
                             seed = 123, model_label = "Testlet", ...) {
  set.seed(seed)
  R_avail <- nrow(ext$a)
  draw_ids <- sample.int(R_avail, size = min(R_eval, R_avail), replace = FALSE)
  
  N <- nrow(Y_new)
  theta_new <- rnorm(N)
  a <- build_params_fun(draw_ids[1], ext, ...)$a
  b <- build_params_fun(draw_ids[1], ext, ...)$b
  
  # Q3 observados
  MU <- outer(theta_new, a, "*") - matrix(rep(a*b, each = N), nrow = N)
  P  <- pnorm(MU)
  Q3_obs_list <- vector("list", length(idx_testlets))
  for (k in seq_along(idx_testlets)) {
    items_k <- idx_testlets[[k]]
    if (length(items_k) < 2) next
    Q3_obs_list[[k]] <- .q3_matrix(Y_new[, items_k, drop = FALSE], P[, items_k, drop = FALSE])
  }
  
  # Q3 replicados
  Q3_rep_accum <- lapply(Q3_obs_list, function(M) if (!is.null(M)) matrix(0, nrow(M), ncol(M)) else NULL)
  
  for (dr in draw_ids) {
    params <- build_params_fun(dr, ext, ...)
    Y_rep <- simulate_Yrep(params, theta_new, idx_testlets, ind_items)
    MU_r  <- outer(theta_new, params$a, "*") -
      matrix(rep(params$a * params$b, each = N), nrow = N)
    P_r <- pnorm(MU_r)
    
    for (k in seq_along(idx_testlets)) {
      items_k <- idx_testlets[[k]]
      if (length(items_k) < 2) next
      Q3_rep_accum[[k]] <- Q3_rep_accum[[k]] +
        .q3_matrix(Y_rep[, items_k, drop = FALSE], P_r[, items_k, drop = FALSE])
    }
  }
  
  # ΔQ3 = obs - média replicada
  Delta_list <- vector("list", length(idx_testlets))
  for (k in seq_along(idx_testlets)) {
    if (is.null(Q3_obs_list[[k]])) next
    Delta_list[[k]] <- Q3_obs_list[[k]] - (Q3_rep_accum[[k]] / R_eval)
  }
  
  # empacota em tibble longo (índices numéricos, sem fatores)
  df_list <- lapply(seq_along(Delta_list), function(k) {
    M <- Delta_list[[k]]
    if (is.null(M)) return(NULL)
    idx <- which(!is.na(M), arr.ind = TRUE)  # pega só entradas válidas
    tibble::tibble(
      Testlet = paste0("T", k),
      Item_i  = as.integer(idx[, 1]),
      Item_j  = as.integer(idx[, 2]),
      DeltaQ3 = as.numeric(M[idx])
    )
  })
  dplyr::bind_rows(df_list) %>%
    dplyr::mutate(model = model_label)
  
}

# ---- combina modelos (ou diferença entre eles) e plota heatmap ----
plot_delta_q3_heatmap <- function(df_testlet, df_2pp = NULL,
                                  diff_only = FALSE,
                                  gray = TRUE,
                                  limits = c(-.5, .5),
                                  base_size = 12,
                                  scales = "free") {
  stopifnot(is.data.frame(df_testlet),
            all(c("Testlet","Item_i","Item_j","DeltaQ3") %in% names(df_testlet)))
  if (!is.null(df_2pp)) stopifnot(all(c("Testlet","Item_i","Item_j","DeltaQ3") %in% names(df_2pp)))
  
  # força tipos
  fix_num <- function(d) {
    d %>%
      dplyr::mutate(
        Item_i = as.integer(Item_i),
        Item_j = as.integer(Item_j),
        Testlet = as.character(Testlet),
        model = if ("model" %in% names(d)) as.character(model) else "Model"
      )
  }
  df_testlet <- fix_num(df_testlet)
  if (!is.null(df_2pp)) df_2pp <- fix_num(df_2pp)
  
  # combina/diferença
  if (!is.null(df_2pp) && diff_only) {
    df_plot <- dplyr::inner_join(
      df_testlet, df_2pp,
      by = c("Testlet","Item_i","Item_j"),
      suffix = c("_test","_2pp")
    ) %>%
      dplyr::mutate(
        DeltaQ3 = DeltaQ3_test - DeltaQ3_2pp,
        model = "Δ(Testlet−2PP)"
      ) %>%
      dplyr::select(Testlet, Item_i, Item_j, DeltaQ3, model)
  } else {
    df_plot <- dplyr::bind_rows(df_testlet, df_2pp)
  }
  
  # triângulo superior e finitos
  df_plot <- df_plot %>%
    dplyr::filter(is.finite(DeltaQ3), Item_i < Item_j)
  
  if (nrow(df_plot) == 0) {
    stop("Sem pares válidos para plotar (verifique se há testlets com >= 2 itens).")
  }
  
  # fatores ordenados
  df_plot <- df_plot %>%
    dplyr::mutate(
      Item_i  = factor(Item_i, levels = sort(unique(Item_i))),
      Item_j  = factor(Item_j, levels = sort(unique(Item_j))),
      Testlet = factor(Testlet, levels = sort(unique(Testlet))),
      model   = factor(model, levels = sort(unique(model)))
    )
  
  p <- ggplot(df_plot, aes(Item_i, Item_j, fill = DeltaQ3)) +
    geom_tile(color = "white", linewidth = 0.2) +
    facet_grid(~ model + Testlet, scales = scales, switch = "x") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(limits = rev, expand = c(0, 0)) +
    scale_fill_gradient2(
      low = if (gray) "gray40" else "blue",
      mid = "white",
      high = if (gray) "gray10" else "red",
      midpoint = 0, limits = limits,
      name = expression(Delta*Q[3])
    ) +
    theme_minimal(base_size = base_size) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid = element_blank(),
      aspect.ratio = 1   # ✅ mantém “quadrado” mesmo com scales livres
    ) +
    labs(x = "Item i", y = "Item j",
         title = '',
         subtitle = if (diff_only) "" else NULL)
  
  p
}



# ---- salvar em PDF ----
save_delta_q3_heatmap_pdf <- function(df_testlet, df_2pp = NULL,
                                      diff_only = FALSE,
                                      file = "Fig_DeltaQ3_heatmap.pdf",
                                      width_cm = 18, height_cm = 10,
                                      embed_fonts = TRUE) {
  p <- plot_delta_q3_heatmap(df_testlet, df_2pp, diff_only = diff_only)
  dev_fun <- if (embed_fonts && capabilities("cairo")) cairo_pdf else "pdf"
  ggsave(filename = file, plot = p,
         width = width_cm, height = height_cm, units = "cm",
         device = dev_fun, dpi = 300)
  message(sprintf("Mapa de calor salvo em '%s'.", file))
  invisible(p)
}



# ---- média do Q3 por testlet (e global) para um Y dado e um conjunto P ----
.mean_upper <- function(M) mean(M[upper.tri(M)], na.rm = TRUE)

q3bar_per_testlet <- function(Y, P, idx_testlets, eps = 1e-9) {
  K <- length(idx_testlets)
  out <- rep(NA_real_, K)
  for (k in seq_len(K)) {
    items_k <- idx_testlets[[k]]
    if (length(items_k) < 2) next
    Qk <- .q3_matrix(Y[, items_k, drop = FALSE], P[, items_k, drop = FALSE], eps)
    out[k] <- .mean_upper(Qk)
  }
  out
}

q3bar_global <- function(Y, P, idx_testlets, eps = 1e-9) {
  vals <- c()
  for (k in seq_along(idx_testlets)) {
    items_k <- idx_testlets[[k]]
    if (length(items_k) < 2) next
    Qk <- .q3_matrix(Y[, items_k, drop = FALSE], P[, items_k, drop = FALSE], eps)
    vals <- c(vals, Qk[upper.tri(Qk)])
  }
  if (length(vals)) mean(vals, na.rm = TRUE) else NA_real_
}

# ---- gera R réplicas de \bar Q3 por testlet e global para UM modelo ----
# ext/build_params_fun como nas suas rotinas; usa o mesmo theta_new para todas as réplicas
q3bar_replicates <- function(ext, build_params_fun, Y_new, idx_testlets,
                             ind_items = integer(0), R_eval = 200, seed = 123, ...) {
  set.seed(seed)
  R_avail <- nrow(ext$a); stopifnot(R_avail >= 1L)
  draw_ids <- sample.int(R_avail, size = min(R_eval, R_avail), replace = FALSE)
  
  N <- nrow(Y_new)
  theta_new <- rnorm(N)
  
  # usamos P baseado em cada draw do próprio modelo (para os resíduos PPC)
  K <- length(idx_testlets)
  Qbar_mat <- matrix(NA_real_, nrow = length(draw_ids), ncol = K)  # por testlet
  Qbar_glb <- numeric(length(draw_ids))                             # global
  
  for (t in seq_along(draw_ids)) {
    dr <- draw_ids[t]
    params <- build_params_fun(dr, ext, ...)
    # réplicas
    Yrep <- simulate_Yrep(params, theta_new, idx_testlets, ind_items)
    MU   <- outer(theta_new, params$a, "*") -
      matrix(rep(params$a * params$b, each = N), nrow = N)
    P    <- pnorm(MU)
    
    Qbar_mat[t, ] <- q3bar_per_testlet(Yrep, P, idx_testlets)
    Qbar_glb[t]   <- q3bar_global(Yrep, P, idx_testlets)
  }
  
  list(Qbar_by_testlet = Qbar_mat,
       Qbar_global     = Qbar_glb,
       draw_ids        = draw_ids)
}

# ---- Probabilidade de superioridade: Pr( \bar Q3_Testlet < \bar Q3_2PP ) ----
# pareamos por índice t (mesmo comprimento = min(R_T, R_2PP))
model_superiority_probability <- function(rep_testlet, rep_2pp, idx_testlets) {
  # alinhar número de réplicas
  m <- min(nrow(rep_testlet$Qbar_by_testlet), nrow(rep_2pp$Qbar_by_testlet))
  Q_T  <- rep_testlet$Qbar_by_testlet[seq_len(m), , drop = FALSE]
  Q_2P <- rep_2pp$Qbar_by_testlet[seq_len(m), , drop = FALSE]
  G_T  <- rep_testlet$Qbar_global[seq_len(m)]
  G_2P <- rep_2pp$Qbar_global[seq_len(m)]
  
  # por testlet
  p_by <- colMeans(Q_T < Q_2P, na.rm = TRUE)
  n_by <- colSums(is.finite(Q_T) & is.finite(Q_2P))
  se_by <- sqrt(p_by * (1 - p_by) / pmax(n_by, 1))
  lo_by <- pmax(0, p_by - 1.96 * se_by)
  hi_by <- pmin(1, p_by + 1.96 * se_by)
  
  # global
  p_g  <- mean(G_T < G_2P, na.rm = TRUE)
  n_g  <- sum(is.finite(G_T) & is.finite(G_2P))
  se_g <- sqrt(p_g * (1 - p_g) / pmax(n_g, 1))
  lo_g <- max(0, p_g - 1.96 * se_g)
  hi_g <- min(1, p_g + 1.96 * se_g)
  
  tibble(
    Scope = c("Global", paste0("T", seq_along(idx_testlets))),
    n     = c(n_g, n_by),
    p_sup = c(p_g, p_by),
    se    = c(se_g, se_by),
    lo95  = c(lo_g, lo_by),
    hi95  = c(hi_g, hi_by)
  )
}

# ---- Plot (tons de cinza) e salvar em PDF ----
plot_superiority_gray <- function(df, base_size = 12) {
  df2 <- df %>% filter(is.finite(p_sup), n > 0)
  ggplot(df2, aes(x = Scope, y = p_sup)) +
    geom_hline(yintercept = 0.5, linetype = "dotted") +
    geom_point(color = "black", size = 2) +
    geom_errorbar(aes(ymin = lo95, ymax = hi95), width = 0.15, color = "black") +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = NULL, y = "Probabilidade de superioridade\nPr(  Q̄3(Testlet) < Q̄3(2PP) )",
         title = "Probabilidade de que o modelo Testlet reduza LID (vs 2PP)") +
    theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold")
    )
}

save_superiority_pdf <- function(df, file = "Prob_Superioridade_Testlet.pdf",
                                 width_cm = 16, height_cm = 9,
                                 embed_fonts = TRUE) {
  p <- plot_superiority_gray(df)
  dev_fun <- if (embed_fonts && capabilities("cairo")) cairo_pdf else "pdf"
  ggsave(file, p, width = width_cm, height = height_cm, units = "cm", device = dev_fun, dpi = 300)
  message(sprintf("Figura salva em '%s'.", file))
  invisible(p)
}








