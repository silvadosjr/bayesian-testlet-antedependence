functions{
  int sum2d(int[,] a){
    int s = 0;
    for (i in 1:size(a))
      s += sum(a[i]);
    return s;
  }
// Constrói matriz de correlação L×L a partir de "rhovec":
  //  - size(rhovec) == 1  -> HU (compound symmetry)
  //  - size(rhovec) == L-1 -> HT (Toeplitz por lag)
  matrix corr_from_rho(int L, vector rhovec) {
    matrix[L,L] R;
    int nr = rows(rhovec);
    for (i in 1:L) {
      for (j in 1:L) {
        if (i == j) {
          R[i,j] = 1;
        } else {
          int d = abs(i - j);
          if (nr == 1) {
            // HU: uma única correlação para todo par
            R[i,j] = rhovec[1];
          } else {
            // HT: correlação depende do lag d (1..L-1)
            R[i,j] = rhovec[d];
          }
        }
      }
    }
    return R;
  }
}
data{
  int<lower=1> I;                 // nº total de itens
  int<lower=1> K;                 // nº de testlets
  int<lower=1> N;                 // nº de respondentes
  int<lower=1> dk[K];             // início de cada testlet (1-based)
  int<lower=1> nk[K];             // tamanho de cada testlet
  int<lower=0> n_ind;             // nº de itens independentes = I - sum(nk)
  int<lower=1,upper=I> ind_items[n_ind]; // índices dos itens independentes
  int<lower=0,upper=1> Y[N,I];

  // “ragged”: offsets/comprimentos dentro do vetor único de rho
  int<lower=0> rho_len[K];        // qtde de rhos que o testlet k usa (1 se HU, nk[k]-1 se HT)
  int<lower=1> rho_start[K];      // posição inicial (1-based) no vetor rho_global

  int<lower=1> S_mc;              // nº de amostras MC p/ log_lik (ex.: 200)
  
//  real<lower=0> sigma_a;
  real<lower=0> sigma_b;
  real<lower=0> sigma_rho;
}
transformed data {
  // For latent variable representation
  int<lower=0> N_pos;
  int<lower=1,upper=N> n_pos[sum2d(Y)];
  int<lower=1,upper=I> d_pos[size(n_pos)];
  int<lower=0> N_neg;
  int<lower=1,upper=N> n_neg[(N*I) - size(n_pos)];
  int<lower=1,upper=I> d_neg[size(n_neg)];

  N_pos = size(n_pos);
  N_neg = size(n_neg);
  {
    int i;
    int j;
    i = 1;
    j = 1;
    for (n in 1:N) {
      for (d in 1:I) {
        if (Y[n,d] == 1) {
          n_pos[i] = n;
          d_pos[i] = d;
          i += 1;
        } else {
          n_neg[j] = n;
          d_neg[j] = d;
          j += 1;
        }
      }
    }
  }
}
parameters{
  vector[N] theta;
  vector<lower=0>[I] a;
  vector[I] b;

  // dados aumentados:
  vector<lower=0>[N_pos] Z_pos;
  vector<upper=0>[N_neg] Z_neg;

  // todos os rhos em um único vetor
  vector<lower=-1,upper=1>[sum(rho_len)] rho_global;
  real<lower=0> sigma_a;
  
 // vector[sum(rho_len)] z_rho;   // sem limites
 
}
transformed parameters{
  matrix[N, I] eta;
  matrix[N, I] Z;
  
//  vector[sum(rho_len)] rho_global;
  
//  rho_global = 2 * inv_logit(z_rho) - 1;  // mapeia R -> (-1,1)

  // buffers gerais
  matrix[I, I] Phi_mat;
  matrix[I, I] D;

  // inicializa
  Phi_mat = rep_matrix(0, I, I);
  D       = rep_matrix(0, I, I);

  // linear predictor
  for (n in 1:N)
    for (i in 1:I)
      eta[n,i] = a[i] * (theta[n] - b[i]);

  for (n in 1:N_pos)
    Z[n_pos[n], d_pos[n]] = Z_pos[n];
  for (n in 1:N_neg)
    Z[n_neg[n], d_neg[n]] = Z_neg[n];

  // monta cada bloco de testlet
  {
    int off; off = 0; // não é usado aqui, mas deixo para clareza
    for (k in 1:K) {
      int L  = nk[k];
      int i0 = dk[k];

      // pega o pedaço de rho deste testlet
      int rs = rho_start[k];
      int rl = rho_len[k];

      // monta Sigma_k
     vector[rl] rh = segment(rho_global, rs, rl);  // len = 1 (HU) ou L-1 (HT)
     matrix[L,L] Sk = corr_from_rho(L, rh);

      // decomposição de Cholesky e conversão para T e Phi 
      {
        matrix[L,L] Lk  = cholesky_decompose(Sk);
        matrix[L,L] Dk  = diag_matrix(diagonal(Lk));
        Lk = Lk * inverse(Dk);
        Dk = Dk * Dk;
        matrix[L,L] Tk  = inverse(Lk);
        matrix[L,L] Phik = diag_matrix(rep_vector(1, L)) - Tk;

        // cola nos blocos
        for (i in 1:L) {
          for (j in 1:L) {
            Phi_mat[i0 + i - 1, i0 + j - 1] = Phik[i, j];
            D[i0 + i - 1,    i0 + j - 1]    = Dk[i, j];
          }
        }
      }
    }
  }
}
model{
  // Priors
   theta ~ normal(0, 1);
  a     ~ gamma(1 / square(sigma_a), 1 / square(sigma_a));
  sigma_a ~ cauchy(0, 5);
//   a~lognormal(log(1),sigma_a);
  b     ~ normal(0, sigma_b);
  //rho_global ~ normal(0, sigma_rho);      
  for (k in 1:K){
  segment(rho_global, rho_start[k], rho_len[k]) ~ normal(0, sigma_rho);
 }


  // Likelihood via encadeamento (igual ao seu — agora em 2 passos):
  // (1) todos os itens dos testlets, seguindo a estrutura Phi_mat/D
  for (n in 1:N){
    for (k in 1:K){
      int L  = nk[k];
      int i0 = dk[k];

      // primeiro do bloco
      target += normal_lpdf(Z[n, i0] | eta[n, i0], sqrt(D[i0, i0]));
      // demais do bloco (regressões condicionais)
      if (L >= 2){
        for (ii in 2:L){
          int i = i0 + ii - 1;
          real pred = eta[n, i];
          for (kk in 0:(ii-2)) {
            int j = i0 + kk;
            pred += Phi_mat[i, j] * (Z[n, j] - eta[n, j]);
          }
          target += normal_lpdf(Z[n, i] | pred, sqrt(D[i, i]));
        }
      }
    }
    // (2) itens independentes (variância 1)
    for (h in 1:n_ind) {
      int i = ind_items[h];
      target += normal_lpdf(Z[n, i] | eta[n, i], 1);
    }
  }
}
generated quantities {
  vector[N] log_lik;

  int max_nk;
  max_nk = nk[1];
  for (k in 2:K) if (nk[k] > max_nk) max_nk = nk[k];

  for (n in 1:N) {
    real ll = 0;

    // 1) independentes: Bernoulli-probit marginal
    for (h in 1:n_ind) {
      int i = ind_items[h];
      real eta_ni = a[i] * (theta[n] - b[i]);
      if (Y[n, i] == 1) ll += log(Phi(eta_ni));
      else              ll += log1m(Phi(eta_ni));
    }

    // 2) testlets: prob. retângulo MVN via MC
    for (k in 1:K) {
      int L  = nk[k];
      int i0 = dk[k];

      // pega rho do k
      int rs = rho_start[k];
      int rl = rho_len[k];

      // buffers de tamanho fixo
      vector[max_nk] mu = rep_vector(0, max_nk);
      matrix[max_nk, max_nk] Sigma = diag_matrix(rep_vector(1.0, max_nk));

      // média (somente 1..L)
      for (t in 1:L) {
        int i = i0 + t - 1;
        mu[t] = a[i] * (theta[n] - b[i]);
      }

      // Sigma_k
      vector[rl] rh = segment(rho_global, rs, rl);
      matrix[L,L] Sk = corr_from_rho(L, rh);

    for (i in 1:L) for (j in 1:L) Sigma[i,j] = Sk[i,j];  // copia para buffer de tamanho fixo

      // Monte Carlo
      {
        int hit = 0;
        for (s in 1:S_mc) {
          vector[max_nk] z_full = multi_normal_rng(mu, Sigma);
          int ok = 1;
          for (t in 1:L) {
            int i = i0 + t - 1;
            if (Y[n, i] == 1) { if (!(z_full[t] > 0)) { ok = 0; break; } }
            else              { if (!(z_full[t] <= 0)) { ok = 0; break; } }
          }
          hit += ok;
        }
        real p_hat = (hit + 0.5) / (S_mc + 1.0); // Laplace
        ll += log(p_hat);
      }
    }

    log_lik[n] = ll;
  }
}

// generated quantities{
//   vector[N] log_lik;
// //  vector<lower=0,upper=1>[vI] prob_pred[N];
//   for (n in 1:N) {
//       real ll = 0;
//       for (i in 1:I) {
//        real pred = a[i] * (theta[n] - b[i]);
//         if (Y[n, i] == 1) {
//           ll += log(Phi(pred));
//         } else {
//           ll += log1m(Phi(pred));  // numericamente estável
//         }
//       }
// //      prob_pred[n,i] =Phi_approx(a[i]*(theta[n] - b[i]));
//       log_lik[n] = ll;    //bernoulli_lpmf(Y[n,i]|prob_pred[n,i]);
//     }
// }


