functions{
  // Sum of all entries in a 2D integer array.
  // Used in transformed data to pre-allocate the positions of Y==1 (positives).
  int sum2d(int[,] a){
    int s = 0;
    for (i in 1:size(a))
      s += sum(a[i]);
    return s;
  }

  // Build an L×L correlation matrix from a vector of correlations "rhovec":
  //  - if size(rhovec) == 1     : HU (compound symmetry / exchangeable)
  //      all off-diagonal correlations are equal to rhovec[1]
  //  - if size(rhovec) == L - 1 : HT (Toeplitz-by-lag)
  //      correlation depends on lag d = |i-j|, using rhovec[d], d=1..L-1
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
            // HU: single correlation for every pair
            R[i,j] = rhovec[1];
          } else {
            // HT: lag-dependent correlation
            R[i,j] = rhovec[d];
          }
        }
      }
    }
    return R;
  }
}

data{
  // Dimensions
  int<lower=1> I;                 // total number of items
  int<lower=1> K;                 // number of testlets
  int<lower=1> N;                 // number of respondents

  // Testlet layout (1-based indexing)
  int<lower=1> dk[K];             // starting index of each testlet block in the item vector
  int<lower=1> nk[K];             // number of items within each testlet

  // Independent (non-testlet) items
  int<lower=0> n_ind;             // number of independent items = I - sum(nk)
  int<lower=1,upper=I> ind_items[n_ind]; // item indices for independent items

  // Binary responses (0/1). Missing must be handled before Stan (e.g., imputed or set to 0).
  int<lower=0,upper=1> Y[N,I];

  // "Ragged" indexing into a single rho vector:
  // each testlet k uses rho_len[k] correlations starting at rho_start[k].
  int<lower=0> rho_len[K];        // number of rho's used by testlet k:
                                 //   1 if HU, or nk[k]-1 if HT
  int<lower=1> rho_start[K];      // starting position (1-based) in rho_global


  // Prior scales (sigma_a is treated as a parameter below; sigma_b and sigma_rho are data)
  real<lower=0> sigma_b;
  real<lower=0> sigma_rho;
}

transformed data {
  // Latent-variable augmentation for probit:
  // For each response Y[n,i], introduce a latent Z[n,i] s.t.
  //   Y=1 <=> Z>0 and Y=0 <=> Z<=0
  //
  // We store indices of positive and negative responses to declare truncated vectors
  // Z_pos (lower=0) and Z_neg (upper=0), which speeds up sampling.
  
  // Adapted from: https://mc-stan.org/docs/2_18/stan-users-guide/multivariate-outcomes.html

  int<lower=0> N_pos;
  int<lower=1,upper=N> n_pos[sum2d(Y)];  // respondent indices for Y==1
  int<lower=1,upper=I> d_pos[size(n_pos)]; // item indices for Y==1

  int<lower=0> N_neg;
  int<lower=1,upper=N> n_neg[(N*I) - size(n_pos)]; // respondent indices for Y==0
  int<lower=1,upper=I> d_neg[size(n_neg)];         // item indices for Y==0

  N_pos = size(n_pos);
  N_neg = size(n_neg);

  // Populate (n_pos, d_pos) and (n_neg, d_neg) with the coordinates of Y==1 and Y==0
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
  // Person parameters
  vector[N] theta;             // latent traits

  // Item parameters
  vector<lower=0>[I] a;        // discriminations (positive)
  vector[I] b;                 // difficulties

  // Data augmentation: latent Z values (truncated to enforce Y)
  vector<lower=0>[N_pos] Z_pos;
  vector<upper=0>[N_neg] Z_neg;

  // Global rho vector, concatenating all testlet correlation parameters
  // (each testlet uses segment(rho_global, rho_start[k], rho_len[k]))
  vector<lower=-1,upper=1>[sum(rho_len)] rho_global;

  // Prior scale for discriminations (treated hierarchically)
  real<lower=0> sigma_a;

  // Alternative unconstrained parameterization (commented out):
  // vector[sum(rho_len)] z_rho;  // unconstrained, then map to (-1,1) via logistic
}

transformed parameters{
  // Linear predictors and augmented latent matrix
  matrix[N, I] eta;  // eta[n,i] = a[i] * (theta[n] - b[i])
  matrix[N, I] Z;    // latent response matrix

  // Matrices defining the antedependence / conditional regression representation
  // for the testlet blocks:
  //   Phi_mat: regression coefficients linking earlier Z's to later Z's (within each testlet)
  //   D      : conditional variances (diagonal in this construction)
  matrix[I, I] Phi_mat;
  matrix[I, I] D;

  // Initialize buffers
  Phi_mat = rep_matrix(0, I, I);
  D       = rep_matrix(0, I, I);

  // Compute linear predictor for all (n,i)
  for (n in 1:N)
    for (i in 1:I)
      eta[n,i] = a[i] * (theta[n] - b[i]);

  // Fill the latent Z matrix using the truncated vectors
  for (n in 1:N_pos)
    Z[n_pos[n], d_pos[n]] = Z_pos[n];
  for (n in 1:N_neg)
    Z[n_neg[n], d_neg[n]] = Z_neg[n];

  // Build block-specific structures for each testlet and place them into
  // the appropriate submatrices of Phi_mat and D.
  {
    int off; off = 0; // not used (left for clarity / possible future refactoring)
    for (k in 1:K) {
      int L  = nk[k];   // testlet size
      int i0 = dk[k];   // starting item index for this testlet block

      // Extract rho segment for testlet k
      int rs = rho_start[k];
      int rl = rho_len[k];

      // Build testlet correlation matrix Sigma_k from rho parameters
      vector[rl] rh = segment(rho_global, rs, rl); // length = 1 (HU) or L-1 (HT)
      matrix[L,L] Sk = corr_from_rho(L, rh);

      // Convert correlation matrix to the antedependence form used in the likelihood.
      // Steps:
      //  1) Cholesky: Sk = Lk * Lk'
      //  2) Normalize Lk to have unit diagonal, move scaling into Dk
      //  3) Tk = inv(Lk) and Phik = I - Tk  (regression coefficients)
      {
        matrix[L,L] Lk  = cholesky_decompose(Sk);
        matrix[L,L] Dk  = diag_matrix(diagonal(Lk));
        Lk = Lk * inverse(Dk);   // now diag(Lk) = 1
        Dk = Dk * Dk;            // conditional variance scaling
        matrix[L,L] Tk  = inverse(Lk);
        matrix[L,L] Phik = diag_matrix(rep_vector(1, L)) - Tk;

        // Paste block matrices into the global I×I containers
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
  // ------------------------------- Priors ------------------------------------
  theta ~ normal(0, 1);

  // Discrimination prior: gamma with mean ~ 1 and variance controlled by sigma_a
  // Using shape=1/sigma_a^2, rate=1/sigma_a^2 implies E[a]=1 and Var[a]=sigma_a^2.
  a ~ gamma(1 / square(sigma_a), 1 / square(sigma_a));
  sigma_a ~ cauchy(0, 5);

  // Difficulty prior
  b ~ normal(0, sigma_b);

  // Rho priors: independent normals within each testlet segment
  // (segment-wise, allowing different lengths per testlet)
  for (k in 1:K){
    segment(rho_global, rho_start[k], rho_len[k]) ~ normal(0, sigma_rho);
  }

  // ----------------------------- Likelihood ----------------------------------
  // Latent-variable likelihood using the conditional (antedependence) factorization.
  // For each respondent n:
  //  (1) For each testlet block, evaluate the sequential conditional normals defined
  //      by Phi_mat and D (within-testlet dependence).
  //  (2) For independent items, use N(eta, 1) for the latent Z.

  for (n in 1:N){
    // (1) Testlet blocks
    for (k in 1:K){
      int L  = nk[k];
      int i0 = dk[k];

      // First item in block: marginal N(eta, sqrt(Dii))
      target += normal_lpdf(Z[n, i0] | eta[n, i0], sqrt(D[i0, i0]));

      // Remaining items: conditional regressions on previous Z's within the testlet
      if (L >= 2){
        for (ii in 2:L){
          int i = i0 + ii - 1;

          // Conditional mean starts at eta[n,i] then adds regression terms
          real pred = eta[n, i];

          // Regress on previous positions in this testlet block
          // (kk indexes earlier items in the block)
          for (kk in 0:(ii-2)) {
            int j = i0 + kk;
            pred += Phi_mat[i, j] * (Z[n, j] - eta[n, j]);
          }

          target += normal_lpdf(Z[n, i] | pred, sqrt(D[i, i]));
        }
      }
    }

    // (2) Independent items: latent Z ~ N(eta, 1)
    for (h in 1:n_ind) {
      int i = ind_items[h];
      target += normal_lpdf(Z[n, i] | eta[n, i], 1);
    }
  }
}


