// SimMult2PP_Chol.stan
data {
  int<lower=0> N;
  int<lower=0> vI;
}
parameters {
  vector<lower=0>[vI] sim_a;
  vector[vI] sim_b;
  cholesky_factor_corr[vI] L_Sigma;  // fator de Cholesky da correlação
}
transformed parameters {
  matrix[vI, vI] Sigma = multiply_lower_tri_self_transpose(L_Sigma); // opcional, se quiser salvar
}
generated quantities {
  matrix[N, vI] epsilon;
  vector[N] sim_theta;
  matrix[N, vI] sim_Z;
  int<lower=0, upper=1> sim_Y[N, vI];

  for (n in 1:N) {
    epsilon[n,] = (multi_normal_cholesky_rng(rep_vector(0, vI), L_Sigma))'; // mais eficiente
    sim_theta[n] = normal_rng(0, 1);
    for (i in 1:vI) {
      sim_Z[n, i] = sim_a[i] * (sim_theta[n] - sim_b[i]) + epsilon[n, i];
      sim_Y[n, i] = (sim_Z[n, i] > 0);
    }
  }
}



