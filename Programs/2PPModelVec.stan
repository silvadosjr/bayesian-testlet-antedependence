data {
 int<lower=0> N;
 int<lower=0> vI;
 int<lower=0,upper=1> Y[N,vI];
 real<lower=0> sigma_a;
 real<lower=0> sigma_b;
 }
 parameters {
 vector[N] theta;
 vector<lower=0>[vI] a;
 vector[vI] b;
 }
 transformed parameters{
 vector<lower=0,upper=1>[vI] prob[N];
  for(n in 1:N){
   for(i in 1:vI){
   prob[n,i] = Phi_approx(a[i]*(theta[n] - b[i])); 
   }
 }
 }
model{
 theta ~ normal(0,1);
 a     ~ gamma(1 / square(sigma_a), 1 / square(sigma_a));
 b ~ normal(0,sigma_b);
 
  for(n in 1:N){
    Y[n,1:vI] ~ bernoulli(prob[n,1:vI]);
   }
}



