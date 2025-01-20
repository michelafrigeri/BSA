
data {
  int T;                     // num union observations (known)
  int M;                     // num basis (fixed)
  int C;                     // num of pollutants (known) 
  int N;                     // num of sites (known)
  int p;                     // num covariates (known)
  int K;                     // num sources (unknown)

  array[N] int T_i;          // num observations x staz (known)
  array[N,T] int idx_ti;     // observations' index x staz (known)

  array[N] matrix[C,T] Y;    // observed values (with NAs)
  matrix[N,p] X;             // covariates
  matrix[M,T] B;             // splines matrix
  matrix[N,N] D;             // distance matrix for sites
  real eps;                  // shrinkage threshold
}



parameters {
  vector<lower=0>[C] sigma;
  
  array[K] vector[p] beta;
  vector[p] m0;
  real<lower=0> s0, r0;

  matrix[K,N] coeff_g;
  matrix<lower=0>[K,C] H;
  vector<lower=0>[K] r;

  matrix[K,M] L;
  matrix<lower=0>[M,K] phi;
  vector<lower=0>[K] delta;
}



transformed parameters {
  matrix[K,T] f = L*B;

  array[K] matrix[N,T] g;
  vector<lower=0>[K] eta;

  real count=0;
  for(k in 1:K){
      eta[k] = prod(delta[1:k]);
      if(max(abs(L[k, ])) > eps) count += 1;
      for (i in 1:N) { g[k][i] = exp(coeff_g[k,i]) * f[k]; }
  }

  array[N] matrix[C,T] mu_y;
  for (i in 1:N){
    for (c in 1:C){
      mu_y[i][c] = rep_row_vector(0,T);
      for (k in 1:K){ mu_y[i][c] += H[k,c]*g[k][i]; }
    }}
}



model {
  
  //error noise
  sigma ~ cauchy(0, 1);

  //spatial range
  r0 ~ gamma(80,0.1);
  r ~ inv_gamma(3,r0);
  
  //regression coefficients hyperparameters
  m0 ~ normal(0,1);
  s0 ~ inv_gamma(3,2);

  //MGPS hyperparameters
  to_vector(phi) ~ gamma(1.5, 1.5);
  delta[1] ~ gamma(10, 1);
  delta[2:K] ~ gamma(20, 1);

  for(k in 1:K){
    beta[k] ~ normal(m0, s0);
    L[k, ] ~ normal(0, 1/sqrt((phi[,k] * eta[k])));
    coeff_g[k] ~ multi_normal((X*beta[k])', exp(-D/r[k]));
    H[k] ~ dirichlet(rep_vector(1, C));
  }

  for (i in 1:N){
    for (c in 1:C){
      for (t in 1:T_i[i]){
        Y[i,c,idx_ti[i,t]] ~ normal(mu_y[i,c,idx_ti[i,t]], sigma[c]);
  }}}

}



generated quantities  {
  //array[N] matrix[C,T] log_lik;
  array[N] matrix[C,T] Y_pred;

  for (i in 1:N){
    Y_pred[i] = rep_matrix(0, C, T);
    for (c in 1:C){
      for (t in 1:T){
        //log_lik[i,c,idx_ti[i,t]] = normal_lpdf(Y[i,c,idx_ti[i,t]] | mu_y[i,c,idx_ti[i,t]], sigma);
        Y_pred[i,c,t] = normal_rng(mu_y[i,c,t], sigma[c]);
      }}}
}




