
data {
  int T;                      // num union observations (known)
  int M;                      // num basis (fixed)
  int C;                      // num of pollutants (known) 
  int N;                      // num of sites (known)
  int p;                      // num covariates (known)
  int K;                      // num sources (NOT known here)

  array[N] int T_i;           // num observations x staz (known)
  array[N,T] int idx_ti;      // observations' index x staz (known)

  array[N,T] int C_ti;        // observations' index x staz (known)
  array[N,T,C] int idx_cti;   // observations' index x staz (known)

  array[N] matrix[C,T] Y;     // observed values (with NAs)
  matrix[N,p] X;              // covariates
  matrix[M,T] B;              // splines matrix
  matrix[N,N] D;              // distance matrix for sites
  vector[C] eps;                   // shrinkage threshold
}



parameters {
  vector<lower=0>[C] sigma;

  array[K] vector[p] beta;

  matrix[K,N] coeff_g;
  array[K] simplex[C] H;
  vector<lower=0>[K] r;

  matrix<lower=0>[K,M] L;
  matrix<lower=0>[M,K] phi;
  vector<lower=0>[K] delta;
}



transformed parameters {
  matrix[K,T] f = L*B;

  array[K] matrix[N,T] g;
  vector<lower=0>[K] eta;

  real count=0;
  vector[K] k_idx = rep_vector(0,K);
  for(k in 1:K){
      eta[k] = prod(delta[1:k]);
      for (i in 1:N) { g[k][i] = exp(coeff_g[k,i]) * f[k]; }
      if((max(H[k, ]) * max(g[k])) > 1) count += 1; 
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
  r ~ inv_gamma(3,1000);

  //MGPS hyperparameters
  to_vector(phi) ~ gamma(1.5, 1.5);
  delta[1] ~ gamma(10, 1);
  delta[2:K] ~ gamma(20, 1);

  for(k in 1:K){
    //beta[k] ~ normal(m0, s0);
    beta[k] ~ normal(0, 1);
    L[k, ] ~ normal(0, sqrt(1/(phi[,k] * eta[k]))) T[0, ];
    coeff_g[k] ~ multi_normal((X*beta[k])', exp(-D/r[k]));
    H[k] ~ dirichlet(rep_vector(1, C));
  }

  for (i in 1:N){
    for (t in 1:T_i[i]){
      for (c in 1:C_ti[i,idx_ti[i,t]]){
        Y[i, idx_cti[i,idx_ti[i,t],c], idx_ti[i,t]] ~ normal(log(mu_y[i, idx_cti[i,idx_ti[i,t],c], idx_ti[i,t]]), sigma[idx_cti[i,idx_ti[i,t],c]]);
  }}}

}



generated quantities  {
  array[N] matrix[C,T] log_lik;
  array[N] matrix[C,T] Y_pred;

  for (i in 1:N){
    Y_pred[i] = rep_matrix(0, C, T);
    for (t in 1:T){
      for (c in 1:C){
        log_lik[i,idx_cti[i,t,c],idx_ti[i,t]] = normal_lpdf(Y[i,idx_cti[i,t,c],idx_ti[i,t]] | log(mu_y[i, idx_cti[i,idx_ti[i,t],c], idx_ti[i,t]]), sigma[idx_cti[i,idx_ti[i,t],c]]);
        Y_pred[i,c,t] = normal_rng(log(mu_y[i,c,t]), sigma[c]);
      }}}
}



