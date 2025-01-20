#-----------------------------------------------------------------#
##################       Applied PM2.5 MCMC       #################
#-----------------------------------------------------------------#


# R packages --------------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(rstan)


# Loading simulation data -------------------------------------------------
load("final_PM25_data.RData")


# Iterative BSA with Stan ------------------------------------------------------
Y = Y_obs
Y[is.na(Y)] = -100
K_star = 10
K_curr = 0
M = numbasis

# Model
bsa_model = stan_model("bsa_PM25_dirichlet.stan") 

bsa_data = list(
  T = T,          # num observations (known)
  M = M,          # num basis (fixed)
  C = C,          # num of pollutants (known) 
  N = N,          # num of sites (known)
  p = ncol(X),    # num covariates (known)
  K = K_star,      
  T_i = T_i,
  idx_ti = t_index,
  C_ti = C_ti,
  idx_cti = idx_cti,
  Y = Y,          # observed values
  X = X,          # covariates
  B = B,          # splines matrix
  D = D,          # distance matrix for sites
  eps = 0.05
)


init_list = "random"
while (K_star != K_curr) {
  bsa_data$K = K_star
  bsa_fit = sampling(
    bsa_model,
    data = bsa_data,
    seed = 1997,
    chains = 1,
    warmup = 500,
    iter = 550,
    refresh = 50
  )
  fit = bsa_fit
  K_curr = K_star
  K_star = as.integer(median(as.array(fit, pars = c("count"))))
  cat(sprintf("K_curr: %g, K_star: %g\n", K_curr, K_star))
}
#2: Longer MCMC with the number of sources K* we just retrieved
bsa_data$K = K_star
bsa_fit = sampling(
  bsa_model,
  data = bsa_data,
  seed = 1997,
  chains = 1,
  warmup = 6000,
  iter = 10000,
  refresh = 1000
)

# Save MCMC
saveRDS(bsa_fit, file = "PM25_fit_dirichlet.RDS")
