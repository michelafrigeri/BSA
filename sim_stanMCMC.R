#-----------------------------------------------------------------#
##################         SIMULATION MCMC        #################
#-----------------------------------------------------------------#


# R packages --------------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(rstan)


# Loading simulation data -------------------------------------------------
load("final_simulation_data.RData")


# Iterative BSA with Stan ----------------------------------------------------------
Y = Y_obs
Y[is.na(Y)] = -100
K_star = 10
K_curr = 0

# Model
bsa_model = stan_model("stan/bsa_dirichlet.stan") 
#bsa_model = stan_model("stan/bsa_model.stan") 

# Data
bsa_data = list(
  T=T,          # num total observations (known)
  M=numbasis,   # num functional basis (fixed)
  C=C,          # num of pollutants (known) 
  N=N,          # num of locations (known)
  p=ncol(X),    # num of covariates (known)
  K=K_star,     # num sources (unknown)
  
  T_i=T_i,
  idx_ti=t_index,
  
  Y=Y,          # observed values
  X=X,          # covariates
  B=B,          # basis splines matrix
  D=D,          # sites' distance matrix
  eps=0.015      # shrinkage threshold
)

# MCMC sampling
#1: Few iterations during sources' shrinkage to find K*
while (K_star != K_curr) {
  bsa_data[[K]] = K_star
  bsa_fit = sampling(
    bsa_model,
    data = bsa_data,
    seed = 2024,
    chains = 1,   
    warmup = 900,          
    iter = 1000,
    refresh = 100
  )
  fit = bsa_fit
  K_curr = K_star
  K_star = as.integer(median(rstan::extract(fit)$count))
}
#2: Longer MCMC with the number of sources K* we just retrieved
bsa_data[[K]] = K_star
bsa_fit = sampling(
  bsa_model,
  data = bsa_data,
  seed = 2024,
  chains = 2,   
  warmup = 3000,          
  iter = 5000,
  refresh = 500
)

# Save MCMC
saveRDS(bsa_fit, file = "bsa_fit_simulation.RDS")



