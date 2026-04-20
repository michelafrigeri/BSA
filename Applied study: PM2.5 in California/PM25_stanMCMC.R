#-----------------------------------------------------------------#
##################       Applied PM2.5 MCMC       #################
#-----------------------------------------------------------------#


# R packages --------------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(rstan)


# Loading simulation data -------------------------------------------------
load("final_PM25_data.RData")


# Iterative BSA with Stan ----------------------------------------------------------
Y = Y_obs
Y[is.na(Y)] = -100
K_star = 10
K_curr = 0
M = numbasis

# Model
bsa_model = stan_model("bsa_PM25.stan") 

bsa_data = list(
  T= T,          # num observations (known)
  M= M,          # num basis (fixed)
  C= C,          # num of pollutants (known) 
  N= N,          # num of sites (known)
  p= ncol(X),    # num covariates (known)
  K= K_star,      
  
  T_i = T_i,
  idx_ti = t_index,
  C_ti = C_ti,
  idx_cti = idx_cti,
  
  Y= Y,          # observed values
  X= X,          # covariates
  B= B,          # splines matrix
  D= D           # distance matrix for sites
)


init_list = "random"
while (K_star != K_curr) {
  bsa_data$K = K_star
  bsa_fit = sampling(
    bsa_model,
    data = bsa_data,
    #init = init_list,
    seed = 1997,
    chains = 1,
    warmup = 1000,
    iter = 1050,
    refresh = 350
  )
  fit = bsa_fit
  K_curr = K_star
  K_star = as.integer(median(as.array(fit, pars = c("count"))))
  cat(sprintf("K_curr: %g, K_star: %g\n", K_curr, K_star))

  init_list = list(
    list(sigma = extract(fit, pars = c("sigma"), permuted=F)[50, , ],
    beta = matrix(extract(fit, pars = c("beta"), permuted=F)[50, , ], nrow=K_curr, ncol=4)[k_sel, ],
    coeff_g = matrix(extract(fit, pars = c("coeff_g"), permuted=F)[50, , ], nrow=K_curr, ncol=N)[k_sel, ],
    H = matrix(extract(fit, pars = c("H"), permuted=F)[50, , ], nrow=K_curr, ncol=C)[k_sel, ],
    r = extract(fit, pars = c("r"), permuted=F)[50, ,k_sel],
    L = matrix(extract(fit, pars = c("L"), permuted=F)[50, , ], nrow=K_curr, ncol=M)[k_sel, ],
    phi = matrix(extract(fit, pars = c("phi"), permuted=F)[50, , ], nrow=M, ncol=K_curr)[ ,k_sel],
    delta = extract(fit, pars = c("delta"), permuted=F)[50, ,k_sel])
  )
}
#2: Longer MCMC with the number of sources K* we just retrieved
bsa_data$K = K_star
bsa_fit = sampling(
  bsa_model,
  data = bsa_data,
  init = init_list,
  seed = 1997,
  chains = 1,
  warmup = 6000,
  iter = 10000,
  refresh = 1000
)

# Save MCMC
saveRDS(bsa_fit, file = "PM25_fit.RDS")
