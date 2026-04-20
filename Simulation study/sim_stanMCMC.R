#-----------------------------------------------------------------#
##################         SIMULATION MCMC        #################
#-----------------------------------------------------------------#


# R packages --------------------------------------------------------------
pacman::p_load(rstan, bayesplot)


# Loading simulation data -------------------------------------------------
load("sim_data10_new.RData")   
#load("sim_data20_new.RData")
#load("sim_data30_new.RData")
#load("sim_data50_new.RData")

# Iterative BSA with Stan ----------------------------------------------------------
Y = Y_obs
Y[is.na(Y)] = -1000
K_star = 10
K_curr = 0
M=numbasis

# Model
bsa_model = stan_model("bsa_sim.stan") 


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
  eps= apply(exp(Y_obs),2,min,na.rm=T)      # shrinkage threshold
  #eps = round(min( exp( - (apply(Y_obs,2,sd,na.rm=T)^2)/2 ))/10,3)
)




# MCMC sampling
#1: Few iterations during sources' shrinkage to find K*
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
  
  L = mcmc_intervals_data(as.matrix(fit, pars = c("L")))
  L_median = matrix(L$m, nrow = K_curr, ncol = numbasis)
  K_star_prova = sum(apply(L_median, 1, max)>0.05)
    
  cat(sprintf("K_curr: %g, K_star: %g, K_star_new: %g\n", K_curr, K_star, K_star_prova))
  init_list = list(
    list(sigma = extract(fit, pars = c("sigma"), permuted=F)[50, , ],
         beta = matrix(extract(fit, pars = c("beta"), permuted=F)[50, , ], nrow=K_curr, ncol=4)[1:K_star, ],
         coeff_g = matrix(extract(fit, pars = c("coeff_g"), permuted=F)[50, , ], nrow=K_curr, ncol=N)[1:K_star, ],
         H = matrix(extract(fit, pars = c("H"), permuted=F)[50, , ], nrow=K_curr, ncol=C)[1:K_star, ],
         r = extract(fit, pars = c("r"), permuted=F)[50, , 1:K_star],
         L = matrix(extract(fit, pars = c("L"), permuted=F)[50, , ], nrow=K_curr, ncol=M)[1:K_star, ],
         phi = matrix(extract(fit, pars = c("phi"), permuted=F)[50, , ], nrow=M, ncol=K_curr)[ ,1:K_star],
         delta = extract(fit, pars = c("delta"), permuted=F)[50, ,1:K_star])
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
saveRDS(bsa_fit, file = "sim_fit10_new.RDS")
#saveRDS(bsa_fit, file = "sim_fit20_new.RDS")
#saveRDS(bsa_fit, file = "sim_fit30_new.RDS")
#saveRDS(bsa_fit, file = "sim_fit50_new.RDS")


