#-----------------------------------------------------------------#
##################         SIMULATION DATA        #################
#-----------------------------------------------------------------#

# R packages --------------------------------------------------------------
library(readr)
library(splines)
library(MASS)
library(ggplot2)
library(gridExtra)
library(fastDummies)
library(geosphere)
library(LDATS)
library(reshape2)
library(extraDistr)


# Simulation study: REALISTIC DATA ------------------------------

## Using EPA AQS monitoring network: real covariates and locations
staz = read_csv("FinalStations.csv")
staz$LandUse = as.factor(staz$LandUse)
staz$LandSetting = as.factor(staz$LandSetting)
staz = staz[order(staz$`Local Site Name`), ]

# Considering only "Elevation" and "Land setting" as regressors
# -> expressing "Land setting" as dummy variable
dummies = fastDummies::dummy_columns(staz[ ,"LandSetting"])
# baseline: Rural
dummies = dummies[ , -c(1:2)] 
# covariates matrix
X = cbind(dummies, std_altitude=scale(staz$Elevation), intercept=1)
names(X)

# Dimensionality of the synthetic data
N = nrow(staz)  #same as real data
C = 6           #same as real data
K = 2
# for reproducible analysis
seed = 1604

# Simulating cubic B-splines for longitudinal data over T_i
t_obs = list(NULL)
T_i = rep(0,N)
set.seed(seed)
n_obs = sample(100:150, N)
for (i in 1:N) {
  t_obs[[i]] = sort(sample(seq(1,365,by=2), n_obs[i], replace=F))
  T_i[i] = length(t_obs[[i]])
}

t_union = unlist(t_obs)
t_basis = sort(unique(t_union))
T = length(t_basis)

# cubic B-splines
M_knots = quantile(t_union, p=seq(0,1,by=0.1))
B = t(bs(t_basis, knots=M_knots, degree=3, intercept = T))
matplot(t(B), type='l')
numbasis = nrow(B)
numbasis

t_index = matrix(rep(0,N*T), nrow=N, ncol=T)
for (i in 1:N) { t_index[i, 1:T_i[i]] = match(t_obs[[i]], t_basis) }


### Generating global sources' profiles: f_k(t) ----------------------
set.seed(seed)
#w = rbind(mvrnorm(1, rep(.8,numbasis), diag(runif(numbasis,0,.5))), mvrnorm(1, rep(.5,numbasis), diag(runif(numbasis,0,1))))
w = mvrnorm(2, rep(1,numbasis), diag(rep(1,numbasis)))
f = t(w%*%B)

matplot(t_basis, f, type='l', lwd=3, lty=1,ylab="", xlab="")



### Generating local sources' profiles: g_ki(t) ----------------------
latlong = staz[ ,2:1]
D = distm(latlong, latlong, fun = distHaversine)/1000 
D = round(D,3)
max(D)/3
r1 = 600 
r2 = 300
H1 = round(exp(-D/r1), 3)
H2 = round(exp(-D/r2), 3)

names(X)
betas = rbind(c(0.1, 0.2, -0.3, 0.4),
              c(-0.5, -0.4, -0.2, -0.2))
set.seed(seed)
coeff_g = t(rbind(mvrnorm(1, as.vector(betas[1,]%*%t(X)), H1), 
                  mvrnorm(1, as.vector(betas[2,]%*%t(X)), H2)))

g = array(rep(0, T*N*K), dim=c(N,K,T))
for (k in 1:K) {
  for (i in 1:N) {
    g[i,k, ] = exp(coeff_g[i, k]) * f[ ,k]
  }}

#plot
df_g = as.data.frame(t(g[ ,2, ]))
df_f = data.frame(values = f[ ,2], variable=factor(1))
df_g$id = df_f$id = t_basis
plot_data = melt(df_g,id.var="id")
g1 = ggplot(plot_data, aes(x=id,y=value,group=variable,colour=variable)) +
  geom_line(aes(lty=factor(as.numeric(variable)%%4)), lwd=.65) + ggtitle("Source 1") + 
  scale_colour_grey(start=0.35, end=0.8) + theme_light() + labs(x = NULL, y = NULL) +
  geom_line(data = df_f, aes(x=id, y=values, group=variable), lty=1, col=1, lwd=1.5) +
  theme(plot.title=element_text(hjust=.5), text=element_text(size=20), legend.position="none") +
  coord_cartesian(xlim=c(-1,366), ylim=c(-0.9,4.5), expand=F)
df_g = as.data.frame(t(g[ ,1, ]))
df_f = data.frame(values = f[ ,1], variable=factor(1))
df_g$id = df_f$id = t_basis
plot_data = melt(df_g,id.var="id")
g2 = ggplot(plot_data, aes(x=id,y=value,group=variable,colour=variable)) +
  geom_line(aes(lty=factor(as.numeric(variable)%%4)), lwd=.65) + ggtitle("Source 2") + 
  scale_colour_grey(start=0.35, end=0.8) + theme_light() + labs(x = NULL, y = NULL) +
  geom_line(data = df_f, aes(x=id, y=values, group=variable), lty=1, col=1, lwd=1.5) +
  theme(plot.title=element_text(hjust=.5), text=element_text(size=20), legend.position="none")+
  coord_cartesian(xlim=c(-1,366), ylim=c(-0.05,4.5), expand=F)
grid.arrange(g1, g2, ncol=2)




### Generating the pollutants' concentrations: y_ic(t) ------------------
set.seed(seed+1)
h = rdirichlet(2,rep(1,C))
coeff_c = round(h,3)
coeff_c
rowSums(coeff_c)

y_true = array(data=rep(0,N*C*T), dim=c(N,C,T))
for (i in 1:N) {
  for (c in 1:C) {
    y_true[i,c, ] = g[i,1, ]*coeff_c[1,c] + g[i,2, ]*coeff_c[2,c]
  }}
range(y_true)

par(mfrow=c(3,2))
for (c in 1:C) { matplot(t(y_true[ ,c, ]),type='l', lwd=2) }



### Adding measurement noise ----------------------------------------------
sigma_c = round(apply(y_true, 2, sd)/sqrt(10), 3)
sigma_c
set.seed(seed)
Y_final = y_true 
for (c in 1:C) { Y_final[,c,] = Y_final[ ,c, ] + rnorm(T*N, 0, sigma_c[c]) }

par(mfrow=c(3,2))
for (c in 1:C) { matplot(t(Y_final[ ,c, ]),type='l', lwd=.7) }
range(Y_final)



### Synthetic longitudinal data: noisy observations with missing values ----
Y_obs = array(data=rep(NA, N*C*T), dim=c(N,C,T))
for (i in 1:N) { Y_obs[i, ,t_index[i, 1:T_i[i]]] = Y_final[i, ,t_index[i, 1:T_i[i]]] }



# Saving RData ------------------------------------------------------------

save.image(file="final_simulation_data.RData")




