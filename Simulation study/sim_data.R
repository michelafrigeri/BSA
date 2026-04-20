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
library(truncnorm)


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
seed = 416

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

#w = rbind(mvrnorm(1, rep(.8,numbasis), diag(runif(numbasis,0,.5))), mvrnorm(1, rep(.5,numbasis), diag(runif(numbasis,0,1))))
#w = mvrnorm(2, rep(1,numbasis), diag(rep(1,numbasis)))
set.seed(seed)
w = rbind(rtruncnorm(n=numbasis,a=0,mean=0,sd=10), rtruncnorm(n=numbasis,a=0,mean=0,sd=5))

f = t(w%*%B)

matplot(t_basis, f, type='l', lwd=3, lty=1,ylab="", xlab="")



### Generating local sources' profiles: g_ki(t) ----------------------
latlong = staz[ ,2:1]
D = distm(latlong, latlong, fun = distHaversine)/1000 
D = round(D,3)
max(D)/3
r1 = 300 
r2 = 900
H1 = round(exp(-D/r1), 3)
H2 = round(exp(-D/r2), 3)

names(X)
#betas = rbind(c(0.1, 0.2, -0.3, 0.4),
#              c(-0.5, -0.4, -0.2, -0.2))
betas = rbind(c(1.25, 1, -0.1, -0.5),
              c(-1.5, -1, 0.2, 0.25))
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
  #coord_cartesian(xlim=c(-1,366), ylim=c(-0.9,4.5), expand=F)
  coord_cartesian(xlim=c(-1,366), expand=F)
df_g = as.data.frame(t(g[ ,1, ]))
df_f = data.frame(values = f[ ,1], variable=factor(1))
df_g$id = df_f$id = t_basis
plot_data = melt(df_g,id.var="id")
g2 = ggplot(plot_data, aes(x=id,y=value,group=variable,colour=variable)) +
  geom_line(aes(lty=factor(as.numeric(variable)%%4)), lwd=.65) + ggtitle("Source 2") + 
  scale_colour_grey(start=0.35, end=0.8) + theme_light() + labs(x = NULL, y = NULL) +
  geom_line(data = df_f, aes(x=id, y=values, group=variable), lty=1, col=1, lwd=1.5) +
  theme(plot.title=element_text(hjust=.5), text=element_text(size=20), legend.position="none")+
  #coord_cartesian(xlim=c(-1,366), ylim=c(-0.05,4.5), expand=F)
  coord_cartesian(xlim=c(-1,366),  expand=F)
grid.arrange(g1, g2, ncol=2)




### Generating the pollutants' concentrations: y_ic(t) ------------------
set.seed(seed)
h = rdirichlet(2,rep(10.1,C))
coeff_c = round(h,3)
coeff_c
coeff_c[2,4] = 0.06
coeff_c[2,3] = 0.272
rowSums(coeff_c)

#coeff_c= matrix(c(0.05, 0.25, 0.15, 0.3, 0.05, 0.2,
#                  0.03, 0.15, 0.2, 0.07, 0.25, 0.3), nrow=2, ncol=6, byrow=T)
rowSums(coeff_c)
matplot(t(coeff_c))


y_true = array(data=rep(0,N*C*T), dim=c(N,C,T))
for (i in 1:N) {
  for (c in 1:C) {
    y_true[i,c, ] = log(g[i,1, ]*coeff_c[1,c] + g[i,2, ]*coeff_c[2,c])
  }}
range(y_true)

par(mfrow=c(3,2))
for (c in 1:C) { matplot(t(y_true[ ,c, ]),type='l', lwd=2) }



### Adding measurement noise ----------------------------------------------
sigma_c = apply(y_true, 2, sd)*sqrt(0.3)
sigma_c
sigma_c^2

set.seed(seed)
Y_final = y_true 
for (c in 1:C) { Y_final[,c,] = Y_final[ ,c, ] + rnorm(T*N, 0, sigma_c[c]) }

par(mfrow=c(3,2))
for (c in 1:C) { matplot(t(Y_final[ ,c, ]),type='l', lwd=.7) }
apply(Y_final, 2, range)
apply(exp(Y_final), 2, range)


### Synthetic longitudinal data: noisy observations with missing values ----
Y_obs = array(data=rep(NA, N*C*T), dim=c(N,C,T))
for (i in 1:N) { Y_obs[i, ,t_index[i, 1:T_i[i]]] = Y_final[i, ,t_index[i, 1:T_i[i]]] }



# Saving RData ------------------------------------------------------------
save.image(file="sim_data30_new.RData")

load(file="sim_data10_new.RData")
apply(exp(Y_obs),2,sd,na.rm=T)^2/diff(apply(exp(Y_obs),2,range,na.rm=T))
#cdv
min(apply(exp(Y_obs),2,sd,na.rm=T)/abs(apply(exp(Y_obs),2,median,na.rm=T)))

quantile(exp(Y_obs), p=0.25, na.rm=T) - 1.5*IQR(exp(Y_obs),na.rm=T)

round( min( exp( apply(Y_obs,2,quantile,p=0.25,na.rm=T) - 1.5*apply(Y_obs,2,IQR,na.rm=T) - (apply(Y_obs,2,sd,na.rm=T)^2)/2 ))/2,3)

#epsilon approx zero
round(min( exp( apply(Y_obs,2,min,na.rm=T) - (apply(Y_obs,2,sd,na.rm=T)^2)/2 ))/2,3)
#real data
round(min( exp( apply(log(exp(Y_obs)-0.4),2,min,na.rm=T) - (apply(log(exp(Y_obs)-0.4),2,sd,na.rm=T)^2)/2 ))/2,3)


# Color plots -------------------------------------------------------------

#plot
df_g = as.data.frame(t(g[ ,2, ]))
df_f = data.frame(values = f[ ,2], variable=factor(1))
df_g$id = df_f$id = t_basis
plot_data = melt(df_g,id.var="id")
g1 = ggplot(plot_data, aes(x=id,y=value,group=variable,colour=variable)) +
  geom_line(aes(lty=factor(as.numeric(variable)%%4)), lwd=.65) + ggtitle("Source 1") + 
  scale_colour_hue() + theme_light() + labs(x = NULL, y = NULL) +
  geom_line(data = df_f, aes(x=id, y=values, group=variable), lty=1, col=1, lwd=2) +
  theme(plot.title=element_text(hjust=.5), text=element_text(size=20), legend.position="none") +
  coord_cartesian(xlim=c(-1,366), ylim=c(-0.9,4.5), expand=F)
df_g = as.data.frame(t(g[ ,1, ]))
df_f = data.frame(values = f[ ,1], variable=factor(1))
df_g$id = df_f$id = t_basis
plot_data = melt(df_g,id.var="id")
g2 = ggplot(plot_data, aes(x=id,y=value,group=variable,colour=variable)) +
  geom_line(aes(lty=factor(as.numeric(variable)%%4)), lwd=.65) + ggtitle("Source 2") + 
  scale_colour_hue() + theme_light() + labs(x = NULL, y = NULL) +
  geom_line(data = df_f, aes(x=id, y=values, group=variable), lty=1, col="red4", lwd=2) +
  theme(plot.title=element_text(hjust=.5), text=element_text(size=20), legend.position="none")+
  coord_cartesian(xlim=c(-1,366), ylim=c(-0.05,4.5), expand=F)
grid.arrange(g1, g2, ncol=2)
g1
g2

ggplot(data = data.frame(values = f[ ,1], variable=factor(1), id = t_basis), aes(x=id, y=values, group=variable)) +
  geom_line(lty=1, col="red3", lwd=3) + theme_light() + labs(x = NULL, y = NULL) +
  geom_line(data = data.frame(values = f[ ,2], variable=factor(1), id = t_basis), aes(x=id, y=values, group=variable), lty=1, col=1, lwd=3) +
  theme(plot.title=element_text(hjust=.5), text=element_text(size=20), legend.position="none")




#load(file="final_simulation_data.RData")
#load(file="simulation_data_var20.RData")
#load(file="simulation_data_var30.RData")
load(file="sim_data30_new.RData")
plot_list <- list() 
lims = list(c(-0.6,1.9), c(-0.6,1.5), c(-0.6,1.5), c(-0.05,0.18), c(-0.6,1.5), c(-0.35,1.3))
ymaxs=ymins=ymeds= rep(0,T)
for (c in 1:C){ 
  for (tt in 1:T){
    curr_data = Y_final[ ,c, tt]
    ymaxs[tt]=max(curr_data, na.rm=T)
    ymeds[tt]=median(curr_data, na.rm=T)
    ymins[tt]=min(curr_data, na.rm=T)}
    maxmin_data = data.frame(ymaxs,ymeds,ymins)
    plot_list[[c]] = ggplot(maxmin_data, aes(x=t_basis,y=ymeds)) +
      geom_ribbon(aes(ymin = ymins, ymax = ymaxs), fill='gray', alpha = 0.3) +
      geom_line(data=data.frame(y=Y_final[23,c, ], x=t_basis), aes(y=y, x=x), color="gray60") +
      #geom_point(data= data.frame(y=Y_final[4,c, ], x=1:T), aes(y=y, x=x), color="gray60") +
      geom_line(data= data.frame(y=Y_final[31,c, ], x=t_basis), aes(y=y, x=x), color="gray40") +
      #geom_point(data= data.frame(y=Y_final[9,c, ], x=1:T), aes(y=y, x=x), color="gray40") +
      geom_line(lty=2, lwd=.8) + ggtitle(paste0("Pollutant ",c)) + theme_light() + labs(x = NULL, y = NULL) + 
      theme(text=element_text(size=20), legend.position="none", plot.title = element_text(hjust = 0.5)) +
      ylim(lims[[c]])
  }
grid.arrange(grobs=plot_list, nrow=2)


