# R packages --------------------------------------------------------------
library(splines)
library(MASS)
library(ggplot2)
library(LDATS)
library(readr)
library(fastDummies)
library(geosphere)
library(gridExtra)

# Reduced dataset ---------------------------------------------------------
# only one year, 32 stations and 6 pollutants
data = read_csv("realData_subset.csv")
data = data[order(data$SiteName, data$t), ]
names(data)


# alcuni negativi?? per ora metto a zero ma CHIEDERE
min(data[,3:8], na.rm=T)
matplot(data[,3:8], type='l')
prova = data[,3:8] +0.5
matplot(log(prova), type='l')

# Trasformiamo i dati in scala log
data[ ,3:8] = log(data[ ,3:8]+0.5)

C = 6
stations = sort(unique(data$SiteName))
N = length(stations)

t_basis = sort(unique(data$t))
T = length(t_basis)
T_i = tapply(data$t, sort(data$SiteName), length)

# Bsplines
set.seed(2024)
B = t(bs(t_basis, knots=quantile(data$t, p=seq(0,1,length=16)), degree=3, intercept = TRUE))
matplot(t(B), type='l')
numbasis = nrow(B)
numbasis

t_obs = tapply(data$t, data$SiteName, sort)
t_index = matrix(rep(0,N*T), nrow=N, ncol=T)
for (i in 1:N) {
  t_index[i, 1:T_i[i]] = match(t_obs[[i]], t_basis)
}


sum(rowSums(is.na(data)))


# creating data in the right format
Y_obs = array(data=rep(NA, N*C*T), dim=c(N,C,T))
for (i in 1:N) {
  temp_staz = data[data$SiteName==stations[i], 3:8]
  for (c in 1:C) {
    for (t in 1:T_i[i]) {
      Y_obs[i, c, t_index[i,t]] = as.numeric(temp_staz[t,c])
    }}
}


C_ti = matrix(rep(0,N*T), nrow=N, ncol=T)
for (i in 1:N) {
  for (t in 1:T_i[i]) {
    C_ti[i,t_index[i,t]] = sum(!is.na(Y_obs[i, ,t_index[i,t]]))
  }}
idx_cti = array(data=rep(0, N*C*T), dim=c(N,T,C))
for (i in 1:N) {
  for (t in 1:T_i[i]) {
    temp_yi = Y_obs[i, ,t_index[i,t]]
    idx_cti[i, t_index[i,t], 1:C_ti[i,t_index[i,t]]] = which(!is.na(temp_yi))
  }
}



# Covariates
## Real covariates and locations
staz = read_csv("FinalStations.csv")
staz$LandSetting = as.factor(staz$LandSetting)
staz = staz[order(staz$`Local Site Name`), ]
dummies = fastDummies::dummy_columns(staz[,"LandSetting"])
# baseline: Rural and Mobile
dummies = dummies[ , -c(1:2)] 
X = cbind(dummies, std_altitude=scale(staz$Elevation), intercept=1)
names(X)


# Spatial structure
# Simulating sources at i=1,..,10 different locations: g_i(t)
latlong = staz[ ,2:1]
D = distm(latlong, latlong, fun = distHaversine)/1000 
D = round(D,2)
max(D)/3

save.image(file="final_PM25_data.RData")

# Check on non-NA indexing for Stan
Y = Y_obs
Y[is.na(Y)] = -1000
yobs = array(data=rep(NA, N*C*T), dim=c(N,T,C))
for (i in 1:N) {
  for (t in 1:T_i[i]) {
    for (c in 1:C_ti[i,t_index[i,t]]) {
      yobs[i,t,c] = Y[i, idx_cti[i,t_index[i,t],c],t_index[i,t]]
    }
  }
}
range(yobs, na.rm = T)

