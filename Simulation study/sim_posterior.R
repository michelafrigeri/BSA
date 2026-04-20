#-----------------------------------------------------------------#
##################       POSTERIOR INFERENCE      #################
#-----------------------------------------------------------------#


# R packages --------------------------------------------------------------
library(MASS)
library(rstan)
library(bayesplot)
library(ggplot2)
library(gridExtra)
library(readr)
library(xtable)


# Loading data and MCMC ---------------------------------------------------
load("sim_data10_new.RData")
# load("sim_data20_new.RData")
# load("sim_data30_new.RData")
# load("sim_data50_new.RData")

fit10 = read_rds("sim_fit10_new.RDS")
fit20 = read_rds("sim_fit20_new.RDS")
fit30 = read_rds("sim_fit30_new.RDS")
fit50 = read_rds("sim_fit50_new.RDS")

col_light = c("deepskyblue", "yellowgreen", "orange", "deeppink")
col_dark = c("deepskyblue3", "forestgreen", "darkorange3", "deeppink2")
col_points = c("blue3", "darkgreen", "darkorange4", "deeppink3")


# Shrinkage check: L and K* ------------------------------------------------------
fit = fit10
K_post = as.array(fit, pars = c("count"))
mcmc_trace(K_post)
K_star = median(as.array(fit, pars = c("count")))
mean(as.array(fit, pars = c("count")))

L = mcmc_intervals_data(as.matrix(fit, pars = c("L")))
L_mean = matrix(L$m, nrow = K_star, ncol = numbasis)
matplot(t(L_mean))
abline(h=c(0,0.02,0.05), lty=2)

L_sup = matrix(L$h, nrow = K_star, ncol = numbasis)
matplot(t(L_sup))
#matplot(abs(t(L_median)))
abline(h=c(0.02), lty=2)

mcmc_intervals(as.array(fit, pars = c("delta")))


# Global source contributions: f_k --------------------------------------------
f_post = as.array(fit, pars = c("f"))
f_ci = mcmc_intervals_data(f_post, prob=0.95, prob_outer = 1)


f_1 = f_ci[seq(1,T*K_star,by=K_star), ]
f_1$t = t_basis
f_2 = f_ci[seq(1,T*K_star,by=K_star)+1, ]
f_2$t = t_basis

ggplot(data=f_1, aes(x=t)) +
  theme_light() + geom_line(aes(y=m), color = 'blue', lwd=.8)  +
  geom_ribbon(data=f_1, aes(ymin = l, ymax = h), fill='deepskyblue', alpha = 0.3) +
  geom_line(aes(x=t,y=f[ ,1]), color = 'red', lwd=.8)  +
  geom_line(data=f_2, aes(y=m), color = 'forestgreen', lwd=.8)  +
  geom_ribbon(data=f_2, aes(ymin = l, ymax = h), fill='yellowgreen', alpha = 0.3) +
  geom_line(aes(x=t,y=f[ ,2]), color = 'red', lwd=.8)  +
  labs(x = NULL, y = NULL)

# Standardized f_k
L = mcmc_intervals_data(as.matrix(fit, pars = c("L")))
L_mean = matrix(L$m, nrow = K_star, ncol = numbasis)
f1_std = f_1[,c("l","m","h")]/sqrt(sum(L_mean[1,]^2))
f2_std = f_2[,c("l","m","h")]/sqrt(sum(L_mean[2,]^2))
fstd = cbind(f[,1]/sqrt(sum(w[1,]^2)), f[,2]/sqrt(sum(w[2,]^2)))
f1_std$t = f2_std$t = t_basis

ggplot(data=f1_std, aes(x=t)) +
  theme_light() + geom_line(aes(y=m), color = 'blue', lwd=.8)  +
  geom_ribbon(data=f1_std, aes(ymin = l, ymax = h), fill='deepskyblue', alpha = 0.3) +
  geom_line(data=f2_std, aes(y=m), color = 'forestgreen', lwd=.8)  +
  geom_ribbon(data=f2_std, aes(ymin = l, ymax = h), fill='yellowgreen', alpha = 0.3) +
  geom_line(aes(x=t,y=fstd[ ,1]), color = 'red3', lwd=.8)  +
  geom_line(aes(x=t,y=fstd[ ,2]), color = 'red3', lwd=.8)  +
  labs(x = NULL, y = NULL)




# Regression coefficients: beta_k -----------------------------------------
fit = fit10
betas_post = as.array(fit, pars = c("beta"))
betas_ci = mcmc_intervals_data(betas_post, prob=0.95, prob_outer = 1)
beta1 = betas_ci[seq(1,4*K_star,by=K_star), c(1,6:8)]
beta2 = betas_ci[seq(1,4*K_star,by=K_star)+1, c(1,6:8)]
mcmc_trace(betas_post)


beta_table = data.frame(
  suburb = rbind(sprintf('%s (%s)', toString(round(beta1[1,3],2)), toString(round(beta1[1,c(2,4)], 2))),
                 sprintf('%s (%s)', toString(round(beta2[1,3],2)), toString(round(beta2[1,c(2,4)], 2)))),
  urban = rbind( sprintf('%s (%s)', toString(round(beta1[2,3],2)), toString(round(beta1[2,c(2,4)], 2))),
                 sprintf('%s (%s)', toString(round(beta2[2,3],2)), toString(round(beta2[2,c(2,4)], 2)))),
  altitude = rbind(sprintf('%s (%s)', toString(round(beta1[3,3],2)), toString(round(beta1[3,c(2,4)], 2))),
                   sprintf('%s (%s)', toString(round(beta2[3,3],2)), toString(round(beta2[3,c(2,4)], 2)))),
  intercept = rbind(sprintf('%s (%s)', toString(round(beta1[4,3],2)), toString(round(beta1[4,c(2,4)], 2))),
                   sprintf('%s (%s)', toString(round(beta2[4,3],2)), toString(round(beta2[4,c(2,4)], 2)))))
xtable(beta_table)


# Standardized betas
beta1_std = rbind(beta1[1:3,c("l","m","h")], 
                  beta1[4,c("l","m","h")] + log(sqrt(sum(L_mean[1,]^2))))
beta2_std = rbind(beta2[1:3,c("l","m","h")], 
                  beta2[4,c("l","m","h")] + log(sqrt(sum(L_mean[2,]^2))))
beta1_std$parameter = beta2_std$parameter = factor(beta1$parameter)
betas_std = cbind(c(betas[1,1:3],betas[1,4]+log(sqrt(sum(w[1,]^2)))), 
                  c(betas[2,1:3],betas[2,4]+log(sqrt(sum(w[2,]^2)))))

ggplot(beta2_std, aes(x=parameter, y=m)) + #color=parameter
  geom_errorbar(aes(ymin=l, ymax=h), width=.1, lwd=0.8, colour='black') +
  geom_point(size=5) + geom_point(size=5, pch=21, color='black', fill="deepskyblue") + 
  ylab(" ") + xlab(" ") +
  geom_hline(yintercept = 0.0 , colour='red4', lwd=0.7) +  
  geom_point(aes(y=betas_std[,2]), size=3, pch=23, color='black', fill='red2') +
  theme_light() + legend_none() +
  theme(text = element_text(size=20)) 
ggplot(beta1_std, aes(x=parameter, y=m)) + #color=parameter
  geom_errorbar(aes(ymin=l, ymax=h), width=.1, lwd=0.8, colour='black') +
  geom_point(size=5) + geom_point(size=5, pch=21, color='black', fill="deepskyblue") + 
  ylab(" ") + xlab(" ") +
  geom_hline(yintercept = 0.0 , colour='red4', lwd=0.7) +  
  geom_point(aes(y=betas_std[,1]), size=3, pch=23, color='black', fill='red2') +
  theme_light() + legend_none() +
  theme(text = element_text(size=20)) 




# Spatial GP: rho_k ---------------------------------------------------------
rho_post = as.array(fit, pars = c("r"))
mcmc_combo(rho_post)
mcmc_intervals(rho_post, prob = .95, prob_outer = 1) + 
  geom_point(aes(x=c(r2,r1), y=1:2), fill=2, size=3, pch=21)
rho_ci=mcmc_intervals_data(rho_post, prob = .95, prob_outer = 1)
rho_ci$l
rho_ci$m
rho_ci$h

# Pollution profiles: H ---------------------------------------
H_post = as.array(fit, pars = c("H"))
#mcmc_trace(H_post)
#H_true = c(rev(coeff_c[,1]), rev(coeff_c[,2]), rev(coeff_c[,3]), rev(coeff_c[,4]), rev(coeff_c[,5]), rev(coeff_c[,6]))
H_true = c(coeff_c[,1], coeff_c[,2], coeff_c[,3], coeff_c[,4], coeff_c[,5], coeff_c[,6])
H_ci = mcmc_intervals_data(H_post, prob_outer = 1, prob = .95)
ggplot(H_ci, aes(x=parameter, y=m, color='black')) + 
  geom_errorbar(aes(ymin=l, ymax=h), width=.3, lwd=0.8, colour='blue4') +
  geom_point(size=3) + geom_point(size=4, pch=21, color='black', fill='deepskyblue1') + 
  geom_point(aes(y=H_true), size=1.8, pch=23, color='black', fill='red2') + 
  ylab(" ") + xlab(" ") +  theme_light() + legend_none() +
  theme(text = element_text(size=20)) +
  #scale_x_discrete(labels=names(cov_selection[-9])) 
  scale_x_discrete(labels=c(expression(h[1]^1), expression(h[2]^1), expression(h[1]^2),
                            expression(h[2]^2), expression(h[1]^3), expression(h[2]^3),
                            expression(h[1]^4), expression(h[1]^4), expression(h[1]^5),
                            expression(h[2]^5), expression(h[1]^6), expression(h[2]^6)))

hci_string = paste(round(H_ci$m,3), "(", round(H_ci$l,3),", ", round(H_ci$h,3), ")", sep = "")
h_table = data.frame(
  true_h1 = coeff_c[1,],
  est_h1 = hci_string[seq(2,12,by=2)],
  true_h2 = coeff_c[2,],
  est_h2 = hci_string[seq(1,12,by=2)]
)
xtable(h_table)



### Marginal posterior 95% CIs 
plot_list <- list() 
for (k in 1:K_star) {
  h_ci = H_ci[seq(1,C*K_star,by=K_star)+k-1, ]
  h_ci$true = H_true[seq(1,C*K_star,by=K_star)+k-1]
  if(k==1) plot_list[[k]] <- ggplot(h_ci, aes(y=parameter, x=m, color='black')) + 
    geom_errorbar(aes(xmin=l, xmax=h), width=.5, lwd=0.6, colour='blue4') +
    geom_point(size=4, pch=21, color='blue4', fill='deepskyblue') + 
    geom_point(aes(x=true ), size=2, pch=23, color='black', fill='black') + 
    ylab(" ") + xlab(" ") +  theme_light() + legend_none() +
    theme(text = element_text(size=15)) + ggtitle(paste0("Source ",k)) +
    scale_y_discrete(labels=c(paste0("pollutant ", 1:6)))
  else plot_list[[k]] <- ggplot(h_ci, aes(y=parameter, x=m, color='black')) + 
    geom_errorbar(aes(xmin=l, xmax=h), width=.3, lwd=0.7, colour='red3') +
    geom_point(size=4, pch=21, color='red3', fill='orange', alpha=.8) + 
    geom_point(aes(x=true ), size=6, pch="x", color='black', fill='black') + 
    ylab(" ") + xlab(" ") +  theme_light() + legend_none() +
    theme(text = element_text(size=15)) + ggtitle(paste0("Source ",k)) #+
    #theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
}
#grid.arrange(grobs=plot_list, widths=c(0.55, 0.45), ncol=K_star)
grid.arrange(grobs=plot_list, nrow=K_star)



#NEW H PLOT x BIOM
H_ci10 = mcmc_intervals_data(as.array(fit10, pars = c("H")), prob_outer = 1, prob = .95)
H_ci20 = mcmc_intervals_data(as.array(fit20, pars = c("H")), prob_outer = 1, prob = .95)
H_ci30 = mcmc_intervals_data(as.array(fit30, pars = c("H")), prob_outer = 1, prob = .95)
H_ci50 = mcmc_intervals_data(as.array(fit50, pars = c("H")), prob_outer = 1, prob = .95)

k = 2 #2
yy = 1:6
h_ci = H_ci10[seq(1,C*K_star,by=K_star)+k-1, ]
h_ci$true = H_true[seq(1,C*K_star,by=K_star)+k-1]
colH = hcl.colors(5,"BluGrn")[4:1]
#colH = c('grey70', 'grey60', 'grey45', 'grey30')

ggplot(h_ci, aes(y=yy+0.3, x=m, color='black')) + 
  geom_errorbar(aes(xmin=l, xmax=h), width=.175, lwd=1, colour=colH[1]) +
  geom_point(size=2.5, pch=21, color=colH[1], fill=colH[1], alpha=.8) + 
  geom_point(aes(x=true ), size=6, pch="x", color='black', fill='black') + 
  ylab(" ") + xlab(" ") +  theme_light() + legend_none() +
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) + 
  scale_y_continuous(breaks=yy,labels=c(paste0("pollutant ", 1:6))) +
  geom_hline(yintercept=yy[-1]-0.5, lwd=.5, lty=2) + ggtitle(paste0("Source ",k)) +
  # var20
  geom_errorbar(data=H_ci20[seq(1,C*K_star,by=K_star)+k-1, ], 
                aes(y=yy+0.1, xmin=l, xmax=h), width=.175, lwd=1, colour=colH[2]) +
  geom_point(data=H_ci20[seq(1,C*K_star,by=K_star)+k-1, ],
             aes(y=yy+0.1, x=m, color='black'), size=2.5, pch=21, color=colH[2], fill=colH[2], alpha=.8) +
  geom_point(aes(y=yy+0.1, x=true), size=6, pch="x", color='black', fill='black') +
  # var30
  geom_errorbar(data=H_ci30[seq(1,C*K_star,by=K_star)+k-1, ], 
                aes(y=yy-0.1, xmin=l, xmax=h), width=.175, lwd=1, colour=colH[3]) +
  geom_point(data=H_ci20[seq(1,C*K_star,by=K_star)+k-1, ],
             aes(y=yy-0.1, x=m, color='black'), size=2.5, pch=21, color=colH[3], fill=colH[3], alpha=.6) +
  geom_point(aes(y=yy-0.1, x=true), size=6, pch="x", color='black', fill='black') +
  # var50
  geom_errorbar(data=H_ci50[seq(1,C*K_star,by=K_star)+k-1, ],
                aes(y=yy-0.3, xmin=l, xmax=h), width=.175, lwd=1, colour=colH[4]) +
  geom_point(data=H_ci20[seq(1,C*K_star,by=K_star)+k-1, ],
             aes(y=yy-0.3, x=m, color='black'), size=2.5, pch=21, color=colH[4], fill=colH[4], alpha=.6) +
  geom_point(aes(y=yy-0.3, x=true), size=6, pch="x", color='black', fill='black')




### Barplot
data_bar = H_ci
species = NULL
for (c in 1:C) { species = c(species, rep(paste0("pollutant ", c), K_star)) }
data_bar$specie = species
data_bar$source = factor(rep(1:K_star, C))
data_bar$xx = c(rbind(1:6+0.2, 1:6-0.2)) 
trueH = data.frame(value=c(coeff_c[1,], coeff_c[2,]), source=factor(c(rep(1,C), rep(2,C))), specie=rep(unique(species),2))
# with colors
ggplot(trueH, aes(fill=source, y=value, x=specie)) + 
  geom_bar(position="dodge", stat="identity", width=.8, alpha=.8) +
  xlab("") + ylab("") + theme_light() +
  scale_fill_manual(values=col_light[c(1,3)]) +
  theme(text=element_text(size=20)) + #theme(legend.position="bottom") + 
  geom_errorbar(data=data_bar, aes(x=xx, y=m, ymin=l, ymax=h), lwd=.6, width=.4) +
  geom_point(data=data_bar, aes(x=xx, y=m), pch=19, size=2.5)
# B&W
ggplot(trueH, aes(fill=source, y=value, x=specie)) + 
  geom_bar(position="dodge", stat="identity", width=.8) +
  xlab("") + ylab("") + theme_light() +
  scale_fill_manual(values=c("grey70","grey85")) +
  theme(text=element_text(size=20), legend.position=" ") + 
  geom_errorbar(data=data_bar, aes(x=xx, y=m, ymin=l, ymax=h), lwd=.6, width=.4) +
  geom_point(data=data_bar, aes(x=xx, y=m), pch=19, size=2)





# Local source contributions: g_ki ---------------------------------------------
L_mean = L_median
### Local coefficients: gamma_ki
coeff_g_i = matrix(mcmc_intervals_data(as.array(fit, pars = c("coeff_g")))$m, nrow=N, ncol = K_star, byrow = T)
# standardized gamma_ki
plot(exp(coeff_g[,1])*sqrt(sum(w[1,]^2)), pch=19)
points(exp(coeff_g_i[,1])*sqrt(sum(L_mean[1,]^2)), col=2, lwd=2)
plot(exp(coeff_g[,2])*sqrt(sum(w[2,]^2)), pch=19)
points(exp(coeff_g_i[,2])*sqrt(sum(L_mean[2,]^2)), col=2, lwd=2)

plot(exp(coeff_g_i[,3]), col=2, lwd=2)

plot(max(exp(coeff_g_i[,3]))*f_ci[seq(1,T*K_star,by=K_star)+2, ]$m, type='l', ylim=c(0,0.1))
lines(max(exp(coeff_g_i[,3]))*f_ci[seq(1,T*K_star,by=K_star)+2, ]$h)
lines(max(exp(coeff_g_i[,3]))*f_ci[seq(1,T*K_star,by=K_star)+2, ]$l)
### Functional estimates: G_i
# Multiplot checking all locations for k=1,2
fit = fit10
G10_ci = mcmc_intervals_data(as.array(fit10, pars = c("g")), prob=0.95, prob_outer = 1)
G20_ci = mcmc_intervals_data(as.array(fit20, pars = c("g")), prob=0.95, prob_outer = 1)
G30_ci = mcmc_intervals_data(as.array(fit30, pars = c("g")), prob=0.95, prob_outer = 1)
G50_ci = mcmc_intervals_data(as.array(fit50, pars = c("g")), prob=0.95, prob_outer = 1)

# plotting ALL NxK functional profiles
G_ci = G20_ci
k = 2 #k=1
g_k = G_ci[seq(1,N*T*K_star,by=K_star)+k-1, ]
plot_list <- list() 
for(i in 1:N){
  g_ki = g_k[seq(1,N*T,by=N)+i-1, ]
  g_ki$t = t_basis
  g_ki$g_true = g[i,1, ] # g[i,2, ]
  plot_list[[i]] <- ggplot(data=g_ki, aes(x=t)) +
    geom_line(aes(y=m), color = 'blue', lwd=.8)  +
    geom_ribbon(data=g_ki, aes(ymin = l, ymax = h), fill='deepskyblue', alpha = 0.3) +
    theme_light() + labs(x = NULL, y = NULL) + ggtitle(paste0("location ",i)) +
    geom_line(aes(y=g_true), color = 'red', lwd=.6)
}
grid.arrange(grobs=plot_list[1:16], ncol=4)
grid.arrange(grobs=plot_list[17:32], ncol=4)

matplot(t(apply(g[ ,2, ], 1, range)), pch=c(1,1))
matplot(t(apply(g[ ,1, ], 1, range)), add=T, pch=c(2,2))
abline(h=30)

### Sample locations: i = 11,9,18
i_list = c(5,7,17)
# B&W
plot_list <- list() 
for (k in 1:K_star) {
g_k = G_ci[seq(1,N*T*K_star,by=K_star)+k-1, ]
for(i in i_list){
  g_ki = g_k[seq(1,N*T,by=N)+i-1, ]
  g_ki$t = t_basis
  g_ki$g_true = g[i,as.numeric(k==2)+1, ]
  if(i==i_list[1])
    plot_list[[k]] <- ggplot(data=g_ki, aes(x=t)) +
    theme_light() + labs(x = NULL, y = NULL) + 
    geom_line(aes(y=m), color = 'grey50', lwd=.8)  +
    geom_ribbon(aes(ymin = l, ymax = h), fill='grey70', alpha = 0.4) +
    geom_line(aes(y=g_true), color = 'black', lwd=.7, lty=2) + ggtitle(paste0("Source ",k)) +
    theme(text=element_text(size=20), plot.title = element_text(hjust = 0.5)) +
    coord_cartesian(xlim=c(-1,366), expand=F)
  else plot_list[[k]] <- plot_list[[k]] +  
    geom_line(data=g_ki, aes(y=m), color = 'grey45', lwd=.8)  +
    geom_ribbon(data=g_ki, aes(ymin = l, ymax = h), fill='grey70', alpha = 0.4) +
    geom_line(data=g_ki, aes(y=g_true), color = 'black', lwd=.7, lty=2)
}}
grid.arrange(grobs=plot_list, ncol=2)

# with colors
G_ci = G10_ci
col_list = c("darkolivegreen3", "goldenrod1", "deeppink2")
i_list = c(5,7,17)
plot_list <- list() 
for (k in 1:K_star) {
  g_k = G_ci[seq(1,N*T*K_star,by=K_star)+k-1, ]
  for(i in i_list){
    g_ki = g_k[seq(1,N*T,by=N)+i-1, ]
    g_ki$t = t_basis
    g_ki$g_true = g[i,as.numeric(k==2) +1, ]
    if(i==i_list[1])
      plot_list[[k]] <- ggplot(data=g_ki, aes(x=t)) +
      theme_light() + labs(x = NULL, y = NULL) + 
      geom_line(aes(y=m), color = col_list[i_list==i], lwd=.8)  +
      geom_ribbon(aes(ymin = l, ymax = h), fill=col_list[i_list==i], alpha = 0.375) +
      geom_line(aes(y=g_true), color = 'black', lwd=1, lty=2) + ggtitle(paste0("Source ",k)) +
      theme(text=element_text(size=20), plot.title = element_text(hjust = 0.5)) +
      coord_cartesian(xlim=c(-1,366), expand=F, ylim=c(ifelse(k==1,0.2,0), ifelse(k==1,12.6,7.8)))
    else plot_list[[k]] <- plot_list[[k]] +  
      geom_line(data=g_ki, aes(y=m), color = col_list[which(i_list==i)], lwd=.8)  +
      geom_ribbon(data=g_ki, aes(ymin = l, ymax = h), fill=col_list[i_list==i], alpha = 0.375) +
      geom_line(data=g_ki, aes(y=g_true), color = 'black', lwd=1, lty=2)
  }}
grid.arrange(grobs=plot_list, ncol=2)
grid.arrange(grobs=plot_list, nrow=2)



# Pure noise variance: sigma_c --------------------------------------------
sigma_post = as.array(fit, pars = c("sigma"))
mcmc_combo(sigma_post)
mcmc_intervals(sigma_post, prob = .95, prob_outer = 1) + 
  geom_point(aes(x=as.numeric(sigma_c), y=6:1), fill="red", size=2, pch=21)




# Posterior predictive distribution: y_ic ---------------------------------
yci = mcmc_intervals_data(as.array(fit, pars = c("Y_pred")), prob=0.95, prob_outer = 1)
i = 7
y_i = yci[seq(1,N*T*C,by=N)+i-1, ]
plot_list <- list() 
for(c in 1:C){
  y_ic = y_i[seq(1,C*T,by=C)+c-1, ]
  y_ic$t = t_basis
  y_ic$yobs = Y_obs[i,c, ]
  plot_list[[c]] <- ggplot(data=y_ic, aes(x=t)) +
    geom_line(aes(y=m), color = 'blue', lwd=.8)  +
    geom_ribbon(data=y_ic, aes(ymin = l, ymax = h), fill='deepskyblue', alpha = 0.3) +
    theme_light() + labs(x = NULL, y = NULL) + ggtitle(paste0("location ",i)) +
    geom_point(aes(y=yobs), color = 'red', size=1)
}
grid.arrange(grobs=plot_list, ncol=3)




# Additional graphics: g_ki ----------------------------------------------

# B&W 
i = 11
k = 1
g_k = G_ci[seq(1,N*T*K_star,by=K_star)+k-1, ]
g_ki = g_k[seq(1,N*T,by=N)+i-1, ]
g_ki$t = t_basis
g_ki$g_true = g[i,2, ]
plot_1 <- ggplot(data=g_ki, aes(x=t)) +
  geom_line(aes(y=m), color = 'grey50', lwd=.8)  +
  geom_ribbon(data=g_ki, aes(ymin = l, ymax = h), fill='grey70', alpha = 0.4) +
  theme_light() + labs(x = NULL, y = NULL) + ggtitle(paste0("Source ",1)) +
  geom_line(aes(y=g_true), color = 'black', lwd=.7, lty=2) +
  theme(text=element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(xlim=c(-1,366), ylim=c(-0.19,0.9), expand=F)
k = 2
g_k = G_ci[seq(1,N*T*K_star,by=K_star)+k-1, ]
g_ki = g_k[seq(1,N*T,by=N)+i-1, ]
g_ki$t = t_basis
g_ki$g_true = g[i,1, ]
plot_2 <- ggplot(data=g_ki, aes(x=t)) +
  geom_line(aes(y=m), color = 'grey45', lwd=.8)  +
  geom_ribbon(data=g_ki, aes(ymin = l, ymax = h), fill='grey70', alpha = 0.4) +
  theme_light() + labs(x = NULL, y = NULL) + ggtitle(paste0("Source ",2)) +
  geom_line(aes(y=g_true), color = 'black', lwd=.7, lty=2) +
  theme(text=element_text(size=20), plot.title = element_text(hjust = 0.5))+
  coord_cartesian(xlim=c(-1,366), ylim=c(-0.05,0.9), expand=F)

i = 9 
k = 1
g_k = G_ci[seq(1,N*T*K_star,by=K_star)+k-1, ]
g_ki = g_k[seq(1,N*T,by=N)+i-1, ]
g_ki$t = t_basis
g_ki$g_true = g[i,2, ]
plot_1 = plot_1 + 
  geom_line(data=g_ki, aes(y=m), color = 'grey45', lwd=.8)  +
  geom_ribbon(data=g_ki, aes(ymin = l, ymax = h), fill='grey70', alpha = 0.4) +
  geom_line(data=g_ki, aes(y=g_true), color = 'black', lwd=.7, lty=2)
k = 2
g_k = G_ci[seq(1,N*T*K_star,by=K_star)+k-1, ]
g_ki = g_k[seq(1,N*T,by=N)+i-1, ]
g_ki$t = t_basis
g_ki$g_true = g[i,1, ]
plot_2 = plot_2 + 
  geom_line(data=g_ki, aes(y=m), color = 'grey45', lwd=.8)  +
  geom_ribbon(data=g_ki, aes(ymin = l, ymax = h), fill='grey70', alpha = 0.4) +
  geom_line(data=g_ki, aes(y=g_true), color = 'black', lwd=.7, lty=2)

i = 18
k = 1
g_k = G_ci[seq(1,N*T*K_star,by=K_star)+k-1, ]
g_ki = g_k[seq(1,N*T,by=N)+i-1, ]
g_ki$t = t_basis
g_ki$g_true = g[i,2, ]
plot_1 = plot_1 + 
  geom_line(data=g_ki, aes(y=m), color = 'grey45', lwd=.8)  +
  geom_ribbon(data=g_ki, aes(ymin = l, ymax = h), fill='grey70', alpha = 0.4) +
  geom_line(data=g_ki, aes(y=g_true), color = 'black', lwd=.7, lty=2)
k = 2
g_k = G_ci[seq(1,N*T*K_star,by=K_star)+k-1, ]
g_ki = g_k[seq(1,N*T,by=N)+i-1, ]
g_ki$t = t_basis
g_ki$g_true = g[i,1, ]
plot_2 = plot_2 + 
  geom_line(data=g_ki, aes(y=m), color = 'grey45', lwd=.8)  +
  geom_ribbon(data=g_ki, aes(ymin = l, ymax = h), fill='grey70', alpha = 0.4) +
  geom_line(data=g_ki, aes(y=g_true), color = 'black', lwd=.7, lty=2)

grid.arrange(plot_1,plot_2, ncol=2)


# with colors 
i = 11
k = 1
g_k = G_ci[seq(1,N*T*K_star,by=K_star)+k-1, ]
g_ki = g_k[seq(1,N*T,by=N)+i-1, ]
g_ki$t = t_basis
g_ki$g_true = g[i,2, ]
plot_1 <- ggplot(data=g_ki, aes(x=t)) +
  geom_line(aes(y=m), color = 'deepskyblue3', lwd=1.1)  +
  geom_ribbon(data=g_ki, aes(ymin = l, ymax = h), fill='deepskyblue', alpha = 0.3) +
  theme_light() + labs(x = NULL, y = NULL) + ggtitle(paste0("Source ",1)) +
  geom_line(aes(y=g_true), color = 'black', lwd=1, lty=2) +
  theme(text=element_text(size=20), plot.title = element_text(hjust = 0.5))+
  coord_cartesian(xlim=c(-1,366), ylim=c(-0.19,0.9), expand=F)
k = 2
g_k = G_ci[seq(1,N*T*K_star,by=K_star)+k-1, ]
g_ki = g_k[seq(1,N*T,by=N)+i-1, ]
g_ki$t = t_basis
g_ki$g_true = g[i,1, ]
plot_2 <- ggplot(data=g_ki, aes(x=t)) +
  geom_line(aes(y=m), color = 'deepskyblue3', lwd=1.1)  +
  geom_ribbon(data=g_ki, aes(ymin = l, ymax = h), fill='deepskyblue', alpha = 0.3) +
  theme_light() + labs(x = NULL, y = NULL) + ggtitle(paste0("Source ",2)) +
  geom_line(aes(y=g_true), color = 'black', lwd=1, lty=2) +
  theme(text=element_text(size=20), plot.title = element_text(hjust = 0.5))+
  coord_cartesian(xlim=c(-1,366), ylim=c(-0.05,0.9), expand=F)

i = 9
k = 1
g_k = G_ci[seq(1,N*T*K_star,by=K_star)+k-1, ]
g_ki = g_k[seq(1,N*T,by=N)+i-1, ]
g_ki$t = t_basis
g_ki$g_true = g[i,2, ]
plot_1_new = plot_1 + 
  geom_line(data=g_ki, aes(y=m), color = 'deeppink3', lwd=1.1)  +
  geom_ribbon(data=g_ki, aes(ymin = l, ymax = h), fill='deeppink', alpha = 0.3) +
  geom_line(data=g_ki, aes(y=g_true), color = 'black', lwd=1, lty=2)
k = 2
g_k = G_ci[seq(1,N*T*K_star,by=K_star)+k-1, ]
g_ki = g_k[seq(1,N*T,by=N)+i-1, ]
g_ki$t = t_basis
g_ki$g_true = g[i,1, ]
plot_2_new = plot_2 + 
  geom_line(data=g_ki, aes(y=m), color = 'deeppink3', lwd=1.1)  +
  geom_ribbon(data=g_ki, aes(ymin = l, ymax = h), fill='deeppink', alpha = 0.3) +
  geom_line(data=g_ki, aes(y=g_true), color = 'black', lwd=1, lty=2)

i = 18 
k = 1
g_k = G_ci[seq(1,N*T*K_star,by=K_star)+k-1, ]
g_ki = g_k[seq(1,N*T,by=N)+i-1, ]
g_ki$t = t_basis
g_ki$g_true = g[i,2, ]
plot_1_new = plot_1_new + 
  geom_line(data=g_ki, aes(y=m), color = 'darkorange', lwd=1.1)  +
  geom_ribbon(data=g_ki, aes(ymin = l, ymax = h), fill='orange', alpha = 0.3) +
  geom_line(data=g_ki, aes(y=g_true), color = 'black', lwd=1, lty=2)
k = 2
g_k = G_ci[seq(1,N*T*K_star,by=K_star)+k-1, ]
g_ki = g_k[seq(1,N*T,by=N)+i-1, ]
g_ki$t = t_basis
g_ki$g_true = g[i,1, ]
plot_2_new = plot_2_new + 
  geom_line(data=g_ki, aes(y=m), color = 'darkorange', lwd=1.1)  +
  geom_ribbon(data=g_ki, aes(ymin = l, ymax = h), fill='orange', alpha = 0.3) +
  geom_line(data=g_ki, aes(y=g_true), color = 'black', lwd=1, lty=2)

grid.arrange(plot_1_new,plot_2_new, ncol=2)


# additional plots for simulated data -------------------------------------
max_per_c = c(40,35,35,20,35,35)
library(reshape2)
plot_list <- list() 
load("sim_data10_new.RData")
for(c in 1:6){
  df_y = as.data.frame(t(exp(Y_obs[ ,c, ])))
  df_med = data.frame(values=apply(exp(Y_obs[ ,c, ]),2,median,na.rm=T), variable=1)
  df_y$id = df_med$id = t_basis
  plot_data = melt(df_y,id.var="id")
  plot_list[[c]] <- ggplot(plot_data, aes(x=id,y=value,group=variable,colour=variable)) +
    geom_line(aes(lty=factor(as.numeric(variable)%%4)), lwd=.65) + ggtitle(paste0("pollutant ",c))+ 
    scale_colour_viridis_d(begin=0,end=1) + theme_light() + labs(x = NULL, y = NULL) +
    geom_line(data = df_med, aes(x=id, y=values, group=variable), lty=1, col=1, lwd=1.5) +
    theme(plot.title=element_text(hjust=.5), text=element_text(size=20), legend.position="none") +
    #coord_cartesian(xlim=c(-1,366), ylim=c(-0.9,4.5), expand=F)
    coord_cartesian(xlim=c(-1,366), expand=F, ylim=c(0, max_per_c[c]))
}
load("sim_data20_new.RData")
for(c in 1:6){
  df_y = as.data.frame(t(exp(Y_obs[ ,c, ])))
  df_med = data.frame(values=apply(exp(Y_obs[ ,c, ]),2,median,na.rm=T), variable=1)
  df_y$id = df_med$id = t_basis
  plot_data = melt(df_y,id.var="id")
  plot_list[[c+6]] <- ggplot(plot_data, aes(x=id,y=value,group=variable,colour=variable)) +
    geom_line(aes(lty=factor(as.numeric(variable)%%4)), lwd=.65) + ggtitle(paste0("pollutant ",c))+ 
    scale_colour_viridis_d(begin=0,end=1) + theme_light() + labs(x = NULL, y = NULL) +
    geom_line(data = df_med, aes(x=id, y=values, group=variable), lty=1, col=1, lwd=1.5) +
    theme(plot.title=element_text(hjust=.5), text=element_text(size=20), legend.position="none") +
    #coord_cartesian(xlim=c(-1,366), ylim=c(-0.9,4.5), expand=F)
    coord_cartesian(xlim=c(-1,366), expand=F, ylim=c(0, max_per_c[c]))
}
load("sim_data30_new.RData")
for(c in 1:6){
  df_y = as.data.frame(t(exp(Y_obs[ ,c, ])))
  df_med = data.frame(values=apply(exp(Y_obs[ ,c, ]),2,median,na.rm=T), variable=1)
  df_y$id = df_med$id = t_basis
  plot_data = melt(df_y,id.var="id")
  plot_list[[c+12]] <- ggplot(plot_data, aes(x=id,y=value,group=variable,colour=variable)) +
    geom_line(aes(lty=factor(as.numeric(variable)%%4)), lwd=.65) + ggtitle(paste0("pollutant ",c))+ 
    scale_colour_viridis_d(begin=0,end=1) + theme_light() + labs(x = NULL, y = NULL) +
    geom_line(data = df_med, aes(x=id, y=values, group=variable), lty=1, col=1, lwd=1.5) +
    theme(plot.title=element_text(hjust=.5), text=element_text(size=20), legend.position="none") +
    #coord_cartesian(xlim=c(-1,366), ylim=c(-0.9,4.5), expand=F)
    coord_cartesian(xlim=c(-1,366), expand=F, ylim=c(0, max_per_c[c]))
}
load("sim_data50_new.RData")
for(c in 1:6){
  df_y = as.data.frame(t(exp(Y_obs[ ,c, ])))
  df_med = data.frame(values=apply(exp(Y_obs[ ,c, ]),2,median,na.rm=T), variable=1)
  df_y$id = df_med$id = t_basis
  plot_data = melt(df_y,id.var="id")
  plot_list[[c+18]] <- ggplot(plot_data, aes(x=id,y=value,group=variable,colour=variable)) +
    geom_line(aes(lty=factor(as.numeric(variable)%%4)), lwd=.65) + ggtitle(paste0("pollutant ",c))+ 
    scale_colour_viridis_d(begin=0,end=1) + theme_light() + labs(x = NULL, y = NULL) +
    geom_line(data = df_med, aes(x=id, y=values, group=variable), lty=1, col=1, lwd=1.5) +
    theme(plot.title=element_text(hjust=.5), text=element_text(size=20), legend.position="none") +
    #coord_cartesian(xlim=c(-1,366), ylim=c(-0.9,4.5), expand=F)
    coord_cartesian(xlim=c(-1,366), expand=F, ylim=c(0, max_per_c[c]))
}
new_order = c(1:3,1:3+6,1:3+12,1:3+18)
new_order = c(4:6,4:6+6,4:6+12,4:6+18)
grid.arrange(grobs=plot_list[new_order], nrow=4, ncol=3)



# plot g_ki f_k
df_gf = as.data.frame(t(g[,1,]))
df_med = data.frame(values=f[,1], variable=1)
df_gf$id = df_med$id= t_basis
plot_data = melt(df_gf,id.var="id")
plot_g1=ggplot(plot_data, aes(x=id,y=value,group=variable,colour=variable)) +
  geom_line(aes(lty=factor(as.numeric(variable)%%4)), lwd=.6) + ggtitle("Source 1")+ 
  scale_colour_viridis_d(begin=0,end=.95) + theme_light() + labs(x = NULL, y = NULL) +
  geom_line(data = df_med, aes(x=id, y=values, group=variable), lty=1, col=1, lwd=3) +
  theme(plot.title=element_text(hjust=.5), text=element_text(size=20), legend.position="none") +
  #coord_cartesian(xlim=c(-1,366), ylim=c(-0.9,4.5), expand=F)
  coord_cartesian(xlim=c(-1,366), expand=F)
df_gf = as.data.frame(t(g[,2,]))
df_med = data.frame(values=f[,2], variable=1)
df_gf$id = df_med$id= t_basis
plot_data = melt(df_gf,id.var="id")
plot_g2=ggplot(plot_data, aes(x=id,y=value,group=variable,colour=variable)) +
  geom_line(aes(lty=factor(as.numeric(variable)%%4)), lwd=.6) + ggtitle("Source 2")+ 
  scale_colour_viridis_d(begin=0,end=.95) + theme_light() + labs(x = NULL, y = NULL) +
  geom_line(data = df_med, aes(x=id, y=values, group=variable), lty=1, col=1, lwd=3) +
  theme(plot.title=element_text(hjust=.5), text=element_text(size=20), legend.position="none") +
  #coord_cartesian(xlim=c(-1,366), ylim=c(-0.9,4.5), expand=F)
  coord_cartesian(xlim=c(-1,366), expand=F)
grid.arrange(grobs=list(plot_g1,plot_g2), ncol=2)



