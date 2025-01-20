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
library(dplyr)
library(xtable)

# Loading data and MCMC ---------------------------------------------------
load("final_PM25_data.RData")
fit = read_rds("PM25_fit_dirichlet.RDS")

col_light = c("deepskyblue", "deeppink", "orange") #"yellowgreen"
col_dark = c("blue", "deeppink3", "darkorange3") # "forestgreen",
col_points = c("blue3", "deeppink3", "darkorange4") #"darkgreen", # nolint

# Shrinkage check: L and K* ----------------------------------------------------
K_post = as.array(fit, pars = c("count"))
K_star = median(K_post) # nolint: object_name_linter.
mcmc_trace(K_post)

L = mcmc_intervals_data(as.matrix(fit, pars = c("L")))
L_mean = matrix(L$m, nrow = K_star, ncol = numbasis)

matplot(t(L_mean))
matplot(abs(t(L_mean)))
abline(h=c(0,0.05), lty=2)
mcmc_trace(as.array(fit, pars = c("L")))


# Regression coefficients: beta_k -----------------------------------------
betas_post = as.array(fit, pars = c("beta"))
mcmc_trace(betas_post)
mcmc_intervals(betas_post, prob=0.95, prob_outer = 1)

betas_ci = mcmc_intervals_data(betas_post, prob=0.95, prob_outer = 1)
beta1 = betas_ci[seq(1, 4*K_star,by=K_star), c(1,6:8)]
beta2 = betas_ci[seq(1, 4*K_star,by=K_star)+1, c(1,6:8)]
beta3 = betas_ci[seq(1, 4*K_star,by=K_star)+2, c(1,6:8)]

# table of values to be included in the paper
beta_table = data.frame(
  suburb = rbind(sprintf('%s (%s)', toString(round(beta1[1,3],2)), toString(round(beta1[1,c(2,4)], 2))),
                 sprintf('%s (%s)', toString(round(beta2[1,3],2)), toString(round(beta2[1,c(2,4)], 2))),
                 sprintf('%s (%s)', toString(round(beta3[1,3],2)), toString(round(beta3[1,c(2,4)], 2)))),
  urban = rbind( sprintf('%s (%s)', toString(round(beta1[2,3],2)), toString(round(beta1[2,c(2,4)], 2))),
                 sprintf('%s (%s)', toString(round(beta2[2,3],2)), toString(round(beta2[2,c(2,4)], 2))),
                 sprintf('%s (%s)', toString(round(beta3[2,3],2)), toString(round(beta3[2,c(2,4)], 2)))),
  altitude = rbind(sprintf('%s (%s)', toString(round(beta1[3,3],2)), toString(round(beta1[3,c(2,4)], 2))),
                   sprintf('%s (%s)', toString(round(beta2[3,3],2)), toString(round(beta2[3,c(2,4)], 2))),
                   sprintf('%s (%s)', toString(round(beta3[3,3],2)), toString(round(beta3[3,c(2,4)], 2)))))
xtable(beta_table)


# Spatial GP: rho_k ---------------------------------------------------------
rho_post = as.array(fit, pars = c("r"))
mcmc_combo(rho_post)
rho_ci = mcmc_intervals_data(rho_post, prob = .95, prob_outer = 1)
rho_ci



# Posterior predictive: Y_ic -------------
yCI = mcmc_intervals_data(as.array(fit, pars=c("Y_pred")), prob=.95, prob_outer=1)
i_list = c(17, 19, 31)
i_names = staz$`Local Site Name`[i_list]
i_names
pollutants = names(data)[3:8]

plot_list <- list() 
for (i in i_list) {
  y_i = yCI[seq(1,N*T*C,by=N)+i-1, ]
  for (c in 1:C){
    y_ci = y_i[seq(1,T*C,by=C)+c-1, ]
    y_ci$t = t_basis
    y_ci$ytrue = Y_obs[i,c, ]
    if(i==i_list[1])
      plot_list[[c]] <- ggplot(data=y_ci, aes(x=t)) +
      theme_light() + geom_line(aes(y=m), color = col_dark[1], lwd=.8)  +
      geom_ribbon(data=y_ci, aes(ymin=l, ymax=h), fill=col_light[1], alpha=.3) +
      geom_point(aes(x=t,y=ytrue), color = col_points[1], size=.8)  +
      labs(x = NULL, y = NULL) + ggtitle(pollutants[c])
    else plot_list[[c]] = plot_list[[c]] + geom_line(data=y_ci, aes(y=m), color = col_dark[i_list==i], lwd=.8)  +
      geom_ribbon(data=y_ci, aes(ymin = l, ymax = h), fill=col_light[i_list==i], alpha = 0.3) +
      geom_point(data=y_ci, aes(x=t,y=ytrue), color = col_points[i_list==i], size=.8)
  }}
grid.arrange(grobs=plot_list, ncol=2)



# Local functional contributions: g_ki ---------------
gpred = as.array(fit, pars = c("g"))
g_ci = mcmc_intervals_data(gpred, prob=0.95, prob_outer = 1)
i_list = c(17,19,31) # 3 locations
ltypes = c(1:2,4)    # 3 line types

# plot with colors
plot_list <- list() 
for (k in 1:K_star) {
  g_k = g_ci[seq(1,N*T*K_star,by=K_star)+k-1, ]
  for(i in i_list){
    g_ki = g_k[seq(1,N*T,by=N)+i-1, ]
    g_ki$t = sort(unique(data$Date))
    if(i==i_list[1])
      plot_list[[k]] <- ggplot(data=g_ki, aes(x=t)) +
      theme_light() + labs(x = NULL, y = NULL) + 
      geom_line(aes(y=m), color = col_dark[1], lwd=.9, lty=ltypes[i_list==i])  +
      geom_ribbon(aes(ymin = l, ymax = h), fill=col_light[1], alpha = 0.3,
                  color= col_light[i_list==i], lty=ltypes[i_list==i], linewidth=.3) +
      ggtitle(paste0("k = ", k)) + theme(text=element_text(size=20), plot.title = element_text(hjust = 0.5))
    else plot_list[[k]] <- plot_list[[k]] +  
      geom_line(data=g_ki, aes(y=m), color = col_dark[i_list==i], lwd=.9, lty=ltypes[i_list==i])  +
      geom_ribbon(data=g_ki, aes(ymin = l, ymax = h), fill=col_light[i_list==i], alpha = 0.3,
                  color= col_light[i_list==i], lty=ltypes[i_list==i], linewidth=.3)
  }}
#grid.arrange(grobs=plot_list, ncol=2)
grid.arrange(grobs=plot_list, ncol=3)

# B&W plot for the paper
plot_list <- list() 
for (k in 1:K_star) {
  g_k = g_ci[seq(1,N*T*K_star,by=K_star)+k-1, ]
  for(i in i_list){
    g_ki = g_k[seq(1,N*T,by=N)+i-1, ]
    g_ki$t = sort(unique(data$Date))
    if(i==i_list[1])
      plot_list[[k]] <- ggplot(data=g_ki, aes(x=t)) +
      theme_light() + labs(x = NULL, y = NULL) + 
      geom_ribbon(aes(ymin = l, ymax = h), fill="grey70", alpha = 0.3,
                  color= "grey50", lty=ltypes[i_list==i], linewidth=.35) +
      geom_line(aes(y=m), color = "black", lwd=.85, lty=ltypes[i_list==i]) +
      #coord_cartesian(xlim=c(min(g_ki$t)-1, max(g_ki$t)+1), expand=F) +
      ggtitle(paste0("k = ", k)) + theme(text=element_text(size=15), plot.title = element_text(hjust = 0.5))
    else plot_list[[k]] <- plot_list[[k]] +  
      geom_ribbon(data=g_ki, aes(ymin = l, ymax = h), fill="grey70", alpha = 0.3,
                  color= "grey50", lty=ltypes[i_list==i], linewidth=.35) +
      geom_line(data=g_ki, aes(y=m), color = "black", lwd=.85, lty=ltypes[i_list==i])
  }}
#grid.arrange(grobs=plot_list, ncol=2)
grid.arrange(grobs=plot_list, ncol=3)




# PLOT per gli H -----------------
CI_h = mcmc_intervals_data(as.array(fit, pars = c("H")), prob=.95, prob_outer = 1)
pollutants = names(data)[3:8]
pollutants[5] = "Nitrate"

plot_list <- list() 
for (k in 1:K_star) {
  h_ci = CI_h[seq(1,C*K_star,by=K_star)+k-1, ]
  if(k==1) plot_list[[k]] <- ggplot(h_ci, aes(y=parameter, x=m, color='black')) + 
      geom_errorbar(aes(xmin=l, xmax=h), width=.5, lwd=0.6, colour=col_dark[k]) +
      geom_point(size=4, pch=21, color=col_dark[k], fill=col_light[k]) + 
      ylab(" ") + xlab(" ") +  theme_light() + legend_none() +
      theme(text = element_text(size=15)) + ggtitle(paste0("Source ",k)) +
      scale_y_discrete(labels=pollutants)
  else plot_list[[k]] <- ggplot(h_ci, aes(y=parameter, x=m, color='black')) + 
      geom_errorbar(aes(xmin=l, xmax=h), width=.5, lwd=0.6, colour=col_dark[k]) +
      geom_point(size=4, pch=21, color=col_dark[k], fill=col_light[k]) + 
      ylab(" ") + xlab(" ") +  theme_light() + legend_none() +
      theme(text = element_text(size=15)) + ggtitle(paste0("Source ",k)) +
      theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
}
grid.arrange(grobs=plot_list, ncol=K_star)


yy = 1:6
plot_tot = plot_list[[1]]
for (k in 2:K_star) {
  h_ci = CI_h[seq(1,C*K_star,by=K_star)+k-1, ]
  plot_tot = plot_tot  +
    geom_errorbar(data=h_ci, aes(y=yy, xmin=l, xmax=h), width=.3, lwd=0.8, colour=col_dark[k]) +
    geom_point(data=h_ci, aes(y=yy, x=m), size=4, pch=21, color=col_dark[k], fill=col_light[k]) + 
    ylab(" ") + xlab(" ") +  theme_light() + legend_none() +
    theme(text = element_text(size=15)) + ggtitle("")
  
}
plot_tot


# prova barplot
data_bar = CI_h
species = NULL
for (c in 1:C) {
  species = c(species, rep(pollutants[c], K_star))
}
data_bar$specie = species
data_bar$source = factor(rep(1:K_star, C))

ggplot(data_bar, aes(fill=source, y=specie, x=m)) + 
  geom_bar(position="dodge", stat="identity", width=.7) +
  xlab("") + ylab("") + theme_light() +
#  scale_fill_manual(values=col_light) +
  theme(text=element_text(size=20))


#final plot H
data_bar |>
  #mutate(year = factor(year), type = factor(type, c("inter", "VIIT","HIIT"))) |> 
  group_by(specie, source) |> 
  summarise(n = m) |>
  mutate(pct = prop.table(n)) |>
  ggplot(aes(pct, specie, fill = source)) + geom_col(alpha=0.8) + 
  geom_text(aes(label = scales::percent(pct, accuracy = 1)), 
            position = position_stack(vjust = .5), size=7) +
  xlab("") + ylab("") + theme_light() +
  scale_fill_manual(values=c("mediumpurple1", "lightskyblue1", "darkolivegreen2", "gold")) +
  theme(text=element_text(size=25)) + coord_flip() +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())


# B&W plot for the paper
data_bar |>
  #mutate(year = factor(year), type = factor(type, c("inter", "VIIT","HIIT"))) |> 
  group_by(specie, source) |> 
  summarise(n = m) |>
  mutate(pct = prop.table(n)) |>
  ggplot(aes(pct, specie, fill = source)) + geom_col(alpha=0.8) + 
  geom_text(aes(label = scales::percent(pct, accuracy = 1)), 
            position = position_stack(vjust = .5), size=7) +
  xlab("") + ylab("") + theme_light() +
  scale_fill_manual(values=c("grey50", "grey75", "grey90")) +
  theme(text=element_text(size=25)) + coord_flip() +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())


# Coverage ----------------------------------------------------------------

yCI =  mcmc_intervals_data(as.array(fit, pars = c("Y_pred")), prob=0.95, prob_outer = 1)
cvg = matrix(rep(0,N*C), nrow=N, ncol= C)

for (i in 1:N) {
  y_i = yCI[seq(1,N*T*C,by=N)+i-1, ]
  for(c in 1:C){
    y_ci = y_i[seq(1,T*C,by=C)+c-1, ]
    y_ci$t = t_basis
    y_ci$ytrue = Y_obs[i,c, ]
    sum_temp = sum(y_ci$ytrue>y_ci$l & y_ci$ytrue<y_ci$h, na.rm = T)
    cvg[i,c] = sum_temp/sum(!is.na(y_ci$ytrue))
  }
}
boxplot(cvg) #ok
species = NULL
for (c in 1:C) {
  species = c(species, rep(pollutants[c], N))
}

# B&W boxplot per paper
box_gg = data.frame(value=as.vector(cvg), pollutants=species)
ggplot(data=box_gg, aes(x=pollutants, y=value, fill=pollutants)) +
  stat_boxplot(geom='errorbar', linetype=1, width=0.5, linewidth=.7)+
  geom_boxplot() + theme_light() + scale_fill_ordinal() + legend_none() +
  theme(text=element_text(size=25)) + xlab("") + ylab("") 

# B&W boxplot per paper
box_gg = data.frame(value=as.vector(cvg), pollutants=species)
ggplot(data=box_gg, aes(x=pollutants, y=value, fill=pollutants)) +
  stat_boxplot(geom='errorbar', linetype=1, width=0.5, linewidth=.7)+
  geom_boxplot() + theme_light() + scale_fill_grey(start = 0.8, end = 0.8) + 
  legend_none() + theme(text=element_text(size=25)) + xlab("") + ylab("")




# Posterior predictive: Y_ic on the exp scale
i = 2
#17ok, 19ok, 31ok, 12ok
y_i = yCI[seq(1,N*T*C,by=N)+i-1, ]

plot_list <- list() 
for(c in 1:C){
  y_ci = y_i[seq(1,T*C,by=C)+c-1, ]
  y_ci$t = sort(unique(data$Date))
  y_ci$ytrue = Y_obs[i,c, ]
  plot_list[[c]] <- ggplot(data=y_ci, aes(x=t)) +
    theme_light() + geom_line(aes(y=exp(m)), color = 'purple4', lwd=1.5)  +
    geom_ribbon(data=y_ci, aes(ymin = exp(l), ymax = exp(h)), fill='orchid2', alpha = 0.4) +
    geom_point(aes(x=t,y=exp(ytrue)), color = 'black', size=1.5)  +
    labs(x = NULL, y = NULL) + ggtitle(pollutants[c]) + 
    theme(text=element_text(size=18), plot.title = element_text(hjust = 0.5))
}
grid.arrange(grobs=plot_list, ncol=2)



# Posterior predictive: Y_ic on the exp scale -- BW plot
i = 12
#0k: 9,12,14,19,22,23,25,31,32
#best: 12, 25
y_i = yCI[seq(1,N*T*C,by=N)+i-1, ]

plot_list <- list() 
for(c in 1:C){
  y_ci = y_i[seq(1,T*C,by=C)+c-1, ]
  y_ci$t = sort(unique(data$Date))
  y_ci$ytrue = Y_obs[i,c, ]
  plot_list[[c]] <- ggplot(data=y_ci, aes(x=t)) +
    theme_light() + geom_line(aes(y=exp(m)), color = 'grey35', lwd=1.5)  +
    geom_ribbon(data=y_ci, aes(ymin = exp(l), ymax = exp(h)), fill='grey75', alpha = 0.4) +
    geom_point(aes(x=t,y=exp(ytrue)), color = 'black', size=1.5)  +
    labs(x = NULL, y = NULL) + ggtitle(pollutants[c]) + 
    theme(text=element_text(size=18), plot.title = element_text(hjust = 0.5))
}
grid.arrange(grobs=plot_list, ncol=3)

unique(data$SiteName)[i]
