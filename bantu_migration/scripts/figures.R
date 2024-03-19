# Load Libraries and spatial data ----
library(here)
library(dplyr)
library(stringr)
library(truncnorm)
library(cascsim)
library(corrplot)
library(ggplot2)
library(ggridges)
library(rnaturalearth)
library(nimbleCarbon)
library(rcarbon)
library(sf)
library(viridis)
library(cowplot)
library(wesanderson)
library(latex2exp)
library(gridExtra)
library(grid)
library(gridBase)
library(diagram)
library(quantreg)
library(coda)
library(graphics)
library(ggthemes)


source(here('src','orderPPlot.R'))
source(here('src','gpqrSim.R'))

`%!in%` <- Negate(`%in%`)

#===============================================================================
## Load Data

# Load Observed Data
load(here('data','c14.RData'))

# Load quantile regression results
load(here('output','quantreg_res.RData'))

#===============================================================================
# Generate Spatial Window for Analyses: Sub-Saharan Africa ----
#Sampling window ---- Background Map
win <- ne_countries(continent = "Africa", scale = 10, returnclass = "sf") %>%
  filter_all(., any_vars(str_detect(., "Sub-Saharan"))) %>% 
  filter(name_en %in% subSahara_countries) %>% 
  filter(name_en != "Madagascar") #We focus on mainland sub-Saharan Africa

#===============================================================================
#Sites per date plot ---- FIGURE 1.1

date_freq  <- dateInfo %>% 
  count(siteID, sort=TRUE) %>%
  rename(n_dates = n)  
  #count(n_dates, sort=TRUE) %>% 
  #rename(n_sites = n)
  
pdf(file=here('output','figures','figure1.1.pdf'), width=8.5, height=7)
ggplot(date_freq, aes(x=n_dates)) +
  geom_bar() +
  scale_y_continuous(name="Number of sites", breaks=seq(0, 350, 25)) +
  scale_x_continuous(name="Number of dates per site", breaks=seq(1, 29, 2)) 
dev.off()

#===============================================================================
##Bayesian Quantile Regression

## quantreg.R Figures ----

# Bayesian Quantile Regression Model with Measurement Error Plot ---- FIGURE 1
## Compute Fitted Model Confidence Intervals

# rq and median calibrated date
rq.ci <- predict.rq(fit_rq, newdata=data.frame(dist_org=0:4600), interval='confidence') #Furthest site from proposed origin, Ngoume, is 4561km

# Bayesian model 
qr.ch1 <- do.call(rbind,quantreg_sample)
post.alpha.quantreg <- qr.ch1[,'alpha']
post.beta.quantreg <- qr.ch1[,'beta']
post.alpha.beta.quantreg  <- data.frame(alpha=post.alpha.quantreg, beta=post.beta.quantreg)
post.theta.quantreg  <- qr.ch1[,grep('theta', colnames(qr.ch1))]
post.theta.med <- apply(post.theta.quantreg, 2, median)
post.quantreg <- apply(post.alpha.beta.quantreg, 1, function(x){x[1]-x[2]*0:4600})
post.ci <- t(apply(post.quantreg, 1, quantile, c(0.025,0.5,0.975)))

## Plot and compare
# Transparency color utility function
col.alpha <- function(x,a=1){xx=col2rgb(x)/255;return(rgb(xx[1],xx[2],xx[2],a))}

pdf(file=here('output','figures','figure1.pdf'), width=8.5, height=7)
plot(NULL, xlim=c(0,4600), ylim=c(4100,100), axes=F, xlab='Distance from Ngoume Site (in km)', ylab='Cal BP') 
rect(xleft=-200, xright=4600, ybottom=2720, ytop=2350, col=col.alpha('grey',0.2), border=NA) #Demarcating Calibration Plateau Region
abline(h=2720,lty=4)
abline(h=2350,lty=4)
axis(1, at=c(0,500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500)) #X-axis
axis(2, at=c(0,500, 1000, 1500, 2000, 2500, 3000, 3500, 4000)) #Y-axis left
axis(4, at=BCADtoBP(c(-1800, -1400, -1000, -600, -200, 200, 600, 1000, 1400, 1800)), labels=c('1800BC', '1400BC', '1000BC', '600BC', '200BC', '200AD', '600AD', '1000AD', '1400AD','1800AD'), cex.axis=0.6) #Y-axis right
points(constants$dist_org, siteInfo$earliest) #Median Calibrated Date #Option: col=siteInfo$dataorigin
points(constants$dist_org, post.theta.med, pch=20) #Median Posterior theta
for (i in 1:nrow(siteInfo))
{
  lines(rep(constants$dist_org[i],2), c(siteInfo$earliest[i],post.theta.med[i]),lty=2)
}
lines(0:4600, rq.ci[,1], lty=1, lwd=2, col='blue') #Quantile Regression on Median Dates
polygon(x=c(0:4600, 4600:0), c(rq.ci[,2], rev(rq.ci[,3])), col=col.alpha('lightblue', 0.4), border=NA)

lines(0:4600, post.ci[,2], lty=1, lwd=2, col='indianred') #Bayesian Quantile Regression with Measurement Error
polygon(x=c(0:4600, 4600:0), c(post.ci[,1], rev(post.ci[,3])), col=col.alpha('indianred', 0.2), border=NA)

text(x=4000, y=2450, labels='Calibration Plateau')
legend('bottomright', legend=c('Median Calibrated Date', TeX('Median Posterior $\\theta$'), 'Quantile Regression on Median Dates', 'Bayesian Quantile Regression with Measurement Error'), pch=c(1,20,NA,NA), lwd=c(NA,NA,2,2), col=c(1,1,'blue','indianred'), cex=0.8)
box()
dev.off()


#-------------------------------------------------------------------------------
# Posterior Dispersal Rate of non-spatial quantile regression Plot ---- FIGURE 2
pdf(here('output','figures','figure2.pdf'),height=5,width=5.5)
postHPDplot(1/post.beta.quantreg, xlab='km/year', ylab='Probability Density', prob=.90, main=TeX('Posterior of $1/\\beta_1$'))
dev.off()

#-------------------------------------------------------------------------------
##Prior predictive checks 

#Prior Predictive Check beta0, beta1 ---- FIGURE 2.1

nsim <- 5000
beta0.prior <- rnorm(nsim, mean=3300, sd=200)
beta1.prior  <- rexp(nsim, rate=1)
slope  <-  beta1.prior
beta0.prior  <- beta0.prior[which((1/slope)>0)] #Ensuring beta0 is positive
slope  <- slope[which((1/slope)>0)] #Ensuring dispersal rate is always positive
nsim2  <- length(slope)
dists  <- -100:4600
slope.mat = matrix(NA, nrow=nsim2, ncol=length(dists))
for (i in 1:nsim2)
{
  slope.mat[i,] <- beta0.prior[i] - slope[i]*c(dists)	
}

pdf(file=here('output','figures','figure2.1.pdf'), width=6, height=6)
plot(NULL, xlim=c(0,4500), ylim=c(3400,1300), type='n', xlab='Distance (km)', ylab='Cal BP', axes=F)
axis(1, at=c(0,500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500), cex.axis=0.9) #X-axis
axis(2, at=seq(3400, 1400, -400))
axis(4, at=BCADtoBP(c(-1400, -1000, -600, -200, 200, 600)), labels=c('1400BC','1000BC','600BC','200BC','200AD','600AD'), cex.axis=0.9)
box()
polygon(x=c(dists, rev(dists)), y=c(apply(slope.mat, 2, quantile,prob=0.025), rev(apply(slope.mat, 2, quantile, prob=0.975))), border=NA, col=rgb(0.67,0.84,0.9,0.5))
polygon(x=c(dists, rev(dists)), y=c(apply(slope.mat, 2, quantile,prob=0.25), rev(apply(slope.mat, 2, quantile, prob=0.75))), border=NA, col=rgb(0.25,0.41,0.88,0.5))

abline(a=3300, b=-1/0.5, lty=2)
text(x=500, y=1600, label='0.5km/yr')

abline(a=3300, b=-1, lty=2)
text(x=1700, y=1900, label='1km/yr')

abline(a=3300, b=-1/3, lty=2)
text(x=3500, y=2000, label='3km/yr')

abline(a=3300, b=-1/5, lty=2)
text(x=2250, y=2950, label='5km/yr')

legend('bottomright', legend=c('50% percentile range', '95% percentile range'), fill=c(rgb(0.67,0.84,0.9,0.5), rgb(0.25,0.41,0.88,0.5)))
dev.off()




#===============================================================================
###Gaussian Process Quantile Regression
#===============================================================================

##Effect of parameter variation in covariance matrix
# Variance Covariance Function Parameters ---- FIGURE 3
pdf(here('output','figures','figure3.pdf'), width=10, height=5)

layout(matrix(c(1,2), nrow = 1, ncol = 2, byrow = TRUE))

plot(NULL, xlim=c(0,500), ylim=c(0.00, 0.05), xlab=TeX('$d_{i,j}$'), ylab=TeX('$cov(i,j)$'))

#Default parameter values
etasq_def = 0.05
rho_def = 150
#Varying parameter values
rho = c(50, 100, 150, 200, 250)
etasq = c(0.025, 0.05, 0.075, 0.1)

cl <- rainbow(5) #Assigning random colours
for (i in 1:length(rho))
{
  curve(etasq_def*exp(-0.5*(x/rho[i])^2), from=0, to=500, col=cl[i], add = TRUE)
}
abline(h=0.05,lty=2)
lines(x=c(rho_def, rho_def), y=c(-0.01, etasq_def*exp(-0.5)),lty=2)
lines(x=c(-20, rho_def), y=c(etasq_def*exp(-0.5), etasq_def*exp(-0.5)),lty=2)

text(400, 0.046, label=TeX('$\\eta^2 = 0.05$'))
text(75, etasq_def*exp(-0.5)-0.004, label=TeX('$\\eta^2e^{-0.5}$'))
text(75, 0.006, label=TeX('$\\rho = 50$'), cex=0.7)
text(168, 0.006, label=TeX('$\\rho = 100$'), cex=0.7)
text(260, 0.006, label=TeX('$\\rho = 150$'), cex=0.85)
text(365, 0.006, label=TeX('$\\rho = 200$'), cex=0.7)
text(465, 0.006, label=TeX('$\\rho = 250$'), cex=0.7)

plot(NULL, xlim=c(0,500), ylim=c(0.00, 0.1), xlab=TeX('$d_{i,j}$'), ylab=TeX('$cov(i,j)$'), cex.sub = 0.6) #sub=TeX("The quadratic exponential model and its relationship to the parameters $\\rho$, $\\eta$, and $d_{i,j}$, assuming ${i \\neq j}$.")
for (i in 1:length(etasq))
{
  curve(etasq[i]*exp(-0.5*(x/rho_def)^2), from=0, to=500, col=cl[i], add = TRUE)
}
text(450, 0.096, label=TeX('$\\rho = 150$'))
text(40, 0.02, label=TeX('$\\eta^2 = 0.025$'), cex=0.7)
text(40, 0.045, label=TeX('$\\eta^2 = 0.05$'), cex=0.7)
text(40, 0.07, label=TeX('$\\eta^2 = 0.075$'), cex=0.7)
text(40, 0.095, label=TeX('$\\eta^2 = 0.1$'), cex=0.7)

dev.off()


#-------------------------------------------------------------------------------
# Impact of rho and etasq on variability in dispersal rate ---- FIGURE 4

# Fixed parameters
n = 600
origin.point = c(11.40, 5.48)
beta0 = 3000
beta1 = 0.5
sigma = 100
seed = 144231

# Sweep parameters
etasq = c(0.01, 0.05, 0.15)
rho = c(10, 200, 600)
cov_param = expand.grid(etasq=etasq, rho=rho)
tmp  <- vector('list',length=nrow(cov_param))
out  <- vector('list',length=nrow(cov_param))


# Simulate sites and plot
for (i in 1:nrow(cov_param))
{
  tmp[[i]]  <- gpqrSim(win = win,
                       n = n,
                       beta0 = beta0,
                       beta1 = beta1,
                       sigma = sigma,
                       origin.point = origin.point,
                       etasq = cov_param$etasq[i],
                       rho =cov_param$rho[i],
                       seed = seed)
  
  out[[i]] <- ggplot() +
    geom_sf(data=win, aes(), fill='grey66', show.legend=FALSE, lwd=0) +
    geom_sf(data=tmp[[i]], mapping = aes(fill=rate), pch=21, col='darkgrey', size=1.5) +
    xlim(-15,50) + #Center the frame on sub-Saharan Africa
    ylim(-35,30) +
    labs(title=TeX(sprintf("$\\eta^2 = %g \\, \\rho = %g$", cov_param$etasq[i], cov_param$rho[i])), fill='Dispersal Rate \n (km/yr)') +
    scale_fill_viridis(option = 'turbo', limits = c(0.7, 4), oob = scales::squish) +
    theme(plot.title = element_text(hjust = 0.5, size=11), panel.background = element_rect(fill='lightblue'), panel.grid.major = element_line(size = 0.1), legend.position=c(0.2, 0.2), legend.text = element_text(size=7), legend.key.width= unit(0.1, 'in'), legend.key.size = unit(0.08, "in"), legend.background=element_rect(fill = alpha("white", 0.5)), legend.title=element_text(size=7), axis.text=element_blank(), axis.ticks=element_blank(), plot.margin = unit(c(0,0,0,0), "in"))

}

pdf(file=here('output','figures','figure4.pdf'), width=8, height=8)
grid.arrange(out[[1]], out[[2]], out[[3]], out[[1]], out[[4]], out[[7]], nrow=2, ncol=3)
dev.off()

# pdf(file=here('output','figures','figure4.pdf'), width=8, height=8)
# grid.arrange(grobs=out, nrow=3, ncol=3)
# dev.off()

#-------------------------------------------------------------------------------
# Relationship between slope and rate of dispersal ---- FIGURE 5
pdf(here('output','figures','figure5.pdf'), width=5, height=5)
par(mar=c(5,6,2,2))
curve(-1/x, from=-2.5, to=-0.05, xlab=TeX('$Slope \\, s - \\beta_1$'), ylab=TeX('$Speed \\, \\frac{-1}{s-\\beta_1}$'), axes=TRUE, sub=TeX("The relationship between regression slope and its negative reciprocal, i.e. the rate of dispersal."), cex.sub = 0.6)
dev.off()


#-------------------------------------------------------------------------------
##Prior predictive checks 

#Prior Predictive Check beta0, beta1, s ---- FIGURE 6
nsim <- 5000
beta0.prior <- rnorm(nsim, mean=3300, sd=200)
beta1.prior  <- rexp(nsim, rate=2)
s.prior  <- rnorm(nsim, mean=0, sd=sqrt(rexp(nsim, rate=20)))
slope  <- s.prior - beta1.prior
beta0.prior  <- beta0.prior[which((-1/slope)>0)] #Ensuring beta0 is positive
slope  <- slope[which((-1/slope)>0)] #Ensuring dispersal rate is always positive
nsim2  <- length(slope)
dists  <- -100:4600
slope.mat = matrix(NA, nrow=nsim2, ncol=length(dists))
for (i in 1:nsim2)
{
  slope.mat[i,] <- beta0.prior[i] + slope[i]*c(dists)	
}

pdf(file=here('output','figures','figure6.pdf'), width=6, height=6)
plot(NULL, xlim=c(0,4500), ylim=c(3400,1300), type='n', xlab='Distance (km)', ylab='Cal BP', axes=F)
axis(1, at=c(0,500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500), cex.axis=0.9) #X-axis
axis(2, at=seq(3400, 1400, -400))
axis(4, at=BCADtoBP(c(-1400, -1000, -600, -200, 200, 600)), labels=c('1400BC','1000BC','600BC','200BC','200AD','600AD'), cex.axis=0.9)
box()
polygon(x=c(dists, rev(dists)), y=c(apply(slope.mat, 2, quantile,prob=0.025), rev(apply(slope.mat, 2, quantile, prob=0.975))), border=NA, col=rgb(0.67,0.84,0.9,0.5))
polygon(x=c(dists, rev(dists)), y=c(apply(slope.mat, 2, quantile,prob=0.25), rev(apply(slope.mat, 2, quantile, prob=0.75))), border=NA, col=rgb(0.25,0.41,0.88,0.5))

abline(a=3300, b=-1/0.5, lty=2)
text(x=500, y=1600, label='0.5km/yr')

abline(a=3300, b=-1, lty=2)
text(x=1700, y=1900, label='1km/yr')

abline(a=3300, b=-1/3, lty=2)
text(x=3500, y=2000, label='3km/yr')

abline(a=3300, b=-1/5, lty=2)
text(x=2250, y=2950, label='5km/yr')

legend('bottomright', legend=c('50% percentile range', '95% percentile range'), fill=c(rgb(0.67,0.84,0.9,0.5), rgb(0.25,0.41,0.88,0.5)))
dev.off()

#-------------------------------------------------------------------------------
# Prior Predictive Check etasq and rho ---- FIGURE 7
nsim  <- 1000
etasq.prior  <- rexp(nsim, 20)
rho.prior  <- rtgamma(nsim,10,(10-1)/200, 1, 4600)
cov.mat = matrix(NA, nrow=nsim, ncol=length(0:1000))
for (i in 1:nsim)
{
  cov.mat[i,] = etasq.prior[i]*exp(-0.5*(0:1000/rho.prior[i])^2)
}

pdf(file=here('output', 'figures', 'figure7.pdf'), width=6, height=5)
plot(NULL, xlab='Distance (km)', ylab='Covariance', xlim=c(0,1000), ylim=c(0,0.2))
polygon(c(0:1000, 1000:0), c(apply(cov.mat, 2, quantile, 0.025), rev(apply(cov.mat, 2, quantile,0.975))), border=NA, col=rgb(0.67,0.84,0.9,0.5))
polygon(c(0:1000, 1000:0), c(apply(cov.mat, 2, quantile, 0.5), rev(apply(cov.mat, 2, quantile,0.75))), border=NA, col=rgb(0.25,0.41,0.88,0.5))

legend('topright', legend=c('50% percentile range','95% percentile range'), fill=c(rgb(0.67,0.84,0.9,0.5), rgb(0.25,0.41,0.88,0.5)))
dev.off()

#===============================================================================
## gpqr_tactical.R Figures

#--------------
# Tactical Simulation ---- FIGURE 8

# Load data
load(here('data', 'tactical_sim_gpqr.RData'))
load(here('output', 'gpqr_tactsim.RData'))

gpqr_tactsim_post  <- do.call(rbind, gpqr_tactsim)
tactsim_post_s  <- gpqr_tactsim_post[ , paste0('s[',1:nrow(sim_sites),']')]
tactsim_post_beta1  <- gpqr_tactsim_post[ ,'beta1']
tactsim_post_rate  <-  -1/(tactsim_post_s-tactsim_post_beta1)
sim_sites$pred.rate  <- apply(tactsim_post_rate,2,median)


s8a  <- ggplot() +
  geom_sf(data = win, aes(), fill='grey66', show.legend=FALSE, lwd=0) +
  geom_sf(data = sim_sites, mapping = aes(fill=rate), pch=21, col='darkgrey', size=3) + 
  xlim(-15,50) + #Center the frame on sub-Saharan Africa
  ylim(-35,30) +
  labs(fill='Dispersal Rate \n (km/yr)') + 
  scale_fill_viridis(option="turbo", limits=c(1,4), oob = scales::squish) +
  annotate("text", x = 45, y = 28, label = TeX('$\\beta_0 = 3300$')) +
  annotate("text", x = 45, y = 25, label = TeX('$\\beta_1 = 0.4$')) +
  annotate("text", x = 45, y = 22, label = TeX('$\\eta^2 = 0.02$')) +
  annotate("text", x = 45, y = 19, label = TeX('$\\rho = 350$')) +
  ggtitle('Simulated Dispersal Rates') +
  theme(legend.position = c(0.2, 0.3), legend.background=element_rect(fill = alpha("white",0.5)), axis.title.x=element_blank(), axis.title.y=element_blank())


s8b  <- ggplot() +
  geom_sf(data = win, aes(), fill='grey66', show.legend=FALSE, lwd=0) +
  geom_sf(data = sim_sites, mapping = aes(fill=pred.rate), pch=21, col='darkgrey', size=3) + 
  xlim(-15,50) + #Center the frame on sub-Saharan Africa
  ylim(-35,30) +
  labs(fill='Dispersal Rate \n (km/yr)') + 
  scale_fill_viridis(option="turbo", limits=c(1,4), oob = scales::squish) +
  ggtitle('Predicted Dispersal Rates') +
  theme(legend.position = c(0.2, 0.3), legend.background=element_rect(fill = alpha("white",0.5)), axis.title.x=element_blank(), axis.title.y=element_blank())


pdf(here('output', 'figures','figure8.pdf'), width=10, height=7)
grid.arrange(s8a, s8b, ncol=2)
dev.off()

#-------------------------------------------------------------------------------
# Posterior vs True values of s for Tactical Simulation ---- FIGURE 9

load(here('output','gpqr_tactsim.RData'))
gpqr_tactsim_post  <- do.call(rbind,gpqr_tactsim)
tactsim_post_s  <- gpqr_tactsim_post[,paste0('s[',1:nrow(sim_sites),']')]
tactsim_post_s_med  <- apply(tactsim_post_s, 2, median)
tactsim_post_s_lo  <- apply(tactsim_post_s, 2, quantile, 0.025)
tactsim_post_s_hi  <- apply(tactsim_post_s, 2, quantile, 0.975)
rr = c(min(c(tactsim_post_s_lo, sim_sites$s)), max(c(tactsim_post_s_hi, sim_sites$s)))

pdf(here('output', 'figures', 'figure9.pdf'), height=6, width=6)
plot(NULL, xlim=rr, ylim=rr, xlab='Simulated s', ylab='Predicted s')
points(sim_sites$s, tactsim_post_s_med, pch=20)
for (i in 1:nrow(sim_sites))
{
  lines(x=c(sim_sites$s[i], sim_sites$s[i]), y=c(tactsim_post_s_lo[i], tactsim_post_s_hi[i]))
}
abline(a=0, b=1, lty=2, col='red')
dev.off()

#-------------------------------------------------------------------------------
# Posterior vs True values of beta0,beta1,rho,etasq for Tactical Simulation ---- FIGURE 10

tactsim_post_beta0  <- gpqr_tactsim_post[,'beta0']
tactsim_post_beta1 <- gpqr_tactsim_post[,'beta1']
tactsim_post_etasq  <- gpqr_tactsim_post[,'etasq']
tactsim_post_rho  <- gpqr_tactsim_post[,'rho']
true_beta0_with_tau09  <- qnorm(0.9, true_param$beta0, true_param$sigma) 


pdf(here('output', 'figures', 'figure10.pdf'), height=8, width=8)
par(mfrow=c(2,2))
postHPDplot(tactsim_post_beta0, xlab='Cal BP', ylab='Posterior Probability', main=TeX('$\\beta_0$'), prob = 0.95)
abline(v=true_beta0_with_tau09, lty=2)
postHPDplot(tactsim_post_beta1, xlab='', ylab='Posterior Probability', main=TeX('$\\beta_1$'), prob=0.95)
abline(v=true_param$beta1, lty=2)
postHPDplot(tactsim_post_etasq, xlab='', ylab='Posterior Probability', main=TeX('$\\eta^2$'), prob=0.95)
abline(v=true_param$etasq, lty=2)
postHPDplot(tactsim_post_rho, xlab='km', ylab='Posterior Probability', main=TeX('$\\rho$'), prob=0.95)
abline(v=true_param$rho, lty=2)
dev.off()



#===============================================================================
## qpqr_tau90.R and qpqr_tau99.R Figures ----

# Load data
load(here('output','gpqr_tau90.RData'))
load(here('output','gpqr_tau99.RData'))

#-------------------------------------------------------------------------------
# Posterior Mean of dispersal rate deviations ---- FIGURE 11

post.gpqr.tau90  <- do.call(rbind,gpqr_tau90)
post.gpqr.tau99  <- do.call(rbind,gpqr_tau99)

nmcmc  <- nrow(post.gpqr.tau90) #number of MCMC samples (same for tau90 and tau99)
post.s.tau90  <- post.gpqr.tau90[,grep('s\\[',colnames(post.gpqr.tau90))]
post.s.tau99  <- post.gpqr.tau99[,grep('s\\[',colnames(post.gpqr.tau99))]

#post.arrival <- matrix(NA,nmcmc,constants$N.sites)
post.rate.tau90 <- post.rate.tau99  <-  matrix(NA,nmcmc,constants$n_sites)
for (i in 1:nmcmc)
{
  post.rate.tau90[i,] = -1 / (post.s.tau90[i,]-post.gpqr.tau90[i,'beta1'])
  post.rate.tau99[i,] = -1 / (post.s.tau99[i,]-post.gpqr.tau99[i,'beta1'])
}

sites@data$s.m.tau90 <- apply(post.s.tau90,2,median)
sites@data$s.m.tau99 <- apply(post.s.tau99,2,median)
sites@data$rate.m.tau90  <- apply(post.rate.tau90,2,median)
sites@data$rate.m.tau99  <- apply(post.rate.tau99,2,median)

sites.sf <- as(sites,'sf')

f2a <- ggplot() +
  geom_text(data=data.frame(x=49, y=32, label='A'), aes(x=x, y=y, label=label)) +
  geom_sf(data=win, aes(), fill='grey66', show.legend=FALSE, lwd=0) +
  geom_sf(data=sites.sf, mapping = aes(fill=s.m.tau90), pch=21, col='black', size=2) + 
  xlim(-15,50) + #Center the frame on sub-Saharan Africa
  ylim(-35,35) +
  labs(title=TeX(r"(Posterior Median of s with $\tau = 0.90$)"), fill='s') + 
  scale_fill_gradient2(low='blue', high='red', mid='white') +
  theme(plot.title = element_text(hjust = 0.5, vjust = -1.5, size=10), 
        panel.background = element_rect(fill='lightblue'), 
        panel.grid.major = element_line(size = 0.1), 
        legend.position=c(0.2,0.2), 
        legend.text = element_text(size=6), 
        legend.key.width= unit(0.1, 'in'), 
        legend.key.size = unit(0.08, "in"), 
        legend.background=element_rect(fill = alpha("white", 0.5)), 
        legend.title=element_text(size=7), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.margin = unit(c(0,0,0,0), "in"), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())

f2b <- ggplot() +
  geom_text(data=data.frame(x=49, y=32, label='B'), aes(x=x, y=y, label=label)) +
  geom_sf(data=win, aes(), fill='grey66', show.legend=FALSE, lwd=0) +
  geom_sf(data=sites.sf, mapping = aes(fill=s.m.tau99), pch=21, col='black', size=2) +
  xlim(-15,50) + #Center the frame on sub-Saharan Africa
  ylim(-35,35) +
  labs(title=TeX(r"(Posterior Median of s with $\tau = 0.99$)"), fill='s') +
  scale_fill_gradient2(low='blue', high='red', mid='white') +
  theme(plot.title = element_text(hjust = 0.5, vjust = -1.5, size=10),
        panel.background = element_rect(fill='lightblue'),
        panel.grid.major = element_line(size = 0.1),
        legend.position=c(0.2,0.2),
        legend.text = element_text(size=6),
        legend.key.width= unit(0.1, 'in'),
        legend.key.size = unit(0.08, "in"),
        legend.background=element_rect(fill = alpha("white", 0.5)),
        legend.title=element_text(size=7),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        plot.margin = unit(c(0,0,0,0), "in"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

f2c <- ggplot() +
  geom_text(data=data.frame(x=49, y=32, label='C'), aes(x=x, y=y, label=label)) +
  geom_sf(data=win, aes(), fill='grey66', show.legend=FALSE, lwd=0) +
  geom_sf(data=sites.sf, mapping = aes(fill=rate.m.tau90), pch=21, col='black', size=2) + 
  xlim(-15,50) + #Center the frame on sub-Saharan Africa
  ylim(-35,35) +
  labs(title=TeX(r"(Posterior median of dispersal rate with $\tau = 0.90$)"), fill='Dispersal Rate (km/year)') +
  scale_fill_viridis(option="turbo", limits=c(0,4)) +
  theme(plot.title = element_text(hjust = 0.5, vjust = -1.5, size=10), 
        panel.background = element_rect(fill='lightblue'), 
        panel.grid.major = element_line(size = 0.1), 
        legend.position=c(0.21,0.2), 
        legend.text = element_text(size=6), 
        legend.key.width= unit(0.1, 'in'), 
        legend.key.size = unit(0.08, "in"), 
        legend.background=element_rect(fill = alpha("white", 0.5)), 
        legend.title=element_text(size=7), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.margin = unit(c(0,0,0,0), "in"), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())

f2d <- ggplot() +
  geom_text(data=data.frame(x=49, y=32, label='D'), aes(x=x, y=y, label=label)) +
  geom_sf(data=win, aes(), fill='grey66', show.legend=FALSE, lwd=0) +
  geom_sf(data=sites.sf, mapping = aes(fill=rate.m.tau99), pch=21, col='black', size=2) +
  xlim(-15,50) + #Center the frame on sub-Saharan Africa
  ylim(-35,35) +
  labs(title=TeX(r"(Posterior median of dispersal rate with $\tau = 0.99$)"), fill='Dispersal Rate (km/year)') +
  scale_fill_viridis(option="turbo", limits=c(0,4)) +
  theme(plot.title = element_text(hjust = 0.5, vjust = -1.5, size=10),
        panel.background = element_rect(fill='lightblue'),
        panel.grid.major = element_line(size = 0.1),
        legend.position=c(0.21,0.2),
        legend.text = element_text(size=6),
        legend.key.width= unit(0.1, 'in'),
        legend.key.size = unit(0.08, "in"),
        legend.background=element_rect(fill = alpha("white", 0.5)),
        legend.title=element_text(size=7),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        plot.margin = unit(c(0,0,0,0), "in"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())


pdf(file=here('output','figures','figure11.pdf'), width=7, height=7)
grid.arrange(f2a, f2b, f2c, f2d, ncol=2, padding=0)
dev.off()

#-------------------------------------------------------------------------------
# Traceplot of beta0, beta1, rhosq, and etasq for tau=0.9) ---- FIGURE 12

pdf(file=here('output','figures','figure12.pdf'), width=8, height=8)
par(mfrow=c(2,2))
traceplot(gpqr_tau90[,'beta0'], main=TeX('$\\beta_0$'), smooth=TRUE)
traceplot(gpqr_tau90[,'beta1'], main=TeX('$\\beta_1$'), smooth=TRUE)
traceplot(gpqr_tau90[,'rho'], main=TeX('$\\rho$'), smooth=TRUE)
traceplot(gpqr_tau90[,'etasq'], main=TeX('$\\eta^2$'), smooth=TRUE)
dev.off()

#-------------------------------------------------------------------------------
# Traceplot of beta0, beta1, rhosq, and etasq for tau=0.99 ---- FIGURE 13

pdf(file=here('output','figures','figure13.pdf'), width=8, height=8)
par(mfrow=c(2,2))
traceplot(gpqr_tau99[,'beta0'], main=TeX('$\\beta_0$'), smooth=TRUE)
traceplot(gpqr_tau99[,'beta1'], main=TeX('$\\beta_1$'), smooth=TRUE)
traceplot(gpqr_tau99[,'rho'], main=TeX('$\\rho$'), smooth=TRUE)
traceplot(gpqr_tau99[,'etasq'], main=TeX('$\\eta^2$'), smooth=TRUE)
dev.off()

#-------------------------------------------------------------------------------
# Marginal posteriors of beta0, beta1, rho, etasq for tau = 0.9 ---- FIGURE 14

gpqr.tau90.comb  <- do.call(rbind, gpqr_tau90)

pdf(file=here('output','figures','figure14.pdf'), width=8, height=8)
par(mfrow=c(2,2))
postHPDplot(gpqr.tau90.comb[,'beta0'], main=TeX('$\\beta_0$'), xlab='Cal BP', ylab='')
postHPDplot(gpqr.tau90.comb[,'beta1'], main=TeX('$\\beta_1$'), xlab='', ylab='')
postHPDplot(gpqr.tau90.comb[,'rho'], main=TeX('$\\rho$'), xlab='km', ylab='')
postHPDplot(gpqr.tau90.comb[,'etasq'], main=TeX('$\\eta^2$'), xlab='', ylab='')
dev.off()

#-------------------------------------------------------------------------------
# Marginal posteriors of beta0, beta1, rho, etasq for tau = 0.99 ---- FIGURE 15

gpqr.tau99.comb  <- do.call(rbind, gpqr_tau99)

pdf(file=here('output','figures','figure15.pdf'), width=8, height=8)
par(mfrow=c(2,2))
postHPDplot(gpqr.tau99.comb[,'beta0'], main=TeX('$\\beta_0$'), xlab='Cal BP', ylab='')
postHPDplot(gpqr.tau99.comb[,'beta1'], main=TeX('$\\beta_1$'), xlab='', ylab='')
postHPDplot(gpqr.tau99.comb[,'rho'], main=TeX('$\\rho$'), xlab='km', ylab='')
postHPDplot(gpqr.tau99.comb[,'etasq'], main=TeX('$\\eta^2$'), xlab='', ylab='')
dev.off()


#===============================================================================
###Bayesian Hierarchical Phase Model
#===============================================================================
##Tactical Simulation of Bayesian Hierarchical Phase Model

#Load Data ----
load(here("output", "phasemodel_tactsim.RData"))

#-------------------------------------------------------------------------------
# Tactical Simulation Posterior Predictive Check for nu and upsilon ---- FIGURE 16

#For models (i) and (ii) select parameters a and b (i.e. start and end date of occupation in the region)
post.model.i  <- do.call(rbind, mcmc.samples1)[ , c(1,2)]
post.model.ii  <- do.call(rbind, mcmc.samples2)[ , c(1,27)]

dens.i.nu  <- density(post.model.i[,1],bw = 5)
dens.i.upsilon  <- density(post.model.i[,2],bw=5)
dens.ii.nu  <- density(post.model.ii[,1],bw=5)
dens.ii.upsilon  <- density(post.model.ii[,2],bw=5)

pdf(file=here('output','figures','figure16.pdf'), width=8, height=8)

plot(NULL, xlim=c(4100,2700), ylim=c(0,0.022), xlab='Cal BP', ylab='Posterior Probability') 
polygon(c(dens.i.nu$x, rev(dens.i.nu$x)), c(rep(0,length(dens.i.nu$x)), rev(dens.i.nu$y)), border=NA, col=rgb(0,0.4,0,0.5))
polygon(c(dens.i.upsilon$x, rev(dens.i.upsilon$x)), c(rep(0,length(dens.i.upsilon$x)), rev(dens.i.upsilon$y)), border=NA, col=rgb(0,0.4,0,0.5))
polygon(c(dens.ii.nu$x, rev(dens.ii.nu$x)), c(rep(0,length(dens.ii.nu$x)), rev(dens.ii.nu$y)), border=NA, col=rgb(1,0.55,0,0.5))
polygon(c(dens.ii.upsilon$x, rev(dens.ii.upsilon$x)), c(rep(0,length(dens.ii.upsilon$x)), rev(dens.ii.upsilon$y)), border=NA, col=rgb(1,0.55,0,0.5))
abline(v=c(3700, 3200),lty=2)
axis(3,at=c(3700, 3200),labels=c(TeX('$a$'),TeX('$b$')))
legend('topright', legend=c('Non hierarchichal','Hierarchichal'), fill=c('darkgreen','darkorange'))

dev.off()

#-------------------------------------------------------------------------------
# Traceplot of start and end of occupation (a, b) ---- FIGURE 16.2

pdf(file=here('output','figures','figure16.2.pdf'), width=8, height=8)
par(mfrow=c(2,2))
traceplot(mcmc.samples1[,'a'], main=TeX('$a$'), smooth=TRUE)
traceplot(mcmc.samples1[,'b'], main=TeX('$b$'), smooth=TRUE)
traceplot(mcmc.samples2[,'a'], main=TeX('$a$'), smooth=TRUE)
traceplot(mcmc.samples2[,'b'], main=TeX('$b$'), smooth=TRUE)
dev.off()


#===============================================================================
##Bayesian Hierarchical Phase Models without constraints

#Classify hex areas
East_hex <- c(35, 31, 30, 29, 28, 25, 24, 19, 18, 14, 10, 9, 6, 5, 3)
West_hex <- c(26, 22, 21, 17, 16, 12)
No_site_hex <- c(1, 2, 4, 7, 8, 11, 15, 20, 32, 33, 37, 38 ,39, 40, 41, 42, 44, 45, 46, 47)
Origin_hex <- 34

#Hex areas with and without out sites
Hex_without_sites <- c(8, 11, 15, 20, 32, 33, 37:42)
Hex_with_sites <- which(rep(1:43) %!in% c(8, 11, 15, 20, 32, 33, 37:42))

#Load Data ----
load(here("output", "phase_model_a.RData"))
load(here("output","phase_model_b.RData"))

#-------------------------------------------------------------------------------
# Prior Predictive check for duration parameter, delta ---- FIGURE 17
nsim  <- 5000

set.seed(123)

gamma1  <- runif(nsim,1,20)
gamma2  <- rtruncnorm(nsim, mean=200, sd=100, 1, 500)
delta.mat = matrix(NA, ncol=1000, nrow=nsim) #Initialise

for (i in 1:nsim) {
  delta.mat[i,] = dgamma(1:1000, gamma1[i], (gamma1[i]-1)/gamma2[i])
  }

pdf(file=here('output','figures','figure17.pdf'), height=6, width=6)

plot(NULL,xlab=TeX('$\\delta$'),ylab='Probability Density',xlim=c(1,1000),ylim=c(0,0.02))
polygon(x=c(1:1000, 1000:1), y=c(apply(delta.mat,2,quantile,prob=0.025), rev(apply(delta.mat,2,quantile,prob=0.975))), border=NA, col=rgb(0.67,0.84,0.9,0.5))
polygon(x=c(1:1000, 1000:1), y=c(apply(delta.mat,2,quantile,prob=0.25), rev(apply(delta.mat,2,quantile,prob=0.75))), border=NA, col=rgb(0.25,0.41,0.88,0.5))
legend('topright', legend=c('50% percentile range', '95% percentile range'), fill=c(rgb(0.67,0.84,0.9,0.5), rgb(0.25,0.41,0.88,0.5)))

dev.off()


#-------------------------------------------------------------------------------
# Marginal Posterior Distribution of nu, model a ---- FIGURE 18

out.comb.unif.model.a  <- do.call(rbind, out_unif_model_a)
post.nu.model.a  <- out.comb.unif.model.a[,paste0('a[',1:57,']')] %>%  round() #57 hex areas
model.a.long  <- data.frame(value = as.numeric(post.nu.model.a),
                            area = rep(1:57, each=nrow(post.nu.model.a)))

for(i in 1:nrow(model.a.long)){
  if(model.a.long$area[i] %in% East_hex){
    model.a.long$stream[i] = "East"
  } else if(model.a.long$area[i] %in% West_hex){
    model.a.long$stream[i] = "West"
  } else if(model.a.long$area[i] %in% No_site_hex){
    model.a.long$stream[i] = "No sites"
  } else if(model.a.long$area[i] == Origin_hex){
    model.a.long$stream[i] = "Origin"
  } else {
    model.a.long$stream[i] = "Neither"
  }
}

model.a.long  <- model.a.long %>%
  mutate(area = factor(area, levels=paste0(1:57), ordered=TRUE),
         stream = factor(stream, levels=c("Origin", "East", "West", "Neither", "No sites"))) %>% 
  filter(area %in% c(3,5,6,8:43)) #Don't plot areas where we know no Bantu Expansion took place


#Plot
pdf(file=here('output','figures','figure18.pdf'), height=10, width=8)

ggplot(model.a.long, aes(x = value, y = area, fill=stream)) + 
  geom_density_ridges() +
  scale_x_reverse(limits=c(5000, 150), breaks=BCADtoBP(c(-3000, -2600, -2200, -1800, -1400, -1000, -600, -200, 200, 600, 1000, 1400, 1800)), labels=c('3000BC', '2600BC', '2200BC', '1800BC', '1400BC', '1000BC', '600BC', '200BC', '200AD', '600AD', '1000AD', '1400AD','1800AD')) +
  scale_fill_viridis_d(name = "Stream") +
  xlab(paste('Arrival time,', TeX('$a_k$'))) +
  ylab(paste('Area,', TeX('$k$')))

dev.off()

#-------------------------------------------------------------------------------
# Probability Matrix of nu, model a ---- FIGURE 19
source(here('src','orderPPlot.R'))
post.nu.model.a_rel  <- out.comb.unif.model.a[,paste0('a[',c(3,5,6,8:43),']')] %>%  round() #Keep 47 relevant hex areas

pdf(file=here('output','figures','figure19.pdf'), width=10, height=10.5)
orderPPlot(post.nu.model.a_rel, name.vec=paste("Area", c(3,5,6,8:43)))
dev.off()


#-------------------------------------------------------------------------------
##Bayesian Hierarchical Phase Models with wave-of-advance constraints

# Marginal Posterior Distribution of nu and upsilon, model b ---- FIGURE 21

out.comb.unif.model.b  <- do.call(rbind, out.unif.model_b)
post.nu.model.b  <- out.comb.unif.model.b[,paste0('a[',1:57,']')] %>%  round()
model.b.long  <- data.frame(value = as.numeric(post.nu.model.b),
                            area = rep(1:57, each=nrow(post.nu.model.b)))

for(i in 1:nrow(model.b.long)){
  if(model.b.long$area[i] %in% East_hex){
    model.b.long$stream[i] = "East"
  } else if(model.b.long$area[i] %in% West_hex){
    model.b.long$stream[i] = "West"
  } else if(model.b.long$area[i] %in% No_site_hex){
    model.b.long$stream[i] = "No sites"
  } else if(model.b.long$area[i] == Origin_hex){
    model.b.long$stream[i] = "Origin"
  } else {
    model.b.long$stream[i] = "Neither"
  }
}

model.b.long  <- model.b.long %>%
  mutate(area = factor(area, levels=paste0(1:57), ordered=TRUE),
         stream = factor(stream, levels=c("Origin", "East", "West", "Neither", "No sites"))) %>% 
  filter(area %in% c(3,5,6,8:43)) #Don't plot areas where we know no Bantu Expansion took place


#Plot
pdf(file=here('output','figures','figure21.pdf'), height=10, width=8)

ggplot(model.b.long, aes(x = value, y = area, fill=stream)) + 
  geom_density_ridges() +
  scale_x_reverse(limits=c(5000, 150), breaks=BCADtoBP(c(-3000, -2600, -2200, -1800, -1400, -1000, -600, -200, 200, 600, 1000, 1400, 1800)), labels=c('3000BC', '2600BC', '2200BC', '1800BC', '1400BC', '1000BC', '600BC', '200BC', '200AD', '600AD', '1000AD', '1400AD','1800AD')) +
  scale_fill_viridis_d(name = "Stream") +
  xlab(paste('Arrival time,', TeX('$a_k$'))) +
  ylab(paste('Area,', TeX('$k$')))

dev.off()

#-------------------------------------------------------------------------------
# Estimated Arrival Date ---- FIGURE 22
# Setup Functions and Variables
extract <- function(x)
{
  tmp = do.call(rbind, x)
  tmp2 = tmp[ , grep('^a\\[',colnames(tmp))]
  qta = apply(tmp2, 2, quantile, prob=c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))
  return(qta)
}

post.bar <- function(x, i, h, col)
{
  # 	lines(c(x[1],x[7]),c(i,i),col=col)
  rect(xleft=x[2], xright=x[6], ybottom=i-h/5, ytop=i+h/5, border=NA, col=col)
  rect(xleft=x[3], xright=x[5], ybottom=i-h/3, ytop=i+h/3, border=NA, col=col)
  lines(c(x[4],x[4]), c(i-h/2,i+h/2), lwd=2, col='grey44')
}


main.col <- c(rgb(68, 1, 84, maxColorValue=255), #Origin
              rgb(59, 82, 139, maxColorValue=255), #East
              rgb(33, 144, 140, maxColorValue=255), #West
              rgb(93, 200, 99, maxColorValue=255), #Neither
              rgb(253, 231, 37, maxColorValue=255)) #No sites


pdf(file=here('output','figures','figure22.pdf'), width=10, height=18, pointsize=4)

# Posterior Arrival Times
plot(NULL, xlim=c(5000,150), ylim=c(3,125), xlab=paste('Arrival time,', TeX('$a_k$')), ylab=paste('Area,', TeX('$k$')), cex.lab = 2, axes=F)
tmp.a = extract(out_unif_model_a)
tmp.b = extract(out.unif.model_b)
iseq.a = seq(2,by=3,length.out=43)
iseq.b = seq(1,by=3,length.out=43)
abline(h=seq(3,by=3,length.out=42), col='darkgrey',lty=2)

for (i in c(3,5,6,8:43)) #Don't plot areas where we know no Bantu Expansion took place
{
  #Asign colour index
  if(i %in% East_hex){ci <- 2
  } else if(i %in% West_hex){ci <- 3
  } else if(i %in% No_site_hex){ci <- 5
  } else if(i == Origin_hex){ci <- 1
  } else {ci <- 4
  }
  #Plot bar in area i
  post.bar(tmp.a[,i], i=iseq.a[i], h=0.9, col=main.col[ci])
  post.bar(tmp.b[,i], i=iseq.b[i], h=0.9, col=main.col[ci])
}


axis(2, at=iseq.a+0.5, labels = paste0(1:43), las=2, cex.axis=1.7)
axis(1, at = BCADtoBP(c(-3000, -2200, -1400, -600, 200, 1000, 1800)), labels=c('3000BC', '2200BC', '1400BC', '600BC', '200AD', '1000AD', '1800AD'), tck=-0.01, cex.axis=1.7)
axis(3, at = seq(5000, 150, -600), labels=paste0(seq(5000,150,-600),'BP'), tck=-0.01, cex.axis=1.7)
axis(1, at = BCADtoBP(c(-3000, -2600, -2200, -1800, -1400, -1000, -600, -200, 200, 600, 1000, 1400, 1800)), labels=NA, tck=-0.01) #Minor tick marks
axis(3, at = seq(5000, 150, -300), labels=NA, tck=-0.01) #Minor tick marks
box()

post.bar(c(800,700,600,500,400,300,200), i=1.5, h=0.9, col='lightgrey')
arrows(x0=700, x1=300, y0=0.3, y1=0.3, angle = 90, code = 3, length = 0.01)
arrows(x0=600, x1=400, y0=0.8, y1=0.8, angle = 90, code = 3, length = 0.01)
text(x=785, y=0.9, "50% HPDI", cex=1.5)
text(x=875, y=0.07,"90% HPDI", cex=1.5)
text(x=820, y=2.3, "Median Posterior", cex=1.5)
lines(x=c(550,500), y=c(2.35, 2.28))
text(x=5000, y=2, 'Model A', cex=1.5)
text(x=5000, y=1, 'Model B', cex=1.5)
legend(x = 500, y = 122,        
       legend = c("Origin", "East", "West", "Neither", "No sites"),
       pch=15,        
       col = main.col,
       title = "Stream",
       cex=1.9,
       bty = "n") 
dev.off()


#-------------------------------------------------------------------------------
## Plot HEX areas with median arrival times ---- FIGURE map_figure3

#Extract arrival times for model A
out.comb.unif.modela  <- do.call(rbind, out_unif_model_a)
post.nu.modela  <- out.comb.unif.modela[,paste0('a[',1:43,']')]  %>% round() 
hpdi.modela  <- apply(post.nu.modela, 2, function(x){HPDinterval(as.mcmc(x), prob = .90)}) 
med.modela  <- apply(post.nu.modela, 2, median)
hi90_modA  <- hpdi.modela[1,]
lo90_modA  <- hpdi.modela[2,]

median_hex_dates_modA <- hex_area_win %>% 
  filter(area_ID %in% 1:43) %>% 
  mutate(median_date = med.modela,
         hpdi_high = hi90_modA,
         hpdi_low = lo90_modA,
         contains_sites = as.factor(case_when(area_ID %in% Hex_with_sites ~ 1, area_ID %in% Hex_without_sites ~ 0))) %>% 
  filter(area_ID %!in% c(1, 2, 3, 4, 7))

#Extract arrival times for model B
out.comb.unif.modelb  <- do.call(rbind, out.unif.model_b)
post.nu.modelb  <- out.comb.unif.modelb[,paste0('a[',1:43,']')]  %>% round()
hpdi.modelb  <- apply(post.nu.modelb, 2, function(x){HPDinterval(as.mcmc(x), prob = .90)}) 
med.modelb  <- apply(post.nu.modelb, 2, median)
hi90_modB  <- hpdi.modelb[1,]
lo90_modB  <- hpdi.modelb[2,]

median_hex_dates_modB <- hex_area_win %>% 
  filter(area_ID %in% 1:43) %>% 
  mutate(median_date = med.modelb,
         hpdi_high = hi90_modB,
         hpdi_low = lo90_modB) %>% 
  filter(area_ID %!in% c(1, 2, 3, 4, 7))


#Plot
#-----MODEL A
modA <- ggplot(data = median_hex_dates_modA) +
          geom_sf(data = st_buffer(st_as_sf(sampling_win, crs = 4326), 40000), aes(color = "grey50")) + #sampling window with coastal buffer
          geom_sf(aes(fill = median_date, alpha=contains_sites)) + #hex grid
          scale_fill_viridis_c(option="F", direction=-1) +
          scale_alpha_manual(values=c(0.45, 1)) +
          xlab('Longitude') +
          ylab('Latitude') +
          geom_sf_label(aes(label = ifelse(contains_sites==0, NA, paste0(median_date, "BP"))), label.size  = NA, alpha = 0.4, size=3.5) + #hex grid labels
          theme(panel.background = element_rect(fill = "lightblue",
                                                colour = "lightblue",
                                                size = 0.5,
                                                linetype = "solid"),
                legend.position = "none")




modAHPDIlow <- ggplot(data = median_hex_dates_modA) +
  geom_sf(data = st_buffer(st_as_sf(sampling_win, crs = 4326), 40000), aes(color = "grey50")) + #sampling window with coastal buffer
  geom_sf(aes(fill = hpdi_low, alpha=contains_sites)) + 
  ggtitle('90% HPDI (low)') +
  scale_fill_viridis_c(option="F", direction=-1) +
  scale_alpha_manual(values=c(0.45, 1)) +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(size=12),
        legend.position = "none")

modAHPDIhigh <- ggplot(data = median_hex_dates_modA) +
  geom_sf(data = st_buffer(st_as_sf(sampling_win, crs = 4326), 40000), aes(color = "grey50")) + #sampling window with coastal buffer
  geom_sf(aes(fill = hpdi_high, alpha=contains_sites)) + 
  ggtitle('90% HPDI (high)') +
  scale_fill_viridis_c(option="F", direction=-1) +
  scale_alpha_manual(values=c(0.45, 1)) +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(size=12),
        legend.position = "none")


A <- cowplot::ggdraw() +
  draw_plot(modA) +
  draw_plot(modAHPDIlow, 
            x = .73, y = .285, width = .25, height = .25) +
  draw_plot(modAHPDIhigh, 
            x = .73, y = .06, width = .25, height = .25)

#-----MODEL B
modB <- ggplot(data = median_hex_dates_modB) +
          geom_sf(data = st_buffer(st_as_sf(sampling_win, crs = 4326), 40000), aes(color = "grey50")) + #sampling window with coastal buffer
          geom_sf(aes(fill = median_date)) + #hex grid
          scale_fill_viridis_c(option="F", direction=-1) +
          geom_sf_label(aes(label = paste0(median_date, "BP")), label.size  = NA, alpha = 0.4, size=3.5) + #hex grid labels
          xlab('Longitude') +
          ylab('Latitude') +
          theme(panel.background = element_rect(fill = "lightblue",
                                                colour = "lightblue",
                                                size = 0.5,
                                                linetype = "solid"),
                legend.position = "none")

modBHPDIlow <- ggplot(data = median_hex_dates_modB) +
  geom_sf(data = st_buffer(st_as_sf(sampling_win, crs = 4326), 40000), aes(color = "grey50")) + #sampling window with coastal buffer
  geom_sf(aes(fill = hpdi_low)) + 
  ggtitle('90% HPDI (low)') +
  scale_fill_viridis_c(option="F", direction=-1) +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(size=12),
        legend.position = "none")

modBHPDIhigh <- ggplot(data = median_hex_dates_modB) +
  geom_sf(data = st_buffer(st_as_sf(sampling_win, crs = 4326), 40000), aes(color = "grey50")) + #sampling window with coastal buffer
  geom_sf(aes(fill = hpdi_high)) + 
  ggtitle('90% HPDI (high)') +
  scale_fill_viridis_c(option="F", direction=-1) +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(size=12),
        legend.position = "none")


B <- cowplot::ggdraw() +
  draw_plot(modB) +
  draw_plot(modBHPDIlow, 
            x = .73, y = .285, width = .25, height = .25) +
  draw_plot(modBHPDIhigh, 
            x = .73, y = .06, width = .25, height = .25)

#Output
pdf(file=here('output','figures','map_figure_arrival.pdf'), width=15, height=8)
grid.arrange(A, B, ncol=2, padding=0)
dev.off()



#===============================================================================
#===============================================================================
#Plot Tactical simulation of ICAR model

##Load Data ----
load(here("output", "phasemodel_tactsim_TOY.RData"))
load(here('data', 'tactical_sim_phase_TOY.RData'))
load(here('data','trig.RData')) #nodes and edges between hex area centroids
#Combine constants
constants <- c(constants, constants_trig)

#-------------
##Traceplot of start and end of occupation (a, b) 

pdf(file=here('output', 'figures','figure30.pdf'), width=8, height=8)
par(mfrow=c(2,2))
traceplot(mcmc.samples1[,'a[18]'], main=TeX('$a[18]$'),smooth=TRUE) #region area 18, as an example
traceplot(mcmc.samples1[,'b[18]'], main=TeX('$b[18]$'),smooth=TRUE)
traceplot(mcmc.samples1[,'nabla[18]'], main=TeX('$nabla[18]$'), smooth=TRUE)
dev.off()

pdf(file=here('output', 'figures','figure31.pdf'), width=8, height=8)
par(mfrow=c(2,2))
traceplot(mcmc.samples2[,'a[18]'], main=TeX('$a[18]$'),smooth=TRUE) #region area 18, as an example
traceplot(mcmc.samples2[,'b[18]'], main=TeX('$b[18]$'),smooth=TRUE)
traceplot(mcmc.samples2[,'nabla[18]'], main=TeX('$nabla[18]$'), smooth=TRUE)
dev.off()

#-------------------------------------------------------------------------------
## Tactical Simulation Posterior Predictive Check for a and b in a given region

#For model (i) and (ii) select parameters a and b (i.e. start and end date of occupation in the region)
post.model.i  <- do.call(rbind, mcmc.samples1)[ , c(18, 59)] #area 18 (selecting a[18] and b[18])
post.model.ii  <- do.call(rbind, mcmc.samples2)[ , c(18, 159)] #area 18 (selecting a[18] and b[18]) #to find b[18] index: which(colnames(as.data.frame(mcmc.samples2$chain1)) == 'b[18]')

dens.i.a  <- density(post.model.i[,1],bw = 5)
dens.i.b  <- density(post.model.i[,2],bw=5)
dens.ii.a  <- density(post.model.ii[,1],bw = 5)
dens.ii.b  <- density(post.model.ii[,2],bw=5)

pdf(file=here('output', 'figures','figure32.pdf'), width=8, height=8)
plot(NULL, xlim=c(4000,2400), ylim=c(0,0.022), xlab='Cal BP', ylab='Posterior Probability')
polygon(c(dens.i.a$x, rev(dens.i.a$x)), c(rep(0,length(dens.i.a$x)), rev(dens.i.a$y)), border=NA, col=rgb(0,0.4,0,0.5))
polygon(c(dens.i.b$x, rev(dens.i.b$x)), c(rep(0,length(dens.i.b$x)), rev(dens.i.b$y)), border=NA, col=rgb(0,0.4,0,0.5))
polygon(c(dens.ii.a$x, rev(dens.i.a$x)), c(rep(0,length(dens.ii.a$x)), rev(dens.ii.a$y)), border=NA, col=rgb(1,0.55,0,0.5))
polygon(c(dens.ii.b$x, rev(dens.i.b$x)), c(rep(0,length(dens.ii.b$x)), rev(dens.ii.b$y)), border=NA, col=rgb(1,0.55,0,0.5))
abline(v=c(3700, 3200),lty=2)
axis(3,at=c(3700, 3200),labels=c(TeX('$a$'),TeX('$b$')))
legend('topright', legend=c('Non hierarchichal','Hierarchichal'), fill=c('darkgreen','darkorange'))
dev.off()


#-------------------------------------------------------------------------------
# ## Tactical Simulation Posterior Predictive Check for gradient in a given region
 
# #For model (i) and (ii) select parameters a and b (i.e. start and end date of occupation in the region)
# post.grad.model.i  <- do.call(rbind, mcmc.samples1)[ , 100] #area 18 (selecting nabla[18])
# post.grad.model.ii  <- do.call(rbind, mcmc.samples2)[ , 300]
# 
# dens.i.nabla  <- density(post.grad.model.i,bw = 5)
# dens.ii.nabla  <- density(post.grad.model.ii,bw = 5)
# 
# pdf(file=here('output', 'figures','figure34.pdf'), width=8, height=8)
# plot(NULL, xlim=c(-20,20), ylim=c(0,0.1), xlab='Gradient', ylab='Posterior Probability')
# polygon(c(dens.i.nabla$x, rev(dens.i.nabla$x)), c(rep(0,length(dens.i.nabla$x)), rev(dens.i.nabla$y)), border=NA, col=rgb(0,0.4,0,0.5))
# polygon(c(dens.ii.nabla$x, rev(dens.ii.nabla$x)), c(rep(0,length(dens.ii.nabla$x)), rev(dens.ii.nabla$y)), border=NA, col=rgb(1,0.55,0,0.5))
# legend('topright', legend=c('Non hierarchichal','Hierarchichal'), fill=c('darkgreen','darkorange'))
# dev.off()


#-------------------------------------------------------------------------------
##Plot magnitude and direction of gradients for model (i)

##Setup Functions and Variables
#Extract gradient distribution information
extract_gradinfo <- function(x)
{
  tmp = do.call(rbind, x)
  tmp2 = tmp[ , grep('^nabla\\[',colnames(tmp))]
  return(tmp2)
}

#Determine proportion of distribution which is positive
prop_gthan_zero <- function(data) {
  # Count the number of elements greater than zero
  count_positive <- sum(data > 0)
  # Calculate the proportion
  proportion <- count_positive / length(data)
  return(proportion)
}

# Define a custom function to plot arrows
plot_arrows <- function(edges, ...) {
  for (i in 1:nrow(edges)) {
    
    ##Determine which is the start node and which is the end node
    if(edges$mean_gradient[i] >= 0){
      x_start = edges$region1_x[i]
      y_start = edges$region1_y[i]
      x_end = edges$region2_x[i]
      y_end = edges$region2_y[i] } else {
        x_start = edges$region2_x[i]
        y_start = edges$region2_y[i]
        x_end = edges$region1_x[i]
        y_end = edges$region1_y[i]
      }
    
    ##Determine the angle of the edge
    #Calculate the differences in x and y coordinates
    delta_x <- x_end - x_start
    delta_y <- y_end - y_start
    
    # Calculate the angle in radians using the arctangent function (atan2)
    angle_rad <- atan2(delta_y, delta_x)
    
    # Define the length of the arrow
    arrow_length <- abs(edges$mean_gradient[i])*40 #TODO: relative magnitude of gradient is what is important, '40' is merely a scaling parameter
    
    # Calculate the coordinates of the arrow head
    arrow_head_x <- x_end - arrow_length * cos(angle_rad)
    arrow_head_y <- y_end - arrow_length * sin(angle_rad)
    
    # Define alpha transparency value (0 to 1) depending on uncertainty in gradient
    alpha <- edges$uncertainty[i] #the proportion of the distribution in this direction
    
    arrows(x0 = x_start, y0 = y_start, x1 = arrow_head_x, y1 = arrow_head_y, col = rgb(0, 0, 1, alpha), ...)
  }
}

#--------
#Extract quantile information for models
tmp = extract_gradinfo(mcmc.samples1) 
#tmp.ii = extract_gradinfo(mcmc.samples2) 
qta = apply(tmp, 2, quantile, prob=c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))

#Extract uncertainty information
uncert = sapply(as.data.frame(tmp), prop_gthan_zero)

#Add info to edges dataframe
edge_info <- edge_info %>% 
  mutate(mean_gradient = qta[4,], #50% quantile
         uncertainty = uncert) #% of distribution > zero
#Remove extremely long edges which are 'artificially' created along the internal window boundary or require extreme coastal movement along external window boundary
edge_info <- edge_info %>% 
  filter(distance <= mean(distance)*1.3)

#Create nodes
nodes <- st_coordinates(hex_area_win$area_center)

#--------
#PLOT
# Get the bounding box of the sample window
sample_win_buff <- st_buffer(st_as_sf(sampling_win, crs = 4326), 40000)
bbox <- sf::st_bbox(sample_win_buff)

# Calculate the aspect ratio of the bounding box
aspect_ratio <- diff(range(c(bbox["ymin"], bbox["ymax"]))) / diff(range(c(bbox["xmin"], bbox["xmax"])))
gridsize <- 1 # Adjust the denominator to change grid density


pdf(file=here('output', 'figures','figure33.pdf'), width=8, height=8)
# Create an empty plot with the appropriate range
plot(x = c(bbox["xmin"], bbox["xmax"]), y = c(bbox["ymin"], bbox["ymax"]), type = "n",
     xlab = "Latitude", ylab = "Longitude", main = "Gradient surface of arrival times", asp = 1/aspect_ratio, axes = F)
# Adding theme elements
rect(bbox["xmin"], bbox["ymin"], bbox["xmax"], bbox["ymax"], col = "lightblue", border=NA) #a blue-colored bounding box
# Add background grid
abline(v = seq(ceiling(bbox["xmin"]), floor(bbox["xmax"]), by = gridsize), col = "white")
abline(h = seq(ceiling(bbox["ymin"]), floor(bbox["ymax"]), by = gridsize), col = "white")
# Add labeled axes for the background grid
axis(side = 1, at = seq(ceiling(bbox["xmin"]), floor(bbox["xmax"]), by = gridsize*2),
     labels = format(seq(ceiling(bbox["xmin"]), floor(bbox["xmax"]), by = gridsize*2), nsmall = 1))
axis(side = 2, at = seq(ceiling(bbox["ymin"]), floor(bbox["ymax"]), by = gridsize*2),
     labels = format(seq(ceiling(bbox["ymin"]), floor(bbox["ymax"]), by = gridsize*2), nsmall = 1), las = 1, add = T)
# Add north arrow in the bottom right corner
north_arrow_length <- 4
north_arrow_x <- bbox["xmax"] - 0.05 * diff(c(bbox["xmin"], bbox["xmax"]))
north_arrow_y <- bbox["ymin"] + 0.05 * diff(c(bbox["ymin"], bbox["ymax"]))
arrows(x0 = north_arrow_x, y0 = north_arrow_y, x1 = north_arrow_x, y1 = north_arrow_y + north_arrow_length, 
       length = 0.1, angle = 30, col = "black", lwd=2)
# Plotting the sampling window with coastal buffer
plot(sample_win_buff, col = "grey", border = rgb(0, 0, 0, 0.2), add = TRUE)
# Plotting the hex grid
plot(hex_area_win$geometry, add = TRUE, border = rgb(0, 0, 0, 0.2))
# Plotting the hex origins
points(nodes[,"X"], nodes[,"Y"], pch = 20, col = rgb(0, 0, 0, 0.8), cex = 2)
# Plot gradient arrows using the custom function
plot_arrows(edge_info, lwd=4, length = 0.1) #length parameter defines arrowhead size
dev.off()






