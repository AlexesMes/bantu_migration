# Load Libraries and spatial data ----
library(here)
library(truncnorm)
library(cascsim)
library(corrplot)
library(ggplot2)
library(ggridges)
library(rnaturalearth)
library(nimbleCarbon)
library(rcarbon)
library(maptools)
library(sf)
library(rgeos)
library(viridis)
library(latex2exp)
library(gridExtra)
library(diagram)
library(quantreg)
library(coda)

source(here('src','gpqrSim.R'))

#===============================================================================
## Load Data

# Load Observed Data
load(here('data','c14.RData'))
# Load quantile regression results
load(here('output','quantreg_res.RData'))

#===============================================================================
##Bayesian Quantile Regression

## quantreg.R Figures ----

# Bayesian Quantile Regression Model with Measurement Error Plot ---- FIGURE 1
## Compute Fitted Model Confidence Intervals

# rq and median calibrated date
rq.ci <- predict.rq(fit_rq, newdata=data.frame(dist_org=0:4400), interval='confidence') #Furthest site from proposed origin, Obobogo, is 4400km

# Bayesian model 
qr.ch1 <- do.call(rbind,quantreg_sample)
post.alpha.quantreg <- qr.ch1[,'alpha']
post.beta.quantreg <- qr.ch1[,'beta']
post.alpha.beta.quantreg  <- data.frame(alpha=post.alpha.quantreg, beta=post.beta.quantreg)
post.theta.quantreg  <- qr.ch1[,grep('theta', colnames(qr.ch1))]
post.theta.med <- apply(post.theta.quantreg, 2, median)
post.quantreg <- apply(post.alpha.beta.quantreg, 1, function(x){x[1]-x[2]*0:4400})
post.ci <- t(apply(post.quantreg, 1, quantile, c(0.025,0.5,0.975)))

## Plot and compare
# Transparency color utility function
col.alpha <- function(x,a=1){xx=col2rgb(x)/255;return(rgb(xx[1],xx[2],xx[2],a))}

pdf(file=here('output','figures','figure1.pdf'), width=8.5, height=7)
plot(NULL, xlim=c(0,4400), ylim=c(3600,100), axes=F, xlab='Distance from Obobogo Site (in km)', ylab='Cal BP') 
rect(xleft=-200, xright=4600, ybottom=2720, ytop=2350, col=col.alpha('grey',0.2), border=NA) #Demarcating Calibration Plateau Region
abline(h=2720,lty=4)
abline(h=2350,lty=4)
axis(1, at=c(0,500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500)) #X-axis
axis(2) #Y-axis left
axis(4, at=BCADtoBP(c(-1400, -1000, -600, -200, 200, 600, 1000, 1400, 1800)), labels=c('1400BC', '1000BC', '600BC', '200BC', '200AD', '600AD', '1000AD', '1400BC','1800BC'), cex.axis=0.7) #Y-axis right
points(constants$dist_org, siteInfo$earliest) #Median Calibrated Date
points(constants$dist_org, post.theta.med, pch=20) #Median Posterior theta
for (i in 1:nrow(siteInfo))
{
  lines(rep(constants$dist_org[i],2),c(siteInfo$earliest[i],post.theta.med[i]),lty=2)
}
lines(0:4400, rq.ci[,1], lty=1, lwd=2, col='blue') #Quantile Regression on Median Dates
polygon(x=c(0:4400,4400:0),c(rq.ci[,2],rev(rq.ci[,3])), col=col.alpha('lightblue', 0.4), border=NA)

lines(0:4400, post.ci[,2], lty=1, lwd=2, col='indianred') #Bayesian Quantile Regression with Measurement Error
polygon(x=c(0:4400,4400:0), c(post.ci[,1], rev(post.ci[,3])), col=col.alpha('indianred', 0.2), border=NA)

text(x=4000, y=2450, labels='Calibration Plateau')
legend('bottomright', legend=c('Median Calibrated Date', TeX('Median Posterior $\\theta$'), 'Quantile Regression on Median Dates', 'Bayesian Quantile Regression with Measurement Error'), pch=c(1,20,NA,NA), lwd=c(NA,NA,2,2), col=c(1,1,'blue','indianred'), cex=0.8)
box()
dev.off()


#-------------------------------------------------------------------------------
# Posterior Dispersal Rate of non-spatial quantile regression Plot ---- FIGURE 2
pdf(here('output','figures','figure2.pdf'),height=5,width=5.5)
postHPDplot(1/post.beta.quantreg, xlab='km/year', ylab='Probability Density', prob=.90, main=TeX('Posterior of $1/\\beta_1$'))
dev.off()



#===============================================================================
##Gaussian Process Quantile Regression

##Effect of parameter variation in covariance matrix

# Variance Covariance Function Parameters ---- FIGURE 3
pdf(here('output','figures','figure3.pdf'), width=5, height=5)
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

plot(NULL, xlim=c(0,500), ylim=c(0.00, 0.1), xlab=TeX('$d_{i,j}$'), ylab=TeX('$cov(i,j)$'), sub=TeX("The quadratic exponential model and its relationship to the parameters $\\rho$, $\\eta$, and $d_{i,j}$, assuming ${i \\neq j}$."), cex.sub = 0.6)
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

#Come back to this figur once the tactical simulations using gpqrSim.R is working ....
#Try: df <- ne_countries(continent = "Africa", returnclass = "sf") #See: https://github.com/ropensci/rnaturalearth/issues/34 ... 








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
beta0.prior <- rnorm(nsim, mean=3000, sd=200)
beta1.prior  <- rexp(nsim, rate=1)
s.prior  <- rnorm(nsim, mean=0, sd=sqrt(rexp(nsim, rate=20)))
slope  <- s.prior - beta1.prior
beta0.prior  <- beta0.prior[which((-1/slope)>0)] #Ensuring beta0 is positive
slope  <- slope[which((-1/slope)>0)] #Ensuring dispersal rate is always positive
nsim2  <- length(slope)
dists  <- -100:4400
slope.mat = matrix(NA, nrow=nsim2, ncol=length(dists))
for (i in 1:nsim2)
{
  slope.mat[i,] <- beta0.prior[i] + slope[i]*c(dists)	
}

pdf(file=here('output','figures','figure6.pdf'), width=6, height=6)
plot(NULL, xlim=c(0,4400), ylim=c(3400,1300), type='n', xlab='Distance (km)', ylab='Cal BP', axes=F)
axis(1, at=c(0,500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500), cex.axis=0.9) #X-axis
axis(2, at=seq(3400, 1400, -400))
axis(4, at=BCADtoBP(c(-1400, -1000, -600, -200, 200, 600)), labels=c('1400BC','1000BC','600BC','200BC','200AD','600AD'), cex.axis=0.9)
box()
polygon(x=c(dists, rev(dists)), y=c(apply(slope.mat, 2, quantile,prob=0.025), rev(apply(slope.mat, 2, quantile, prob=0.975))), border=NA, col=rgb(0.67,0.84,0.9,0.5))
polygon(x=c(dists, rev(dists)), y=c(apply(slope.mat, 2, quantile,prob=0.25), rev(apply(slope.mat, 2, quantile, prob=0.75))), border=NA, col=rgb(0.25,0.41,0.88,0.5))

abline(a=3000, b=-1/0.5, lty=2)
text(x=1050, y=1600, label='0.5km/yr')

abline(a=3000, b=-1, lty=2)
text(x=1400, y=1900, label='1km/yr')

abline(a=3000, b=-1/3, lty=2)
text(x=3500, y=2000, label='3km/yr')

abline(a=3000, b=-1/5, lty=2)
text(x=1350, y=2850, label='5km/yr')

legend('bottomright', legend=c('50% percentile range', '95% percentile range'), fill=c(rgb(0.67,0.84,0.9,0.5), rgb(0.25,0.41,0.88,0.5)))
dev.off()

#-------------------------------------------------------------------------------
# Prior Predictive Check etasq and rho ---- FIGURE 7
nsim  <- 1000
etasq.prior  <- rexp(nsim,20)
rho.prior  <- rtgamma(nsim,10,(10-1)/150, 1, 4400)
cov.mat = matrix(NA, nrow=nsim, ncol=length(0:1000))
for (i in 1:nsim)
{
  cov.mat[i,] = etasq.prior[i]*exp(-0.5*(0:1000/rho.prior[i])^2)
}

pdf(file=here('output','figures','figure7.pdf'), width=6, height=5)
plot(NULL, xlab='Distance (km)', ylab='Covariance', xlim=c(0,1000), ylim=c(0,0.2))
polygon(c(0:1000, 1000:0), c(apply(cov.mat, 2, quantile, 0.025), rev(apply(cov.mat, 2, quantile,0.975))), border=NA, col=rgb(0.67,0.84,0.9,0.5))
polygon(c(0:1000, 1000:0), c(apply(cov.mat, 2, quantile, 0.5), rev(apply(cov.mat, 2, quantile,0.75))), border=NA, col=rgb(0.25,0.41,0.88,0.5))

legend('topright', legend=c('50% percentile range','95% percentile range'), fill=c(rgb(0.67,0.84,0.9,0.5), rgb(0.25,0.41,0.88,0.5)))
dev.off()


#--------------
## qpqr_tau90.R Figures ----












