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
#library(quantreg)
library(coda)
library(graphics)
library(ggthemes)
library(gganimate) #animate ggplot
library(gifski) #save output as gif


source(here('src','orderPPlot.R'))
source(here('src','gpqrSim.R'))
source(here('src', 'grad_funcs.R'))

`%!in%` <- Negate(`%in%`)

#===============================================================================
## Load Data

# Load Observed Data
load(here('data','eastc14.RData')) #East and Southern Africa

#Load nodes and edges between hex area centroids
load(here('data','trig.RData'))

# Load quantile regression results
load(here('output','quantreg_res.RData'))
#===============================================================================
#Sites per date plot ---- FIGURE 1

date_freq  <- dateInfo %>% 
  count(siteID, sort=TRUE) %>%
  rename(n_dates = n)  
  #count(n_dates, sort=TRUE) %>% 
  #rename(n_sites = n)
  
pdf(file=here('output','figures','figure1.pdf'), width=8.5, height=7)
ggplot(date_freq, aes(x=n_dates)) +
  geom_bar() +
  scale_y_continuous(name="Number of sites", breaks=seq(0, 350, 25)) +
  scale_x_continuous(name="Number of dates per site", breaks=seq(1, 29, 2)) 
dev.off()

#===============================================================================
##Bayesian Quantile Regression

## quantreg.R Figures ----

# Bayesian Quantile Regression Model with Measurement Error Plot ---- FIGURE 2
## Compute Fitted Model Confidence Intervals

# rq and median calibrated date
rq.ci <- predict.rq(fit_rq, newdata=data.frame(dist_org=0:3500), interval='confidence') #Furthest site from proposed origin, Ngoume, is 4561km

# Bayesian model 
qr.ch1 <- do.call(rbind,quantreg_sample)
post.alpha.quantreg <- qr.ch1[,'alpha']
post.beta.quantreg <- qr.ch1[,'beta']
post.alpha.beta.quantreg  <- data.frame(alpha=post.alpha.quantreg, beta=post.beta.quantreg)
post.theta.quantreg  <- qr.ch1[,grep('theta', colnames(qr.ch1))]
post.theta.med <- apply(post.theta.quantreg, 2, median)
post.quantreg <- apply(post.alpha.beta.quantreg, 1, function(x){x[1]-x[2]*0:3500})
post.ci <- t(apply(post.quantreg, 1, quantile, c(0.025,0.5,0.975)))

## Plot and compare
# Transparency color utility function
col.alpha <- function(x,a=1){xx=col2rgb(x)/255;return(rgb(xx[1],xx[2],xx[2],a))}

pdf(file=here('output','figures','figure2.pdf'), width=8.5, height=7)
plot(NULL, xlim=c(0,3500), ylim=c(3600,100), axes=F, xlab='Distance from Katuruka Site (in km)', ylab='Cal BP') 
rect(xleft=-200, xright=3500, ybottom=2720, ytop=2350, col=col.alpha('grey',0.2), border=NA) #Demarcating Calibration Plateau Region
abline(h=2720,lty=4)
abline(h=2350,lty=4)
axis(1, at=c(0,500, 1000, 1500, 2000, 2500, 3000, 3500)) #X-axis
axis(2, at=c(0,500, 1000, 1500, 2000, 2500, 3000, 3500)) #Y-axis left
axis(4, at=BCADtoBP(c(-1400, -1000, -600, -200, 200, 600, 1000, 1400, 1800)), labels=c('1400BC', '1000BC', '600BC', '200BC', '200AD', '600AD', '1000AD', '1400AD','1800AD'), cex.axis=0.6) #Y-axis right
points(constants$dist_org, siteInfo$earliest) #Median Calibrated Date #Option: col=siteInfo$dataorigin
points(constants$dist_org, post.theta.med, pch=20) #Median Posterior theta
for (i in 1:nrow(siteInfo))
{
  lines(rep(constants$dist_org[i],2), c(siteInfo$earliest[i],post.theta.med[i]),lty=2)
}
lines(0:3500, rq.ci[,1], lty=1, lwd=2, col='blue') #Quantile Regression on Median Dates
polygon(x=c(0:3500, 3500:0), c(rq.ci[,2], rev(rq.ci[,3])), col=col.alpha('lightblue', 0.4), border=NA)

lines(0:3500, post.ci[,2], lty=1, lwd=2, col='indianred') #Bayesian Quantile Regression with Measurement Error
polygon(x=c(0:3500, 3500:0), c(post.ci[,1], rev(post.ci[,3])), col=col.alpha('indianred', 0.2), border=NA)

text(x=3000, y=2450, labels='Calibration Plateau')
legend('bottomright', legend=c('Median Calibrated Date', TeX('Median Posterior $\\theta$'), 'Quantile Regression on Median Dates', 'Bayesian Quantile Regression with Measurement Error'), pch=c(1,20,NA,NA), lwd=c(NA,NA,2,2), col=c(1,1,'blue','indianred'), cex=0.8)
box()
dev.off()


#-------------------------------------------------------------------------------
# Posterior Dispersal Rate of non-spatial quantile regression Plot ---- FIGURE 3
pdf(here('output','figures','figure3.pdf'),height=5,width=5.5)
postHPDplot(1/post.beta.quantreg, xlab='km/year', ylab='Probability Density', xlim = c(0,25), prob=.90, main=TeX('Posterior of $1/\\beta_1$'))
dev.off()

#-------------------------------------------------------------------------------
##Prior predictive checks 

#Prior Predictive Check beta0, beta1 ---- FIGURE 4
nsim <- 5000
beta0.prior <- rnorm(nsim, mean=2500, sd=200)
beta1.prior  <- rexp(nsim, rate=1)
slope  <-  beta1.prior
beta0.prior  <- beta0.prior[which((1/slope)>0)] #Ensuring beta0 is positive
slope  <- slope[which((1/slope)>0)] #Ensuring dispersal rate is always positive
nsim2  <- length(slope)
dists  <- -100:3500
slope.mat = matrix(NA, nrow=nsim2, ncol=length(dists))
for (i in 1:nsim2)
{
  slope.mat[i,] <- beta0.prior[i] - slope[i]*c(dists)	
}

pdf(file=here('output','figures','figure4.pdf'), width=6, height=6)
plot(NULL, xlim=c(0,3500), ylim=c(3000, 1300), type='n', xlab='Distance (km)', ylab='Cal BP', axes=F)
axis(1, at=c(0,500, 1000, 1500, 2000, 2500, 3000, 3500), cex.axis=0.9) #X-axis
axis(2, at=seq(3000, 1400, -400))
axis(4, at=BCADtoBP(c(-1000, -600, -200, 200, 600)), labels=c('1000BC','600BC','200BC','200AD','600AD'), cex.axis=0.9)
box()
polygon(x=c(dists, rev(dists)), y=c(apply(slope.mat, 2, quantile,prob=0.025), rev(apply(slope.mat, 2, quantile, prob=0.975))), border=NA, col=rgb(0.67,0.84,0.9,0.5))
polygon(x=c(dists, rev(dists)), y=c(apply(slope.mat, 2, quantile,prob=0.25), rev(apply(slope.mat, 2, quantile, prob=0.75))), border=NA, col=rgb(0.25,0.41,0.88,0.5))

abline(a=2500, b=-1/0.5, lty=2)
text(x=500, y=1600, label='0.5km/yr')

abline(a=2500, b=-1, lty=2)
text(x=600, y=1800, label='1km/yr')

abline(a=2500, b=-1/3, lty=2)
text(x=1200, y=2000, label='3km/yr')

abline(a=2500, b=-1/5, lty=2)
text(x=2500, y=2150, label='5km/yr')

legend('bottomright', legend=c('50% percentile range', '95% percentile range'), fill=c(rgb(0.67,0.84,0.9,0.5), rgb(0.25,0.41,0.88,0.5)))
dev.off()

#-------------------------------------------------------------------------------
#Traceplots
pdf(file=here('output', 'figures','figure4.1.pdf'), width=8, height=8)
par(mfrow=c(2,2))
traceplot(quantreg_sample[, "alpha"], main=TeX('$alpha$'),smooth=TRUE)
traceplot(quantreg_sample[,'beta'], main=TeX('$beta$'),smooth=TRUE)
traceplot(quantreg_sample[,'sigma'], main=TeX('$sigma$'), smooth=TRUE)
traceplot(quantreg_sample[,'theta[14]'], main=TeX('$theta[14]$'), smooth=TRUE) #example theta
dev.off()


#===============================================================================
###Bayesian Hierarchical Phase Model
#===============================================================================
##Tactical Simulation of Bayesian Hierarchical Phase Model

#Load Data ----
load(here("output", "phasemodel_tactsim.RData"))

#-------------------------------------------------------------------------------
# Tactical Simulation Posterior Predictive Check for nu and upsilon ---- FIGURE 5

#For models (i) and (ii) select parameters a and b (i.e. start and end date of occupation in the region)
post.model.i  <- do.call(rbind, mcmc.samples1)[ , c(1,2)]
post.model.ii  <- do.call(rbind, mcmc.samples2)[ , c(1,27)]

dens.i.nu  <- density(post.model.i[,1],bw = 5)
dens.i.upsilon  <- density(post.model.i[,2],bw=5)
dens.ii.nu  <- density(post.model.ii[,1],bw=5)
dens.ii.upsilon  <- density(post.model.ii[,2],bw=5)

pdf(file=here('output','figures','figure5.pdf'), width=8, height=8)

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
# Traceplot of start and end of occupation (a, b) ---- FIGURE 6

pdf(file=here('output','figures','figure6.pdf'), width=8, height=8)
par(mfrow=c(2,2))
traceplot(mcmc.samples1[,'a'], main=TeX('$a$'), smooth=TRUE)
traceplot(mcmc.samples1[,'b'], main=TeX('$b$'), smooth=TRUE)
traceplot(mcmc.samples2[,'a'], main=TeX('$a$'), smooth=TRUE)
traceplot(mcmc.samples2[,'b'], main=TeX('$b$'), smooth=TRUE)
dev.off()


#===============================================================================
##Bayesian Hierarchical Phase Models without constraints

#Hex areas with and without out sites
Hex_with_sites <- unique(siteInfo$area_id)
Hex_without_sites <- which(rep(1:41) %!in% Hex_with_sites)

#Load Data ----
load(here("output", "phase_model_a.RData"))

#-------------------------------------------------------------------------------
# Prior Predictive check for duration parameter, delta ---- FIGURE 7
nsim  <- 5000

set.seed(123)

gamma1  <- runif(nsim,1,20)
gamma2  <- rtruncnorm(nsim, mean=200, sd=100, 1, 500)
delta.mat = matrix(NA, ncol=1000, nrow=nsim) #Initialise

for (i in 1:nsim) {
  delta.mat[i,] = dgamma(1:1000, gamma1[i], (gamma1[i]-1)/gamma2[i])
  }

pdf(file=here('output','figures','figure7.pdf'), height=6, width=6)

plot(NULL,xlab=TeX('$\\delta$'),ylab='Probability Density',xlim=c(1,1000),ylim=c(0,0.02))
polygon(x=c(1:1000, 1000:1), y=c(apply(delta.mat,2,quantile,prob=0.025), rev(apply(delta.mat,2,quantile,prob=0.975))), border=NA, col=rgb(0.67,0.84,0.9,0.5))
polygon(x=c(1:1000, 1000:1), y=c(apply(delta.mat,2,quantile,prob=0.25), rev(apply(delta.mat,2,quantile,prob=0.75))), border=NA, col=rgb(0.25,0.41,0.88,0.5))
legend('topright', legend=c('50% percentile range', '95% percentile range'), fill=c(rgb(0.67,0.84,0.9,0.5), rgb(0.25,0.41,0.88,0.5)))

dev.off()


#-------------------------------------------------------------------------------
# Marginal Posterior Distribution of nu, model a ---- FIGURE 8

out.comb.unif.model.a  <- do.call(rbind, out_unif_model_a)
post.nu.model.a  <- out.comb.unif.model.a[,paste0('a[',1:41,']')] %>%  round() #41 hex areas
model.a.long  <- data.frame(value = as.numeric(post.nu.model.a),
                            area = rep(1:41, each=nrow(post.nu.model.a)))

model.a.long  <- model.a.long %>%
  mutate(area = factor(area, levels=paste0(1:41), ordered=TRUE)) %>% 
  filter(area %!in% c(1, 2, 3, 4, 6)) #The Bantu hadn't settled in this area by the time the dutch arrived in the Cape (~1600AD). To back this up there are no EIA sites in these regions.


#Plot
pdf(file=here('output','figures','figure8.pdf'), height=10, width=8)

ggplot(model.a.long, aes(x = value, y = area, fill='orange')) + 
  geom_density_ridges() +
  scale_x_reverse(limits=c(2670, 1200), breaks=BCADtoBP(c(-700, -500, -300, -100, 100, 300, 500, 700)), labels=c('700BC', '500BC', '300BC', '100BC', '100AD', '300AD', '500AD', '700AD')) +
  scale_fill_manual(values='orange') +
  xlab(paste('Arrival time,', TeX('$a_k$'))) +
  ylab(paste('Area,', TeX('$k$'))) +
  theme(legend.position = "none")

dev.off()

#-------------------------------------------------------------------------------
# Probability Matrix of nu, model a ---- FIGURE 9
source(here('src','orderPPlot.R'))

post.nu.model.a_rel  <- out.comb.unif.model.a[,paste0('a[',c(5,7:41),']')] %>%  round() #Keep relevant hex areas

pdf(file=here('output','figures','figure9.pdf'), width=10, height=10.5)
orderPPlot(post.nu.model.a_rel, name.vec=paste("Area", c(5,7:41)))
dev.off()

#-------------------------------------------------------------------------------
# Estimated Arrival Date ---- FIGURE 10
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
  rect(xleft=x[2], xright=x[6], ybottom=i-h/5, ytop=i+h/5, border=NA, col=col)
  rect(xleft=x[3], xright=x[5], ybottom=i-h/3, ytop=i+h/3, border=NA, col=col)
  lines(c(x[4],x[4]), c(i-h/2,i+h/2), lwd=2, col='grey44')
}


# Posterior Arrival Times
pdf(file=here('output','figures','figure10.pdf'), width=10, height=18, pointsize=4)
plot(NULL, xlim=c(2670, 1200), ylim=c(3,70), xlab=paste('Arrival time,', TeX('$a_k$')), ylab=paste('Area,', TeX('$k$')), cex.lab = 2, axes=F)
tmp.a = extract(out_unif_model_a)
iseq.a = seq(2,by=2,length.out=36)
abline(h=seq(1,by=2,length.out=36), col='darkgrey',lty=2)

counter <- 1 #indexing counter
for (i in c(5,7:41)) #all relevant hex areas
{
  #Plot bar in area i
  post.bar(tmp.a[,i], i=iseq.a[counter], h=0.9, col='lightblue')
  counter <- counter + 1
}

axis(2, at=iseq.a+0.5, labels = paste0(c(5, 7:41)), las=2, cex.axis=1.7)
axis(1, at = BCADtoBP(c(-700, -500, -300, -100, 100, 300, 500, 700)), labels=c('700BC', '500BC', '300BC', '100BC', '100AD', '300AD', '500AD', '700AD'), tck=-0.01, cex.axis=1.7)
axis(3, at = seq(2700, 1200, -200), labels=paste0(seq(2700, 1200, -200),'BP'), tck=-0.01, cex.axis=1.7)
axis(1, at = BCADtoBP(c(-600, -400, -200, 1, 200, 400, 600)), labels=NA, tck=-0.01) #Minor tick marks
axis(3, at = seq(2700, 1200, -100), labels=NA, tck=-0.01) #Minor tick marks
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
theme(legend.position = "none")
dev.off()


#-------------------------------------------------------------------------------
## Plot HEX areas with median arrival times ---- FIGURE 11

#Extract arrival times for model A
out.comb.unif.modela  <- do.call(rbind, out_unif_model_a)
post.nu.modela  <- out.comb.unif.modela[,paste0('a[',c(5,7:41),']')]  %>% round() 
hpdi.modela  <- apply(post.nu.modela, 2, function(x){HPDinterval(as.mcmc(x), prob = .90)}) 
med.modela  <- apply(post.nu.modela, 2, median)
hi90_modA  <- hpdi.modela[1,]
lo90_modA  <- hpdi.modela[2,]

median_hex_dates_modA <- hex_area_win %>% 
  filter(area_ID %in% c(5,7:41)) %>% 
  mutate(median_date = med.modela,
         hpdi_high = hi90_modA,
         hpdi_low = lo90_modA,
         contains_sites = as.factor(case_when(area_ID %in% Hex_with_sites ~ 1, area_ID %in% Hex_without_sites ~ 0))) 


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

#Output
pdf(file=here('output','figures','figure11.pdf'), width=15, height=8)
cowplot::ggdraw() +
  draw_plot(modA) +
  draw_plot(modAHPDIlow, 
            x = .73, y = .285, width = .25, height = .25) +
  draw_plot(modAHPDIhigh, 
            x = .73, y = .06, width = .25, height = .25)
dev.off()


#===============================================================================
###Bayesian Hierarchical ICAR Model
#===============================================================================
#Tactical simulation of ICAR model 

##Load Data ----
load(here("output", "ICARmodel_tactsim.RData")) #load(here('output', 'Womblemodel_tactsim.RData'))
load(here('data', 'tactical_sim_ICAR.RData')) 

#Combine constants
constants <- c(constants, constants_trig)

#-------------
##Traceplot of start and end of occupation (a, b) -- FIGURE 12 and 13

pdf(file=here('output', 'figures','figure12.pdf'), width=8, height=8)
par(mfrow=c(2,2))
traceplot(mcmc.samples1[,'a[18]'], main=TeX('$a[18]$'),smooth=TRUE) #region area 18, as an example
traceplot(mcmc.samples1[,'b[18]'], main=TeX('$b[18]$'),smooth=TRUE)
#traceplot(mcmc.samples1[,'nabla[18]'], main=TeX('$nabla[18]$'), smooth=TRUE)
dev.off()

pdf(file=here('output', 'figures','figure13.pdf'), width=8, height=8)
par(mfrow=c(2,2))
traceplot(mcmc.samples2[,'a[18]'], main=TeX('$a[18]$'),smooth=TRUE) #region area 18, as an example
traceplot(mcmc.samples2[,'b[18]'], main=TeX('$b[18]$'),smooth=TRUE)
#traceplot(mcmc.samples2[,'nabla[18]'], main=TeX('$nabla[18]$'), smooth=TRUE)
dev.off()

#-------------
#Plotting Priors for ICAR -- exploratory (not currently saved as a fig)

nsim  <- 5000

gamma1  <- runif(nsim, 0.5, 0.5) #runif(nsim, 0.01, 3) 
gamma2  <- runif(nsim, 0.5, 0.55)  #runif(nsim, 0.1, 2) 
tau.mat = matrix(NA, ncol=100, nrow=nsim) #Initialise

for (i in 1:nsim) {
  tau.mat[i,] = dgamma(1:100, shape = gamma1[i], rate = gamma2[i])
}

plot(NULL,xlab=TeX('$\\tau$ (ICAR precision parameter)'),ylab='Probability Density',xlim=c(1,15),ylim=c(0,0.5))
polygon(x=c(1:100, 100:1), y=c(apply(tau.mat,2,quantile,prob=0.025), rev(apply(tau.mat,2,quantile,prob=0.975))), border=NA, col=rgb(0.67,0.84,0.9,0.5))
polygon(x=c(1:100, 100:1), y=c(apply(tau.mat,2,quantile,prob=0.25), rev(apply(tau.mat,2,quantile,prob=0.75))), border=NA, col=rgb(0.25,0.41,0.88,0.5))
legend('topright', legend=c('50% percentile range', '95% percentile range'), fill=c(rgb(0.67,0.84,0.9,0.5), rgb(0.25,0.41,0.88,0.5)))


#gamma1  <- runif(nsim,1,20)
# gamma2  <- rtruncnorm(nsim, mean=200, sd=100, 1, 500)
# delta.mat = matrix(NA, ncol=1000, nrow=nsim) #Initialise
# 
# for (i in 1:nsim) {
#   delta.mat[i,] = dgamma(1:1000, gamma1[i], (gamma1[i]-1)/gamma2[i])
# }

#tau1 <- 1/sigma1^2
#sigma1 ~ dexp(1) #dunif(0,100)
#sigma1= rexp(1,1)) #runif(1,0,100))

#tau1 ~ dgamma(50, 50)

#-------------------------------------------------------------------------------
## Tactical Simulation Posterior Predictive Check for a and b in a given region -- FIGURE 14

#For model (i) and (ii) select parameters a and b (i.e. start and end date of occupation in the region)
sim_a <- constants$true_a
sim_b <- constants$true_b

pdf(file=here('output', 'figures','figure14.pdf'), width=14, height=18)
# Define the layout for the plots
par(mfrow = c(7, 6))

for (k in 1:41) #all hex areas
{
  post.model.i <- do.call(rbind, mcmc.samples0)[ , c(k, k+41)] #selecting a[k] and b[k]
  post.model.ii <- do.call(rbind, mcmc.samples1)[ , c(k, k+141)] #selecting a[k] and b[k]
  post.model.iii <- do.call(rbind, mcmc.samples2)[ , c(k, k+41)] #selecting a[k] and b[k]
  post.model.iv <- do.call(rbind, mcmc.samples3)[ , c(k, k+141)] #selecting a[k] and b[k] #e.g.to find b[18] index: which(colnames(as.data.frame(mcmc.samples2$chain1)) == 'b[18]')
  
  dens.i.a <- density(post.model.i[,1],bw = 5)
  dens.i.b <- density(post.model.i[,2],bw=5)
  dens.ii.a <- density(post.model.ii[,1],bw = 5)
  dens.ii.b <- density(post.model.ii[,2],bw=5)
  dens.iii.a <- density(post.model.iii[,1],bw = 5)
  dens.iii.b <- density(post.model.iii[,2],bw=5)
  dens.iv.a <- density(post.model.iv[,1],bw = 5)
  dens.iv.b <- density(post.model.iv[,2],bw=5)
  
  # Plot
  plot(NULL, xlim=c(sim_a[[k]]+1000, sim_b[[k]]-1000), ylim=c(0,0.022), xlab='Cal BP', ylab='Posterior Probability')
  polygon(c(dens.i.a$x, rev(dens.i.a$x)), c(rep(0,length(dens.i.a$x)), rev(dens.i.a$y)), border=NA, col=rgb(0,0.4,0,0.5))
  #polygon(c(dens.i.b$x, rev(dens.i.b$x)), c(rep(0,length(dens.i.b$x)), rev(dens.i.b$y)), border=NA, col=rgb(0,0.4,0,0.5))
  polygon(c(dens.ii.a$x, rev(dens.ii.a$x)), c(rep(0,length(dens.ii.a$x)), rev(dens.ii.a$y)), border=NA, col=rgb(1,0.55,0,0.5))
  #polygon(c(dens.ii.b$x, rev(dens.ii.b$x)), c(rep(0,length(dens.ii.b$x)), rev(dens.ii.b$y)), border=NA, col=rgb(1,0.55,0,0.5))
  polygon(c(dens.iii.a$x, rev(dens.iii.a$x)), c(rep(0,length(dens.iii.a$x)), rev(dens.iii.a$y)), border=NA, col=rgb(0.82,0.086,0,0.5))
  #polygon(c(dens.iii.b$x, rev(dens.iii.b$x)), c(rep(0,length(dens.iii.b$x)), rev(dens.iii.b$y)), border=NA, col=rgb(0.82,0.086,0,0.5))
  polygon(c(dens.iv.a$x, rev(dens.iv.a$x)), c(rep(0,length(dens.iv.a$x)), rev(dens.iv.a$y)), border=NA, col=rgb(0.004,0,0.82,0.5))
  #polygon(c(dens.iv.b$x, rev(dens.iv.b$x)), c(rep(0,length(dens.iv.b$x)), rev(dens.iv.b$y)), border=NA, col=rgb(0.004,0,0.82,0.5))
  abline(v=sim_a[[k]],lty=2) #abline(v=c(sim_a[[k]], sim_b[[k]]),lty=2)
  axis(3,at=sim_a[[k]],labels=TeX('$a$')) #axis(3,at=c(sim_a[[k]], sim_b[[k]]),labels=c(TeX('$a$'),TeX('$b$')))
  legend('topright', legend=c('Non-hierarchical Phase',
                              'Hierarchical Phase',
                              'Non-hierarchical ICAR',
                              'Hierarchical ICAR'), 
         fill= c('darkgreen','darkorange','darkred','darkblue')) 
  title(main = paste("Area", k))
}

dev.off()


# pdf(file=here('output', 'figures','figure14W.pdf'), width=14, height=18)
# # Define the layout for the plots
# par(mfrow = c(7, 6))
# 
# for (k in 1:41) #all hex areas
# {
#   post.model.i <- do.call(rbind, mcmc.samplesW)[ , c(k, k+41)] #selecting a[k] and b[k]
# 
#   dens.i.a <- density(post.model.i[,1],bw = 5)
#   dens.i.b <- density(post.model.i[,2],bw=5)
#   
#   # Plot
#   plot(NULL, xlim=c(sim_a[[k]]+1000, sim_b[[k]]-1000), ylim=c(0,0.022), xlab='Cal BP', ylab='Posterior Probability')
#   polygon(c(dens.i.a$x, rev(dens.i.a$x)), c(rep(0,length(dens.i.a$x)), rev(dens.i.a$y)), border=NA, col=rgb(0,0.4,0,0.5))
#   #polygon(c(dens.i.b$x, rev(dens.i.b$x)), c(rep(0,length(dens.i.b$x)), rev(dens.i.b$y)), border=NA, col=rgb(0,0.4,0,0.5))
#   abline(v=sim_a[[k]],lty=2) #abline(v=c(sim_a[[k]], sim_b[[k]]),lty=2)
#   axis(3,at=sim_a[[k]],labels=TeX('$a$')) #axis(3,at=c(sim_a[[k]], sim_b[[k]]),labels=c(TeX('$a$'),TeX('$b$')))
#   legend('topright', legend=c('Womble ICAR'), 
#          fill= c('darkgreen')) 
#   title(main = paste("Area", k))
# }
# 
# dev.off()

#-------------
##Show alpha values are recovered at sites for hierarchical ICAR model -- FIGURE 14.2
sim_alpha <- constants$true_alpha

pdf(file=here('output', 'figures','figure14.2.pdf'), width=14, height=18)
# Define the layout for the plots
par(mfrow = c(10, 10))

for (j in 1:100){ #all sites
  post.model.iii <- do.call(rbind, mcmc.samples3)[ , c(j+41)]
  dens.iii.a <- density(post.model.iii,bw = 5)

  # Plot
  plot(NULL, xlim=c(sim_alpha[[j]]+300, sim_alpha[[j]]-300), ylim=c(0,0.06), xlab='Cal BP', ylab='Posterior Probability')
  polygon(c(dens.iii.a$x, rev(dens.iii.a$x)), c(rep(0,length(dens.iii.a$x)), rev(dens.iii.a$y)), border=NA, col=rgb(1,0.55,0,0.5))
  abline(v=sim_alpha[[j]],lty=2)
  axis(3,at=sim_alpha[[j]],labels=TeX('$alpha$'))
  title(main = paste("Site", j))
}

dev.off()

#-------------------------------------------------------------------------------
## Tactical Simulation Posterior Predictive Check for gradient in a given region -- FIGURE 15
## For hierarchical ICAR model

#For model (i) and (ii) select parameters a and b (i.e. start and end date of occupation in the region)
pdf(file=here('output', 'figures','figure15.pdf'), width=14, height=18)
# Define the layout for the plots
par(mfrow = c(6, 6))

for (k in c(5, 7:41)) #all hex areas
{
  post.grad.model.iii <- do.call(rbind, mcmc.samples3)[ , k+282] #selecting nabla[k] #selecting nabla[k] #e.g.to find nabla[18] index: which(colnames(as.data.frame(mcmc.samples3$chain1)) == 'nabla[18]')
  dens.iii.nabla  <- density(post.grad.model.iii, bw = 5)

  # Plot
  plot(NULL, xlim=c(-20,20), ylim=c(0,0.1), xlab='Gradient', ylab='Posterior Probability')
  polygon(c(dens.iii.nabla$x, rev(dens.iii.nabla$x)), c(rep(0,length(dens.iii.nabla$x)), rev(dens.iii.nabla$y)), border=NA, col=rgb(0,0.4,0,0.5))
  legend('topright', legend=c('Hierarchichal ICAR'), fill=c('darkgreen'))
  title(main = paste("Area", k))
}

dev.off()


#-------------------------------------------------------------------------------
##Plot magnitude and direction of gradients for model (i) -- FIGURE 16
## For hierarchical ICAR model

#Extract quantile information for models
tmp = extract_gradinfo(mcmc.samples3) 
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

pdf(file=here('output', 'figures','figure16.pdf'), width=8, height=8)
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
plot_arrows(edge_info, lwd=4, scale_par=8, length = 0.1) #length parameter defines arrowhead size
dev.off()

#-------------------------------------------------------------------------------
##Display proportion of distribution before specified time slices -- FIGURE 25

#With the tactical simulation data from hierarchical ICAR model
out.comb.tac_icar.model  <- do.call(rbind, mcmc.samples3)
post.model.tac_icar  <- out.comb.tac_icar.model[,paste0('a[',1:41,']')]  %>% round()

#Extract arrival times for tactical icar model
med.model.tac_icar  <- apply(post.model.tac_icar, 2, median)

time_slices <- seq(3200, 1400, -200)

#Extract proportion of MCMC samples occurring before the specified time threshold
prop_model_tac_icar  <- lapply(time_slices,
                               function(t) data.frame(x = 1:41,
                                                      y = sapply(as.data.frame(post.model.tac_icar), 
                                                                 prop_gthan_threshold, 
                                                                 threshold = t)))

#Save figure for each time slice
plot_list <- list() #Create a list to store the plots

for (k in 1:length(time_slices)) #all time slices
{
  #Save in data structure
  prop_threshold_tac_icar <- hex_area_win %>%
    filter(area_ID %in% 1:41) %>%
    mutate(median_date = med.model.tac_icar,
           prop_threshold = prop_model_tac_icar[[k]]$y,
           contains_sites = as.factor(case_when(area_ID %in% Hex_with_sites ~ 1, area_ID %in% Hex_without_sites ~ 0))) %>%
    filter(area_ID %!in% c(1, 2, 3, 4, 6)) #The Bantu hadn't settled in this area by the time the dutch arrived in the Cape (~1600AD). To back this up there are no EIA sites in these regions.
  
  #Plot
  p <- ggplot(data = prop_threshold_tac_icar) +
    geom_sf(data = st_buffer(st_as_sf(sampling_win, crs = 4326), 40000), aes(color = "grey50")) + #sampling window with coastal buffer
    geom_sf(aes(fill = median_date, alpha=prop_threshold)) + #hex grid #alpha=contains_sites
    scale_fill_viridis_c(option="F", direction=-1) +
    scale_alpha_continuous(range = c(0, 1)) +  # Use for continuous alpha values
    xlab('Longitude') +
    ylab('Latitude') +
    ggtitle(paste0('t = ', time_slices[k], ' BP')) +
    #geom_sf_label(aes(label = paste0(median_date, "BP")), label.size  = NA, alpha = 0.4, size=3.5) + #hex grid labels #label = ifelse(contains_sites==0, NA, paste0(median_date, "BP")))
    theme(panel.background = element_rect(fill = "lightblue",
                                          colour = "lightblue",
                                          size = 0.5,
                                          linetype = "solid"),
          legend.position = "none")
  
  plot_list[[k]] <- p
}

#Output
pdf(file=here('output','figures','figure25.pdf'), width=15, height=8)
grid.arrange(grobs=plot_list, ncol=5, nrow=2, padding=0) # Define the layout for the plots
dev.off()

#-------------------------------------------------------------------------------
##Animate through time slices to display proportion of distribution before specified time -- FIGURE 25_animate

#Extract proportion of MCMC samples occurring before the specified time threshold (note more time divisions for animation)
time_slices2 <- seq(4000, 1000, -50) 

prop_model_tac_icar2  <- lapply(time_slices2, 
                            function(t) data.frame(x = 1:41, 
                                                   y = sapply(as.data.frame(post.model.tac_icar), 
                                                            prop_gthan_threshold, 
                                                            threshold = t)))

# Combine all time slices into one dataset
all_data <- lapply(seq_along(time_slices2), function(k) {
  hex_area_win %>%
    filter(area_ID %in% 1:41) %>%
    mutate(median_date = med.model.tac_icar,
           prop_threshold = prop_model_tac_icar2[[k]]$y,
           time_slice = -time_slices2[k]) %>% #minus sign ensures animation time proceeds in the correct direction
    filter(area_ID %!in% c(1, 2, 3, 4, 6)) %>% # Exclude specific areas
    st_as_sf()})

# Verify that all elements are sf objects
if (!all(sapply(all_data, inherits, "sf"))) {
  stop("One or more elements of all_data are not sf objects!")
}

# Combine all sf objects safely
all_data <- do.call(rbind, all_data)

#co-ordinates of area centroids used for label positions
label_data <- st_coordinates(all_data$area_center)

# Animated ggplot
p <- ggplot(data = all_data) +
  geom_sf(data = st_buffer(st_as_sf(sampling_win, crs = 4326), 40000), aes(color = "grey50")) +
  geom_sf(aes(fill = median_date, alpha = prop_threshold)) +
  scale_fill_viridis_c(option = "F", direction = -1) +
  scale_alpha_continuous(range = c(0, 1)) +
  xlab('Longitude') +
  ylab('Latitude') +
  geom_text(data = label_data,
            aes(x = X, y = Y, label = paste0(all_data$median_date, " BP")),
            size = 4.5, alpha = 0.4) +  #Use geom_text instead of geom_sf_label
  theme(
    panel.background = element_rect(
      fill = "lightblue",
      colour = "lightblue",
      size = 0.5,
      linetype = "solid"
    ),
    legend.position = "none",
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 15, face = "bold"),  # X-axis label
    axis.title.y = element_text(size = 15, face = "bold")   # Y-axis label
  ) +
  labs(title = "Time Slice: {closest_state} BP") #Note: use '{closest_state}' with 'transition_states()' and '{frame_time}' with 'transition_time()'

anim <- p +
  transition_states(time_slice, transition_length = 2, state_length = 1) +
  enter_fade() +
  exit_fade()

#Save the animation
anim_save("figure25_animation.gif",
          animation = animate(anim, width = 1200, height = 900, fps = 10))

#===============================================================================
##Bayesian ICAR Models

`#Load Data ----
load(here("output", "ICAR_model_a.RData")) #model (i) -- no sample interdependence
load(here("output","ICAR_model_b.RData")) #model (ii) -- hierarchical structure
`
#-------------------------------------------------------------------------------
# Marginal Posterior Distribution of a[k], model i ---- FIGURE 17

out.comb.icar.model.a  <- do.call(rbind, out_icar_model_a)
post.a.model.i  <- out.comb.icar.model.a[,paste0('a[',c(5,7:41),']')] %>%  round() #all relevant hex areas
model.i.long  <- data.frame(value = as.numeric(post.a.model.i),
                            area = rep(c(5,7:41), each=nrow(post.a.model.i)))

model.i.long  <- model.i.long %>%
  mutate(area = factor(area, levels=paste0(c(5,7:41)), ordered=TRUE))

#Plot
pdf(file=here('output','figures','figure17.pdf'), height=10, width=8)
ggplot(model.i.long, aes(x = value, y = area, fill='orange')) + 
  geom_density_ridges() +
  scale_x_reverse(limits=c(3470, 1000), breaks=BCADtoBP(c(-1500, -1300, -1100, -900, -700, -500, -300, -100, 100, 300, 500, 700, 900)), labels=c('1500BC','1300BC', '1100BC','900BC','700BC', '500BC', '300BC', '100BC', '100AD', '300AD', '500AD', '700AD', '900AD')) +
  scale_fill_manual(values='orange') +
  xlab(paste('Arrival time,', TeX('$a_k$'))) +
  ylab(paste('Area,', TeX('$k$'))) +
  theme(legend.position = "none")
dev.off()

#------
# Marginal Posterior Distribution of a[k], model ii ---- FIGURE 18

out.comb.icar.model.b  <- do.call(rbind, out_icar_model_b)
post.a.model.ii  <- out.comb.icar.model.b[,paste0('a[',c(5,7:41),']')] %>%  round() #all relevant hex areas
model.ii.long  <- data.frame(value = as.numeric(post.a.model.ii),
                            area = rep(c(5,7:41), each=nrow(post.a.model.ii)))

model.ii.long  <- model.ii.long %>%
  mutate(area = factor(area, levels=paste0(c(5,7:41)), ordered=TRUE))

#Plot
pdf(file=here('output','figures','figure18.pdf'), height=10, width=8)
ggplot(model.ii.long, aes(x = value, y = area, fill='orange')) + 
  geom_density_ridges() +
  scale_x_reverse(limits=c(3470, 1000), breaks=BCADtoBP(c(-1500, -1300, -1100, -900, -700, -500, -300, -100, 100, 300, 500, 700, 900)), labels=c('1500BC','1300BC', '1100BC','900BC','700BC', '500BC', '300BC', '100BC', '100AD', '300AD', '500AD', '700AD', '900AD')) +
  scale_fill_manual(values='orange') +
  xlab(paste('Arrival time,', TeX('$a_k$'))) +
  ylab(paste('Area,', TeX('$k$'))) +
  theme(legend.position = "none")
dev.off()

#-------------------------------------------------------------------------------
# Probability Matrix of a[k], model i ---- FIGURE 19
source(here('src','orderPPlot.R'))

pdf(file=here('output','figures','figure19.pdf'), width=10, height=10.5)
orderPPlot(post.a.model.i, name.vec=paste("Area", c(5,7:41)))
dev.off()

#------
# Probability Matrix of a[k], model ii ---- FIGURE 20
pdf(file=here('output','figures','figure20.pdf'), width=10, height=10.5)
orderPPlot(post.a.model.ii, name.vec=paste("Area", c(5,7:41)))
dev.off()

#-------------------------------------------------------------------------------
# Estimated Arrival Date ---- FIGURE 21
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
  rect(xleft=x[2], xright=x[6], ybottom=i-h/5, ytop=i+h/5, border=NA, col=col)
  rect(xleft=x[3], xright=x[5], ybottom=i-h/3, ytop=i+h/3, border=NA, col=col)
  lines(c(x[4],x[4]), c(i-h/2,i+h/2), lwd=2, col='grey44')
}

# Posterior Arrival Times
pdf(file=here('output','figures','figure21.pdf'), width=10, height=18, pointsize=4)
plot(NULL, xlim=c(3470, 1000), ylim=c(3,105), xlab=paste('Arrival time,', TeX('$a_k$')), ylab=paste('Area,', TeX('$k$')), cex.lab = 2, axes=F)
tmp.a = extract(out_icar_model_a)
tmp.b = extract(out_icar_model_b)
iseq.a = seq(2,by=3,length.out=36)
iseq.b = seq(1,by=3,length.out=36)
abline(h=seq(3,by=3,length.out=36), col='darkgrey',lty=2)

counter <- 1 #indexing counter
for (i in c(5,7:41)) #all relevant hex areas
{
  #Plot bar in area i
  post.bar(tmp.a[,i], i=iseq.a[counter], h=0.9, col='lightblue')
  post.bar(tmp.b[,i], i=iseq.b[counter], h=0.9, col='lightgreen')
  counter <- counter + 1
}

axis(2, at=iseq.a+0.5, labels = paste0(c(5,7:41)), las=2, cex.axis=1.7)
axis(1, at = BCADtoBP(c(-1500, -1300, -1100, -900, -700, -500, -300, -100, 100, 300, 500, 700, 900)), labels=c('1500BC','1300BC', '1100BC','900BC','700BC', '500BC', '300BC', '100BC', '100AD', '300AD', '500AD', '700AD', '900AD'), tck=-0.01, cex.axis=1.7)
axis(3, at = seq(3400, 1000, -200), labels=paste0(seq(3400, 1000, -200),'BP'), tck=-0.01, cex.axis=1.7)
axis(1, at = BCADtoBP(c(-600, -400, -200, 1, 200, 400, 600)), labels=NA, tck=-0.01) #Minor tick marks
axis(3, at = seq(3400, 1000, -100), labels=NA, tck=-0.01) #Minor tick marks
box()

# post.bar(c(800,700,600,500,400,300,200), i=1.5, h=0.9, col='lightgrey')
# arrows(x0=700, x1=300, y0=0.3, y1=0.3, angle = 90, code = 3, length = 0.01)
# arrows(x0=600, x1=400, y0=0.8, y1=0.8, angle = 90, code = 3, length = 0.01)
# text(x=1085, y=0.9, "50% HPDI", cex=1.5)
# text(x=1275, y=0.07,"90% HPDI", cex=1.5)
# text(x=1100, y=2.3, "Median Posterior", cex=1.5)
# lines(x=c(550,500), y=c(2.35, 2.28))
# text(x=4000, y=2, 'Model A', cex=1.5)
# text(x=4000, y=1, 'Model B', cex=1.5)
# theme(legend.position = "none")
dev.off()

#-------------------------------------------------------------------------------
## Plot HEX areas with median arrival times ---- FIGURE 22

#Extract arrival times for model (i)
out.comb.icar.modela  <- do.call(rbind, out_icar_model_a)
post.a.model.i  <- out.comb.icar.modela[,paste0('a[',1:41,']')]  %>% round() 
hpdi.model.i  <- apply(post.a.model.i, 2, function(x){HPDinterval(as.mcmc(x), prob = .90)}) 
med.model.i  <- apply(post.a.model.i, 2, median)
hi90_mod.i  <- hpdi.model.i[1,]
lo90_mod.i  <- hpdi.model.i[2,]

median_hex_dates_mod.i <- hex_area_win %>% 
  filter(area_ID %in% 1:41) %>% 
  mutate(median_date = med.model.i,
         hpdi_high = hi90_mod.i,
         hpdi_low = lo90_mod.i,
         contains_sites = as.factor(case_when(area_ID %in% Hex_with_sites ~ 1, area_ID %in% Hex_without_sites ~ 0))) %>% 
  filter(area_ID %!in% c(1, 2, 3, 4, 6)) #The Bantu hadn't settled in this area by the time the dutch arrived in the Cape (~1600AD). To back this up there are no EIA sites in these regions.

#Extract arrival times for model (ii)
out.comb.icar.modelb  <- do.call(rbind, out_icar_model_b)
post.a.model.ii  <- out.comb.icar.modelb[,paste0('a[',1:41,']')]  %>% round()
hpdi.model.ii  <- apply(post.a.model.ii, 2, function(x){HPDinterval(as.mcmc(x), prob = .90)}) 
med.model.ii  <- apply(post.a.model.ii, 2, median)
hi90_mod.ii  <- hpdi.model.ii[1,]
lo90_mod.ii  <- hpdi.model.ii[2,]

median_hex_dates_mod.ii <- hex_area_win %>% 
  filter(area_ID %in% 1:41) %>% 
  mutate(median_date = med.model.ii,
         hpdi_high = hi90_mod.ii,
         hpdi_low = lo90_mod.ii,
         contains_sites = as.factor(case_when(area_ID %in% Hex_with_sites ~ 1, area_ID %in% Hex_without_sites ~ 0))) %>% 
  filter(area_ID %!in% c(1, 2, 3, 4, 6)) #The Bantu hadn't settled in this area by the time the dutch arrived in the Cape (~1600AD). To back this up there are no EIA sites in these regions.

#Plot
#-----MODEL (i)
modi <- ggplot(data = median_hex_dates_mod.i) +
  geom_sf(data = st_buffer(st_as_sf(sampling_win, crs = 4326), 40000), aes(color = "grey50")) + #sampling window with coastal buffer
  geom_sf(aes(fill = median_date)) + #hex grid #alpha=contains_sites
  scale_fill_viridis_c(option="F", direction=-1) +
  scale_alpha_manual(values=c(0.45, 1)) +
  xlab('Longitude') +
  ylab('Latitude') +
  geom_sf_label(aes(label = paste0(median_date, "BP")), label.size  = NA, alpha = 0.4, size=3.5) + #hex grid labels #label = ifelse(contains_sites==0, NA, paste0(median_date, "BP")))
  theme(panel.background = element_rect(fill = "lightblue",
                                        colour = "lightblue",
                                        size = 0.5,
                                        linetype = "solid"),
        legend.position = "none")


modiHPDIlow <- ggplot(data = median_hex_dates_mod.i) +
  geom_sf(data = st_buffer(st_as_sf(sampling_win, crs = 4326), 40000), aes(color = "grey50")) + #sampling window with coastal buffer
  geom_sf(aes(fill = hpdi_low)) + #alpha=contains_sites
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

modiHPDIhigh <- ggplot(data = median_hex_dates_mod.i) +
  geom_sf(data = st_buffer(st_as_sf(sampling_win, crs = 4326), 40000), aes(color = "grey50")) + #sampling window with coastal buffer
  geom_sf(aes(fill = hpdi_high)) + #alpha=contains_sites
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
  draw_plot(modi) +
  draw_plot(modiHPDIlow, 
            x = .73, y = .285, width = .25, height = .25) +
  draw_plot(modiHPDIhigh, 
            x = .73, y = .06, width = .25, height = .25)

#-----MODEL (ii)
modii <- ggplot(data = median_hex_dates_mod.ii) +
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

modiiHPDIlow <- ggplot(data = median_hex_dates_mod.ii) +
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

modiiHPDIhigh <- ggplot(data = median_hex_dates_mod.ii) +
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
  draw_plot(modii) +
  draw_plot(modiiHPDIlow, 
            x = .73, y = .285, width = .25, height = .25) +
  draw_plot(modiiHPDIhigh, 
            x = .73, y = .06, width = .25, height = .25)

#Output
pdf(file=here('output','figures','figure22.pdf'), width=15, height=8)
grid.arrange(A, B, ncol=2, padding=0)
dev.off()


#-------------------------------------------------------------------------------
##Plot magnitude and direction of gradients for model (i) and (ii) -- FIGURE 23

#Extract quantile information for models
tmp.i = extract_gradinfo(out_icar_model_a) 
tmp.ii = extract_gradinfo(out_icar_model_b) 

qta.i = apply(tmp.i, 2, quantile, prob=c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))
qta.ii = apply(tmp.ii, 2, quantile, prob=c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))

#Extract uncertainty information
uncert.i = sapply(as.data.frame(tmp.i), prop_gthan_zero)
uncert.ii = sapply(as.data.frame(tmp.ii), prop_gthan_zero)

#Add info to edges dataframe
edge_info.i <- edge_info %>% 
  mutate(mean_gradient = qta.i[4,], #50% quantile
         uncertainty = uncert.i) %>% #% of distribution > zero
  filter(distance <= mean(distance)*1.3, #Remove extremely long edges which are 'artificially' created along the internal window boundary or require extreme coastal movement along external window boundary
         region1_id %!in% c(1,2,3,4,6),
         region2_id %!in% c(1,2,3,4,6))

edge_info.ii <- edge_info %>% 
  mutate(mean_gradient = qta.ii[4,], #50% quantile
         uncertainty = uncert.ii) %>% #% of distribution > zero
  filter(distance <= mean(distance)*1.3, #Remove extremely long edges which are 'artificially' created along the internal window boundary or require extreme coastal movement along external window boundary
         region1_id %!in% c(1,2,3,4,6),
         region2_id %!in% c(1,2,3,4,6))

#Create nodes
rel_hex_win <- hex_area_win %>% filter(area_ID %!in% c(1,2,3,4,6))
nodes <- st_coordinates(rel_hex_win$area_center)

#--------
#PLOT

plot_grad <- function(edges_info, scale_par, sampling_window){
  # Get the bounding box of the sample window
  sample_win_buff <- st_buffer(st_as_sf(sampling_window, crs = 4326), 40000)
  bbox <- sf::st_bbox(sample_win_buff)
  
  # Calculate the aspect ratio of the bounding box
  aspect_ratio <- diff(range(c(bbox["ymin"], bbox["ymax"]))) / diff(range(c(bbox["xmin"], bbox["xmax"])))
  gridsize <- 1 # Adjust the denominator to change grid density
  
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
  plot_arrows(edges_info, scale_par, lwd=4, length = 0.1) #length parameter defines arrowhead size
}


#Output
pdf(file=here('output','figures','figure23.pdf'))
A.i <- plot_grad(edge_info.i, 3, sampling_win)
B.ii <- plot_grad(edge_info.ii, 3, sampling_win)
grid.arrange(A.i, B.ii, ncol=2, padding=0)
dev.off()

#-------------------------------------------------------------------------------
# Traceplot of start and end of occupation (a, b) ---- FIGURE 24

pdf(file=here('output', 'figures','figure24.pdf'), width=8, height=8)
par(mfrow=c(2,2))
traceplot(out_icar_model_a[,'a[18]'], main=TeX('$a[18]$'),smooth=TRUE) #region area 18, as an example
traceplot(out_icar_model_a[,'b[18]'], main=TeX('$b[18]$'),smooth=TRUE)
traceplot(out_icar_model_a[,'nabla[18]'], main=TeX('$nabla[18]$'), smooth=TRUE)
dev.off()

pdf(file=here('output', 'figures','figure26.pdf'), width=8, height=8)
par(mfrow=c(2,2))
traceplot(out_icar_model_b[,'a[18]'], main=TeX('$a[18]$'),smooth=TRUE) #region area 18, as an example
traceplot(out_icar_model_b[,'b[18]'], main=TeX('$b[18]$'),smooth=TRUE)
traceplot(out_icar_model_b[,'nabla[18]'], main=TeX('$nabla[18]$'), smooth=TRUE)
dev.off()

#-------------------------------------------------------------------------------
##Display proportion of distribution before specified time slices -- FIGURE 27

#Hierarchical ICAR model
out.comb.icar.modelb  <- do.call(rbind, out_icar_model_b)
post.modelb.icar  <- out.comb.icar.modelb[,paste0('a[',1:41,']')]  %>% round()

#Extract arrival times for tactical icar model
med.modelb.icar  <- apply(post.modelb.icar, 2, median)

time_slices <- seq(3200, 1400, -200)

#Extract proportion of MCMC samples occurring before the specified time threshold
prop_modelb_icar  <- lapply(time_slices,
                            function(t) data.frame(x = 1:41,
                                                   y = sapply(as.data.frame(post.modelb.icar), 
                                                              prop_gthan_threshold, 
                                                              threshold = t)))

#Save figure for each time slice
plot_list <- list() #Create a list to store the plots

for (k in 1:length(time_slices)) #all time slices
{
  #Save in data structure
  prop_thresholdb_icar <- hex_area_win %>%
    filter(area_ID %in% 1:41) %>%
    mutate(median_date = med.modelb.icar,
           prop_threshold = prop_modelb_icar[[k]]$y,
           contains_sites = as.factor(case_when(area_ID %in% Hex_with_sites ~ 1, area_ID %in% Hex_without_sites ~ 0))) %>%
    filter(area_ID %!in% c(1, 2, 3, 4, 6)) #The Bantu hadn't settled in this area by the time the dutch arrived in the Cape (~1600AD). To back this up there are no EIA sites in these regions.
  
  #Plot
  p <- ggplot(data = prop_thresholdb_icar) +
    geom_sf(data = st_buffer(st_as_sf(sampling_win, crs = 4326), 40000), aes(colour = "grey50"), color=NA) + #sampling window with coastal buffer
    geom_sf(aes(fill = median_date, alpha=prop_threshold)) + #hex grid #alpha=contains_sites
    scale_fill_viridis_c(option="F", direction=-1) +
    scale_alpha_continuous(range = c(0, 1)) +  # Use for continuous alpha values
    xlab('Longitude') +
    ylab('Latitude') +
    ggtitle(paste0('t = ', time_slices[k], ' BP')) +
    #geom_sf_label(aes(label = paste0(median_date, "BP")), label.size  = NA, alpha = 0.4, size=3.5) + #hex grid labels #label = ifelse(contains_sites==0, NA, paste0(median_date, "BP")))
    theme(panel.background = element_rect(fill = "lightblue",
                                          colour = "lightblue",
                                          size = 0.5,
                                          linetype = "solid"),
          legend.position = "none")
  
  plot_list[[k]] <- p
}

#Output
pdf(file=here('output','figures','figure27.pdf'), width=15, height=8)
grid.arrange(grobs=plot_list, ncol=5, nrow=2, padding=0) # Define the layout for the plots
dev.off()

#-------------------------------------------------------------------------------
##Animate through time slices to display proportion of distribution before specified time -- FIGURE 27_animate

#Extract proportion of MCMC samples occurring before the specified time threshold (note more time divisions for animation)
time_slices2 <- seq(4000, 1000, -50) 

prop_modelb_icar2  <- lapply(time_slices2, 
                             function(t) data.frame(x = 1:41, 
                                                    y = sapply(as.data.frame(post.modelb.icar), 
                                                               prop_gthan_threshold, 
                                                               threshold = t)))

# Combine all time slices into one dataset
all_data <- lapply(seq_along(time_slices2), function(k) {
  hex_area_win %>%
    filter(area_ID %in% 1:41) %>%
    mutate(median_date = med.modelb.icar,
           prop_threshold = prop_modelb_icar2[[k]]$y,
           time_slice = -time_slices2[k]) %>% #minus sign ensures animation time proceeds in the correct direction
    filter(area_ID %!in% c(1, 2, 3, 4, 6)) %>% # Exclude specific areas
    st_as_sf()})

# Verify that all elements are sf objects
if (!all(sapply(all_data, inherits, "sf"))) {
  stop("One or more elements of all_data are not sf objects!")
}

# Combine all sf objects safely
all_data <- do.call(rbind, all_data)

#co-ordinates of area centroids used for label positions
label_data <- st_coordinates(all_data$area_center)

# Animated ggplot
p <- ggplot(data = all_data) +
  geom_sf(data = st_buffer(st_as_sf(sampling_win, crs = 4326), 40000), aes(color = "grey50"), color=NA) +
  geom_sf(aes(fill = median_date, alpha = prop_threshold)) +
  scale_fill_viridis_c(option = "F", direction = -1) +
  scale_alpha_continuous(range = c(0, 1)) +
  xlab('Longitude') +
  ylab('Latitude') +
  geom_text(data = label_data,
            aes(x = X, y = Y, label = paste0(all_data$median_date, " BP")),
            size = 4.5, alpha = 0.4) +  #Use geom_text instead of geom_sf_label
  theme(
    panel.background = element_rect(
      fill = "lightblue",
      colour = "lightblue",
      size = 0.5,
      linetype = "solid"
    ),
    legend.position = "none",
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 15, face = "bold"),  # X-axis label
    axis.title.y = element_text(size = 15, face = "bold")   # Y-axis label
  ) +
  labs(title = "Time Slice: {closest_state} BP") #Note: use '{closest_state}' with 'transition_states()' and '{frame_time}' with 'transition_time()'

anim <- p +
  transition_states(time_slice, transition_length = 2, state_length = 1) +
  enter_fade() +
  exit_fade()

#Save the animation
anim_save("figure27_animation.gif",
          animation = animate(anim, width = 1200, height = 900, fps = 10))


# #-------------------------------------------------------------------------------
# # Marginal posteriors of beta0, beta1, rho, etasq for tau = 0.9
# 
# gpqr.tau90.comb  <- do.call(rbind, gpqr_tau90)
# 
# pdf(file=here('output','figures','figure14.pdf'), width=8, height=8)
# par(mfrow=c(2,2))
# postHPDplot(gpqr.tau90.comb[,'beta0'], main=TeX('$\\beta_0$'), xlab='Cal BP', ylab='')
# postHPDplot(gpqr.tau90.comb[,'beta1'], main=TeX('$\\beta_1$'), xlab='', ylab='')
# postHPDplot(gpqr.tau90.comb[,'rho'], main=TeX('$\\rho$'), xlab='km', ylab='')
# postHPDplot(gpqr.tau90.comb[,'etasq'], main=TeX('$\\eta^2$'), xlab='', ylab='')
# dev.off()


#===============================================================================
##Bayesian Hierarchical ICAR Model -- over entire sub-Saharan region

## Load Data

# Load Observed Data
load(here('data','c14.RData')) #sub-Saharan Africa
#Load nodes and edges between hex area centroids
load(here('data','trig_cont.RData'))
#Load model data
load(here("output", "ICAR_model_a_cont.RData"))
load(here("output", "ICAR_model_b_cont.RData"))

#----

#Hex areas with and without out sites
Hex_with_sites <- unique(siteInfo$area_id)
Hex_without_sites <- which(rep(1:47) %!in% Hex_with_sites)


#-------------------------------------------------------------------------------
# Marginal Posterior Distribution of a[k], model i ---- FIGURE 28

out.comb.icar.model.a  <- do.call(rbind, out_icar_model_a)
post.a.model.i  <- out.comb.icar.model.a[,paste0('a[',c(1:47),']')] %>%  round() #all relevant hex areas
model.i.long  <- data.frame(value = as.numeric(post.a.model.i),
                            area = rep(c(1:47), each=nrow(post.a.model.i)))

model.i.long  <- model.i.long %>%
  mutate(area = factor(area, levels=paste0(c(1:47)), ordered=TRUE))

#Plot
pdf(file=here('output','figures','figure28.pdf'), height=10, width=8)
ggplot(model.i.long, aes(x = value, y = area, fill='orange')) + 
  geom_density_ridges() +
  scale_x_reverse(limits=c(4470, 1000), breaks=BCADtoBP(c(-2300, -1900, -1500, -1100, -700, -300, 100, 500, 900)), labels=c('2300BC', '1900BC', '1500BC', '1100BC','700BC',  '300BC',  '100AD', '500AD', '900AD')) +
  scale_fill_manual(values='orange') +
  xlab(paste('Arrival time,', TeX('$a_k$'))) +
  ylab(paste('Area,', TeX('$k$'))) +
  theme(legend.position = "none")
dev.off()

#------
# Marginal Posterior Distribution of a[k], model ii ---- FIGURE 29

out.comb.icar.model.b  <- do.call(rbind, out_icar_model_b)
post.a.model.ii  <- out.comb.icar.model.b[,paste0('a[',c(1:47),']')] %>%  round() #all relevant hex areas
model.ii.long  <- data.frame(value = as.numeric(post.a.model.ii),
                             area = rep(c(1:47), each=nrow(post.a.model.ii)))

model.ii.long  <- model.ii.long %>%
  mutate(area = factor(area, levels=paste0(c(1:47)), ordered=TRUE))

#Plot
pdf(file=here('output','figures','figure29.pdf'), height=10, width=8)
ggplot(model.ii.long, aes(x = value, y = area, fill='orange')) + 
  geom_density_ridges() +
  scale_x_reverse(limits=c(4470, 1000), breaks=BCADtoBP(c(-2300, -1900, -1500, -1100, -700, -300, 100, 500, 900)), labels=c('2300BC', '1900BC', '1500BC', '1100BC','700BC',  '300BC',  '100AD', '500AD', '900AD')) +
  scale_fill_manual(values='orange') +
  xlab(paste('Arrival time,', TeX('$a_k$'))) +
  ylab(paste('Area,', TeX('$k$'))) +
  theme(legend.position = "none")
dev.off()

#-------------------------------------------------------------------------------
# Probability Matrix of a[k], model i ---- FIGURE 30
source(here('src','orderPPlot.R'))

pdf(file=here('output','figures','figure30.pdf'), width=10, height=10.5)
orderPPlot(post.a.model.i, name.vec=paste("Area", c(1:47)))
dev.off()

#------
# Probability Matrix of a[k], model ii ---- FIGURE 31
pdf(file=here('output','figures','figure31.pdf'), width=10, height=10.5)
orderPPlot(post.a.model.ii, name.vec=paste("Area", c(1:47)))
dev.off()

#-------------------------------------------------------------------------------
## Plot HEX areas with median arrival times ---- FIGURE 32

#Extract arrival times for model (i)
out.comb.icar.modela  <- do.call(rbind, out_icar_model_a)
post.a.model.i  <- out.comb.icar.modela[,paste0('a[',1:47,']')]  %>% round() 
hpdi.model.i  <- apply(post.a.model.i, 2, function(x){HPDinterval(as.mcmc(x), prob = .90)}) 
med.model.i  <- apply(post.a.model.i, 2, median)
hi90_mod.i  <- hpdi.model.i[1,]
lo90_mod.i  <- hpdi.model.i[2,]

median_hex_dates_mod.i <- hex_area_win %>% 
  filter(area_ID %in% 1:47) %>% 
  mutate(median_date = med.model.i,
         hpdi_high = hi90_mod.i,
         hpdi_low = lo90_mod.i,
         contains_sites = as.factor(case_when(area_ID %in% Hex_with_sites ~ 1, area_ID %in% Hex_without_sites ~ 0))) 

#Extract arrival times for model (ii)
out.comb.icar.modelb  <- do.call(rbind, out_icar_model_b)
post.a.model.ii  <- out.comb.icar.modelb[,paste0('a[',1:47,']')]  %>% round()
hpdi.model.ii  <- apply(post.a.model.ii, 2, function(x){HPDinterval(as.mcmc(x), prob = .90)}) 
med.model.ii  <- apply(post.a.model.ii, 2, median)
hi90_mod.ii  <- hpdi.model.ii[1,]
lo90_mod.ii  <- hpdi.model.ii[2,]

median_hex_dates_mod.ii <- hex_area_win %>% 
  filter(area_ID %in% 1:47) %>% 
  mutate(median_date = med.model.ii,
         hpdi_high = hi90_mod.ii,
         hpdi_low = lo90_mod.ii,
         contains_sites = as.factor(case_when(area_ID %in% Hex_with_sites ~ 1, area_ID %in% Hex_without_sites ~ 0))) 

#Plot
#-----MODEL (i)
modi <- ggplot(data = median_hex_dates_mod.i) +
  geom_sf(data = st_buffer(st_as_sf(sampling_win, crs = 4326), 40000), aes(color = "grey50")) + #sampling window with coastal buffer
  geom_sf(aes(fill = median_date)) + #hex grid #alpha=contains_sites
  scale_fill_viridis_c(option="F", direction=-1) +
  scale_alpha_manual(values=c(0.45, 1)) +
  xlab('Longitude') +
  ylab('Latitude') +
  geom_sf_label(aes(label = paste0(median_date, "BP")), label.size  = NA, alpha = 0.4, size=3.5) + #hex grid labels #label = ifelse(contains_sites==0, NA, paste0(median_date, "BP")))
  theme(panel.background = element_rect(fill = "lightblue",
                                        colour = "lightblue",
                                        size = 0.5,
                                        linetype = "solid"),
        legend.position = "none")


modiHPDIlow <- ggplot(data = median_hex_dates_mod.i) +
  geom_sf(data = st_buffer(st_as_sf(sampling_win, crs = 4326), 40000), aes(color = "grey50")) + #sampling window with coastal buffer
  geom_sf(aes(fill = hpdi_low)) + #alpha=contains_sites
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

modiHPDIhigh <- ggplot(data = median_hex_dates_mod.i) +
  geom_sf(data = st_buffer(st_as_sf(sampling_win, crs = 4326), 40000), aes(color = "grey50")) + #sampling window with coastal buffer
  geom_sf(aes(fill = hpdi_high)) + #alpha=contains_sites
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
  draw_plot(modi) +
  draw_plot(modiHPDIlow, 
            x = .73, y = .285, width = .25, height = .25) +
  draw_plot(modiHPDIhigh, 
            x = .73, y = .06, width = .25, height = .25)

#-----MODEL (ii)
modii <- ggplot(data = median_hex_dates_mod.ii) +
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

modiiHPDIlow <- ggplot(data = median_hex_dates_mod.ii) +
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

modiiHPDIhigh <- ggplot(data = median_hex_dates_mod.ii) +
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
  draw_plot(modii) +
  draw_plot(modiiHPDIlow, 
            x = .73, y = .285, width = .25, height = .25) +
  draw_plot(modiiHPDIhigh, 
            x = .73, y = .06, width = .25, height = .25)

#Output
pdf(file=here('output','figures','figure32.pdf'), width=15, height=8)
grid.arrange(A, B, ncol=2, padding=0)
dev.off()


#-------------------------------------------------------------------------------
##Plot magnitude and direction of gradients for model (i) and (ii) -- FIGURE 33

#Extract quantile information for models
tmp.i = extract_gradinfo(out_icar_model_a) 
tmp.ii = extract_gradinfo(out_icar_model_b) 

qta.i = apply(tmp.i, 2, quantile, prob=c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))
qta.ii = apply(tmp.ii, 2, quantile, prob=c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))

#Extract uncertainty information
uncert.i = sapply(as.data.frame(tmp.i), prop_gthan_zero)
uncert.ii = sapply(as.data.frame(tmp.ii), prop_gthan_zero)

#Add info to edges dataframe
edge_info.i <- edge_info %>% 
  mutate(mean_gradient = qta.i[4,], #50% quantile
         uncertainty = uncert.i) %>% #% of distribution > zero
  filter(distance <= mean(distance)*1.3) #Remove extremely long edges which are 'artificially' created along the internal window boundary or require extreme coastal movement along external window boundary

edge_info.ii <- edge_info %>% 
  mutate(mean_gradient = qta.ii[4,], #50% quantile
         uncertainty = uncert.ii) %>% #% of distribution > zero
  filter(distance <= mean(distance)*1.3) #Remove extremely long edges which are 'artificially' created along the internal window boundary or require extreme coastal movement along external window boundary

#Create nodes
nodes <- st_coordinates(hex_area_win$area_center)

#--------
#PLOT

plot_grad <- function(edges_info, scale_par, sampling_window){
  # Get the bounding box of the sample window
  sample_win_buff <- st_buffer(st_as_sf(sampling_window, crs = 4326), 40000)
  bbox <- sf::st_bbox(sample_win_buff)
  
  # Calculate the aspect ratio of the bounding box
  aspect_ratio <- diff(range(c(bbox["ymin"], bbox["ymax"]))) / diff(range(c(bbox["xmin"], bbox["xmax"])))
  gridsize <- 1 # Adjust the denominator to change grid density
  
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
  plot_arrows(edges_info, scale_par, lwd=4, length = 0.1) #length parameter defines arrowhead size
}


#Output
pdf(file=here('output','figures','figure33.pdf'))
A.i <- plot_grad(edge_info.i, 3, sampling_win)
B.ii <- plot_grad(edge_info.ii, 3, sampling_win)
grid.arrange(A.i, B.ii, ncol=2, padding=0)
dev.off()

#-------------------------------------------------------------------------------
# Traceplot of start and end of occupation (a, b) ---- FIGURE 34

pdf(file=here('output', 'figures','figure34.pdf'), width=8, height=8)
par(mfrow=c(2,2))
traceplot(out_icar_model_a[,'a[18]'], main=TeX('$a[18]$'),smooth=TRUE) #region area 18, as an example
traceplot(out_icar_model_a[,'b[18]'], main=TeX('$b[18]$'),smooth=TRUE)
traceplot(out_icar_model_a[,'nabla[18]'], main=TeX('$nabla[18]$'), smooth=TRUE)
dev.off()

pdf(file=here('output', 'figures','figure35.pdf'), width=8, height=8)
par(mfrow=c(2,2))
traceplot(out_icar_model_b[,'a[18]'], main=TeX('$a[18]$'),smooth=TRUE) #region area 18, as an example
traceplot(out_icar_model_b[,'b[18]'], main=TeX('$b[18]$'),smooth=TRUE)
traceplot(out_icar_model_b[,'nabla[18]'], main=TeX('$nabla[18]$'), smooth=TRUE)
dev.off()

#-------------------------------------------------------------------------------
##Display proportion of distribution before specified time slices -- FIGURE 36

#Hierarchical ICAR model
out.comb.icar.modelb  <- do.call(rbind, out_icar_model_b)
post.modelb.icar  <- out.comb.icar.modelb[,paste0('a[',1:47,']')]  %>% round()

#Extract arrival times for tactical icar model
med.modelb.icar  <- apply(post.modelb.icar, 2, median)

time_slices <- seq(3200, 1400, -200)

#Extract proportion of MCMC samples occurring before the specified time threshold
prop_modelb_icar  <- lapply(time_slices,
                            function(t) data.frame(x = 1:47,
                                                   y = sapply(as.data.frame(post.modelb.icar), 
                                                              prop_gthan_threshold, 
                                                              threshold = t)))

#Save figure for each time slice
plot_list <- list() #Create a list to store the plots

for (k in 1:length(time_slices)) #all time slices
{
  #Save in data structure
  prop_thresholdb_icar <- hex_area_win %>%
    filter(area_ID %in% 1:47) %>%
    mutate(median_date = med.modelb.icar,
           prop_threshold = prop_modelb_icar[[k]]$y,
           contains_sites = as.factor(case_when(area_ID %in% Hex_with_sites ~ 1, area_ID %in% Hex_without_sites ~ 0))) 
  #Plot
  p <- ggplot(data = prop_thresholdb_icar) +
    geom_sf(data = st_buffer(st_as_sf(sampling_win, crs = 4326), 40000), aes(colour = "grey50"), color=NA) + #sampling window with coastal buffer
    geom_sf(aes(fill = median_date, alpha=prop_threshold)) + #hex grid #alpha=contains_sites
    scale_fill_viridis_c(option="F", direction=-1) +
    scale_alpha_continuous(range = c(0, 1)) +  # Use for continuous alpha values
    xlab('Longitude') +
    ylab('Latitude') +
    ggtitle(paste0('t = ', time_slices[k], ' BP')) +
    #geom_sf_label(aes(label = paste0(median_date, "BP")), label.size  = NA, alpha = 0.4, size=3.5) + #hex grid labels #label = ifelse(contains_sites==0, NA, paste0(median_date, "BP")))
    theme(panel.background = element_rect(fill = "lightblue",
                                          colour = "lightblue",
                                          size = 0.5,
                                          linetype = "solid"),
          legend.position = "none")
  
  plot_list[[k]] <- p
}

#Output
pdf(file=here('output','figures','figure36.pdf'), width=15, height=8)
grid.arrange(grobs=plot_list, ncol=5, nrow=2, padding=0) # Define the layout for the plots
dev.off()

#-------------------------------------------------------------------------------
##Animate through time slices to display proportion of distribution before specified time -- FIGURE 36_animate

#Extract proportion of MCMC samples occurring before the specified time threshold (note more time divisions for animation)
time_slices2 <- seq(6000, 500, -50) 

prop_modelb_icar2  <- lapply(time_slices2, 
                             function(t) data.frame(x = 1:47, 
                                                    y = sapply(as.data.frame(post.modelb.icar), 
                                                               prop_gthan_threshold, 
                                                               threshold = t)))

# Combine all time slices into one dataset
all_data <- lapply(seq_along(time_slices2), function(k) {
  hex_area_win %>%
    filter(area_ID %in% 1:47) %>%
    mutate(median_date = med.modelb.icar,
           prop_threshold = prop_modelb_icar2[[k]]$y,
           time_slice = -time_slices2[k]) %>% #minus sign ensures animation time proceeds in the correct direction
    st_as_sf()})

# Verify that all elements are sf objects
if (!all(sapply(all_data, inherits, "sf"))) {
  stop("One or more elements of all_data are not sf objects!")
}

# Combine all sf objects safely
all_data <- do.call(rbind, all_data)

#co-ordinates of area centroids used for label positions
label_data <- st_coordinates(all_data$area_center)

# Animated ggplot
p <- ggplot(data = all_data) +
  geom_sf(data = st_buffer(st_as_sf(sampling_win, crs = 4326), 40000), aes(color = "grey50"), color=NA) +
  geom_sf(aes(fill = median_date, alpha = prop_threshold)) +
  scale_fill_viridis_c(option = "F", direction = -1) +
  scale_alpha_continuous(range = c(0, 1)) +
  xlab('Longitude') +
  ylab('Latitude') +
  geom_text(data = label_data,
            aes(x = X, y = Y, label = paste0(all_data$median_date, " BP")),
            size = 4.5, alpha = 0.4) +  #Use geom_text instead of geom_sf_label
  theme(
    panel.background = element_rect(
      fill = "lightblue",
      colour = "lightblue",
      size = 0.5,
      linetype = "solid"
    ),
    legend.position = "none",
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 15, face = "bold"),  # X-axis label
    axis.title.y = element_text(size = 15, face = "bold")   # Y-axis label
  ) +
  labs(title = "Time Slice: {closest_state} BP") #Note: use '{closest_state}' with 'transition_states()' and '{frame_time}' with 'transition_time()'

anim <- p +
  transition_states(time_slice, transition_length = 2, state_length = 1) +
  enter_fade() +
  exit_fade()

#Save the animation
anim_save("figure36_animation.gif",
          animation = animate(anim, width = 1200, height = 900, fps = 5))



