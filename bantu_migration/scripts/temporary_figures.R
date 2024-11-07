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

`%!in%` <- Negate(`%in%`)


#===============================================================================
#Tactical simulation of ICAR model 

##Load Data ----
load(here("output", "ICARmodel_tactsim2.RData"))
load(here('data', 'tactical_sim_ICAR.RData'))
load(here('data','trig.RData')) #nodes and edges between hex area centroids

#Combine constants
constants <- c(constants, constants_trig)

## Tactical Simulation Posterior Predictive Check for a and b in a given region -- FIGURE 14

#For model (i) and (ii) select parameters a and b (i.e. start and end date of occupation in the region)

sim_a <- constants$true_a
sim_b <- constants$true_b

pdf(file=here('output', 'figures','figure14.3.pdf'), width=14, height=18)
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
  plot(NULL, xlim=c(sim_a[[k]]+300, sim_b[[k]]-300), ylim=c(0,0.022), xlab='Cal BP', ylab='Posterior Probability')
  polygon(c(dens.i.a$x, rev(dens.i.a$x)), c(rep(0,length(dens.i.a$x)), rev(dens.i.a$y)), border=NA, col=rgb(0,0.4,0,0.5))
  #polygon(c(dens.i.b$x, rev(dens.i.b$x)), c(rep(0,length(dens.i.b$x)), rev(dens.i.b$y)), border=NA, col=rgb(0,0.4,0,0.5))
  polygon(c(dens.ii.a$x, rev(dens.ii.a$x)), c(rep(0,length(dens.ii.a$x)), rev(dens.ii.a$y)), border=NA, col=rgb(1,0.55,0,0.5))
  #polygon(c(dens.ii.b$x, rev(dens.ii.b$x)), c(rep(0,length(dens.ii.b$x)), rev(dens.ii.b$y)), border=NA, col=rgb(1,0.55,0,0.5))
  polygon(c(dens.iii.a$x, rev(dens.iii.a$x)), c(rep(0,length(dens.iii.a$x)), rev(dens.iii.a$y)), border=NA, col=rgb(0.82,0.086,0,0.5))
  #polygon(c(dens.iii.b$x, rev(dens.iii.b$x)), c(rep(0,length(dens.iii.b$x)), rev(dens.iii.b$y)), border=NA, col=rgb(0.82,0.086,0,0.5))
  #polygon(c(dens.iv.a$x, rev(dens.iv.a$x)), c(rep(0,length(dens.iv.a$x)), rev(dens.iv.a$y)), border=NA, col=rgb(0.004,0,0.82,0.5))
  #polygon(c(dens.iv.b$x, rev(dens.iv.b$x)), c(rep(0,length(dens.iv.b$x)), rev(dens.iv.b$y)), border=NA, col=rgb(0.004,0,0.82,0.5))
  abline(v=sim_a[[k]],lty=2) #abline(v=c(sim_a[[k]], sim_b[[k]]),lty=2)
  axis(3,at=sim_a[[k]],labels=TeX('$a$')) #axis(3,at=c(sim_a[[k]], sim_b[[k]]),labels=c(TeX('$a$'),TeX('$b$')))
  legend('topright', legend=c('Non-hierarchical Phase',
                              'Hierarchical Phase',
                              'Non-hierarchical ICAR'),
                              #'Non-hierarchical ICAR',
                              #'Hierarchical ICAR'), 
         fill=c('darkgreen','darkorange','darkred')) #c('darkgreen','darkorange','darkred','darkblue'))
  title(main = paste("Area", k))
}

dev.off()


# ##Show alpha values are recovered at sites -- FIGURE 14.2
# sim_alpha <- constants$true_alpha
# 
# pdf(file=here('output', 'figures','figure14.2.pdf'), width=14, height=18)
# # Define the layout for the plots
# par(mfrow = c(10, 10))
# 
# for (j in 1:100){ #all sites
#   post.model.ii <- do.call(rbind, mcmc.samples2)[ , c(j+41)]
#   dens.ii.a <- density(post.model.ii,bw = 5)
#   
#   # Plot
#   plot(NULL, xlim=c(sim_alpha[[j]]+300, sim_alpha[[j]]-300), ylim=c(0,0.06), xlab='Cal BP', ylab='Posterior Probability')
#   polygon(c(dens.ii.a$x, rev(dens.ii.a$x)), c(rep(0,length(dens.ii.a$x)), rev(dens.ii.a$y)), border=NA, col=rgb(1,0.55,0,0.5))
#   abline(v=sim_alpha[[j]],lty=2)
#   axis(3,at=sim_alpha[[j]],labels=TeX('$alpha$'))
#   legend('topright', legend='Hierarchichal', fill='darkorange')
#   title(main = paste("Site", j))
# }
# 
# dev.off()