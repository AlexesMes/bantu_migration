#Load Libraries ----
library(here)
library(coda)
library(nimbleCarbon)
library(rcarbon)
library(dplyr)

rm(list = ls())
`%!in%` <- Negate(`%in%`)

set.seed(123)

#-------------------------------------------------------------------------------
## Data Setup ----
load(here('data', 'tactical_sim_phase_TOY.RData'))

#-------------------------------------------------------------------------------
## Initialise Parameters ----

# Initialise regional parameters ----
buffer <- 100

# Initialise hex areas which contain sites
init_a  <- aggregate(earliest~area_id, FUN=max, data=siteInfo) #find earliest date in each region k
init_b  <- aggregate(latest~area_id, FUN=min, data=siteInfo) #find latest date in each area k

#Initialise hex areas which do not contain sites
init_empty_area <- function(init_df) {
  for(i in 1:constants$n_area){
    
    area_ids <- init_df$area_id #List of hex areas ids with sites
    
    if (i %!in% area_ids){
      empty_hex_id <- i #Id of empty hex
      neighbour_hex <- which.min(abs(i - area_ids)) #Determine closest hex neighbor which has sites. If there are more than one neighbour hex with sites, it selects the first observation (i.e. hex with the smallest id, since ids are in ascending order). 
      neighbour_hex_id <- area_ids[neighbour_hex] #Determine area id of closest hex neighbor
      neighbour_date <- init_df[neighbour_hex , 2] #Select the date associated with the neighbour hex 
      init_df <- rbind(init_df, c(i, neighbour_date)) #Assign this date to the empty hex
    }
  }
  return(init_df)
}
init_a <- init_a %>% init_empty_area() %>%  arrange(area_id)
init_b <- init_b %>% init_empty_area() %>%  arrange(area_id)

#Add buffer
init_a  <- init_a[ ,2] + buffer
init_b  <- init_b[ ,2] - buffer

#-------------------------------------------------------------------------------
#Spatial data

library(spdep)
nb_areas <- poly2nb(as(hex_area_win, 'Spatial'), queen=FALSE, row.names = hex_area_win$area_ID) #neighboring areas using sp library 
#nb_areas <- st_intersects(hex_area_win, hex_area_win, remove_self = TRUE) #neighboring areas using sf library

nbInfo <- nb2WB(nb_areas) #transform into iCAR inputs: adjacent matrix, weights, number of neighbors (for WinBUGS)

#-------------------------------------------------------------------------------
## Model assuming independence of samples ----

model1 <- nimbleCode({
  for (i in 1:n_dates){
    theta[i] ~ dunif(min = b[id_area[i]], max = a[id_area[i]]);
  }
  
  # Set Prior for Each Region
  for (k in 1:n_areas){
    b[k] ~ dunif(50,5000);
    constraint_uniform[k] ~ dconstraint(a[k]>b[k]) #In each area, start date of occupation, a_k, must be greater than the end date of occupation, b_k (note: BP dates in the positive direction)
  }

  # ICAR Model Prior
  a[1:n_areas] ~ dcar_normal(adj[1:L], weights[1:L], num[1:n_areas], tau1, zero_mean =0)
  tau1 <- 1/sigma1^2
  sigma1 ~ dunif(0,100)
  
})

#Constants ----
constants$n_dates <- sim_constants$n_dates 
constants$n_areas <- constants$n_areas
constants$adj <- nbInfo$adj
constants$weights <- nbInfo$weights
constants$num <- nbInfo$num
constants$L <- length(nbInfo$adj)
constants$id_area <- sim_constants$id_areas
constants$id_sites <- sim_constants$id_sites

#Define initial values ---- 
d1 <- list(theta=sim_df$cra, 
           constraint_uniform = rep(1, constants$n_areas)) 


inits1 <- list(a=init_a,
               b=init_b,
               sigma1=runif(1,0,100))


#Run MCMC ----
mcmc.samples1 <- nimbleMCMC(code = model1,
                           constants = constants,
                           data = d1,
                           niter = 200, 
                           nchains = 3, 
                           thin = 5, 
                           nburnin = 100,
                           monitors = c('a','b','theta'), 
                           inits = inits1, 
                           samplesAsCodaMCMC=TRUE)

#Diagnostics ----
#rhat1  <- gelman.diag(mcmc.samples1, multivariate = FALSE)
#ess1  <- effectiveSize(mcmc.samples1)

#-------------------------------------------------------------------------------
# Save output ----
save(mcmc.samples1, file=here('output','phasemodel_tactsim_TOY.RData'))
#rhat1, ess1


#--------------------------
##Plot ----
#
##Inputs
# nsim  <- 1000
# sigma1.prior  <- runif(nsim, 0, 100)
# tau1.prior  <- 1/sigma1.prior^2
# s1.mat = matrix(NA, nrow=nsim, ncol=length(0:1000))
# for (i in 1:nsim)
# {
#   s1.mat[i,] = 
#     #s1[1:n_areas] ~ dcar_normal(adj[1:L], weights[1:L], num[1:n_areas], tau1, zero_mean =0)
# }
# 
# plot(NULL, xlab='Distance (km)', ylab='Covariance', xlim=c(0,1000), ylim=c(0,0.2))
# polygon(c(0:1000, 1000:0), c(apply(tau1.mat, 2, quantile, 0.025), rev(apply(tau1.mat, 2, quantile,0.975))), border=NA, col=rgb(0.67,0.84,0.9,0.5))
# 
# plot(tau1.prior)


## Traceplot of start and end of occupation (a, b) 
#Load Data ----
load(here("output", "phasemodel_tactsim_TOY.RData"))

#par(mfrow=c(2,2))
traceplot(mcmc.samples1[,'a[18]'], smooth=TRUE) #region area 18
traceplot(mcmc.samples1[,'b[18]'], smooth=TRUE)

#load(here("output", "phasemodel_tactsim.RData"))
#traceplot(mcmc.samples2[,'a'], smooth=TRUE)
#traceplot(mcmc.samples2[,'b'], smooth=TRUE)
#-------------------------
