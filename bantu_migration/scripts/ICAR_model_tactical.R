#Load Libraries ----
library(here)
library(coda)
library(nimbleCarbon)
library(rcarbon)
library(dplyr)

rm(list = ls())
`%!in%` <- Negate(`%in%`)

set.seed(123)

##ICAR model tactical simulation without errors (i.e. using dates that have no associated error and calibration)

#-------------------------------------------------------------------------------
## Data Setup ----
load(here('data', 'tactical_sim_ICAR.RData')) #'tactical_sim_ICAR.RData' 'tactical_sim_ICAR_spatial_auto.RData'
load(here('data','trig.RData')) #nodes and edges between hex area centroids

#Combine constants
constants <- c(constants, constants_trig)
#-------------------------------------------------------------------------------
## Initialise Parameters ----

# Initialise regional parameters ----
buffer <- 100
delta_init <- siteInfo$diff + buffer
alpha_init <- siteInfo$earliest + buffer/2

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

# #Gradient parameter initialisation
init_nabla <- 0
for (t in 1:constants$n_trans){
  m <- constants$edge_id1[t] #transition/edge t, select first area
  n <- constants$edge_id2[t] #transition/edge t, select second area

  init_nabla[t] <- (init_a$earliest[init_a$area_id==m] - init_a$earliest[init_a$area_id==n])/constants$edge_dist[t]
}

#Add buffer
init_a  <- init_a[ ,2] + buffer
init_b  <- init_b[ ,2] - buffer

#-------------------------------------------------------------------------------
#Spatial data ----

library(spdep)
nb_areas <- poly2nb(as(hex_area_win, 'Spatial'), queen=FALSE, row.names = hex_area_win$area_ID) #neighboring areas using sp library 
#nb_areas <- st_intersects(hex_area_win, hex_area_win, remove_self = TRUE) #neighboring areas using sf library

nbInfo <- nb2WB(nb_areas) #transform into iCAR inputs: adjacent matrix, weights, number of neighbors (for WinBUGS)

#-------------------------------------------------------------------------------
#Constants ----
constants$adj <- nbInfo$adj
constants$weights <- nbInfo$weights
constants$num <- nbInfo$num
constants$L <- length(nbInfo$adj)
constants <- constants[names(constants) %!in% c("dist_mat", "dist_org", "center_coords")] #remove constants which aren't used


#-------------------------------------------------------------------------------
## Model 0: Phase model assuming independence of samples (i.e. no hierarchical structure) ----

model0 <- nimbleCode({
  #For Each Date
  for (i in 1:n_dates){
    theta[i] ~ dunif(min = b[id_areas[id_sites[i]]], max = a[id_areas[id_sites[i]]]); 
  }
  
  #For Each Region
  for (k in 1:n_areas){
    a[k] ~ dunif(50,5000);
    b[k] ~ dunif(50,5000);
    constraint_uniform[k] ~ dconstraint(a[k]>b[k]) #In each area, start date of occupation, a_k, must be greater than the end date of occupation, b_k (note: BP dates in the positive direction)
  }
})

#Define initial values ---- 
d0 <- list(theta=sim_df$cra, ##Theta should be provided as data -- it should not be initalised (this changes when radiocarbon dates -- with associated error -- are used)
           constraint_uniform = rep(1, constants$n_areas)) #unif.const=1


inits0 <- list(a=init_a,
               b=init_b) 

#Run MCMC ----
mcmc.samples0 <- nimbleMCMC(code = model0,
                            constants = constants,
                            data = d0,
                            niter = 2000000, #50000, 
                            nchains = 4, 
                            thin= 100, #50,  
                            nburnin = 1000000, #10000,
                            monitors = c('a', 'b', 'theta'),
                            inits = inits0, 
                            samplesAsCodaMCMC=TRUE)

#Diagnostics ----
rhat0  <- gelman.diag(mcmc.samples0, multivariate = FALSE)
ess0  <- effectiveSize(mcmc.samples0)


#-------------------------------------------------------------------------------
## Model 1: Phase model adding hierarchical structure ----

model1 <- nimbleCode({
  #For Each Site
  for (j in 1:n_sites)
  {
    delta[j] ~ dgamma(gamma1, (gamma1-1)/gamma2)
    alpha[j] ~ dunif(max = a[id_areas[j]], min = b[id_areas[j]]);
  }
  
  #For Each Date
  for (i in 1:n_dates){
    theta[i] ~ dunif(min = (alpha[id_sites[i]] - (delta[id_sites[i]]+1)), max = alpha[id_sites[i]]);
  }
  
  #For Each Region
  for (k in 1:n_areas){
    a[k] ~ dunif(50,5000);
    b[k] ~ dunif(50,5000);
    constraint_uniform[k] ~ dconstraint(a[k]>b[k]) #In each area, start date of occupation, a_k, must be greater than the end date of occupation, b_k (note: BP dates in the positive direction)
  }
  # Hyperprior for duration
  gamma1 ~ dunif(1,20) #Hyperprior for rate
  gamma2 ~ T(dnorm(mean=200,sd=100), 1, 500) #Hyperprior for mode
})

#Define initial values ---- 
d1 <- list(theta=sim_df$cra, 
  constraint_uniform = rep(1, constants$n_areas)) #unif.const=1


inits1 <- list(a=init_a,
               b=init_b,
               alpha=alpha_init, 
               delta=delta_init,
               gamma1=10,
               gamma2=200) 

#Run MCMC ----
mcmc.samples1 <- nimbleMCMC(code = model1,
                            constants = constants,
                            data = d1,
                            niter = 2000000, 
                            nchains = 4, 
                            thin= 100, 
                            nburnin = 1000000,
                            monitors = c('a', 'b', 'theta', 'gamma1', 'gamma2', 'delta', 'alpha'),
                            inits = inits1, 
                            samplesAsCodaMCMC=TRUE)

#Diagnostics ----
rhat1  <- gelman.diag(mcmc.samples1, multivariate = FALSE)
ess1  <- effectiveSize(mcmc.samples1)

#-------------------------------------------------------------------------------
## Model 2: ICAR assuming independence of samples (i.e. no hierarchical structure) ----
model2 <- nimbleCode({
  #For Each Date
  for (i in 1:n_dates){
    theta[i] ~ dunif(min = b[id_areas[id_sites[i]]], max = a[id_areas[id_sites[i]]]); 
  }
  
  #For Each Region
  for (k in 1:n_areas){
    b[k] ~ dunif(50,5000);
    constraint_uniform[k] ~ dconstraint(a[k]>b[k]) #In each area, start date of occupation, a_k, must be greater than the end date of occupation, b_k (note: BP dates in the positive direction)
  }
  
  # ICAR Model Prior
  a[1:n_areas] ~ dcar_normal(adj[1:L], weights[1:L], num[1:n_areas], tau1, zero_mean =0)
  tau1 ~ dgamma(2, 0.5) #dunif(0.001, 10) #dgamma(50, 50)
  #tau1 <- 1/sigma1^2
  #sigma1 ~ dexp(1) #dunif(0,100)
  
})

#Define initial values ---- 
d2 <- list(theta=sim_df$cra, 
           constraint_uniform = rep(1, constants$n_areas)) 


inits2 <- list(a=init_a,
               b=init_b,
               #sigma1= rexp(1,1)) #runif(1,0,100))
               tau1=rgamma(1, shape = 2, rate = 0.5)) #runif(1,0,10))


#Run MCMC ----
mcmc.samples2 <- nimbleMCMC(code = model2,
                           constants = constants,
                           data = d2,
                           niter = 2000000, 
                           nchains = 4, 
                           thin= 100, 
                           nburnin = 1000000,
                           monitors = c('a', 'b', 'theta'),
                           inits = inits2, 
                           samplesAsCodaMCMC=TRUE)

#Diagnostics ----
rhat2  <- gelman.diag(mcmc.samples2, multivariate = FALSE)
ess2  <- effectiveSize(mcmc.samples2)

#-------------------------------------------------------------------------------
# Model 3: ICAR integrating sample interdependence, i.e. the addition of a hierarchical model ----

model3 <- nimbleCode({
  #For Each Site
  for (j in 1:n_sites)
  {
    delta[j] ~ dgamma(gamma1, (gamma1-1)/gamma2)
    alpha[j] ~ dunif(max = a[id_areas[j]], min = b[id_areas[j]]);
  }
  
  #For Each Date
  for (i in 1:n_dates){
    theta[i] ~ dunif(min = (alpha[id_sites[i]] - (delta[id_sites[i]]+1)), max = alpha[id_sites[i]]);
  }
  
  #For Each Region
  for (k in 1:n_areas){
    b[k] ~ dunif(50,5000);
    constraint_uniform[k] ~ dconstraint(a[k]>b[k]) #In each area, start date of occupation, a_k, must be greater than the end date of occupation, b_k (note: BP dates in the positive direction)
  }
  
  # ICAR Model Prior
  a[1:n_areas] ~ dcar_normal(adj[1:L], weights[1:L], num[1:n_areas], tau1, zero_mean =0)
  tau1 ~ dgamma(1, 0.1)  #dgamma(2, 0.5) 
  #tau1 <- 1/sigma1^2
  #sigma1 ~ dexp(1) #dunif(0,100)
  # Hyperprior for duration
  gamma1 ~ dunif(1,20) #Hyperprior for rate
  gamma2 ~ T(dnorm(mean=200, sd=100), 1, 500) #Hyperprior for mode
  
})

#Define initial values ---- 
d3 <- list(theta=sim_df$cra, 
           constraint_uniform = rep(1, constants$n_areas)) 


inits3 <- list(a=init_a,
               b=init_b,
               alpha=alpha_init, 
               delta=delta_init,
               tau1=rgamma(1, shape = 1, rate = 0.1), #runif(1,0,20)
               #sigma1= rexp(1,1), #runif(1,0,100),
               gamma1=10,
               gamma2=200)


#Run MCMC ----
mcmc.samples3 <- nimbleMCMC(code = model3,
                            constants = constants,
                            data = d3,
                            niter = 2000000, 
                            nchains = 4, 
                            thin= 100, 
                            nburnin = 1000000,
                            monitors = c('a', 'b', 'theta', 'delta','alpha'),
                            inits = inits3, 
                            samplesAsCodaMCMC=TRUE)

#Diagnostics ----
rhat3  <- gelman.diag(mcmc.samples3, multivariate = FALSE)
ess3  <- effectiveSize(mcmc.samples3)

#-------------------------------------------------------------------------------
# Save output ----
save(mcmc.samples0, rhat0, ess0,
     mcmc.samples1, rhat1, ess1,
     mcmc.samples2, rhat2, ess2, 
     mcmc.samples3, rhat3, ess3,
     file=here('output','ICARmodel_tactsim.RData'))