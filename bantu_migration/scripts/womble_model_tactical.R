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
load(here('data', 'tactical_sim_womble.RData'))
load(here('data','trig_cont.RData')) #nodes and edges between hex area centroids


#Combine constants
constants <- c(constants, constants_trig)

#-------------------------------------------------------------------------------
## Environmental Data ----
hex_area_win <- hex_area_win %>% 
  mutate(forest_present = case_when(area_ID %in% c(16,24,19) ~ +400, TRUE ~ 0)) #assume the forest provides a 250 year delay to expansion


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

##Delta parameter initialisation
init_nabla <- 0
for (t in 1:constants$n_trans){
  m <- constants$edge_id1[t] #transition/edge t, select first area
  n <- constants$edge_id2[t] #transition/edge t, select second area
  
  init_nabla[t] <- (init_a$earliest[init_a$area_id==m] - init_a$earliest[init_a$area_id==n])
}

#Add buffer
init_a  <- init_a[ ,2] + buffer
init_b  <- init_b[ ,2] - buffer

# Initialise spatial residues
init_phi <- init_a*0.5 #rep(0, constants$n_areas)

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
constants$beta1 <- 1 #binary: turn on/off the effect of the forest covariate
constants$beta2 <- -1 #binary: turn on/off the effect of the forest covariate
constants$beta3 <- 0 #binary: turn on/off the effect of the forest covariate
constants <- constants[names(constants) %!in% c("dist_mat", "dist_org", "center_coords")] #remove constants which aren't used

# #-------------------------------------------------------------------------------
# # Model W: ICAR integrating sample interdependence, i.e. the addition of a hierarchical model ----

##Model 1: Beta1==1

modelW1 <- nimbleCode({
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

    #a[k] ~ dnorm(phi[k], tau.err)
    a[k] <- beta0 + (x[k]*beta1) + phi[k] #mu[k] <- beta0 + (x1[k]*beta1[k]) + phi[k]
  }

  # ICAR Model prior to capture spatial random effects
  phi[1:n_areas] ~ dcar_normal(adj[1:L], weights[1:L], num[1:n_areas], tau1, zero_mean =0)

  #For Each Boundary
  for (t in 1:n_trans){
    #nabla defines the difference in arrival time across a boundary
    nabla[t] <- abs(a[edge_id1[t]] - a[edge_id2[t]]) #edge t: select first area, m, and second area, n
  }
  
  #Priors
  beta0 ~ dnorm(2500, sd=500); #Intercept
  tau1 ~ dgamma(1, 0.1);  #weak prior for ICAR model -- spatial autocorrelation precision parameter

  tau.err <- 1/sigma^2;
  sigma ~ dunif(0,100);

  # Hyperprior for duration
  gamma1 ~ dunif(1,20); #Hyperprior for rate
  gamma2 ~ T(dnorm(mean=200, sd=100), 1, 500) #Hyperprior for mode

})

#Define initial values ----
dW1 <- list(theta=sim_df$cra,
            constraint_uniform = rep(1, constants$n_areas),
            x = hex_area_win$forest_present)


initsW1 <- list(a=init_a,
                b=init_b,
                alpha=alpha_init,
                delta=delta_init,
                phi=init_phi,
                tau1=rgamma(1, shape = 1, rate = 0.1),
                sigma= runif(1,0,100),
                beta0=rnorm(1, 2500, 500),
                gamma1=10,
                gamma2=200)


#Run MCMC ----
mcmc.samplesW1 <- nimbleMCMC(code = modelW1,
                             constants = constants,
                             data = dW1,
                             niter = 2000000,
                             nchains = 4,
                             thin= 100,
                             nburnin = 1000000,
                             monitors = c('a', 'b', 'theta', 'delta', 'nabla', 'alpha','phi'),
                             inits = initsW1,
                             samplesAsCodaMCMC=TRUE)

#Diagnostics ----
rhatW1  <- gelman.diag(mcmc.samplesW1, multivariate = FALSE)
essW1  <- effectiveSize(mcmc.samplesW1)

#-------------------------------------------------------------------------------
##Model 2: Beta2==-1 

modelW2 <- nimbleCode({
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
    
    a[k] <- beta0 + (x[k]*beta2) + phi[k]
  }
  
  # ICAR Model prior to capture spatial random effects
  phi[1:n_areas] ~ dcar_normal(adj[1:L], weights[1:L], num[1:n_areas], tau1, zero_mean =0)
  
  #For Each Boundary
  for (t in 1:n_trans){
    #nabla defines the difference in arrival time across a boundary
    nabla[t] <- abs(a[edge_id1[t]] - a[edge_id2[t]]) #edge t: select first area, m, and second area, n
  }
  
  #Priors
  beta0 ~ dnorm(2500, sd=500); #Intercept
  tau1 ~ dgamma(1, 0.1);  #weak prior for ICAR model -- spatial autocorrelation precision parameter
  
  tau.err <- 1/sigma^2;
  sigma ~ dunif(0,100);
  
  # Hyperprior for duration
  gamma1 ~ dunif(1,20); #Hyperprior for rate
  gamma2 ~ T(dnorm(mean=200, sd=100), 1, 500) #Hyperprior for mode
  
})

#Define initial values ----
dW2 <- list(theta=sim_df$cra,
            constraint_uniform = rep(1, constants$n_areas),
            x = hex_area_win$forest_present)


initsW2 <- list(a=init_a,
                b=init_b,
                alpha=alpha_init,
                delta=delta_init,
                nabla = init_nabla,
                phi=init_phi,
                tau1=rgamma(1, shape = 1, rate = 0.1),
                sigma= runif(1,0,100),
                beta0=rnorm(1, 2500, 500),
                gamma1=10,
                gamma2=200)


#Run MCMC ----
mcmc.samplesW2 <- nimbleMCMC(code = modelW2,
                             constants = constants,
                             data = dW2,
                             niter = 2000000,
                             nchains = 4,
                             thin= 100,
                             nburnin = 1000000,
                             monitors = c('a', 'b', 'theta', 'delta', 'nabla', 'alpha','phi'),
                             inits = initsW2,
                             samplesAsCodaMCMC=TRUE)

#Diagnostics ----
rhatW2  <- gelman.diag(mcmc.samplesW2, multivariate = FALSE)
essW2  <- effectiveSize(mcmc.samplesW2)
 
#-------------------------------------------------------------------------------
##Model 3: Beta3==0
 
modelW3 <- nimbleCode({
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
    
    a[k] <- beta0 + (x[k]*beta3) + phi[k]
  }
  
  # ICAR Model prior to capture spatial random effects
  phi[1:n_areas] ~ dcar_normal(adj[1:L], weights[1:L], num[1:n_areas], tau1, zero_mean =0)
  
  #For Each Boundary
  for (t in 1:n_trans){
    #nabla defines the difference in arrival time across a boundary
    nabla[t] <- abs(a[edge_id1[t]] - a[edge_id2[t]]) #edge t: select first area, m, and second area, n
  }
  
  #Priors
  beta0 ~ dnorm(2500, sd=500); #Intercept
  tau1 ~ dgamma(1, 0.1);  #weak prior for ICAR model -- spatial autocorrelation precision parameter
  
  tau.err <- 1/sigma^2;
  sigma ~ dunif(0,100);
  
  # Hyperprior for duration
  gamma1 ~ dunif(1,20); #Hyperprior for rate
  gamma2 ~ T(dnorm(mean=200, sd=100), 1, 500) #Hyperprior for mode
  
})

#Define initial values ----
dW3 <- list(theta=sim_df$cra,
            constraint_uniform = rep(1, constants$n_areas),
            x = hex_area_win$forest_present)


initsW3 <- list(a=init_a,
                b=init_b,
                alpha=alpha_init,
                delta=delta_init,
                phi=init_phi,
                tau1=rgamma(1, shape = 1, rate = 0.1),
                sigma= runif(1,0,100),
                beta0=rnorm(1, 2500, 500),
                gamma1=10,
                gamma2=200)


#Run MCMC ----
mcmc.samplesW3 <- nimbleMCMC(code = modelW3,
                             constants = constants,
                             data = dW3,
                             niter = 2000000,
                             nchains = 4,
                             thin= 100,
                             nburnin = 1000000,
                             monitors = c('a', 'b', 'theta', 'nabla', 'delta','alpha','phi'),
                             inits = initsW3,
                             samplesAsCodaMCMC=TRUE)

#Diagnostics ----
rhatW3  <- gelman.diag(mcmc.samplesW3, multivariate = FALSE)
essW3  <- effectiveSize(mcmc.samplesW3)

#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Save output ----
save(mcmc.samplesW1, rhatW1, essW1,
     mcmc.samplesW2, rhatW2, essW2,
     mcmc.samplesW3, rhatW3, essW3,
     file=here('output','Womblemodel_tactsim.RData'))

