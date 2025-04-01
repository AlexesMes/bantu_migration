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
  mutate(forest_present = case_when(area_ID %in% c(16,24,19) ~ 1, TRUE ~ 0)) #assume the forest provides a 250 year delay to expansion


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
init_phi <- init_a #init_a*0.5 (if using an intercept, beta0, in model)

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

# #-------------------------------------------------------------------------------
## Model W: ICAR and wombling integrating sample interdependence, i.e. the addition of a hierarchical model ----

modelW <- nimbleCode({
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
    a[k] <- (x[k]*beta4[k]) + phi[k] #beta0 + (x[k]*beta4[k]) + phi[k]
    beta4[k] ~ dnorm(0, sd=200); #dunif(-1, 1); #determining strength of forest covariate in each area
  }
  
  # ICAR Model prior to capture spatial random effects
  phi[1:n_areas] ~ dcar_normal(adj[1:L], weights[1:L], num[1:n_areas], tau1, zero_mean =0)
  
  #For Each Boundary
  for (t in 1:n_trans){
    #nabla defines the difference in arrival time across a boundary
    nabla[t] <- abs(a[edge_id1[t]] - a[edge_id2[t]]) #edge t: select first area, m, and second area, n
  }
  
  #Priors
  #beta0 ~ dnorm(2500, sd=500); #Intercept
  tau1 ~ dgamma(1, 0.1);  #weak prior for ICAR model -- spatial autocorrelation precision parameter
  
  #tau.err <- 1/sigma^2;
  #sigma ~ dunif(0,100);
  
  # Hyperprior for duration
  gamma1 ~ dunif(1,20); #Hyperprior for rate
  gamma2 ~ T(dnorm(mean=200, sd=100), 1, 500) #Hyperprior for mode
  
})

#Define initial values ----
dW <- list(theta=sim_df$cra,
            constraint_uniform = rep(1, constants$n_areas),
            x = hex_area_win$forest_present)


initsW <- list(a=init_a,
                b=init_b,
                alpha=alpha_init,
                delta=delta_init,
                phi=init_phi,
                tau1=rgamma(1, shape = 1, rate = 0.1),
                #sigma= runif(1,0,100),
                #beta0=rnorm(1, 2500, 500),
                beta4=rnorm(1:constants$n_areas, 0, sd=200), #runif(1:constants$n_areas, min = -1, max = 1)
                gamma1=10,
                gamma2=200)


#Run MCMC ----
mcmc.samplesW <- nimbleMCMC(code = modelW,
                             constants = constants,
                             data = dW,
                             niter = 200000,
                             nchains = 4,
                             thin= 100,
                             nburnin = 100000,
                             monitors = c('a', 'b', 'theta', 'nabla', 'delta', 'alpha', 'phi', 'beta4'), #beta0
                             inits = initsW,
                             samplesAsCodaMCMC=TRUE)

#Diagnostics ----
rhatW  <- gelman.diag(mcmc.samplesW, multivariate = FALSE)
essW  <- effectiveSize(mcmc.samplesW)

#-------------------------------------------------------------------------------
# Save output ----
save(mcmc.samplesW, rhatW, essW,
     file=here('output','Womblemodel_tactsim_covariate.RData'))
