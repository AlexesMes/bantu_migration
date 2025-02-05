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
load(here('data', 'tactical_sim_ICAR.RData'))
load(here('data','boundary_edges.RData')) #nodes and edges between hex area centroids

load(here('data','environ_vars.RData')) #load areal environmental variables
load(here('data','elevation.RData')) #load areal elevation data

#Combine constants
constants <- constants[names(constants) %!in% names(constants_trig)]
constants <- c(constants, constants_trig)

#-------------------------------------------------------------------------------
## Environmental data

#Exceptions: Assign hex areas which are too small to have elevation data that of their nearest neighbor.
mean_hex_elv[1,1] <- mean_hex_elv[3,1] 
mean_hex_elv[16,1] <- mean_hex_elv[20,1] 
  
environ_df <- cbind(mean_hex_clim_df, mean_hex_elv) 


#Select 3 covariates (as a first example from Cologne PCA -- later determine which covariates are important!)
environ_df <- environ_df %>% 
  select("BWA_elv_msk", #Elevation
         "bio07", #temperature annual range
         "bio14") #precipitation of the driest month
         #"rugosity") #measure of surface changes 

#-------------------------------------------------------------------------------
# ## Select relevant hexagons ----
# 
# sites <- sites %>% 
#   filter(area_id %in% c(29, 24, 19))
# 
# siteInfo <- siteInfo %>% 
#   filter(area_id %in% c(29, 24, 19))  
# 
# tiles <- c(tiles[29], tiles[24], tiles[19])
# 
# hex_area_win <- hex_area_win %>% 
#   filter(area_ID %in% c(29, 24, 19)) 
# 
# edge_info <- edge_info %>% 
#   filter(region1_id %in% c(29, 24, 19) & region2_id %in% c(29, 24, 19))

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
# Model Womble: Hierarchical ICAR with boundary identified ----

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
    
    a[k] ~ dnorm(mu[k], tau.err)
    mu[k] <- beta0 + (x1[k]*beta1 + x2[k]*beta2 + x3[k]*beta3) + phi[k]
    #a[k] <- phi[k]
    }
  
  # ICAR Model prior to capture spatial random effects
  phi[1:n_areas] ~ dcar_normal(adj[1:L], weights[1:L], num[1:n_areas], tau1, zero_mean =0)
  #a[1:n_areas] ~ dcar_normal(adj[1:L], weights[1:L], num[1:n_areas], tau1, zero_mean =0)
  
  #For each boundary
  #for (t in 1:n_trans){
  #  #nabla defines the difference in arrival time across a boundary
  #  nabla[t] <- abs(a[edge_id1[t]] - a[edge_id2[t]]) #edge t: select first area, m, and second area, n
  #}
  
  #Priors
  beta0 ~ dnorm(3300, sd=200); #Intercept
  tau1 ~ dgamma(1, 0.1);  #weak prior for ICAR model -- spatial autocorrelation precision parameter 
  
  beta1 ~ dnorm(0, 0.000001);
  beta2 ~ dnorm(0, 0.000001);
  beta3 ~ dnorm(0, 0.000001);
  
  tau.err <- 1/sigma^2;
  sigma ~ dunif(0,100);
  
  # Hyperprior for duration
  gamma1 ~ dunif(1,20); #Hyperprior for rate
  gamma2 ~ T(dnorm(mean=200, sd=100), 1, 500) #Hyperprior for mode
})

#Define initial values and data ---- 
dW <- list(theta=sim_df$cra, 
           constraint_uniform = rep(1, constants$n_areas),
           x1 = environ_df[,1],
           x2 = environ_df[,2],
           x3 = environ_df[,3])


initsW <- list(a=init_a,
               b=init_b,
               alpha=alpha_init, 
               delta=delta_init,
               nabla = init_nabla,
               tau1=rgamma(1, shape = 1, rate = 0.1),
               sigma= runif(1,0,100),
               beta0=rnorm(1, 3300, 200),
               beta1=rnorm(1, 0, 0.000001),
               beta2=rnorm(1, 0, 0.000001),
               beta3=rnorm(1, 0, 0.000001),
               gamma1=10,
               gamma2=200)


#Run MCMC ----
mcmc.samplesW <- nimbleMCMC(code = modelW,
                            constants = constants,
                            data = dW,
                            niter = 5, #2000000, 
                            nchains = 4, 
                            thin= 1, #100, 
                            nburnin = 2, #1000000,
                            monitors = c('a', 'b'), #c('a', 'b', 'theta', 'delta', 'alpha', 'nabla'),
                            inits = initsW, 
                            samplesAsCodaMCMC=TRUE)

#Diagnostics ----
rhatW  <- gelman.diag(mcmc.samplesW, multivariate = FALSE)
essW  <- effectiveSize(mcmc.samplesW)


