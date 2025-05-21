#Load Libraries ----
library(here)
library(coda)
library(nimbleCarbon)
library(rcarbon)
library(dplyr)
library(msm)

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
# ##Simulated data summary
# #Number of dates in area
# dates_in_areas_summarise <- as.data.frame(table(sites$area_id))
# #Number of sites in area
# sites_in_areas_summarise <- sites %>% group_by(area_id) %>% summarize(n_sites =n_distinct(site_id))

#-------------------------------------------------------------------------------
#Covariates -- center, standardise and save in a matrix

covariates <- c("distance")
cov_df <- as.data.frame(hex_area_win$dist_from_origin)
colnames(cov_df) <- covariates
cov_scaled <- sapply(cov_df, function (x) x / max(x)) #sapply(cov_df, function(x) (x - mean(x)) / sd(x))


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
init_nabla_phi <- init_nabla*0.5

#Add buffer
init_a  <- init_a[ ,2] + buffer
init_b  <- init_b[ ,2] - buffer

# Initialise spatial residues
init_phi <- init_a*1.5 #(use init_a if no intercept (beta0) is included and zero_mean=0)

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
constants <- constants[names(constants) %!in% c("dist_mat", 
                                                "dist_org", 
                                                "center_coords",
                                                "beta0",
                                                "beta1",
                                                "beta2",
                                                "x1",
                                                "x2",
                                                "mu",
                                                "mu2",
                                                "tau.err",
                                                "tau")] #remove constants which aren't used

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
    b[k] ~ dunif(0, 3000);
    constraint_uniform[k] ~ dconstraint(b[k]<a[k]); #In each area, start date of occupation, a_k, must be greater than the end date of occupation, b_k (note: BP dates in the positive direction)
    
    a[k] <- phi[k] + x1[k]*beta1; #beta0[k] + phi[k] + X[k]*beta1; # + x1[k]*beta1[k]; #+ x2[k]*beta2; 
    #beta1[k] ~ dnorm(0, sd=600);
    #beta0[k] ~ dunif(1000, 3000);
  }
  
  # ICAR Model prior to capture spatial random effects
  phi[1:n_areas] ~ dcar_normal(adj[1:L], weights[1:L], num[1:n_areas], tau1, zero_mean=0)
  #beta1[1:n_areas] ~ dcar_normal(adj[1:L], weights[1:L], num[1:n_areas], tau2, zero_mean=0)
  
  #For Each Boundary
  for (t in 1:n_trans){
    #nabla defines the difference in arrival time across a boundary
    nabla[t] <- abs(a[edge_id1[t]] - a[edge_id2[t]]) #edge t: select first area, m, and second area, n
    #nabla_phi defines the difference in spatial residues across a boundary
    nabla_phi[t] <- abs(phi[edge_id1[t]] - phi[edge_id2[t]])
  }
  
  #Priors
  #beta0 ~ dunif(1000, 3000); #dnorm(2000, sd=300); #Intercept
  beta1 ~ dnorm(0, sd=300); #dunif(-1, 1); #determining strength of forest covariate in each area
  #beta2 ~ dnorm(0, sd=300);
  tau1 ~ dgamma(0.8, 0.1);  #weak prior for ICAR model -- spatial autocorrelation precision parameter
  #tau2 ~ dgamma(0.8, 0.1);  #weak prior for ICAR model -- spatial autocorrelation precision parameter
  
  #tau.err <- 1/sigma^2;
  #sigma ~ dunif(0,100);
  
  # Hyperprior for duration
  gamma1 ~ dunif(1,20); #Hyperprior for rate
  gamma2 ~ T(dnorm(mean=200, sd=100), 1, 500) #Hyperprior for mode
  
})


# # Generate a sequence of values for beta0
# beta0 <- seq(0, 4000, length.out = 4000)
# # Compute the density
# density0 <- dnorm(beta0, mean = 2000, sd = 300)
# # Plot the density
# plot(beta0, density0, type = "l",
#      main = "Density of beta0 ~ N(2000, 400²)", xlab = "beta0", ylab = "Density", col = "blue", lwd = 2)
# grid()

# # Generate a sequence of values for beta0
# beta1 <- seq(-2000, 2000, length.out = 4000)
# # Compute the density
# density1 <- dnorm(beta1, mean = 0, sd =  600) #dtruncnorm(beta1, a=-2000, b=0, mean = 0, sd = 300)
# # Plot the density
# plot(beta1, density1, type = "l",
#      main = "Density of beta1 ~ N(0, 200²)", xlab = "beta1", ylab = "Density", col = "blue", lwd = 2)
# grid()




#Define initial values ----
dW <- list(theta=sim_df$cra,
           constraint_uniform = rep(1, constants$n_areas),
           #X = as.data.frame(cov_scaled)$distance)
           x1 = hex_area_win$forest_present)
           #x2 = hex_area_win$water_present)


initsW <- list(#a=init_a,
                b=init_b,
                alpha=alpha_init,
                delta=delta_init,
                phi=init_phi,
                tau1=rgamma(1, shape = 0.8, rate = 0.1),
                #tau2=rgamma(1, shape = 0.8, rate = 0.1),
                #beta0=runif(1:constants$n_areas, 1000, 3000), #rnorm(1, 2000, 300),
                beta1=rnorm(1, 0, sd=300),  #rnorm(1:constants$n_areas, 0, sd=200), #runif(1:constants$n_areas, min = -1, max = 1)
                #beta2=rnorm(1, 0, sd=300),
                gamma1=10,
                gamma2=200)


#Run MCMC ----
mcmc.samplesW <- nimbleMCMC(code = modelW,
                             constants = constants,
                             data = dW,
                             niter = 2000000,
                             nchains = 4,
                             thin= 100,
                             nburnin = 1000000,
                             monitors = c('a', 'b', 'theta', 'nabla', 'nabla_phi', 'delta', 'alpha', 'phi', 'beta1'), 
                             inits = initsW,
                             samplesAsCodaMCMC=TRUE)

#Diagnostics ----
rhatW  <- gelman.diag(mcmc.samplesW, multivariate = FALSE)
essW  <- effectiveSize(mcmc.samplesW)

#-------------------------------------------------------------------------------
# Save output ----
save(mcmc.samplesW, rhatW, essW,
     file=here('output','Womblemodel_tactsim_south_forest_fixed.RData'))
