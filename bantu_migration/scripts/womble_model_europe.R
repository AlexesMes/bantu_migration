#Load Libraries ----
library(here)
library(coda)
library(nimbleCarbon)
library(rcarbon)
library(dplyr)
library(parallel)
library(spdep)

rm(list = ls())
`%!in%` <- Negate(`%in%`)

#-------------------------------------------------------------------------------
## Data Setup ----
load(here('data', 'europe.RData'))

load(here('data','trig_europe.RData')) #nodes and edges between hex area centroids

#Combine constants
constants <- c(constants, constants_trig)
#-------------------------------------------------------------------------------
# General Setup ----
# Data --
dat <- list(cra = europe_sites_df$c14date, 
            cra_error = europe_sites_df$c14std,
            constraint_uniform = rep(1, constants$n_areas),
            cra_constraint = rep(1, constants$n_sites)) # Set-up constraint for ignoring inference outside calibration range

#Calibration curve
constants$cc <- as.numeric(as.factor(europe_sites_df$calCurve)) #intcal20==1

# Dummy extension of the calibration curve
constants$calBP <- c(1000000, constants$calBP, -1000000)
constants$C14BP <- c(1000000, constants$C14BP, -1000000)
constants$C14err <- c(1000, constants$C14err, 1000)

## Initialise Parameters ----
theta_init <- europe_sites_df$median_dates

#===============================================================================
# Initialise regional parameters
buffer <- 300
delta_init <-  rep(buffer, constants$n_sites)
alpha_init <- europe_sites_df$median_dates + buffer/2

# Initialise hex areas which contain sites
init_a  <- aggregate(median_dates~area_id, FUN=max, data=europe_sites_df) #find earliest date in each region k
init_b  <- aggregate(median_dates~area_id, FUN=min, data=europe_sites_df) #find latest date in each area k

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

# Initialise spatial residues
init_phi <- init_a*1.5

#-------------------------------------------------------------------------------
#Spatial data
nb_areas <- poly2nb(as(hex_area_win, 'Spatial'), queen=FALSE, row.names = hex_area_win$area_ID) #neighboring areas using sp library 
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
                                                "countries",
                                                "origin_point")] #remove constants which aren't used








# model <- nimbleCode({
#   #For Each Site
#   for (j in 1:n_sites)
#   {
#     delta[j] ~ dgamma(gamma1, (gamma1-1)/gamma2)
#     alpha[j] ~ dunif(max = a[id_areas[j]], min = b[id_areas[j]]);
#     theta[j] ~ dunif(min = (alpha[j] - (delta[j]+1)), max = alpha[j]);
#   }
#   #For Each Region
#   for (k in 1:n_areas){
#     b[k] ~ dunif(0, 3000);
#     constraint_uniform[k] ~ dconstraint(b[k]<a[k]); #In each area, start date of occupation, a_k, must be greater than the end date of occupation, b_k (note: BP dates in the positive direction)
#     a[k] <- phi[k];
#   }
#   
#   # ICAR Model prior to capture spatial random effects
#   phi[1:n_areas] ~ dcar_normal(adj[1:L], weights[1:L], num[1:n_areas], tau1, zero_mean=0)
#   
#   #Priors
#   tau1 ~ dgamma(0.8, 0.1);  #weak prior for ICAR model -- spatial autocorrelation precision parameter
#   
#   # Hyperprior for duration
#   gamma1 ~ dunif(1,20); #Hyperprior for rate
#   gamma2 ~ T(dnorm(mean=200, sd=100), 1, 500) #Hyperprior for mode
# })
# 
# 
# #Define initial values ----
# dW <- list(theta=europe_sites_df$c14date,
#            constraint_uniform = rep(1, constants$n_areas))
# 
# 
# initsW <- list(a=init_a,
#   b=init_b,
#   alpha=alpha_init,
#   delta=delta_init,
#   phi=init_phi,
#   tau1=rgamma(1, shape = 0.8, rate = 0.1),
#   gamma1=10,
#   gamma2=200)
# 
# 
# #Run MCMC ----
# mcmc.samplesW <- nimbleMCMC(code = model,
#                             constants = constants,
#                             data = dW,
#                             niter = 200,
#                             nchains = 4,
#                             thin= 10,
#                             nburnin = 100,
#                             monitors = c('a', 'b', 'theta', 'delta', 'alpha', 'phi'), 
#                             inits = initsW,
#                             samplesAsCodaMCMC=TRUE)

















#===============================================================================
#-------------------------------------------------------------------------------
# Model Womble: Hierarchical ICAR with boundary identified ----

modelW <- function(seed, d, theta_init, alpha_init, delta_init, init_a, init_b, init_phi, constants, nburnin, thin, niter)
{
  #Load Library
  library(nimbleCarbon)
  #Define Core Model
  model <- nimbleCode({
    #For Each Site
    for (j in 1:n_sites)
    {
      delta[j] ~ dgamma(gamma1, (gamma1-1)/gamma2)
      alpha[j] ~ dunif(max = a[id_areas[j]], min = b[id_areas[j]]);
    }
    
    #For Each Date
    for (i in 1:n_sites){ #n_dates=n_sites
      theta[i] ~ dunif(min = (alpha[i] - (delta[i]+1)), max = alpha[i]);
      #Calibration
      mu[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]); #c14age 
      cra_constraint[i] ~ dconstraint(mu[i] < 50193 & mu[i] > 95) #C14 age must be within the calibration range
      sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]);
      sd[i] <- (cra_error[i]^2 + sigmaCurve[i]^2)^(1/2);
      cra[i] ~ dnorm(mean=mu[i], sd=sd[i]);
    }
    
    #For Each Region
    for (k in 1:n_areas){
      b[k] ~ dunif(50,5000);
      constraint_uniform[k] ~ dconstraint(a[k]>b[k]) #In each area, start date of occupation, a_k, must be greater than the end date of occupation, b_k (note: BP dates in the positive direction)
      
      a[k] <- phi[k]; #x3[k]*beta3 + phi[k] + beta0 #x1[k]*beta1 + x2[k]*beta2; 
    }
    
    # ICAR Model prior to capture spatial random effects
    phi[1:n_areas] ~ dcar_normal(adj[1:L], weights[1:L], num[1:n_areas], tau1, zero_mean =0)
    
    #For Each Edge Transition
    for (t in 1:n_trans){
      #nabla defines the difference in arrival time across a boundary
      nabla[t] <- abs(a[edge_id1[t]] - a[edge_id2[t]]) #edge t: select first area, m, and second area, n
      #nabla_phi defines the difference in spatial residues across a boundary
      nabla_phi[t] <- abs(phi[edge_id1[t]] - phi[edge_id2[t]])
    }
    
    #Priors
    tau1 ~ dgamma(1, 0.1);  #weak prior for ICAR model -- spatial autocorrelation precision parameter
    
    # Hyperprior for duration
    gamma1 ~ dunif(1,20) #Hyperprior for rate
    gamma2 ~ T(dnorm(mean=200,sd=100), 1, 500) #Hyperprior for mode
  })
  
  #Define initial values ---- 
  inits <- list(theta=theta_init,
                a=init_a,
                b=init_b,
                phi=init_phi,
                alpha=alpha_init, 
                delta=delta_init, 
                tau=rgamma(1, shape = 1, rate = 0.1),
                gamma1=10,
                gamma2=200)
  
  # Compile and Run model	----
  model <- nimbleModel(model, constants=constants, data=d, inits=inits)
  cModel <- compileNimble(model)
  conf <- configureMCMC(model, control=list(adaptInterval=20000, adaptFactorExponent=0.1))
  conf$addMonitors(c('a','b','nabla', 'nabla_phi','theta','delta','alpha'))
  MCMC <- buildMCMC(conf)
  cMCMC <- compileNimble(MCMC)
  results <- runMCMC(cMCMC, niter = niter, thin = thin, nburnin = nburnin, samplesAsCodaMCMC = T, setSeed = seed) 
}

#-------------------------------------------------------------------------------
## Run MCMCs ----

# MCMC Setup
ncores  <-  4
cl <- makeCluster(ncores)
seeds <- c(12, 34, 56, 78)
niter  <- 2000000
nburnin  <- 1000000
thin  <- 100

#Hierarchical Womble Model 
out_womble_model <-  parLapply(cl = cl, 
                               X = seeds, 
                               fun = modelW, 
                               d = dat, 
                               constants = constants, 
                               theta_init = theta_init, 
                               init_a = init_a, 
                               init_b = init_b, 
                               alpha_init = alpha_init, 
                               delta_init = delta_init,
                               init_phi = init_phi,
                               niter = niter, 
                               nburnin = nburnin,
                               thin = thin)
out_womble_model <- mcmc.list(out_womble_model)

## Diagnostics ----
rhat_womble_model  <- gelman.diag(out_womble_model, multivariate = FALSE)
ess_womble_model  <- effectiveSize(out_womble_model)
agg_womble_model <- agreementIndex(dat$cra,
                                   dat$cra_error,
                                   calCurve = dateInfo$calCurve,
                                   theta = out_womble_model[[1]][ , grep("theta", colnames(out_womble_model[[1]]))],
                                   verbose = F)

#-------------------------------------------------------------------------------
# Save output ----
save(out_womble_model, 
     rhat_womble_model, 
     ess_womble_model, 
     agg_womble_model, 
     file=here('output','Womble_europe.RData'))
