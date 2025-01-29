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
#load(here('data', 'eastc14.RData')) #East and Southern Africa
#load(here('data','trig.RData')) #nodes and edges between hex area centroids
load(here('data', 'c14.RData')) #sub-Saharan Africa
load(here('data','trig_cont.RData'))

#Combine constants
constants <- c(constants, constants_trig)
#-------------------------------------------------------------------------------
# General Setup ----
# Data --
dat <- list(cra = dateInfo$cra, #theta=sim_df$cra
            cra_error = dateInfo$cra_error,
            constraint_uniform = rep(1, constants$n_areas),
            cra_constraint = rep(1, constants$n_dates)) # Set-up constraint for ignoring inference outside calibration range

#Calibration curve
constants$cc <- as.numeric(as.factor(dateInfo$calCurve)) #intcal20==1 and shcal20==2

# Dummy extension of the calibration curve
constants$calBP <- c(1000000, constants$calBP, -1000000)
constants$C14BP <- rbind(c(1000000,1000000), constants$C14BP, c(-1000000,-1000000))
constants$C14err <- rbind(c(1000,1000), constants$C14err, c(1000,1000))

## Initialise Parameters ----
theta_init <- dateInfo$median_dates

# Initialise regional parameters
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

#Gradient parameter initialisation
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
#Spatial data
nb_areas <- poly2nb(as(hex_area_win, 'Spatial'), queen=FALSE, row.names = hex_area_win$area_ID) #neighboring areas using sp library 
#nb_areas <- st_intersects(hex_area_win, hex_area_win, remove_self = TRUE) #neighboring areas using sf library

nbInfo <- nb2WB(nb_areas) #transform into iCAR inputs: adjacent matrix, weights, number of neighbors (for WinBUGS)
#help('CAR-Normal') for details of these parameters
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
                                                #"eastEIAcountries",
                                                "origin_point")] #remove constants which aren't used

#-------------------------------------------------------------------------------
## ICAR Model assuming independence of samples ----
icar_model  <- function(seed, d, theta_init, init_nabla, init_a, init_b, constants, nburnin, thin, niter)
{
  #Load Library
  library(nimbleCarbon)
  #Define Core Model
  model <- nimbleCode({
    for (i in 1:n_dates){
      theta[i] ~ dunif(min = b[id_areas[id_sites[i]]], max = a[id_areas[id_sites[i]]]);
      #Calibration
      mu[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[ , cc[i]]); #c14age #Index cc selects the correct calibration curve
      cra_constraint[i] ~ dconstraint(mu[i] < 50193 & mu[i] > 95) #C14 age must be within the calibration range
      sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[ , cc[i]]);
      sd[i] <- (cra_error[i]^2 + sigmaCurve[i]^2)^(1/2);
      cra[i] ~ dnorm(mean=mu[i], sd=sd[i]);
    }

    #For Each Region
    for (k in 1:n_areas){
      b[k] ~ dunif(50,5000);
      constraint_uniform[k] ~ dconstraint(a[k]>b[k]) #In each area, start date of occupation, a_k, must be greater than the end date of occupation, b_k (note: BP dates in the positive direction)
    }

    #For each edge transition
    for (t in 1:n_trans){
      #nabla defines the gradient along the edge
      nabla[t] <- (a[edge_id1[t]] - a[edge_id2[t]])/edge_dist[t] #edge t: select first area, m, and second area, n
    }

    # ICAR Model Prior
    a[1:n_areas] ~ dcar_normal(adj[1:L], weights[1:L], num[1:n_areas], tau, zero_mean =0)
    tau ~ dgamma(1, 0.1)
  })

  #Define initial values ----
  inits <- list(theta=theta_init,
                 a=init_a,
                 b=init_b,
                 nabla=init_nabla,
                 tau=rgamma(1, shape = 1, rate = 0.1))

  # Compile and Run model	----
  model <- nimbleModel(model, constants=constants, data=d, inits=inits)
  cModel <- compileNimble(model)
  conf <- configureMCMC(model, control=list(adaptInterval=20000, adaptFactorExponent=0.1))
  conf$addMonitors(c('a','b','nabla','theta'))
  MCMC <- buildMCMC(conf)
  cMCMC <- compileNimble(MCMC)
  results <- runMCMC(cMCMC, niter = niter, thin = thin, nburnin = nburnin, samplesAsCodaMCMC = T, setSeed = seed)
}

#-------------------------------------------------------------------------------
## Hierarchical ICAR Model assuming interdependence of samples ----
icar_model_b  <- function(seed, d, theta_init, init_nabla, alpha_init, delta_init, init_a, init_b, constants, nburnin, thin, niter)
{
  #Load Library
  library(nimbleCarbon)
  #Define Core Model
  model <- nimbleCode({
    #For Each Date
    for (j in 1:n_sites)
    {
      delta[j] ~ dgamma(gamma1, (gamma1-1)/gamma2)
      alpha[j] ~ dunif(max = a[id_areas[j]], min = b[id_areas[j]]);
    }
    
    #For Each Site
    for (i in 1:n_dates){
      theta[i] ~ dunif(min = (alpha[id_sites[i]] - (delta[id_sites[i]]+1)), max = alpha[id_sites[i]]);
      #Calibration
      mu[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[ , cc[i]]); #c14age #Index cc selects the correct calibration curve
      cra_constraint[i] ~ dconstraint(mu[i] < 50193 & mu[i] > 95) #C14 age must be within the calibration range
      sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[ , cc[i]]);
      sd[i] <- (cra_error[i]^2 + sigmaCurve[i]^2)^(1/2);
      cra[i] ~ dnorm(mean=mu[i], sd=sd[i]);
    }
    
    #For Each Region
    for (k in 1:n_areas){
      b[k] ~ dunif(50,5000);
      constraint_uniform[k] ~ dconstraint(a[k]>b[k]) #In each area, start date of occupation, a_k, must be greater than the end date of occupation, b_k (note: BP dates in the positive direction)
    }
    
    #For each edge transition
    for (t in 1:n_trans){
      #nabla defines the gradient along the edge
      nabla[t] <- (a[edge_id1[t]] - a[edge_id2[t]])/edge_dist[t] #edge t: select first area, m, and second area, n
    }
    
    # ICAR Model Prior
    a[1:n_areas] ~ dcar_normal(adj[1:L], weights[1:L], num[1:n_areas], tau, zero_mean =0)
    tau ~ dgamma(1, 0.1)
    
    # Hyperprior for duration
    gamma1 ~ dunif(1,20) #Hyperprior for rate
    gamma2 ~ T(dnorm(mean=200,sd=100), 1, 500) #Hyperprior for mode
  })
  
  #Define initial values ---- 
  inits <- list(theta=theta_init,
                a=init_a,
                b=init_b,
                alpha=alpha_init, 
                delta=delta_init, 
                nabla=init_nabla,
                tau=rgamma(1, shape = 1, rate = 0.1),
                gamma1=10,
                gamma2=200)
  
  # Compile and Run model	----
  model <- nimbleModel(model, constants=constants, data=d, inits=inits)
  cModel <- compileNimble(model)
  conf <- configureMCMC(model, control=list(adaptInterval=20000, adaptFactorExponent=0.1))
  conf$addMonitors(c('a','b','nabla','theta','delta','alpha'))
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
niter  <- 1000000
nburnin  <- 500000
thin  <- 100

#Model A -- ICAR Model
out_icar_model_a  <-  parLapply(cl = cl, 
                                X = seeds, 
                                fun = icar_model,
                                d = dat, 
                                constants = constants, 
                                theta_init = theta_init, 
                                init_a = init_a, 
                                init_b = init_b, 
                                init_nabla = init_nabla,
                                niter = niter, 
                                nburnin = nburnin,
                                thin = thin)
out_icar_model_a <- mcmc.list(out_icar_model_a)

#Model B -- Hierarchical ICAR Model 
out_icar_model_b  <-  parLapply(cl = cl, 
                                X = seeds, 
                                fun = icar_model_b, 
                                d = dat, 
                                constants = constants, 
                                theta_init = theta_init, 
                                init_a = init_a, 
                                init_b = init_b, 
                                init_nabla = init_nabla,
                                alpha_init = alpha_init, 
                                delta_init = delta_init,
                                niter = niter, 
                                nburnin = nburnin,
                                thin = thin)
out_icar_model_b <- mcmc.list(out_icar_model_b)

## Diagnostics ----
#Model A
rhat_icar_model_a  <- gelman.diag(out_icar_model_a, multivariate = FALSE)
ess_icar_model_a  <- effectiveSize(out_icar_model_a)
agg_icar_model_a <- agreementIndex(dat$cra,
                                   dat$cra_error,
                                   calCurve = dateInfo$calCurve,
                                   theta = out_icar_model_a[[1]][ , grep("theta", colnames(out_icar_model_a[[1]]))],
                                   verbose = F)
#Model B
rhat_icar_model_b  <- gelman.diag(out_icar_model_b, multivariate = FALSE)
ess_icar_model_b  <- effectiveSize(out_icar_model_b)
agg_icar_model_b <- agreementIndex(dat$cra,
                                   dat$cra_error,
                                   calCurve = dateInfo$calCurve,
                                   theta = out_icar_model_b[[1]][ , grep("theta", colnames(out_icar_model_b[[1]]))],
                                   verbose = F)

#-------------------------------------------------------------------------------
# Save output ----
save(out_icar_model_a, 
     rhat_icar_model_a, 
     ess_icar_model_a, 
     agg_icar_model_a, 
     file=here('output','ICAR_model_a_cont.RData')) #'ICAR_model_a.RData'
save(out_icar_model_b, 
     rhat_icar_model_b, 
     ess_icar_model_b, 
     agg_icar_model_b, 
     file=here('output','ICAR_model_b_cont.RData')) #'ICAR_model_b.RData'
