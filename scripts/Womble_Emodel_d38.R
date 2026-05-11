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

#===============================================================================
###ICAR MODEL WITH ENVIRONMENTAL COVARIATES -- npp and rugosity
#===============================================================================
#-------------------------------------------------------------------------------
## Data Setup ----
#Chose model and spatial scale
df_dat <- c('eastc14.RData','trig_d38.RData','Krapp_enviro_variables.RData') #ICAR_d38: ICAR mode with scale at d=3.8
#df_dat <- c('eastc14_wanAB.RData','trig_d38.RData','Krapp_enviro_variables.RData') #ICAR_d38_wanAB: dataset with Wanyika grade C and D dates filtered out (~15% of the original data)
#df_dat <- c('eastc14_d29.RData','trig_d29.RData','Krapp_enviro_variables_d29.RData') #ICAR_d29: model running with a sample window with a finer spatial scale (i.e. 83 hexagonal units, d=2.9)

#Load data
load(here('data', df_dat[1])) #East and Southern Africa data
load(here('data', df_dat[2])) #nodes and edges between hex area centroids
load(here('data', df_dat[3])) #Environmental data

#Combine constants
constants <- c(constants, constants_trig)
#-------------------------------------------------------------------------------
# General Setup ----
# Data --
dat <- list(cra = dateInfo$cra, #theta=sim_df$cra
            cra_error = dateInfo$cra_error,
            constraint_uniform = rep(1, constants$n_areas),
            cra_constraint = rep(1, constants$n_dates), # Set-up constraint for ignoring inference outside calibration range
            x1 = scaled_hex_clim_df[,"npp"],
            x2 = scaled_hex_clim_df[,"rugosity"]) 

# #Center covariates
# dat$x1_centered <- dat$x1 - mean(dat$x1)

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
  
  init_nabla[t] <- (init_a$earliest[init_a$area_id==m] - init_a$earliest[init_a$area_id==n])
}
init_nabla_phi <- init_nabla*0.5

#Add buffer
init_a  <- init_a[ ,2] + buffer
init_b  <- init_b[ ,2] - buffer

# Initialise spatial residues
init_phi <- init_a #(use init_a if no intercept (beta0) is included and zero_mean=0)
init_phi_raw <- rep(0, constants$n_areas)

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
                                                "eastEIAcountries",
                                                "origin_point")] #remove constants which aren't used

#-------------------------------------------------------------------------------
# Model Womble: Hierarchical ICAR with boundary identified ----

modelW <- function(seed, d, theta_init, alpha_init, delta_init, init_a, init_b, init_phi_raw, constants, nburnin, thin, niter)
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
      phi[k] <- phi_raw[k] - phi_mean; #manually implement the zero_mean=1 constraint
      a[k] <- phi[k] + x1[k]*beta1 + x2[k]*beta2 + beta0;
      b[k] ~ T(dunif(50, 5000), 50, a[k]); #b[k] ~ dunif(50, 2300) #In each area, start date of occupation, a_k, must be greater than the end date of occupation, b_k (note: BP dates in the positive direction). Truncating in this way avoids dconstraint.
    }
    
    # ICAR Model prior to capture spatial random effects
    phi_raw[1:n_areas] ~ dcar_normal(adj[1:L], weights[1:L], num[1:n_areas], tau1, zero_mean=0)
    phi_mean <- mean(phi_raw[1:n_areas])
    
    
    #For Each Edge Transition
    for (t in 1:n_trans){
      #nabla defines the difference in arrival time across a boundary
      nabla[t] <- abs(a[edge_id1[t]] - a[edge_id2[t]]) #edge t: select first area, m, and second area, n
      #nabla_phi defines the difference in spatial residues across a boundary
      nabla_phi[t] <- abs(phi[edge_id1[t]] - phi[edge_id2[t]])
    }
    
    #Priors
    beta0 ~ dnorm(2500, sd=100); #Intercept
    beta1 ~ dnorm(0, sd=150); #determining strength of the covariate
    beta2 ~ dnorm(0, sd=150);
    tau1 ~dgamma(0.8, 0.1);  #weak prior for ICAR model -- spatial autocorrelation precision parameter #Also try: dgamma(1, 1)
    
    # Hyperprior for duration
    gamma1 ~ dunif(1,20) #Hyperprior for rate
    gamma2 ~ T(dnorm(mean=200,sd=100), 1, 500) #Hyperprior for mode
  })
  
  #Define initial values ---- 
  inits <- list(theta=theta_init,
                a=init_a,
                b=init_b,
                phi_raw=init_phi_raw,
                alpha=alpha_init, 
                delta=delta_init, 
                tau1=rgamma(1, shape = 0.8, rate = 0.1),
                beta0=rnorm(1, 2500, 100),
                beta1=rnorm(1, 0, sd=150), #rnorm(1:constants$n_areas, 0, sd=200), #runif(1:constants$n_areas, min = -1, max = 1)
                beta2=rnorm(1, 0, sd=150),
                gamma1=10,
                gamma2=200)
  
  # Compile and Run model	----
  model <- nimbleModel(model, constants=constants, data=d, inits=inits)
  cModel <- compileNimble(model)
  
  #Configure MCMC with conjugacy where possible
  conf <- configureMCMC(model, useConjugacy = TRUE, control = list(adaptInterval=5000, adaptFactorExponent=0.1))
  
  #Add monitors of importance
  conf$addMonitors(c('a','b','theta','delta','alpha','phi','nabla', 'nabla_phi', 'beta1','beta2','beta0'))
  #Build, compile, and run MCMC
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
     file=here('output','Womble_Emodel_d38.RData')) 
#'Womble_Emodel_d38_reduce.RData' for model run on dataset with Wanyika grade C and D dates filtered out
#'Womble_Emodel_d29.RData' for model running with a sample window with a finer spatial scale (i.e. 83 hexagonal units, d=2.9)