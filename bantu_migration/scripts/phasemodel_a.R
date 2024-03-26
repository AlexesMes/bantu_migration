# Load Library and Data ----
library(here)
library(dplyr)
library(nimbleCarbon)
library(parallel)
library(coda)
library(rcarbon)

`%!in%` <- Negate(`%in%`)

#-------------------------------------------------------------------------------
## Data Setup ----
# Read 14C dates
load(here('data','eastc14.RData'))

# General Setup ----
# Data --
dat <- list(cra = dateInfo$cra,
            cra_error = dateInfo$cra_error,
            constraint_uniform = rep(1, constants$n_areas),
            cra_constraint = rep(1, constants$n_dates)) # Set-up constraint for ignoring inference outside calibration range

# Initial parameters --
buffer <- 100
theta_init <- dateInfo$median_dates
delta_init <- siteInfo$diff + buffer
alpha_init <- siteInfo$earliest + buffer/2


#Calibration curve
constants$cc <- as.numeric(as.factor(dateInfo$calCurve)) #intcal20==1 and shcal20==2

# Dummy extension of the calibration curve
constants$calBP <- c(1000000, constants$calBP, -1000000)
constants$C14BP <- rbind(c(1000000,1000000), constants$C14BP, c(-1000000,-1000000))
constants$C14err <- rbind(c(1000,1000), constants$C14err, c(1000,1000))

# Initialise regional parameters ----
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

#===============================================================================
# MCMC RunScript (Uniform Model a) ----
unif_model_a  <- function(seed, d, theta_init, alpha_init, delta_init, init_a, init_b, constants, nburnin, thin, niter)
{
  #Load Library
  library(nimbleCarbon)
  #Define Core Model
  model <- nimbleCode({
    for (j in 1:n_sites)
    {
      delta[j] ~ dgamma(gamma1, (gamma1-1)/gamma2)
      alpha[j] ~ dunif(max = a[id_area[j]], min = b[id_area[j]]);
    }
    
    for (i in 1:n_dates){
      theta[i] ~ dunif(min = (alpha[id_sites[i]] - (delta[id_sites[i]]+1)), max = alpha[id_sites[i]]);
      #Calibration
      mu[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[ , cc[i]]); #c14age #Index cc selects the correct calibration curve
      cra_constraint[i] ~ dconstraint(mu[i] < 50193 & mu[i] > 95) #C14 age must be within the calibration range
      sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[ , cc[i]]);
      sd[i] <- (cra_error[i]^2 + sigmaCurve[i]^2)^(1/2);
      cra[i] ~ dnorm(mean=mu[i], sd=sd[i]);
    }
    
    # Set Prior for Each Region
    for (k in 1:n_areas){
      a[k] ~ dunif(50,5000);
      b[k] ~ dunif(50,5000);
      constraint_uniform[k] ~ dconstraint(a[k]>b[k]) #In each area, start date of occupation, a_k, must be greater than the end date of occupation, b_k (note: BP dates in the positive direction)
    }
    # Hyperprior for duration
    gamma1 ~ dunif(1,20) #Hyperprior for rate
    gamma2 ~ T(dnorm(mean=200,sd=100), 1, 500) #Hyperprior for mode
  })
  
  # Define Initial values ----
  inits <- list(theta=theta_init, 
                alpha=alpha_init, 
                delta=delta_init, 
                a=init_a, 
                b=init_b)
  inits$gamma1  <- 10
  inits$gamma2  <- 200
  
  # Compile and Run model	----
  model <- nimbleModel(model, constants=constants, data=d, inits=inits)
  cModel <- compileNimble(model)
  conf <- configureMCMC(model, control=list(adaptInterval=20000, adaptFactorExponent=0.1))
  conf$addMonitors(c('theta','delta','alpha'))
  MCMC <- buildMCMC(conf)
  cMCMC <- compileNimble(MCMC)
  results <- runMCMC(cMCMC, niter = niter, thin = thin, nburnin = nburnin, samplesAsCodaMCMC = T, setSeed = seed) 
}

# Run MCMCs ----

# MCMC Setup
ncores  <-  4
cl <- makeCluster(ncores)
seeds <- c(12, 34, 56, 78)
niter  <- 2000000
nburnin  <- 1000000
thin  <- 100

out_unif_model_a  <-  parLapply(cl = cl, 
                                X = seeds, 
                                fun = unif_model_a, 
                                d = dat, 
                                constants = constants, 
                                theta_init = theta_init, 
                                alpha_init = alpha_init, 
                                delta_init = delta_init,  
                                init_a = init_a, 
                                init_b = init_b, 
                                niter = niter, 
                                nburnin = nburnin,
                                thin = thin)

out_unif_model_a <- mcmc.list(out_unif_model_a)


# Diagnostics ----
rhat_unif_model_a <- gelman.diag(out_unif_model_a, multivariate = FALSE)
ess_unif_model_a <- effectiveSize(out_unif_model_a)
agg_unif_model_a <- agreementIndex(dat$cra,
                                   dat$cra_error,
                                   calCurve = dateInfo$calCurve,
                                   theta = out_unif_model_a[[1]][ , grep("theta", colnames(out_unif_model_a[[1]]))],
                                   verbose = F)


#-------------------------------------------------------------------------------
# Save output ----
save(out_unif_model_a, 
     rhat_unif_model_a, 
     ess_unif_model_a, 
     agg_unif_model_a, 
     file=here("output","phase_model_a.RData"))
