# Load Library and Data ----
library(here)
library(dplyr)
library(nimbleCarbon)
library(parallel)
library(coda)
library(rcarbon)

`%!in%` <- Negate(`%in%`)

#PHASE MODEL WITH CONSTRAINTS
#-------------------------------------------------------------------------------
## Data Setup ----
# Read 14C dates
load(here('data','c14.RData'))

# General Setup ----
# Data --
dat <- list(cra = dateInfo$cra,
            cra_error = dateInfo$cra_error,
            constraint_uniform = rep(1, constants$n_areas),
            cra_constraint = rep(1, constants$n_dates), # Set-up constraint for ignoring inference outside calibration range
            constraint_dispersal=1)

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
# MCMC RunScript (Uniform Model b) ----

unif.model.b <- function(seed, d, theta_init, alpha_init, delta_init, constants, init_a, init_b, nburnin, thin, niter)
{
  #Load Library
  library(nimbleCarbon)
  #Define Core Model
  model <- nimbleCode({
    # Model Start/End Dates at Individual Sites	
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
      constraint_uniform[k] ~ dconstraint(a[k]>b[k])
    }
    # Hyperprior for duration
    gamma1 ~ dunif(1,20)
    gamma2 ~ T(dnorm(mean=200,sd=100),1,500)
    
    # Define Dispersal Constraint
    constraint_dispersal ~ dconstraint(a[34]>a[27] & a[34]>a[26] & 
                                         a[34]>a[21] & a[34]>a[22] &
                                         a[27]>a[16] & a[27]>a[17] &
                                         a[22]>a[12] &
                                         a[34]>a[35] & a[34]>a[43] &
                                         a[34]>a[28] &
                                         a[28]>a[23] &  
                                         a[14]>a[9] &
                                         a[9]>a[5] &
                                         a[29]>a[25] & a[29]>a[30] &
                                         a[30]>a[31] &
                                         a[28]>a[18] & a[28]>a[19] &
                                         a[18]>a[9] & a[18]>a[14]) 
    # Constraints which are not implemented: a[28]>a[29] & a[29]>a[24] #to allow for potential leapfrogging ... 
    #                                        a[6]>a[3] & a[5]>a[3] & 
    #                                        a[14]>a[10] & 
    #                                        a[9]>a[6] &
    #                                        a[18]>a[13]
    
  })
  # Define Inits
  inits <- list(theta=theta_init, 
                alpha=alpha_init, 
                delta=delta_init, 
                a=init_a, 
                b=init_b)
  inits$gamma1  <- 10
  inits$gamma2  <- 200
  
  #Add constraint
  d$constraint_dispersal  <- 1
  
  #Compile and Run
  model <- nimbleModel(model, constants = constants, data=d, inits=inits)
  cModel <- compileNimble(model)
  conf <- configureMCMC(model, control = list(adaptInterval=20000, adaptFactorExponent=0.1))
  conf$addMonitors(c('theta','delta','alpha'))
  MCMC <- buildMCMC(conf)
  cMCMC <- compileNimble(MCMC)
  results <- runMCMC(cMCMC, niter = niter, thin=thin, nburnin = nburnin, samplesAsCodaMCMC = T, setSeed=seed) 
  return(results)
}

# Run MCMCs ----

# MCMC Setup
ncores  <-  4
cl <- makeCluster(ncores)
seeds <- c(12,34,56,78)
niter  <- 6000000
nburnin  <- 3000000
thin  <- 300

out.unif.model_b  <-  parLapply(cl = cl, 
                                X = seeds, 
                                fun = unif.model.b, 
                                d = dat,
                                constants = constants, 
                                theta_init = theta_init, 
                                alpha_init = alpha_init, 
                                delta_init = delta_init,  
                                niter = niter, 
                                init_a = init_a, 
                                init_b = init_b, 
                                nburnin = nburnin,
                                thin = thin)

out.unif.model_b <- mcmc.list(out.unif.model_b)

# Diagnostics ----

rhat.unif.model_b <- gelman.diag(out.unif.model_b, multivariate = FALSE)
ess.unif.model_b <- effectiveSize(out.unif.model_b)
a.unif.model_b <- agreementIndex(dat$cra,
                                 dat$cra_error,
                                 calCurve = dateInfo$calCurve,
                                 theta = out.unif.model_b[[1]][ , grep("theta", colnames(out.unif.model_b[[1]]))], verbose = F)

#-------------------------------------------------------------------------------
# Save output ----
save(out.unif.model_b,
     rhat.unif.model_b,
     ess.unif.model_b,
     a.unif.model_b,file = here("output", "phase_model_b.RData"))

