# Load Library and Data ----
library(here)
library(dplyr)
library(nimbleCarbon)
library(parallel)
library(coda)
library(rcarbon)
library(deldir)

`%!in%` <- Negate(`%in%`)

#-------------------------------------------------------------------------------
## Data Setup ----
# Read 14C dates
load(here('data','c14.RData'))
load(here('data','trig.RData'))


##General Setup ----

#Filter for test case --
siteInfo <- siteInfo %>% filter(area_id %in% c(13, 14, 18))
dateInfo <- dateInfo %>% filter(siteID %in% siteInfo$siteID)
n_areas <- 3
n_dates <- nrow(dateInfo)

#Data --
dat <- list(cra = dateInfo$cra,
            cra_error = dateInfo$cra_error,
            constraint_uniform = rep(1, n_areas),
            cra_constraint = rep(1, n_dates)) # Set-up constraint for ignoring inference outside calibration range

#Calibration curve --
constants$cc <- as.numeric(as.factor(dateInfo$calCurve)) #intcal20==1 and shcal20==2

# Dummy extension of the calibration curve
constants$calBP <- c(1000000, constants$calBP, -1000000)
constants$C14BP <- rbind(c(1000000,1000000), constants$C14BP, c(-1000000,-1000000))
constants$C14err <- rbind(c(1000,1000), constants$C14err, c(1000,1000))


#Initial parameters --
buffer <- 100
theta_init <- dateInfo$median_dates

#Initialise regional parameters
#Initialise hex areas which contain sites
init_a  <- init_b <- aggregate(earliest~area_id, FUN=max, data=siteInfo) #two parameters of earliest date in each region k

#Add buffer
init_a  <- init_a[ ,2] + buffer
init_b  <- init_b[ ,2] - buffer


#===============================================================================
# MCMC RunScript (Model assuming independence of samples, no site hierarchy) ----

#Define Core Model
model <- nimbleCode({
  for (i in 1:n_dates){
    theta[i] ~ dunif(min = b[id_area[id_sites[i]]], max = a[id_area[id_sites[i]]]);
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
  }
  
  for (t in 1:nrow(transitions)){  
    m <- transitions[t,1] #transition in row t, select first area
    n <- transitions[t,2]  #transition in row t, select second area
    
    #Duration parameter 
    delta[t] <- dist_org[m,n]/(kappa*lambda) #time = distance/speed #TODO: make speed regional -- add indices
  }
  # Hyperpriors
  lambda ~ dexp(1) #speed -- from pilot study estimate is between 3.5 and 4 km/year. speed is always positive
  kappa[1:nrow(transitions)] ~ dbinom();
})

# Define Initial values ----
inits <- list(theta=theta_init, 
              a=init_a, 
              b=init_b)
inits$lambda <- rexp(1, rate=1/3)
inits$kappa <- rbinom(1, 3.5, 0.5)


# Compile and Run model	----
model <- nimbleModel(model, constants=constants, data=dat, inits=inits)
cModel <- compileNimble(model)
conf <- configureMCMC(model, control=list(adaptInterval=20000, adaptFactorExponent=0.1))
conf$addMonitors(c('theta'))
MCMC <- buildMCMC(conf)
cMCMC <- compileNimble(MCMC)

niter  <- 20
nburnin  <- 10
thin  <- 2

results <- runMCMC(cMCMC, niter = niter, thin = thin, nburnin = nburnin, samplesAsCodaMCMC = T, setSeed = 12) 

out_move_model <- mcmc.list(results)


# Diagnostics ----
rhat_move_model <- gelman.diag(out_move_model, multivariate = FALSE)
ess_move_model <- effectiveSize(out_move_model)
agg_move_model <- agreementIndex(dat$cra,
                                 dat$cra_error,
                                 calCurve = dateInfo$calCurve,
                                 theta = out_move_model[[1]][ , grep("theta", colnames(out_move_model[[1]]))],
                                 verbose = F)


#-------------------------------------------------------------------------------
# Save output ----
save(out_move_model, 
     rhat_move_model, 
     ess_move_model, 
     agg_move_model, 
     file=here("output","move_model.RData"))