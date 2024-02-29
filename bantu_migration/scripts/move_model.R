# Load Library and Data ----
library(here)
library(dplyr)
library(nimbleCarbon)
library(parallel)
library(coda)
library(rcarbon)
library(deldir)

rm(list = ls())

`%!in%` <- Negate(`%in%`)

#-------------------------------------------------------------------------------
## Data Setup ----
# Read 14C dates
load(here('data','c14.RData'))
load(here('data','trig.RData'))


##General Setup ----

#Filter for test case -- #TODO: remove when considering whole dataset
siteInfo <- siteInfo %>% filter(area_id %in% c(13, 14, 18)) %>% mutate(area_id = case_when(area_id == 13 ~ 1, area_id == 14 ~ 2, area_id == 18 ~3))
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

#Dummy extension of the calibration curve
constants$calBP <- c(1000000, constants$calBP, -1000000)
constants$C14BP <- rbind(c(1000000,1000000), constants$C14BP, c(-1000000,-1000000))
constants$C14err <- rbind(c(1000,1000), constants$C14err, c(1000,1000))


#Initial parameters --
buffer <- 100
theta_init <- dateInfo$median_dates

#Initialise regional parameters
#Initialise hex areas which contain sites
init_a  <- aggregate(earliest~area_id, FUN=max, data=siteInfo) #parameter of earliest date in each region k

#Duration parameter initialisation
delta_init <- 0 
for (t in 1:constants$n_trans){ 
  m <- constants$transitions[[t,'region1_id']] #transition in row t, select first area
  n <- constants$transitions[[t,'region2_id']] #transition in row t, select second area
  
  delta_init[t] <- abs(init_a$earliest[init_a$area_id==m] - init_a$earliest[init_a$area_id==n]) + buffer
}

#Add buffer
init_a  <- init_a[ ,2] + buffer

#===============================================================================
# MCMC RunScript (Model assuming independence of samples, no site hierarchy) ----

#Define Core Model
model <- nimbleCode({
  for (i in 1:n_dates){
    theta[i] ~ dunif(min = (a[id_area[id_sites[i]]] - (delta[id_area[id_sites[i]]]+1)), max = a[id_area[id_sites[i]]]);
    #Calibration
    mu[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[ , cc[i]]); #c14age #Index cc selects the correct calibration curve
    cra_constraint[i] ~ dconstraint(mu[i] < 50193 & mu[i] > 95) #C14 age must be within the calibration range
    sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[ , cc[i]]);
    sd[i] <- (cra_error[i]^2 + sigmaCurve[i]^2)^(1/2);
    cra[i] ~ dnorm(mean=mu[i], sd=sd[i]);
  }

  for (k in 1:n_areas)
  {
	  a[k] ~ dunif(50,5000)
  }
  
  for (t in 1:constants$n_trans){ 
    
    #Duration parameter 
    kappa[t] ~ dbinom(size=1,prob=eta[t]); #kappa is a binomial with probability eta[5] of being 1  
    eta[t] ~ dunif(0,1)
    lambda[t] ~ dexp(1) #possibly replace with gamma with hyperprior at later stage
    delta[t] <- constants$transitions[t,7]/((-1 + 2 *kappa[t])*lambda[t]) #time = distance/speed, and transitions[t,'distance'] is the distance between the two areas #Normalised
  }
  # Hyperpriors
  # iota ~ dexp(1); #speed -- from pilot study estimate is between 3.5 and 4 km/year. speed is always positive
})

# Define Initial values ----
inits <- list(theta=theta_init, 
              a=init_a,
              delta=delta_init)
inits$lambda <- rexp(1, rate=1/3)
inits$kappa <- rbeta(seq(0,1, by=1/constants$n_trans), 5, 5)


# Compile and Run model	----
model <- nimbleModel(model, constants=constants, data=dat, inits=inits)
cModel <- compileNimble(model)
conf <- configureMCMC(model, control=list(adaptInterval=20000, adaptFactorExponent=0.1))
conf$addMonitors(c('theta', 'lambda', 'kappa', 'delta'))
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
