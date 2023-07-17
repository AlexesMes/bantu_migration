# Load Libraries ----
library(nimbleCarbon)
library(rcarbon)
library(here)
library(quantreg)
library(parallel)

# Load and prepare data ----
load(here('data','c14.RData'))


#===============================================================================
## Quantile Regression based on median earliest date at each site ----

# Add distance from origin to siteInfo
siteInfo$dist_org <- constants$dist_org

# Compute Quantile Regression (95th percentile)
fit_rq <- rq(earliest ~ dist_org, tau = 0.99, data=siteInfo, alpha=0.95)

# Derive overall estimated rate of dispersal
-1/summary(fit_rq)$coefficients[2,] 
#Note: Between 2.83 to 18.67 km/year, with estimate at 3.32 km/year. 'Lower bd' and 'upper bd' represent the endpoints of confidence intervals for the model coefficients
#Additional note: when removing Russell_EIA data and only including SARD and ADRAC we find a much smaller confidence interval: between 2.53 to 3.77 km/year with an estimate of 3.33 km/year

#===============================================================================
## Bayesian Quantile Regression ----

# Setup Data, Constants, and Inits ----

# Data #Note: consider only the earliest date at each site
subset_dateInfo  <- subset(dateInfo, earliestAtSite == TRUE)
subset_dateInfo  <- subset_dateInfo[order(subset_dateInfo$siteID, decreasing = F),] #Arrange in ascending site ID order
# Generate list of observed data
dat  <- list(cra = subset_dateInfo$cra, 
             cra_error = subset_dateInfo$cra_error, 
             calCurve = as.numeric(as.factor(subset_dateInfo$calCurve))) #intcal20==1 and shcal20==2

# Constants
# Remove constants defined in prepare_data.R that we aren't going to need right now
constants$n_sites <- NULL
constants$id_sites  <- NULL
constants$dist_mat  <- NULL
# Update number of dates, since we now only have the earliest dates in our dataset
constants$n_dates  <- nrow(subset_dateInfo)
# Define Quantile
constants$tau <- 0.99

# Dummy extension of the calibration curves -- 'bookend' values to ensure the regression algorithm never falls out of bounds
constants$calBP <- c(1000000, constants$calBP, -1000000)
constants$C14BP <- rbind(c(1000000,1000000), constants$C14BP, c(-1000000,-1000000))
constants$C14err <- rbind(c(1000,1000), constants$C14err, c(1000,1000))

# constants$iC14BP <- c(1000000,constants$C14BP[,1],-1000000)
# constants$sC14BP <- c(1000000,constants$C14BP[,2],-1000000)
# constants$iC14err <- c(1000000,constants$C14err[,1],-1000000)
# constants$sC14err <- c(1000000,constants$C14err[,2],-1000000)

# Constraint for ignoring inference outside calibration range. Creates a list in the dat called 'constraint' which is filled with 1's used later as an indicator that data is within calibration curve
dat$cra_constraint = rep(1, constants$n_dates) 

# Initilise theta
theta_init  <- subset_dateInfo$median_dates


#-------------------------------------------------------------------------------
## Main Script ----
runFun <- function(seed, dat, theta_init, constants, nburnin, thin, niter)
{
  library(nimbleCarbon)
  model <- nimbleCode({
    for (i in 1:n_dates){
      cc <- dat$calCurve[i] #index selecting correct calibration curve for the site
      # C14BP_cc <- C14BP[ ,cc] #select calibration curve at beginning, instead of using cc index directly on constants$C14BP, because dynamic indexing on constants not working so well with nimble
      # C14err_cc <- C14err[ ,cc]
      #Way less elegant version (without dynamic indexing) which also isn't working...
      # C14BP_cc <- if(cc==1){C14BP_cc <- iC14BP} else {C14BP_cc <- sC14BP}
      # C14err_cc <- if(cc==1){C14err_cc <- iC14BP} else {C14err_cc <- sC14err}
      
      # Model
      mu[i] <- alpha - beta*dist_org[i]
      theta[i] ~ dAsymLaplace(mu=mu[i], sigma=sigma, tau=tau) 
      c14age[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[ ,cc]); #Index cc selects the correct calibration curve
      cra_constraint[i] ~ dconstraint(c14age[i] < 50193 & c14age[i] > 95) #C14 age must be within the calibration range
      sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[ ,cc]);
      sigmaDate[i] <- (cra_error[i]^2+sigmaCurve[i]^2)^(1/2);
      cra[i] ~ dnorm(mean=c14age[i],sd=sigmaDate[i]);
    }
    #priors #TODO: adjust priors... 
    alpha ~ dnorm(3000, sd=200); #beta_0 #Assume the first migration to be somewhere between 3500BP and 2500BP. Note: age of approximate origin, Obobogo, 3070 +- 95BP
    beta ~ dexp(1) #beta_1 #If we were focused on the introduction of farming, a sensible prior can be based on known archaeological examples of farming dispersal rates
    sigma ~ dexp(0.01) #lambda
  })
  set.seed(seed)
  inits  <- list(alpha=rnorm(1,3000,200), beta=rexp(1,1), sigma=rexp(1,0.01), theta=theta_init)
  model_asymlap <- nimbleModel(model, constants = constants, data=dat, inits=inits)

  while(any(model_asymlap$logProb_theta==-Inf))
  {
    inits  <- list(alpha=rnorm(1,3000,200), beta=rexp(1,1), sigma=rexp(1,0.01), theta=theta_init)
    model_asymlap <- nimbleModel(model, constants = constants, data=dat, inits=inits)
  }

  cModel_asymlap <- compileNimble(model_asymlap)
  conf_asymlap <- configureMCMC(model_asymlap, control=list(adaptInterval=20000, adaptFactorExponent=0.1)) #TODO: explore different parameters
  conf_asymlap$removeSampler(c('alpha','beta','sigma'))
  conf_asymlap$addSampler('alpha', control=list(adaptInterval=200, adaptFactorExponent=0.8))
  conf_asymlap$addSampler('beta', control=list(adaptInterval=200, adaptFactorExponent=0.8))
  conf_asymlap$addSampler('sigma', control=list(adaptInterval=200, adaptFactorExponent=0.8))
  conf_asymlap$addMonitors('theta')
  MCMC_asymlap <- buildMCMC(conf_asymlap) #Create an MCMC function from a NIMBLE model
  cMCMC_asymlap <- compileNimble(MCMC_asymlap)
  results <- runMCMC(cMCMC_asymlap, niter = niter, thin=thin, nburnin = nburnin, samplesAsCodaMCMC = T, setSeed = seed) #Can specify number of chains here
  return(results)
}

# Setup and Execution of MCMC in Parallel ----

# Setup parallel processing
ncores  <-  4
cl <- makeCluster(ncores)
seeds  <-  c(12, 34, 56, 78)
niter  <- 20000 #Number of iterations #Working 1000000 and 500000
nburnin  <- 10000 #Burn-in iterations (after this assuming chain has reached a good approximation of stationary)
thin  <- 100 #Parameters sampled every 300 steps

#Run the model in parallel
chain_output = parLapply(cl = cl, X = seeds, fun = runFun, d = dat, constants = constants, theta = theta_init, niter = niter, nburnin = nburnin, thin = thin)
stopCluster(cl)

# Convert into a mcmc_list object for diagnostic (see below)
quantreg_sample <- coda::mcmc.list(chain_output)
qrhat <- coda::gelman.diag(quantreg_sample, multivariate = FALSE)
qess <- coda::effectiveSize(quantreg_sample)


#-------------------------------------------------------------------------------
## Store Output ----
save(qrhat, qess, fit_rq, quantreg_sample, file=here('output','quantreg_res.RData'))

