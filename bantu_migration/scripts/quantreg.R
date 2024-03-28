# Load Libraries ----
library(nimbleCarbon)
library(rcarbon)
library(here)
library(quantreg)
library(parallel)

# Load and prepare data ----
load(here('data', 'eastc14.RData'))


#===============================================================================
## Quantile Regression based on median earliest date at each site ----

# Add distance from origin to siteInfo
siteInfo$dist_org <- constants$dist_org

# Compute Quantile Regression (95th percentile)
fit_rq <- rq(earliest ~ dist_org, tau = 0.99, data=siteInfo, alpha=0.95)

# Derive overall estimated rate of dispersal
-1/summary(fit_rq)$coefficients[2,] 
#Note: Between 1.16 to 11.7 km/year, with estimate at 6.1 km/year. 'Lower bd' and 'upper bd' represent the endpoints of confidence intervals for the model coefficients

#===============================================================================
## Bayesian Quantile Regression ----

# Setup Data, Constants, and Inits ----

# Data #Note: consider only the earliest date at each site
subset_dateInfo  <- subset(dateInfo, earliestAtSite == TRUE)
subset_dateInfo  <- subset_dateInfo[order(subset_dateInfo$siteID, decreasing = F),] #Arrange in ascending site ID order
# Generate list of observed data
dat  <- list(cra = subset_dateInfo$cra, 
             cra_error = subset_dateInfo$cra_error)

# Constants
# Remove constants defined in prepare_data.R that we aren't going to need right now
constants$n_sites <- NULL
constants$id_sites  <- NULL
constants$dist_mat  <- NULL
# Update number of dates, since we now only have the earliest dates in our dataset
constants$n_dates  <- nrow(subset_dateInfo)
# Define Quantile
constants$tau <- 0.99

#Calibration curve
constants$cc <- as.numeric(as.factor(subset_dateInfo$calCurve)) #intcal20==1 and shcal20==2

# Dummy extension of the calibration curves -- 'bookend' values to ensure the regression algorithm never falls out of bounds
constants$calBP <- c(1000000, constants$calBP, -1000000)
constants$C14BP <- rbind(c(1000000,1000000), constants$C14BP, c(-1000000,-1000000))
constants$C14err <- rbind(c(1000,1000), constants$C14err, c(1000,1000))

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
      # Model
      mu[i] <- alpha - beta*dist_org[i]
      theta[i] ~ dAsymLaplace(mu=mu[i], sigma=sigma, tau=tau) 
      c14age[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[ ,cc[i]]); #Index cc selects the correct calibration curve
      cra_constraint[i] ~ dconstraint(c14age[i] < 50193 & c14age[i] > 95) #C14 age must be within the calibration range
      sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[ ,cc[i]]);
      sigmaDate[i] <- (cra_error[i]^2 + sigmaCurve[i]^2)^(1/2);
      cra[i] ~ dnorm(mean=c14age[i],sd=sigmaDate[i]);
    }
    #priors 
    alpha ~ dnorm(2500, sd=200); #beta_0 #Assume the first migration to be somewhere between 2400BP and 2600BP. Note: age of approximate origin, Katuruka, 2570 BP
    beta ~ dexp(1) #beta_1 #If we were focused on the introduction of farming, a sensible prior can be based on known archaeological examples of farming dispersal rates
    sigma ~ dexp(0.01) #lambda
  })
  set.seed(seed)
  inits  <- list(alpha=rnorm(1,3500,200), beta=rexp(1,1), sigma=rexp(1,0.01), theta=theta_init)
  model_asymlap <- nimbleModel(model, constants = constants, data=dat, inits=inits)

  while(any(model_asymlap$logProb_theta==-Inf))
  {
    inits  <- list(alpha=rnorm(1,3500,200), beta=rexp(1,1), sigma=rexp(1,0.01), theta=theta_init)
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
niter  <- 2000000 #Number of iterations
nburnin  <- 1000000 #Burn-in iterations (after this assuming chain has reached a good approximation of stationary)
thin  <- 100 #Parameters sampled every 100 steps

#Run the model in parallel
chain_output = parLapply(cl = cl, X = seeds, fun = runFun, d = dat, constants = constants, theta = theta_init, niter = niter, nburnin = nburnin, thin = thin)
stopCluster(cl)

# Convert into a mcmc_list object for diagnostic (see below)
quantreg_sample <- coda::mcmc.list(chain_output)
qrhat <- coda::gelman.diag(quantreg_sample, multivariate = FALSE)
qess <- coda::effectiveSize(quantreg_sample)

#Agreement index
agree_quantreg <- agreementIndex(dat$cra,
                                 dat$cra_error,
                                 calCurve = subset_dateInfo$calCurve, 
                                 theta = quantreg_sample[[1]][ , grep("theta", colnames(quantreg_sample[[1]]))] , #a Matrix containing the posterior samples of each date
                                 verbose = F) #For chain 1 -- can do this for other chains too
#mean(agree_quantreg$agreement)
##Agreement index of Mubuga (siteID=389), Kabacusi (167), and Mucucu (390)
#agree_quantreg$agreement[389]
#agree_quantreg$agreement[167]
#agree_quantreg$agreement[390]

#-------------------------------------------------------------------------------
## Store Output ----
save(qrhat, qess, fit_rq, quantreg_sample, file=here('output','quantreg_res.RData'))

