# Load Libraries and Data ----
library(nimbleCarbon)
library(rcarbon)
library(here)
library(parallel)
library(sf)
library(maptools)
library(coda)

#-------------------------------------------------------------------------------
## Data Setup ----
# Load Simulated Data
load(here('data','tactical_sim_gpqr_TOY.RData'))

# Constants ---
constants  <- list()
constants$N <- true_param$n
constants$dist_mat  <- spDists(as_Spatial(sim_sites), longlat = TRUE)
constants$dist_org  <- spDistsN1(as_Spatial(sim_sites), true_param$origin_point, longlat = TRUE)

# Data ---
dat  <- list()
dat$cra  <- sim_sites$cra

# Theta Init ---
theta_init  <-  sim_sites$cra


#-------------------------------------------------------------------------------
# MCMC Function ----
runFun  <- function(seed, dat, theta_init, constants, niter, nburnin, thin)
{
  library(nimbleCarbon)
  library(truncnorm)
  library(cascsim)
  
  # Variance Covariance Function ----
  cov_ExpQ <- nimbleFunction(run = function(dists = double(2), rho = double(0), etasq = double(0), sigmasq = double(0)) 
  {
    returnType(double(2))
    n <- dim(dists)[1]
    result <- matrix(nrow = n, ncol = n, init = FALSE)
    deltaij <- matrix(nrow = n, ncol = n,init = TRUE)
    diag(deltaij) <- 1
    for(i in 1:n)
      for(j in 1:n)
        result[i, j] <- etasq*exp(-0.5*(dists[i,j]/rho)^2) + sigmasq*deltaij[i,j]
    return(result)
  })
  Ccov_ExpQ <- compileNimble(cov_ExpQ)
  assign('cov_ExpQ', cov_ExpQ, envir=.GlobalEnv)
  
  # Handle constraints on dipsersal rate
  dat$lim  <- rep(1, constants$N)
  
  
  # Core model ----
  model <- nimbleCode({
    for (i in 1:N){
      # Model
      rate[i] <- -1/(s[i]-beta1)
      lim[i] ~ dconstraint(rate[i]>0)
      mu[i] <- beta0 + (s[i]-beta1)*dist_org[i]
      theta[i] ~ dAsymLaplace(mu=mu[i], sigma=sigma, tau=tau)
    }
    #priors
    beta0 ~ dnorm(3300, sd=200);
    beta1 ~ dexp(1)
    sigma ~ dexp(0.01)
    etasq ~ dexp(20);
    rho ~ T(dgamma(10, (10-1)/150), 1, 4600); 
    mu_s[1:N] <- 0;
    cov_s[1:N, 1:N] <- cov_ExpQ(dist_mat[1:N, 1:N], rho, etasq, 0.000001)
    s[1:N] ~ dmnorm(mu_s[1:N], cov = cov_s[1:N, 1:N])
  }) 
  
  # MCMC initialisation
  set.seed(seed)
  inits  <-  list()
  inits$theta  <- theta_init
  inits$beta0 <- rnorm(1, 3300, 200)
  inits$beta1 <- rexp(1, rate=2)
  inits$sigma  <- rexp(1, 0.01)
  inits$rho  <- rtgamma(1, shape=10, scale=(10-1)/200, min=1, max=4600)
  inits$etasq  <- rexp(1,20)
  inits$s  <- rep(0, constants$N)
  inits$cov_s <- Ccov_ExpQ(constants$dist_mat, inits$rho, inits$etasq, 0.000001)
  inits$s <-  t(chol(inits$cov_s)) %*% rnorm(constants$N)
  inits$s <- inits$s[ , 1]  # so can give nimble a vector rather than one-column matrix
  
  # Model Compilation
  model.gpqr <- nimbleModel(model, constants = constants, data=dat, inits=inits)
  cModel.gpqr <- compileNimble(model.gpqr)
  
  # MCMC configuration
  conf.gpqr <- configureMCMC(model.gpqr)
  conf.gpqr$addMonitors('s')
  conf.gpqr$addMonitors('rho')
  conf.gpqr$addMonitors('etasq')
  conf.gpqr$removeSamplers('s[1:600]')
  conf.gpqr$removeSamplers('beta1')
  conf.gpqr$addSampler(c('beta1','s[1:600]'), type='AF_slice') 
  MCMC.gpqr <- buildMCMC(conf.gpqr)
  cMCMC.gpqr <- compileNimble(MCMC.gpqr)
  
  # MCMC execution
  results <- runMCMC(cMCMC.gpqr, nchain=1, niter = niter, thin=thin, nburnin = nburnin, samplesAsCodaMCMC = T, setSeed = seed) 
  return(results)
}


## Setup and Execution of MCMC in Parallel ----
ncores <- 8
cl <- makeCluster(ncores)
# Run the model in parallel:
seeds <- c(12, 45, 67, 89, 21, 54, 76, 98)
niter = 8
nburnin = 4
thin = 2
chain_output <- parLapply(cl = cl, X = seeds, fun = runFun, dat = dat, constants = constants, theta = theta_init, niter = niter, nburnin = nburnin,thin = thin)
stopCluster(cl)

# Convert into a mcmc.list object for diagnostic (see below)
gpqr_tactsim <- coda::mcmc.list(chain_output)
rhat <- coda::gelman.diag(gpqr_tactsim, multivariate = F)
range(rhat$psrf[,1])
ii  <- which(rhat$psrf[,1]>1.01)  
ess  <- coda::effectiveSize(gpqr_tactsim)
range(ess)


#-------------------------------------------------------------------------------
## Store Output ----
save(gpqr_tactsim, file=here('output', 'gpqr_tactsim_TOY.RData'))