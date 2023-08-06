#Load Libraries ----
library(here)
library(coda)
library(nimbleCarbon)
library(rcarbon)

#-------------------------------------------------------------------------------
## Data Setup ----
load(here('data', 'tactical_sim_phase.RData'))

#-------------------------------------------------------------------------------
# Model assuming independence of samples ----

model1 <- nimbleCode({
  for (i in 1:Ndates)
  {
    theta[i] ~ dunif(max=a, min=b); #where a and b are the start and end dates of the focus area (switched around because BP dates go in the positive direction for nimbleCarbon)
    # Calibration
    mu[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]); #C14 age on the relevant calibration curve
    sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]); #error on the calibration curve
    error[i] <- (cra_error[i]^2 + sigmaCurve[i]^2)^(1/2); #the samples' C14 error + error on calibration curve
    cra[i] ~ dnorm(mean=mu[i], sd=error[i]); #observed radiocarbon age of the sample
  }
  a ~ dunif(0,10000);
  b ~ dunif(0,10000);
  unif.const ~ dconstraint(b<a); #date a is earlier than date b
})

#Constants ----
data("intcal20") 
constants1 <- list(Ndates = sim_constants$Ndates,
                   calBP = intcal20$CalBP,
                   C14BP = intcal20$C14Age,
                   C14err = intcal20$C14Age.sigma)

#Initialise parameters ---- 
d1 <- list(cra=sim_df$cra, cra_error=sim_df$cra_error, unif.const=1)
theta.init = medCal(calibrate(d1$cra, 
                              d1$cra_error, 
                              verbose = FALSE))

inits1 <- list(a=5000,
               b=500,
               theta=theta.init)


#Run MCMC ----
mcmc.samples1<- nimbleMCMC(code = model1,
                           constants = constants1,
                           data = d1,
                           niter = 2000000, 
                           nchains = 3, 
                           thin=100, 
                           nburnin = 1000000,
                           monitors=c('a','b','theta'), 
                           inits=inits1, 
                           samplesAsCodaMCMC=TRUE)

#Diagnostics ----
rhat1  <- gelman.diag(mcmc.samples1, multivariate = FALSE)
ess1  <- effectiveSize(mcmc.samples1)



#-------------------------------------------------------------------------------
# Model integrating sample interdependence, i.e. the addition of a hierarchical model ----

model2 <- nimbleCode({
  for (j in 1:Nsites)
  {
    delta[j] ~ dgamma(gamma1,(gamma1-1)/gamma2)
    alpha[j] ~ dunif(max=a, min=b);
  }
  
  for (i in 1:Ndates){
    theta[i] ~ dunif(min=(alpha[id_sites[i]] - (delta[id_sites[i]]+1)), max=alpha[id_sites[i]]); 
    #Calibration
    mu[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]);
    sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]);
    error[i] <- (cra_error[i]^2 + sigmaCurve[i]^2)^(1/2);
    cra[i] ~ dnorm(mean=mu[i], sd=error[i]);
  }
  
  # Set Prior for Duration
  a ~ dunif(0,10000);
  b ~ dunif(0,10000);	
  unif.const ~ dconstraint(b<a);
  gamma1 ~ dunif(1,20); #A uniform hyper-prior for the rate (bounded between 1 and 20) 
  gamma2 ~ T(dnorm(mean=200,sd=100), 1, 500); #Hyper-prior for the mode
})

#Constants ----
data("intcal20") 
constants2 <- list(Ndates=sim_constants$Ndates,
                   Nsites=sim_constants$Nsites,
                   calBP=intcal20$CalBP,
                   C14BP=intcal20$C14Age,
                   C14err=intcal20$C14Age.sigma,
                   id_sites=sim_constants$id_sites)

#Initialise parameters ----
d2 <- list(cra=sim_df$cra, cra_error=sim_df$cra_error, unif.const=1)
theta.init  <-  medCal(calibrate(d1$cra,
                                 d1$cra_error,
                                 verbose = FALSE))
dd  <-  data.frame(theta=theta.init, id=constants2$id_sites)
buffer  <- 100
earliest  <- aggregate(theta~id, max, data=dd) #Earliest date at each site
latest  <- aggregate(theta~id, min, data=dd) #Latest date at each site
diff.age  <- earliest$theta - latest$theta
delta.init  <- diff.age + buffer
alpha.init  <- earliest$theta + buffer/2

inits2 <- list(a = 5000,
               b = 500,
               theta = theta.init,
               alpha = alpha.init,
               delta = delta.init,
               gamma1 = 5,
               gamma2 = 200)


#Run MCMC -----
mcmc.samples2<- nimbleMCMC(code = model2, 
                           constants = constants2, 
                           data = d2, 
                           niter = 2000000, 
                           nchains = 3, 
                           thin=100, 
                           nburnin = 1000000, 
                           monitors=c('a','b','theta','delta','alpha','gamma1','gamma2'), 
                           inits=inits2, 
                           samplesAsCodaMCMC=TRUE)

#Diagnostics ----
rhat2  <- gelman.diag(mcmc.samples2, multivariate = FALSE)
ess2  <- effectiveSize(mcmc.samples2)



#-------------------------------------------------------------------------------
# Save output ----
save(mcmc.samples1, rhat1, ess1, mcmc.samples2, rhat2, ess2, file=here('output','phasemodel_tactsim.RData'))
