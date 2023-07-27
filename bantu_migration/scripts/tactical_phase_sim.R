library(rcarbon)
library(nimbleCarbon)
library(here)

set.seed(123)
Nsites <- 10
Ndates <- 30
id.sites <- c(1:Nsites, 
              sample(1:Nsites,
                     size=Ndates-Nsites,
                     replace=TRUE,
                     prob=dexp(1:Nsites,rate=1)/sum(dexp(1:Nsites,rate=1))))

#Model ----
sim.model <- nimbleCode({
  for (k in 1:Nsites)
  {
    delta[k] ~ dgamma(5,(5-1)/200); #Site duration parameter.
    alpha[k] ~ dunif(max=a,min=b);
    beta[k] <- alpha[k] - (delta[k] + 1); #The minus sign here is just convention -- NimbleCarbon was written so that the BP dates go in the positive direction. The minus sign here accounts for that. The +1 ensures a minimum where there are two dates at a site with 1 year between them.
  }
  
  for (i in 1:Ndates){
    theta[i] ~ dunif(beta[id.sites[i]], alpha[id.sites[i]]);
  }
})

#Define constants ----
sim.constants <- list()
sim.constants$Nsites <- Nsites
sim.constants$Ndates  <- Ndates
sim.constants$id.sites  <- id.sites
sim.constants$a <- 3500
sim.constants$b <- 3000

#Simulate ----
set.seed(123)
simModel <- nimbleModel(code = sim.model, constants = sim.constants)
simModel$simulate('delta')
simModel$simulate('alpha')
simModel$simulate('beta')
simModel$simulate('theta')

# Combine data ----
cra = uncalibrate(round(simModel$theta))$rCRA #Can only specify one curve currently for all dates (therefore can't specify calCurves=out.df$calCurve). Default intcal20. TODO: Check if this is problem?
cra_error = rep(20,length(cra))
d.sim <- list(cra = cra,
              cra_error = cra_error,
              id.site = id.sites)

#Store output ----
save(d.sim, sim.constants, file=here('data','tactical_sim_phase.RData'))

