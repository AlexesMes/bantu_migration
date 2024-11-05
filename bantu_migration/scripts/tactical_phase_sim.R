library(rcarbon)
library(nimbleCarbon)
library(here)

##SCRIPT TO SIMULATE DATA FOR PHASEMODEL WITH ERRORS (i.e. uncalibrated radiocarbon dates with associated errors)

set.seed(1223)
Nsites <- 25
Ndates <- 60
id_sites <- c(1:Nsites, 
              sample(1:Nsites,
                     size=Ndates-Nsites,
                     replace=TRUE,
                     prob=dexp(1:Nsites,rate=1)/sum(dexp(1:Nsites,rate=1))))

#Model ----
sim_model <- nimbleCode({
  for (j in 1:Nsites)
  {
    delta[j] ~ dgamma(5,(5-1)/200); #Site duration parameter.
    alpha[j] ~ dunif(max=a, min=b);
    beta[j] <- alpha[j] - (delta[j] + 1); #The minus sign here is just convention -- NimbleCarbon was written so that the BP dates go in the positive direction. The minus sign here accounts for that. The +1 ensures a minimum where there are two dates at a site with 1 year between them.
  }
  
  for (i in 1:Ndates){
    theta[i] ~ dunif(min=beta[id_sites[i]], max=alpha[id_sites[i]]);
  }
})

#Define constants ----
sim_constants <- list()
sim_constants$Nsites <- Nsites
sim_constants$Ndates  <- Ndates
sim_constants$id_sites  <- id_sites
sim_constants$a <- 3700
sim_constants$b <- 3200

#Simulate ----
simModel <- nimbleModel(code = sim_model, constants = sim_constants)
simModel$simulate('delta')
simModel$simulate('alpha')
simModel$simulate('beta')
simModel$simulate('theta')

# Combine data ----
cra = uncalibrate(round(simModel$theta))$rCRA 
cra_error = rep(20,length(cra))
sim_df <- list(cra = cra,
              cra_error = cra_error,
              site_id = id_sites)

#Store output ----
save(sim_df, sim_constants, file=here('data','tactical_sim_phase.RData'))

