#Load Libraries ----
library(here)
library(coda)
library(nimbleCarbon)
library(rcarbon)

rm(list = ls())
`%!in%` <- Negate(`%in%`)

set.seed(123)

#-------------------------------------------------------------------------------
## Data Setup ----
load(here('data', 'tactical_sim_phase_TOY.RData'))

#-------------------------------------------------------------------------------
## Initialise Parameters ----

# Initialise regional parameters ----
buffer <- 100

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

#-------------------------------------------------------------------------------
#Spatial data

library(spdep)
nb_areas <- poly2nb(as(hex_area_win, 'Spatial'), queen=FALSE, row.names = hex_area_win$area_ID) #neighboring areas using sp library 
#nb_areas <- st_intersects(hex_area_win, hex_area_win, remove_self = TRUE) #neighboring areas using sf library

nbInfo <- nb2WB(nb_areas) #transform into iCAR inputs: adjacent matrix, weights, number of neighbors (for WinBUGS)

#-------------------------------------------------------------------------------
## Model assuming independence of samples ----

model1 <- nimbleCode({
  for (i in 1:n_dates){
    theta[i] ~ dunif(min = b[id_area[i]], max = a[id_area[i]]);
  }
  
  # Set Prior for Each Region
  for (k in 1:n_areas){
    a[k] <- beta_0 - (beta_1 + s1[k])*dist_org[k]; #dist_org[k] can be replaced with some sort of cumulative friction c[k]
    b[k] ~ dunif(50,5000);
    constraint_uniform[k] ~ dconstraint(a[k]>b[k]) #In each area, start date of occupation, a_k, must be greater than the end date of occupation, b_k (note: BP dates in the positive direction)
  }
  #Priors
  beta_0 ~ dnorm(2500, sd=200); #beta_0 #Assume the first migration to be somewhere between 2300BP and 2700BP. Note: age of approximate origin, Katuruka, 2549BP
  beta_1 ~ dexp(1) #beta_1 #If we were focused on the introduction of farming, a sensible prior can be based on known archaeological examples of farming dispersal rates
  
  # ICAR Model Prior
  s1[1:n_areas] ~ dcar_normal(adj[1:L], weights[1:L], num[1:n_areas], tau1, zero_mean =0)
  tau1 <- 1/sigma1^2
  sigma1 ~ dunif(0,100)
  
})

#Constants ----
constants$n_dates <- sim_constants$n_dates
constants$n_areas <- constants$n_areas
constants$adj <- nbInfo$adj
constants$weights <- nbInfo$weights
constants$num <- nbInfo$num
constants$L <- length(nbInfo$adj)

#Define initial values ---- 
d1 <- list(cra=sim_df$cra, unif.const=1)
theta.init = d1$cra

inits1 <- list(a=init_a,
               b=init_b,
               theta=theta.init,
               beta_0=rnorm(1,2500,200),
               beta_1=rexp(1,1),
               s1 = rnorm(constants$n_areas, sd=0.001),
               sigma1=runif(1,0,100))


#Run MCMC ----
mcmc.samples1 <- nimbleMCMC(code = model1,
                           constants = constants,
                           data = d1,
                           niter = 20, 
                           nchains = 3, 
                           thin = 2, 
                           nburnin = 10,
                           monitors = c('a','b','theta'), 
                           inits = inits1, 
                           samplesAsCodaMCMC=TRUE)

#Diagnostics ----
rhat1  <- gelman.diag(mcmc.samples1, multivariate = FALSE)
ess1  <- effectiveSize(mcmc.samples1)

#-------------------------------------------------------------------------------
# Save output ----
save(mcmc.samples1, rhat1, ess1, file=here('output','phasemodel_tactsim_TOY.RData'))

