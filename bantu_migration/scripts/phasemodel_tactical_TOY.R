#Load Libraries ----
library(here)
library(coda)
library(nimbleCarbon)
library(rcarbon)
library(dplyr)

rm(list = ls())
`%!in%` <- Negate(`%in%`)

set.seed(123)

#-------------------------------------------------------------------------------
## Data Setup ----
load(here('data', 'tactical_sim_phase_TOY.RData'))
load(here('data','trig.RData')) #nodes and edges between hex area centroids

#Combine constants
constants <- c(constants, constants_trig)
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

#Gradient parameter initialisation
init_nabla <- 0 
for (t in 1:constants$n_trans){ 
  m <- constants$transitions[[t,'region1_id']] #transition in row t, select first area
  n <- constants$transitions[[t,'region2_id']] #transition in row t, select second area
  
  init_nabla[t] <- abs(init_a$earliest[init_a$area_id==m] - init_a$earliest[init_a$area_id==n])/constants$transitions[[t, 'distance']]
}

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
    theta[i] ~ dunif(min = b[id_areas[i]], max = a[id_areas[i]]);
  }
  
  #For Each Region
  for (k in 1:n_areas){
    b[k] ~ dunif(50,5000);
    constraint_uniform[k] ~ dconstraint(a[k]>b[k]) #In each area, start date of occupation, a_k, must be greater than the end date of occupation, b_k (note: BP dates in the positive direction)
  }
  
  #For each edge transition
  for (t in 1:n_trans){
    #nabla defines the gradient along the edge
    nabla[t] <- abs(a[edge_id1[t]] - a[edge_id2[t]])/edge_dist[t] #edge t: select first area, m, and second area, n
  }

  # ICAR Model Prior
  a[1:n_areas] ~ dcar_normal(adj[1:L], weights[1:L], num[1:n_areas], tau1, zero_mean =0)
  tau1 <- 1/sigma1^2
  sigma1 ~ dunif(0,100)
  
})

#Constants ----
constants$adj <- nbInfo$adj
constants$weights <- nbInfo$weights
constants$num <- nbInfo$num
constants$L <- length(nbInfo$adj)
constants <- constants[names(constants) %!in% c("dist_mat", "dist_org", "center_coords")] #remove constants which aren't used
#transform transitions into usable format
edge_info <- as.data.frame(constants$transitions)
constants$edge_id1 <- edge_info$region1_id
constants$edge_id2 <- edge_info$region2_id 
constants$edge_dist <- edge_info$distance
constants <- constants[names(constants) %!in% c("transitions")]

#Define initial values ---- 
d1 <- list(theta=sim_df$cra, 
           constraint_uniform = rep(1, constants$n_areas)) 


inits1 <- list(a=init_a,
               b=init_b,
               nabla=init_nabla,
               sigma1=runif(1,0,100))


#Run MCMC ----
mcmc.samples1 <- nimbleMCMC(code = model1,
                           constants = constants,
                           data = d1,
                           niter = 2000000, 
                           nchains = 4, 
                           thin=100, 
                           nburnin = 1000000,
                           monitors = c('a', 'b', 'nabla', 'theta'),
                           inits = inits1, 
                           samplesAsCodaMCMC=TRUE)

#Diagnostics ----
rhat1  <- gelman.diag(mcmc.samples1, multivariate = FALSE)
ess1  <- effectiveSize(mcmc.samples1)


#-------------------------------------------------------------------------------
# Save output ----
save(mcmc.samples1, rhat1, ess1, file=here('output','phasemodel_tactsim_TOY.RData'))

#===============================================================================
#===============================================================================
##Plot ----
library(graphics)
library(latex2exp)

# ## Traceplot of start and end of occupation (a, b) 
# #Load Data ----
load(here("output", "phasemodel_tactsim_TOY.RData"))

par(mfrow=c(2,2))
traceplot(mcmc.samples1[,'a[18]'], smooth=TRUE) #region area 18
traceplot(mcmc.samples1[,'b[18]'], smooth=TRUE)
dev.off()

#-------------------------------------------------------------------------------
## Tactical Simulation Posterior Predictive Check for a and b in a given region

#For model (i) select parameters a and b (i.e. start and end date of occupation in the region)
post.model.i  <- do.call(rbind, mcmc.samples1)[ , c(18, 59)] #area 18 (selecting a[18] and b[18])

dens.i.a  <- density(post.model.i[,1],bw = 5)
dens.i.b  <- density(post.model.i[,2],bw=5)

plot(NULL, xlim=c(4100,2700), ylim=c(0,0.022), xlab='Cal BP', ylab='Posterior Probability')
polygon(c(dens.i.a$x, rev(dens.i.a$x)), c(rep(0,length(dens.i.a$x)), rev(dens.i.a$y)), border=NA, col=rgb(0,0.4,0,0.5))
polygon(c(dens.i.b$x, rev(dens.i.b$x)), c(rep(0,length(dens.i.b$x)), rev(dens.i.b$y)), border=NA, col=rgb(0,0.4,0,0.5))
abline(v=c(3700, 3200),lty=2)
axis(3,at=c(3700, 3200),labels=c(TeX('$a$'),TeX('$b$')))
legend('topright', legend=c('Non hierarchichal'), fill=c('darkgreen'))

#-------------------------------------------------------------------------------

