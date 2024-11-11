# Load Libraries ----
library(here)
library(dplyr)
library(stringr)
library(nimbleCarbon)
library(rnaturalearth)
library(sf)
library(ggplot2)
library(viridis)
library(rcarbon)
library(units)

rm(list = ls())
`%!in%` <- Negate(`%in%`)
set.seed(123)

##SCRIPT TO SIMULATE DATA (WITH UNDERLYING SPATIAL AUTOCORRELATED STRUCTURE) FOR ICAR MODEL WITHOUT ERRORS

# Load sample window data ----
load(here('data','sample_window.RData'))
load(here('data','trig.RData')) #nodes and edges between hex area centroids

#-------------------------------------------------------------------------------
## List of countries ----
subSahara_countries <- constants_sw$countries #sub-Saharan Africa
eastEIA_countries <- constants_sw$eastEIAcountries #Eastern Sub-Saharan Africa 

#-------------------------------------------------------------------------------
#Spatial Window for Analyses ----
sf::sf_use_s2(FALSE) #turn off spherical co-ordinates
sampling_win <-  sampling_win %>%
  st_make_valid() %>%
  st_union()
sf::sf_use_s2(TRUE) #turn on spherical co-ordinates

#-------------------------------------------------------------------------------
# Target Parameters ----
n_sites  <- 100
n_dates  <- 500 
origin_point <- st_sfc(st_point(c(-1.45, 31.77))) #dispersal origin point -- approximately at Katuruka

#-------------------------------------------------------------------------------
#Spatial data ----

library(spdep)
nb_areas <- poly2nb(as(hex_area_win, 'Spatial'), queen=FALSE, row.names = hex_area_win$area_ID) #neighboring areas using sp library 
nbInfo <- nb2WB(nb_areas) #transform into iCAR inputs: adjacent matrix, weights, number of neighbors (for WinBUGS)

#-------------------------------------------------------------------------------
#Constants ----
#Combine constants
constants <- constants_trig

constants$adj <- nbInfo$adj
constants$weights <- nbInfo$weights
constants$num <- nbInfo$num
constants$L <- length(nbInfo$adj)

#-------------------------------------------------------------------------------
#Simulate Data ----

#Generate sites and calculate distances from origin ---
sites <- st_sample(sampling_win,  size = n_sites, type = 'random')
st_crs(sites) <- 4326
st_crs(origin_point) <- 4326
dist_mat  <- set_units(st_distance(sites), 'km')
dist_org  <- as.vector(set_units(st_distance(x=sites, y=origin_point), 'km'))

#Assign hex area id to each site ----
sites <- st_as_sf(sites) %>% rename(geometry = x)
sites$site_id <- row_number(sites)
sites$area_id <- as.integer(st_within(sites$geometry, hex_area_win$geometry))

#Simulate multiple observations at each site ----
id_sites <- c(1:n_sites, 
              sample(1:n_sites,
                     size = n_dates-n_sites,
                     replace=TRUE,
                     prob = dnorm(1:n_sites, mean = mean(1:n_sites), sd = sd(1:n_sites)) / sum(dnorm(1:n_sites, mean = mean(1:n_sites), sd = sd(1:n_sites))))) #generate a site id for each date

dates <- sites[1, ] #initialise dates_sf 
for(i in 2:n_dates){
  current_site_id <- id_sites[i]
  dates[i, ] <- sites[current_site_id, ]
}

##CHECK ---
site_freq  <- plyr::count(dates, 'site_id') ##See how many observations at each site
area_freq  <- plyr::count(dates, 'area_id') ##See how many observations in each hex area

#-------------------------------------------------------------------------------
##MODEL ---
sim_model <- nimbleCode({
  for (k in 1:n_areas){
    b[k] ~ dunif(100,5000);
    constraint_uniform[k] ~ dconstraint(a[k]>b[k]);
  }
  
  a[1:n_areas] ~ dcar_normal(adj[1:L], weights[1:L], num[1:n_areas], tau1, zero_mean =0)
  tau1 ~ dgamma(2, 0.5)
  
  for (j in 1:n_sites)
  {
    delta[j] ~ dgamma(5,(5-1)/100); #Site duration parameter.
    alpha[j] ~ dunif(max=a[id_areas[j]], min=b[id_areas[j]]);
    beta[j] <- alpha[j] - (delta[j] + 1); #The +1 ensures at a minimum where there are two dates at a site there will be 1 year between them.
    constraint_duration[j] ~ dconstraint(alpha[j]>(delta[j]+1)); #Site can't have have a duration longer than its time of first arrival
  }
  
  for (i in 1:n_dates){
    theta[i] ~ dunif(min=beta[id_sites[i]], max=alpha[id_sites[i]]);
    cra_constraint[i] ~ dconstraint(theta[i] > 0);
  }
})


#Define constants ----
sim_constants <- constants
sim_constants$n_sites <- n_sites
sim_constants$n_dates  <- n_dates
sim_constants$n_areas  <- constants_sw$n_areas
sim_constants$id_sites  <- dates$site_id
sim_constants$id_areas <- sites$area_id

#Define constraints, data, and initial values ----
dat <- list(constraint_uniform = rep(1, sim_constants$n_areas),
            constraint_duration = rep(1, n_sites),
            cra_constraint = rep(1, n_dates))

init_a <- runif(1:sim_constants$n_areas, min = 600, max = 3500)
init_b <- init_a - runif(1:sim_constants$n_areas, min = 50, max = 600)
inits <- list(a = init_a,
              b = init_b,
              tau1 = 2)

#Simulate ----
set.seed(1223)
simModel <- nimbleModel(code = sim_model, constants = sim_constants, data = dat, inits = inits)
simModel$simulate('delta')
simModel$simulate('alpha')
simModel$simulate('beta')
simModel$simulate('theta')

# Combine data ---- ##Note: no model uncertainty has yet been added in... generates dates, not uncalibrated radiocarbon dates with associated error (for this see tactical_ICAR_sim_uncert.R)
cra = round(simModel$theta)
sim_df <- list(cra = cra,
               site_id = id_sites)

# Collect site level information ----
sites_sf <- sites #save a copy
sites <- sites %>% 
  st_drop_geometry() %>% 
  left_join(as.data.frame(sim_df), by='site_id')

earliest_dates <- aggregate(cra ~ site_id, data = sites, FUN = max) #Earliest Date for Each Site 
latest_dates <- aggregate(cra ~ site_id, data=sites, FUN=min) #Latest Date for Each Site
n_dates_persite <- aggregate(cra ~ site_id, data=sites, FUN=length) #Number of Dates for Each Site

siteInfo <- data.frame(site_id = earliest_dates$site_id,
                       earliest = earliest_dates$cra,
                       latest = latest_dates$cra,
                       diff = earliest_dates$cra - latest_dates$cra,
                       n_dates = n_dates_persite$cra) %>% unique()

#Assign hex area id to each site ----
siteInfo$area_id <- as.integer(st_within(sites_sf$geometry, hex_area_win$geometry))

#Save constants ----
constants <- sim_constants[names(sim_constants) %!in% c("a", "b")]
constants$dist_mat  <- dist_mat
constants$dist_org  <- dist_org
constants$n_areas  <- nrow(hex_area_win)
constants$true_a <- init_a
constants$true_b <- init_b
constants$true_alpha <- simModel$alpha
constants$true_beta <- simModel$beta

#Store output ----
save(sites, siteInfo, sim_df, constants, sampling_win, hex_area_win, file=here('data','tactical_sim_ICAR_spatial_auto.RData'))

