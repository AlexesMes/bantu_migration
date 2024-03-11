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
set.seed(123)

##SCRIPT TO SIMULATE DATA FOR PHASEMODEL WITHOUT ERRORS 

# Load sample window data ----
load(here('data','sample_window.RData'))
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
#sites$area_origin <- hex_area_win$area_center[sites$area_id]


#Simulate multiple observations at each site ----
id_sites <- c(1:n_sites, 
              sample(1:n_sites,
                     size = n_dates-n_sites,
                     replace=TRUE,
                     prob = dexp(1:n_sites,rate=1)/sum(dexp(1:n_sites,rate=1)))) #generate a site id for each date

dates <- sites[1, ] #initialise dates_sf 
for(i in 2:n_dates){
  current_site_id <- id_sites[i]
  dates[i, ] <- sites[current_site_id, ]
}

##CHECK ---
#site_freq  <- plyr::count(dates, 'site_id') ##See how many observations at each site
#table(id_sites)
#area_freq  <- plyr::count(dates, 'area_id') ##See how many observations in each hex area

# #Check that this lines up visually with how many sites are in each hex area
# ggplot(data = hex_area_win) +
#   geom_sf(data = st_buffer(st_as_sf(sampling_win, crs = 4326), 40000), aes(color = "grey50")) + #sampling window with coastal buffer
#   geom_sf() + #hex grid
#   geom_sf_label(aes(label = area_ID)) +
#   geom_sf(data = sites, size=2, alpha=0.5) + #sites
#   #geom_sf(data = hex_area_win$area_center, size=2, alpha=1, aes(color = "purple")) + #hex-origins
#   theme(panel.background = element_rect(fill = "lightblue",
#                                         colour = "lightblue",
#                                         size = 0.5,
#                                         linetype = "solid"),
#         legend.position = "none")


#-------------------------------------------------------------------------------
#Model ----
sim_model <- nimbleCode({
  for (j in 1:n_sites)
  {
    delta[j] ~ dgamma(5,(5-1)/200); #Site duration parameter.
    alpha[j] ~ dunif(max=a, min=b);
    beta[j] <- alpha[j] - (delta[j] + 1); 
  }
  
  for (i in 1:n_dates){
    theta[i] ~ dunif(min=beta[id_sites[i]], max=alpha[id_sites[i]]);
  }
})

#Define constants ----
sim_constants <- list()
sim_constants$n_sites <- n_sites
sim_constants$n_dates  <- n_dates
sim_constants$id_sites  <- id_sites
sim_constants$a <- 3700
sim_constants$b <- 3200

#Simulate ----
set.seed(1223)
simModel <- nimbleModel(code = sim_model, constants = sim_constants)
simModel$simulate('delta')
simModel$simulate('alpha')
simModel$simulate('beta')
simModel$simulate('theta')

# Combine data ----
cra = round(simModel$theta)
sim_df <- list(cra = cra,
               site_id = id_sites)

# Collect site level information ----
sites <- sites %>% 
  st_drop_geometry() %>% 
  left_join(as.data.frame(sim_df), by='site_id')

earliest_dates <- aggregate(cra ~ site_id, data = sites, FUN = max) #Earliest Date for Each Site 
latest_dates <- aggregate(cra ~ site_id, data=sites, FUN=min) #Latest Date for Each Site
n_dates_persite <- aggregate(cra ~ site_id, data=sites, FUN=length) #Number of Dates for Each Site

siteInfo <- data.frame(site_id = earliest_dates$site_id,
                       area_id = sites$area_id,
                       earliest = earliest_dates$cra,
                       latest = latest_dates$cra,
                       diff = earliest_dates$cra - latest_dates$cra,
                       n_dates = n_dates_persite$cra) %>% unique()

#Save constants ----
constants <- list()
constants$dist_mat  <- dist_mat
constants$dist_org  <- dist_org
constants$n_areas  <- nrow(hex_area_win)
constants$n_dates  <- n_dates
constants$n_sites  <- n_sites

#Store output ----
save(sites, siteInfo, sim_df, sim_constants, constants, sampling_win, hex_area_win, file=here('data','tactical_sim_phase_TOY.RData'))