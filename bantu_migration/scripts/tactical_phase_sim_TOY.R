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

set.seed(123)

source(here('src','hex_areas.R'))

##SCRIPT TO SIMULATE DATA WITHOUT ERRORS 

#-------------------------------------------------------------------------------
## List of countries in sub-Saharan Africa ----
subSahara_countries <- c("South Africa", 
                         "Lesotho", 
                         "eSwatini", 
                         "Botswana",
                         "Zimbabwe",
                         "Namibia",
                         "Angola",
                         "Zambia",
                         "Mozambique",
                         "Malawi",
                         "Madagascar",
                         "Tanzania",
                         "Rwanda",
                         "Burundi",
                         "Kenya",
                         "Uganda",
                         "Somalia",
                         "Ethiopia",
                         "Central African Republic",
                         "Cameroon",
                         "Democratic Republic of the Congo",
                         "Republic of the Congo",
                         "Gabon",
                         "Cameroon",
                         "Nigeria",
                         "Equatorial Guinea",
                         "Sudan",
                         "South Sudan",
                         "Chad")
#-------------------------------------------------------------------------------
# Generate Spatial Window for Analyses: Sub-Saharan Africa ----
sampling_win <- ne_countries(continent = "Africa", returnclass = "sf") %>%
  filter_all(., any_vars(str_detect(., "Sub-Saharan"))) %>% 
  filter(name_en %in% subSahara_countries) %>%
  filter(name_en != "Madagascar") #We focus on mainland sub-Saharan Africa

#Generate Hex Areas over Spatial Window ----
hex_area_win <- hex_areas(sampling_win, cell_d = 7.5)

#Spatial Window for Analyses ----
sf::sf_use_s2(FALSE) #turn off spherical co-ordinates
sampling_win <-  sampling_win %>%
  st_make_valid() %>%
  st_union()
sf::sf_use_s2(TRUE) #turn on spherical co-ordinates

#-------------------------------------------------------------------------------
# Target Parameters ----
n_sites  <- 10
n_dates  <- 30 
origin_point <- st_point(c(11.40, 5.48)) #dispersal origin point -- approximately at Ngoume
# true_param$beta0 <- 3300 #approximate mean date at origin point
# true_param$beta1 <- 0.4 #reciprocal of dispersal rate 
# true_param$sigma <- 100 


#-------------------------------------------------------------------------------
#Simulate Data ----

#Generate sites and calculate distances from origin ---
sites <- st_sample(sampling_win,  size = n_sites, type = 'random')
st_crs(sites) <- 4326
dist_mat  <- st_distance(sites)/1000
#dist_org  <-  st_distance(x=sites, y=origin_point)/1000 


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
# site_freq  <- plyr::count(dates_sf, 'site_id') ##See how many observations at each site
# table(id_sites)
# area_freq  <- plyr::count(dates_sf, 'area_id') ##See how many observations in each hex area

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
sim_constants$a <- 3500
sim_constants$b <- 3000

#Simulate ----
set.seed(123)
simModel <- nimbleModel(code = sim_model, constants = sim_constants)
simModel$simulate('delta')
simModel$simulate('alpha')
simModel$simulate('beta')
simModel$simulate('theta')

# Combine data ----
cra = round(simModel$theta)
sim_df <- list(cra = cra,
               site_id = id_sites)

#Store output ----
save(sim_df, sim_constants, file=here('data','tactical_sim_phase_TOY.RData'))