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

##SCRIPT TO SIMULATE DATA (WITH UNDERLYING SPATIAL AUTOCORRELATED STRUCTURE AND COVARIATE EFFECT) FOR WOMBLING MODEL WITHOUT ERRORS

# Load sample window data ----
load(here('data','sample_window_cont.RData'))
load(here('data','trig_cont.RData')) #nodes and edges between hex area centroids

#-------------------------------------------------------------------------------
## List of countries ----
subSahara_countries <- constants_sw$countries #sub-Saharan Africa

#-------------------------------------------------------------------------------
#Spatial Window for Analyses ----
sf::sf_use_s2(FALSE) #turn off spherical co-ordinates
sampling_win_outline <- sampling_win %>%
  st_make_valid() %>%
  st_union()

sampling_win <- hex_area_win %>%
  st_make_valid() %>%
  st_union()
sf::sf_use_s2(TRUE) #turn on spherical co-ordinates


#-------------------------------------------------------------------------------
# Target Parameters ----
n_sites  <- 100
n_dates  <- 500 
origin_point <- st_sfc(st_point(c(11.4, 5.483))) #dispersal origin point -- approximately at Katuruka, st_point(c(-1.45, 31.77)) (east Africa) or Ngoume (sub-Saharan Africa)

#-------------------------------------------------------------------------------
#Spatial data ----
library(spdep)
nb_areas <- poly2nb(as(hex_area_win, 'Spatial'), queen=FALSE, row.names = hex_area_win$area_ID) #neighboring areas using sp library 
nbInfo <- nb2WB(nb_areas) #transform into iCAR inputs: adjacent matrix, weights, number of neighbors (for WinBUGS)

adj <- nbInfo$adj
weights <- nbInfo$weights
num <- nbInfo$num
L <- length(nbInfo$adj)

#-------
#Define full adjacent matrix from sparse representation 
n_areas <- length(num) #number of regions. Alternatively, nrow(adj_mat)
#adj_mat <- nb2mat(nb_areas, style="B") #generate adjacency matrix with binary weights, no standardisation (style=B)

#-------------------------------------------------------------------------------
#Constants ----
#Combine constants
constants <- constants_trig

constants$adj <- adj
#constants$adj_mat <- adj_mat
constants$n_areas <- n_areas
constants$weights <- weights
constants$num <- num
constants$L <- L

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
##Create binary environmental variable --
#For example: A forest in central Africa

hex_area_win <- hex_area_win %>% 
  mutate(forest_present = case_when(area_ID %in% c(16,24,19) ~ -400, TRUE ~ 0)) #assume the forest provides a 250 year delay to expansion

##Visulaise the presence/absence of forests
ggplot(data = hex_area_win) +
  geom_sf(data = st_buffer(st_as_sf(sampling_win_outline, crs = 4326), 40000), aes(color = "grey50")) + #sampling window with coastal buffer
  geom_sf(aes(fill = factor(forest_present))) +  #color hex grid by binary variable
  scale_fill_manual(values = c("0" = "grey90", "1" = "green")) + # Define color
  geom_sf(data = as(sites, 'sf'), size=2, alpha=0.5) + #sites
  geom_sf_label(aes(label = area_ID)) + #hex grid labels
  theme(panel.background = element_rect(fill = "lightblue",
                                        colour = "lightblue",
                                        size = 0.5,
                                        linetype = "solid"),
        legend.position = "none")


#-------------------------------------------------------------------------------
##MODEL ---

sim_model <- nimbleCode({
  # Simulate spatially correlated data for all k in 1:n_areas
  for (k in 1:n_areas){
    a[k] ~ dnorm(nu[k], tau.err);
    nu[k] <- x1[k]*beta1[k] + phi[k];
  }
  phi[1:n_areas] ~ dcar_proper(mu = mu[1:n_areas], adj=adj[1:L], num=num[1:n_areas], tau=tau, gamma=gamma) # ICAR prior to capture spatial random effects
  d[1:n_areas] ~ dcar_proper(mu = mu2[1:n_areas], adj=adj[1:L], num=num[1:n_areas], tau=tau, gamma=gamma) 
  b[1:n_areas] <- a[1:n_areas] - abs(d[1:n_areas]) #duration must be positive
  #tau ~ dgamma(2, 0.5)
  
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
sim_constants$x1 <- hex_area_win$forest_present
sim_constants$beta1 <- rep(1, n_areas) #binary: turn on the effect of the forest covariate
sim_constants$mu <- rep(2000, n_areas) #runif(1:sim_constants$n_areas, min = 600, max = 3500) #rep(0, n_areas)
sim_constants$mu2 <- rep(500, n_areas) #runif(1:sim_constants$n_areas, min = 50, max = 600)
sim_constants$tau <- 0.000005
sim_constants$tau.err <- 0.5
sim_constants$gamma <- 0.99

#Define constraints, data, and initial values ----
dat <- list(constraint_uniform = rep(1, sim_constants$n_areas),
            constraint_duration = rep(1, n_sites),
            cra_constraint = rep(1, n_dates))

init_a <- runif(1:sim_constants$n_areas, min = 600, max = 3500)
init_b <- init_a - runif(1:sim_constants$n_areas, min = 50, max = 600)
inits <- list(a = init_a,
              b = init_b)
#              tau = 2)

#Simulate ----
set.seed(1223)
simModel <- nimbleModel(code = sim_model, constants = sim_constants, data = dat, inits = inits)

nodesToSim <- simModel$getDependencies(c("a", "phi", "d", "b", "delta", "alpha", "beta", "theta"), self = T, downstream = T)

simModel$simulate(nodesToSim)
# simModel$a #check variables
# simModel$b
# simModel$theta


##Check spatial autocorrelation with Moran's statistic
#nbw <- nb2listw(nb_areas)
#map <- as(hex_area_win, 'Spatial')
#map$a <- simModel$a
#gmoran <- moran.test(map$a, nbw, alternative = "greater")

##Bounds on gamma
# carMinBound(
#   C= CAR_calcC(adj, num), 
#   adj = adj, 
#   num = num, 
#   M = rep(2000, 41))
# 
# carMaxBound(
#   C= CAR_calcC(adj, num), 
#   adj = adj, 
#   num = num, 
#   M = rep(2000, 41))

#-------------------------------------------------------------------------------
##Visualize arrival times

#Extract simulated arrival times
true_hex_dates <- hex_area_win %>% 
  filter(area_ID %in% 1:47) %>% 
  mutate(true_a = simModel$a)

#Plot
ggplot(data = true_hex_dates) +
  geom_sf(data = st_buffer(st_as_sf(sampling_win, crs = 4326), 40000), aes(color = "grey50")) + #sampling window with coastal buffer
  geom_sf(aes(fill = true_a)) + 
  scale_fill_viridis_c(option="F", direction=-1) +
  scale_alpha_manual(values=c(0.45, 1)) +
  xlab('Longitude') +
  ylab('Latitude') +
  geom_sf_label(aes(label = paste0(round(true_a), "BP")), label.size  = NA, alpha = 0.4, size=3.5) + 
  theme(panel.background = element_rect(fill = "lightblue",
                                        colour = "lightblue",
                                        size = 0.5,
                                        linetype = "solid"),
        legend.position = "none")

#-----------------------------------------------
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
constants$true_a <- simModel$a
constants$true_b <- simModel$b
constants$true_alpha <- simModel$alpha
constants$true_beta <- simModel$beta

#Store output ----
save(sites, siteInfo, sim_df, constants, sampling_win, hex_area_win, file=here('data','tactical_sim_womble.RData'))

