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

source(here('src','gpqrSim.R'))

# Load data (to access constants) ----
load(here('data','eastc14.RData'))

#-------------------------------------------------------------------------------
## List of countries ----
subSahara_countries <- constants$countries #sub-Saharan Africa
eastEIA_countries <- constants$eastEIAcountries #Eastern Sub-Saharan Africa 

#-------------------------------------------------------------------------------
# Generate Spatial Window for Analyses----
sf::sf_use_s2(FALSE) #turn off spherical co-ordinates
sampling_win <-  sampling_win %>%
  st_make_valid() %>%
  st_union()
sf::sf_use_s2(TRUE) #turn on spherical co-ordinates

# Target Parameters ----
true_param  <- list()
true_param$n  <- 400 #number of sites & dates
true_param$origin_point <- st_sfc(st_point(c(-1.45, 31.77))) #dispersal origin point -- approximately at Katuruka
true_param$beta0 <- 2500 #approximate mean date at origin point
true_param$beta1 <- 1 #reciprocal of dispersal rate 
true_param$sigma <- 100 
true_param$etasq <- 0.02 #variability of dispersal rate
true_param$rho <- 350 #range of spatial autocorrelation
true_param$seed <- 1233 #random seed

#Projection ----
st_crs(sampling_win) <- 4326
st_crs(true_param$origin_point) <- 4326

# Simulate Data ----
sim_sites  <- gpqrSim(win = sampling_win,
                      n = true_param$n,
                      beta0 = true_param$beta0,
                      beta1 = true_param$beta1,
                      sigma = true_param$sigma,
                      origin.point = true_param$origin_point,
                      etasq = true_param$etasq,
                      rho = true_param$rho,
                      seed = true_param$seed)

# Save simulation output ----
save(true_param, sim_sites, file=here('data', 'tactical_sim_gpqr.RData'))