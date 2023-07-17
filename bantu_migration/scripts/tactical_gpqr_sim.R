# Load Libraries ----
library(here)
library(dplyr)
library(stringr)
library(nimbleCarbon)
library(rnaturalearth)
library(sp)
library(maptools)
library(sf)
library(ggplot2)
library(viridis)
library(rgeos)
library(rcarbon)

source(here('src','gpqrSim.R'))

# Generate Spatial Window for Analyses: Sub-Saharan Africa ----
sf_subsah_africa <- ne_countries(continent = "Africa", returnclass = "sf") %>%
  filter_all(., any_vars(str_detect(., "Sub-Saharan"))) %>% 
  filter_all(name_en != "Madagascar") #We focus on mainland sub-Saharan Africa

sp_ss_africa <- sf_subsah_africa %>% as("Spatial") #convert sf to sp object

sampling_win <- as(sp_ss_africa, "SpatialPolygons") |>  unionSpatialPolygons(IDs = rep(1, nrow(sp_ss_africa)))
sampling_win <- disaggregate(sampling_win) #create new raster layer with higher resolution (smaller cells)
sampling_win  <- sampling_win[order(raster::area(sampling_win), decreasing=TRUE)]



# Target Parameters ----
true_param  <- list()
true_param$n  <- 1600 #number of sites & dates
true_param$origin_point <- c(11.50, 3.82) #dispersal origin point -- approximately at Obobogo
true_param$beta0 <- 3070 #mean date at origin point
true_param$beta1 <- 0.15 #reciprocal of dispersal rate ##TODO: This should be increase to around 0.3, but other parameters need changing such that theta (out.df$theta in gpqrSim()) remains positive
true_param$sigma <- 100 
true_param$etasq <- 0.06 #variability of dispersal rate
true_param$rho <- 150 #range of spatial autocorrelation
true_param$seed <- 1233 #random seed

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