# Load Libraries and Data ----
library(rcarbon)
library(nimbleCarbon)
library(sf)
library(rnaturalearth)
library(stringr)
library(dplyr)
library(tidyr)
library(readxl)
library(here)
library(ggplot2)
library(ggthemes)
library(parallel)
library(units)
library(geodata)
library(terra)
library(stars)
library(raster)
library(readxl)

rm(list = ls())

ncores = (detectCores() - 1)

`%!in%` <- Negate(`%in%`)

source(here('src','hex_areas.R'))

#-------------------------------------------------------------------------------
##List of countries in europe (excluding france and norway which need to be added separately) ---
europe_countries <- c("Albania", "Latvia", "Andorra", "Liechtenstein",
                      "Armenia", "Lithuania", "Austria", "Luxembourg",
                      "Azerbaijan", "Malta", "Belarus", "Moldova",
                      "Belgium", "Monaco", "Bosnia and Herzegovina", "Montenegro",
                      "Bulgaria", "Netherlands", "Croatia", 
                      "Cyprus", "Poland", "Czechia", "Portugal",
                      "Denmark", "Romania", "Estonia","Finland", "San Marino",
                      "Macedonia", "Republic of Serbia", "Slovakia",
                      "Georgia", "Slovenia", "Germany", "Spain", "Greece", "Sweden",
                      "Hungary", "Switzerland", "Ireland", "Turkey",
                      "Italy", "Ukraine", "Kosovo", "United Kingdom")

#-------------------------------------------------------------------------------
## Read data sets  ----
betti_dat <- read_excel(here("data", "neolithic_europe_dates_betti.xlsx"))

#-------------------------------------------------------------------------------
## Data filtering and cleaning ----

#Filter dates
betti_sum_df <- betti_dat %>%
  dplyr::select("Lab Number", "Site/Centre name", "Lat", "Longitude", "Uncal C14 BP", "Uncal C14_SD", "Country") %>%
  rename(labCode="Lab Number", siteName="Site/Centre name", lat="Lat", long="Longitude", c14date="Uncal C14 BP", c14std="Uncal C14_SD", country="Country") %>%
  filter(!is.na(lat) & !is.na(long) & !is.na(c14std) & !is.na(c14date) & !is.na(labCode) & labCode!="?") %>% 
  filter(country %in% c(europe_countries, "France", "Norway"))

# Assign ID ----
europe_sites_df <- betti_sum_df %>% mutate(siteID = row_number())

# Assign Site ID ----
#europe_sites_df$siteID  <- as.numeric(factor(europe_sites_df$siteName)) #unnecessary, in this case we only have 1 date per site... 

##-----------
##Save Output as csv
write.csv(europe_sites_df, here('data','europe_sites_df.csv'), row.names = FALSE)

#===============================================================================
## Determining which calibration curve should be used---- (should all be intcal20 since we afre in the northern hemisphere)

#Remove unnecessary information
europe_sites_df <- europe_sites_df %>% dplyr::select(-country)

#Assign calibration curve ----
europe_sites_df$calCurve <- ifelse((europe_sites_df$lat>=0), 'intcal20', 'shcal20') #Assign the calibration curve to use based on the site's position relative to the equator #TODO: refine this -- weird to have a hard step-change between calibration curves at the equator -- maybe use a gradient change function? Mixed Curve? Or is there perhaps better regional calibration curves to use?

#-------------------------------------------------------------------------------
## Restructure Data for Bayesian Analyses ----

# Compute median calibrated dates ----
europe_sites_df$median_dates = medCal(calibrate(europe_sites_df$c14date,
                                                europe_sites_df$c14std,
                                                calCurve = europe_sites_df$calCurve,
                                                ncores = ncores))

#-------------------------------------------------------------------------------
## Designating approximate origin ----

# Possible start-point (oldest date) in easter_EIA dataset
possible_origin_dat <- europe_sites_df %>%
  slice_max(c14date, n=1)

## Compute Great-Arc Distances in km ----
sites <- st_as_sf(europe_sites_df, coords = c('long','lat'))
st_crs(sites)  <- 4326 
dist_mat  <- set_units(st_distance(sites), 'km') #inter-site distance matrix in km: each site's distance from every other site (i.e. with n sites, this matrix is n^2)
origin_point  <- sites %>% filter(siteName == possible_origin_dat$siteName)
dist_org  <-  as.vector(set_units(st_distance(x=sites, y=origin_point), 'km')) #distance from origin site

#-------------------------------------------------------------------------------
##Generate Spatial Window for Analyses
#Sampling window: Europe ----

sampling_win <- st_union(st_make_valid(ne_countries(country=europe_countries, returnclass = "sf")),
                         st_make_valid(ne_countries(geounit = c("france", "norway"), type = "map_units", returnclass = "sf"))) # France filter map_units by geounit to exclude French Guiana. Similar problem for Norway

##Plot sample window
# ggplot(data = hex_area_win) +
#   geom_sf(data = st_buffer(st_as_sf(sampling_win, crs = 4326), 40000), aes(color = "grey50"), lwd=1) + #internal country borders
#   geom_sf_label(data = sampling_win, aes(label = admin, alpha=0.6), color="darkred", size=4) + #country labels
#   theme(panel.background = element_rect(fill = "lightblue",
#                                         colour = "lightblue",
#                                         size = 0.5,
#                                         linetype = "solid"),
#         legend.position = "none")

#Generate Spatial Hexagons ---- ##see code block below to determine hex diameter, cell_d 
hex_area_win <- hex_areas(sampling_win, cell_d = 4.4) 

#Remove spatial hexagons where the neolithic farmers didn't reach
hex_area_win <- hex_area_win %>%
  filter(area_ID %!in% c(35, 44, 54, 64, 75, 86, 81, 92, 98, 87, 76, 65, 70, 59, 40, 28, 49)) %>%
  mutate(area_ID = row_number())

##CHECK -- plot hexs and sites
# ggplot(data = hex_area_win) +
# geom_sf(data = st_buffer(st_as_sf(sampling_win, crs = 4326), 40000), aes(color = "grey50")) + #sampling window with coastal buffer
# geom_sf() + #hex grid
# geom_sf(data = as(sites, 'sf'), size=2, alpha=0.5) + #sites
# geom_sf_label(aes(label = area_ID)) + #hex grid labels
# geom_sf(data = hex_area_win$area_center, size=2, alpha=1, aes(color = "purple")) + #hex-origins
# theme(panel.background = element_rect(fill = "lightblue",
#                                       colour = "lightblue",
#                                       size = 0.5,
#                                       linetype = "solid"),
#       legend.position = "none")


#Assign hex area id to each site ----
europe_sites_df$area_id <- as.integer(st_within(sites$geometry, hex_area_win$geometry)) 
#Remove NAs
europe_sites_df <- europe_sites_df %>% filter(!is.na(area_id))
sites <- sites %>% filter(siteID %in% europe_sites_df$siteID)
  
##CHECK ---
#area_freq  <- plyr::count(europe_sites_df, 'area_id') ##See how many sites fall in each hex area. 

#--------------------------------
## Determining hex size ---
#Under changing hex size, determine the proportion of areal hex units in the sampling window with sites
# prop_units_df <- data.frame(d = numeric(), prop_with_sites = numeric())
# 
# for (d in seq(1, 15, 0.1)){
#   hex_area_win <- hex_areas(sampling_win, cell_d = d)
#   europe_sites_df$area_id <- as.integer(st_within(sites$geometry, hex_area_win$geometry))
# 
#   hex_with_sites <- length(unique(europe_sites_df$area_id))
#   all_hex <- length(hex_area_win$area_ID)
# 
#   prop_with_sites <- hex_with_sites/all_hex
# 
#   prop_units_df <- rbind(prop_units_df, data.frame(d = d, prop_with_sites = prop_with_sites))
# }

## Plot results
# ggplot(prop_units_df, aes(x = d, y = prop_with_sites)) +
#   geom_line() +
#   geom_point() +
#   scale_x_continuous(breaks=seq(0,15,by=1))+
#   labs(x = "Hexagon Size (d)", y = "Proportion of Hexagons with Sites", title = "Effect of Hexagon Size on Site Coverage") +
#   theme_minimal()

#-------------------------------------------------------------------------------
## Create list with constants and data ----

# Data
europe_dat <- list(cra=europe_sites_df$c14date,
                   cra_error=europe_sites_df$c14std) #Creating europe_dat with only date information

## Constants
data(intcal20)
constants <- list()
constants$countries <- c(europe_countries, "France", "Norway") 
constants$n_sites <- nrow(europe_sites_df) #n_sites = n_dates
constants$n_areas  <- nrow(hex_area_win) #All areas (even empty ones) are included #Only occupied areas: length(unique(siteInfo$area_id))
constants$id_sites <- europe_sites_df$siteID
constants$id_areas  <- europe_sites_df$area_id 
constants$dist_mat  <- dist_mat
constants$dist_org  <- dist_org
constants$origin_point <- st_coordinates(origin_point)

#Calibration curves
constants$calBP <- intcal20$CalBP 
constants$C14BP  <- intcal20$C14Age #Northern hemisphere calibration curves
constants$C14err  <- intcal20$C14Age.sigma

#-------------------------------------------------------------------------------
## Save everything on a R image file ----
save(sites, constants, europe_dat, europe_sites_df, sampling_win, hex_area_win, file=here('data','europe.RData'))