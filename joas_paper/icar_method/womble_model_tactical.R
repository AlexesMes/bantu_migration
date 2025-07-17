#Load Libraries ----
library(here)
library(coda)
library(nimbleCarbon)
library(rcarbon)
library(dplyr)
library(msm)

rm(list = ls())
`%!in%` <- Negate(`%in%`)

set.seed(123)

##ICAR model tactical simulation without errors (i.e. using dates that have no associated error and calibration)

#-------------------------------------------------------------------------------
## Data Setup ----
load(here('data', 'tactical_sim_womble.RData'))
load(here('data','trig.RData')) #nodes and edges between hex area centroids

#Combine constants
constants <- c(constants, constants_trig)
#-------------------------------------------------------------------------------
# ##Simulated data summary
# #Number of dates in area
#dates_in_areas_summarise <- as.data.frame(table(sites$area_id))
# #Number of sites in area
#sites_in_areas_summarise <- sites %>% group_by(area_id) %>% summarize(n_sites =n_distinct(site_id))

#-------------------------------------------------------------------------------





