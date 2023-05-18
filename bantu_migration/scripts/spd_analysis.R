# Load Libraries and Data ----

library(rcarbon)
library(nimbleCarbon)
library(maptools)
library(sf)
library(stringr)
library(dplyr)
library(parallel)


ncores = (detectCores() - 1)

#-------------------------------------------------------------------------------
## NOTES ----
#Exploratory Analysis is marked ExA

#-------------------------------------------------------------------------------
## Read input files -----

source("scripts/prepare_data.R")


#-------------------------------------------------------------------------------
## Set constants -----

min_14Cage <- 0
n_bins <- 100 #time window for binning dates (in years) within the same site
n_simulations <- 100 #num of simulations
timeRange <- c(4000, 0) #set the timerange of analysis in calBP, older date first
breaks <- seq(4000, 0, -200) #200 year blocks

cal_norm <- FALSE #TODO: Explore how results differ when normalised vs unnormalised dates are used
calCurve <- 'intcal20'



#-------------------------------------------------------------------------------
## Calibrate, bin, spd ----


## Calibrate ----
bantu_caldates <- rcarbon::calibrate(x = bantu_sites_df$c14date, 			
                                     errors = bantu_sites_df$c14std,
                                     calCurves = calCurve, 
                                     normalised = cal_norm,
                                     ncores = ncores)

#Basic plot
plot(bantu_caldates)
#summary(bantu_sites_df.caldates) #QQ: Why doesn't this work, but calling summary on cal_dates_restricted (below) works?

#ExA: Extract specific dates (all dates with a CDF of 0.5 or over between 2000 and 500 cal BP)
which.CalDates(bantu_caldates, BP<=2000 & BP>=500, p=0.5)
cal_dates_restricted <- subset(bantu_caldates, BP<=2000 & BP>=500, p=0.5)
summary(cal_dates_restricted)


## Binning ----
#The archaeological record in West Africa seems particularly biased with several sites having a potentially outsized number of dates -- most likely along the river networks
bins <- binPrep(sites = bantu_sites_df$site,
                ages = bantu_sites_df$c14date,
                h = n_bins)

#QQ: Should I rather use the median calibrated dates to bin?
#bins_cal <- binPrep(sites = bantu_sites_df$site, ages = bantu_caldates, h = n_bins) ##QQ: Why does this line of code not work?



print(paste("Analysis based on", nrow(bantu_sites_df), "dates and", length(unique(bins)), "bins"))



## Summed Probability Distribution ----

bantu_spd =spd(bantu_caldates, bins=bins, timeRange = timeRange)
plot(bantu_spd)
#plot(bantu_spd,runm=200,add=TRUE,type="simple",col="darkorange",lwd=1.5,lty=2) #using a rolling average of 200 years for smoothing
#plot(bantu_spd, calendar='BCAD') #Using a BC/AD timescale


#Sensitivity Analysis -- visual assessment of how different cut-off values modify the shape of the spd
binsense(x=bantu_caldates, y=bantu_sites_df$site, h=seq(0,500,100), timeRange = timeRange)

#Visualising bins
binsmed = binMed(x = bantu_caldates, bins=bins)
plot(bantu_spd,runm=200)
barCodes(binsmed, yrng = c(0,0.01))





