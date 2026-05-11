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
library(readr)

rm(list = ls())

ncores = (detectCores() - 1)

`%!in%` <- Negate(`%in%`)

source(here('src','hex_areas.R'))

#-------------------------------------------------------------------------------
## List of countries in eastern and southern sub-Saharan Africa ----
eastEIA_countries <- c("South Africa", 
                       "Lesotho", 
                       "Eswatini", "eSwatini", "Swaziland", 
                       "Botswana",
                       "Zimbabwe",
                       "Zambia",
                       "Mozambique",
                       "Malawi",
                       "Tanzania", "United Republic of Tanzania",
                       "Rwanda",
                       "Burundi",
                       "Kenya",
                       "Uganda") 
                       #"Madagascar",
                       #"Comoros")

#-------------------------------------------------------------------------------
## Data sources ----

##SARD
#SARD_dat contains data from the Southern African Radiocarbon Database
#Downloaded from: https://github.com/emmaloftus/Southern-African-Radiocarbon-Database
#See Loftus et al., 2019, "An archaeological radiocarbon database for southern Africa" Antiquity. (https://doi.org/10.15184/aqy.2019.75).

##Collected dates 
#Compiled dates covering Eastern and Southern Africa

##Wanyika database
#Eastern africa radiocarbon dates from Kenya, Tanzania, the Comoros Islands, and Madagascar
#https://pandoradata.earth/dataset/wanyika

#-------------------------------------------------------------------------------
## Read data sets  ----
SARD_dat <- read.csv(here("data", "SARD_database.csv"))
Collected_EIA_dat <- read_excel(here("data", "collected_dates.xlsx"), sheet = "Collected_dates_final") 
Wanyika_dat <- read_csv("data/wanyika_database.csv")

#-------------------------------------------------------------------------------
## Data filtering and cleaning ----

SARD_dat <- SARD_dat %>% 
  mutate(Date = case_when(Lab.ID =="Pta-2360" ~ '1760', .default = Date),
         Uncertainty = case_when(Lab.ID =="Pta-2360" ~ '40', .default = Uncertainty)) %>% #see "Excavations at Silver Leaves: A Final Report", M. Klapwijk and T. N. Huffman
  mutate(Archaeological.Period = case_when(Lab.ID %in% c("Pta-1818", "Pta-1959", "Pta-1961") ~ "Iron Age", .default = Archaeological.Period)) #see "SOME RECENT RADIOCARBON DATES FROM SOUTHERN AFRICA", M. HALL AND J. C. VOGEL, 1980

#Filter dates
SARD_df <- SARD_dat %>%
  filter(Archaeological.Period=="Iron Age") %>% 
  dplyr::select(Lab.ID, X.Site, DecdegE, DecdegS, Date, Uncertainty, Material.dated, Country) %>%
  rename(labCode=Lab.ID, siteName=X.Site, lat=DecdegS, long=DecdegE, c14date=Date, c14std=Uncertainty, material=Material.dated, country=Country) %>%
  mutate(c14date = as.numeric(c14date), c14std=as.numeric(c14std), dataorigin="SARD")  %>%
  filter(!is.na(lat) & !is.na(long) & !is.na(c14std) & !is.na(c14date)) %>% 
  filter(siteName !="Bambata Cave") %>% #Bambata has a compicated hisotry with being included or excluded as EFC ware
  filter(labCode != "Beta-11112") #Marine shell -- several centuries before EIA communities arrived in region

collected_df <- Collected_EIA_dat %>%
  dplyr::select(LabID, Age, Error, Long, Lat, SiteName,  Material, Country) %>%
  rename(labCode=LabID, siteName=SiteName, lat=Lat, long=Long, c14date=Age, c14std=Error, material=Material, country=Country) %>%
  mutate(c14date = as.numeric(c14date), c14std=as.numeric(c14std), dataorigin="Collected")  %>%
  filter(!is.na(lat) & !is.na(long) & !is.na(c14std) & !is.na(c14date)) 

wanyika_df <- Wanyika_dat %>%
  dplyr::select("Labcode", "Date BP", "Date BP SD", "Longitude", "Latitude", "Site Name", "Dated Material", "Country", "Regional Cultural Phase (Eastern Africa)", "MEAN GRADE: Chrono Hygiene 1 and 2 (Stratigraphic integrity  +  SD scores)") %>%
  rename(labCode="Labcode", siteName="Site Name", lat="Latitude", long="Longitude", c14date="Date BP", c14std="Date BP SD", material="Dated Material", country="Country", grade="MEAN GRADE: Chrono Hygiene 1 and 2 (Stratigraphic integrity  +  SD scores)", phase="Regional Cultural Phase (Eastern Africa)") %>%
  mutate(c14date = as.numeric(c14date), c14std=as.numeric(c14std), dataorigin="Wanyika")  %>%
  filter(phase %in% c("EIA","MIA", "LIA","EIA/MIA","MIA/LIA","EIA/MIA/LIA")) %>% #grade %in% c("A","B") #(UNCOMMENT if generating 'eastc14_wanAB.RData')
  filter(!is.na(lat) & !is.na(long) & !is.na(c14std) & !is.na(c14date)) %>% 
  filter(labCode %!in% c("Pta-8527", "Pta-8522", "Ua-38476","OxA-18870")) %>% #These dates from "Ukunju Cave","Makangale","Wei wei Valley" are questioned
  dplyr::select(-grade, -phase)
  

##Filter out dates in collected_df which already exist in the wanyika and SARD databases ----
overlap_labID_wc <- merge(wanyika_df, collected_df, by="labCode")$labCode
collected_df <- collected_df %>% filter(labCode %!in% overlap_labID_wc)

overlap_labID_sc <- merge(SARD_df, collected_df, by="labCode")$labCode
collected_df <- collected_df %>% filter(labCode %!in% overlap_labID_sc)

##Dates to check
#collected_df_prefilter <- collected_df %>% filter(labCode %!in% overlap_labID_wc)
#dates_to_check <- anti_join(collected_df, collected_df_prefilter, by="labCode")
#write.csv(dates_to_check, here('dates_to_check.csv'), row.names = FALSE)

##------------
## Combine datasets ----
bantu_sites_df <- bind_rows(SARD_df, wanyika_df, collected_df) %>% 
  mutate(dataorigin=as.factor(dataorigin)) %>% 
  filter((c14date != 0) & (c14std != 0)) %>% #Some of these datasets had modern dates (indicated with c14 date and error of 0), we remove these
  filter((c14date <=7000) & (c14date >=246)) #We assume 1652 AD as the cut-off.
  
##-----------
##Filter for eastern EIA stream
eastEIA_sites_df <- bantu_sites_df %>% 
  filter(country %in% eastEIA_countries) %>% 
  mutate(ID = row_number()) %>%  #Assign ID
  mutate(siteID = as.numeric(factor(siteName))) #Assign Site ID 

##-----------
##Save Output as csv
write.csv(eastEIA_sites_df, here('data','eastEIA_dataset.csv'), row.names = FALSE)

#===============================================================================
## Determining which calibration curve should be used----

#Remove unnecessary information
eastEIA_sites_df <- eastEIA_sites_df %>% 
  dplyr::select(-material, -country)

#Assign calibration curve ----
eastEIA_sites_df$calCurve <- ifelse((eastEIA_sites_df$lat>=0), 'intcal20', 'shcal20') #Assign the calibration curve to use based on the site's position relative to the equator #TODO: refine this -- weird to have a hard step-change between calibration curves at the equator -- maybe use a gradient change function? Mixed Curve? Or is there perhaps better regional calibration curves to use?


#-------------------------------------------------------------------------------
## Restructure Data for Bayesian Analyses ----

# Compute median calibrated dates ----
eastEIA_sites_df$median_dates = medCal(calibrate(eastEIA_sites_df$c14date,
                                               eastEIA_sites_df$c14std,
                                               calCurve = eastEIA_sites_df$calCurve,
                                               ncores = ncores))


# Collect site level information ----
earliest_dates <- aggregate(median_dates ~ siteID, data=eastEIA_sites_df, FUN = max) #Earliest medCal Date for Each Site 
latest_dates <- aggregate(median_dates ~ siteID, data=eastEIA_sites_df, FUN=min) #Latest medCal Date for Each Site
n_dates <- aggregate(median_dates ~ siteID, data=eastEIA_sites_df, FUN=length) #Number of medCal Date for Each Site

siteInfo <- data.frame(siteID = earliest_dates$siteID,
                       earliest = earliest_dates$median_dates,
                       latest = latest_dates$median_dates,
                       diff = earliest_dates$median_dates - latest_dates$median_dates,
                       n_dates = n_dates$median_dates) %>% unique()

siteInfo <- siteInfo %>% 
  left_join(unique(dplyr::select(eastEIA_sites_df,
                                 lat,
                                 long,
                                 siteID,
                                 siteName,
                                 dataorigin,
                                 calCurve))) %>% 
  distinct(siteID, .keep_all = TRUE)  #Some site duplicates introduced by slight differences in co-ordinate decimal places

# Collect date level information ----
dateInfo <- unique(dplyr::select(eastEIA_sites_df,
                                  ID,
                                  labCode,
                                  siteID,
                                  siteName,
                                  cra=c14date,
                                  cra_error=c14std,
                                  median_dates=median_dates,
                                  calCurve=calCurve)) %>% arrange(ID) 
dateInfo$earliestAtSite  <- FALSE #initialize 

for (i in unique(siteInfo$siteID))
{
  tmp_index  <-  which(dateInfo$siteID==i) #Identify the row indexes of all sites that share siteID x
  ii  <- tmp_index[which.max(dateInfo$median_dates[tmp_index])] #Select the index of the row with the maximum median date (earliest date) within the sub-group of rows selected in the previous line (given by tmp.index)
  dateInfo$earliestAtSite[ii]  <- TRUE #Designate this observation as the earliest date
} 
#table(dateInfo$earliestAtSite) #Check

#-------------------------------------------------------------------------------
## Designating approximate origin ---- (Note: not used ICAR model)

# Possible start-point (oldest date) in eastEIA dataset
possible_origin_dat <- eastEIA_sites_df %>%
  slice_max(c14date, n=1)

## Compute Great-Arc Distances in km ----
sites <- st_as_sf(siteInfo, coords = c('long','lat'))
st_crs(sites)  <- 4326 
dist_mat  <- set_units(st_distance(sites), 'km') #inter-site distance matrix in km: each site's distance from every other site (i.e. with n sites, this matrix is n^2)
origin_point  <- sites %>% filter(siteName == possible_origin_dat$siteName)
dist_org  <-  as.vector(set_units(st_distance(x=sites, y=origin_point), 'km')) #distance from origin site

#-------------------------------------------------------------------------------
# Generate Spatial Window for Analyses: Sub-Saharan Africa ----

#Sampling window: Eastern Sub-Saharan Africa ----
sampling_win <- ne_countries(continent = "Africa", country = eastEIA_countries, returnclass = "sf", scale="large") #the detailed resolution (scale='large') ensures the Comoros islands are included

#Generate Spatial Hexagons ---- ##see code block below to determine hex diameter, cell_d 
hex_area_win <- hex_areas(sampling_win, cell_d = 3.8)

#Remove spatial hexagons which are sample window edges or where EIA didn't reach ---
##East Africa, cell-diameter 3.8
hex_area_win <- hex_area_win %>%
  filter(area_ID %!in% c(1,2,3,4,5,8,9,12,13,25,38,44,50,55,59,60,64)) %>%
  mutate(area_ID = row_number())

##East Africa, cell-diameter 2.9 (i.e. finer spatial resolution -- UNCOMMENT if generating 'eastc14_d29.RData')
# hex_area_win <- hex_area_win %>%
#   filter(area_ID %!in% c(1,2,3,4,5,6,8,9,10,13,14,18,19,24,77,78,90,100,102)) %>%
#   mutate(area_ID = row_number())

#---------
##CHECK -- plot hexs and sites
# #interest_sites <- eastEIA_sites_df %>% filter(siteName %in% c("University Campus","Kamukombe","Zitundo","Caimane"))
# #interest_sites$points <- st_geometry(st_as_sf(interest_sites, coords = c("long","lat"), crs=4326))
# 
# ggplot(data = hex_area_win) +
#   geom_sf(data = st_buffer(st_as_sf(sampling_win, crs = 4326), 2000), aes(color = "grey50"),lwd=1.5) + #sampling window with coastal buffer
#   geom_sf(aes(alpha=0.2)) + #hex grid
#   geom_sf(data = as(sites, 'sf'), size=2, alpha=0.5) + #sites
#   geom_sf_label(aes(label = area_ID),size=3) + #hex grid labels
#   #geom_sf(data = interest_sites$points, size=2, alpha=1, color = "purple") +
#   #geom_sf_label(data = interest_sites$points, aes(label = interest_sites$siteName),size=3) + #hex grid labels
#   #geom_sf_label(data = sampling_win, aes(label = admin, alpha=0.6), color="darkred", size=4) + #country labels
#   theme(panel.background = element_rect(fill = "lightblue",
#                                         colour = "lightblue",
#                                         size = 0.5,
#                                         linetype = "solid"),
#         legend.position = "none")

#Assign hex area id to each site ----
siteInfo$area_id <- as.integer(st_within(sites$geometry, hex_area_win$geometry))

#Assign hex area id to each date ----
dateInfo$area_id <- siteInfo$area_id[match(dateInfo$siteID, siteInfo$siteID)]

##CHECK ---
#area_freq  <- plyr::count(siteInfo, 'area_id') ##See how many sites fall in each hex area. Also make sure there are no 'NA' entries

##-------------------------------------------------------------------------------
## Create list with constants and data ----

# Data
eastEIA_dat <- list(cra=dateInfo$cra,
                    cra_error=dateInfo$cra_error) #Creating eastEIA_dat with only date information

## Constants
data(intcal20)
data(shcal20)
constants <- list()
constants$eastEIAcountries <- eastEIA_countries
constants$n_sites <- nrow(siteInfo)
constants$n_dates  <- nrow(dateInfo)
constants$n_areas  <- nrow(hex_area_win) #All areas (even empty ones) are included #Only occupied areas: length(unique(siteInfo$area_id))
constants$id_sites <- dateInfo$siteID
constants$id_areas  <- siteInfo$area_id 
constants$dist_mat  <- dist_mat
constants$dist_org  <- dist_org
constants$origin_point <- st_coordinates(origin_point)
#Calibration curves
constants$calBP <- intcal20$CalBP #Same for intcal20 and shcal20 
constants$C14BP  <- cbind(intcal20$C14Age, shcal20$C14Age) #Northern and southern hemisphere calibration curves
constants$C14err  <- cbind(intcal20$C14Age.sigma, shcal20$C14Age.sigma)

#-------------------------------------------------------------------------------
## Save everything on a R image file ----
save(sites, constants, eastEIA_dat, siteInfo, dateInfo, sampling_win, hex_area_win, file=here('data','eastc14.RData')) #'eastc14_d29.RData' #'eastc14_wanAB.RData'

#-------------------------------------------------------------------------------
## Save sampling window specific information separately ----
 
constants_sw <- list()
constants_sw$eastEIAcountries <- eastEIA_countries
constants_sw$n_areas  <- constants$n_areas
constants_sw$calBP <- constants$calBP
constants_sw$C14BP  <- constants$C14BP 
constants_sw$C14err  <- constants$C14err

save(constants_sw, sampling_win, hex_area_win, file=here('data','sample_window.RData')) #'sample_window_d29.RData'