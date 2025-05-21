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
## List of countries in sub-Saharan Africa ----
subSahara_countries <- c("South Africa", 
                         "Lesotho", 
                         "Eswatini", "eSwatini", "Swaziland", 
                         "Botswana",
                         "Zimbabwe",
                         "Namibia",
                         "Angola",
                         "Zambia",
                         "Mozambique",
                         "Malawi",
                         "Madagascar",
                         "Comoros",
                         "Tanzania", "United Republic of Tanzania",
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

## List of countries associated with eastern EIA stream ----
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
                       "Uganda", 
                       "Madagascar",
                       "Comoros")

#-------------------------------------------------------------------------------
## Data sources ----

##SARD
#SARD_dat contains data from the Southern African Radiocarbon Database
#Downloaded from: https://github.com/emmaloftus/Southern-African-Radiocarbon-Database
#See Loftus et al., 2019, "An archaeological radiocarbon database for southern Africa" Antiquity. (https://doi.org/10.15184/aqy.2019.75).

##aDRAC
#HumActCA_dat contains data from the Human Activity in Central Africa (Archive des datations radiocarbones d’Afrique centrale) Database 
#Downloaded from: https://zenodo.org/record/4394894 or https://github.com/dirkseidensticker/aDRAC [Latter appears to be more up to date]
#See Seidensticker et al., 2021, "Population collapse in Congo rainforest from 400 CE urges reassessment of the Bantu Expansion" (https://www.science.org/doi/full/10.1126/sciadv.abd8352)

##Collected dates 
#Compiled dates covering Eastern and Southern Africa

##Wanyika database
#Eastern africa radiocarbon dates from Kenya, Tanzania, the Comoros Islands, and Madagascar
#https://pandoradata.earth/dataset/wanyika

#-------------------------------------------------------------------------------
## Read data sets  ----
SARD_dat <- read.csv(here("data", "SARD_Mar2021_14C.csv"))
HumActCA_dat <- read.csv(here("data", "HumActCA_Dec2020_14C.csv"))
aDRAC_dat <- read.csv(here("data", "aDRAC_Feb2024_14C.csv"))#This is a more up-to-date version of HumActCA_dat, but not all the dates in aDRAC_dat are iron age, whereas all dates in the HumActCA are
Collected_EIA_dat <- read_excel(here("data", "collected_dates.xlsx"), sheet = "Collected_dates_final") 
Wanyika_dat <- read_csv("data/wanyika_chronological_database.csv")

#-------------------------------------------------------------------------------
## Data filtering and cleaning ----

#Cleaning up aDRAC database
aDRAC_dat <- aDRAC_dat %>%
  mutate(IRON = case_when(IRON %in% c("", "-") ~ NA, .default = IRON)) %>% 
  mutate(COUNTRY = case_when(COUNTRY == "AGO" ~ "Angola",
                             COUNTRY == "BDI" ~ "Burundi",
                             COUNTRY == "CAF" ~ "Central African Republic",
                             COUNTRY == "CMR" ~ "Cameroon",
                             COUNTRY == "COD" ~ "Democratic Republic of the Congo",
                             COUNTRY == "COG" ~ "Republic of the Congo",
                             COUNTRY == "GAB" ~ "Gabon",
                             COUNTRY == "GEQ" ~ "Equatorial Guinea",
                             COUNTRY == "GNQ" ~ "Equatorial Guinea",
                             COUNTRY == "RWA" ~ "Rwanda",
                             COUNTRY == "TCD" ~ "Chad"),
                             .default = COUNTRY)


#Correct incorrect data entries
aDRAC_dat <- aDRAC_dat %>%
  mutate(SITE = case_when(SITE == "Ngoma III" ~ "Ngoma", .default = SITE)) %>% 
  mutate(PHASE = case_when(LABNR %in% c("OxTL-209a", "OxTL-209c", "OxTL-209d") ~ "Iron Age", .default = PHASE))
  
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
  filter(siteName !="Bambata Cave") %>% #Bambata designated pre-bantu (references given in Isern and Fort 2019, Suplementary Material S1, https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0215573) 
  filter(labCode != "Beta-11112") #Marine shell -- several centuries before EIA communities arrived in region

collected_df <- Collected_EIA_dat %>%
  dplyr::select(LabID, Age, Error, Long, Lat, SiteName,  Material, Country) %>%
  rename(labCode=LabID, siteName=SiteName, lat=Lat, long=Long, c14date=Age, c14std=Error, material=Material, country=Country) %>%
  mutate(c14date = as.numeric(c14date), c14std=as.numeric(c14std), dataorigin="Collected")  %>%
  filter(!is.na(lat) & !is.na(long) & !is.na(c14std) & !is.na(c14date)) 

aDRAC_sum_df <- aDRAC_dat %>% #To find extra dates not in HumActCA
  dplyr::select(LABNR, SITE, LAT, LONG, C14AGE, C14STD, PHASE, CLASS, IRON, MATERIAL, COUNTRY) %>%
  rename(labCode= LABNR, siteName=SITE, lat=LAT, long=LONG, c14date=C14AGE, c14std=C14STD, material=MATERIAL, country=COUNTRY) %>%
  filter(!is.na(lat) & !is.na(long) & !is.na(c14std) & !is.na(c14date)) %>%
  mutate(c14date = as.numeric(c14date), c14std=as.numeric(c14std), dataorigin="aDRAC") %>%
  filter(siteName %in% HumActCA_dat$SITE | PHASE %in% c("N; EIA", "LIA", "EIA", "Iron Age") | !is.na(IRON)) %>% #All dates in HumActCA_dat should be associated with Bantu pottery finds
  filter(CLASS %!in% c('IIa', 'IIb', 'IIc', 'IIIa', 'IIIb', 'IIIc')) %>% #Filter for unreliable class of dates III and irrelevant (according to Seidensticker et al.) dates II
  dplyr::select(-PHASE, -CLASS, -IRON)

wanyika_df <- Wanyika_dat %>%
  dplyr::select("Labcode", "Date BP", "Date BP SD", "Longitude", "Latitude", "Site Name", "Dated Material", "Country") %>%
  rename(labCode="Labcode", siteName="Site Name", lat="Latitude", long="Longitude", c14date="Date BP", c14std="Date BP SD", material="Dated Material", country="Country") %>%
  mutate(c14date = as.numeric(c14date), c14std=as.numeric(c14std), dataorigin="Wanyika")  %>%
  filter(!is.na(lat) & !is.na(long) & !is.na(c14std) & !is.na(c14date)) 

##Filter out dates in collected_df which already exist in the wanyika and SARD databases ----
overlap_labID_wc <- merge(wanyika_df, collected_df, by="labCode")$labCode
collected_df <- collected_df %>% filter(labCode %!in% overlap_labID_wc)

overlap_labID_sc <- merge(SARD_df, collected_df, by="labCode")$labCode
collected_df <- collected_df %>% filter(labCode %!in% overlap_labID_sc)

##------------
## Combine datasets ----
bantu_sites_df <- bind_rows(SARD_df, wanyika_df, aDRAC_sum_df, collected_df) %>% 
  mutate(dataorigin=as.factor(dataorigin)) %>% 
  filter((c14date != 0) & (c14std != 0)) %>% #Apparently some of these datasets had modern dates (indicated with c14 date and error of 0), we remove these
  filter((c14date <=7000) & (c14date >=246)) #Dates earlier than this are assumed to not be of Bantu origin #TODO: build some flexibility into this... #Further, we assume the Dutch arrival in the Cape (1652) as the cut-off. Dates after this point of colonial contact are not considered.
  
# Assign ID ----
bantu_sites_df <- bantu_sites_df %>% mutate(ID = row_number())

# Assign Site ID ----
bantu_sites_df$siteID  <- as.numeric(factor(bantu_sites_df$siteName))

##-----------
##Filter for eastern EIA stream
eastEIA_sites_df <- bantu_sites_df %>% 
  filter(country %in% eastEIA_countries) %>% 
  mutate(ID = row_number()) %>%  #reassign ID for new dataset
  mutate(siteID = as.numeric(factor(siteName))) #reassign site ID for new dataset

##-----------
##Save Output as csv
write.csv(bantu_sites_df, here('data','bantu_dataset.csv'), row.names = FALSE)
write.csv(eastEIA_sites_df, here('data','eastEIA_dataset.csv'), row.names = FALSE)

#===============================================================================
##FROM HERE ONWARDS WE USE THE EASTERN EIA STREAM
#But can easily add back in western EIA and rainforest pottery-using (aDRAC) 
#dates by using bantu_sites_df instead of eastEIA_sites_df below
#===============================================================================
## Determining which calibration curve should be used----

#Remove unnecessary information
eastEIA_sites_df <- eastEIA_sites_df %>% #bantu_sites_df %>% 
  dplyr::select(-material, -country)

#Assign calibration curve ----
eastEIA_sites_df$calCurve <- ifelse((eastEIA_sites_df$lat>=0), 'intcal20', 'shcal20') #Assign the calibration curve to use based on the site's position relative to the equator #TODO: refine this -- weird to have a hard step-change between calibration curves at the equator -- maybe use a gradient change function? Mixed Curve? Or is there perhaps better regional calibration curves to use?
#table(eastEIA_sites_df$calCurve)

#-------------------------------------------------------------------------------
## Restructure Data for Bayesian Analyses ----

# Compute median calibrated dates ----
eastEIA_sites_df$median_dates = medCal(calibrate(eastEIA_sites_df$c14date,
                                               eastEIA_sites_df$c14std,
                                               calCurve = eastEIA_sites_df$calCurve,
                                               ncores = ncores))


# Collect site level information ----
earliest_dates <- aggregate(median_dates ~ siteID, data = eastEIA_sites_df, FUN = max) #Earliest medCal Date for Each Site 
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
## Designating approximate origin ----

# #A generally acknowledged 'homeland region on the border between Nigeria and Cameroon' (Seidensticket,et al. 2020)
# possible_origin_dat <- HumActCA_dat %>%
#   filter(COUNTRY %in% c('CMR', 'NGA')) %>% #I don't think there are any sites in Nigeria in this dataset, but still ...
#   arrange(desc(C14AGE))
# 
# #Note: oldest class I date is Shum Laka -- but see recent genetic analysis discussing why this is a hunter-gather site ('Ancient West African foragers in the context of African population history', Lipson, 2020, https://www.nature.com/articles/s41586-020-1929-1). Also, see Piere de Maret's thoughts on the site. #TODO: Decide whether to designate pre-Bantu and exclude Shum Laka from analyses entirely...
# #Therefore, as an approximate origin we select the next oldest class I date in the region: Ngoume
# possible_origin_dat <- possible_origin_dat %>%
#   filter(SITE=='Ngoume' & (CLASS %in% c('Ia', 'Ib', 'Ic'))) %>% #Alternatively select Obobogo
#   slice(1L)
# #Additional note: Recent linguistic origin used ('Exploring the relationships between genetic, linguistic and geographic distances in Bantu-speaking populations', Gonzalez-Santos, 2022) was that of the Lemande population: (lat, long) = (4.50, 11.08)

# Possible start-point (oldest date) in easter_EIA dataset
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

#Sampling window: Sub-Saharan Africa ----
#sampling_win <- st_union(ne_countries(country = subSahara_countries, returnclass = "sf") %>%
# filter(name_en %!in% c("Madagascar","Sudan"))) #We focus on mainland sub-Saharan Africa (also, there is something wrong with the geometry of Sudan -- remove country since we have no iron age dates there anyway)

#Sampling window: Eastern Sub-Saharan Africa ----
sampling_win <- ne_countries(continent = "Africa", country = eastEIA_countries, returnclass = "sf", scale="large") #the detailed resolution (scale='large') ensures the Comoros islands are included

#Generate Spatial Hexagons ---- ##see code block below to determine hex diameter, cell_d 
hex_area_win <- hex_areas(sampling_win, cell_d = 3.8)

#Remove spatial hexagons where the Bantu Expansion didn't reach
#Sub-Saharan Africa, cell-diameter 6.2
# hex_area_win <- hex_area_win %>%
#   filter(area_ID %!in% c(13, 25, 12, 24, 18, 30, 41, 3, 11, 23, 35, 46, 54, 61, 65, 66, 64, 57, 63)) %>%
#   mutate(area_ID = row_number())
#Sub-Saharan Africa, cell-diameter 5.2
# hex_area_win <- hex_area_win %>%
#  filter(area_ID %!in% c(18, 25, 39, 24, 17, 38, 31, 45, 51, 23, 37, 4, 11, 88, 86, 87, 84, 83, 81, 73, 77, 68, 63)) %>%
#  mutate(area_ID = row_number())
#Sub-Saharan Africa, cell-diameter 7.2
# hex_area_win <- hex_area_win %>% 
#  filter(area_ID %!in% c(1, 3, 9, 19, 30, 35, 25, 14, 20, 40, 48, 53, 55, 54, 51, 44, 10, 15)) %>%
#  mutate(area_ID = row_number())
#East Africa, cell-diameter 2.7
# hex_area_win <- hex_area_win %>%
#   filter(area_ID %!in% c(1,2,3,4,5,9,6,10,15,20,97,122,127,121,116)) %>%
#   mutate(area_ID = row_number())
# #East Africa, cell-diameter 3.8
hex_area_win <- hex_area_win %>%
  filter(area_ID %!in% c(1,2,3,4,5,7,11,15,47,58)) %>%
  mutate(area_ID = row_number())
# #East Africa, cell-diameter 5.7
# hex_area_win <- hex_area_win %>%
#   filter(area_ID %!in% c(1,2,3,6,30)) %>%
#   mutate(area_ID = row_number())

##CHECK -- plot hexs and sites
# ggplot(data = hex_area_win) +
#   geom_sf(data = st_buffer(st_as_sf(sampling_win, crs = 4326), 40000), aes(color = "grey50")) + #sampling window with coastal buffer
#   geom_sf() + #hex grid
#   geom_sf(data = as(sites, 'sf'), size=2, alpha=0.5) + #sites
#   geom_sf_label(aes(label = area_ID),size=3) + #hex grid labels
#   #geom_sf(data = hex_area_win$area_center, size=2, alpha=1, aes(color = "purple")) + #hex-origins
#   theme(panel.background = element_rect(fill = "lightblue",
#                                         colour = "lightblue",
#                                         size = 0.5,
#                                         linetype = "solid"),
#         legend.position = "none")

#Assign hex area id to each site ----
siteInfo$area_id <- as.integer(st_within(sites$geometry, hex_area_win$geometry))
siteInfo$area_id[siteInfo$siteName=="Dembeni"] <- 65 #Impute area_ID for site=Dembeni on the far east of Comoros islands


#Assign hex area id to each date ----
dateInfo$area_id <- siteInfo$area_id[match(dateInfo$siteID, siteInfo$siteID)]

# #CHECK ---
area_freq  <- plyr::count(siteInfo, 'area_id') ##See how many sites fall in each hex area. Also make sure there are no 'NA' entries
#In order to check that this lines up visually with how many sites are in each hex area see map_figure2

#--------------------------------
## Determining hex size ---
#Under changing hex size, determine the proportion of areal hex units in the sampling window with sites
# prop_units_df <- data.frame(d = numeric(), prop_with_sites = numeric())
# 
# for (d in seq(1, 15, 0.1)){
#   hex_area_win <- hex_areas(sampling_win, cell_d = d)
#   siteInfo$area_id <- as.integer(st_within(sites$geometry, hex_area_win$geometry))
# 
#   hex_with_sites <- length(unique(siteInfo$area_id))
#   all_hex <- length(hex_area_win$area_ID)
# 
#   prop_with_sites <- hex_with_sites/all_hex
# 
#   prop_units_df <- rbind(prop_units_df, data.frame(d = d, prop_with_sites = prop_with_sites))
# }
# 
# # Plot results
# pdf(here('output','figures','figure_hexsize2.pdf'),height=5,width=5.5)
# ggplot(prop_units_df, aes(x = d, y = prop_with_sites)) +
#   geom_line() +
#   geom_point() +
#   scale_x_continuous(breaks=seq(0,15,by=1))+
#   labs(x = "Hexagon Size (d)", y = "Proportion of Hexagons with Sites", title = "Effect of Hexagon Size on Site Coverage") +
#   theme_minimal()
# dev.off()
# #-------------------------------------------------------------------------------
## Create list with constants and data ----

# Data
eastEIA_dat <- list(cra=dateInfo$cra,
                    cra_error=dateInfo$cra_error) #Creating eastEIA_dat with only date information

## Constants
data(intcal20)
data(shcal20)
constants <- list()
constants$countries <- subSahara_countries 
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
save(sites, constants, eastEIA_dat, siteInfo, dateInfo, sampling_win, hex_area_win, file=here('data','eastc14.RData')) #c14.RData

#-------------------------------------------------------------------------------
## Save sampling window specific information separately ----
 
constants_sw <- list()
constants_sw$countries <- subSahara_countries
constants_sw$eastEIAcountries <- eastEIA_countries
constants_sw$n_areas  <- constants$n_areas
constants_sw$calBP <- constants$calBP
constants_sw$C14BP  <- constants$C14BP 
constants_sw$C14err  <- constants$C14err

save(constants_sw, sampling_win, hex_area_win, file=here('data','sample_window.RData'))

#-------------------------------------------------------------------------------
## Elevation Data

country_codes <- country_codes() %>% filter(NAME %in% eastEIA_countries) #obtain country codes 

#Import elevation data
SRTM90m <- elevation_30s(country_codes$ISO3[1], path=here('input'), mask=TRUE)
for (i in 2:nrow(country_codes)){
  SRTM90m <- merge(SRTM90m, elevation_30s(country_codes$ISO3[i], path=here('input'), mask=TRUE))
}

plot(SRTM90m)
#plot(hex_area_win$geometry, add = T)

#Add area IDs
mean_hex_elv <- data.frame(area_ID = hex_area_win$area_ID,
                           mean_elevation = terra::zonal(SRTM90m, terra::vect(hex_area_win), fun = "mean", na.rm = TRUE)) #calculate mean elevation in each hexagon

#Impute missing values from nearest neighbors 
mean_hex_elv[mean_hex_elv$area_ID==17, ] = c(17, mean_hex_elv[mean_hex_elv$area_ID==8, ]$BWA_elv_msk)
mean_hex_elv[mean_hex_elv$area_ID==47, ] = c(47, mean_hex_elv[mean_hex_elv$area_ID==44, ]$BWA_elv_msk)
mean_hex_elv[mean_hex_elv$area_ID==45, ] = c(45, mean_hex_elv[mean_hex_elv$area_ID==48, ]$BWA_elv_msk)
mean_hex_elv[mean_hex_elv$area_ID==37, ] = c(37, mean_hex_elv[mean_hex_elv$area_ID==32, ]$BWA_elv_msk)

#Normalise elevation to [0,1] scale
mean_hex_elv <- mean_hex_elv %>% 
  mutate(norm_mean_elv = BWA_elv_msk/max(BWA_elv_msk))

#Save elevation output
save(mean_hex_elv, SRTM90m, file=here('data','elevation.RData'))

#-------------------------------------------------------------------------------
## Crop suitability data

#Read in crop data ----
pearl_millet_sf <- read_stars("data/environment/Soil_suitability_Chemura/Cereals/Suit/pmillet_curr.tif") %>% 
  st_as_sf() %>% 
  rename("pmillet_suit" = "pmillet_curr.tif")

sorghum_sf <- read_stars("data/environment/Soil_suitability_Chemura/Cereals/Suit/sorghum_curr.tif") %>% 
  st_as_sf() %>% 
  rename("sorghum_suit" = "sorghum_curr.tif")

#--------------  
#Aggregate suitability values for each hexagonal area ----

#Assign hex area id
pearl_millet_sf$area_id <- as.integer(st_within(pearl_millet_sf$geometry, hex_area_win$geometry))

sorghum_sf$area_id <- as.integer(st_within(sorghum_sf$geometry, hex_area_win$geometry))


#Filter for suitability values within the sample window
pearl_millet_sf <- pearl_millet_sf %>% 
  filter(!is.na(area_id))

sorghum_sf <- sorghum_sf %>% 
  filter(!is.na(area_id))

#------
#Aggregate
pearl_millet_df <- pearl_millet_sf %>% 
  group_by(area_id) %>% 
  summarize(mean_pmillet_suit = mean(pmillet_suit)) %>% 
  as.data.frame() %>% 
  dplyr::select(area_id, mean_pmillet_suit) %>% 
  na.omit()

sorghum_df <- sorghum_sf %>% 
  group_by(area_id) %>% 
  summarize(mean_sorghum_suit = mean(sorghum_suit)) %>% 
  as.data.frame() %>% 
  dplyr::select(area_id, mean_sorghum_suit) %>% 
  na.omit()

#--------------  
#Impute values from neighbors for hex areas that are too small to have suitability values ----
#NB: This needs to be changed for different sample windows!
pearl_millet_df[nrow(pearl_millet_df) + 1,] = c(16, pearl_millet_df[pearl_millet_df$area_id==27, ]$mean_pmillet_suit)
pearl_millet_df[nrow(pearl_millet_df) + 1,] = c(17, pearl_millet_df[pearl_millet_df$area_id==12, ]$mean_pmillet_suit)
pearl_millet_df[nrow(pearl_millet_df) + 1,] = c(38, pearl_millet_df[pearl_millet_df$area_id==29, ]$mean_pmillet_suit)
pearl_millet_df[nrow(pearl_millet_df) + 1,] = c(45, pearl_millet_df[pearl_millet_df$area_id==41, ]$mean_pmillet_suit)
pearl_millet_df[nrow(pearl_millet_df) + 1,] = c(46, pearl_millet_df[pearl_millet_df$area_id==39, ]$mean_pmillet_suit)
pearl_millet_df[nrow(pearl_millet_df) + 1,] = c(47, pearl_millet_df[pearl_millet_df$area_id==44, ]$mean_pmillet_suit)

sorghum_df[nrow(sorghum_df) + 1,] = c(16, sorghum_df[sorghum_df$area_id==27, ]$mean_sorghum_suit)
sorghum_df[nrow(sorghum_df) + 1,] = c(17, sorghum_df[sorghum_df$area_id==12, ]$mean_sorghum_suit)
sorghum_df[nrow(sorghum_df) + 1,] = c(38, sorghum_df[sorghum_df$area_id==29, ]$mean_sorghum_suit)
sorghum_df[nrow(sorghum_df) + 1,] = c(45, sorghum_df[sorghum_df$area_id==41, ]$mean_sorghum_suit)
sorghum_df[nrow(sorghum_df) + 1,] = c(46, sorghum_df[sorghum_df$area_id==39, ]$mean_sorghum_suit)
sorghum_df[nrow(sorghum_df) + 1,] = c(47, sorghum_df[sorghum_df$area_id==44, ]$mean_sorghum_suit)


pearl_millet_df <- pearl_millet_df %>% arrange(area_id) 
sorghum_df <- sorghum_df %>% arrange(area_id) 

#Plot
#Pearl Millet
ggplot(data = hex_area_win) +
  geom_sf(data = st_buffer(st_as_sf(sampling_win, crs = 4326), 40000), color = "grey50") + #sampling window with coastal buffer
  geom_sf(aes(fill = pearl_millet_df$mean_pmillet_suit)) +
  scale_fill_gradientn(colours = rev(terrain.colors(7)), name = "Pearl Millet Suitability") +
  theme(panel.background = element_rect(fill = "lightblue",
                                        colour = "lightblue",
                                        size = 0.5,
                                        linetype = "solid"))
# ggplot(data = pearl_millet_sf) +
#   geom_sf(data = st_buffer(st_as_sf(sampling_win, crs = 4326), 40000), color = "grey50") + #sampling window with coastal buffer
#   geom_sf(aes(fill = pmillet_suit)) +
#   scale_fill_gradientn(colours = rev(terrain.colors(7)), name = "Pearl Millet Suitability") +
#   theme(panel.background = element_rect(fill = "lightblue",
#                                         colour = "lightblue",
#                                         size = 0.5,
#                                         linetype = "solid"))

#Sorghum
ggplot(data = hex_area_win) +
  geom_sf(data = st_buffer(st_as_sf(sampling_win, crs = 4326), 40000), color = "grey50") + #sampling window with coastal buffer
  geom_sf(aes(fill = sorghum_df$mean_sorghum_suit)) +
  scale_fill_gradientn(colours = rev(terrain.colors(7)), name = "Sorghum Suitability") +
  theme(panel.background = element_rect(fill = "lightblue",
                                        colour = "lightblue",
                                        size = 0.5,
                                        linetype = "solid"))
# ggplot(data = sorghum_sf) +
#   geom_sf(data = st_buffer(st_as_sf(sampling_win, crs = 4326), 40000), color = "grey50") + #sampling window with coastal buffer
#   geom_sf(aes(fill = sorghum_suit)) +
#   scale_fill_gradientn(colours = rev(terrain.colors(7)), name = "Sorghum Suitability") +
#   theme(panel.background = element_rect(fill = "lightblue",
#                                         colour = "lightblue",
#                                         size = 0.5,
#                                         linetype = "solid"))

#--------------
#Join crop data ----
agg_crop_suitability <- pearl_millet_df %>% 
  left_join(sorghum_df, by=join_by(area_id)) %>% 
  rowwise() %>% 
  mutate(max_crop_suit = as.numeric(max(mean_pmillet_suit, mean_sorghum_suit))) %>% 
  select(area_id, max_crop_suit)

#--------------
#Save crop suitability output ----
save(agg_crop_suitability, file=here('data','crop_suitability.RData'))


#----------------- ALTERNATIVE AGRICULTURAL SUITABILITY ------------------------

## Sedentary animal husbandry suitability data

#Read in animal husbandry data ----
agr_suit_sf <- read_stars("data/environment/Agriculture_suitability_Beck/agr_suit.asc") %>% 
  st_as_sf() %>% 
  rename("agr_suit" = "agr_suit.asc")

st_crs(agr_suit_sf)  <- 4326 
#--------------  
#Aggregate suitability values for each hexagonal area ----

#Assign hex area id
agr_suit_sf$area_id <- as.integer(st_within(agr_suit_sf$geometry, hex_area_win$geometry))

#Filter for suitability values within the sample window
agr_suit_sf <- agr_suit_sf %>% 
  filter(!is.na(area_id))

#------
#Aggregate
agr_suit_df <- agr_suit_sf %>% 
  group_by(area_id) %>% 
  summarize(mean_agr_suit = mean(agr_suit)) %>% 
  as.data.frame() %>% 
  dplyr::select(area_id, mean_agr_suit) %>% 
  na.omit()

#--------------  
#Impute values from neighbors for hex areas that are too small to have suitability values ----
#NB: This needs to be changed for different sample windows!
agr_suit_df[nrow(agr_suit_df) + 1,] = c(17, agr_suit_df[agr_suit_df$area_id==12, ]$mean_agr_suit)
agr_suit_df[nrow(agr_suit_df) + 1,] = c(47, agr_suit_df[agr_suit_df$area_id==44, ]$mean_agr_suit)

agr_suit_df <- agr_suit_df %>% arrange(area_id) 

#Plot
#Agricultural suitability
ggplot(data = hex_area_win) +
  geom_sf(data = st_buffer(st_as_sf(sampling_win, crs = 4326), 40000), color = "grey50") + #sampling window with coastal buffer
  geom_sf(aes(fill = agr_suit_df$mean_agr_suit)) +
  scale_fill_gradientn(colours = rev(terrain.colors(7)), name = "Agriculture Suitability") +
  theme(panel.background = element_rect(fill = "lightblue",
                                        colour = "lightblue",
                                        size = 0.5,
                                        linetype = "solid"))
# ggplot(data = agr_suit_sf) +
#   geom_sf(data = st_buffer(st_as_sf(sampling_win, crs = 4326), 40000), color = "grey50") + #sampling window with coastal buffer
#   geom_sf(aes(fill = agr_suit)) +
#   scale_fill_gradientn(colours = rev(terrain.colors(7)), name = "Agriculture Suitability") +
#   theme(panel.background = element_rect(fill = "lightblue",
#                                         colour = "lightblue",
#                                         size = 0.5,
#                                         linetype = "solid"))


#Save alternative agriculutural suitability output ----
save(agr_suit_df, file=here('data','crop_suitability2.RData'))

#-------------------------------------------------------------------------------
## Sedentary animal husbandry suitability data

#Read in animal husbandry data ----
animal_hus_sf <- read_stars("data/environment/AnimalHusbandry_suitability_Beck/anim_suit.asc") %>% 
  st_as_sf() %>% 
  rename("animal_hus_suit" = "anim_suit.asc")

st_crs(animal_hus_sf)  <- 4326 
#--------------  
#Aggregate suitability values for each hexagonal area ----

#Assign hex area id
animal_hus_sf$area_id <- as.integer(st_within(animal_hus_sf$geometry, hex_area_win$geometry))

#Filter for suitability values within the sample window
animal_hus_sf <- animal_hus_sf %>% 
  filter(!is.na(area_id))

#------
#Aggregate
amimal_hus_df <- animal_hus_sf %>% 
  group_by(area_id) %>% 
  summarize(mean_animal_hus_suit = mean(animal_hus_suit)) %>% 
  as.data.frame() %>% 
  dplyr::select(area_id, mean_animal_hus_suit) %>% 
  na.omit()

#--------------  
#Impute values from neighbors for hex areas that are too small to have suitability values ----
#NB: This needs to be changed for different sample windows!
amimal_hus_df[nrow(amimal_hus_df) + 1,] = c(17, amimal_hus_df[amimal_hus_df$area_id==12, ]$mean_animal_hus_suit)
amimal_hus_df[nrow(amimal_hus_df) + 1,] = c(47, amimal_hus_df[amimal_hus_df$area_id==44, ]$mean_animal_hus_suit)

amimal_hus_df <- amimal_hus_df %>% arrange(area_id) 

#Plot
#Animal Husbandry
ggplot(data = hex_area_win) +
  geom_sf(data = st_buffer(st_as_sf(sampling_win, crs = 4326), 40000), color = "grey50") + #sampling window with coastal buffer
  geom_sf(aes(fill = amimal_hus_df$mean_animal_hus_suit)) +
  scale_fill_gradientn(colours = rev(terrain.colors(7)), name = "Animal Husbandry Suitability") +
  theme(panel.background = element_rect(fill = "lightblue",
                                        colour = "lightblue",
                                        size = 0.5,
                                        linetype = "solid"))
# ggplot(data = animal_hus_sf) +
#   geom_sf(data = st_buffer(st_as_sf(sampling_win, crs = 4326), 40000), color = "grey50") + #sampling window with coastal buffer
#   geom_sf(aes(fill = animal_hus_suit)) +
#   scale_fill_gradientn(colours = rev(terrain.colors(7)), name = "Animal Husbandry Suitability") +
#   theme(panel.background = element_rect(fill = "lightblue",
#                                         colour = "lightblue",
#                                         size = 0.5,
#                                         linetype = "solid"))


#Save animal husbandry suitability output ----
save(amimal_hus_df, file=here('data','animal_hus_suitability.RData'))


