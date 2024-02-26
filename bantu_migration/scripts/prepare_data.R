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
library(p3k14c)

rm(list = ls())

ncores = (detectCores() - 1)

`%!in%` <- Negate(`%in%`)

source(here('src','hex_areas.R'))

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

#-------------------------------------------------------------------------------
## Read data sets  ----
SARD_dat <- read.csv(here("data", "SARD_Mar2021_14C.csv"))
HumActCA_dat <- read.csv(here("data", "HumActCA_Dec2020_14C.csv"))
aDRAC_dat <- read.csv(here("data", "aDRAC_Feb2024_14C.csv"))#This is a more up-to-date version of HumActCA_dat, but not all the dates in aDRAC_dat are iron age, whereas all dates in the HumActCA are
Collected_EIA_dat <- read_excel(here("data", "collected_dates.xlsx"), sheet = "Collected_dates_final") 

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
SARD_sum_df <- SARD_dat %>%
  filter(Archaeological.Period=="Iron Age") %>% 
  dplyr::select(Lab.ID, X.Site, DecdegE, DecdegS, Date, Uncertainty, refcode, Material.dated, Country) %>%
  rename(labCode=Lab.ID, siteName=X.Site, lat=DecdegS, long=DecdegE, c14date=Date, c14std=Uncertainty, reference=refcode, material=Material.dated, country=Country) %>%
  mutate(c14date = as.numeric(c14date), c14std=as.numeric(c14std), dataorigin="SARD")  %>%
  filter(!is.na(lat) & !is.na(long) & !is.na(c14std) & !is.na(c14date)) %>% 
  filter(siteName !="Bambata Cave") #Designated pre-bantu (references given in Isern and Fort 2019, Suplementary Material S1, https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0215573) 

Collected_sum_df <- Collected_EIA_dat %>%
  dplyr::select(LabID, Age, Error, Long, Lat, SiteName, Reference, Material, Country) %>%
  rename(labCode=LabID, siteName=SiteName, lat=Lat, long=Long, c14date=Age, c14std=Error, reference=Reference, material=Material, country=Country) %>%
  mutate(c14date = as.numeric(c14date), c14std=as.numeric(c14std), dataorigin="Collected")  %>%
  filter(!is.na(lat) & !is.na(long) & !is.na(c14std) & !is.na(c14date)) 

aDRAC_sum_df <- aDRAC_dat %>% #To find extra dates not in HumActCA
  dplyr::select(LABNR, SITE, LAT, LONG, C14AGE, C14STD, PHASE, CLASS, IRON, SOURCE, MATERIAL, COUNTRY) %>%
  rename(labCode= LABNR, siteName=SITE, lat=LAT, long=LONG, c14date=C14AGE, c14std=C14STD, reference=SOURCE, material=MATERIAL, country=COUNTRY) %>%
  filter(!is.na(lat) & !is.na(long) & !is.na(c14std) & !is.na(c14date)) %>%
  mutate(c14date = as.numeric(c14date), c14std=as.numeric(c14std), dataorigin="aDRAC") %>%
  filter(siteName %in% HumActCA_dat$SITE | PHASE %in% c("N; EIA", "LIA", "EIA", "Iron Age") | !is.na(IRON)) %>% #All dates in HumActCA_dat should be associated with Bantu pottery finds
  filter(CLASS %!in% c('IIa', 'IIb', 'IIc', 'IIIa', 'IIIb', 'IIIc')) %>% #Filter for unreliable class of dates III and irrelevant (according to Seidensticker et al.) dates II
  dplyr::select(-PHASE, -CLASS, -IRON)

##------------
## Combine datasets ----
bantu_sites_df <- bind_rows(SARD_sum_df, aDRAC_sum_df, Collected_sum_df) %>% 
  mutate(dataorigin=as.factor(dataorigin)) %>% 
  filter((c14date != 0) & (c14std != 0)) %>% #Apparently some of these datasets had modern dates (indicated with c14 date and error of 0), we remove these
  filter((c14date <=3357) & (c14date >=246)) #We assume an approximate origin at Ngoume, date 3357 +- 95BP (see later in this script). Dates earlier than this are assumed to not be of Bantu origin #TODO: build some flexibility into this... #Further, we assume the Dutch arrival in the Cape (1652) as the cut-off. Dates after this point of colonial contact are not considered.
  
# Assign ID ----
bantu_sites_df <- bantu_sites_df %>% mutate(ID = row_number())

# Assign Site ID ----
bantu_sites_df$siteID  <- as.numeric(factor(bantu_sites_df$siteName))

##-----------
##Save Output as csv
write.csv(bantu_sites_df, here('data','bantu_dataset.csv'), row.names = FALSE)

#===============================================================================
## Determining which calibration curve should be used----

#Remove unnecessary information
bantu_sites_df <- bantu_sites_df %>%
  dplyr::select(-reference, -material, -country)

#Assign calibration curve ----
bantu_sites_df$calCurve <- ifelse((bantu_sites_df$lat>=0), 'intcal20', 'shcal20') #Assign the calibration curve to use based on the site's position relative to the equator #TODO: refine this -- weird to have a hard step-change between calibration curves at the equator -- maybe use a gradient change function? Mixed Curve? Or is there perhaps better regional calibration curves to use?
#table(bantu_sites_df$calCurve)

#-------------------------------------------------------------------------------
## Restructure Data for Bayesian Analyses ----

# Compute median calibrated dates ----
bantu_sites_df$median_dates = medCal(calibrate(bantu_sites_df$c14date,
                                               bantu_sites_df$c14std,
                                               calCurve = bantu_sites_df$calCurve,
                                               ncores = ncores))


# Collect site level information ----
earliest_dates <- aggregate(median_dates ~ siteID, data = bantu_sites_df, FUN = max) #Earliest medCal Date for Each Site 
latest_dates <- aggregate(median_dates ~ siteID, data=bantu_sites_df, FUN=min) #Latest medCal Date for Each Site
n_dates <- aggregate(median_dates ~ siteID, data=bantu_sites_df, FUN=length) #Number of medCal Date for Each Site

siteInfo <- data.frame(siteID = earliest_dates$siteID,
                       earliest = earliest_dates$median_dates,
                       latest = latest_dates$median_dates,
                       diff = earliest_dates$median_dates - latest_dates$median_dates,
                       n_dates = n_dates$median_dates) %>% unique()

siteInfo <- siteInfo %>% 
  left_join(unique(dplyr::select(bantu_sites_df,
                                 lat,
                                 long,
                                 siteID,
                                 siteName,
                                 dataorigin,
                                 calCurve))) %>% 
  distinct(siteID, .keep_all = TRUE)  #Some site duplicates introduced by slight differences in co-ordinate decimal places

# Collect date level information ----
dateInfo <- unique(dplyr::select(bantu_sites_df,
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

#A generally acknowledged 'homeland region on the border between Nigeria and Cameroon' (Seidensticket,et al. 2020)
possible_origin_dat <- HumActCA_dat %>%
  filter(COUNTRY %in% c('CMR', 'NGA')) %>% #I don't think there are any sites in Nigeria in this dataset, but still ...
  arrange(desc(C14AGE))

#Note: oldest class I date is Shum Laka -- but see recent genetic analysis discussing why this is a hunter-gather site ('Ancient West African foragers in the context of African population history', Lipson, 2020, https://www.nature.com/articles/s41586-020-1929-1). Also, see Piere de Maret's thoughts on the site. #TODO: Decide whether to designate pre-Bantu and exclude Shum Laka from analyses entirely...
#Therefore, as an approximate origin we select the next oldest class I date in the region: Ngoume
possible_origin_dat <- possible_origin_dat %>%
  filter(SITE=='Ngoume' & (CLASS %in% c('Ia', 'Ib', 'Ic'))) %>% #Alternatively select Obobogo
  slice(1L)
#Additional note: Recent linguistic origin used ('Exploring the relationships between genetic, linguistic and geographic distances in Bantu-speaking populations', Gonzalez-Santos, 2022) was that of the Lemande population: (lat, long) = (4.50, 11.08)

## Possible Southern Africa start-point (oldest date) in SARD dataset
# possible_origin_dat <- SARD_sum_df %>%
#   slice_max(c14date, n=1)


## Compute Great-Arc Distances in km ----
sites <- st_as_sf(siteInfo, coords = c('long','lat'))
st_crs(sites)  <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
dist_mat  <- st_distance(sites)/1000 #inter-site distance matrix in km: each site's distance from every other site (i.e. with n sites, this matrix is n^2)
origin_point  <- sites %>% filter(siteName == possible_origin_dat$SITE)
dist_org  <-  st_distance(x=sites, y=origin_point)/1000 #distance from origin site

#-------------------------------------------------------------------------------
# Generate Spatial Window for Analyses: Sub-Saharan Africa ----

#Sampling window ----
sampling_win <- ne_countries(continent = "Africa", returnclass = "sf") %>%
  filter_all(., any_vars(str_detect(., "Sub-Saharan"))) %>% 
  filter(name_en %in% subSahara_countries) %>% 
  filter(name_en != "Madagascar") #We focus on mainland sub-Saharan Africa

#Generate Spatial Hexagons ----
hex_area_win <- hex_areas(sampling_win, cell_d = 10)

#Assign hex area id to each site ----
siteInfo$area_id <- as.integer(st_within(sites$geometry, hex_area_win$geometry))

# #CHECK ---
area_freq  <- plyr::count(siteInfo, 'area_id') ##See how many sites fall in each hex area. Also make sure there are no 'NA' entries
#In order to check that this lines up visually with how many sites are in each hex area see map_figure2


#-------------------------------------------------------------------------------
## Create list with constants and data ----

# Data
bantu_dat <- list(cra=dateInfo$cra,
                  cra_error=dateInfo$cra_error) #Creating bantu_dat with only date information

## Constants
data(intcal20)
data(shcal20)
constants <- list()
constants$countries <- subSahara_countries 
constants$n_sites <- nrow(siteInfo)
constants$n_dates  <- nrow(dateInfo)
constants$n_areas  <- nrow(hex_area_win) #All areas (even empty ones) are included #Only occupied areas: length(unique(siteInfo$area_id))
constants$id_sites <- dateInfo$siteID
constants$id_area  <- siteInfo$area_id 
constants$dist_mat  <- dist_mat
constants$dist_org  <- dist_org
constants$origin_point <- origin_point
#Calibration curves
constants$calBP <- intcal20$CalBP #Same for intcal20 and shcal20 
constants$C14BP  <- cbind(intcal20$C14Age, shcal20$C14Age) #Northern and southern hemisphere calibration curves
constants$C14err  <- cbind(intcal20$C14Age.sigma, shcal20$C14Age.sigma)

#-------------------------------------------------------------------------------
## Save everything on a R image file ----
save(sites, constants, bantu_dat, siteInfo, dateInfo, sampling_win, hex_area_win, file=here('data','c14.RData'))
