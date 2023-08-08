#install.packages("remotes")
#remotes::install_github("people3k/p3k14c@2022.06")

# Load Libraries and Data ----
library(rcarbon)
library(nimbleCarbon)
library(maptools)
library(sf)
library(rnaturalearth)
library(rgeos)
library(stringr)
library(dplyr)
library(tidyr)
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

##East African data
#Note: There is potential to use these dates, but most seem to be geological (and not archaeological) in origin. Also seems unclear where this database has been used in research subsequent to its creation
#EA_dat contains data from the York Institute for Tropical Ecosystems dataverse
#Downloaded from: https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/NJLNRJ 
#See Mustaphi et al., 2016, "Radiocarbon dates from eastern Africa in the CARD2.0 format"

##Russell Earliest Bantu dates
#Note: This database has already been agregated at the site level (keeping the earliest dates at each site) #TODO: Obtain unagregated site from Thembi Russell 
#Downloaded from: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0087854 (Can also be downloaded from: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0215573)
#See Russell et al., 2014, "Modelling the spread of farming in the Bantu-Speaking regions of Africa"

##Kay West Africa data
#KayWA_dat contains data collected in West Africa associated with LandUse6K
#Summary data can be downloaded from: https://link.springer.com/article/10.1007/s10963-019-09131-2#Sec32 
#See Kay et al., 2019, "Diversification, Intensification and Specialization: Changing Land Use in Western Africa from 1800 BC to AD 1500"
#Radiocarbon dates are included in p3k14c database, but complete raw data isn't seperately available yet.


##p3k14c
#p3k14c_dat contains data from the p13k14c Global Archaeological Radiocarbon Database
#Download from: https://github.com/people3k/p3k14c/ [Raw data is accessed through tDAR: https://core.tdar.org/dataset/459173/p3k14c-version-202201-raw-data]
#Note this is an overview database containing the following African databases: Kay's West Africa database, aDRAC, CalPal.... 
#See Bird et al., 2022, "p3k14c, a synthetic global database of archaeological radiocarbon dates" (https://www.nature.com/articles/s41597-022-01118-7)


#-------------------------------------------------------------------------------
## Read data sets  ----
SARD_dat <- read.csv(here("data", "SARD_Mar2021_14C.csv"))
HumActCA_dat <- read.csv(here("data", "HumActCA_Dec2020_14C.csv")) #Alternatively read in "aDRAC_April2023_14C.csv" -- this is a more up to date version from the aDRAC GitHub repository, but then we need to also filter for Bantu associated dates
Russell_EIA_dat <- read.csv(here("data", "Russell_EarliestBantu_Jan2014_14C.csv")) 
#EA_dat <- read.csv(here("data", "EastAfrica_April2016_14C.csv"))
#KayWA_dat <- read.csv(here("data", "KayWA_SUMMARY_May2019_14C.csv")) #TODO: Replace this summary data with actual radiocarbon data when available --- for now, p3k14c contains these dates

p3k14c_dat <- read.csv(here("data", "p3k14c_raw.csv")) #p3k14c::p3k14c_data will give the scrubbed and fuzzed p3k14c dataset

#-------------------------------------------------------------------------------
## Data filtering and cleaning ----

SARD_sum_df <- SARD_dat %>%
  filter(Archaeological.Period=="Iron Age") %>% 
  dplyr::select(Lab.ID, X.Site, DecdegE, DecdegS, Date, Uncertainty) %>%
  rename(labCode=Lab.ID, siteName=X.Site, lat=DecdegS, long=DecdegE, c14date=Date, c14std=Uncertainty) %>%
  mutate(c14date = as.numeric(c14date), c14std=as.numeric(c14std), dataorigin="SARD")  %>%
  filter(!is.na(lat) & !is.na(long) & !is.na(c14std) & !is.na(c14date)) %>% 
  filter(siteName !="Bambata Cave") #Designated pre-bantu (references given in Isern and Fort 2019, Suplementary Material S1, https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0215573) #TODO: decide if Shongweni North and South dates are to be removed? Isern and Fort also mention removing Shongweni Waterworks Park (2030 uncal BP) as being pre-bantu

Russell_sum_df <- Russell_EIA_dat %>%
  dplyr::select(Site, Latitude, Longitude, Uncal.BP, St.Dev) %>%
  rename(siteName=Site, lat=Latitude, long=Longitude, c14date=Uncal.BP, c14std=St.Dev) %>%
  mutate(c14date = as.numeric(c14date), c14std=as.numeric(c14std), dataorigin="RussellEIA")  %>%
  filter(!is.na(lat) & !is.na(long) & !is.na(c14std) & !is.na(c14date)) %>% 
  filter(siteName %!in% c("Bambata Cave Series", "Nakapapula", "Kwelikwiji", "Kumadzulo", "Ndonde", "Shongweni Waterworks Park")) %>% #Designated pre-bantu (references given in Isern and Fort 2019, Suplementary Material S1, https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0215573)
  filter(siteName %!in% c('Bisoli', 'Bokele', 'Boma', 'Broadhurst', 'Broederstroom, 24/73 K', 'Campo', 'Chibuene', 'Ditouba', 'Djambala', 'Enkwazini', 'Hola hola', 'Imbonga', 'Kinsale Farm', 'Klein Afrika', 'Langubhela', 'Mabveni', 'Makokou', 'Malemba-Nkulu', 'Maluba', 'Massingir', 'Matola IV', 'Mbilap_ 4', 'Mpassa farm', 'Munda', 'Ntadi ntadi Cave', 'Ntsitsana, Pit 1', 'Oyem 2', 'Pikunda', 'Silver Leaves', 'Taukome (lower levels)', 'Toteng I', 'Tchissanga Ouest')) #Duplicate site entries already present in aDRAC and SARD databases. We prefer the entries in those databases, since they have labcodes attached and are unagregated
Russell_sum_df[,"labCode"] <- NA #No Lab Code present in this site level dataset. Perhaps in the unagregated df?

HumActCA_sum_df <- HumActCA_dat %>% #I think all the sites in this database are associated with Bantu pottery finds
  filter(CLASS %!in% c('IIa', 'IIb', 'IIc', 'IIIa', 'IIIb', 'IIIc')) %>% #Filter for unreliable class of dates III and irrelevant (according to Seidensticker et al.) dates II (TODO: check whether class II dates should actually be excluded)
  dplyr::select(LABNR, SITE, LAT, LONG, C14AGE, C14STD) %>%
  rename(labCode= LABNR, siteName=SITE, lat=LAT, long=LONG, c14date=C14AGE, c14std=C14STD) %>%
  mutate(c14date = as.numeric(c14date), c14std=as.numeric(c14std), dataorigin="HumActCA") %>%
  filter(!is.na(lat) & !is.na(long) & !is.na(c14std) & !is.na(c14date)) 

# EA_sum_df <- EA_dat %>%  #TODO: filter for Bantu dates
#   filter(Type.of.Date == "Archaeological") %>% 
#   select(Lab.Number, Site.Identifier, Latitude, Longitude, Normalized.Age, NA.Sigma) %>%
#   rename(labCode= Lab.Number, siteName=Site.Identifier, lat=Latitude, long=Longitude, c14date=Normalized.Age, c14std=NA.Sigma) %>%
#   mutate(c14date = as.numeric(c14date), c14std=as.numeric(c14std), dataorigin="EastAfrican") %>%
#   filter(!is.na(lat) & !is.na(long) & !is.na(c14std) & !is.na(c14date))

# KayWA_sum_df <- KayWA_dat %>%  #TODO: filter for Bantu dates
#   select(Site.Name, X, Y) %>%
#   rename(site=Site.Name, lat=X, long=Y) %>%
#   mutate(c14date = as.numeric(c14date), c14std=as.numeric(c14std), dataorigin="KayWestAfrican") 


p3k14c_sum_df <- p3k14c_dat %>%
  filter(Continent=="Africa" & 
           (Country %in% c(subSahara_countries, "CAR", "Equat.Guinea", "DRC")) & #TODO: refine the determination of whether the site is Bantu at a later stage -- specifically filter for the correct period
           (Age != "Bomb C14")) %>% #Just cleaning up...
  dplyr::select(LabID, SiteName, Lat, Long, Age, Error, Source) %>%
  rename(labCode = LabID, siteName=SiteName, lat=Lat, long=Long, c14date=Age, c14std=Error, source=Source) %>%
  mutate(lat=as.numeric(lat), long=as.numeric(long), c14date = as.numeric(c14date), c14std=as.numeric(c14std), dataorigin="p3k14c") %>%
  filter(!is.na(lat) & !is.na(long) & !is.na(c14std) & !is.na(c14date)) 
#The sites without co-ordinates come from the following sources: aDRAC, CALPAL, Kay_WestAfrica, Vermeersch2019. There are 95 of these sites in total. We remove them. #TODO:Circle back to this and see if these dates can be found in original references.
#The sites without C14 Errors come from the following sources: Kay_WestAfrica and SARD. There are 46 of these sites in total. We remove them. #TODO:Circle back to this and see if these dates can be found in original references. 
#The sites without C14 Ages all come from the SARD database. There are 33 of these sites in total. We remove them. #TODO:Circle back to this and see if these dates can be found in original references. 


#Note, since p3k14c is a database of databases it contains a fair amount of information already present in the other databases, such as SARD. We want to use the original databases as far as possible. 
p3k14c_fil_df <- p3k14c_sum_df %>%
  filter(!(source %in% c("aDRAC", "SARD", "RussellEIA"))) %>%  # This represents the HumActCA_dat, SARD_dat, and Russell_sum_df databases respectively
  dplyr::select(-source)
#Note there are 7 sites which are excluded (above) from aDRAC (being class II or class III) dates, which p3k14c adds back into the dataset here. TODO: take a look at this more closely... 

# Combine datasets ----
bantu_sites_df <- bind_rows(SARD_sum_df, HumActCA_sum_df, Russell_sum_df) %>% #TODO: Might need more careful filtering of Bantu dates... 
  #bind_rows(p3k14c_fil_df, SARD_sum_df, HumActCA_sum_df, Russell_sum_df) %>%
  mutate(dataorigin=as.factor(dataorigin)) %>% 
  filter((c14date != 0) & (c14std != 0)) %>% #Apparently some of these datasets had modern dates (indicated with c14 date and error of 0), we remove these
  filter((c14date <=3357) & (c14date >=246)) #We assume an approximate origin at Ngoume, date 3357 +- 95BP (see later in this script). Dates earlier than this are assumed to not be of Bantu origin #TODO: build some flexibility into this... #Further, we assume the Dutch arrival in the Cape (1652) as the cut-off. Dates after this point of colonial contact are not considered.
  #filter((c14date <=3070) & (c14date >=298)) #We assume an approximate origin at Obobogo, date 3070 +- 95BP (see later in this script). Dates earlier than this are assumed to not be of Bantu origin #TODO: build some flexibility into this... #Further, we assume the Dutch arrival in the Cape (1652) as the cut-off. Dates after this point of colonial contact are not considered.

  
# Assign ID ----
bantu_sites_df <- bantu_sites_df %>% mutate(ID = row_number())

# Assign Site ID ----
bantu_sites_df$siteID  <- as.numeric(factor(bantu_sites_df$siteName))

#-------------------------------------------------------------------------------
## Determining which calibration curve should be used----

# Assign calibration curve ----
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
sites <- siteInfo
coordinates(sites) <- c('long','lat')
proj4string(sites)  <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
dist_mat  <- spDists(sites,longlat=TRUE) #inter-site distance matrix: each site's distance from every other site (i.e. with n sites, this matrix is n^2)
origin_point  <- c(possible_origin_dat$LONG, possible_origin_dat$LAT) #Ngoume, HumActCA Origin
#origin_point  <- c(possible_origin_dat$long, possible_origin_dat$lat) #SARD Origin
dist_org  <-  spDistsN1(sites, origin_point, longlat=TRUE) #distance from origin site


#-------------------------------------------------------------------------------
# Generate Spatial Window for Analyses: Sub-Saharan Africa ----

#Sampling window ----
sf_subsah_africa <- ne_countries(continent = "Africa", returnclass = "sf") %>%
  filter_all(., any_vars(str_detect(., "Sub-Saharan"))) %>% 
  filter(name_en %in% subSahara_countries) %>% 
  filter(name_en != "Madagascar") #We focus on mainland sub-Saharan Africa

sampling_win <- sf_subsah_africa %>% as("Spatial") #convert sf to sp object


#Generate Hex Areas over Spatial Window ----
hex_area_win <- hex_areas(sampling_win, cell_d = 7)
  
#Assign hex area id to each site ----
sites_sf <- as(sites, 'sf')
siteInfo$area_id <- as.integer(st_within(sites_sf$geometry, hex_area_win$geometry))

# #CHECK ---
# area_freq  <- plyr::count(siteInfo, 'area_id') ##See how many sites fall in each hex area. Also make sure there are no 'NA' entries
# To check that this lines up visually with how many sites are in each hex area -- see map_figure2

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
