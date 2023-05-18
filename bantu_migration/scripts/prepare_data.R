#install.packages("remotes")
#remotes::install_github("people3k/p3k14c@2022.06")

# Load Libraries and Data ----
library(rcarbon)
library(nimbleCarbon)
library(maptools)
library(sf)
library(stringr)
library(dplyr)
library(tidyr)
library(here)
library(ggplot2)
library(ggthemes)
library(parallel)
library(p3k14c)

#-------------------------------------------------------------------------------
##List of countries in sub-Saharan Africa
subSahara_countries <- c("South Africa", 
                         "Lesotho", 
                         "Eswatini", 
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
                         "Congo",
                         "Gabon",
                         "Cameroon",
                         "Nigeria",
                         "Equatorial Guinea",
                         "Sudan",
                         "Chad")
  
  
#-------------------------------------------------------------------------------
## Data sources  ----

##SARD
#SARD_dat contains data from the Southern African Radiocarbon Database
#Downloaded from: https://github.com/emmaloftus/Southern-African-Radiocarbon-Database
#See Loftus et al., 2019, "An archaeological radiocarbon database for southern Africa" Antiquity. (https://doi.org/10.15184/aqy.2019.75).


##aDRAC
#HumActCA_dat contains data from the Human Activity in Central Africa (Archive des datations radiocarbones d’Afrique centrale) Database 
#Downloaded from: https://zenodo.org/record/4394894 or https://github.com/dirkseidensticker/aDRAC [Latter appears to be more up to date]
#See Seidensticker et al., 2021, "Population collapse in Congo rainforest from 400 CE urges reassessment of the Bantu Expansion" (https://www.science.org/doi/full/10.1126/sciadv.abd8352)


##East African data
#EA_dat contains data from the York Institute for Tropical Ecosystems dataverse
#Downloaded from: https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/NJLNRJ 
#See Mustaphi et al., 2016, "Radiocarbon dates from eastern Africa in the CARD2.0 format"


##Kay West Africa data
#KayWA_dat contains data collected in West Africa associated with LandUse6K
#Summary data can be downloaded from: https://link.springer.com/article/10.1007/s10963-019-09131-2#Sec32 
#See Kay et al., 2019, "Diversification, Intensification and Specialization: Changing Land Use in Western Africa from 1800 BC to AD 1500"


##p3k14c
#p3k14c_dat contains data from the p13k14c Global Archaeological Radiocarbon Database
#Download from: https://github.com/people3k/p3k14c/ [If access to the raw data is required: https://core.tdar.org/dataset/459173/p3k14c-version-202201-raw-data]
#Note this is an overview database containing the following African databases: Kay's West Africa database, aDRAC, CalPal.... 
#See Bird et al., 2022, "p3k14c, a synthetic global database of archaeological radiocarbon dates" (https://www.nature.com/articles/s41597-022-01118-7)


#-------------------------------------------------------------------------------
## Read data sets  ----
SARD_dat <- read.csv(here("data", "SARD_Mar2021_14C.csv"))
HumActCA_dat <- read.csv(here("data", "HumActCA_Dec2020_14C.csv")) #Alternatively read in "aDRAC_April2023_14C.csv" -- this is a more up to date version from the aDRAC GitHub repository, but then we need to also filter for Bantu associated dates
EA_dat <- read.csv(here("data", "EastAfrica_April2016_14C.csv"))
#KayWA_dat <- read.csv(here("data", "KayWA_SUMMARY_May2019_14C.csv")) #TODO: Replace this summary data with actual radiocarbon data

p3k14c_dat <- p3k14c::p3k14c_data #TODO: replace with the raw data -- this will include periods

#-------------------------------------------------------------------------------
## Data filtering and cleaning ----

SARD_sum_df <- SARD_dat %>%
  filter(Archaeological.Period=="Iron Age") %>% #TODO: refine the determination of whether the site is Bantu at a later stage...
  select(X.Site, DecdegE, DecdegS, Date, Uncertainty) %>%
  rename(site=X.Site, lat=DecdegS, long=DecdegE, c14date=Date, c14std=Uncertainty) %>%
  mutate(c14date = as.numeric(c14date), c14std=as.numeric(c14std), dataorigin="SARD") 

HumActCA_sum_df <- HumActCA_dat %>% #TODO: filter for unreliable classes of dates, eg. III (and maybe II dates) #I think all the sites in this database are associated with Bantu pottery finds
  select(SITE, LAT, LONG, C14AGE, C14STD) %>%
  rename(site=SITE, lat=LAT, long=LONG, c14date=C14AGE, c14std=C14STD) %>%
  mutate(c14date = as.numeric(c14date), c14std=as.numeric(c14std), dataorigin="HumActCA")

EA_sum_df <- EA_dat %>%  #TODO: filter for Bantu dates
  select(Site.Identifier, Latitude, Longitude, Normalized.Age, NA.Sigma) %>%
  rename(site=Site.Identifier, lat=Latitude, long=Longitude, c14date=Normalized.Age, c14std=NA.Sigma) %>%
  mutate(c14date = as.numeric(c14date), c14std=as.numeric(c14std), dataorigin="EastAfrican") 

# KayWA_sum_df <- KayWA_dat %>%  #TODO: filter for Bantu dates
#   select(Site.Name, X, Y) %>%
#   rename(site=Site.Name, lat=X, long=Y) %>%
#   mutate(c14date = as.numeric(c14date), c14std=as.numeric(c14std), dataorigin="KayWestAfrican") 


p3k14c_sum_df <- p3k14c_dat %>%
  filter(Continent=="Africa" & 
           (Country %in% c(subSahara_countries, "CAR", "Equat.Guinea", "DRC"))) %>% #TODO: refine the determination of whether the site is Bantu at a later stage -- specifically filter for the correct period (waiting for raw data)
  select(SiteName, Lat, Long, Age, Error) %>%
  rename(site=SiteName, lat=Lat, long=Long, c14date=Age, c14std=Error) %>%
  mutate(c14date = as.numeric(c14date), c14std=as.numeric(c14std), dataorigin="p3k14c") %>%
  filter(!is.na(lat) | !is.na(long))


#Combine datasets
bantu_sites_df <- rbind(HumActCA_sum_df, 
                        SARD_sum_df, 
                        EA_sum_df) %>%  #TODO: Not all these sites are Bantu -- need more careful filtering... #TODO: Add p3k14c_sum_df once it's been replaced with a filtered Bantu version
  mutate(dataorigin=as.factor(dataorigin)) %>% 
  drop_na()


