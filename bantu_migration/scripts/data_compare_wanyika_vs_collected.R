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

##Collected dates 
#Compiled dates covering Eastern and Southern Africa

##Wanyika database
#Eastern africa radiocarbon dates from Kenya, Tanzania, the Comoros Islands, and Madagascar
#https://pandoradata.earth/dataset/wanyika

#-------------------------------------------------------------------------------
## Read data sets  ----
SARD_dat <- read.csv(here("data", "SARD_Mar2021_14C.csv"))
Collected_EIA_dat <- read_excel(here("data", "collected_dates.xlsx"), sheet = "Collected_dates_final") 
Wanyika_dat <- read_csv("data/wanyika_chronological_database.csv")

#-------------------------------------------------------------------------------
## Clean data ----

##Collected
collected_df <- Collected_EIA_dat %>%
  dplyr::select(LabID, Age, Error, Long, Lat, SiteName, Material, Country) %>%
  rename(labCode=LabID, siteName=SiteName, lat=Lat, long=Long, c14date=Age, c14std=Error, material=Material, country=Country) %>%
  mutate(c14date = as.numeric(c14date), c14std=as.numeric(c14std), dataorigin="Collected")  %>%
  filter(!is.na(lat) & !is.na(long) & !is.na(c14std) & !is.na(c14date)) 

##Wanyika
wanyika_df <- Wanyika_dat %>%
  dplyr::select("Labcode", "Date BP", "Date BP SD", "Longitude", "Latitude", "Site Name", "Dated Material", "Country") %>%
  rename(labCode="Labcode", siteName="Site Name", lat="Latitude", long="Longitude", c14date="Date BP", c14std="Date BP SD", material="Dated Material", country="Country") %>%
  mutate(c14date = as.numeric(c14date), c14std=as.numeric(c14std), dataorigin="Wanyika")  %>%
  filter(!is.na(lat) & !is.na(long) & !is.na(c14std) & !is.na(c14date)) 

##SARD
#Correct incorrect data entries
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


##Left join
#X_df <- merge(wanyika_df, collected_df, by="labCode")

## Combine datasets ----
bantu_sites_df <- bind_rows(SARD_df, wanyika_df, collected_df) %>% 
  mutate(dataorigin=as.factor(dataorigin)) %>% 
  filter((c14date != 0) & (c14std != 0)) %>% #Apparently some of these datasets had modern dates (indicated with c14 date and error of 0), we remove these
  filter((c14date <=7000) & (c14date >=246)) #Dates earlier than this are assumed to not be of Bantu origin #TODO: build some flexibility into this... #Further, we assume the Dutch arrival in the Cape (1652) as the cut-off. Dates after this point of colonial contact are not considered.

# Assign ID ----
bantu_sites_df <- bantu_sites_df %>% mutate(ID = row_number())

# Assign Site ID ----
bantu_sites_df$siteID  <- as.numeric(factor(bantu_sites_df$siteName))









## Data preparation ----

#Sampling window without internal boundaries
countries <- constants$countries #eastEIAcountries
#cntry_sampling_win <- ne_countries(continent = "Africa", country = countries, returnclass = "sf")
cntry_sampling_win <- ne_countries(country = countries, returnclass = "sf") %>%
  filter(name_en %!in% c("Madagascar","Sudan")) #We focus on mainland sub-Saharan Africa (also, there is something wrong with the geometry of Sudan -- remove country since we have no iron age dates there anyway)

EA_cntry_sampling_win <- ne_countries(country = countries, returnclass = "sf") %>%
  filter(name_en %in% constants$eastEIAcountries) 

#-------------------------------------------------------------------------------
## Plot Data  ---- FIGURE figure_map

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

minimap <- ggplot(data = world) +
  geom_sf(color = NA, fill = "grey") + 
  geom_rect(xmin = 7, xmax = 50, 
            ymin = -35, ymax = 6.5, 
            fill = NA, color = "black") + 
  coord_sf(xlim = c(-15, 50), 
           ylim = c(-35, 35)) + 
  theme_void() + 
  theme(panel.border = element_rect(colour = "darkgrey", 
                                    fill = NA, size = .5))

#Basemap taken from HumActCA project, White's vegetation descriptions. Assuming vegetation areas need to be extended to rest of sub-Saharan Africa...
basemap <- function(){
  white <- sf::st_read("input/Whites vegetation.shp") %>%
    st_set_crs(4326) %>%
    dplyr::filter(DESCRIPTIO %in% c("Anthropic landscapes",
                                    "Dry forest and thicket",
                                    "Swamp forest and mangrove",
                                    "Tropical lowland rainforest"))  
  
  # Vector layers ----
  rivers10 <- ne_download(scale = 10, type = "rivers_lake_centerlines", category = "physical", returnclass="sf")
  lakes10 <- ne_download(scale = 10, type = "lakes", category = "physical", returnclass="sf")
  coast10 <- ne_download(scale = 10, type = "coastline", category = "physical", returnclass="sf")
  land10 <- ne_download(scale = 10, type = "land", category = "physical", returnclass="sf")
  boundary_lines_land10 <- ne_download(scale = 10, type = "boundary_lines_land", category = "cultural", returnclass="sf")
  
  # Base map plot ----
  plt <- ggplot() + 
    geom_sf(data = white, fill = "grey", color = NA) + 
    geom_sf(data = coast10, size = .5, color = '#808080') + 
    geom_sf(data = rivers10, size = .5, color = '#808080') + 
    geom_sf(data = lakes10, fill = '#808080', color = NA) + 
    geom_sf(data = boundary_lines_land10, size = .1, color = 'black') 
  
  return(plt)
}

#sites <- sites %>% filter(dataorigin %in% c("aDRAC", "SARD"))

#Ploting sites with basemap ----
plt.main <- basemap() +
  geom_sf(data = sites,
          aes(colour=dataorigin),
          size = 2,
          alpha=0.5) +
  geom_point() +
  #geom_point(aes(x=constants$origin_point[1], y=constants$origin_point[2]), colour="purple", size=3) +
  ggsn::north(data = sites, location="bottomright", anchor = c(x = 43, y = -31)) + 
  ggsn::scalebar(sites,
                 location  = "bottomright",
                 anchor = c(x = 46, y = -33),
                 dist = 500, 
                 dist_unit = "km",
                 transform = TRUE, 
                 model = "WGS84",
                 height = .01, 
                 st.dist = .025,
                 border.size = .1, 
                 st.size = 3) +
  coord_sf(xlim = c(7, 50),
           ylim = c(-35, 6.5)) +
  scale_x_continuous(breaks = seq(8, 50, 2)) +
  labs(colour="Original dataset") +
  scale_colour_discrete(labels = c("aDRAC", "Collected", "SARD")) + #c(Collected", "SARD")
  theme_few() +
  theme(axis.title = element_blank(),
        plot.background = element_rect(color = NA,
                                       fill = NA))


pdf(file=here('output','figures','figure_map_cont.pdf'), width=8.5, height=7)
cowplot::ggdraw() +
  draw_plot(plt.main) +
  draw_plot(minimap, 
            x = .05, y = .275, width = .15, height = .15)
dev.off()


