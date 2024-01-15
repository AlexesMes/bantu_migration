# Load Libraries and Data ----
library(maptools)
library(sf)
library(stringr)
library(dplyr)
library(here)
library(ggplot2)
library(ggthemes)
library(ggsn) #for scale bar and north arrow
library(rnaturalearth)
library(rnaturalearthdata)
library(parallel)
library(RColorBrewer)
library(cowplot)
library(rgeos)
#library(RCDT) #constrained delaunay tesselation
library(ggforce)
library(deldir)

# Load and prepare data ----
load(here('data','c14.RData'))


#-------------------------------------------------------------------------------
## Data preparation ----

#Convert to sf objects
bantu_sites_sf <- sf::st_as_sf(siteInfo, 
                               coords = c("long", "lat"), 
                               remove = F, 
                               crs = 4326, 
                               na.fail = F)


#-------------------------------------------------------------------------------
## Compute Great-Arc Distances in km between area centers ---
#We take the center of each area k to be a point representing that whole area
hex_area_centers <- hex_area_win$area_center %>% as('Spatial')
hex_dist_mat <- spDists(hex_area_centers, longlat=TRUE) #Inter-area distance matrix: each area's distance from every other area.


#Delaunay triangulation between hex centers
del <- deldir(hex_area_centers@coords, id=hex_area_win$area_ID)
tiles <- tile.list(del)



##Plot delaunay triangulation
ggplot(data = hex_area_win) +
  geom_sf(data = st_buffer(as(gUnaryUnion(sampling_win), 'sf'), 40000), aes(color = "grey50")) + #sampling window with coastal buffer
  geom_sf() + #hex grid
  #geom_sf_text(aes(label = area_ID), size=2, alpha=0.5) + #hex grid labels
  geom_sf(data = hex_area_win$area_center, size=2, alpha=1, aes(color = "purple")) + #hex-origins
  geom_delaunay_segment(aes(x=hex_area_centers@coords[,1], y=hex_area_centers@coords[,2]), 
                        alpha=0.5, 
                        colour='purple',
                        size=0.8) +
  labs(x = "Longitude", y = "Latitude") +
  theme(panel.background = element_rect(fill = "lightblue",
                                        colour = "lightblue",
                                        size = 0.5,
                                        linetype = "solid"),
        legend.position = "none")




