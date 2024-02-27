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
proj4string(hex_area_centers)  <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
hex_dist_mat <- spDists(hex_area_centers, longlat=TRUE) #Inter-area distance matrix: each area's distance from every other area.


##Delaunay triangulation between hex centers
#center_coords <- hex_area_centers@coords #TODO: Uncomment
#del <- deldir(center_coords, id=hex_area_win$area_ID) 
#tiles <- tile.list(del)
##Selecting 3 triangles of interest -- TEST CASE
center_coords <- hex_area_centers@coords[c(13, 14, 18), ]
del <- deldir(center_coords, id=c(1, 2, 3)) #Relabel ids for nimble (must be sequential from 1) 13 -> 1, 14 -> 2, 18 -> 3
tiles <- tile.list(del)

# Add center_coords in constants
constants$center_coords <- center_coords

##Plot delaunay triangulation
ggplot(data = hex_area_win[c(13,14, 18),]) + #ggplot(data = hex_area_win) + #TODO: Uncomment
  geom_sf(data = st_buffer(as(gUnaryUnion(sampling_win), 'sf'), 40000), aes(color = "grey50")) + #sampling window with coastal buffer
  geom_sf() + #hex grid
  geom_sf_text(aes(label = c('1','2','3')), size=4, alpha=0.8) + #hex grid labels #aes(label = area_ID)
  geom_sf(data = hex_area_win$area_center, size=2, alpha=1, aes(color = "purple")) + #hex-origins
  geom_delaunay_segment(aes(x=center_coords[,1], y=center_coords[,2]), 
                        alpha=0.5, 
                        colour='purple',
                        size=0.8) +
  labs(x = "Longitude", y = "Latitude") +
  theme(panel.background = element_rect(fill = "lightblue",
                                        colour = "lightblue",
                                        size = 0.5,
                                        linetype = "solid"),
        legend.position = "none")

#-------------------------------------------------------------------------------
##Transitions dataframe --
transitions <- del$delsgs %>% 
  dplyr::select(region1_id = ind1, 
                region2_id = ind2,
                region1_x = x1,
                region1_y = y1,
                region2_x = x2,
                region2_y = y1) %>% 
  mutate(region1_id = as.integer(region1_id), 
         region2_id = as.integer(region2_id)) %>% 
  rowwise() %>% 
  mutate(distance = hex_dist_mat[region1_id, region2_id]) #Great-arc distance between transitions in km


# Add transitions as matrix in constants
constants$transitions  <- as.matrix(transitions)
constants$n_trans <- nrow(transitions) #Number of transitions

#-------------------------------------------------------------------------------
## Save everything on a R image file ----
save(del, 
     tiles, 
     constants,
     file=here('data','trig.RData'))


