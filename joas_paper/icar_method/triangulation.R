# Load Libraries and Data ----
library(sf)
library(stringr)
library(dplyr)
library(here)
library(ggplot2)
library(ggthemes)
library(parallel)
library(RColorBrewer)
library(cowplot)
library(ggforce)
library(deldir)
library(units)
library(igraph)
library(rnaturalearth)

source(here('src','hex_areas.R'))

#-------------------------------------------------------------------------------
# Generate Spatial Window for Analyses: Box same size as Europe ----

#-------------------------------------------------------------------------------
##Sampling window: Europe ----
#List of countries in europe (excluding france and norway which need to be added separately if needed) ---
europe_countries <- c("Albania", "Latvia", "Andorra", "Liechtenstein",
                      "Armenia", "Lithuania", "Austria", "Luxembourg",
                      "Azerbaijan", "Malta", "Belarus", "Moldova",
                      "Belgium", "Monaco", "Bosnia and Herzegovina", "Montenegro",
                      "Bulgaria", "Netherlands", "Croatia", 
                      "Cyprus", "Poland", "Czechia", "Portugal",
                      "Denmark", "Romania", "Estonia","Finland", "San Marino",
                      "Macedonia", "Republic of Serbia", "Slovakia",
                      "Georgia", "Slovenia", "Germany", "Spain", "Greece", "Sweden",
                      "Hungary", "Switzerland", "Ireland", "Turkey",
                      "Italy", "Ukraine", "Kosovo", "United Kingdom")

sampling_europe_win <- ne_countries(country=europe_countries, returnclass = "sf")

#Generate Spatial Hexagons --
hex_areas <- function(win, cell_d = 6){
  
  #Spatial Window for Analyses ----
  sf::sf_use_s2(FALSE) #turn off spherical co-ordinates
  sampling_win <-  win %>%
    st_make_valid() %>%
    st_union() %>% 
    st_bbox() %>% 
    st_as_sfc(crs=4326)
  sf::sf_use_s2(TRUE) #turn on spherical co-ordinates
  
  #Generate Spatial Hexagons ----
  cell_diameter <- cell_d
  #hex_points <- st_make_grid(buffer_win, square=FALSE,  cellsize = cell_diameter, what = "centers")
  hex_grid <- st_make_grid(sampling_win, square=FALSE,  cellsize = cell_diameter) #makes an hexagonal grid (default: what = "polygons)
  
  #Projection ----
  st_crs(hex_grid)  <- 4326
  
  #Assign hex IDs ----
  hex_grid <- st_as_sf(hex_grid) %>%
    rename(geometry = x) %>%
    mutate(area_ID = row_number(),
           area_center = st_centroid(hex_grid))
  
  #Return Output ----
  return(hex_grid)
  
}  

hex_area_win <- hex_areas(sampling_europe_win, cell_d = 5)
sampling_win_box <- st_as_sfc(st_bbox(hex_area_win), crs=4326)

##Plot sample window and hex areal units
ggplot(data = hex_area_win) +
  geom_sf(data = st_as_sfc(st_bbox(sampling_win_box), crs=4326), fill = "green", alpha = 0.2, lwd=0) + #internal country borders
  geom_sf() + #hex grid
  geom_point(aes(x=-11,y=72.5), colour="purple", size=3) + #origin
  theme(panel.background = element_rect(fill = "lightblue",
                                        colour = "lightblue",
                                        size = 0.5,
                                        linetype = "solid"),
        legend.position = "none")

#-------------------------------------------------------------------------------
## Compute Great-Arc Distances in km between area centers ---
#We take the center of each area k to be a point representing that whole area
hex_area_centers <- st_as_sf(hex_area_win$area_center)
st_crs(hex_area_centers)  <- 4326
hex_dist_mat <- set_units(st_distance(hex_area_centers), 'km') #Inter-area distance matrix in km: each area's distance from every other area.


##Delaunay triangulation between hex centers
center_coords <- st_coordinates(hex_area_centers) #TODO: Uncomment
del <- deldir(center_coords, id=hex_area_win$area_ID)
tiles <- tile.list(del)

# Add center_coords in constants
constants_trig <- list()
constants_trig$center_coords <- center_coords

##Plot delaunay triangulation
# ggplot(data = hex_area_win) + 
#   geom_sf(data = sampling_win_box, fill = "green", alpha = 0.2, lwd=0) + #sampling window
#   geom_sf() + #hex grid
#   geom_sf_text(aes(label = area_ID), size=4, alpha=0.8) + #hex grid labels #aes(label =  c('1','2','3'))
#   geom_sf(data = hex_area_win$area_center, size=2, alpha=1, aes(color = "purple")) + #hex-origins
#   geom_delaunay_segment(aes(x=center_coords[,1], y=center_coords[,2]),
#                         alpha=0.5,
#                         colour='purple',
#                         size=0.8) +
#   labs(x = "Longitude", y = "Latitude") +
#   theme(panel.background = element_rect(fill = "lightblue",
#                                         colour = "lightblue",
#                                         size = 0.5,
#                                         linetype = "solid"),
#         legend.position = "none")

#-------------------------------------------------------------------------------
##Transitions dataframe --
transitions <- del$delsgs %>% 
  dplyr::select(region1_id = ind1, 
                region2_id = ind2,
                region1_x = x1,
                region1_y = y1,
                region2_x = x2,
                region2_y = y2) %>% 
  mutate(region1_id = as.integer(region1_id), 
         region2_id = as.integer(region2_id)) %>% 
  rowwise() %>% 
  mutate(distance = hex_dist_mat[region1_id, region2_id]) #Great-arc distance between transitions in km

#Transform transitions into usable format to save in constants
edge_info <- as.data.frame(transitions)
constants_trig$n_trans <- nrow(transitions) #Number of transitions
constants_trig$edge_id1 <- edge_info$region1_id
constants_trig$edge_id2 <- edge_info$region2_id 
constants_trig$edge_dist <- edge_info$distance


#-------------------------------------------------------------------------------
## Save transition data on a R image file ----
save(del, tiles, edge_info, constants_trig, file=here('data','trig.RData'))

## Save sampling window specific information separately ----
save(sampling_win_box, hex_area_win, file=here('data','sample_window.RData'))

