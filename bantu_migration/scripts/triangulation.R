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
#library(RCDT) #constrained delaunay tesselation
library(ggforce)
library(deldir)
library(units)
library(igraph)

# Load and prepare data ----
load(here('data','eastc14.RData'))


#-------------------------------------------------------------------------------
## Data preparation ----

#Convert to sf objects
eastEIA_sites_sf <- sf::st_as_sf(siteInfo, 
                               coords = c("long", "lat"), 
                               remove = F, 
                               crs = 4326, 
                               na.fail = F)

#Sampling window
sampling_win <- st_as_sf(sampling_win, crs = 4326)
#Sampling window without internal boundaries
sf::sf_use_s2(FALSE) #turn off spherical co-ordinates
sampling_win_ext <-  sampling_win %>%
  st_make_valid() %>%
  st_union() 
sf::sf_use_s2(TRUE) #turn on spherical co-ordinates

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
##Selecting 3 triangles of interest -- TEST CASE
# center_coords <- st_coordinates(hex_area_centers[c(13, 18, 22), ])
# del <- deldir(center_coords, id=c(1, 2, 3)) #Relabel ids for nimble (must be sequential from 1) 13 -> 1, 14 -> 2, 18 -> 3
# tiles <- tile.list(del)

##Remove external edges that are outside sampling window
#Create edges as linestrings
st_segment <- function(r){st_linestring(t(matrix(unlist(r), 2, 2)))}
edges_coords <- del$delsgs[ , 1:4]
edges_coords$geom <- st_sfc(sapply(1:nrow(edges_coords), 
                           function(i){st_segment(edges_coords[i,])},simplify=FALSE))
st_crs(edges_coords$geom)  <- 4326
#Check if edges are external
check_edges_ext <- st_within(edges_coords$geom, sampling_win_ext)
ind_ext_edges <- which(is.na(as.numeric(check_edges_ext))) # indices of edges which are external -- to be removed if we aren't considering ocean movement
##Check
#internal_edges <- edges_coords[-ind_ext_edges, ]
#plot(sampling_win_ext)
#plot(internal_edges$geom, add=T) #compare to plot(edges_coords$geom, add=T)

# Add center_coords in constants
constants_trig <- list()
constants_trig$center_coords <- center_coords

##Plot delaunay triangulation
# #ggplot(data = hex_area_win[c(13, 18, 22),]) +
# ggplot(data = hex_area_win) + #TODO: Uncomment
#   geom_sf(data = st_buffer(sampling_win, 40000), aes(color = "grey50")) + #sampling window with coastal buffer
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
                region2_y = y1) %>% 
  mutate(region1_id = as.integer(region1_id), 
         region2_id = as.integer(region2_id)) %>% 
  rowwise() %>% 
  mutate(distance = hex_dist_mat[region1_id, region2_id]) #Great-arc distance between transitions in km


# Add transitions as matrix in constants
constants_trig$transitions  <- as.matrix(transitions)
constants_trig$n_trans <- nrow(transitions) #Number of transitions

#-------------------------------------------------------------------------------
## Save everything on a R image file ----
save(del, 
     tiles, 
     constants_trig,
     file=here('data','trig.RData'))


#===============================================================================
##Simple paths

relations <- transitions %>% 
  dplyr::select(from = region1_id, to = region2_id) 


vertices <- st_drop_geometry(hex_area_win) 

#vertices <- vertices[c(13, 18, 22),] %>% 
#  mutate(area_ID = case_when(area_ID == 13 ~ 1, area_ID == 18 ~ 2, area_ID == 22 ~ 3))  #test_case

hex_centers_graph <- graph_from_data_frame(relations, directed=FALSE, vertices = NULL)
plot(hex_centers_graph)

create_paths <- all_simple_paths(hex_centers_graph, from = 25, mode = "out") #area_ID = 25 contains 'Katuruka' our putative origin


#-------------------------------------------------------------------------------
## Save everything on a R image file ----
save(create_paths, file=here('data','simplepaths.RData'))






