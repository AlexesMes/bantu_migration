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
##List of countries in europe (excluding france and norway which need to be added separately) ---
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
#-------------------------------------------------------------------------------
##Generate Spatial Window for Analyses
#Sampling window: Europe ----

sf::sf_use_s2(FALSE)
sampling_europe_win <- st_union(st_make_valid(ne_countries(country=europe_countries, returnclass = "sf")),
                                st_make_valid(ne_countries(geounit = c("france", "norway"), type = "map_units", returnclass = "sf"))) # France filter map_units by geounit to exclude French Guiana. Similar problem for Norway
sf::sf_use_s2(TRUE)

sampling_win <- st_as_sfc(st_bbox(sampling_europe_win, crs = 4326))

#Generate Spatial Hexagons --
hex_areas <- function(win, cell_d = 6){
  cell_diameter <- cell_d
  hex_grid <- st_make_grid(win, square=FALSE,  cellsize = cell_diameter) #makes an hexagonal grid (default: what = "polygons)
  
  #Projection ----
  st_crs(hex_grid)  <- 4326
  
  #Clip to sampling window with coastal buffer ---
  #coastal_buffer_win <- st_as_sf(win, crs = 4326)
  #hex_grid_clipped <- st_intersection(st_as_sf(hex_grid), coastal_buffer_win)
  
  #Assign hex IDs ----
  hex_grid <- st_as_sf(hex_grid) %>%
    rename(geometry = x) %>%
    mutate(area_ID = row_number(),
           area_center = st_centroid(hex_grid))
  
  #Return Output ----
  return(hex_grid)
}  

hex_area_win <- hex_areas(sampling_win, cell_d = 5) 

##Plot sample window and hex areal units
ggplot(data = hex_area_win) +
  geom_sf(data = sampling_win, aes(color = "grey50"), lwd=1) + #internal country borders
  geom_sf() + #hex grid
  theme(panel.background = element_rect(fill = "lightblue",
                                        colour = "lightblue",
                                        size = 0.5,
                                        linetype = "solid"),
        legend.position = "none")























#Assign hex area id to each site ----
siteInfo$area_id <- as.integer(st_within(sites$geometry, hex_area_win$geometry))
#Assign hex area id to each date ----
dateInfo$area_id <- siteInfo$area_id[match(dateInfo$siteID, siteInfo$siteID)]

# #CHECK ---
area_freq  <- plyr::count(siteInfo, 'area_id') ##See how many sites fall in each hex area. Also make sure there are no 'NA' entries
#In order to check that this lines up visually with how many sites are in each hex area see map_figure2

#--------------------------------
# ## Determining hex size ---
# #Under changing hex size, determine the proportion of areal hex units in the sampling window with sites
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
# pdf(here('output','figures','figure_hexsize_cont.pdf'),height=5,width=5.5)
# ggplot(prop_units_df, aes(x = d, y = prop_with_sites)) +
#   geom_line() +
#   geom_point() +
#   scale_x_continuous(breaks=seq(0,15,by=1))+
#   labs(x = "Hexagon Size (d)", y = "Proportion of Hexagons with Sites", title = "Effect of Hexagon Size on Site Coverage") +
#   theme_minimal()
# dev.off()







#Convert to sf objects
eastEIA_sites_sf <- sf::st_as_sf(siteInfo, 
                                 coords = c("long", "lat"), 
                                 remove = F, 
                                 crs = 4326, 
                                 na.fail = F)

#Sampling window
sampling_win <- st_as_sf(sampling_win, crs = 4326)
#Sampling window without internal boundaries (and added coastal buffer)
sf::sf_use_s2(FALSE) #turn off spherical co-ordinates
sampling_win_ext <-  sampling_win %>%
  st_make_valid() %>%
  st_union() 
sf::sf_use_s2(TRUE) #turn on spherical co-ordinates

sampling_win_ext <- st_buffer(st_as_sf(sampling_win_ext, crs = 4326), 40000)

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


##Remove external edges that are outside sampling window ---

#EITHER...
# #Remove boundary transitions
# average_trans <- mean(transitions$distance)
# transitions <- transitions %>% filter(distance < average_trans)

#OR....

#Create edges as linestrings
st_segment <- function(r){st_linestring(t(matrix(unlist(r), 2, 2)))}
edges_coords <- del$delsgs[ , 1:4]
edges_coords$geom <- st_sfc(sapply(1:nrow(edges_coords), 
                                   function(i){st_segment(edges_coords[i,])}, simplify=FALSE))
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

# ##Plot delaunay triangulation
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
                region2_y = y2) %>% 
  mutate(region1_id = as.integer(region1_id), 
         region2_id = as.integer(region2_id)) %>% 
  rowwise() %>% 
  mutate(distance = hex_dist_mat[region1_id, region2_id]) #Great-arc distance between transitions in km

transitions <- transitions[-ind_ext_edges,] #remove transitions outside the sampling window

#Transform transitions into usable format to save in constants
edge_info <- as.data.frame(transitions)
constants_trig$n_trans <- nrow(transitions) #Number of transitions
constants_trig$edge_id1 <- edge_info$region1_id
constants_trig$edge_id2 <- edge_info$region2_id 
constants_trig$edge_dist <- edge_info$distance



#-------------------------------------------------------------------------------
## Save everything on a R image file ----
save(del, 
     tiles,
     edge_info,
     constants_trig,
     file=here('data','trig_cont.RData')) #'trig.RData'

save(edge_info,
     constants_trig,
     file=here('data','boundary_edges.RData'))