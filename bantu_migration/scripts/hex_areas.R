# Load Libraries ----
library(here)
library(nimbleCarbon)
library(parallel)
library(coda)
library(rcarbon)
library(dplyr)
library(stringr)
library(rnaturalearth)
library(sp)
library(maptools)
library(sf)
library(ggplot2)
library(viridis)
library(rgeos)
library(raster)
library(gridExtra)
library(rasterVis)
library(rgbif)
library(rgdal)

#-------------------------------------------------------------------------------
# Load and prepare data ----
load(here('data','c14.RData'))

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
# Generate Spatial Window for Analyses: Sub-Saharan Africa ----

sf_subsah_africa <- ne_countries(continent = "Africa", returnclass = "sf") %>%
  filter_all(., any_vars(str_detect(., "Sub-Saharan"))) %>% 
  filter(name_en %in% subSahara_countries) %>% 
  filter(name_en != "Madagascar") #We focus on mainland sub-Saharan Africa

sp_ss_africa <- sf_subsah_africa %>% as("Spatial") #convert sf to sp object

sampling_win <- sp_ss_africa %>%
  disaggregate %>% #create new raster layer with higher resolution (smaller cells)
  geometry %>% 
  gUnaryUnion() #Remove internal country borders


#Generate Spatial Hexagons (without projection) ----

cell_diameter <- 5
buffer_win <- as(extent(sampling_win) + cell_diameter, "SpatialPolygons")

hex_points <- spsample(buffer_win, type = "hexagonal", cellsize = cell_diameter, offset = c(0.5, 0.5)) #offset ensures same grid of points are generated each time
hex_grid <- HexPoints2SpatialPolygons(hex_points, dx = cell_diameter)
hex_grid_clipped <- rgeos::gIntersection(hex_grid, sampling_win, byid = TRUE)

row.names(hex_grid_clipped) <- as.character(1:length(hex_grid_clipped))

#Plot ----
plot(sampling_win, col = "grey50", bg = "light blue", axes = TRUE, cex = 20) 
plot(hex_points, col = "black", pch = 20, cex = 0.5, add = T)
plot(hex_grid_clipped, border = "orange", add = T)









##Example (with projection, where cellsize is important)
# make_grid <- function(x, cell_diameter, cell_area, clip = FALSE) {
#   if (missing(cell_diameter)) {
#     if (missing(cell_area)) {
#       stop("Must provide cell_diameter or cell_area")
#     } else {
#       cell_diameter <- sqrt(2 * cell_area / sqrt(3))
#     }
#   }
#   ext <- as(extent(x) + cell_diameter, "SpatialPolygons")
#   projection(ext) <- projection(x)
#   # generate array of hexagon centers
#   g <- spsample(ext, type = "hexagonal", cellsize = cell_diameter, 
#                 offset = c(0.5, 0.5))
#   # convert center points to hexagons
#   g <- HexPoints2SpatialPolygons(g, dx = cell_diameter)
#   # clip to boundary of study area
#   if (clip) {
#     g <- gIntersection(g, x, byid = TRUE)
#   } else {
#     g <- g[x, ]
#   }
#   # clean up feature IDs
#   row.names(g) <- as.character(1:length(g))
#   return(g)
# }
# 
# study_area <- sampling_win
# plot(study_area, col = "grey50", bg = "light blue", axes = TRUE, cex = 20)
# 
# study_area_utm <- CRS("+proj=utm +zone=44 +datum=WGS84 +units=km +no_defs") %>% 
#   spTransform(study_area, .)
# # without clipping
# hex_grid <- make_grid(study_area_utm, cell_area = 200000, clip = FALSE)
# plot(study_area_utm, col = "grey50", bg = "light blue", axes = FALSE)
# plot(hex_grid_clipped, border = "orange", add = TRUE)
# box()
# # with clipping
# hex_grid <- make_grid(study_area_utm, cell_area = 200000, clip = TRUE)
# plot(study_area_utm, col = "grey50", bg = "light blue", axes = FALSE)
# plot(hex_grid, border = "orange", add = TRUE)
# box()







