hex_areas <- function(win, cell_d = 6){
  require(dplyr)
  require(stringr)
  require(rnaturalearth)
  require(sp)
  require(sf)
  require(raster)
  require(rgeos)
  
  #Spatial Window for Analyses ----
  sampling_win <- win %>%
    disaggregate %>% #create new raster layer with higher resolution (smaller cells)
    geometry %>% 
    gUnaryUnion() #Remove internal country borders
  
  
  #Generate Spatial Hexagons ----
  cell_diameter <- cell_d
  buffer_win <- as(extent(sampling_win) + cell_diameter*2, "SpatialPolygons")
  
  hex_points <- spsample(buffer_win, type = "hexagonal", cellsize = cell_diameter, offset = c(0.5, 0.5)) #offset ensures same grid of points are generated each time
  hex_grid <- HexPoints2SpatialPolygons(hex_points, dx = cell_diameter)

  #Projection ----
  proj4string(hex_points)  <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
  proj4string(hex_grid)  <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
  
  
  #Clip to sampling window with coastal buffer ---
  coastal_buffer_win <- st_buffer(as(sampling_win, 'sf'), 40000) #40km buffer to ensure all sites fall in sample_win
  hex_points_clipped <- st_intersection(as(hex_points,'sf'), coastal_buffer_win)
  hex_grid_clipped <- st_intersection(as(hex_grid,'sf'), coastal_buffer_win)
  
  #Assign hex IDs ----
  hex_grid_clipped <- hex_grid_clipped %>% 
    mutate(area_ID = row_number())  
  
  
  #Return Output ----
  return(hex_grid_clipped)
  
}



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