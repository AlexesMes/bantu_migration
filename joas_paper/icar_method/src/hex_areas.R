hex_areas <- function(win, cell_d = 6){
  require(dplyr)
  require(stringr)
  require(rnaturalearth)
  require(sf)
  
  #Spatial Window for Analyses ----
  sf::sf_use_s2(FALSE) #turn off spherical co-ordinates
  sampling_win <-  win %>%
    st_make_valid() %>%
    st_union()
  sf::sf_use_s2(TRUE) #turn on spherical co-ordinates
  
  #Generate Spatial Hexagons ----
  cell_diameter <- cell_d
  buffer_win <- st_as_sfc(st_bbox(sampling_win), crs = 4326)
  
  #hex_points <- st_make_grid(buffer_win, square=FALSE,  cellsize = cell_diameter, what = "centers")
  hex_grid <- st_make_grid(buffer_win, square=FALSE,  cellsize = cell_diameter) #makes an hexagonal grid (default: what = "polygons)
  
  #Projection ----
  st_crs(hex_grid)  <- 4326
  
  #Clip to sampling window with coastal buffer ---
  coastal_buffer_win <- st_buffer(st_as_sf(sampling_win, crs = 4326), 40000) #40km buffer to ensure all sites fall in sample_win
  hex_grid_clipped <- st_intersection(st_as_sf(hex_grid), coastal_buffer_win)
  
  #Assign hex IDs ----
  hex_grid_clipped <- hex_grid_clipped %>%
    rename(geometry = x) %>%
    mutate(area_ID = row_number(),
           area_center = st_centroid(hex_grid_clipped))
  
  #Return Output ----
  return(hex_grid_clipped)
  
}  