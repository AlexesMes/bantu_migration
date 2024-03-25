# Load Libraries and Data ----
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

# Load and prepare data ----
load(here('data','eastc14.RData'))


#-------------------------------------------------------------------------------
## Data preparation ----

#Convert to sf objects
bantu_sites_sf <- sf::st_as_sf(siteInfo, 
                               coords = c("long", "lat"), 
                               remove = F, 
                               crs = 4326, 
                               na.fail = F)

#Sampling window without internal boundaries
sf::sf_use_s2(FALSE) #turn off spherical co-ordinates
sampling_win <-  sampling_win %>%
  st_make_valid() %>%
  st_union()
sf::sf_use_s2(TRUE) #turn on spherical co-ordinates

#-------------------------------------------------------------------------------
## Plot Data  ---- FIGURE map_figure

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

#Ploting sites with basemap ----
plt.main <- basemap() +
  geom_sf(data = bantu_sites_sf,
          aes(colour=dataorigin),
          size = 2,
          alpha=0.5) +
  geom_point() +
  geom_point(aes(x=st_coordinates(constants$origin_point)[1], y=st_coordinates(constants$origin_point)[2]), colour="purple", size=3) +
  ggsn::north(data = bantu_sites_sf, location="bottomright", anchor = c(x = 43, y = -31)) + 
  ggsn::scalebar(bantu_sites_sf,
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
  scale_colour_discrete(labels = c("aDRAC", "Collected", "SARD")) +
  theme_few() +
  theme(axis.title = element_blank(),
        plot.background = element_rect(color = NA,
                                       fill = NA))


pdf(file=here('output','figures','map_figure_EAST.pdf'), width=8.5, height=7)
cowplot::ggdraw() +
  draw_plot(plt.main) +
  draw_plot(minimap, 
            x = .05, y = .275, width = .15, height = .15)
dev.off()


#-------------------------------------------------------------------------------
## Plot HEX areas  ---- FIGURE map_figure2

pdf(file=here('output','figures','map_figure2_EAST.pdf'), width=8.5, height=7)

ggplot(data = hex_area_win) +
  geom_sf(data = st_buffer(st_as_sf(sampling_win, crs = 4326), 40000), aes(color = "grey50")) + #sampling window with coastal buffer
  geom_sf() + #hex grid
  geom_sf(data = as(sites, 'sf'), size=2, alpha=0.5) + #sites
  geom_sf_label(aes(label = area_ID)) + #hex grid labels
  #geom_sf(data = hex_area_win$area_center, size=2, alpha=1, aes(color = "purple")) + #hex-origins
  theme(panel.background = element_rect(fill = "lightblue",
                                        colour = "lightblue",
                                        size = 0.5,
                                        linetype = "solid"),
        legend.position = "none")

dev.off()


#Plotting country borders
cntry_sampling_win <- ne_countries(continent = "Africa", country = eastEIA_countries, returnclass = "sf")

ggplot(data = hex_area_win) +
  geom_sf(data = st_buffer(st_as_sf(cntry_sampling_win, crs = 4326), 40000), aes(color = "grey50"), lwd=2) + #internal country borders
  geom_sf(data = as(sites, 'sf'), size=2, alpha=0.5) + #sites
  geom_sf(aes(alpha=0.01)) + #hex grid
  geom_sf_label(aes(label = area_ID)) + #hex grid labels
  geom_sf_text(data = cntry_sampling_win, aes(label = admin), color="darkred", size=4) + #country labels
  theme(panel.background = element_rect(fill = "lightblue",
                                        colour = "lightblue",
                                        size = 0.5,
                                        linetype = "solid"),
        legend.position = "none")