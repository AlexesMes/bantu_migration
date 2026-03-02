# Load Libraries and Data ----
library(sf)
library(stringr)
library(dplyr)
library(here)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(grid)
library(gridBase)
library(ggspatial)
library(rnaturalearth)
library(rnaturalearthdata)
library(parallel)
library(RColorBrewer)
library(cowplot)
library(viridis)

# Load and prepare data ----
load(here('data','eastc14.RData'))

`%!in%` <- Negate(`%in%`)
#-------------------------------------------------------------------------------
## Data preparation ----

#Sampling window without internal boundaries
countries <- constants$countries
cntry_sampling_win <- ne_countries(country = countries, returnclass = "sf") %>%
                            filter(name_en %!in% c("Madagascar","Sudan")) #We focus on mainland sub-Saharan Africa (also, there is something wrong with the geometry of Sudan -- remove country since we have no iron age dates there anyway)

EA_cntry_sampling_win <- ne_countries(country = countries, returnclass = "sf") %>%
  filter(name_en %in% constants$eastEIAcountries) 

#-------------------------------------------------------------------------------
##FIGURE -- study region and sites used in the analysis

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

#Plotting sites with basemap ----
plt.main <- basemap() +
  geom_sf(data = sites,
          aes(colour=dataorigin),
          size = 2,
          alpha=0.6) +
  scale_colour_manual(
    values = c("Collected" = "#E69F00","SARD" = "#0072B2","Wanyika" = "#009E73" ),
    name = "Data Origin",
    labels = c(
      "Collected" = "Newly digitalised",
      "SARD"      = "SARD",
      "Wanyika"   = "Wanyika")) +
  geom_sf(data = st_buffer(st_as_sf(st_union(EA_cntry_sampling_win), crs = 4326), 60000), color = "grey30" , alpha=0, lwd=1.5) + #demarcate study area
  annotation_scale(
    location = "br",
    width_hint = 0.3,
    text_cex = 0.8,
    pad_x = unit(0.5, "cm"),
    pad_y = unit(0.5, "cm")
  ) +
  annotation_north_arrow(
    location = "br",
    which_north = "true",
    style = north_arrow_fancy_orienteering,
    pad_x = unit(1.5, "cm"),
    pad_y = unit(1, "cm")
  ) +
  coord_sf(xlim = c(7, 50),
           ylim = c(-35, 6.5)) +
  scale_x_continuous(breaks = seq(8, 50, 4)) +
  scale_y_continuous(breaks = seq(-35, 6.5, 5)) +
  theme_few() +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.text  = element_text(size = 11), 
        legend.title = element_text(size = 12), 
        panel.grid.major = element_line(color = "grey95", linewidth = 0.2))


pdf(file=here('output','figures','figure_sites_map.pdf'), width=8.5, height=7)
cowplot::ggdraw() +
  draw_plot(plt.main) +
  draw_plot(minimap, 
            x = .05, y = .275, width = .15, height = .15)
dev.off()


#-------------------------------------------------------------------------------
##SUPPLEMENTARY FIGURE -- Labeled hex areas and country borders

# Load Africa country boundaries as sf
africa <- ne_countries(scale = "medium", continent = "Africa", returnclass = "sf")

#Plot
hex_area_plot <- ggplot(data = hex_area_win) +
  geom_sf(data = st_union(africa), fill = "#ECE6DD", color = "black") +
  geom_sf() + #hex grid
  geom_sf(data = st_buffer(st_as_sf(st_union(sampling_win), crs = 4326), 60000), color = "grey30", alpha=0, lwd=1.5) + #sampling window with coastal buffer
  geom_sf(data = as(sites, 'sf'), size=2, alpha=0.5) + #sites
  geom_sf_label(aes(label = area_ID, alpha=0.7)) + #hex grid labels
  coord_sf(xlim = c(7, 45),
           ylim = c(-35, 6.5)) +
  scale_x_continuous(breaks = seq(8, 45, 4)) +
  scale_y_continuous(breaks = seq(-35, 6.5, 5)) +
  xlab("Longitude") + 
  ylab("Latitude") +
  theme(panel.background = element_rect(fill = "lightblue",
                                        colour = NA,
                                        size = 0.5,
                                        linetype = "solid"),
        panel.border = element_rect(color = "grey50", fill = NA, size = 0.5),
        legend.position = "none")

pdf(file=here('output','figures_supplementary','figure_map_hex_areas.pdf'), width=7, height=7)
hex_area_plot
dev.off()

# cntry_plot <- ggplot(data = hex_area_win) +
#   geom_sf(data = st_buffer(st_as_sf(EA_cntry_sampling_win, crs = 4326), 40000), aes(color = "grey50"), lwd=2) + #internal country borders
#   geom_sf(data = as(sites, 'sf'), size=2, alpha=0.5) + #sites
#   geom_sf_label(data = EA_cntry_sampling_win, aes(label = admin, alpha=0.6), color="darkred", size=4) + #country labels
#   theme(panel.background = element_rect(fill = "lightblue",
#                                         colour = "lightblue",
#                                         size = 0.5,
#                                         linetype = "solid"),
#         axis.text.x = element_text(size = 10),
#         axis.text.y = element_text(size = 10),
#         legend.position = "none")

##Both graphs overlaid -- helpful as a visual aid
# ggplot(data = hex_area_win) +
#   geom_sf(data = st_buffer(st_as_sf(cntry_sampling_win, crs = 4326), 40000), aes(color = "grey50"), lwd=2) + #internal country borders
#   geom_sf(data = as(sites, 'sf'), size=2, alpha=0.5) + #sites
#   geom_sf(aes(alpha=0.01)) + #hex grid
#   geom_sf_label(aes(label = area_ID)) + #hex grid labels
#   geom_sf_text(data = cntry_sampling_win, aes(label = admin), color="darkred", size=4) + #country labels
#   theme(panel.background = element_rect(fill = "lightblue",
#                                         colour = "lightblue",
#                                         size = 0.5,
#                                         linetype = "solid"),
#         legend.position = "none")

#-------------------------------------------------------------------------------
##SUPPLEMENTARY FIGURE -- Labeled hex areas and country borders with regions 23 and neighbours highlighted

# Load Africa country boundaries as sf
africa <- ne_countries(scale = "medium", continent = "Africa", returnclass = "sf")

#Select Area 23 and neighbours specifically (to examine likelihood of leapfrog transmission)
interest_areas <- c(33,24,14,28,18,9,13,5,17,8,12,23)
hex_area_win$neighbours <- hex_area_win$area_ID %in% interest_areas

#Plot
hex_area23_plot <- ggplot(data = hex_area_win) +
  geom_sf(data = st_union(africa), fill = "#ECE6DD", color = "black") +
  geom_sf(aes(fill = neighbours), color = "grey40") + #hex grid
  geom_sf_label(aes(label = area_ID, alpha=0.7)) + #hex grid labels
  coord_sf(xlim = c(7, 45),
           ylim = c(-35, 6.5)) +
  scale_x_continuous(breaks = seq(8, 45, 4)) +
  scale_y_continuous(breaks = seq(-35, 6.5, 5)) +
  scale_fill_manual(values = c(`TRUE` = "#7CCD7C", `FALSE` = "grey80"),guide = "none") +
  xlab("Longitude") + 
  ylab("Latitude") +
  theme(panel.background = element_rect(fill = "lightblue",
                                        colour = NA,
                                        size = 0.5,
                                        linetype = "solid"),
        panel.border = element_rect(color = "grey50", fill = NA, size = 0.5),
        legend.position = "none")

pdf(file=here('output','figures_supplementary','figure_map_hex_areas23.pdf'), width=7, height=7)
hex_area23_plot
dev.off()