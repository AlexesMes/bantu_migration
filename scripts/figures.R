# Load Libraries and spatial data ----
library(here)
library(dplyr)
library(tidyr)
library(stringr)
library(truncnorm)
library(cascsim)
library(corrplot)
library(ggplot2)
library(ggridges)
library(rnaturalearth)
library(nimbleCarbon)
library(rcarbon)
library(sf)
library(viridis)
library(cowplot)
library(wesanderson)
library(latex2exp)
library(gridExtra)
library(grid)
library(gridBase)
library(diagram)
library(quantreg)
library(coda)
library(graphics)
library(ggthemes)
library(gganimate)
library(gifski) 
library(ggrepel)
library(patchwork)
library(geodata)
library(colorspace)
library(RColorBrewer)
library(ggspatial)
library(rnaturalearthdata)
library(parallel)

source(here('src', 'grad_funcs.R'))
source(here('src', 'accuracy_precision.R'))
source(here('src','hex_areas.R'))
source(here('src','plotting_funcs.R'))

`%!in%` <- Negate(`%in%`)

#===============================================================================
## Load Data

#Load Observed Data
load(here('data','eastc14.RData'))
#Load nodes and edges between hex area centroids
load(here('data','trig_d38.RData'))
#Load environmental data
load(here('data','Krapp_enviro_variables.RData'))

#===============================================================================
## Set up
#Combine constants
constants <- c(constants, constants_trig)

#Africa country boundaries
africa <- ne_countries(scale = "medium", continent = "Africa", returnclass = "sf")

#Topographical features to map
lakes <- ne_download(scale = 10, type = "lakes", category = "physical", returnclass = "sf")
rivers <- ne_download(scale = 10,type = "rivers_lake_centerlines",category = "physical",returnclass = "sf")

#Obtain country codes
country_codes <- as.data.frame(africa) %>% dplyr::select(name, iso_a3) %>% filter(name %in% constants$eastEIAcountries)

#Import elevation data
SRTM90m <- elevation_30s(country_codes$iso_a3[1], path=here('input'), mask=TRUE)
for (i in 2:nrow(country_codes)){
  SRTM90m <- merge(SRTM90m, elevation_30s(country_codes$iso_a3[i], path=here('input'), mask=TRUE))
}

#Convert DEM to data frame 
SRTM_small <- aggregate(SRTM90m, fact = 5) #aggregate to ~5 km resolution
dem_df <- as.data.frame(SRTM_small, xy = TRUE)
colnames(dem_df) <- c("lon", "lat", "elevation")

#Plotting colours
barcolours1 = c("skyblue","dodgerblue","darkblue","darkgreen")
barcolours2 = c("plum","orchid","purple","darkgreen")

#-------------------------------------------------------------------------------
##FIGURE 1 -- study region and sites used in the analysis

## Data preparation ----

#Sampling window without internal boundaries
countries <- constants$countries
cntry_sampling_win <- ne_countries(country = countries, returnclass = "sf") %>%
  filter(name_en %!in% c("Madagascar","Sudan")) #We focus on mainland sub-Saharan Africa (also, there is something wrong with the geometry of Sudan -- remove country since we have no iron age dates there anyway)

EA_cntry_sampling_win <- ne_countries(country = countries, returnclass = "sf") %>%
  filter(name_en %in% constants$eastEIAcountries) %>% 
  mutate(lon = st_coordinates(st_centroid(geometry))[,1],
         lat = st_coordinates(st_centroid(geometry))[,2]) %>% #country_id assign spatially 
  arrange(desc(lat), lon) %>% 
  mutate(country_id = row_number())

#-------
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
    geom_sf(data = white, fill = "#C7BFBF", color = NA) + 
    geom_sf(data = coast10, size = .5, color = '#808080') + 
    geom_sf(data = rivers10, size = .5, color = '#8691BD') + 
    geom_sf(data = lakes10, fill = '#8691BD', color = NA) + 
    geom_sf(data = boundary_lines_land10, size = .5, color = 'black') 
  return(plt)
}

#Labels ----
countries_labels <- EA_cntry_sampling_win %>%
  st_point_on_surface() %>%
  select(country_id, name)

legend_df <- EA_cntry_sampling_win %>%
  st_drop_geometry() %>%
  select(country_id, name) %>%
  arrange(country_id)


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
  geom_sf_text(data = countries_labels,
               aes(label = country_id),
               size = 6,
               color = "#6E2727") +
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
  theme(legend.justification = c("right", "top"),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.text  = element_text(size = 11), 
        legend.title = element_text(size = 12), 
        panel.grid.major = element_line(color = "grey95", linewidth = 0.2))


legend_plot <- ggplot(legend_df, aes(x = 1, y = -country_id)) +
  geom_text(aes(label = paste0(country_id, ": ", name)),
            hjust = 0, size = 3) +
  theme_void() +
  xlim(1, 5)

pdf(file=here('output','figures','figure_sites_map.pdf'), width=8.5, height=7)
cowplot::ggdraw() +
  draw_plot(plt.main) +
  draw_plot(legend_plot, x = 0.8, y = 0.1, width = 0.2, height = 0.4) +
  draw_plot(minimap, 
            x = .05, y = .275, width = .15, height = .15)
dev.off()

#===============================================================================
##WOMBLE MODEL WITH ENVIRONEMENTAL COVARIATES ---

##Load Data ----
load(here("output","Womble_Emodel_d38.RData"))

#Extract
out.comb.womble.model  <- do.call(rbind, out_womble_model)
post.model.womble <- out.comb.womble.model[,paste0('a[',1:47,']')]  %>% round()  

#Hex areas with and without out sites
Hex_with_sites <- unique(siteInfo$area_id)
Hex_without_sites <- which(rep(1:47) %!in% Hex_with_sites)

#-----
#Calculate credible interval
ci_95 = credible_interval(out_womble_model, 0.95)

#Posterior summaries HDPI
low.model.womble  <- apply(post.model.womble, 2, quantile, probs = 0.025)
high.model.womble <- apply(post.model.womble, 2, quantile, probs = 0.975)
med.model.womble  <- apply(post.model.womble, 2, median)

arrival_df <- hex_area_win %>% 
  filter(area_ID %in% 1:47) %>% 
  mutate(median_date = med.model.womble,
         lower_ci    = low.model.womble,
         upper_ci    = high.model.womble,
         contains_sites = as.factor(case_when(area_ID %in% Hex_with_sites ~ 1, area_ID %in% Hex_without_sites ~ 0))) 


#-------------------------------------------------------------------------------
##FIGURE 2 -- Arrival by time panels 

time_slices <-  seq(2600, 1400, -150)
time_slices_BCAD <- BPtoBCAD(time_slices)
time_slices_BCAD_ch <- ifelse(
  time_slices_BCAD < 0,
  paste0(abs(time_slices_BCAD), " BCE"),
  paste0(time_slices_BCAD, " CE")
)


#Extract proportion of MCMC samples occurring before the specified time threshold
prop_model_womble  <- lapply(time_slices,
                             function(t) data.frame(x = 1:47,
                                                    y = sapply(as.data.frame(post.model.womble), 
                                                               prop_gthan_threshold, 
                                                               threshold = t)))
#Plot
plot_list <- list() 

for (k in 1:length(time_slices)) #all time slices
{
  #Save in data structure
  prop_threshold_womble <- hex_area_win %>%
    filter(area_ID %in% 1:47) %>%
    mutate(median_date = med.model.womble,
           prop_threshold = prop_model_womble[[k]]$y) 
  #Plot
  p <- ggplot(data = prop_threshold_womble) +
    geom_sf(data = st_union(africa), fill = "grey70", color = "black") +
    geom_sf(data = st_buffer(st_as_sf(st_union(sampling_win), crs = 4326), 2000), color = "grey30",fill = "white", lwd=0.3) + #sampling window with coastal buffer
    geom_sf(aes(fill = prop_threshold),color=NA) + #hex grid 
    coord_sf(xlim = c(7, 45),
             ylim = c(-35, 6.5)) +
    scale_fill_gradient(high="darkblue", low="white", name = "Probabilty of arrival", limits = c(0,1), guide = guide_colorbar(direction = "horizontal", barwidth = 13)) + #horizontal legend
    #scale_fill_viridis_c(option = "viridis",direction = -1, limits = c(0, 1), name = "Probability of arrival", guide = guide_colorbar(direction = "horizontal", barwidth = 13)) +
    scale_x_continuous(breaks = seq(8, 45, 8)) +
    scale_y_continuous(breaks = seq(-35, 6.5, 5)) +
    xlab('Longitude') +
    ylab('Latitude') +
    ggtitle(paste0('t = ', time_slices[k], ' BP / ', time_slices_BCAD_ch[k])) +
    theme(panel.background = element_rect(fill = "lightblue",
                                          colour = NA,
                                          size = 0.5,
                                          linetype = "solid"),
          panel.border = element_rect(color = "grey50", fill = NA, size = 0.5))
  
  plot_list[[k]] <- p
}

#Output
pdf(file=here('output','figures','fig_arrival_by_time.pdf'), width=8, height=10)
wrap_plots(plot_list, ncol = 3, nrow = 3, guides = "collect") +
  plot_annotation(theme = theme(legend.position = "bottom",
                                legend.box = "vertical"))
dev.off()

#-------------------------------------------------------------------------------
#FIGURE 3 -- Wombling to display significant boundaries

median_hex_dates_mod.i <- hex_area_win %>% 
  filter(area_ID %in% 1:47) %>% 
  mutate(median_date = med.model.womble,
         contains_sites = as.factor(case_when(area_ID %in% Hex_with_sites ~ 1, area_ID %in% Hex_without_sites ~ 0))) 

#With the data from hierarchical wombling model
post.model.womble_nab  <- out.comb.womble.model[,paste0('nabla[',1:117,']')]  %>% round()

#Extract differences in arrival times for wombling model
med.model.womble_nab  <- apply(post.model.womble_nab, 2, median)

#Extract proportion of MCMC sample differences which are significant over a specified time difference
prop_model_womble_nab  <- data.frame(x = 1:117,
                                     y = sapply(as.data.frame(post.model.womble_nab),
                                                prop_gthan_threshold,
                                                threshold = 600))

#Add info to edges dataframe
edge_info.i <- edge_info %>%
  mutate(mean_gradient = med.model.womble_nab, #50% quantile
         prob_BLV = prop_model_womble_nab$y, #% of distribution > specified threshold
         boundary = mapply(function(a, b) {intersection <- st_intersection(hex_area_win$geometry[[a]], hex_area_win$geometry[[b]])
         if (st_is_empty(intersection) || st_is(intersection, "MULTILINESTRING")) return(st_linestring()) else return(intersection)},
         edge_info$region1_id,
         edge_info$region2_id)) #shared boundary between two subareas

#If less than 15% chance of boundary -> 0 (to make boundary distinctions clearer)
edge_info.i$prob_BLV[edge_info.i$prob_BLV < 0.15] <- 0 

#Create nodes
nodes <- st_coordinates(hex_area_win$area_center)

#Create boundary segments
boundaries <- st_sf(prob_BLV = edge_info.i$prob_BLV,
                    geometry = st_sfc(edge_info.i$boundary)) 
st_crs(boundaries) <- 4326  # Set CRS to EPSG:4326 (WGS 84)


#Plot
nabla_plot <- ggplot(data = median_hex_dates_mod.i) +
  geom_sf(data = st_union(africa), fill = "#ECE6DD", color = "black", show.legend = FALSE) + 
  geom_sf(data = st_buffer(st_as_sf(st_union(sampling_win), crs = 4326), 20000), fill = "#ECE6DD", color="black", show.legend = FALSE) +
  geom_sf(data = st_as_sf(st_union(median_hex_dates_mod.i), crs = 4326), color="grey50", fill = "#ECE6DD", show.legend = FALSE) +#sampling window with coastal buffer
  geom_sf() + #hex grid
  geom_raster(data = dem_df, aes(x = lon, y = lat, fill = elevation)) +
  scale_fill_gradient(low = "white", high = "black", name = "Elevation (m)",limits = c(0, 3000)) +
  geom_sf(data = lakes, fill = "#4660dd", colour = "#4660dd", show.legend = FALSE) +
  geom_sf(data = boundaries, lwd=2.5, aes(alpha=prob_BLV), color = "red") +
  scale_alpha_continuous(range = c(0, 1), name="Probability of BLV", breaks=c(0, 0.25, 0.5, 0.75), labels=c("<0.15","0.25", "0.5", "0.75")) +
  coord_sf(xlim = c(7, 45),
           ylim = c(-35, 6.5)) +
  xlab('Longitude') +
  ylab('Latitude') +
  theme(panel.background = element_rect(fill = "lightblue",
                                        colour = NA,
                                        size = 0.5,
                                        linetype = "solid"),
        panel.border = element_rect(color = "grey50", fill = NA, size = 0.5),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13))
#---
##Repeat wombling but with spatial residues, phi[k]

#With the tactical simulation data from hierarchical wombling model
post.model.womble_nab_phi  <- out.comb.womble.model[,paste0('nabla_phi[',1:117,']')]  %>% round()

#Extract differences in spatial residues times for tactical wombling model
med.model.womble_nab_phi  <- apply(post.model.womble_nab_phi, 2, median)

#Extract proportion of MCMC sample differences which are significant over a specified time difference
prop_model_womble_nab_phi  <- data.frame(x = 1:117,
                                         y = sapply(as.data.frame(post.model.womble_nab_phi), 
                                                    prop_gthan_threshold, 
                                                    threshold = 700)) #change to 500 or 700 years for sensitivity analysis (supplementary figures)

#Add info to edges dataframe
edge_info_phi.i <- edge_info %>%
  mutate(mean_gradient = med.model.womble_nab_phi, #50% quantile
         prob_BLV = prop_model_womble_nab_phi$y, #% of distribution > specified threshold
         boundary = mapply(function(a, b) {intersection <- st_intersection(hex_area_win$geometry[[a]], hex_area_win$geometry[[b]])
         if (st_is_empty(intersection) || st_is(intersection, "MULTILINESTRING")) return(st_linestring()) else return(intersection)}, 
         edge_info$region1_id, 
         edge_info$region2_id)) #shared boundary between two subareas 

#If less than 15% chance of boundary -> 0 (to make boundary distinctions clearer)
edge_info_phi.i$prob_BLV[edge_info_phi.i$prob_BLV < 0.15] <- 0 

#Create boundary segments
boundaries_phi <- st_sf(prob_BLV = edge_info_phi.i$prob_BLV,
                        geometry = st_sfc(edge_info_phi.i$boundary))
st_crs(boundaries_phi) <- 4326  # Set CRS to EPSG:4326 (WGS 84)


#Plot
nabla_phiplot <- ggplot(data = median_hex_dates_mod.i) +
  geom_sf(data = st_union(africa), fill = "#ECE6DD", color = "black", show.legend = FALSE) + 
  geom_sf(data = st_buffer(st_as_sf(st_union(sampling_win), crs = 4326), 20000), fill = "#ECE6DD", color="black", show.legend = FALSE) +
  geom_sf(data = st_as_sf(st_union(median_hex_dates_mod.i), crs = 4326), color="grey50", fill = "#ECE6DD", show.legend = FALSE) +#sampling window with coastal buffer
  geom_sf() + #hex grid 
  geom_raster(data = dem_df, aes(x = lon, y = lat, fill = elevation)) +
  scale_fill_gradient(low = "white", high = "black", name = "Elevation (m)",limits = c(0, 3000)) +
  geom_sf(data = lakes, fill = "#4660dd", colour = "#4660dd", show.legend = FALSE) +
  geom_sf(data = boundaries_phi, lwd=2.5, aes(alpha=prob_BLV), color = "red",show.legend = FALSE) +
  scale_alpha_continuous(range = c(0, 1), name="Probability of BLV", breaks=c(0, 0.25, 0.5, 0.75), labels=c("<0.15","0.25", "0.5", "0.75")) +
  coord_sf(xlim = c(7, 45),
           ylim = c(-35, 6.5)) +
  xlab('Longitude') +
  ylab('Latitude') +
  theme(panel.background = element_rect(fill = "lightblue",
                                        colour = NA,
                                        size = 0.5,
                                        linetype = "solid"),
        panel.border = element_rect(color = "grey50", fill = NA, size = 0.5),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13))

# Save as PNG
combined_plot <- wrap_plots(nabla_plot, nabla_phiplot, ncol = 2, nrow = 1, guides = "collect") &
  theme(legend.position = "bottom",legend.box = "horizontal") &
  guides(fill=guide_colorbar(direction="horizontal"), alpha=guide_legend(direction="horizontal"))
ggsave(filename = here("output", "figures", "fig_wombling_boundaries.png"),plot = combined_plot,width = 13,height = 8,dpi = 300)