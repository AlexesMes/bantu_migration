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

source(here('src', 'grad_funcs.R'))
source(here('src', 'accuracy_precision.R'))
source(here('src','hex_areas.R'))
source(here('src','plotting_funcs.R'))

`%!in%` <- Negate(`%in%`)

#===============================================================================
## Load Data

#Load Observed Data
load(here('data','eastc14.RData')) #East and Southern Africa

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

#===============================================================================
##SUPPLEMENTARY FIGURES -- Exploring the data

#Dates per site plot
date_freq  <- dateInfo %>% 
  count(siteID, sort=TRUE) %>%
  rename(n_dates = n)
  
pdf(file=here('output','figures_supplementary','fig_num_dates.pdf'), width=6, height=5.5)
ggplot(date_freq, aes(x=n_dates)) +
  geom_bar() +
  scale_y_continuous(name="Number of sites", breaks=seq(0, 200, 25)) +
  scale_x_continuous(name="Number of dates per site", breaks=seq(1, 53, 4)) +
  theme_minimal() +
  theme(axis.text = element_text(size=14),
        text = element_text(size=14))
dev.off()

#------
#Effect of Hexagon Size on Site Coverage figure
#Under changing hex size, the proportion of areal hex units in the sampling window with sites changes
prop_units_df <- data.frame(d = numeric(), prop_with_sites = numeric())

for (d in seq(1, 15, 0.1)){
  hex_area_win <- hex_areas(sampling_win, cell_d = d)
  siteInfo$area_id <- as.integer(st_within(sites$geometry, hex_area_win$geometry))
  hex_with_sites <- length(unique(siteInfo$area_id))
  all_hex <- length(hex_area_win$area_ID)
  prop_with_sites <- hex_with_sites/all_hex
  prop_units_df <- rbind(prop_units_df, data.frame(d = d, prop_with_sites = prop_with_sites))
}

#Highlight the chosen d at which analysis occurs 
chosen_scale_d <- c(2.9, 3.8, 5.2)

#Select points at chosen d scale
prop_units_df <- prop_units_df %>% 
  dplyr::mutate(highlight = purrr::map_lgl(d,~ any(dplyr::near(.x, chosen_scale_d, tol = 1e-6))))

#Plot
pdf(here('output','figures_supplementary','fig_hexsize.pdf'),height=5,width=5.5)
ggplot(prop_units_df, aes(x = d, y = prop_with_sites)) +
  geom_line(color = "black") +
  geom_point(color = "black", size = 1.8) +
  geom_point(data = subset(prop_units_df, highlight),color = "blue",size = 3) +
  geom_vline(xintercept = chosen_scale_d,color = "blue",linewidth = 2, alpha=0.2) +
  geom_hline(data = subset(prop_units_df, highlight),aes(yintercept = prop_with_sites),color = "blue",linewidth = 2, alpha=0.2) +
  scale_x_continuous(breaks=seq(0,15,by=1))+
  scale_y_continuous(breaks=seq(0,0.8,0.1)) +
  labs(x = "Hexagon Size (d)", y = "Proportion of Hexagons with Sites") +
  theme_minimal()+
  theme(axis.text = element_text(size=14),
        text = element_text(size=14))
dev.off()

#---------------
##SUPPLEMENTARY TABLE -- Summary of areas, sites, dates
data_summary_df <- dateInfo %>% 
  group_by(area_id) %>% 
  summarise(
    n_sites = n_distinct(siteID),
    n_dates = n()) %>% 
  right_join(data.frame("area_id"=1:47), by="area_id") %>% 
  replace_na(list(n_sites=0, n_dates=0)) %>% 
  arrange(area_id)

write.csv(data_summary_df,file=here('output','tables','data_summary.csv'), row.names = FALSE)

#===============================================================================
##SUPPLEMENTARY FIGURES -- Environmental covariates

scaled_hex_clim_df <- as.data.frame(scaled_hex_clim_df)

#Plot net primary productivity
npp_plot <- ggplot(data = hex_area_win) +
  geom_sf(data = st_union(africa), fill = "#ECE6DD", color = "black") +
  geom_sf() + #hex grid
  geom_sf(aes(fill = scaled_hex_clim_df$npp)) +
  scale_fill_continuous_sequential(palette = "ag_GrnYl") +
  labs(fill="Scaled npp") +
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
        legend.position = "bottom")

#Plot elevation
elev_plot <- ggplot(data = hex_area_win) +
  geom_sf(data = st_union(africa), fill = "#ECE6DD", color = "black") +
  geom_sf() + #hex grid
  geom_sf(aes(fill = scaled_hex_clim_df$rugosity)) +
  scale_fill_continuous_sequential(palette = "Blues 3") +
  labs(fill="Scaled rugosity") +
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
        legend.position = "bottom")

pdf(file=here('output','figures_supplementary','fig_env_covs.pdf'), width=10, height=6)
grid.arrange(grobs=c(npp_plot, elev_plot), ncol=2, nrow=1, padding=0) # Define the layout for the plots
dev.off()

#===============================================================================
##WOMBLE MODEL WITH ENVIRONEMENTAL COVARIATES ---

##Load Data ----
load(here("output","Womble_Emodel_d38.RData"))

#Hex areas with and without out sites
Hex_with_sites <- unique(siteInfo$area_id)
Hex_without_sites <- which(rep(1:47) %!in% Hex_with_sites)

#Calculate credible interval
ci_95 = credible_interval(out_womble_model, 0.95)

#-------------------------------------------------------------------------------
##SUPPLEMENTARY Diagnostics: Traceplots and metric table 

#Create list to store recorded plots
plots <- vector("list", 47)

#Output
pdf(file = here("output", "figures_supplementary", "traceplots_d38.pdf"), width = 10, height = 15, onefile=TRUE)
par(mfrow = c(10,10), mar = c(2,0,2,0), oma = c(0, 0, 0, 0), mgp = c(1, 0.2, 0))
for (i in 1:47) {
  param_name <- paste0("a[", i, "]")
  plot.new()
  traceplot(out_womble_model[, param_name], main = TeX(paste0("$a[", i, "]$")), smooth = TRUE)
  plots[[i]] <- recordPlot()
}
dev.off()

#Diagnostics table
out.comb.tac_icar.model  <- do.call(rbind, out_womble_model)
post.a.model.i  <- out.comb.tac_icar.model[,paste0('a[',1:47,']')]  %>% round() 
med.model.i  <- apply(post.a.model.i, 2, median)

diagnostic_df <- data.frame(median_posterior = paste(med.model.i, "BP"),
                            HPDI95_low = paste(round(ci_95[1,]), "BP"),
                            HPDI95_high = paste(round(ci_95[2,]), "BP"),
                            rhat = round(rhat_womble_model$psrf[1:47,1],2),
                            ESS = round(ess_womble_model[1:47]))

write.csv(diagnostic_df,file=here("output","tables",'diagnostics_d38.csv'), row.names = TRUE)

#-------------------------------------------------------------------------------
##SUPPLEMENTARY FIGURE -- Beta Posteriors

pdf(file = here("output","figures_supplementary", "beta_posteriors.pdf"),width = 15, height = 5, onefile = TRUE)
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1.2))

plot_beta_density <- function(beta_name, label) {
  post.beta <- do.call(rbind, out_womble_model)[, beta_name] %>% round()
  dens.beta <- density(post.beta, bw = 5)
  med.beta  <- median(post.beta)
  plot(NULL,xlim = c(med.beta - 700, med.beta + 900),ylim = c(0, 0.01),xlab = "Cal BP",ylab = "Posterior Probability",cex.lab=1.5,cex.axis=1)
  polygon(c(dens.beta$x, rev(dens.beta$x)),c(rep(0, length(dens.beta$x)), rev(dens.beta$y)),border = NA,col = rgb(0, 0.4, 0, 0.5))
  abline(v = med.beta, lty = 2)
  axis(3, at = med.beta, labels = TeX(label), cex.axis=2)
}

plot_beta_density("beta0", "$\\beta_0$")
plot_beta_density("beta1", "$\\beta_1$")
plot_beta_density("beta2", "$\\beta_2$")
dev.off()

#Traceplots
pdf(file=here('output', 'figures_supplementary','traceplots_betas.pdf'), width=15, height=5)
par(mfrow=c(1,3))
traceplot(out_womble_model[,'beta0'], main=TeX('$\\beta_1$'),smooth=TRUE,cex.main=2)
traceplot(out_womble_model[,'beta1'], main=TeX('$\\beta_2$'),smooth=TRUE,cex.main=2)
traceplot(out_womble_model[,'beta2'], main=TeX('$\\beta_3$'),smooth=TRUE,cex.main=2)
dev.off()

#-------------------------------------------------------------------------------
##FIGURE -- Arrival time in time slices

out.comb.womble.model  <- do.call(rbind, out_womble_model)
post.model.womble <- out.comb.womble.model[,paste0('a[',1:47,']')]  %>% round() 

#Extract arrival times for tactical icar model
med.model.womble  <- apply(post.model.womble, 2, median)

time_slices <- seq(2200, 1200, -200)

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
    geom_sf(aes(fill = median_date, alpha=prop_threshold),color=NA) + #hex grid 
    coord_sf(xlim = c(7, 45),
             ylim = c(-35, 6.5)) +
    scale_fill_viridis_c(option="F", direction=-1, name = "Arrival time (in BP)   ", limits = c(1200, 2700), guide = guide_colorbar(direction = "horizontal", barwidth = 13)) + #horizontal legend
    scale_alpha_continuous(range = c(0, 1), guide="none", limits=c(0,1)) +  # Use for continuous alpha values
    scale_x_continuous(breaks = seq(8, 45, 8)) +
    scale_y_continuous(breaks = seq(-35, 6.5, 5)) +
    xlab('Longitude') +
    ylab('Latitude') +
    ggtitle(paste0('t = ', time_slices[k], ' BP')) +
    theme(panel.background = element_rect(fill = "lightblue",
                                          colour = NA,
                                          size = 0.5,
                                          linetype = "solid"),
          panel.border = element_rect(color = "grey50", fill = NA, size = 0.5))
  
  plot_list[[k]] <- p
}

#Output
pdf(file=here('output','figures','fig_arrival_time_slices.pdf'), width=8, height=7)
wrap_plots(plot_list, ncol = 3, nrow = 2, guides = "collect") +
  plot_annotation(theme = theme(legend.position = "bottom",
                                legend.box = "vertical"))
dev.off()

#-------------------------------------------------------------------------------
##FIGURE -- Median posterior arrival times with credible interval displayed

#Posterior summaries HDPI
low.model.womble  <- apply(post.model.womble, 2, quantile, probs = 0.025)
high.model.womble <- apply(post.model.womble, 2, quantile, probs = 0.975)

arrival_df <- hex_area_win %>% 
  filter(area_ID %in% 1:47) %>% 
  mutate(median_date = med.model.womble,
         lower_ci    = low.model.womble,
         upper_ci    = high.model.womble,
         contains_sites = as.factor(case_when(area_ID %in% Hex_with_sites ~ 1, area_ID %in% Hex_without_sites ~ 0))) 

#Plot -- median, lower CI, and upper CI dates
plot_list <- list() 
post_percent <- c("median_date", "lower_ci", "upper_ci")

for (k in 1:length(post_percent)) #all time slices
{
  p <- ggplot(data = arrival_df) +
    geom_sf(data = st_union(africa), fill = "grey70", color = "black") +
    geom_sf(data = st_buffer(st_as_sf(st_union(sampling_win), crs = 4326), 2000), color = "grey30",fill = "white", lwd=0.3) + #sampling window with coastal buffer
    geom_sf(aes(fill = .data[[post_percent[k]]]),color=NA) + #hex grid #alpha=contains_sites
    coord_sf(xlim = c(7, 45),
             ylim = c(-35, 6.5)) +
    scale_fill_viridis_c(option="viridis", direction=-1, name = "Arrival time (in BP)   ", limits = c(750, 3550), guide = guide_colorbar(direction = "horizontal", barwidth = 13)) + #horizontal legend
    scale_x_continuous(breaks = seq(8, 45, 8)) +
    scale_y_continuous(breaks = seq(-35, 6.5, 5)) +
    xlab('Longitude') +
    ylab('Latitude') +
    #geom_sf_label(aes(label = paste0(post_percent[k], "BP")), label.size  = NA, alpha = 0.4, size=3.5) +
    theme(panel.background = element_rect(fill = "lightblue",
                                          colour = NA,
                                          size = 0.5,
                                          linetype = "solid"),
          panel.border = element_rect(color = "grey50", fill = NA, size = 0.5))
  
  plot_list[[k]] <- p
}

#Output
pdf(file=here('output','figures','fig_arrival_posterior_summary.pdf'), width=10, height=8)
wrap_plots(plot_list[[1]], wrap_plots(plot_list[[2]], plot_list[[3]], ncol = 1), ncol = 2, widths = c(2.15, 1), guides = "collect") +
  plot_annotation(
    theme = theme(legend.position = "bottom",legend.box = "vertical"))
dev.off()

#-------------------------------------------------------------------------------
##FIGURE -- Arrival by time panels 

time_slices <-  seq(2600, 1200, -250)

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
    #scale_fill_gradient(high="black", low="white", name = "Probabilty of arrival", limits = c(0,1), guide = guide_colorbar(direction = "horizontal", barwidth = 13)) + #horizontal legend
    scale_fill_viridis_c(
      option = "rocket",direction = -1, limits = c(0, 1), name = "Probability of arrival", guide = guide_colorbar(direction = "horizontal", barwidth = 13)
      ) +
    scale_x_continuous(breaks = seq(8, 45, 8)) +
    scale_y_continuous(breaks = seq(-35, 6.5, 5)) +
    xlab('Longitude') +
    ylab('Latitude') +
    ggtitle(paste0('t = ', time_slices[k], ' BP')) +
    theme(panel.background = element_rect(fill = "lightblue",
                                          colour = NA,
                                          size = 0.5,
                                          linetype = "solid"),
          panel.border = element_rect(color = "grey50", fill = NA, size = 0.5))
  
  plot_list[[k]] <- p
}

#Output
pdf(file=here('output','figures','fig_arrival_by_time.pdf'), width=8, height=7)
wrap_plots(plot_list, ncol = 3, nrow = 2, guides = "collect") +
  plot_annotation(theme = theme(legend.position = "bottom",
                                legend.box = "vertical"))
dev.off()

#-------------------------------------------------------------------------------
##SUPPLEMENTARY FIGURE -- Animation: Arrival time in time slices
#Animate through time slices to display proportion of distribution before specified time

#Extract proportion of MCMC samples occurring before the specified time threshold (note more time divisions for animation)
time_slices2 <- seq(2500, 1000, -15) 

prop_model_icar2  <- lapply(time_slices2, 
                            function(t) data.frame(x = 1:47, 
                                                   y = sapply(as.data.frame(post.model.womble), 
                                                              prop_gthan_threshold, 
                                                              threshold = t)))

# Combine all time slices into one dataset
all_data <- lapply(seq_along(time_slices2), function(k) {
  hex_area_win %>%
    filter(area_ID %in% 1:47) %>%
    mutate(median_date = med.model.womble,
           prop_threshold = prop_model_icar2[[k]]$y,
           time_slice = -time_slices2[k]) %>% #minus sign ensures animation time proceeds in the correct direction
    st_as_sf()})

# Verify that all elements are sf objects
if (!all(sapply(all_data, inherits, "sf"))) {
  stop("One or more elements of all_data are not sf objects!")
}

# Combine all sf objects safely
all_data <- do.call(rbind, all_data)

#co-ordinates of area centroids used for label positions
label_data <- st_coordinates(all_data$area_center)

# Animated ggplot
p <- ggplot(data = all_data) +
  geom_sf(data = st_union(africa), fill = "#ECE6DD", color = "black") +
  geom_sf(data = st_buffer(st_as_sf(sampling_win, crs = 4326), 2000), aes(color = "grey30"), color=NA) +
  geom_sf(aes(fill = median_date, alpha = prop_threshold),color=NA) +
  coord_sf(xlim = c(7, 45),
           ylim = c(-35, 6.5)) +
  scale_fill_viridis_c(option = "F", direction = -1) +
  scale_alpha_continuous(range = c(0, 1)) +
  xlab('Longitude') +
  ylab('Latitude') +
  scale_x_continuous(breaks = seq(8, 45, 8)) +
  scale_y_continuous(breaks = seq(-35, 6.5, 5)) +
  # geom_text(data = label_data,
  #           aes(x = X, y = Y, label = paste0(all_data$median_date, " BP")),
  #           size = 4.5, alpha = 0.4) +  #Use geom_text instead of geom_sf_label
  theme(
    panel.background = element_rect(fill = "lightblue",colour = "lightblue",size = 0.5,linetype = "solid"),
    legend.position = "none",
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 15, face = "bold"),  # X-axis label
    axis.title.y = element_text(size = 15, face = "bold")   # Y-axis label
  ) +
  labs(title = "Time Slice: {closest_state} BP") #Note: use '{closest_state}' with 'transition_states()' and '{frame_time}' with 'transition_time()'

anim <- p +
  transition_states(time_slice, transition_length = 2, state_length = 1) +
  enter_fade() +
  exit_fade()

#Save the animation
gganimate::anim_save("arrival_time_slices_animation.gif",
          animation = gganimate::animate(anim, width = 1200, height = 900, fps = 4))

#-------------------------------------------------------------------------------
#FIGURE -- Wombling to display significant boundaries

post.a.model.i  <- out.comb.womble.model[,paste0('a[',1:47,']')]  %>% round() 
med.model.i  <- apply(post.a.model.i, 2, median)

median_hex_dates_mod.i <- hex_area_win %>% 
  filter(area_ID %in% 1:47) %>% 
  mutate(median_date = med.model.i,
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
                                                    threshold = 600))

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

#-------------------------------------------------------------------------------
##SUPPLEMENTARY FIGURE -- Probability Matrix of arrival times

source(here('src','plotting_funcs.R'))

pdf(file=here('output','figures_supplementary','fig_order_areas.pdf'),width=7,height=7.5)
orderPPlot(post.model.womble, name.vec=paste("Area",as.character(1:47)))
dev.off()

#----
#Select Area 23 and neighbours specifically (to examine likelihood of leapfrog transmission)
interest_areas <- c(33,24,14,28,18,9,13,5,17,8,12,23)

out.comb.womble.model.23  <- do.call(rbind, out_womble_model)
post.model.womble.23 <- out.comb.womble.model[,paste0('a[',interest_areas,']')]  %>% round() 

pdf(file=here('output','figures_supplementary','fig_order_area23.pdf'),width=7,height=7.5)
orderPPlot(post.model.womble.23, name.vec=paste("Area",as.character(interest_areas)))
dev.off()

#Difference Matrix plot of arrival times of area 23 and neighbours
make_mat <- function(n) {
  mat <- matrix(0, n, n)
  counter <- n + 1
  
  for (j in 1:n) {
    mat[j, j] <- j
    if (j < n) {
      mat[(j + 1):n, j] <- counter:(counter + n - j - 1)
      counter <- counter + n - j
    }
  }
  mat
}
mat <- make_mat(length(interest_areas))

pdf(file=here('output','figures_supplementary','fig_dif_matrix_area23.pdf'),width=16,height=11)
layout(mat)
par(mar=c(0,0,0,0))
for (i in interest_areas){plot(NULL,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",axes=F);text(0.5,0.5,paste('Area',as.character(i)),cex=3)}
par(mar=c(3,0,0,1))
for (i in seq_along(interest_areas)){
  for (j in seq_along(interest_areas)){
    if (i < j)
    {diffDens(post.model.womble.23[,i],post.model.womble.23[,j],xlim=c(-1400,1400),prob=0.9)}
  }
}
dev.off()

#-------------------------------------------------------------------------------
##SUPPLEMENTARY FIGURE  -- Comparing posterior distributions between models with and without environmental covariates

#Extract values for ICAR model with environmental covariates 
model.i.long <- data.frame(value = as.numeric(post.model.womble),
                            area = rep(c(1:47), each=nrow(post.model.womble))) %>% 
  mutate(area = factor(area, levels=paste0(c(1:47)), ordered=TRUE))
tmp.a.model.i = extract(out_womble_model)

#Load and extract ICAR model without environmental covariates 
load(here('output', 'Womble_Emodel_d38_nocov.RData'))

out.comb.womble.ii.model  <- do.call(rbind, out_womble_model)
post.model.ii.womble <- out.comb.womble.ii.model[,paste0('a[',1:47,']')]  %>% round() 

model.ii.long <- data.frame(value = as.numeric(post.model.ii.womble),
                           area = rep(c(1:47), each=nrow(post.model.ii.womble))) %>% 
  mutate(area = factor(area, levels=paste0(c(1:47)), ordered=TRUE))
tmp.a.model.ii = extract(out_womble_model)

#--------
#Plot
pdf(file=here('output','figures_supplementary','fig_compare_posteriors_ak.pdf'), width=10, height=15, pointsize=4)
par(mar = c(5, 5, 4, 2))   #pad space around plot
plot(NULL, xlim=c(3450, 500), ylim=c(1.5, 46.5), xlab=paste('Arrival time (BP),', TeX('$a_k$')), ylab=paste('Area,', TeX('$k$')), cex.lab = 2, axes=F)

iseq.a = seq(1,by=1,length.out=47)
abline(h=seq(1,by=1,length.out=47), col='lightgrey')

counter <- 1 #indexing counter
for (i in 1:47) #all even hex areas #for odd hex areas use: seq(1,81,2)
{
  #Plot bars from sim 2 and sim 7 in area i
  post.bar(tmp.a.model.i[,i], i=iseq.a[counter], h=0.5, barcolours=barcolours1)
  post.bar(tmp.a.model.ii[,i], i=iseq.a[counter]+0.3, h=0.5, barcolours=barcolours2)
  counter <- counter + 1
}

axis(2, at=iseq.a, labels = paste0(1:47), las=2, cex.axis=1.7) 
axis(1, at = BCADtoBP(c(-1500, -1100, -700, -300, 100, 500, 900, 1300)), labels=c('1500BC','1100BC','700BC','300BC','100AD','500AD','900AD','1300AD'), tck=-0.01, cex.axis=1.7)
axis(3, at = seq(3400, 500, -400), labels=paste0(seq(3400, 500, -400),'BP'), tck=-0.01, cex.axis=1.7)
axis(1, at = BCADtoBP(c(-1400,-1300,-1200,-1000,-900,-800,-600,-500,-400,-200,-100,1,200,300,400,600,800,1000,1100,1200,1400)), labels=NA, tck=-0.01) #Minor tick marks
axis(3, at = seq(3400, 500, -50), labels=NA, tck=-0.01) #Minor tick marks
box()
#Legend
post.bar(c(1300,1200,1100,1000,900,800,700), i=44.5, h=0.9)
text(x=750, y=45, "50% HPDI", cex=2)
text(x=1070, y=45,"95% HPDI", cex=2)
text(x=900, y=44, "Median Posterior", cex=2)
rect(xleft=1450, xright=500, ybottom=42, ytop=45.5, border="darkgrey", col=NA, lwd=2.5)
text(x=900, y=43, "ICAR Model with covariates", col="black", cex=2)
segments(x0=1400, x1=1350, y0=43, col="dodgerblue", lwd=4)
text(x=900, y=42.5, "ICAR Model without covariates", col="black", cex=2)
segments(x0=1400, x1=1350, y0=42.5, col="orchid", lwd=4)
dev.off()

#-------------------------------------------------------------------------------
##SUPPLEMENTARY FIGURE  -- Arrival times when data is reduced

#Areas with Wanyika data
wan_areas <- siteInfo %>%filter(dataorigin=="Wanyika")
wan_areas <- unique(wan_areas$area_id)

#Load and extract ICAR model with 15% of data missing (by filtering Wanyika dataset for A,B graded dates) 
load(here('data','eastc14_wanAB.RData')) #Reduced data
load(here('output', 'Womble_Emodel_d38_reduce.RData'))#Inferred results

out_womble_model_reduc <- out_womble_model
out.comb.womble.model.reduc  <- do.call(rbind, out_womble_model_reduc)
post.model.womble.reduc <- out.comb.womble.model.reduc[,paste0('a[',1:47,']')]  %>% round() 

#Extract arrival times for tactical icar model
med.model.womble  <- apply(post.model.womble.reduc, 2, median)

time_slices <- seq(2000, 1200, -400)

#Extract proportion of MCMC samples occurring before the specified time threshold
prop_model_womble  <- lapply(time_slices,
                             function(t) data.frame(x = 1:47,
                                                    y = sapply(as.data.frame(post.model.womble.reduc), 
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
           prop_threshold = prop_model_womble[[k]]$y,
           contains_sites = as.factor(case_when(area_ID %in% Hex_with_sites ~ 1, area_ID %in% Hex_without_sites ~ 0))) 
  #Plot
  p <- ggplot(data = prop_threshold_womble) +
    geom_sf(data = st_union(africa), fill = "#ECE6DD", color = "black") +
    geom_sf(data = st_buffer(st_as_sf(st_union(sampling_win), crs = 4326), 2000), color = "grey30", lwd=0.3) + #sampling window with coastal buffer
    geom_sf(aes(fill = median_date, alpha=prop_threshold),color=NA) + #hex grid #alpha=contains_sites
    coord_sf(xlim = c(7, 45),
             ylim = c(-35, 6.5)) +
    scale_fill_viridis_c(option="F", direction=-1, name = "Arrival time (in BP)   ", limits = c(1200, 2700), guide = guide_colorbar(direction = "horizontal", barwidth = 13)) + #horizontal legend
    scale_alpha_continuous(range = c(0, 1), guide="none") +  # Use for continuous alpha values
    scale_x_continuous(breaks = seq(8, 45, 8)) +
    scale_y_continuous(breaks = seq(-35, 6.5, 5)) +
    xlab('Longitude') +
    ylab('Latitude') +
    ggtitle(paste0('t = ', time_slices[k], ' BP')) +
    theme(panel.background = element_rect(fill = "lightblue",
                                          colour = NA,
                                          size = 0.5,
                                          linetype = "solid"),
          panel.border = element_rect(color = "grey50", fill = NA, size = 0.5))
  
  plot_list[[k]] <- p
}

#Output
pdf(file=here('output','figures_supplementary','fig_arrival_time_reduced.pdf'), width=8, height=4)
wrap_plots(plot_list, ncol = 3, nrow = 1, guides = "collect") +
  plot_annotation(theme = theme(legend.position = "bottom",
                                legend.box = "vertical"))
dev.off()

#-------------------------------------------------------------------------------
##SUPPLEMENTARY FIGURE  -- Comparing posterior distributions between models with reduced data

#Extract values for ICAR model with 15% of data missing (by filtering Wanyika dataset for A,B graded dates)
model.iii.long <- data.frame(value = as.numeric(post.model.womble.reduc),
                             area = rep(c(1:47), each=nrow(post.model.womble.reduc))) %>% 
  mutate(area = factor(area, levels=paste0(c(1:47)), ordered=TRUE))
tmp.a.model.iii = extract(out_womble_model_reduc)

#--------
#Plot
pdf(file=here('output','figures_supplementary','fig_compare_posteriors_ak_reduced.pdf'), width=10, height=5, pointsize=4)
par(mar = c(5, 5, 4, 2))   #pad space around plot
plot(NULL, xlim=c(3450, 500), ylim=c(0.5, 12.5), xlab=paste('Arrival time (BP),', TeX('$a_k$')), ylab=paste('Area,', TeX('$k$')), cex.lab = 2, axes=F)

iseq.a = seq(1,by=1,length.out=12)
abline(h=seq(1,by=1,length.out=12), col='lightgrey')

counter <- 1 #indexing counter
for (i in wan_areas) #all even hex areas #for odd hex areas use: seq(1,81,2)
{
  #Plot bars from sim 2 and sim 7 in area i
  post.bar(tmp.a.model.i[,i], i=iseq.a[counter], h=0.5, barcolours=barcolours1)
  post.bar(tmp.a.model.iii[,i], i=iseq.a[counter]+0.3, h=0.5, barcolours=barcolours2)
  counter <- counter + 1
}

axis(2, at=iseq.a, labels = paste0(wan_areas), las=2, cex.axis=1.7) 
axis(1, at = BCADtoBP(c(-1500, -1100, -700, -300, 100, 500, 900, 1300)), labels=c('1500BC','1100BC','700BC','300BC','100AD','500AD','900AD','1300AD'), tck=-0.01, cex.axis=1.7)
axis(3, at = seq(3400, 500, -400), labels=paste0(seq(3400, 500, -400),'BP'), tck=-0.01, cex.axis=1.7)
axis(1, at = BCADtoBP(c(-1400,-1300,-1200,-1000,-900,-800,-600,-500,-400,-200,-100,1,200,300,400,600,800,1000,1100,1200,1400)), labels=NA, tck=-0.01) #Minor tick marks
axis(3, at = seq(3400, 500, -50), labels=NA, tck=-0.01) #Minor tick marks
box()
#Legend
post.bar(c(1300,1200,1100,1000,900,800,700), i=10, h=0.9)
text(x=750, y=10.5, "50% HPDI", cex=2)
text(x=1070, y=10.5,"95% HPDI", cex=2)
text(x=900, y=9.5, "Median Posterior", cex=2)
rect(xleft=1450, xright=400, ybottom=7.5, ytop=11, border="darkgrey", col=NA, lwd=2.5)
text(x=850, y=8.5, "ICAR Model with covariates", col="black", cex=1.85)
segments(x0=1400, x1=1350, y0=8.5, col="dodgerblue", lwd=4)
text(x=860, y=8, "ICAR Model with covs, Wanyika filtered", col="black", cex=1.85)
segments(x0=1400, x1=1350, y0=8, col="orchid", lwd=4)
dev.off()

#-------------------------------------------------------------------------------
##SUPPLEMENTARY FIGURE  -- Arrival times when scale is finer: d=2.9

#Load and extract ICAR model when scale is finer, d=2.9 
load(here('data','eastc14_d29.RData')) #Reduced data
load(here('output', 'Womble_Emodel_d29.RData'))#Inferred results

out_womble_model_d29 <- out_womble_model
out.comb.womble.model.d29  <- do.call(rbind, out_womble_model_d29)
post.model.womble.d29 <- out.comb.womble.model.d29[,paste0('a[',1:83,']')]  %>% round() 

#Extract arrival times for tactical icar model
med.model.womble  <- apply(post.model.womble.d29, 2, median)

time_slices <- seq(2200, 1200, -200)

#Extract proportion of MCMC samples occurring before the specified time threshold
prop_model_womble  <- lapply(time_slices,
                             function(t) data.frame(x = 1:83,
                                                    y = sapply(as.data.frame(post.model.womble.d29), 
                                                               prop_gthan_threshold, 
                                                               threshold = t)))
#Plot
plot_list <- list() 

for (k in 1:length(time_slices)) #all time slices
{
  #Save in data structure
  prop_threshold_womble <- hex_area_win %>%
    filter(area_ID %in% 1:83) %>%
    mutate(median_date = med.model.womble,
           prop_threshold = prop_model_womble[[k]]$y,
           contains_sites = as.factor(case_when(area_ID %in% Hex_with_sites ~ 1, area_ID %in% Hex_without_sites ~ 0))) 
  #Plot
  p <- ggplot(data = prop_threshold_womble) +
    geom_sf(data = st_union(africa), fill = "grey70", color = "black") +
    geom_sf(data = st_buffer(st_as_sf(st_union(sampling_win), crs = 4326), 2000), color = "grey30", lwd=0.3) + #sampling window with coastal buffer
    geom_sf(aes(fill = median_date, alpha=prop_threshold),color=NA) + #hex grid #alpha=contains_sites
    coord_sf(xlim = c(7, 45),
             ylim = c(-35, 6.5)) +
    scale_fill_viridis_c(option="F", direction=-1, name = "Arrival time (in BP)   ", limits = c(1200, 2860), guide = guide_colorbar(direction = "horizontal", barwidth = 13)) + #horizontal legend
    scale_alpha_continuous(range = c(0, 1), guide="none") +  # Use for continuous alpha values
    scale_x_continuous(breaks = seq(8, 45, 8)) +
    scale_y_continuous(breaks = seq(-35, 6.5, 5)) +
    xlab('Longitude') +
    ylab('Latitude') +
    ggtitle(paste0('t = ', time_slices[k], ' BP')) +
    theme(panel.background = element_rect(fill = "lightblue",
                                          colour = NA,
                                          size = 0.5,
                                          linetype = "solid"),
          panel.border = element_rect(color = "grey50", fill = NA, size = 0.5))
  
  plot_list[[k]] <- p
}

#Output
pdf(file=here('output','figures_supplementary','fig_arrival_time_d29.pdf'), width=8, height=7)
wrap_plots(plot_list, ncol = 3, nrow = 2, guides = "collect") +
  plot_annotation(theme = theme(legend.position = "bottom",
                                legend.box = "vertical"))
dev.off()
#-------------------------------------------------------------------------------
##SUPPLEMENTARY FIGURE -- Beta Posteriors when scale d=2.9

pdf(file = here("output","figures_supplementary", "beta_posteriors_d29.pdf"),width = 15, height = 5, onefile = TRUE)
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1.2))

plot_beta_density <- function(beta_name, label) {
  post.beta <- do.call(rbind, out_womble_model)[, beta_name] %>% round()
  dens.beta <- density(post.beta, bw = 5)
  med.beta  <- median(post.beta)
  plot(NULL,xlim = c(med.beta - 700, med.beta + 900),ylim = c(0, 0.01),xlab = "Cal BP",ylab = "Posterior Probability",cex.lab=1.5,cex.axis=1)
  polygon(c(dens.beta$x, rev(dens.beta$x)),c(rep(0, length(dens.beta$x)), rev(dens.beta$y)),border = NA,col = rgb(0, 0.4, 0, 0.5))
  abline(v = med.beta, lty = 2)
  axis(3, at = med.beta, labels = TeX(label), cex.axis=2)
}

plot_beta_density("beta0", "$\\beta_0$")
plot_beta_density("beta1", "$\\beta_1$")
plot_beta_density("beta2", "$\\beta_2$")
dev.off()

#-------------------------------------------------------------------------------
##SUPPLEMENTARY FIGURE -- Wombling to display significant boundaries when d=2.9

post.a.model.i  <- out.comb.womble.model[,paste0('a[',1:83,']')]  %>% round() 
med.model.i  <- apply(post.a.model.i, 2, median)

median_hex_dates_mod.i <- hex_area_win %>% 
  filter(area_ID %in% 1:83) %>% 
  mutate(median_date = med.model.i,
         contains_sites = as.factor(case_when(area_ID %in% Hex_with_sites ~ 1, area_ID %in% Hex_without_sites ~ 0))) 

#With the data from hierarchical wombling model
post.model.womble_nab  <- out.comb.womble.model[,paste0('nabla[',1:216,']')]  %>% round()

#Extract differences in arrival times for wombling model
med.model.womble_nab  <- apply(post.model.womble_nab, 2, median)

#Extract proportion of MCMC sample differences which are significant over a specified time difference
prop_model_womble_nab  <- data.frame(x = 1:216,
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
post.model.womble_nab_phi  <- out.comb.womble.model[,paste0('nabla_phi[',1:216,']')]  %>% round()

#Extract differences in spatial residues times for tactical wombling model
med.model.womble_nab_phi  <- apply(post.model.womble_nab_phi, 2, median)

#Extract proportion of MCMC sample differences which are significant over a specified time difference
prop_model_womble_nab_phi  <- data.frame(x = 1:216,
                                         y = sapply(as.data.frame(post.model.womble_nab_phi), 
                                                    prop_gthan_threshold, 
                                                    threshold = 600))

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
ggsave(filename = here("output", "figures_supplementary", "fig_wombling_boundaries_d29.png"),plot = combined_plot,width = 13,height = 8,dpi = 300)
