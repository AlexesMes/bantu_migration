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
library(coda)
library(graphics)
library(ggthemes)

rm(list = ls())
`%!in%` <- Negate(`%in%`)

#===============================================================================
## Load Data ----

source(here('src', 'grad_funcs.R'))
source(here('src', 'hex_areas.R'))
source(here('src', 'accuracy_precision.R'))

load(here('data','trig.RData')) #Load nodes and edges between hex area centroids
load(here('data','sample_window.RData')) 

#===============================================================================
##FIGURE 1 -- Map of sampling window, with delaunay triangulation and hex areas 

#Convert edge_info to sf LINESTRING geometry
edge_info_sf <- edge_info %>%
  # Keep only the filtered transitions you already created
  mutate(
    geometry = purrr::pmap(
      list(region1_x, region1_y, region2_x, region2_y),
      ~ st_linestring(matrix(c(..1, ..2, ..3, ..4), ncol = 2, byrow = TRUE))
    )
  ) %>%
  st_as_sf(crs = 3035)  #projected CRS

#Plot
y <- ggplot() +
      geom_sf(data = st_buffer(sampling_win_proj, 40000), fill=NA, color = "grey50", linewidth=1) +
      geom_sf(data = hex_area_win_proj) +
      geom_sf_text(data = hex_area_win_proj, aes(label = area_ID), size=4, alpha=0.8) +
      geom_sf(data = hex_area_win_proj$area_center, size=2, alpha=0.6, color = "purple") +
      geom_sf(data = edge_info_sf, colour='purple', size=0.8, alpha=0.3) +
      labs(x = "Longitude", y = "Latitude") +
      theme(panel.background = element_rect(fill = "lightblue",
                                            colour = "lightblue",
                                            size = 0.5,
                                            linetype = "solid"),
            legend.position = "none")

#Output
pdf(file=here('output','figures','sample_win_hex.pdf'), width=15, height=8)
grid.arrange(y, ncol=1, padding=0)
dev.off()

#===============================================================================
#SIMULATION 1: TACTICAL WOA SIMULATION (with calibrated radiocarbon dates)

#Load data
load(here('data','tactical_sim_woa2.RData')) #simulated data
load(here('output', 'Womblemodel_tactical_woa.RData')) #inferred data

#------------
#Hex areas with and without out sites
Hex_with_sites <- unique(siteInfo$area_id)
Hex_without_sites <- which(rep(1:81) %!in% Hex_with_sites)

#------------
#Chain health check
traceplot(out_womble_model[,'a[20]'], main=TeX('$a$'), smooth=TRUE) #region 20 as an example
traceplot(out_womble_model[,'theta[100]'], main=TeX('$theta$'), smooth=TRUE)

#------------
#Calculate accuracy and precision
ci_95 = credible_interval(out_womble_model, 0.95)
sim_a <- constants$true_a
sim1_model_accuracy = accuracy(sim_a, ci_95)
sim1_model_precision = precision(sim_a, ci_95)


#-------------------------------------------------------------------------------
##FIGURE 2 -- Map of sites and sampling window
#Plot
site_map <- ggplot(data = hex_area_win_proj) +
    geom_sf(data = sampling_win_proj, color = "grey50") +  # sampling window border
    geom_sf() +
    geom_sf(data = sites_sf, size = 2, alpha = 0.5) +  # sites
    geom_sf_label(aes(label = area_ID)) +                     # area labels
    theme(
      panel.background = element_rect(fill = "lightblue", colour = "lightblue"),
      legend.title = element_blank(),
      legend.position = "bottom")

#Output
pdf(file=here('output','figures','sim1_sites.pdf'), width=15, height=8)
grid.arrange(site_map, ncol=1, padding=0)
dev.off()

#------------
##FIGURE 3 -- Map of Inferred arrival times
out.comb.tac_icar.model  <- do.call(rbind, out_womble_model)
post.a.model.i  <- out.comb.tac_icar.model[,paste0('a[',1:81,']')]  %>% round() 
med.model.i  <- apply(post.a.model.i, 2, median) #Extract arrival times for tactical icar model

median_hex_dates_mod.i <- hex_area_win_proj %>% 
  filter(area_ID %in% 1:81) %>% 
  mutate(median_date = med.model.i,
         contains_sites = as.factor(case_when(area_ID %in% Hex_with_sites ~ 1, area_ID %in% Hex_without_sites ~ 0))) 

#Plot
modi <- ggplot(data = median_hex_dates_mod.i) +
  geom_sf(data = st_buffer(sampling_win_proj, 40000), aes(color = "grey50")) + #sampling window with coastal buffer
  geom_sf(aes(fill = median_date)) + #hex grid #alpha=contains_sites
  scale_fill_viridis_c(option="F", direction=-1) +
  scale_alpha_manual(values=c(0.45, 1)) +
  xlab('Longitude') +
  ylab('Latitude') +
  geom_sf_label(aes(label = paste0(median_date, "BP")), label.size  = NA, alpha = 0.4, size=3.5) + #hex grid labels #label = ifelse(contains_sites==0, NA, paste0(median_date, "BP")))
  theme(panel.background = element_rect(fill = "lightblue",
                                        colour = "lightblue",
                                        size = 0.5,
                                        linetype = "solid"),
        legend.position = "none")


#Output
pdf(file=here('output','figures','sim1_arrivaltime.pdf'), width=15, height=8)
grid.arrange(modi, ncol=1, padding=0)
dev.off()

#------------
##FIGURE 4 -- Posterior distributions of arrival times
#For model (i) and (ii) select parameters a and b (i.e. start and end date of occupation in the region)
sim_a <- constants$true_a
sim_b <- constants$true_b

#Functions
extract <- function(x)
{
  tmp = do.call(rbind, x)
  tmp2 = tmp[ , grep('^a\\[',colnames(tmp))]
  qta = apply(tmp2, 2, quantile, prob=c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1))
  return(qta)
}

post.bar <- function(x, i, h, a)
{
  rect(xleft = x[2], xright = x[6], ybottom = i - h/5, ytop = i + h/5, border = NA, col = "skyblue") # 95% interval rectangle
  segments(x[2], i-h/3.5, x[2], i+h/3.5, lwd = 2, col = "skyblue") # horizontal ticks for 95%
  segments(x[6], i-h/3.5, x[6], i+h/3.5, lwd = 2, col = "skyblue") # horizontal ticks for 95%
  
  rect(xleft=x[3], xright=x[5], ybottom=i-h/3, ytop=i+h/3, border=NA, col="dodgerblue") #50% interval
  segments(x[3], i-h/2.5, x[3], i+h/2.5, lwd = 2, col = "dodgerblue")   # horizontal ticks for 50%
  segments(x[5], i-h/2.5, x[5], i+h/2.5, lwd = 2, col = "dodgerblue")   # horizontal ticks for 50%
  
  points(x[4], i, pch = 16, col = "darkblue", cex = 2) #posterior median 
  points(a, i, pch = 4, col = "darkgreen", cex = 2, lwd = 2) #simulated (true) arrival time 
}

#Plot
pdf(file=here('output','figures','sim1_posteriors.pdf'), width=10, height=15, pointsize=4)
par(mar = c(5, 5, 4, 2))   #pad space around plot
plot(NULL, xlim=c(6800, 4000), ylim=c(3, 79), xlab=paste('Arrival time (BP),', TeX('$a_k$')), ylab=paste('Area,', TeX('$k$')), cex.lab = 2, axes=F)
tmp.a = extract(out_womble_model)
iseq.a = seq(1,by=1,length.out=81)
abline(h=seq(1,by=1,length.out=81), col='lightgrey')

counter <- 1 #indexing counter
for (i in c(1:81)) #all relevant hex areas
{
  #Plot bar in area i
  post.bar(tmp.a[,i], i=iseq.a[counter], h=0.5, a= sim_a[[i]])
  counter <- counter + 1
}

axis(2, at=iseq.a, labels = paste0(c(1:81)), las=2, cex.axis=1.7)
axis(1, at = BCADtoBP(c(-4900, -4700, -4500, -4300, -4100, -3900, -3700, -3500, -3300, -3100, -2900, -2700, -2400, -2200, -2000)), labels=c('4900BC','4700BC','4500BC','4300BC', '4100BC', '3900BC', '3700BC', '3500BC', '3300BC', '3100BC', '2900BC', '2700BC', '2400BC','2200BC', '2000BC'), tck=-0.01, cex.axis=1.7)
axis(3, at = seq(6800, 4000, -200), labels=paste0(seq(6800, 4000, -200),'BP'), tck=-0.01, cex.axis=1.7)
axis(1, at = BCADtoBP(c(-4800, -4600, -4400, -4200, -4000, -3800, -3600, -3400, -3200, -3000, -2800, -2600, -2500, -2300, -2100)), labels=NA, tck=-0.01) #Minor tick marks
axis(3, at = seq(6800, 4000, -50), labels=NA, tck=-0.01) #Minor tick marks
box()
#Legend
post.bar(c(6900,6800,6600,6500,6400,6200,6100), i=77, h=0.9, a=6750)
text(x=6550, y=78, "95% HPDI", cex=1.5)
text(x=6300, y=78,"50% HPDI", cex=1.5)
text(x=6350, y=76, "Median Posterior", cex=1.5)
text(x=6700, y=76, "Simulated value", cex=1.5)
rect(xleft=6850, xright=6150, ybottom=75, ytop=79, border="darkgrey", col=NA, lwd=2)
theme(legend.position = "none")
dev.off()

#------------
##FIGURE  -- Map of Wombling Boundaries (highlighting significant boundaries)
#
# #With the tactical simulation data from hierarchical wombling model
# post.model.tac_womble_nab  <- out.comb.tac_icar.model[,paste0('nabla[',1:208,']')]  %>% round()
# 
# #Extract differences in arrival times for tactical wombling model
# med.model.tac_womble_nab  <- apply(post.model.tac_womble_nab, 2, median)
# 
# #Extract proportion of MCMC sample differences which are significant over a specified time difference
# prop_model_tac_womble_nab  <- data.frame(x = 1:208,
#                                          y = sapply(as.data.frame(post.model.tac_womble_nab),
#                                                     prop_gthan_threshold,
#                                                     threshold = 600))
# 
# #Add info to edges dataframe
# edge_info.i <- edge_info %>%
#   mutate(mean_gradient = med.model.tac_womble_nab, #50% quantile
#          prob_BLV = prop_model_tac_womble_nab$y, #% of distribution > specified threshold
#          boundary = mapply(function(a, b) {intersection <- st_intersection(hex_area_win_proj$geometry[[a]], hex_area_win_proj$geometry[[b]])
#          if (st_is_empty(intersection) || st_is(intersection, "MULTILINESTRING")) return(st_linestring()) else return(intersection)},
#          edge_info$region1_id,
#          edge_info$region2_id, SIMPLIFY = FALSE)) #shared boundary between two subareas
# 
# #Create nodes
# nodes <- st_coordinates(hex_area_win_proj$area_center)
# 
# #Create boundary segments
# boundaries <- st_sf(prob_BLV = edge_info.i$prob_BLV,
#                     geometry = st_sfc(edge_info.i$boundary)) #lapply(edge_info.i$boundary[[a]], st_coordinates(a))
# st_crs(boundaries) <- 3035  # Set CRS for correct Europe projection
# 
# #Plot
# pdf(file=here('output','figures','sim1_womble.pdf'))
# ggplot(data = median_hex_dates_mod.i) +
#   geom_sf(data = st_buffer(sampling_win_proj, 40000), fill = "grey80", color = "grey40") + #sampling window with coastal buffer
#   geom_sf(aes(alpha=0.01), color = "grey60") + scale_alpha(range = c(0, 1)) + #hex grid
#   # geom_segment(data=edge_info.i$boundary, aes(#x= region1_x, y= region1_y, xend= region2_x, yend= region2_y, alpha= prob_BLV), color="red", size=2) +
#   geom_sf(data = boundaries, lwd=3, aes(alpha=prob_BLV), color = "red") +
#   geom_sf(data = hex_area_win_proj$area_center, size=2, alpha=1, color = "grey40") + #hex-centers
#   scale_alpha_continuous(range = c(0, 1)) +  # Use for continuous alpha values
#   xlab('Longitude') +
#   ylab('Latitude') +
#   ggtitle(paste0('c = 400', ' years')) +
#   theme(panel.background = element_rect(fill = "lightblue",
#                                         colour = "lightblue",
#                                         size = 0.5,
#                                         linetype = "solid"),
#         legend.position = "none")
# dev.off()

#===============================================================================
#SIMULATION 2: TACTICAL ICAR SIMULATION (with calibrated radiocarbon dates and two covariates)

#Load data
load(here('data','tactical_sim_icar2.RData')) #simulated data
load(here('output', 'Womblemodel_tactical_icar2.RData')) #inferred data

#------------
#Hex areas with and without out sites
Hex_with_sites <- unique(siteInfo$area_id)
Hex_without_sites <- which(rep(1:81) %!in% Hex_with_sites)

#------------
#Chain health check
traceplot(out_womble_model[,'a[20]'], main=TeX('$a$'), smooth=TRUE) #region 20 as an example
traceplot(out_womble_model[,'theta[100]'], main=TeX('$theta$'), smooth=TRUE)

#------------
#Calculate accuracy and precision
ci_95 = credible_interval(out_womble_model, 0.95)
sim_a <- constants$true_a
sim2_model_accuracy = accuracy(sim_a, ci_95)
sim2_model_precision = precision(sim_a, ci_95)

#-------------------------------------------------------------------------------
##FIGURE 5 -- Map of sites and sampling window
#Plot
site_map <- ggplot(data = hex_area_win_proj) +
  geom_sf(data = sampling_win_proj, color = "grey50") +  # sampling window border
  geom_sf() +
  geom_sf(aes(fill = env_type)) + # color by combined type
  scale_fill_manual(
    values = c("Neither" = "grey90", "Forest Only" = "forestgreen", "Water Only" = "skyblue", "Forest & Water" = "hotpink4")) +
  geom_sf(data = sites_sf, size = 2, alpha = 0.5) +  # sites
  geom_sf_label(aes(label = area_ID)) +                     # area labels
  theme(
    panel.background = element_rect(fill = "lightblue", colour = "lightblue"),
    legend.title = element_blank(),
    legend.position = "bottom")

#Output
pdf(file=here('output','figures','sim2_sites.pdf'), width=15, height=8)
grid.arrange(site_map, ncol=1, padding=0)
dev.off()

#------------
##FIGURE 6 -- Map of Inferred arrival times
out.comb.tac_icar.model  <- do.call(rbind, out_womble_model)
post.a.model.i  <- out.comb.tac_icar.model[,paste0('a[',1:81,']')]  %>% round() 
med.model.i  <- apply(post.a.model.i, 2, median) #Extract arrival times for tactical icar model

median_hex_dates_mod.i <- hex_area_win_proj %>% 
  filter(area_ID %in% 1:81) %>% 
  mutate(median_date = med.model.i,
         contains_sites = as.factor(case_when(area_ID %in% Hex_with_sites ~ 1, area_ID %in% Hex_without_sites ~ 0))) 

#Plot
modi <- ggplot(data = median_hex_dates_mod.i) +
  geom_sf(data = st_buffer(sampling_win_proj, 40000), aes(color = "grey50")) + #sampling window with coastal buffer
  geom_sf(aes(fill = median_date)) + #hex grid #alpha=contains_sites
  scale_fill_viridis_c(option="F", direction=-1) +
  scale_alpha_manual(values=c(0.45, 1)) +
  xlab('Longitude') +
  ylab('Latitude') +
  geom_sf_label(aes(label = paste0(median_date, "BP")), label.size  = NA, alpha = 0.4, size=3.5) + #hex grid labels #label = ifelse(contains_sites==0, NA, paste0(median_date, "BP")))
  theme(panel.background = element_rect(fill = "lightblue",
                                        colour = "lightblue",
                                        size = 0.5,
                                        linetype = "solid"),
        legend.position = "none")


#Output
pdf(file=here('output','figures','sim2_arrivaltime.pdf'), width=15, height=8)
grid.arrange(modi, ncol=1, padding=0)
dev.off()


#------------
##FIGURE 7 -- Posterior distributions of arrival times
#For model (i) and (ii) select parameters a and b (i.e. start and end date of occupation in the region)
sim_a <- constants$true_a
sim_b <- constants$true_b

#Functions
extract <- function(x)
{
  tmp = do.call(rbind, x)
  tmp2 = tmp[ , grep('^a\\[',colnames(tmp))]
  qta = apply(tmp2, 2, quantile, prob=c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1))
  return(qta)
}

post.bar <- function(x, i, h, a)
{
  rect(xleft = x[2], xright = x[6], ybottom = i - h/5, ytop = i + h/5, border = NA, col = "skyblue") # 95% interval rectangle
  segments(x[2], i-h/3.5, x[2], i+h/3.5, lwd = 2, col = "skyblue") # horizontal ticks for 95%
  segments(x[6], i-h/3.5, x[6], i+h/3.5, lwd = 2, col = "skyblue") # horizontal ticks for 95%
  
  rect(xleft=x[3], xright=x[5], ybottom=i-h/3, ytop=i+h/3, border=NA, col="dodgerblue") #50% interval
  segments(x[3], i-h/2.5, x[3], i+h/2.5, lwd = 2, col = "dodgerblue")   # horizontal ticks for 50%
  segments(x[5], i-h/2.5, x[5], i+h/2.5, lwd = 2, col = "dodgerblue")   # horizontal ticks for 50%
  
  points(x[4], i, pch = 16, col = "darkblue", cex = 2) #posterior median 
  points(a, i, pch = 4, col = "darkgreen", cex = 2, lwd = 2) #simulated (true) arrival time 
}

#Plot
pdf(file=here('output','figures','sim2_posteriors.pdf'), width=10, height=15, pointsize=4)
par(mar = c(5, 5, 4, 2))   #pad space around plot
plot(NULL, xlim=c(8200, 6000), ylim=c(3, 79), xlab=paste('Arrival time (BP),', TeX('$a_k$')), ylab=paste('Area,', TeX('$k$')), cex.lab = 2, axes=F)
tmp.a = extract(out_womble_model)
iseq.a = seq(1,by=1,length.out=81)
abline(h=seq(1,by=1,length.out=81), col='lightgrey')

counter <- 1 #indexing counter
for (i in c(1:81)) #all relevant hex areas
{
  #Plot bar in area i
  post.bar(tmp.a[,i], i=iseq.a[counter], h=0.5, a= sim_a[[i]])
  counter <- counter + 1
}

axis(2, at=iseq.a, labels = paste0(c(1:81)), las=2, cex.axis=1.7)
axis(1, at = BCADtoBP(c(-6300, -6100, -5900, -5700, -5500, -5300, -5100, -4900, -4700, -4400, -4200, -4000)), labels=c('6300BC', '6100BC', '5900BC', '5700BC', '5500BC', '5300BC', '5100BC', '4900BC', '4700BC', '4400BC','4200BC', '4000BC'), tck=-0.01, cex.axis=1.7)
axis(3, at = seq(8200, 6000, -200), labels=paste0(seq(8200, 6000, -200),'BP'), tck=-0.01, cex.axis=1.7)
axis(1, at = BCADtoBP(c(-6200, -6000, -5800, -5600, -5400, -5200, -5000, -4800, -4600, -4500, -4300, -4100)), labels=NA, tck=-0.01) #Minor tick marks
axis(3, at = seq(8200, 6000, -50), labels=NA, tck=-0.01) #Minor tick marks
box()
#Legend
post.bar(c(6500,6400,6300,6200,6100,6000,5900), i=77, h=0.9, a=6050)
text(x=6200, y=78, "95% HPDI", cex=1.5)
text(x=6400, y=78,"50% HPDI", cex=1.5)
text(x=6300, y=76, "Median Posterior", cex=1.5)
text(x=6050, y=76, "Simulated value", cex=1.5)
rect(xleft=6550, xright=5920, ybottom=75, ytop=79, border="darkgrey", col=NA, lwd=2)
theme(legend.position = "none")
dev.off()


#------------
##FIGURE 8 -- Map of Wombling Boundaries (highlighting significant boundaries)

#With the tactical simulation data from hierarchical wombling model
post.model.tac_womble_nab  <- out.comb.tac_icar.model[,paste0('nabla[',1:208,']')]  %>% round()

#Extract differences in arrival times for tactical wombling model
med.model.tac_womble_nab  <- apply(post.model.tac_womble_nab, 2, median)

#Extract proportion of MCMC sample differences which are significant over a specified time difference
prop_model_tac_womble_nab  <- data.frame(x = 1:208,
                                         y = sapply(as.data.frame(post.model.tac_womble_nab), 
                                                    prop_gthan_threshold, 
                                                    threshold = 400))

#Add info to edges dataframe
edge_info.i <- edge_info %>%
  mutate(mean_gradient = med.model.tac_womble_nab, #50% quantile
         prob_BLV = prop_model_tac_womble_nab$y, #% of distribution > specified threshold
         boundary = mapply(function(a, b) {intersection <- st_intersection(hex_area_win_proj$geometry[[a]], hex_area_win_proj$geometry[[b]])
         if (st_is_empty(intersection) || st_is(intersection, "MULTILINESTRING")) return(st_linestring()) else return(intersection)},
         edge_info$region1_id,
         edge_info$region2_id, SIMPLIFY = FALSE)) #shared boundary between two subareas

#Create nodes
nodes <- st_coordinates(hex_area_win_proj$area_center)

#Create boundary segments
boundaries <- st_sf(prob_BLV = edge_info.i$prob_BLV,
                    geometry = st_sfc(edge_info.i$boundary)) #lapply(edge_info.i$boundary[[a]], st_coordinates(a))
st_crs(boundaries) <- 3035  # Set CRS for correct Europe projection

#Plot
pdf(file=here('output','figures','sim2_womble.pdf'))
ggplot(data = median_hex_dates_mod.i) +
  geom_sf(data = st_buffer(sampling_win_proj, 40000), fill = "grey80", color = "grey40") + #sampling window with coastal buffer
  geom_sf(aes(alpha=0.01), color = "grey60") + scale_alpha(range = c(0, 1)) + #hex grid 
  # geom_segment(data=edge_info.i$boundary, aes(#x= region1_x, y= region1_y, xend= region2_x, yend= region2_y, alpha= prob_BLV), color="red", size=2) +
  geom_sf(data = boundaries, lwd=3, aes(alpha=prob_BLV), color = "red") +
  geom_sf(data = hex_area_win_proj$area_center, size=2, alpha=1, color = "grey40") + #hex-centers
  scale_alpha_continuous(range = c(0, 1)) +  # Use for continuous alpha values
  xlab('Longitude') +
  ylab('Latitude') +
  ggtitle(paste0('c = 400', ' years')) +
  theme(panel.background = element_rect(fill = "lightblue",
                                        colour = "lightblue",
                                        size = 0.5,
                                        linetype = "solid"),
        legend.position = "none")
dev.off()

#===============================================================================
#SIMULATION 3: TACTICAL ICAR SIMULATION (with calibrated radiocarbon dates not the plateau)

#Load data
load(here('data','tactical_sim_icar_withoutplat.RData')) #simulated data
load(here('output', 'Womblemodel_tactical_withoutplat_errors.RData')) #inferred data

#------------
#Hex areas with and without out sites
Hex_with_sites <- unique(siteInfo$area_id)
Hex_without_sites <- which(rep(1:81) %!in% Hex_with_sites)

#------------
#Chain health check
traceplot(out_womble_model[,'a[20]'], main=TeX('$a$'), smooth=TRUE) #region 20 as an example
traceplot(out_womble_model[,'theta[100]'], main=TeX('$theta$'), smooth=TRUE)

#------------
#Calculate accuracy and precision
ci_95 = credible_interval(out_womble_model, 0.95)
sim_a <- constants$true_a
sim3_model_accuracy = accuracy(sim_a, ci_95)
sim3_model_precision = precision(sim_a, ci_95)

#-------------------------------------------------------------------------------
##FIGURE 9 -- Map of sites and sampling window
#Plot
site_map <- ggplot(data = hex_area_win_proj) +
  geom_sf(data = sampling_win_proj, color = "grey50") +  # sampling window border
  geom_sf() +
  geom_sf(data = sites_sf, size = 2, alpha = 0.5) +  # sites
  geom_sf_label(aes(label = area_ID)) +                     # area labels
  theme(
    panel.background = element_rect(fill = "lightblue", colour = "lightblue"),
    legend.title = element_blank(),
    legend.position = "bottom")

#Output
pdf(file=here('output','figures','sim3_sites.pdf'), width=15, height=8)
grid.arrange(site_map, ncol=1, padding=0)
dev.off()

#------------
##FIGURE 10 -- Map of Inferred arrival times
out.comb.tac_icar.model  <- do.call(rbind, out_womble_model)
post.a.model.i  <- out.comb.tac_icar.model[,paste0('a[',1:81,']')]  %>% round() 
med.model.i  <- apply(post.a.model.i, 2, median) #Extract arrival times for tactical icar model

median_hex_dates_mod.i <- hex_area_win_proj %>% 
  filter(area_ID %in% 1:81) %>% 
  mutate(median_date = med.model.i,
         contains_sites = as.factor(case_when(area_ID %in% Hex_with_sites ~ 1, area_ID %in% Hex_without_sites ~ 0))) 

#Plot
modi <- ggplot(data = median_hex_dates_mod.i) +
  geom_sf(data = st_buffer(sampling_win_proj, 40000), aes(color = "grey50")) + #sampling window with coastal buffer
  geom_sf(aes(fill = median_date)) + #hex grid #alpha=contains_sites
  scale_fill_viridis_c(option="F", direction=-1) +
  scale_alpha_manual(values=c(0.45, 1)) +
  xlab('Longitude') +
  ylab('Latitude') +
  geom_sf_label(aes(label = paste0(median_date, "BP")), label.size  = NA, alpha = 0.4, size=3.5) + #hex grid labels #label = ifelse(contains_sites==0, NA, paste0(median_date, "BP")))
  theme(panel.background = element_rect(fill = "lightblue",
                                        colour = "lightblue",
                                        size = 0.5,
                                        linetype = "solid"),
        legend.position = "none")


#Output
pdf(file=here('output','figures','sim3_arrivaltime.pdf'), width=15, height=8)
grid.arrange(modi, ncol=1, padding=0)
dev.off()

#------------
##FIGURE 11 -- Posterior distributions of arrival times
#For model (i) and (ii) select parameters a and b (i.e. start and end date of occupation in the region)
sim_a <- constants$true_a
sim_b <- constants$true_b

#Functions
extract <- function(x)
{
  tmp = do.call(rbind, x)
  tmp2 = tmp[ , grep('^a\\[',colnames(tmp))]
  qta = apply(tmp2, 2, quantile, prob=c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1))
  return(qta)
}

post.bar <- function(x, i, h, a)
{
  rect(xleft = x[2], xright = x[6], ybottom = i - h/5, ytop = i + h/5, border = NA, col = "skyblue") # 95% interval rectangle
  segments(x[2], i-h/3.5, x[2], i+h/3.5, lwd = 2, col = "skyblue") # horizontal ticks for 95%
  segments(x[6], i-h/3.5, x[6], i+h/3.5, lwd = 2, col = "skyblue") # horizontal ticks for 95%
  
  rect(xleft=x[3], xright=x[5], ybottom=i-h/3, ytop=i+h/3, border=NA, col="dodgerblue") #50% interval
  segments(x[3], i-h/2.5, x[3], i+h/2.5, lwd = 2, col = "dodgerblue")   # horizontal ticks for 50%
  segments(x[5], i-h/2.5, x[5], i+h/2.5, lwd = 2, col = "dodgerblue")   # horizontal ticks for 50%
  
  points(x[4], i, pch = 16, col = "darkblue", cex = 2) #posterior median 
  points(a, i, pch = 4, col = "darkgreen", cex = 2, lwd = 2) #simulated (true) arrival time 
}

#Plot
pdf(file=here('output','figures','sim3_posteriors.pdf'), width=10, height=15, pointsize=4)
par(mar = c(5, 5, 4, 2))   #pad space around plot
plot(NULL, xlim=c(6800, 4000), ylim=c(3, 79), xlab=paste('Arrival time (BP),', TeX('$a_k$')), ylab=paste('Area,', TeX('$k$')), cex.lab = 2, axes=F)
tmp.a = extract(out_womble_model)
iseq.a = seq(1,by=1,length.out=81)
abline(h=seq(1,by=1,length.out=81), col='lightgrey')

counter <- 1 #indexing counter
for (i in c(1:81)) #all relevant hex areas
{
  #Plot bar in area i
  post.bar(tmp.a[,i], i=iseq.a[counter], h=0.5, a= sim_a[[i]])
  counter <- counter + 1
}

axis(2, at=iseq.a, labels = paste0(c(1:81)), las=2, cex.axis=1.7)
axis(1, at = BCADtoBP(c(-4900, -4700, -4500, -4300, -4100, -3900, -3700, -3500, -3300, -3100, -2900, -2700, -2400, -2200, -2000)), labels=c('4900BC','4700BC','4500BC','4300BC', '4100BC', '3900BC', '3700BC', '3500BC', '3300BC', '3100BC', '2900BC', '2700BC', '2400BC','2200BC', '2000BC'), tck=-0.01, cex.axis=1.7)
axis(3, at = seq(6800, 4000, -200), labels=paste0(seq(6800, 4000, -200),'BP'), tck=-0.01, cex.axis=1.7)
axis(1, at = BCADtoBP(c(-4800, -4600, -4400, -4200, -4000, -3800, -3600, -3400, -3200, -3000, -2800, -2600, -2500, -2300, -2100)), labels=NA, tck=-0.01) #Minor tick marks
axis(3, at = seq(6800, 4000, -50), labels=NA, tck=-0.01) #Minor tick marks
box()
#Legend
post.bar(c(6900,6800,6600,6500,6400,6200,6100), i=77, h=0.9, a=6750)
text(x=6550, y=78, "95% HPDI", cex=1.5)
text(x=6300, y=78,"50% HPDI", cex=1.5)
text(x=6350, y=76, "Median Posterior", cex=1.5)
text(x=6700, y=76, "Simulated value", cex=1.5)
rect(xleft=6850, xright=6150, ybottom=75, ytop=79, border="darkgrey", col=NA, lwd=2)
theme(legend.position = "none")
dev.off()

#===============================================================================
#SIMULATION 4: TACTICAL ICAR SIMULATION (with uncalibrated radiocarbon dates not in the plateau)

#Load data
load(here('data','tactical_sim_icar_withoutplat.RData')) #simulated data
load(here('output', 'Womblemodel_tactical_withoutplat_noerror.RData')) #inferred data

#------------
#Hex areas with and without out sites
Hex_with_sites <- unique(siteInfo$area_id)
Hex_without_sites <- which(rep(1:81) %!in% Hex_with_sites)

#------------
#Chain health check
traceplot(out_womble_model[,'a[20]'], main=TeX('$a$'), smooth=TRUE) #region 20 as an example
traceplot(out_womble_model[,'theta[100]'], main=TeX('$theta$'), smooth=TRUE)

#------------
#Calculate accuracy and precision
ci_95 = credible_interval(out_womble_model, 0.95)
sim_a <- constants$true_a
sim4_model_accuracy = accuracy(sim_a, ci_95)
sim4_model_precision = precision(sim_a, ci_95)

#-------------------------------------------------------------------------------
##FIGURE 12 -- Map of sites and sampling window
#Plot
site_map <- ggplot(data = hex_area_win_proj) +
  geom_sf(data = sampling_win_proj, color = "grey50") +  # sampling window border
  geom_sf() +
  geom_sf(data = sites_sf, size = 2, alpha = 0.5) +  # sites
  geom_sf_label(aes(label = area_ID)) +                     # area labels
  theme(
    panel.background = element_rect(fill = "lightblue", colour = "lightblue"),
    legend.title = element_blank(),
    legend.position = "bottom")

#Output
pdf(file=here('output','figures','sim4_sites.pdf'), width=15, height=8)
grid.arrange(site_map, ncol=1, padding=0)
dev.off()

#------------
##FIGURE 13 -- Map of Inferred arrival times
out.comb.tac_icar.model  <- do.call(rbind, out_womble_model)
post.a.model.i  <- out.comb.tac_icar.model[,paste0('a[',1:81,']')]  %>% round() 
med.model.i  <- apply(post.a.model.i, 2, median) #Extract arrival times for tactical icar model

median_hex_dates_mod.i <- hex_area_win_proj %>% 
  filter(area_ID %in% 1:81) %>% 
  mutate(median_date = med.model.i,
         contains_sites = as.factor(case_when(area_ID %in% Hex_with_sites ~ 1, area_ID %in% Hex_without_sites ~ 0))) 

#Plot
modi <- ggplot(data = median_hex_dates_mod.i) +
  geom_sf(data = st_buffer(sampling_win_proj, 40000), aes(color = "grey50")) + #sampling window with coastal buffer
  geom_sf(aes(fill = median_date)) + #hex grid #alpha=contains_sites
  scale_fill_viridis_c(option="F", direction=-1) +
  scale_alpha_manual(values=c(0.45, 1)) +
  xlab('Longitude') +
  ylab('Latitude') +
  geom_sf_label(aes(label = paste0(median_date, "BP")), label.size  = NA, alpha = 0.4, size=3.5) + #hex grid labels #label = ifelse(contains_sites==0, NA, paste0(median_date, "BP")))
  theme(panel.background = element_rect(fill = "lightblue",
                                        colour = "lightblue",
                                        size = 0.5,
                                        linetype = "solid"),
        legend.position = "none")


#Output
pdf(file=here('output','figures','sim4_arrivaltime.pdf'), width=15, height=8)
grid.arrange(modi, ncol=1, padding=0)
dev.off()

#------------
##FIGURE 14 -- Posterior distributions of arrival times
#For model (i) and (ii) select parameters a and b (i.e. start and end date of occupation in the region)
sim_a <- constants$true_a
sim_b <- constants$true_b

#Functions
extract <- function(x)
{
  tmp = do.call(rbind, x)
  tmp2 = tmp[ , grep('^a\\[',colnames(tmp))]
  qta = apply(tmp2, 2, quantile, prob=c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1))
  return(qta)
}

post.bar <- function(x, i, h, a)
{
  rect(xleft = x[2], xright = x[6], ybottom = i - h/5, ytop = i + h/5, border = NA, col = "skyblue") # 95% interval rectangle
  segments(x[2], i-h/3.5, x[2], i+h/3.5, lwd = 2, col = "skyblue") # horizontal ticks for 95%
  segments(x[6], i-h/3.5, x[6], i+h/3.5, lwd = 2, col = "skyblue") # horizontal ticks for 95%
  
  rect(xleft=x[3], xright=x[5], ybottom=i-h/3, ytop=i+h/3, border=NA, col="dodgerblue") #50% interval
  segments(x[3], i-h/2.5, x[3], i+h/2.5, lwd = 2, col = "dodgerblue")   # horizontal ticks for 50%
  segments(x[5], i-h/2.5, x[5], i+h/2.5, lwd = 2, col = "dodgerblue")   # horizontal ticks for 50%
  
  points(x[4], i, pch = 16, col = "darkblue", cex = 2) #posterior median 
  points(a, i, pch = 4, col = "darkgreen", cex = 2, lwd = 2) #simulated (true) arrival time 
}

#Plot
pdf(file=here('output','figures','sim4_posteriors.pdf'), width=10, height=15, pointsize=4)
par(mar = c(5, 5, 4, 2))   #pad space around plot
plot(NULL, xlim=c(6800, 4000), ylim=c(3, 79), xlab=paste('Arrival time (BP),', TeX('$a_k$')), ylab=paste('Area,', TeX('$k$')), cex.lab = 2, axes=F)
tmp.a = extract(out_womble_model)
iseq.a = seq(1,by=1,length.out=81)
abline(h=seq(1,by=1,length.out=81), col='lightgrey')

counter <- 1 #indexing counter
for (i in c(1:81)) #all relevant hex areas
{
  #Plot bar in area i
  post.bar(tmp.a[,i], i=iseq.a[counter], h=0.5, a= sim_a[[i]])
  counter <- counter + 1
}

axis(2, at=iseq.a, labels = paste0(c(1:81)), las=2, cex.axis=1.7)
axis(1, at = BCADtoBP(c(-4900, -4700, -4500, -4300, -4100, -3900, -3700, -3500, -3300, -3100, -2900, -2700, -2400, -2200, -2000)), labels=c('4900BC','4700BC','4500BC','4300BC', '4100BC', '3900BC', '3700BC', '3500BC', '3300BC', '3100BC', '2900BC', '2700BC', '2400BC','2200BC', '2000BC'), tck=-0.01, cex.axis=1.7)
axis(3, at = seq(6800, 4000, -200), labels=paste0(seq(6800, 4000, -200),'BP'), tck=-0.01, cex.axis=1.7)
axis(1, at = BCADtoBP(c(-4800, -4600, -4400, -4200, -4000, -3800, -3600, -3400, -3200, -3000, -2800, -2600, -2500, -2300, -2100)), labels=NA, tck=-0.01) #Minor tick marks
axis(3, at = seq(6800, 4000, -50), labels=NA, tck=-0.01) #Minor tick marks
box()
#Legend
post.bar(c(6900,6800,6600,6500,6400,6200,6100), i=77, h=0.9, a=6750)
text(x=6550, y=78, "95% HPDI", cex=1.5)
text(x=6300, y=78,"50% HPDI", cex=1.5)
text(x=6350, y=76, "Median Posterior", cex=1.5)
text(x=6700, y=76, "Simulated value", cex=1.5)
rect(xleft=6850, xright=6150, ybottom=75, ytop=79, border="darkgrey", col=NA, lwd=2)
theme(legend.position = "none")
dev.off()


#===============================================================================
#SIMULATION 5: TACTICAL ICAR SIMULATION (with calibrated radiocarbon in the plateau)

#Load data
load(here('data','tactical_sim_icar_withinplat.RData')) #simulated data
load(here('output', 'Womblemodel_tactical_withplat_errors.RData')) #inferred data

#------------
#Hex areas with and without out sites
Hex_with_sites <- unique(siteInfo$area_id)
Hex_without_sites <- which(rep(1:81) %!in% Hex_with_sites)

#------------
#Chain health check
traceplot(out_womble_model[,'a[20]'], main=TeX('$a$'), smooth=TRUE) #region 20 as an example
traceplot(out_womble_model[,'theta[100]'], main=TeX('$theta$'), smooth=TRUE)

#------------
#Calculate accuracy and precision
ci_95 = credible_interval(out_womble_model, 0.95)
sim_a <- constants$true_a
sim5_model_accuracy = accuracy(sim_a, ci_95)
sim5_model_precision = precision(sim_a, ci_95)

#-------------------------------------------------------------------------------
##FIGURE 15 -- Map of sites and sampling window
#Plot
site_map <- ggplot(data = hex_area_win_proj) +
  geom_sf(data = sampling_win_proj, color = "grey50") +  # sampling window border
  geom_sf() +
  geom_sf(data = sites_sf, size = 2, alpha = 0.5) +  # sites
  geom_sf_label(aes(label = area_ID)) +                     # area labels
  theme(
    panel.background = element_rect(fill = "lightblue", colour = "lightblue"),
    legend.title = element_blank(),
    legend.position = "bottom")

#Output
pdf(file=here('output','figures','sim5_sites.pdf'), width=15, height=8)
grid.arrange(site_map, ncol=1, padding=0)
dev.off()

#------------
##FIGURE 16 -- Map of Inferred arrival times
out.comb.tac_icar.model  <- do.call(rbind, out_womble_model)
post.a.model.i  <- out.comb.tac_icar.model[,paste0('a[',1:81,']')]  %>% round() 
med.model.i  <- apply(post.a.model.i, 2, median) #Extract arrival times for tactical icar model

median_hex_dates_mod.i <- hex_area_win_proj %>% 
  filter(area_ID %in% 1:81) %>% 
  mutate(median_date = med.model.i,
         contains_sites = as.factor(case_when(area_ID %in% Hex_with_sites ~ 1, area_ID %in% Hex_without_sites ~ 0))) 

#Plot
modi <- ggplot(data = median_hex_dates_mod.i) +
  geom_sf(data = st_buffer(sampling_win_proj, 40000), aes(color = "grey50")) + #sampling window with coastal buffer
  geom_sf(aes(fill = median_date)) + #hex grid #alpha=contains_sites
  scale_fill_viridis_c(option="F", direction=-1) +
  scale_alpha_manual(values=c(0.45, 1)) +
  xlab('Longitude') +
  ylab('Latitude') +
  geom_sf_label(aes(label = paste0(median_date, "BP")), label.size  = NA, alpha = 0.4, size=3.5) + #hex grid labels #label = ifelse(contains_sites==0, NA, paste0(median_date, "BP")))
  theme(panel.background = element_rect(fill = "lightblue",
                                        colour = "lightblue",
                                        size = 0.5,
                                        linetype = "solid"),
        legend.position = "none")


#Output
pdf(file=here('output','figures','sim5_arrivaltime.pdf'), width=15, height=8)
grid.arrange(modi, ncol=1, padding=0)
dev.off()

#------------
##FIGURE 17 -- Posterior distributions of arrival times
#For model (i) and (ii) select parameters a and b (i.e. start and end date of occupation in the region)
sim_a <- constants$true_a
sim_b <- constants$true_b

#Functions
extract <- function(x)
{
  tmp = do.call(rbind, x)
  tmp2 = tmp[ , grep('^a\\[',colnames(tmp))]
  qta = apply(tmp2, 2, quantile, prob=c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1))
  return(qta)
}

post.bar <- function(x, i, h, a)
{
  rect(xleft = x[2], xright = x[6], ybottom = i - h/5, ytop = i + h/5, border = NA, col = "skyblue") # 95% interval rectangle
  segments(x[2], i-h/3.5, x[2], i+h/3.5, lwd = 2, col = "skyblue") # horizontal ticks for 95%
  segments(x[6], i-h/3.5, x[6], i+h/3.5, lwd = 2, col = "skyblue") # horizontal ticks for 95%
  
  rect(xleft=x[3], xright=x[5], ybottom=i-h/3, ytop=i+h/3, border=NA, col="dodgerblue") #50% interval
  segments(x[3], i-h/2.5, x[3], i+h/2.5, lwd = 2, col = "dodgerblue")   # horizontal ticks for 50%
  segments(x[5], i-h/2.5, x[5], i+h/2.5, lwd = 2, col = "dodgerblue")   # horizontal ticks for 50%
  
  points(x[4], i, pch = 16, col = "darkblue", cex = 2) #posterior median 
  points(a, i, pch = 4, col = "darkgreen", cex = 2, lwd = 2) #simulated (true) arrival time 
}

#Plot
pdf(file=here('output','figures','sim5_posteriors.pdf'), width=10, height=15, pointsize=4)
par(mar = c(5, 5, 4, 2))   #pad space around plot
plot(NULL, xlim=c(3800, 1000), ylim=c(3, 79), xlab=paste('Arrival time (BP),', TeX('$a_k$')), ylab=paste('Area,', TeX('$k$')), cex.lab = 2, axes=F)
tmp.a = extract(out_womble_model)
iseq.a = seq(1,by=1,length.out=81)
abline(h=seq(1,by=1,length.out=81), col='lightgrey')

counter <- 1 #indexing counter
for (i in c(1:81)) #all relevant hex areas
{
  #Plot bar in area i
  post.bar(tmp.a[,i], i=iseq.a[counter], h=0.5, a= sim_a[[i]])
  counter <- counter + 1
}

axis(2, at=iseq.a, labels = paste0(c(1:81)), las=2, cex.axis=1.7)
axis(1, at = BCADtoBP(c(-1900, -1700, -1500, -1300, -1100, -900, -700, -500, -300, -100, 100, 300, 500, 700, 900)), labels=c('1900BC','1700BC','1500BC','1300BC', '1100BC', '900BC', '700BC', '500BC', '300BC', '100BC', '100AD', '300AD', '500AD','700AD', '900AD'), tck=-0.01, cex.axis=1.7)
axis(3, at = seq(3800, 1000, -200), labels=paste0(seq(3800, 1000, -200),'BP'), tck=-0.01, cex.axis=1.7)
axis(1, at = BCADtoBP(c(-1800, -1600, -1400, -1200, -1000, -800, -600, -400, -200, 1, 200, 400, 600, 800, 1000)), labels=NA, tck=-0.01) #Minor tick marks
axis(3, at = seq(3800, 1000, -50), labels=NA, tck=-0.01) #Minor tick marks
box()
#Legend
post.bar(c(1600,1500,1400,1300,1200,1100,1000), i=77, h=0.9, a=1450)
text(x=1350, y=78, "95% HPDI", cex=1.5)
text(x=1550, y=78,"50% HPDI", cex=1.5)
text(x=1300, y=76, "Median Posterior", cex=1.5)
text(x=1450, y=76, "Simulated value", cex=1.5)
rect(xleft=1650, xright=950, ybottom=75, ytop=79, border="darkgrey", col=NA, lwd=2)
theme(legend.position = "none")
dev.off()

#===============================================================================
#SIMULATION 6: TACTICAL ICAR SIMULATION (with uncalibrated radiocarbon within the plateau)

#Load data
load(here('data','tactical_sim_icar_withinplat.RData')) #simulated data
load(here('output', 'Womblemodel_tactical_withplat_noerrors.RData')) #inferred data

#------------
#Hex areas with and without out sites
Hex_with_sites <- unique(siteInfo$area_id)
Hex_without_sites <- which(rep(1:81) %!in% Hex_with_sites)

#------------
#Chain health check
traceplot(out_womble_model[,'a[20]'], main=TeX('$a$'), smooth=TRUE) #region 20 as an example
traceplot(out_womble_model[,'theta[100]'], main=TeX('$theta$'), smooth=TRUE)

#------------
#Calculate accuracy and precision
ci_95 = credible_interval(out_womble_model, 0.95)
sim_a <- constants$true_a
sim6_model_accuracy = accuracy(sim_a, ci_95)
sim6_model_precision = precision(sim_a, ci_95)

#-------------------------------------------------------------------------------
##FIGURE 18 -- Map of sites and sampling window
#Plot
site_map <- ggplot(data = hex_area_win_proj) +
  geom_sf(data = sampling_win_proj, color = "grey50") +  # sampling window border
  geom_sf() +
  geom_sf(data = sites_sf, size = 2, alpha = 0.5) +  # sites
  geom_sf_label(aes(label = area_ID)) +                     # area labels
  theme(
    panel.background = element_rect(fill = "lightblue", colour = "lightblue"),
    legend.title = element_blank(),
    legend.position = "bottom")

#Output
pdf(file=here('output','figures','sim6_sites.pdf'), width=15, height=8)
grid.arrange(site_map, ncol=1, padding=0)
dev.off()

#------------
##FIGURE 19 -- Map of Inferred arrival times
out.comb.tac_icar.model  <- do.call(rbind, out_womble_model)
post.a.model.i  <- out.comb.tac_icar.model[,paste0('a[',1:81,']')]  %>% round() 
med.model.i  <- apply(post.a.model.i, 2, median) #Extract arrival times for tactical icar model

median_hex_dates_mod.i <- hex_area_win_proj %>% 
  filter(area_ID %in% 1:81) %>% 
  mutate(median_date = med.model.i,
         contains_sites = as.factor(case_when(area_ID %in% Hex_with_sites ~ 1, area_ID %in% Hex_without_sites ~ 0))) 

#Plot
modi <- ggplot(data = median_hex_dates_mod.i) +
  geom_sf(data = st_buffer(sampling_win_proj, 40000), aes(color = "grey50")) + #sampling window with coastal buffer
  geom_sf(aes(fill = median_date)) + #hex grid #alpha=contains_sites
  scale_fill_viridis_c(option="F", direction=-1) +
  scale_alpha_manual(values=c(0.45, 1)) +
  xlab('Longitude') +
  ylab('Latitude') +
  geom_sf_label(aes(label = paste0(median_date, "BP")), label.size  = NA, alpha = 0.4, size=3.5) + #hex grid labels #label = ifelse(contains_sites==0, NA, paste0(median_date, "BP")))
  theme(panel.background = element_rect(fill = "lightblue",
                                        colour = "lightblue",
                                        size = 0.5,
                                        linetype = "solid"),
        legend.position = "none")


#Output
pdf(file=here('output','figures','sim6_arrivaltime.pdf'), width=15, height=8)
grid.arrange(modi, ncol=1, padding=0)
dev.off()

#------------
##FIGURE 20 -- Posterior distributions of arrival times
#For model (i) and (ii) select parameters a and b (i.e. start and end date of occupation in the region)
sim_a <- constants$true_a
sim_b <- constants$true_b

#Functions
extract <- function(x)
{
  tmp = do.call(rbind, x)
  tmp2 = tmp[ , grep('^a\\[',colnames(tmp))]
  qta = apply(tmp2, 2, quantile, prob=c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1))
  return(qta)
}

post.bar <- function(x, i, h, a)
{
  rect(xleft = x[2], xright = x[6], ybottom = i - h/5, ytop = i + h/5, border = NA, col = "skyblue") # 95% interval rectangle
  segments(x[2], i-h/3.5, x[2], i+h/3.5, lwd = 2, col = "skyblue") # horizontal ticks for 95%
  segments(x[6], i-h/3.5, x[6], i+h/3.5, lwd = 2, col = "skyblue") # horizontal ticks for 95%
  
  rect(xleft=x[3], xright=x[5], ybottom=i-h/3, ytop=i+h/3, border=NA, col="dodgerblue") #50% interval
  segments(x[3], i-h/2.5, x[3], i+h/2.5, lwd = 2, col = "dodgerblue")   # horizontal ticks for 50%
  segments(x[5], i-h/2.5, x[5], i+h/2.5, lwd = 2, col = "dodgerblue")   # horizontal ticks for 50%
  
  points(x[4], i, pch = 16, col = "darkblue", cex = 2) #posterior median 
  points(a, i, pch = 4, col = "darkgreen", cex = 2, lwd = 2) #simulated (true) arrival time 
}

#Plot
pdf(file=here('output','figures','sim6_posteriors.pdf'), width=10, height=15, pointsize=4)
par(mar = c(5, 5, 4, 2))   #pad space around plot
plot(NULL, xlim=c(3800, 1000), ylim=c(3, 79), xlab=paste('Arrival time (BP),', TeX('$a_k$')), ylab=paste('Area,', TeX('$k$')), cex.lab = 2, axes=F)
tmp.a = extract(out_womble_model)
iseq.a = seq(1,by=1,length.out=81)
abline(h=seq(1,by=1,length.out=81), col='lightgrey')

counter <- 1 #indexing counter
for (i in c(1:81)) #all relevant hex areas
{
  #Plot bar in area i
  post.bar(tmp.a[,i], i=iseq.a[counter], h=0.5, a= sim_a[[i]])
  counter <- counter + 1
}

axis(2, at=iseq.a, labels = paste0(c(1:81)), las=2, cex.axis=1.7)
axis(1, at = BCADtoBP(c(-1900, -1700, -1500, -1300, -1100, -900, -700, -500, -300, -100, 100, 300, 500, 700, 900)), labels=c('1900BC','1700BC','1500BC','1300BC', '1100BC', '900BC', '700BC', '500BC', '300BC', '100BC', '100AD', '300AD', '500AD','700AD', '900AD'), tck=-0.01, cex.axis=1.7)
axis(3, at = seq(3800, 1000, -200), labels=paste0(seq(3800, 1000, -200),'BP'), tck=-0.01, cex.axis=1.7)
axis(1, at = BCADtoBP(c(-1800, -1600, -1400, -1200, -1000, -800, -600, -400, -200, 1, 200, 400, 600, 800, 1000)), labels=NA, tck=-0.01) #Minor tick marks
axis(3, at = seq(3800, 1000, -50), labels=NA, tck=-0.01) #Minor tick marks
box()
#Legend
post.bar(c(1600,1500,1400,1300,1200,1100,1000), i=77, h=0.9, a=1450)
text(x=1350, y=78, "95% HPDI", cex=1.5)
text(x=1550, y=78,"50% HPDI", cex=1.5)
text(x=1300, y=76, "Median Posterior", cex=1.5)
text(x=1450, y=76, "Simulated value", cex=1.5)
rect(xleft=1650, xright=950, ybottom=75, ytop=79, border="darkgrey", col=NA, lwd=2)
theme(legend.position = "none")
dev.off()

