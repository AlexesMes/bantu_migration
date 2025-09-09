# Load Libraries and spatial data ----
library(here)
library(dplyr)
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

`%!in%` <- Negate(`%in%`)

#===============================================================================
## Load Data ----

source(here('src', 'grad_funcs.R'))
source(here('src', 'hex_areas.R'))

load(here('data','tactical_sim_woa.RData')) #simulated data
load(here('output', 'Womblemodel_tactical_woa.RData')) #inferred data

load(here('data','trig.RData')) #Load nodes and edges between hex area centroids
load(here('data','sample_window.RData')) 


#------------
#Hex areas with and without out sites
Hex_with_sites <- unique(siteInfo$area_id)
Hex_without_sites <- which(rep(1:81) %!in% Hex_with_sites)

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
##Simulation 1 -- WoA

##FIGURE 2 -- Map of sites and sampling window
#Plot
site_map <- ggplot(data = hex_area_win_proj) +
    geom_sf(data = sampling_win_proj, color = "grey50") +  # sampling window border
    geom_sf() +
    geom_sf(data = as(sites, 'sf'), size = 2, alpha = 0.5) +  # sites
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
out.comb.tac_icar.model  <- do.call(rbind, mcmc.samplesW)
post.a.model.i  <- out.comb.tac_icar.model[,paste0('a[',1:143,']')]  %>% round() 
med.model.i  <- apply(post.a.model.i, 2, median) #Extract arrival times for tactical icar model

median_hex_dates_mod.i <- hex_area_win_proj %>% 
  filter(area_ID %in% 1:143) %>% 
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

post.bar <- function(x, i, h, a, col)
{
  rect(xleft=x[2], xright=x[6], ybottom=i-h/5, ytop=i+h/5, border=NA, col=col)
  rect(xleft=x[3], xright=x[5], ybottom=i-h/3, ytop=i+h/3, border=NA, col=col)
  lines(c(x[4],x[4]), c(i-h/2,i+h/2), lwd=3, col='blue')
  lines(c(a,a), c(i-h/2,i+h/2), lwd=3, col='darkgreen') #so=imulated (true) arrival time
}

#Plot
pdf(file=here('output','figures','sim1_posteriors.pdf'), width=10, height=18, pointsize=4)
plot(NULL, xlim=c(2700, 1000), ylim=c(3,95), xlab=paste('Arrival time,', TeX('$a_k$')), ylab=paste('Area,', TeX('$k$')), cex.lab = 2, axes=F)
tmp.a = extract(mcmc.samplesW)
iseq.a = seq(1,by=2,length.out=143)
abline(h=seq(2,by=2,length.out=143), col='darkgrey',lty=2)

counter <- 1 #indexing counter
for (i in c(1:143)) #all relevant hex areas
{
  #Plot bar in area i
  post.bar(tmp.a[,i], i=iseq.a[counter], h=0.9, a= sim_a[[i]], col='lightblue')
  counter <- counter + 1
}

axis(2, at=iseq.a+0.5, labels = paste0(c(1:143)), las=2, cex.axis=1.7)
axis(1, at = BCADtoBP(c(-900, -700, -500, -300, -100, 100, 300, 500, 700, 900)), labels=c('900BC','700BC', '500BC', '300BC', '100BC', '100AD', '300AD', '500AD', '700AD', '900AD'), tck=-0.01, cex.axis=1.7)
axis(3, at = seq(2700, 1000, -200), labels=paste0(seq(2700, 1000, -200),'BP'), tck=-0.01, cex.axis=1.7)
axis(1, at = BCADtoBP(c(-600, -400, -200, 1, 200, 400, 600)), labels=NA, tck=-0.01) #Minor tick marks
axis(3, at = seq(2700, 1000, -50), labels=NA, tck=-0.01) #Minor tick marks
box()

post.bar(c(1500,1400,1300,1200,1100,1000,900), i=96, h=0.9, a=1050, col='lightgrey')
#arrows(x0=1000, x1=500, y0=0.3, y1=90, angle = 90, code = 3, length = 0.01)
#arrows(x0=900, x1=600, y0=0.8, y1=90, angle = 90, code = 3, length = 0.01)
text(x=1200, y=97.5, "95% HPDI", cex=1.5)
text(x=1400, y=97.5,"50% HPDI", cex=1.5)
text(x=1250, y=95, "Median Posterior", cex=1.5)
text(x=1050, y=95, "Simulated value", cex=1.5)
#lines(x=c(550,500), y=c(2.35, 2.28))
theme(legend.position = "none")
dev.off()


#------------
##FIGURE 5 -- Map of Wombling Boundaries (highlighting significant boundaries)

#With the tactical simulation data from hierarchical wombling model
post.model.tac_womble_nab  <- out.comb.tac_icar.model[,paste0('nabla[',1:383,']')]  %>% round()

#Extract differences in arrival times for tactical wombling model
med.model.tac_womble_nab  <- apply(post.model.tac_womble_nab, 2, median)

#Extract proportion of MCMC sample differences which are significant over a specified time difference
prop_model_tac_womble_nab  <- data.frame(x = 1:383,
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
         edge_info$region2_id)) #shared boundary between two subareas 

#Create nodes
nodes <- st_coordinates(hex_area_win_proj$area_center)

#Create boundary segments
boundaries <- st_sf(prob_BLV = edge_info.i$prob_BLV,
                    geometry = st_sfc(edge_info.i$boundary)) #lapply(edge_info.i$boundary[[a]], st_coordinates(a))
st_crs(boundaries) <- 3035  # Set CRS for correct Europe projection

#Plot
pdf(file=here('output','figures','sim1_womble.pdf'))
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


#-------------------------------------------------------------------------------
##Simulation 2 -- ICAR with 2 covariates with error

##FIGURE 6 -- Map of sites and sampling window
#Plot
site_map <- ggplot(data = hex_area_win_proj) +
  geom_sf(data = sampling_win_proj, color = "grey50") +  # sampling window border
  geom_sf() +
  geom_sf(data = as(sites, 'sf'), size = 2, alpha = 0.5) +  # sites
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
##FIGURE 7 -- Map of Inferred arrival times
out.comb.tac_icar.model  <- do.call(rbind, mcmc.samplesW)
post.a.model.i  <- out.comb.tac_icar.model[,paste0('a[',1:143,']')]  %>% round() 
med.model.i  <- apply(post.a.model.i, 2, median) #Extract arrival times for tactical icar model

median_hex_dates_mod.i <- hex_area_win_proj %>% 
  filter(area_ID %in% 1:143) %>% 
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
##FIGURE 8 -- Posterior distributions of arrival times
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

post.bar <- function(x, i, h, a, col)
{
  rect(xleft=x[2], xright=x[6], ybottom=i-h/5, ytop=i+h/5, border=NA, col=col)
  rect(xleft=x[3], xright=x[5], ybottom=i-h/3, ytop=i+h/3, border=NA, col=col)
  lines(c(x[4],x[4]), c(i-h/2,i+h/2), lwd=3, col='blue')
  lines(c(a,a), c(i-h/2,i+h/2), lwd=3, col='darkgreen') #so=imulated (true) arrival time
}

#Plot
pdf(file=here('output','figures','sim2_posteriors.pdf'), width=10, height=18, pointsize=4)
plot(NULL, xlim=c(2700, 1000), ylim=c(3,95), xlab=paste('Arrival time,', TeX('$a_k$')), ylab=paste('Area,', TeX('$k$')), cex.lab = 2, axes=F)
tmp.a = extract(mcmc.samplesW)
iseq.a = seq(1,by=2,length.out=143)
abline(h=seq(2,by=2,length.out=143), col='darkgrey',lty=2)

counter <- 1 #indexing counter
for (i in c(1:143)) #all relevant hex areas
{
  #Plot bar in area i
  post.bar(tmp.a[,i], i=iseq.a[counter], h=0.9, a= sim_a[[i]], col='lightblue')
  counter <- counter + 1
}

axis(2, at=iseq.a+0.5, labels = paste0(c(1:143)), las=2, cex.axis=1.7)
axis(1, at = BCADtoBP(c(-900, -700, -500, -300, -100, 100, 300, 500, 700, 900)), labels=c('900BC','700BC', '500BC', '300BC', '100BC', '100AD', '300AD', '500AD', '700AD', '900AD'), tck=-0.01, cex.axis=1.7)
axis(3, at = seq(2700, 1000, -200), labels=paste0(seq(2700, 1000, -200),'BP'), tck=-0.01, cex.axis=1.7)
axis(1, at = BCADtoBP(c(-600, -400, -200, 1, 200, 400, 600)), labels=NA, tck=-0.01) #Minor tick marks
axis(3, at = seq(2700, 1000, -50), labels=NA, tck=-0.01) #Minor tick marks
box()

post.bar(c(1500,1400,1300,1200,1100,1000,900), i=96, h=0.9, a=1050, col='lightgrey')
#arrows(x0=1000, x1=500, y0=0.3, y1=90, angle = 90, code = 3, length = 0.01)
#arrows(x0=900, x1=600, y0=0.8, y1=90, angle = 90, code = 3, length = 0.01)
text(x=1200, y=97.5, "95% HPDI", cex=1.5)
text(x=1400, y=97.5,"50% HPDI", cex=1.5)
text(x=1250, y=95, "Median Posterior", cex=1.5)
text(x=1050, y=95, "Simulated value", cex=1.5)
#lines(x=c(550,500), y=c(2.35, 2.28))
theme(legend.position = "none")
dev.off()


#------------
##FIGURE 9 -- Map of Wombling Boundaries (highlighting significant boundaries)

#With the tactical simulation data from hierarchical wombling model
post.model.tac_womble_nab  <- out.comb.tac_icar.model[,paste0('nabla[',1:383,']')]  %>% round()

#Extract differences in arrival times for tactical wombling model
med.model.tac_womble_nab  <- apply(post.model.tac_womble_nab, 2, median)

#Extract proportion of MCMC sample differences which are significant over a specified time difference
prop_model_tac_womble_nab  <- data.frame(x = 1:383,
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
         edge_info$region2_id)) #shared boundary between two subareas 

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



