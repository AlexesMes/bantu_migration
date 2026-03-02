library(pastclim)
library(here)
library(sf)
library(dplyr)
library(stringr)
library(vegan)
library(ggplot2)
library(gridExtra)

rm(list = ls())
`%!in%` <- Negate(`%in%`)

#-------------------------------------------------------------------------------
# Load and prepare data ----
set_data_path(path_to_nc = paste0(here(), '/data/environment'))

load(here('data','trig_d38.RData'))

load(here('data','sample_window.RData'))
#download_dataset(dataset = "Krapp2021") #download pastclim data

#List of environmental variables
enviro_vars <- get_vars_for_dataset(dataset = "Krapp2021") 
enviro_vars <- enviro_vars[enviro_vars != 'biome'] #remove biome (for now) as it is a categorical variable
#get_vars_for_dataset("Krapp2021", details=TRUE)
#-------------------------------------------------------------------------------
##ENTIRE SAMPLING WINDOW ---

#SpatVector of sampling window
sampling_win_sv <- vect(sampling_win)

#Get environmental variable information over the sampling region (Eastern and Southern Africa) at 2000BP
ea_climate_2k <- region_slice(
  time_bp = -2000,
  bio_variables = enviro_vars,
  dataset = "Krapp2021",
  crop = sampling_win_sv
)

#Convert to dataframe
ea_clim_2k_df <- as.data.frame(ea_climate_2k, xy=TRUE)

#-----
#Plot

##Plot environmental variables
# pdf(file=here('output','figures_supplementary','fig_env_vars.pdf'), width=20, height=20)
# terra::plot(ea_climate_2k[[c("npp","rugosity")]],
#             main = var_labels(ea_climate_2k[[c("npp","rugosity")]], dataset = "Krapp2021", abbreviated = TRUE))
# dev.off()

##Plot correlation matrix
# vars <- c("bio01","bio04","bio05","bio06","bio07","bio08","bio09",
#           "bio10","bio11","bio12","bio13","bio14","bio15","bio16",
#           "bio17","bio18","bio19","npp","altitude","rugosity")
# 
# ea_clim_2k_df <- ea_clim_2k_df[vars] 
# ea_clim_2k_df <- na.omit(ea_clim_2k_df)
# cor_mat <- cor(ea_clim_2k_df) #calculate the correlation matrix
# 
# # 5. plot with corrplot
# corrplot(
#   cor_mat,
#   method = "color",
#   type = "upper",           # shows upper triangle
#   tl.col = col_mat,         # variable label colors
#   tl.cex = 0.8,
#   addgrid.col = "grey80",
#   col = colorRampPalette(c("blue", "white", "red"))(200),
#   diag = FALSE
# )

#-------------------------------------------------------------------------------
##HEX SAMPLING REGDIONS ----

#Apply the function to convert sfc_GEOMETRY to SpatVector to each row's geometry
hex_area_win_sv <- lapply(hex_area_win$geometry, function(sfc){
  return(vect(sfc))
  })

#Get environmental variable information in each hex region at 2000BP
hex_area_clim_2k <- lapply(hex_area_win_sv, function(raster) {
  region_slice(time_bp = -2000,
               bio_variables = enviro_vars,
               dataset = "Krapp2021",
               crop = raster)
})

#Summarize the mean values in each hex for each bioclimatic variable
mean_hex_clim_list <- lapply(hex_area_clim_2k, function(env_in_hex){
  #Convert each hex spatRaster to a dataframe
  env_hex_area_df <- as.data.frame(env_in_hex)
  
  #Summarize the mean values
  env_hex_means <- env_hex_area_df %>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))
  return(env_hex_means)
})

#Combine into a single dataframe
mean_hex_clim_df <- do.call(rbind, mean_hex_clim_list) 

# #Impute missing values (for d=2.9)
# mean_hex_clim_df[63,] <- mean_hex_clim_df[51,]

#Scale variables
scaled_hex_clim_df <- scale(mean_hex_clim_df)

#-------------------------------------------------------------------------------
## Save environmental variables on a R image file ----
save(mean_hex_clim_df, scaled_hex_clim_df, file=here('data','Krapp_enviro_variables.RData'))
