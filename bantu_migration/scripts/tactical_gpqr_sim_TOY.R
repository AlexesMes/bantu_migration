# Load Libraries ----
library(here)
library(dplyr)
library(stringr)
library(nimbleCarbon)
library(rnaturalearth)
library(sp)
library(maptools)
library(sf)
library(ggplot2)
library(viridis)
library(rgeos)
library(rcarbon)

##SCRIPT TO SIMULATE DATA WITHOUT ERRORS 

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

sampling_win <- as(sp_ss_africa, "SpatialPolygons") |>  unionSpatialPolygons(IDs = rep(1, nrow(sp_ss_africa)))
sampling_win <- disaggregate(sampling_win) #create new raster layer with higher resolution (smaller cells)
sampling_win  <- sampling_win[order(raster::area(sampling_win), decreasing=TRUE)]

#-------------------------------------------------------------------------------
# Target Parameters ----
true_param  <- list()
true_param$n  <- 600 #number of sites & dates
true_param$origin_point <- c(11.40, 5.48) #dispersal origin point -- approximately at Ngoume
true_param$beta0 <- 3300 #approximate mean date at origin point
true_param$beta1 <- 0.4 #reciprocal of dispersal rate 
true_param$sigma <- 100 
true_param$etasq <- 0.02 #variability of dispersal rate
true_param$rho <- 350 #range of spatial autocorrelation
true_param$seed <- 1233 #random seed

# Simulate Data ----

#Gaussian Process Quantile Regression Function
gpqrSim_TOY  <- function(win, n=600, seed=123, beta0=3000, beta1=0.7, sigma=100, etasq=0.05, rho=100, origin.point=c(11.40,5.48))
{
  require(nimbleCarbon)
  require(rcarbon)
  require(sf)
  require(sp)
  require(maptools)
  set.seed(seed)
  
  #Initialise objects for storing output and site information
  out_fin_df <- data.frame()
  tot_sites <- spsample(win, n = 1, type = 'random') #dummy site
  
  while (nrow(out_fin_df) < n) { #where n is the number of sites
    
    #Generate sites and calculate distances from origin ---
    sites <- spsample(win, n = n, type = 'random') 
    dist_mat  <- spDists(sites, longlat=TRUE)
    dist_org  <-  spDistsN1(sites, origin.point, longlat=TRUE)
    
    #Covariance matrix ---
    cov_ExpQ <- nimbleFunction(run = function(dists = double(2), rho = double(0), etasq = double(0), sigmasq = double(0)) 
    {
      returnType(double(2))
      n <- dim(dists)[1]
      result <- matrix(nrow = n, ncol = n, init = FALSE)
      deltaij <- matrix(nrow = n, ncol = n, init = TRUE)
      diag(deltaij) <- 1
      for(i in 1:n)
        for(j in 1:n)
          result[i, j] <- etasq*exp(-0.5*(dists[i,j]/rho)^2) + sigmasq*deltaij[i,j]
      return(result)
    })
    Ccov_ExpQ <- compileNimble(cov_ExpQ)
    assign('cov_ExpQ', cov_ExpQ, envir=.GlobalEnv)
    
    #Model ---
    dispersalmodel <- nimbleCode({
      for (i in 1:N){
        mu[i] <- beta0 + (s[i]-beta1)*dist_org[i]
        theta[i] ~ T(dnorm(mean=mu[i], sd=sigma), 0, 55000) #Truncates to ensure theta > 0. #Yields infinity on some of the nodes -- maybe because truncation so far away from the bulk of the probability mass that the code misbehaves. Later we filter again for (theta != Inf) to get rid of these instances
      }
      mu_s[1:N] <- 0
      cov_s[1:N, 1:N] <- cov_ExpQ(dist_mat[1:N, 1:N], rho, etasq, 0.000001)
      s[1:N] ~ dmnorm(mu_s[1:N], cov = cov_s[1:N, 1:N])
    })
    
    #Define Parameters ---
    constants  <- list()
    constants$N <- n
    constants$dist_org <- dist_org
    constants$dist_mat <- dist_mat
    constants$beta0  <- beta0
    constants$beta1  <- beta1
    constants$sigma  <- sigma
    constants$etasq  <- etasq
    constants$rho  <- rho
    
    #Simulate ---
    set.seed(seed)
    simModel  <- nimbleModel(code=dispersalmodel, constants=constants)
    simModel$simulate('s')
    simModel$simulate('mu')
    simModel$simulate('theta')
    
    #Results ---
    out_df <- data.frame(ID=1:n, theta=simModel$theta)
    
    out_df$s  <- simModel$s
    out_df$rate  <- -1/(out_df$s - beta1)
    out_df$mu  <- beta0 + (out_df$s - beta1)*dist_org
    
    #Constrain 
    out_df <- out_df %>% filter((out_df$rate > 0) & (out_df$theta != Inf) & (out_df$theta >= 0))  #Constrain the simulated data -- don't allow for negative velocity (rate) or for dates (theta) to be in the future
    
    #Store
    out_fin_df  <- rbind(out_fin_df, out_df)
    tot_sites <- rbind(tot_sites, sites)
    
  }
  
  #Keep only required number of sites
  out.df <- head(out_fin_df, n=n)
  tot_sites <- tot_sites[2:(n+1)] #indices ignore dummy sites which was the first site in this df
  
  
  #Calibrate dates ---
  out.df$cra  <- round(out.df$theta)
  #Save as sf object ---
  out  <- as(sites, 'SpatialPointsDataFrame')
  out@data  <- out.df
  out.sf  <- as(out, 'sf')
  
  #Store Output ---
  return(out.sf)
}
sim_sites  <- gpqrSim_TOY(win = sampling_win,
                          n = true_param$n,
                          beta0 = true_param$beta0,
                          beta1 = true_param$beta1,
                          sigma = true_param$sigma,
                          origin.point = true_param$origin_point,
                          etasq = true_param$etasq,
                          rho = true_param$rho,
                          seed = true_param$seed)

# Save simulation output ----
save(true_param, sim_sites, file=here('data', 'tactical_sim_gpqr_TOY.RData'))