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
library(maptools)
library(sf)
library(rgeos)
library(viridis)
library(latex2exp)
library(gridExtra)
library(diagram)
library(quantreg)
library(coda)

#===============================================================================
# Generate Spatial Window for Analyses: Sub-Saharan Africa ----
sf_subsah_africa <- ne_countries(continent = "Africa", returnclass = "sf") %>%
  filter_all(., any_vars(str_detect(., "Sub-Saharan"))) %>% 
  filter(name_en != "Madagascar") #We focus on mainland sub-Saharan Africa

sp_ss_africa <- sf_subsah_africa %>% as("Spatial") #convert sf to sp object

sampling_win <- as(sp_ss_africa, "SpatialPolygons") |>  unionSpatialPolygons(IDs = rep(1, nrow(sp_ss_africa)))
sampling_win <- disaggregate(sampling_win) #create new raster layer with higher resolution (smaller cells)
sampling_win  <- sampling_win[order(raster::area(sampling_win), decreasing=TRUE)]

win.sf  <- as(sampling_win,'sf')
win = sampling_win


# Fixed params
n = 300
origin.point = c(11.50, 3.82)
beta0 = 3000
beta1 = 1
sigma = 100
seed = 144231

etasq = 0.15
rho = 200


#-------------

# Data ---
dat  <- list()


gpqrSim2 <- function(win, n=600, seed=123, beta0=3000, beta1=0.7, sigma=100, etasq=0.05, rho=100, origin.point=c(11.50,3.82))
{
  require(nimbleCarbon)
  require(rcarbon)
  require(sf)
  require(sp)
  require(maptools)
  set.seed(seed)
  
  sites <- spsample(win, n = n, type = 'random') 
  dist_mat  <- spDists(sites, longlat=TRUE)
  dist_org  <-  spDistsN1(sites, origin.point, longlat=TRUE)
  
  #Covariance matrix 
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
  
  # Handle constraints on dipsersal rate
  dat$lim  <- rep(1, n)
  dat$lim2  <- rep(1, n)
  
  dispersalmodel <- nimbleCode({
    for (i in 1:N){
      # Model
      rate[i] <- -1/(s[i]-beta1)
      lim[i] ~ dconstraint(rate[i] > 0)
      
      # slope[i] <- (s[i]-beta1)
      # lim2[i] ~ dconstraint( slope[i] < (beta0/dist_org[i]) )
      
      mu[i] <- beta0 + (s[i]-beta1)*dist_org[i]
      theta[i] ~ dnorm(mean=mu[i], sd=sigma)
    }
    mu_s[1:N] <- 0
    cov_s[1:N, 1:N] <- cov_ExpQ(dist_mat[1:N, 1:N], rho, etasq, 0.000001)
    s[1:N] ~ dmnorm(mu_s[1:N], cov = cov_s[1:N, 1:N])
  })
  
  #Define Parameters
  constants  <- list()
  constants$N <- n
  constants$dist_org <- dist_org
  constants$dist_mat <- dist_mat
  constants$beta0  <- beta0
  constants$beta1  <- beta1
  constants$sigma  <- sigma
  constants$etasq  <- etasq
  constants$rho  <- rho
  
  #Simulate
  set.seed(seed)
  simModel  <- nimbleModel(code=dispersalmodel, data=dat, constants=constants)
  simModel$simulate('s')
  simModel$simulate('mu')
  simModel$simulate('theta')
  
  #Assign calibration curve
  cc <- ifelse((sites@coords[,2]>=0), 'intcal20', 'shcal20')
  
  #Combine Results
  out.df  <- data.frame(ID=1:n, theta=simModel$theta, calCurve=cc)
  #out.df$cra  <- round(uncalibrate(round(out.df$theta))$ccCRA) #Can only specify one curve currently for all dates (therefore can't specify calCurves=out.df$calCurve). Default intcal20. TODO: Check if this is problem?
  out.df$cra.error  <- 20
  # out.df$med.date  <- medCal(calibrate(out.df$cra,
  #                                      out.df$cra.error,
  #                                      calCurve=out.df$calCurve,
  #                                      verbose=F))
  out.df$s  <- simModel$s
  out.df$rate  <- -1/(out.df$s - beta1)
  out.df$mu  <- beta0 + (out.df$s - beta1)*dist_org
  out  <- as(sites, 'SpatialPointsDataFrame')
  out@data  <- out.df
  out.sf  <- as(out, 'sf')
  
  #Store Output
  return(out.sf)
}


tm3  <- gpqrSim2(win = win,
                 n = n,
                 beta0 = beta0,
                 beta1 = beta1,
                 sigma = sigma,
                 origin.point = origin.point,
                 etasq = etasq,
                 rho = rho,
                 seed = seed)