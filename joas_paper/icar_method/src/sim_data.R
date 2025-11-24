##SCRIPT TO SIMULATE DATA (with underlying ICAR structure)
sim_data <- function(with_calibration = FALSE, #use calibrated or uncalibrated radiocarbon dates
                     #structure = 'ICAR', 
                     seed=123,
                     k = 0.3, #clustering strength: 0 = random, 1 = very clustered
                     n_sites = 500, 
                     n_dates = 1200){
  
  # Load Libraries ----
  library(here)
  library(dplyr)
  library(stringr)
  library(nimbleCarbon)
  library(rnaturalearth)
  library(sf)
  library(spatstat.random)
  library(rcarbon)
  library(units)
  library(patchwork)
  
  set.seed(seed)
  
  # Load sample window data ----
  load(here('data','sample_window.RData'))
  load(here('data','trig.RData')) #nodes and edges between hex area centroids
  #-----------------------------------------------------------------------------
  # Target Parameters ----
  origin_point_proj <- st_transform(st_sfc(st_point(c(-5, 63)),crs=4326), crs=3035) #choose dispersal origin point -- approximately in hex 29
  #-----------------------------------------------------------------------------
  #Spatial data ----
  library(spdep)
  nb_areas <- poly2nb(as(hex_area_win_proj, 'Spatial'), queen=FALSE, row.names = hex_area_win_proj$area_ID) #neighboring areas using sp library 
  nbInfo <- nb2WB(nb_areas) #transform into iCAR inputs: adjacent matrix, weights, number of neighbors (for WinBUGS)
  
  constants <- constants_sw
  constants$adj <- nbInfo$adj
  constants$weights <- nbInfo$weights
  constants$num <- nbInfo$num
  constants$L <-length(nbInfo$adj)
  
  #-----------------------------------------------------------------------------
  #Generate sites and calculate distances from origin ---
  sites <- sample_clustered_points(sampling_win_proj, n_sites, k)
  
  dist_mat  <- set_units(st_distance(sites), 'km')
  dist_org  <- as.vector(set_units(st_distance(x=sites, y=origin_point_proj), 'km'))
  
  #Assign hex area id to each site ----
  sites$area_id <- as.integer(st_within(sites$geometry, hex_area_win_proj$geometry))
  
  #Simulate multiple observations at each site ----
  id_sites <- c(1:n_sites, 
                sample(1:n_sites,
                       size = n_dates - n_sites,
                       replace = TRUE,
                       prob = dnorm(1:n_sites, mean = mean(1:n_sites), sd = sd(1:n_sites)) / sum(dnorm(1:n_sites, mean = mean(1:n_sites), sd = sd(1:n_sites))))) #generate a site id for each date
  
  dates <- sites[1, ] #initialise 'dates' as first site 
  
  for(i in 2:n_dates){ #populate 'dates' with repeated sites 
    current_site_id <- id_sites[i]
    dates[i, ] <- sites[current_site_id, ]
  }
  
  #-----------------------------------------------------------------------------
  #Create environmental variables ----
  hex_area_win_proj <- hex_area_win_proj %>%
    mutate(
      forest_present = if_else(area_ID %in% c(7, 15, 23, 32, 39, 44, 6, 10, 18, 27, 35, 14, 31), 1, 0),
      water_present  = if_else(area_ID %in% c(25, 39, 51, 44, 63, 57, 56, 52), 1, 0), 
      env_type = case_when(
        forest_present == 1 & water_present == 1 ~ "Forest & Water",
        forest_present == 1 & water_present == 0 ~ "Forest Only",
        forest_present == 0 & water_present == 1 ~ "Water Only",
        TRUE ~ "Neither"))
  
  #-----------------------------------------------------------------------------
  #MODEL: ICAR Simulate ----
  sim_model <- nimbleCode({
    # Simulate spatially correlated data for all k in 1:n_areas
    for (k in 1:n_areas){
      a[k] ~ dnorm(nu[k], tau.err);
      nu[k] <- phi[k] + x1[k]*beta1 + x2[k]*beta2;
    }
    phi[1:n_areas] ~ dcar_proper(mu = mu[1:n_areas], adj=adj[1:L], num=num[1:n_areas], tau=tau, gamma=gamma) # ICAR prior to capture spatial random effects
    d[1:n_areas] ~ dcar_proper(mu = mu2[1:n_areas], adj=adj[1:L], num=num[1:n_areas], tau=tau, gamma=gamma)
    b[1:n_areas] <- a[1:n_areas] - abs(d[1:n_areas]) #duration must be positive
    #tau ~ dgamma(2, 0.5)

    for (j in 1:n_sites)
    {
      delta[j] ~ dgamma(5,(5-1)/100); #Site duration parameter.
      alpha[j] ~ dunif(max=a[id_areas[j]], min=b[id_areas[j]]);
      beta[j] <- alpha[j] - (delta[j] + 1); #The +1 ensures at a minimum where there are two dates at a site there will be 1 year between them.
      constraint_duration[j] ~ dconstraint(alpha[j]>(delta[j]+1)); #Site can't have have a duration longer than its time of first arrival
    }

    for (i in 1:n_dates){
      theta[i] ~ dunif(min=beta[id_sites[i]], max=alpha[id_sites[i]]);
      cra_constraint[i] ~ dconstraint(theta[i] > 0);
    }
  })
  
  #-----------------------------------------------------------------------------
  #Define constants ----
  sim_constants <- constants
  sim_constants$n_sites <- n_sites
  sim_constants$n_dates  <- n_dates
  sim_constants$n_areas  <- constants$n_areas
  sim_constants$id_sites  <- dates$site_id
  sim_constants$id_areas <- sites$area_id
  sim_constants$x1 <- hex_area_win_proj$forest_present
  sim_constants$beta1 <- -400 #magnitude of the effect of the forest covariate
  sim_constants$x2 <- hex_area_win_proj$water_present
  sim_constants$beta2 <- +300 #magnitude of the effect of the water covariate
  sim_constants$mu <- rep(5000, constants$n_areas) 
  sim_constants$mu2 <- rep(500, constants$n_areas)
  sim_constants$tau <- 1/(500^2)
  sim_constants$tau.err <- 0.5
  sim_constants$gamma <- 0.95
  
  #Define constraints, data, and initial values ----
  dat <- list(constraint_uniform = rep(1, sim_constants$n_areas),
              constraint_duration = rep(1, n_sites),
              cra_constraint = rep(1, n_dates))
  
  init_a <- runif(1:sim_constants$n_areas, min = 3500, max = 6500)
  init_b <- init_a - runif(1:sim_constants$n_areas, min = 50, max = 600)
  inits <- list(a = init_a,
                b = init_b)
  
  #-----------------------------------------------------------------------------
  #Simulate ----
  simModel <- nimbleModel(code = sim_model, constants = sim_constants, data = dat, inits = inits)
  nodesToSim <- simModel$getDependencies(c("a", "phi", "d", "b", "delta", "alpha", "beta", "theta"), self = T, downstream = T)

  simModel$simulate(nodesToSim)

  #-----------------------------------------------------------------------------
  ##Combine data ---- 
  cra = round(simModel$theta) #uncalibrate(round(simModel$theta))$rCRA  #uncomment if simulating calibrated dates
  #cra_error = rep(20,length(cra)) #uncomment if simulating calibrated dates
  
  sim_df <- list(cra = cra,
                 #cra_error = cra_error,
                 site_id = id_sites)
  
  ##Collect site level information ----
  sites_sf <- sites #save a copy
  sites <- sites %>% 
    st_drop_geometry() %>% 
    left_join(as.data.frame(sim_df), by='site_id')
  
  earliest_dates <- aggregate(cra ~ site_id, data = sites, FUN = max) #Earliest Date for Each Site 
  latest_dates <- aggregate(cra ~ site_id, data=sites, FUN=min) #Latest Date for Each Site
  n_dates_persite <- aggregate(cra ~ site_id, data=sites, FUN=length) #Number of Dates for Each Site
  
  siteInfo <- data.frame(site_id = earliest_dates$site_id,
                         earliest = earliest_dates$cra,
                         latest = latest_dates$cra,
                         diff = earliest_dates$cra - latest_dates$cra,
                         n_dates = n_dates_persite$cra) %>% unique()
  
  #Assign hex area id to each site ----
  siteInfo$area_id <- as.integer(st_within(sites_sf$geometry, hex_area_win_proj$geometry))
  
  #-------------------------------------------------------------------------------
  ##Save constants ----
  constants <- sim_constants[names(sim_constants) %!in% c("a", "b")]
  constants$dist_mat  <- dist_mat
  constants$dist_org  <- dist_org
  constants$n_areas  <- nrow(hex_area_win_proj)
  constants$true_a <- simModel$a
  constants$true_b <- simModel$b
  constants$true_alpha <- simModel$alpha
  constants$true_beta <- simModel$beta
  
  ##Output ----
  tactical_sim_data <- list(sites=sites, 
                            sites_sf=sites_sf, 
                            siteInfo=siteInfo, 
                            sim_df=sim_df, 
                            constants=constants, 
                            sampling_win_proj=sampling_win_proj, 
                            hex_area_win_proj=hex_area_win_proj)
  return(tactical_sim_data)
}

#-------------------------------------------------------------------------------
#Controlling level of clustering in generated sites
sample_clustered_points <- function(sampling_win_proj, n_sites, k,
                                    min_clusters = 2, max_clusters = 20,
                                    sigma_loose = 1700e3, #Loose/almost random: σ ≈ 1700 km -- starts to look like a random spread at the continental Europe scale
                                    sigma_tight = 20e3) { #Very tight clusters: σ ≈ 20 km -- clear, compact clusters, metropolitan-scale
  
  #Convert sampling window to spatstat format
  win_spat <- as.owin(sampling_win_proj)
  
  if (k <= 0) {
    pts <- runifpoint(n_sites, win = win_spat) #random sampling
  } else {
    #Map k to parameters -- k controls number of clusters and the tightness of those clusters
    n_clusters <- round(min_clusters + (max_clusters - min_clusters) * k) #interpolate number of clusters (1 = random, many = clustered)
    sigma <- sigma_loose - (sigma_loose - sigma_tight) * k #interpolate cluster spread (sigma_loose and sigma_tight are in meters)
    
    #Translate to Thomas process parameters
    kappa <- n_clusters / area.owin(win_spat) #cluster center density
    mu <- n_sites / n_clusters #mean points per cluster
    
    #Generate clustered points
    pp <- rThomas(kappa = kappa, sigma = sigma, mu = mu, win = win_spat)
    pts <- as.data.frame(pp)[, c("x", "y")]
    
    #Enforce exactly n_sites
    if (nrow(pts) > n_sites) {
      pts <- pts[sample(nrow(pts), n_sites), ]
    } else if (nrow(pts) < n_sites) {
      extra <- as.data.frame(runifpoint(n_sites - nrow(pts), win = win_spat))[, c("x", "y")] #Top up with uniform random points
      pts <- rbind(pts, extra)
    }
  }
  
  #Convert to sf
  st_as_sf(as.data.frame(pts), coords = c("x", "y"), crs = 3035) %>%
    mutate(site_id = row_number())
}

