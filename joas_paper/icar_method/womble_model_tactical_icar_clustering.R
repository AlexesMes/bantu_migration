#Load Libraries ----
library(here)
library(coda)
library(nimbleCarbon)
library(rcarbon)
library(dplyr)
library(parallel)
library(spdep)
library(rlist)

rm(list = ls())
`%!in%` <- Negate(`%in%`)

set.seed(123)

##ICAR model tactical simulation with uncalibrated dates -- varying sampling intensity 

#-------------------------------------------------------------------------------
## Data and Functions Setup ----
source(here('src', 'sim_data.R'))
load(here('data','trig.RData')) #nodes and edges between hex area centroids

#-------------------------------------------------------------------------------
output_varyclusterdeg_models <- list()

#Varying sample intensity (clustering in simulated data)
for (c in seq(0,1,0.1)){
  #Generate simulated data set
  sim_dataset <- sim_data(with_calbration = FALSE, k = c, n_sites = 550, n_dates = 1200)

  sites <- sim_dataset$sites
  sites_sf <- sim_dataset$sites_sf
  siteInfo <- sim_dataset$siteInfo
  sim_df <- sim_dataset$sim_df
  constants <- sim_dataset$constants
  sampling_win_proj <- sim_dataset$sampling_win_proj
  hex_area_win_proj <-sim_dataset$hex_area_win_proj

  #Combine constants
  constants <- c(constants, constants_trig)
  #-------------------------------------------------------------------------------
  ##Simulated data summary
  #Number of dates in area
  dates_in_areas_summarise <- as.data.frame(table(sites$area_id))
  #Number of sites in area
  sites_in_areas_summarise <- sites %>% group_by(area_id) %>% summarize(n_sites =n_distinct(site_id))

  #-------------------------------------------------------------------------------
  ## Initialise Parameters ----

  #Initialise regional parameters
  buffer <- 20
  delta_init <- siteInfo$diff + buffer
  alpha_init <- siteInfo$earliest + buffer/2

  #Initialise hex areas which contain sites
  init_a  <- aggregate(earliest~area_id, FUN=max, data=siteInfo) #find earliest date in each region k
  init_b  <- aggregate(latest~area_id, FUN=min, data=siteInfo) #find latest date in each area k

  #Initialise hex areas which do not contain sites
  init_empty_area <- function(init_df) {
    for(i in 1:constants$n_area){

      area_ids <- init_df$area_id #List of hex areas ids with sites

      if (i %!in% area_ids){
        empty_hex_id <- i #Id of empty hex
        neighbour_hex <- which.min(abs(i - area_ids)) #Determine closest hex neighbor which has sites. If there are more than one neighbour hex with sites, it selects the first observation (i.e. hex with the smallest id, since ids are in ascending order).
        neighbour_hex_id <- area_ids[neighbour_hex] #Determine area id of closest hex neighbor
        neighbour_date <- init_df[neighbour_hex , 2] #Select the date associated with the neighbour hex
        init_df <- rbind(init_df, c(i, neighbour_date)) #Assign this date to the empty hex
      }
    }
    return(init_df)
  }
  init_a <- init_a %>% init_empty_area() %>%  arrange(area_id)
  init_b <- init_b %>% init_empty_area() %>%  arrange(area_id)

  #Delta parameter initialisation
  init_nabla <- 0
  for (t in 1:constants$n_trans){
    m <- constants$edge_id1[t] #transition/edge t, select first area
    n <- constants$edge_id2[t] #transition/edge t, select second area

    init_nabla[t] <- (init_a$earliest[init_a$area_id==m] - init_a$earliest[init_a$area_id==n])
  }
  init_nabla_phi <- init_nabla*0.5

  #Add buffer
  init_a  <- init_a[ ,2] + buffer
  init_b  <- init_b[ ,2] - buffer

  # Initialise spatial residues
  init_phi <- init_a #(use init_a if no intercept (beta0) is included and zero_mean=0)

  #-------------------------------------------------------------------------------
  #Spatial data ----
  library(spdep)
  nb_areas <- poly2nb(as(hex_area_win_proj, 'Spatial'), queen=FALSE, row.names = hex_area_win_proj$area_ID) #neighboring areas using sp library
  nbInfo <- nb2WB(nb_areas) #transform into iCAR inputs: adjacent matrix, weights, number of neighbors (for WinBUGS)

  #-------------------------------------------------------------------------------
  ##General Setup ----
  #Data
  dat <- list(cra=sim_df$cra,
              cra_error=sim_df$cra_error,
              constraint_uniform = rep(1, constants$n_areas),
              cra_constraint = rep(1, constants$n_dates)) # Set-up constraint for ignoring inference outside calibration range
  #x1 = hex_area_win_proj$forest_present)
  #x2 = hex_area_win_proj$water_present)

  theta_init <- medCal(calibrate(dat$cra, dat$cra_error, verbose=FALSE)) #initialize theta parameter

  #Constants
  constants$adj <- nbInfo$adj
  constants$weights <- nbInfo$weights
  constants$num <- nbInfo$num
  constants$L <- length(nbInfo$adj)
  constants <- constants[names(constants) %!in% c("dist_mat",
                                                  "dist_org",
                                                  "center_coords",
                                                  "beta0",
                                                  "beta1",
                                                  "beta2",
                                                  "x1",
                                                  "x2",
                                                  "mu",
                                                  "mu2",
                                                  "tau.err",
                                                  "tau")] #remove constants which aren't used
  #Calibration constants (with dummy extensions)
  data("intcal20")
  constants$calBP <- c(1000000, intcal20$CalBP, -1000000)
  constants$C14BP <- c(1000000, intcal20$C14Age, -1000000)
  constants$C14err <- c(1000, intcal20$C14Age.sigma, 1000)

  # #-------------------------------------------------------------------------------
  ## Model W: ICAR and wombling integrating sample interdependence, i.e. the addition of a hierarchical model ----

  modelW <- function(seed, d, theta_init, alpha_init, delta_init, init_a, init_b, init_phi, constants, nburnin, thin, niter)
  {
    #Load Library
    library(nimbleCarbon)
    #Define Core Model
    model <- nimbleCode({
      #For Each Site
      for (j in 1:n_sites)
      {
        delta[j] ~ dgamma(gamma1, (gamma1-1)/gamma2)
        alpha[j] ~ dunif(max = a[id_areas[j]], min = b[id_areas[j]]);
      }

      #For Each Date
      for (i in 1:n_dates){
        theta[i] ~ dunif(min = (alpha[id_sites[i]] - (delta[id_sites[i]]+1)), max = alpha[id_sites[i]]);
        # Calibration
        mu[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]); #C14 age on the relevant calibration curve
        cra_constraint[i] ~ dconstraint(mu[i] < 50193 & mu[i] > 95) #C14 age must be within the calibration range
        sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]); #error on the calibration curve
        error[i] <- (cra_error[i]^2 + sigmaCurve[i]^2)^(1/2); #the samples' C14 error + error on calibration curve
        cra[i] ~ dnorm(mean=mu[i], sd=error[i]); #observed radiocarbon age of the sample
      }

      #For Each Region
      for (k in 1:n_areas){
        b[k] ~ dunif(0, 3000);
        constraint_uniform[k] ~ dconstraint(b[k]<a[k]); #In each area, start date of occupation, a_k, must be greater than the end date of occupation, b_k (note: BP dates in the positive direction)

        a[k] <- phi[k] # + x1[k]*beta1[k]; #beta0[k] + phi[k] + X[k]*beta1; # + x1[k]*beta1[k]; #+ x2[k]*beta2;
        #beta1[k] ~ dnorm(0, sd=300);
        #beta0[k] ~ dunif(1000, 3000);
      }

      # ICAR Model prior to capture spatial random effects
      phi[1:n_areas] ~ dcar_normal(adj[1:L], weights[1:L], num[1:n_areas], tau1, zero_mean=0)
      #beta1[1:n_areas] ~ dcar_normal(adj[1:L], weights[1:L], num[1:n_areas], tau2, zero_mean=0)

      #For Each Boundary
      for (t in 1:n_trans){
        #nabla defines the difference in arrival time across a boundary
        nabla[t] <- abs(a[edge_id1[t]] - a[edge_id2[t]]) #edge t: select first area, m, and second area, n
        #nabla_phi defines the difference in spatial residues across a boundary
        nabla_phi[t] <- abs(phi[edge_id1[t]] - phi[edge_id2[t]])
      }

      #Priors
      #beta0 ~ dunif(1000, 3000); #dnorm(2000, sd=300); #Intercept
      #beta1 ~ dnorm(0, sd=300); #dunif(-1, 1); #determining strength of forest covariate in each area
      #beta2 ~ dnorm(0, sd=300);
      tau1 ~ dgamma(0.8, 0.1);  #weak prior for ICAR model -- spatial autocorrelation precision parameter
      #tau2 ~ dgamma(0.8, 0.1);  #weak prior for ICAR model -- spatial autocorrelation precision parameter

      #tau.err <- 1/sigma^2;
      #sigma ~ dunif(0,100);

      # Hyperprior for duration
      gamma1 ~ dunif(1,20); #Hyperprior for rate
      gamma2 ~ T(dnorm(mean=200, sd=100), 1, 500) #Hyperprior for mode

    })

    #Define initial values ----
    inits <- list(a=init_a,
                  b=init_b,
                  alpha=alpha_init,
                  delta=delta_init,
                  theta=d$cra, #theta_init
                  phi=init_phi,
                  tau1=rgamma(1, shape = 0.8, rate = 0.1),
                  #tau2=rgamma(1, shape = 0.8, rate = 0.1),
                  #beta0=runif(1:constants$n_areas, 1000, 3000), #rnorm(1, 2000, 300),
                  #beta1=rnorm(1:constants$n_areas, 0, sd=300),  #rnorm(1:constants$n_areas, 0, sd=200), #runif(1:constants$n_areas, min = -1, max = 1)
                  #beta2=rnorm(1, 0, sd=300),
                  gamma1=10,
                  gamma2=200)


    # Compile and Run model	----
    model <- nimbleModel(model, constants=constants, data=d, inits=inits)
    cModel <- compileNimble(model)
    conf <- configureMCMC(model, control=list(adaptInterval=20000, adaptFactorExponent=0.1))
    conf$addMonitors(c('a','b','nabla','nabla_phi','theta','delta','alpha')) #'beta1', 'beta2','beta3'
    MCMC <- buildMCMC(conf)
    cMCMC <- compileNimble(MCMC)
    results <- runMCMC(cMCMC, niter = niter, thin = thin, nburnin = nburnin, samplesAsCodaMCMC = T, setSeed = seed)
  }

  #-------------------------------------------------------------------------------
  # MCMC Setup
  ncores  <-  4
  cl <- makeCluster(ncores)
  seeds <- c(12, 34, 56, 78)
  niter  <- 1000000
  nburnin  <- 500000
  thin  <- 100

  #Hierarchical Womble Model
  out_womble_model <-  parLapply(cl = cl,
                                 X = seeds,
                                 fun = modelW,
                                 d = dat,
                                 constants = constants,
                                 theta_init=dat$cra, #theta_init
                                 init_a = init_a,
                                 init_b = init_b,
                                 alpha_init = alpha_init,
                                 delta_init = delta_init,
                                 init_phi = init_phi,
                                 niter = niter,
                                 nburnin = nburnin,
                                 thin = thin)
  out_womble_model <- mcmc.list(out_womble_model)

  ## Diagnostics ----
  rhat_womble_model  <- gelman.diag(out_womble_model, multivariate = FALSE)
  ess_womble_model  <- effectiveSize(out_womble_model)
  agg_womble_model <- agreementIndex(dat$cra,
                                     dat$cra_error,
                                     calCurve = "intcal20",
                                     theta = out_womble_model[[1]][ , grep("theta", colnames(out_womble_model[[1]]))],
                                     verbose = F)

  #-----------------------------------------------------------------------------
  #Output ----
  output <- list(out_womble_model,
                 rhat_womble_model,
                 ess_womble_model,
                 agg_womble_model)

  output_varyclusterdeg_models <- list.append(output_varyclusterdeg_models, output)
  
}

#-------------------------------------------------------------------------------
# Save output ----
save(output_varyclusterdeg_models, file=here('output','Womblemodel_tactical_icar_clustering.RData'))
