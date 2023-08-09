# Load Libraries and spatial data ----
library(here)
library(rcarbon) 
library(dplyr)
library(coda)

#-------------------------------------------------------------------------------
# Table 1 (Rhat, ESS, and Posterior Summaries of beta0, beta1, rhosq, and etasq for tau=0.9) ----
load(here('output','gpqr_tau90.RData'))

gpqr.tau90.comb  <- do.call(rbind, gpqr_tau90)
params = c('beta0', 'beta1', 'rho', 'etasq')
gpqr.tau90.comb[,'beta1'] = 1/gpqr.tau90.comb[,'beta1']
meds = apply(gpqr.tau90.comb[ ,params], 2, median) %>% round(3)
lo90 = apply(gpqr.tau90.comb[ ,params], 2, function(x){HPDinterval(as.mcmc(x), prob=0.90)[1]}) %>%  round(3)
hi90 = apply(gpqr.tau90.comb[ ,params], 2, function(x){HPDinterval(as.mcmc(x), prob=0.90)[2]}) %>%  round(3)
rhats = gelman.diag(gpqr_tau90)$psrf[params, 1]
ess = effectiveSize(gpqr_tau90)[params]
params[2]='1/beta1'
table.S1 = data.frame(params, meds, lo90, hi90, rhats, ess)
write.table(table.S1,
            file = here('output','tables','table_1.csv'),
            col.names = c('Parameter', 'Median Posterior', '90% HPDI (low)', '90% HPDI (high)', 'Rhat', 'ESS'),
            sep=',',
            row.names=FALSE)


#-------------------------------------------------------------------------------
# Table 2 (Rhat, ESS, and Posterior Summaries of beta0, beta1, rhosq, and etasq for tau=0.99) ----
load(here('results','gpqr_tau99.RData'))

gpqr.tau99.comb  <- do.call(rbind, gpqr_tau99)
params = c('beta0', 'beta1', 'rho', 'etasq')
gpqr.tau90.comb[,'beta1'] = 1/gpqr.tau90.comb[,'beta1']
meds = apply(gpqr.tau99.comb[ ,params], 2, median) %>% round(3)
lo90 = apply(gpqr.tau99.comb[ ,params], 2, function(x){HPDinterval(as.mcmc(x),prob=0.90)[1]}) %>% round(3)
hi90 = apply(gpqr.tau99.comb[ ,params], 2, function(x){HPDinterval(as.mcmc(x),prob=0.90)[2]}) %>% round(3)
rhats = gelman.diag(gpqr_tau99)$psrf[params,1]
ess = effectiveSize(gpqr_tau99)[params]
params[2]='1/beta1'
table.S2 = data.frame(params,meds,lo90,hi90,rhats,ess)
write.table(table.S2,
            file = here('output','tables','table_2.csv'),
            col.names = c('Parameter', 'Median Posterior', '90% HPDI (low)', '90% HPDI (high)', 'Rhat', 'ESS'),
            sep = ',',
            row.names = FALSE)

#-------------------------------------------------------------------------------
# Table 3 Summary Data Per Region ----
load(here('data','c14.RData'))

area_freq  <- plyr::count(siteInfo, 'area_id') 
n_dates <- as.numeric(aggregate(n_dates~area_id, data=siteInfo, sum)$n_dates)

table.S3  <- data.frame(Area=area_freq$area_id, n_dates=n_dates, n_sites=area_freq$freq)

write.table(table.S3,
            file = here('output','tables','table_3.csv'),
            col.names = c('Area', 'Number of Dates', 'Number of Sites'),
            sep=',',
            row.names=FALSE)


#-------------------------------------------------------------------------------
# Table 4 posterior estimates for nu ----
load(here("output","phase_model_a.RData"))
load(here("output","phase_model_b.RData"))

out.comb.unif.modela  <- do.call(rbind, out_unif_model_a)
out.comb.unif.modelb  <- do.call(rbind, out.unif.model_b)
post.nu.modela  <- out.comb.unif.modela[,paste0('a[',1:47,']')]  %>% round() #57 Hex areas, but don't report hex_id>=48 (assume no Bantu Expansion)
post.nu.modelb  <- out.comb.unif.modelb[,paste0('a[',1:47,']')]  %>% round()
hpdi.modela  <- apply(post.nu.modela, 2, function(x){HPDinterval(as.mcmc(x), prob = .90)}) 
hpdi.modelb  <- apply(post.nu.modelb, 2, function(x){HPDinterval(as.mcmc(x), prob = .90)}) 
med.modela  <- apply(post.nu.modela, 2, median)
med.modelb  <- apply(post.nu.modelb, 2, median)


date_strc  <- function(x)
{
  x = BPtoBCAD(x)
  ifelse(x<0, paste(abs(x),'BC'), paste(x,'AD'))
}

models  <- rep(c('Model a','Model b'), each=47) 
area  <- rep(as.character(1:47), 2)
meds  <- c(date_strc(med.modela), date_strc(med.modelb))
hi90  <- c(date_strc(hpdi.modela[1,]), date_strc(hpdi.modelb[1,]))
lo90  <- c(date_strc(hpdi.modela[2,]), date_strc(hpdi.modelb[2,]))
rhat  <- c(rhat_unif_model_a$psrf[1:47, 1], rhat.unif.model_b$psrf[1:47, 1])  %>% round(digits=3)
ess  <- c(ess_unif_model_a[1:47], ess.unif.model_b[1:47]) %>% round()
table.S4  = data.frame(models, area, meds, lo90, hi90, rhat, ess)
write.table(table.S4,
           file = here('output','tables','table_4.csv'),
           col.names = c('Model', 'Area', 'Median Posterior', '90% HPDI (low)', '90% HPDI (high)', 'Rhat', 'ESS'),
           sep = ',',
           row.names = FALSE)
