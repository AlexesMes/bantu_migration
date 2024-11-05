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

load(here('data','trig.RData'))

load(here('data','sample_window.RData'))
#download_dataset(dataset = "Krapp2021") #download pastclim data

#get_downloaded_datasets()
#get_vars_for_dataset(dataset = "Krapp2021", details= TRUE)
#get_time_bp_steps(dataset = "Krapp2021")

#List of environmental variables
enviro_vars <- get_vars_for_dataset(dataset = "Krapp2021") 
enviro_vars <- enviro_vars[enviro_vars != 'biome'] #remove biome (for now) as it is a categorical variable


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


#Plot environmental variables
pdf(file=here('output','figures','figure_env_vars.pdf'), width=20, height=20)
terra::plot(ea_climate_2k[[1:12]],
            main = var_labels(ea_climate_2k[[1:12]], dataset = "Krapp2021", abbreviated = TRUE))
terra::plot(ea_climate_2k[[13:22]],
            main = var_labels(ea_climate_2k[[13:22]], dataset = "Krapp2021", abbreviated = TRUE))
dev.off()

#-------------------------------------------------------------------------------
##HEX SAMPLING REGDIONS ----

#Apply the function to convert sfc_GEOMETRY to SpatVector to each row's geometry
remove_areaID <- c(1, 2, 3, 4, 6) #The Bantu hadn't settled in this area by the time the dutch arrived in the Cape (~1600AD). To back this up there are no EIA sites in these regions.
new_areaID <- setdiff(1:41, remove_areaID)
  
hex_area_win <- hex_area_win %>% 
  filter(area_ID %in% new_areaID) 

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

#-------------------------------------------------------------------------------
#PCA Analysis for environmental variables over entire sampling window ----
#for extra help see: https://www.datacamp.com/tutorial/pca-analysis-r
library('corrr')
library('ggcorrplot')
library("FactoMineR")
library('factoextra')

#Remove location coloumns
env_vars_df <- ea_clim_2k_df[, -which(names(ea_clim_2k_df) %in% c("x", "y"))]

#Checks
str(env_vars_df) #Check only numerical values
colSums(is.na(env_vars_df)) #Check for missing values (which will bias a PCA)

#Normalising the data -- ensuring that each attribute has the same level of contribution, preventing one variable from dominating others
ea_clim_2k_nom <- scale(env_vars_df)

#Correlation matrix
corr_matrix <- cor(ea_clim_2k_nom)
ggcorrplot(corr_matrix) #The higher the value, the most positively correlated the two variables are. The closer the value to -1, the most negatively correlated they are.

#PCA 
ea_clim_pca <- princomp(corr_matrix)
summary(ea_clim_pca) #Examine cumulative proportion: how much of the total variance does the first principal component explain? 
fviz_eig(ea_clim_pca, addlabels = TRUE) #Scree plot: Check number of PC to retain #We see the first 3 principal components (PCs) explain 96% of the total variance in the data

#What do the principal components mean?
#Exploring how the first three PCs relate to each column using the loadings of each PC
pc_loadingmat <-  ea_clim_pca$loadings[ , 1:3]

#Biplot
fviz_pca_var(ea_clim_pca, col.var = "black") 
#Note: 
#(i) all the variables that are grouped together are positively correlated to each other
#(ii) the higher the distance between the variable and the origin, the better represented that variable is
#(iii) variables that are negatively correlated are displayed to the opposite sides of the biplot’s origin


#Cos2 plot
#How much each variable is represented in a given component
fviz_cos2(ea_clim_pca, choice="var", axes=1:3) #change axes to 1 or 1:2
#Note:
#(i) A low value means that the variable is not perfectly represented by that component. 
#(ii) A high value, on the other hand, means a good representation of the variable on that component.
#The top variables with the highest cos2 contribute the most to PC1, PC2, and PC3


#Biplot and Cos2 plot combined
pdf(file=here('output','figures','figure_env_PCA_biplot.pdf'), width=15, height=10)
fviz_pca_var(ea_clim_pca, 
             col.var = "cos2",
             gradient.cols = c("black", "orange", "green"),
             repel = TRUE)
dev.off()

# Extract the top 3 principal components
ea_clim_pca_scores <- predict(ea_clim_pca, newdata = ea_clim_2k_df)
ea_clim_pca_3 <- ea_clim_pca_scores[, 1:3] #Top 3 PCs

# Combine the top 3 PCA scores with the original coordinates
ea_clim_top3pca <- cbind(ea_clim_2k_df[, c("x", "y")], ea_clim_pca_3)
colnames(ea_clim_top3pca) <- c("x", "y", "PC1", "PC2", "PC3")

# View the transformed dataframe
head(ea_clim_top3pca)


#----
#Transform hex areas into PC co-ordinate space
hex_clim_pca <- predict(ea_clim_pca, newdata=mean_hex_clim_df)[ ,1:3] #top 3 PCs

#Add area_ID 
hex_clim_pca <- as.data.frame(hex_clim_pca) %>% 
  mutate(area_ID = new_areaID)

#----
#Calculate the differences in PC components between the hex-areas
#PC1
hex_clim_pc1_dists <- dist(hex_clim_pca$Comp.1)
hex_clim_pc1_mat <- as.matrix(hex_clim_pc1_dists)

#PC2
hex_clim_pc2_dists <- dist(hex_clim_pca$Comp.2)
hex_clim_pc2_mat <- as.matrix(hex_clim_pc2_dists)

#PC3
hex_clim_pc3_dists <- dist(hex_clim_pca$Comp.3)
hex_clim_pc3_mat <- as.matrix(hex_clim_pc3_dists)



#-------------------------------------------------------------------------------
##Bayesian ICAR Models -- calculating differences in arrival times
#Load Data ----
load(here("output", "ICAR_model_a.RData")) #model (i) -- no sample interdependence
load(here("output","ICAR_model_b.RData")) #model (ii) -- hierarchical structure


#Extract arrival time (ak) distribution information
extract_arrival <- function(x)
{
  tmp = do.call(rbind, x)
  tmp2 = tmp[ , grep('^a\\[',colnames(tmp))]
  qta = apply(tmp2, 2, quantile, prob=c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))
  return(qta)
}

#Extract quantile information for models
tmp.i = extract_arrival(out_icar_model_a)
tmp.ii = extract_arrival(out_icar_model_b) 

mean_arrival.i <- tmp.i[4,]   #50% quantile
mean_arrival.i <- mean_arrival.i[-remove_areaID] #remove area IDs we aren't interested in
mean_arrival.i <- as.data.frame(mean_arrival.i)

mean_arrival.ii <- tmp.ii[4,]   #50% quantile
mean_arrival.ii <- mean_arrival.ii[-remove_areaID] #remove area IDs we aren't interested in
mean_arrival.ii <- as.data.frame(mean_arrival.ii)



#----
#Calculate the differences in mean arrival times between the hex-areas
#Model i
mean_arrival.i_dists <- dist(mean_arrival.i)
mean_arrival.i_mat <- as.matrix(mean_arrival.i_dists)

#Model ii
mean_arrival.ii_dists <- dist(mean_arrival.ii)
mean_arrival.ii_mat <- as.matrix(mean_arrival.ii_dists)

#-------------------------------------------------------------------------------
#Mantel Tests

#Model 1 against ecological PC1, PC2, and PC3
m1_pc1 <- mantel(hex_clim_pc1_mat, mean_arrival.i_mat) #mantel(hex_clim_pc1_mat, mean_arrival.i_mat, method = "spearman", permutations = 9999, na.rm = TRUE)
m1_pc2 <- mantel(hex_clim_pc2_mat, mean_arrival.i_mat)
m1_pc3 <- mantel(hex_clim_pc3_mat, mean_arrival.i_mat)

#Model 2 against ecological PC1, PC2, and PC3
m2_pc1 <- mantel(hex_clim_pc1_mat, mean_arrival.ii_mat)
m2_pc2 <- mantel(hex_clim_pc2_mat, mean_arrival.ii_mat)
m2_pc3 <- mantel(hex_clim_pc3_mat, mean_arrival.ii_mat)


#-------------------------------------------------------------------------------
##Plots ---

pdf(file=here('output','figures','figure_env_PCA_mantel.pdf'), width=15, height=10)
par(mfrow=c(2,3))

#Model 1
m1_pc1_plot <- t(rbind(hex_clim_pc1_dists, mean_arrival.i_dists))
plot(m1_pc1_plot[,2], m1_pc1_plot[,1],
     xlab = "Difference in arrival times",
     ylab = "Difference in PC1",
     main = paste0("Model 1: no sample interdependence; rho = ", round(m1_pc1$statistic,3)))

m1_pc2_plot <- t(rbind(hex_clim_pc2_dists, mean_arrival.i_dists))
plot(m1_pc2_plot[,2], m1_pc2_plot[,1],
     xlab = "Difference in arrival times",
     ylab = "Difference in PC2",
     main = paste0("Model 1: no sample interdependence; rho = ", round(m1_pc2$statistic,3)))

m1_pc3_plot <- t(rbind(hex_clim_pc3_dists, mean_arrival.i_dists))
plot(m1_pc3_plot[,2], m1_pc3_plot[,1],
     xlab = "Difference in arrival times",
     ylab = "Difference in PC3",
     main = paste0("Model 1: no sample interdependence; rho = ", round(m1_pc3$statistic,3)))

#----
#Model 2
m2_pc1_plot <- t(rbind(hex_clim_pc1_dists, mean_arrival.ii_dists))
plot(m2_pc1_plot[,2], m2_pc1_plot[,1],
     xlab = "Difference in arrival times",
     ylab = "Difference in PC1",
     main = paste0("Model 2: hierarchical structure; rho = ", round(m2_pc1$statistic,3)))

m2_pc2_plot <- t(rbind(hex_clim_pc2_dists, mean_arrival.ii_dists))
plot(m2_pc2_plot[,2], m2_pc2_plot[,1],
     xlab = "Difference in arrival times",
     ylab = "Difference in PC2",
     main = paste0("Model 2: hierarchical structure; rho = ", round(m2_pc2$statistic,3)))

m2_pc3_plot <- t(rbind(hex_clim_pc3_dists, mean_arrival.ii_dists))
plot(m2_pc3_plot[,2], m2_pc3_plot[,1],
     xlab = "Difference in arrival times",
     ylab = "Difference in PC3",
     main = paste0("Model 2: hierarchical structure; rho = ", round(m2_pc3$statistic,3)))

dev.off()

#-------------------------------------------------------------------------------
#Repeat Mantel test, but only allow delaunay transitions 

#Filter edge info to remove hexs where Bantu Expansion didn't reach
f_transitions <- edge_info %>% 
  filter(region1_id %in% new_areaID, region2_id %in% new_areaID) %>% 
  select(region1_id, region2_id)

#--------------------------
##Arrival times ----
#MODEL 1 
#Add area ID
mean_arrival.i <- mean_arrival.i %>% 
  mutate(area_ID = new_areaID) %>% 
  rename(ak = mean_arrival.i)

# Calculate distances for the specified pairs
f_mean_arrival.i_dists <- f_transitions %>%
  rowwise() %>%
  mutate(
    ak_dist = abs(mean_arrival.i$ak[mean_arrival.i$area_ID ==region1_id] - mean_arrival.i$ak[mean_arrival.i$area_ID == region2_id])) %>% 
  select(ak_dist)

f_mean_arrival.i_vec <- as.vector(as.matrix(f_mean_arrival.i_dists))

#MODEL 2 
#Add area ID
mean_arrival.ii <- mean_arrival.ii %>% 
  mutate(area_ID = new_areaID) %>% 
  rename(ak = mean_arrival.ii)

# Calculate distances for the specified pairs
f_mean_arrival.ii_dists <- f_transitions %>%
  rowwise() %>%
  mutate(
    ak_dist = abs(mean_arrival.ii$ak[mean_arrival.ii$area_ID ==region1_id] - mean_arrival.ii$ak[mean_arrival.ii$area_ID == region2_id])) %>% 
  select(ak_dist)

f_mean_arrival.ii_vec <- as.vector(as.matrix(f_mean_arrival.ii_dists))

#--------------------------
##Principal Components
#PC1 ----
hex_clim_pc1 <- hex_clim_pca %>% 
  select(area_ID, Comp.1)

# Calculate distances for the specified pairs
f_hex_clim_pc1_dists <- f_transitions %>%
  rowwise() %>%
  mutate(
    pc1_dist = abs(hex_clim_pc1$Comp.1[hex_clim_pc1$area_ID ==region1_id] - hex_clim_pc1$Comp.1[hex_clim_pc1$area_ID == region2_id]))

f_hex_clim_pc1_vec <- as.vector(as.matrix(f_hex_clim_pc1_dists$pc1_dist))

#PC2 ----
hex_clim_pc2 <- hex_clim_pca %>% 
  select(area_ID, Comp.2)

# Calculate distances for the specified pairs
f_hex_clim_pc2_dists <- f_transitions %>%
  rowwise() %>%
  mutate(
    pc2_dist = abs(hex_clim_pc2$Comp.2[hex_clim_pc2$area_ID ==region1_id] - hex_clim_pc2$Comp.2[hex_clim_pc2$area_ID == region2_id]))

f_hex_clim_pc2_vec <- as.vector(as.matrix(f_hex_clim_pc2_dists$pc2_dist))

#PC3 ----
hex_clim_pc3 <- hex_clim_pca %>% 
  select(area_ID, Comp.3)

# Calculate distances for the specified pairs
f_hex_clim_pc3_dists <- f_transitions %>%
  rowwise() %>%
  mutate(
    pc3_dist = abs(hex_clim_pc3$Comp.3[hex_clim_pc3$area_ID ==region1_id] - hex_clim_pc3$Comp.3[hex_clim_pc3$area_ID == region2_id]))

f_hex_clim_pc3_vec <- as.vector(as.matrix(f_hex_clim_pc3_dists$pc3_dist))


#--------------------------
#Calculate correlations (e.g. spearman's rank correlation coefficient) ---
#Model 1 against ecological PC1
cor_m1_pc1 <- cor(f_hex_clim_pc1_vec, f_mean_arrival.i_vec, method="spearman")
cor_m1_pc2 <- cor(f_hex_clim_pc2_vec, f_mean_arrival.i_vec, method="spearman")
cor_m1_pc3 <- cor(f_hex_clim_pc3_vec, f_mean_arrival.i_vec, method="spearman")

cor_m2_pc1 <- cor(f_hex_clim_pc1_vec, f_mean_arrival.ii_vec, method="spearman")
cor_m2_pc2 <- cor(f_hex_clim_pc2_vec, f_mean_arrival.ii_vec, method="spearman")
cor_m2_pc3 <- cor(f_hex_clim_pc3_vec, f_mean_arrival.ii_vec, method="spearman")


#Plot correlations ---
#Model 1
f_m1_pc1 <- as.data.frame(t(rbind(f_hex_clim_pc1_vec, f_mean_arrival.i_vec)))
p1 <- ggplot(f_m1_pc1, aes(y = f_hex_clim_pc1_vec, 
                     x = f_mean_arrival.i_vec)) +
        geom_point() +
        geom_smooth(method = "lm", col = "blue") +
        ggtitle(paste0("Model 1: no sample interdependence; rho = ", round(cor_m1_pc1,3))) +
        xlab("Difference in arrival times") +
        ylab("Difference in PC1")

f_m1_pc2 <- as.data.frame(t(rbind(f_hex_clim_pc2_vec, f_mean_arrival.i_vec)))
p2 <- ggplot(f_m1_pc2, aes(y = f_hex_clim_pc2_vec, 
                     x = f_mean_arrival.i_vec)) +
        geom_point() +
        geom_smooth(method = "lm", col = "blue") +
        ggtitle(paste0("Model 1: no sample interdependence; rho = ", round(cor_m1_pc2,3))) +
        xlab("Difference in arrival times") +
        ylab("Difference in PC2")


f_m1_pc3 <- as.data.frame(t(rbind(f_hex_clim_pc3_vec, f_mean_arrival.i_vec)))
p3 <- ggplot(f_m1_pc3, aes(y = f_hex_clim_pc3_vec, 
                     x = f_mean_arrival.i_vec)) +
        geom_point() +
        geom_smooth(method = "lm", col = "blue") +
        ggtitle(paste0("Model 1: no sample interdependence; rho = ", round(cor_m1_pc3,3))) +
        xlab("Difference in arrival times") +
        ylab("Difference in PC3")

#Model 2
f_m2_pc1 <- as.data.frame(t(rbind(f_hex_clim_pc1_vec, f_mean_arrival.ii_vec)))
p4 <- ggplot(f_m2_pc1, aes(y = f_hex_clim_pc1_vec, 
                     x = f_mean_arrival.ii_vec)) +
        geom_point() +
        geom_smooth(method = "lm", col = "blue") +
        ggtitle(paste0("Model 2: hierarchical structure; rho = ", round(cor_m2_pc1,3))) +
        xlab("Difference in arrival times") +
        ylab("Difference in PC1")

f_m2_pc2 <- as.data.frame(t(rbind(f_hex_clim_pc2_vec, f_mean_arrival.ii_vec)))
p5 <- ggplot(f_m2_pc2, aes(y = f_hex_clim_pc2_vec, 
                     x = f_mean_arrival.ii_vec)) +
        geom_point() +
        geom_smooth(method = "lm", col = "blue") +
        ggtitle(paste0("Model 2: hierarchical structure; rho = ", round(cor_m2_pc2,3))) +
        xlab("Difference in arrival times") +
        ylab("Difference in PC2")

f_m2_pc3 <- as.data.frame(t(rbind(f_hex_clim_pc3_vec, f_mean_arrival.ii_vec)))
p6 <- ggplot(f_m2_pc3, aes(y = f_hex_clim_pc3_vec, 
                     x = f_mean_arrival.ii_vec)) +
        geom_point() +
        geom_smooth(method = "lm", col = "blue") +
        ggtitle(paste0("Model 2: hierarchical structure; rho = ", round(cor_m2_pc3,3))) +
        xlab("Difference in arrival times") +
        ylab("Difference in PC3")

pdf(file=here('output','figures','figure_env_PCA_corr.pdf'), width=30, height=20)
grid.arrange(p1, p2, p3, p4, p5, p6, ncol=3, nrow=2)
dev.off()

#-------------------------------------------------------------------------------
##Permutation Tests for significance

#Difference in delaunay transitions PC means
org_trans_mean <- mean(f_hex_clim_pc1_vec)
org_trans_mean2 <- mean(f_hex_clim_pc2_vec)

#Permutation test
permutation.test <- function(outcome, group, transitions, n){
  distribution=c()
  result=0
  for(i in 1:n){
    group = sample(group, length(group), FALSE)
    
    pair_dists <- transitions %>%
      rowwise() %>%
      mutate(dist = abs(outcome[group==region1_id]-outcome[group==region2_id]))
    
    distribution[i] = mean(pair_dists$dist)
  }
  result=sum(abs(distribution) >= abs(org_trans_mean))/(n)
  return(list(result, distribution))
}


test1 <- permutation.test(hex_clim_pc1$Comp.1, 
                          hex_clim_pc1$area_ID, 
                          f_transitions, 
                          10000)

test2 <- permutation.test(hex_clim_pc2$Comp.2, 
                          hex_clim_pc2$area_ID, 
                          f_transitions, 
                          10000)

#plot
hist(test1[[2]], 
     breaks=50, 
     col='grey', 
     main="Permutation Distribution of PC1", 
     las=1, 
     xlab='Difference in PC1 values',
     xlim=c(150, 350))
abline(v=org_trans_mean, lwd=3, col="red")

hist(test2[[2]], 
     breaks=50, 
     col='grey', 
     main="Permutation Distribution of PC2", 
     las=1, 
     xlab='Difference in PC2 values',
     xlim=c(150, 350))
abline(v=org_trans_mean2, lwd=3, col="red")




#-----
#Select biovaribales and add area_ID 
hex_annualT <- mean_hex_clim_df %>% 
  select(bio07) %>% 
  mutate(area_ID = new_areaID)

# Calculate distances for the specified pairs
hex_annualT_pair_dist <- f_transitions %>%
  rowwise() %>%
  mutate(
    annT_dist = abs(hex_annualT$bio07[hex_annualT$area_ID ==region1_id] 
                    - hex_annualT$bio07[hex_annualT$area_ID == region2_id]))
#Mean difference in pair dists
org_trans_bio7mean <- mean(hex_annualT_pair_dist$annT_dist)


test3 <- permutation.test(hex_annualT$bio07, 
                          hex_annualT$area_ID, 
                          f_transitions, 
                          10000)

#plot
hist(test3[[2]], 
     breaks=50, 
     col='grey', 
     main="Permutation Distribution of Annual Temperature range", 
     las=1, 
     xlab='Difference in bio07 values',
     xlim=c(1.5, 6))
abline(v=org_trans_bio7mean, lwd=3, col="red")

  

#---
#All pairs
hex_annualT_dists <- dist(hex_annualT)
mean(hex_annualT_dists)


#Permutation test for all pairs
permutation.test2 <- function(df_vars, n){
  distribution=c()
  result=0
  group = df_vars$area_ID
  
  for(i in 1:n){
    df_vars$area_ID <- sample(group, length(group), FALSE)
    pair_dists <- dist(df_vars)
    distribution[i] = mean(pair_dists)
  }
  result=sum(abs(distribution) >= abs(org_trans_mean))/(n)
  return(list(result, distribution))
}


test_all1 <- permutation.test2(hex_annualT, 10000)
