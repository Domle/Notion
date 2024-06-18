### In this script we compute a final raster stack object, that contains most of the calculated metrices and statistics afterwards used for plotting.

### NOTE: 	When running the code, please always double check input and output directories, to ensure you know where exactly files get saved so they
###			    can be referred as the correct input files within each consecutive script.

# Author: 	
#     Dominic Eriksson
#			Environmental Physics Group, UP
# 		ETH Zurich
#			Switzerland


# Input files:
#       1.  Several calculated statistics:
#               1. PredictorAnnualMean.grd --> annual averages of environmental predictors

# Output files:
#       1. Raster stack as .grd file, later used for final plotting

# derikssong@ethz.ch; 18th of June 2024


# Liebraries
library(raster)

# Directories
wd_out <- "/net/kryo/work/deriksson/Projects/Notion_DE/Code/Tables/1_Output/"

## Get lon and lat and environmental annual averages
# Load environmental data
coords <- stack("/home/deriksson/Projects/Notion_DE/Data/Predictor_variables/2_Output/PredictorAnnualMean.grd")
names(coords) <- paste0( "EnvPar_", names(coords) )

# Add individual ensemble biogeographies majority vote based
filenames <- list.files("/home/deriksson/Projects/Notion_DE/Code/Ensembles/Annual_individual_ensemble_biogeographies/hsi/")
filenames <- grep( "grd", filenames, value = TRUE )
for(i in seq_along(filenames)){

    # Extract name of taxa
    taxon <- gsub( "_annual_ensemble_mean_sd_individualEnsembleMembers.grd", "", filenames[i] )

    # Open file
    r <- brick( paste0("/home/deriksson/Projects/Notion_DE/Code/Ensembles/Annual_individual_ensemble_biogeographies/hsi/", filenames[i]) )
    r <- r[[1:2]]
    names(r) <- c(  paste0("Species_", taxon, "_AnnualEnsembleMean_hsi"), paste0("Species_", taxon, "_AnnualEnsembleSD_hsi"))
    coords <- stack(coords, r)
}

# Add individual annual ensemble biogeographies PA
filenames <- list.files("/home/deriksson/Projects/Notion_DE/Code/Ensembles/Annual_individual_ensemble_biogeographies/pa/")
filenames <- grep( "grd", filenames, value = TRUE )
for(i in seq_along(filenames)){

    # Extract name of taxa
    taxon <- gsub( "annual_ensemble_mean_sd_individualEnsembleMembers.grd", "", filenames[i] )

    # Open file
    r <- brick( paste0("/home/deriksson/Projects/Notion_DE/Code/Ensembles/Annual_individual_ensemble_biogeographies/hsi/", filenames[i]) )
    r <- r[[1:2]]
    names(r) <- c(  paste0("Species_", taxon, "_AnnualEnsembleMean_pa"), paste0("Species_", taxon, "_AnnualEnsembleSD_pa"))
    coords <- stack(coords, r)
}

## Add annual ensemble richness for total, cynaos and NCDs - PA based normalized richness
list.files("/home/deriksson/Projects/Notion_DE/Code/Diversity/Diversities/pa/")

# Total
r <- brick( "/home/deriksson/Projects/Notion_DE/Code/Diversity/Diversities/pa/RichnessNormalized_annual_totalDiaz_Ensemble_and_EnsembleMembers.grd" )
r <- stack(r$mean, r$sd)
names(r) <- c("Diversity_Richness_ensembleMean_paBased", "Diversity_Richness_ensembleSD_paBased")
coords <- stack(coords, r)
# Cyanos
r <- brick( "/home/deriksson/Projects/Notion_DE/Code/Diversity/Diversities/pa/RichnessNormalized_annual_cyanoDiaz_Ensemble_and_EnsembleMembers.grd" )
r <- stack(r$mean, r$sd)
names(r) <- c("Diversity_Richness_ensembleMean_paBased_cyanos", "Diversity_Richness_ensembleSD_paBased_cyanos")
coords <- stack(coords, r)
# NCDs
r <- brick( "/home/deriksson/Projects/Notion_DE/Code/Diversity/Diversities/pa/RichnessNormalized_annual_nonCyanoDiaz_Ensemble_and_EnsembleMembers.grd" )
r <- stack(r$mean, r$sd)
names(r) <- c("Diversity_Richness_ensembleMean_paBased_ncd", "Diversity_Richness_ensembleSD_paBased_ncd")
coords <- stack(coords, r)


## Add annual ensemble richness for total, cynaos and NCDs - HSI based normalized richness
list.files("/home/deriksson/Projects/Notion_DE/Code/Diversity/Diversities/hsi/")
# Total
r <- brick( "/home/deriksson/Projects/Notion_DE/Code/Diversity/Diversities/hsi/RichnessNormalized_annual_totalDiaz_Ensemble_and_EnsembleMembers.grd" )
r <- stack(r$mean, r$sd)
names(r) <- c("Diversity_Richness_ensembleMean_hsiBased", "Diversity_Richness_ensembleSD_hsiBased")
coords <- stack(coords, r)
# Cyanos
r <- brick( "/home/deriksson/Projects/Notion_DE/Code/Diversity/Diversities/hsi/RichnessNormalized_annual_cynaoDiaz_Ensemble_and_EnsembleMembers.grd" )
r <- stack(r$mean, r$sd)
names(r) <- c("Diversity_Richness_ensembleMean_hsiBased_cyanos", "Diversity_Richness_ensembleSD_hsiBased_cyanos")
coords <- stack(coords, r)
# NCDs
r <- brick( "/home/deriksson/Projects/Notion_DE/Code/Diversity/Diversities/hsi/RichnessNormalized_annual_nonCyanoDiaz_Ensemble_and_EnsembleMembers.grd" )
r <- stack(r$mean, r$sd)
names(r) <- c("Diversity_Richness_ensembleMean_hsiBased_ncd", "Diversity_Richness_ensembleSD_hsiBased_ncd")
coords <- stack(coords, r)


## Add Betadiversity components
# Total betadiversity components
filenames <- list.files("/home/deriksson/Projects/Notion_DE/Code/Ensembles/Diversity/")
filenames <- grep( ".grd", filenames, value = TRUE )

# Stack via loop
for( i in seq_along(filenames) ){

    # Diversity index name
    diversity_index <- gsub("Annual_", "", filenames[i])
    diversity_index <- gsub("_Ensemble_and_ensembleMembers.grd", "", diversity_index)

    # Open file
    r <- brick( paste0( "/home/deriksson/Projects/Notion_DE/Code/Ensembles/Diversity/", filenames[i]) )
    r <- r[[1:2]]
    names(r) <- c( paste0("Diversity_", diversity_index, "_ensembleMean_paBased"), paste0("Diversity_", diversity_index, "_ensembleSD_paBased") )

    # We miss some cells, so we will resample our values into a new created complete (180, 360, 64800) raster object
    new_raster <- raster(nrow = 180, ncol = 360, xmn = -180, xmx = 180, ymn = -90, ymx = 90)
    res(new_raster) <- c(1, 1)
    r <- resample(r, new_raster, method = "ngb")

    # Stack
    coords <- stack(coords, r)
}


## Add nitrogen fixation rates
nfix_luo <- read.csv("/home/deriksson/Projects/Notion_DE/Data/n2_fix_Luo_et_al_2014/Copy_to_work_with/Estimate_global_N2fix_luo_2014_copy.csv", sep = ";")
names(nfix_luo)[3] <- "nfix_luo"
# Convert to raster
r <- rasterFromXYZ(nfix_luo)
# We miss some cells, so we will resample our values into a new created complete (180, 360, 64800) raster object
new_raster <- raster(nrow = 180, ncol = 360, xmn = -180, xmx = 180, ymn = -90, ymx = 90)
res(new_raster) <- c(1, 1)
r <- resample(r, new_raster, method = "ngb")
coords <- stack(coords, r)

#
r <- raster("/home/deriksson/Projects/Notion_DE/Data/n2fix_Wang_et_al_2019/1_Output/depth_integrated_n2_fix_top_2_layers.grd")
names(r) <- "nfix_wang"

# Wang has some negative values, we will replace those with zero --> we account for nitrogen fixation which we don't consider to be negative 
r[r < 0] <- 0.001
coords <- stack(coords, r)


# Add in situ nitrogen fixation rates
## Now we use nitrogen fixation rate observations - Let's see if we find a relationship between measured 
## marine nitrogen fixation rates
# Libraries
library(tidyverse)
library(readr)

# Load nitrogen fixation rates
luo <- readr::read_delim("/home/deriksson/Projects/Notion_DE/Data/n2_fix_Luo_et_al_2014/Copy_to_work_with/dat_N2_fixation.csv", delim = ";")
luo <- luo[, c(5, 6, 12)]
df_nfix <- data.frame(
  lon = as.numeric(unlist(luo$LONGITUDE)),
  lat = as.numeric(unlist(luo$LATITUDE)),
  nfix = as.numeric(unlist(luo[, "Total N2 Fixation Rates (μmol N m-3 d-1)"])),
  unit = "μmol N m-3 d-1",
  set = "Luo"
)
#
bonnet <- read.csv("/home/deriksson/Projects/Notion_DE/Data/n2fix_Bonnet_et_al_2017/1_Output/dat_N2_fixation.csv")
bonnet <- bonnet[, c(4, 5, 7:9)]
bonnet <- data.frame(
  lon = bonnet$Longitude,
  lat = bonnet$Latitude,
  nfix = rowMeans(bonnet[, 3:5]),
  unit = "μmol N m-3 d-1",
  set = "bonnet"
)
df_nfix <- rbind(df_nfix, bonnet)
df_nfix <- df_nfix[, 1:3]
df_nfix$x <- round(df_nfix$lon + 0.5) - 0.5
df_nfix$y <- round(df_nfix$lat + 0.5) - 0.5
#
df_nfix <- data.frame(
  x = df_nfix$x,
  y = df_nfix$y,
  nfix = df_nfix$nfix
)


## Load nitrogen fixation rates Luo et al. 2023 - Although measured for the whole water column, we assume majority of N2 fixation is in surface anyway (Shao 2023 Fig. 1)
df_int <- read.csv("/home/deriksson/Projects/Notion_DE/Data/Zhibo_Shao_2023/Copy_to_work_with/N2Fixation_Integral.csv", sep = ";")
df_vol <- read.csv("/home/deriksson/Projects/Notion_DE/Data/Zhibo_Shao_2023/Copy_to_work_with/N2Fixation_volumetric.csv", sep = ";")

# Plot 
names(df_int)


df1 <- df_int[, c("LATITUDE", "LONGITUDE", "Total.N2.Fixation..μmol.N.m.2.d.1.")]
names(df1) <- c("y", "x", "nfix")
df1 <- df1[, c("x", "y", "nfix")]
#
df2 <- df_vol[, c("LATITUDE", "LONGITUDE", "Total.N2.Fixation..μmol.N.m.3.d.1.")]
names(df2) <- c("y", "x", "nfix")
df2 <- df2[, c("x", "y", "nfix")]

# Merge with nitrogen fixation rates
df_nfix <- rbind(df_nfix, df1, df2)
df_nfix$x <- round(df_nfix$x + 0.5) - 0.5
df_nfix$y <- round(df_nfix$y + 0.5) - 0.5


# If one than one measurement falls in same grid cell use mean
df_nfix$id <- paste0(df_nfix$x, "_", df_nfix$y)
df_nfix$nfix <- as.numeric(df_nfix$nfix)
df_nfix <- aggregate( .~id, data = df_nfix, FUN = mean,  na.rm = TRUE )
# Remove id column
df_nfix <- df_nfix[, -1]

# Convert to raster
r <- rasterFromXYZ( df_nfix )
# Remove negative nitrogen fixation rates
r[r$nfix < 0] <- NA

# Extent to same size
new_raster <- raster(nrow = 180, ncol = 360, xmn = -180, xmx = 180, ymn = -90, ymx = 90)
res(new_raster) <- c(1, 1)
r <- resample(r, new_raster, method = "ngb")

# Merge
coords <- stack(coords, r$nfix)


## Add k means cluster based on PA based individual annual presences
df <- read.csv("/home/deriksson/Projects/Notion_DE/Code/10_Niche_analysis/2_Output/KMeans_clustering_PABased_AnnualIndividualEnsembles.csv")
df <- df[, -1]
names(df)[3:ncol(df)] <- paste0( "Cluster_kMeans", 2:10 )
#
r <- rasterFromXYZ(df)
# Extent to same size
new_raster <- raster(nrow = 180, ncol = 360, xmn = -180, xmx = 180, ymn = -90, ymx = 90)
res(new_raster) <- c(1, 1)
r <- resample(r, new_raster, method = "ngb")
coords <- stack(coords, r)

# Save final table
fln <- paste0( wd_out, "FinalTable_AllMetrics.grd" )
writeRaster(coords, fln, overwrite = TRUE)
