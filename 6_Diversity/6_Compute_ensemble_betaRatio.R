## In this script we compute the ensemble betaRatio across all ensemble members and months. We compute mean
## and standard deviation.

## Input files:
##      1.  CSV files that contain the beta-ratio, longitude and latitude.

## Output files:
##      1.  A grid file containing a raster stack with ensemble betadiversity components as layers.

## Strategy: We open each file and merge into one dataframe. We aggregate based on longitude and latitude.

## Author:  Dominic Eriksson
##          Environmental Physics Group, UP
##          ETH Zürich, Zurich
##          Switzerland

## deriksson@ethz.ch; 9th of August 2023 ---------------------------------------------------------

# Clear workspace
rm(list = ls())

# Libraries
library(raster)
library(ggplot2)

# Directories
wd_in <- "/home/deriksson/Projects/Notion_DE/Code/8_Diversity/8_Output/"
wd_out <- "/home/deriksson/Projects/Notion_DE/Code/8_Diversity/9_Output/"

# Vectors
filenames <- list.files(wd_in)
filenames <- grep(".csv", filenames, value = TRUE)


# Load data
l <- list()
for(i in seq_along(filenames)){
    print(i)
    d <- read.csv(paste0(wd_in, filenames[i]))
    l[[i]] <- d
}
# Merge
df <- do.call("rbind", l)

# Aggregate 
df_mean <- aggregate( . ~ x + y, data = df, FUN = mean, na.rm = TRUE)
names(df_mean)[3:6] <- c("mean_jac", "mean_jtu", "mean_jne", "mean_betaRatio_jne")
df_sd <- aggregate( . ~ x + y, data = df, FUN = sd, na.rm = TRUE)
names(df_sd)[3:6] <- c("sd_jac", "sd_jtu", "sd_jne", "sd_betaRatio_jne")

# Save file as raster
raster_stack <- stack( 
    rasterFromXYZ(df_mean[, c("x", "y", "mean_jac")]),
    rasterFromXYZ(df_mean[, c("x", "y", "mean_jtu")]),
    rasterFromXYZ(df_mean[, c("x", "y", "mean_jne")]),
    rasterFromXYZ(df_mean[, c("x", "y", "mean_betaRatio_jne")]),
    rasterFromXYZ(df_sd[, c("x", "y", "sd_jac")]),
    rasterFromXYZ(df_sd[, c("x", "y", "sd_jtu")]),
    rasterFromXYZ(df_sd[, c("x", "y", "sd_jne")]),
    rasterFromXYZ(df_sd[, c("x", "y", "sd_betaRatio_jne")])
 )
 
# 
fln <- paste0( wd_out, "Ensemble_BetaRatio.grd" )
writeRaster(
    raster_stack,
    filename = fln,
    format = "raster",
    overwrite = TRUE)


###==============================================================
### END
###==============================================================