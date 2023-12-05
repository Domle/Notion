## In this script we calculate the final ensemble for betadiversity (turnover and nestedness) across all models
## and all months to get on annual mean global estimate. We compute the mean and standard deviation.

# Input files:
#       1.  CSV files containing longitude, latitude, jaccard dissimilarity, species turnover and nestedness. One file represents
#           one model and one month.

# Output files:
#       2.  A grid file (containing a raster stack object) containing the mean (across month and models) and standard deviation of jaccard dissimilarity, species turnover and
#           nestedness as raster layers.

# Strategy: We load each file and merge everything into one dataframe. Afterwards we aggregate using the mean function based on
#           longitude and latitude.

## Author:  Dominic Eriksson
##          Environmental Physics Group,  UP
##          ETH, Zurich
##          Switzerland

# deriksson@ethz.ch, 23rd of June 2023 -----------------------------------------------------------------


# Clean workspace
rm(list = ls())

# Libraries
library(raster)
library(ggplot2)

# Directories
wd_in <- "/net/kryo/work/deriksson/Projects/Notion_DE/Code/8_Diversity/5_Output/"
wd_out <- "/net/kryo/work/deriksson/Projects/Notion_DE/Code/8_Diversity/6_Output/"

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
names(df_mean)[3:5] <- c("mean_jac", "mean_jtu", "mean_jne")
df_sd <- aggregate( . ~ x + y, data = df, FUN = sd, na.rm = TRUE)
names(df_sd)[3:5] <- c("sd_jac", "sd_jtu", "sd_jne")

# Save file as raster
raster_stack <- stack( 
    rasterFromXYZ(df_mean[, c("x", "y", "mean_jac")]),
    rasterFromXYZ(df_mean[, c("x", "y", "mean_jtu")]),
    rasterFromXYZ(df_mean[, c("x", "y", "mean_jne")]),
    rasterFromXYZ(df_sd[, c("x", "y", "sd_jac")]),
    rasterFromXYZ(df_sd[, c("x", "y", "sd_jtu")]),
    rasterFromXYZ(df_sd[, c("x", "y", "sd_jne")])
 )
 
# Save as raster
fln <- paste0( wd_out, 'TotalEnsemble_Betadiversity.grd')
writeRaster( raster_stack, filename = fln, format = 'raster', overwrite= TRUE )


###==============================================================
### END
###==============================================================