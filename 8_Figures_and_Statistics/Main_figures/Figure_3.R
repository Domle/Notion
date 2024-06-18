### In this script we plot main figure 3 from the Manuscript "Nitrogen fixation increases with diazotroph richness in the global ocean"

### NOTE: 	When running the code, please always double check input and output directories, to ensure you know where exactly files get saved so they
###			    can be referred as the correct input files within each consecutive script.

# Author: 	
#     Dominic Eriksson
#			Environmental Physics Group, UP
# 		ETH Zurich
#			Switzerland


# Input files:
#       1. shapefile goas

# Output files:
#       1. Raster stack as .grd file, later used for final plotting

# derikssong@ethz.ch; 18th of June 2024




## Load libraries
library(raster)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(betapart)
library(sf)
library(rnaturalearth)
library(terra)
library(tidyterra)
library(ggpubr)
library(hexbin)

## This script plot all plots in the manuscript and additionally used anywhere else. 
dark_green_hex <- "#2ca25f" # Cyanobacterial diazotroph
dark_blue_hex <- "#c51b8a" # Non-cyanobacterial diazotroph community
orange_hex <- "#f1a340" # Total diazotroph community


# Load functions
source("/home/deriksson/Projects/Notion_DE/Code/Functions/Dominic/subset_ocean_regions.R")


# Load Table
r <- brick("/home/deriksson/Projects/Notion_DE/Code/Tables/1_Output/FinalTable_AllMetrics.grd")


### ================================================================
### PREPARATIONS FOR PACIFIC CENTERED MAP
### ================================================================
## Required data & data formatting
# Read shapefile with geometries of ocean basins
sf.ocean <- st_read("/net/kryo/work/deriksson/Projects/Notion_DE/Data/Ocean_regions_shapefile/GOaS_v1_20211214/goas_v01.shp")
# Pacific centered world map: Source: https://stackoverflow.com/questions/56146735/visual-bug-when-changing-robinson-projections-central-meridian-with-ggplot2
worldMap <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_make_valid()
# Set projection pacific centered
lon_pacific <- '200' 
target_crs <- st_crs( paste0("+proj=eqc +x_0=0 +y_0=0 +lat_0=0 +lon_0=", lon_pacific) )

# define a long & slim polygon that overlaps the meridian line & set its CRS to match
# that of world
# Centered in lon 200
offset <- 180 - as.numeric(lon_pacific)


polygon <- st_polygon(x = list(rbind(
  c(-0.0001 - offset, 90),
  c(0 - offset, 90),
  c(0 - offset, -90),
  c(-0.0001 - offset, -90),
  c(-0.0001 - offset, 90)
))) %>%
  st_sfc() %>%
  st_set_crs(4326)


# modify world dataset to remove overlapping portions with world's polygons
world2 <- worldMap %>% st_difference(polygon)
#> Warning: attribute variables are assumed to be spatially constant throughout all
#> geometries
# Transform
sf.ocean_pacific <- world2 %>% st_transform(crs = target_crs)





#################################################################################################
## GLOBAL MAP - Ensemble BetaRatio (PA Based)
#################################################################################################
sort(names(r))
## Total community
dat <- as.data.frame(r$Diversity_BetaRatio_ensembleMean_paBased, xy = TRUE)
dat <- rast(dat)
# Set reference system
crs(dat) <- 'epsg:4326'

# Calculate upper bound 70th percentile of coefficient of variation
df <- cbind(
  as.data.frame(r$Diversity_BetaRatio_ensembleMean_paBased, xy = TRUE),
  as.data.frame(r$Diversity_BetaRatio_ensembleSD_paBased)
)
df$coefficient_of_variation <- df$Diversity_BetaRatio_ensembleSD_paBased/df$Diversity_BetaRatio_ensembleMean_paBased
percentile <- quantile(df$coefficient_of_variation, probs = 0.7, na.rm = TRUE) 
df1 <- df[which(!is.na(df$coefficient_of_variation)), ]
df1 <- df1[which(df1$coefficient_of_variation > percentile  ), ]
hb  <- hexbin::erode(hexbin(df1$x, df1$y, xbins = 50))
df70 <- as.data.frame(hcell2xy(hb))
# Convert stippling points to the Pacific-centered CRS
stipple_points_sf <- st_as_sf(df70, coords = c("x", "y"), crs = 4326)
stipple_points_sf <- st_transform(stipple_points_sf, crs = target_crs)

# Plot using ggplot2
gg_diversity <- 
  
  # Plot world map
  ggplot(sf.ocean_pacific) + 
  geom_sf(fill = "grey") +
  
  # Plot countour plot & and contour line settings
  geom_spatraster_contour_filled(data = dat, breaks = seq(0, 1, 0.1), show.legend = F) + 
  geom_spatraster_contour(
    data = dat,
    color = "grey",
    linewidth = 0.1) +

  # Extend legend values to the break values above, froom 0 to 1
  scale_fill_viridis_d(drop = FALSE) +


  # Add stippling upper 70th percentile bound
  geom_sf(data = stipple_points_sf, color = "white", alpha = 0.5, size = 0.5) +

  # Legend settings
  # guides(fill = guide_legend(nrow = 1, title = "Ensemble species turnover", title.position = "top")) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.title.align = 0.5,
    legend.key.size = unit(10, 'cm'), # change legend key size
    legend.key.height = unit(0.3, 'cm'), # change legend key height
    legend.key.width = unit(0.3, 'cm'),
    legend.title = element_text(size = 10), # change legend title font size
    legend.text = element_text(size = 8),

    # Plot background color
    panel.background = element_rect(fill = "white")) +

    # Axis labels
    xlab("") +
    ylab("") 


# Save plot
fln <- paste0( "/home/deriksson/Projects/Notion_DE/Code/Tables/Visualizations/", "GlobalMap_BetaRatio_totalDiaz_annualEnsembleMean.png" )
ggsave(
  filename = fln,
  plot = gg_diversity,
  # width = 
  # height = 
  dpi = 300
)


