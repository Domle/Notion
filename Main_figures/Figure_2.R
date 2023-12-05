## This script computes Figure 2 from Dominic Eriksson et al.(submitted). Global map on 
## annual ensemble richness from cyanobacterial diazotrophs and non-cyanobacterial diazotrophs.


#   Input files:
#           1.  Grid file with ensemble richness on annual basis from cyanobacterial and non-cyanobacteria
#               diazotrophs.


#   Output files:
#           2. PNG files for each panel.


# Strategy: We create a pacific centered world map and plot using ggplot2 package. We save legends separately.


# Figure Caption:   Annual ensemble mean species richness of A) cyanobacterial and B) non-cyanobacterial 
#                   diazotrophs. The ensemble (n = 18) has been computed from 1° longitudinal and 1° 
#                   latitudinal monthly Species Distribution Model outputs that account for uncertainties 
#                   related to predictors, algorithms and background selection strategies chosen. Ensemble 
#                   richness has been normalized by the number of species modeled with blue colors indicating 
#                   low richness and the yellowish colors increasing the annual ensemble richness.

# Author:   Dominic Eriksson
#		    Environmental Physics Group, UP
#			ETH, Zurich
#			Switzerland

# 9th of October 2023, deriksson@ethz.ch --------------------------------------------------------------------------



### ====================================================================
### Preparations
### ====================================================================
# Clear workspace
rm(list = ls())

# Libraries
library(raster)

# Load functions
source('/home/deriksson/Projects/Notion_DE/Code/Functions/Dominic/global_map.R')

# Load data
r <- brick('/net/kryo/work/deriksson/Projects/Notion_DE/Code/8_Diversity/3_Output/Ensemble_mean/Total_dataset/div_annually/pa/Ensemble_richness_acrossAllModels_CyanosVsNonCyanos.grd')


### ====================================================================
### 1. Pacific Center Worldmap
### ====================================================================
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

### ====================================================================
### Figure 2 panel A: Cyanobacterial diazotroph annual ensemble richness
### ====================================================================

# Convert raster to spatRaster
dat <- rast(r$Total_ensemble_mean)
# Set reference system
crs(dat) <- 'epsg:4326'

## Plot
# Plot using ggplot2
gg_diversity <- 
  
  # Plot world map
  ggplot(sf.ocean_pacific) + 
  geom_sf(fill = "white") +
  
  # Plot countour plot & and contour line settings
  geom_spatraster_contour_filled(data = dat, show.legend = TRUE) +
  geom_spatraster_contour(
    data = dat,
    color = "grey",
    linewidth = 0.1) +

  # Legend settings
  scale_fill_viridis_d(guide = guide_legend()) +
  guides(fill = guide_legend(nrow = 1, title = "Ensemble species turnover", title.position = "top")) +
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
    panel.background = element_rect(fill = "darkgrey")) +

    # Axis labels
    xlab("") +
    ylab("") +
    labs(tag = 'A')

## Save plot
fln <- paste0( wd_out, '2A_Richness_Cyanos_pacificCentered.png' )
ggsave(
  plot = gg_diversity,
  filename = fln ,
  dpi = 300
)


# Create legend separately
d <- as.data.frame(dat, xy = TRUE)
plot <- 
  ggplot(d, aes(x, y, fill = Total_ensemble_mean)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme(
    legend.position = 'bottom',
    legend.key.width = unit(4, 'cm'),
    legend.key.height = unit(0.8, 'cm'),
    legend.text = element_text(size = 22)
  ) + 
  labs(fill = '')

# Extract legend
legend <- get_legend(plot)


# Save legend
fln <- paste0( wd_out, 'Legend_2A.png' )
ggsave(
  plot = legend,
  filename = fln,
  dpi = 300
)  


### ====================================================================
### Figure 2 panel B: Non-cyanobacterial diazotroph annual ensemble richness
### ====================================================================

# Convert raster to spatRaster
dat <- rast(r$layer)
# Set reference system
crs(dat) <- 'epsg:4326'

# Plot using ggplot2
gg_diversity <- 
  
  # Plot world map
  ggplot(sf.ocean_pacific) + 
  geom_sf(fill = "white") +
  
  # Plot countour plot & and contour line settings
  geom_spatraster_contour_filled(data = dat, show.legend = TRUE) +
    geom_spatraster_contour(
    data = dat,
    color = "grey",
    linewidth = 0.1) +

  # Legend settings
  scale_fill_viridis_d(guide = guide_legend()) +
  guides(fill = guide_legend(nrow = 1, title = "Ensemble species turnover", title.position = "top")) +
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
    panel.background = element_rect(fill = "darkgrey")) +

    # Axis labels
    xlab("") +
    ylab("") +
    labs(tag = 'B')

## Save plot
fln <- paste0( wd_out, '2B_Richness_NonCyanos_pacificCentered.png' )
ggsave(
  plot = gg_diversity,
  filename = fln ,
  dpi = 300
)


# Create legend separately
d <- as.data.frame(dat, xy = TRUE)
plot <- 
  ggplot(d, aes(x, y, fill = layer)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme(
    legend.position = 'bottom',
    legend.key.width = unit(4, 'cm'),
    legend.key.height = unit(0.8, 'cm'),
    legend.text = element_text(size = 22)
  ) + 
  labs(fill = '')

# Extract legend
legend <- get_legend(plot)


# Save legend
fln <- paste0( wd_out, 'Legend_2B.png' )
ggsave(
  plot = legend,
  filename = fln,
  dpi = 300
)  


###==============================================================
### END
###==============================================================