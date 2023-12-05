## This script computes Figure 1 from Dominic Eriksson et al.(submitted)


#   Input files:
#       1.  Shapefile of the world
#       2.  Grid file containing the raster stack of the global annual ensemble richness. The layers of the raster stack are different statistics
#           such as mean, standard deviation, etc.
#       3.  RData files containing the annual latitudinal richness gradients for each model.
#       4.  Grid file containing the components of betadiversity ensemble

#   Output files:
#       1. One png file for each panel containing the plot.

# Strategy:         Pacific Center Map: We will reformat a world shapefile to be pacific centered. This reformated worldmap will then be used for plotting. The 
#                   methodology to pacific center the map can be found here: 
#                   https://stackoverflow.com/questions/56146735/visual-bug-when-changing-robinson-projections-central-meridian-with-ggplot2
#      
#                   For Figure 1 panel A (global richness map) we use the annual ensemble richness based on the total dataset as input to the Species Distribution Model, and calculated from presence 
#                   absence converted projection outputs. Therefore, the global map in Panel A shows the average number of annual presences across all ensemble members. We use
#                   a white stippling to highlight marine regions with higher degrees of uncertainties based on the coefficient of variation (stippling based 
#                   on 70th percentile threshold).

#                   For Figure 1 panel B (latitduinal binned mean line plots) we load all files containing latitudinal data and merge it into one dataframe. We aggregate
#                   using the mean function and use the color option from ggplot to color different models.

#                   For Figure 1 panel C and D (Global species turnover and nestedness map) we create legends manually and save them. 

#       

# Figure Caption:   Global diazotroph ensemble diversity. A) Global annual ensemble mean of diazotroph species richness. 
#                   B) Global 1° binned latitudinal richness gradients for each one of the 18 models (colored lines) 
#                   the ensemble has been generated from. The black line is the mean across all 18 models. White stipples 
#                   indicate areas where the coefficient of variation was above the 70th percentile marking greater differences 
#                   between model projections. C) Global annual ensemble species turnover. D) Global ensemble nestedness. 
#                   Species turnover measures the degree of species replacement and nestedness observes species assemblages 
#                   that are smaller subsets of larger sets, based on Jaccard’s dissimilarity index.

# Author: Dominic Eriksson
#					Environmental Physics Group, UP
#					ETH, Zurich
#					Switzerland

# 9th of October 2023, deriksson@ethz.ch --------------------------------------------------------------------------


### ====================================================================
### Preparations
### ====================================================================

# Libraries
lib_vec <- c(
    "sf", # read shapefile
    "raster", # work with grid files and spatial raster
    "openxlsx", # read.xlsx function
    "terra", # Package for spatial data analysis, use the rast function to create SpatRaster object
    "tidyterra", # function geom_spatraster_contour_filled

    "ggplot2", # used for plotting

    "ggpubr", # calculate spearman correlation and plot it

    "gridExtra", # package to align multiple plots
    "metR", # package to do the stippling
    "hexbin", # package to do the stippling
    
    "polynom", # package to print polynomial function on plot
    "rnaturalearth"
)
sapply(lib_vec, library, character.only = TRUE)

# Directories
wd_out <- "/net/kryo/work/deriksson/Projects/Notion_DE/Figures/1_Output/"


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

### =========================================================
### Figure 1, Panel A: Global ensemble richness - Pacific Centered
### =========================================================

# Diazotroph ensemble mean richness & uncterainty (SD) grid file
prj.stack <- brick( "/home/deriksson/Projects/Notion_DE/Code/8_Diversity/3_Output/Final_ensemble_used/Ensemble_AcrossEveryModelandMonth.grd" ) 
df <- as.data.frame(prj.stack$mean_richness, xy = TRUE)
names(df) <- c("lon", "lat", "richness")
#
df <- cbind(df, as.data.frame(prj.stack$sd))
df$coefficient_of_variation <- (df$sd/df$richness) * 100

# Create spatRaster
dat <- rast(df[, c("lon", "lat", "richness")])
# Set reference system
crs(dat) <- 'epsg:4326'

# Calculate upper bound 70th percentile of coefficient of variation
quantile(df$coefficient_of_variation, probs = 0.7, na.rm = TRUE) # 0.33
df1 <- df[which(!is.na(df$richness)), ]
df1 <- df1[which(df1$coefficient_of_variation > 45  ), ]
hb  <- erode(hexbin(df1$lon, df1$lat, xbins = 50))
df70 <- as.data.frame(hcell2xy(hb))
# Convert stippling points to the Pacific-centered CRS
# stipple_points_sf <- st_as_sf(stipple_points, coords = c("lon", "lat"), crs = 4326)
stipple_points_sf <- st_as_sf(df70, coords = c("x", "y"), crs = 4326)
# stipple_points_sf <- st_transform(stipple_points_sf, crs = pacific_crs)
stipple_points_sf <- st_transform(stipple_points_sf, crs = target_crs)


## Plot
gg_diversity <- 
  
  # Plot world map
  ggplot(sf.ocean_pacific) + 
  geom_sf(fill = "white") +
  
  # Plot countour plot & and contour line settings
  geom_spatraster_contour_filled(data = dat, breaks = seq(0, 1, 0.1), show.legend = FALSE) + 
  geom_spatraster_contour(
    data = dat, breaks = seq(0, 1, 0.1),
    # color = "grey",
    linewidth = 0.1) +

  # Add stippling upper 70th percentile bound
  geom_sf(data = stipple_points_sf, color = "white", alpha = 0.5, size = 0.5) +

  # Legend settings
    # labs(fill = 'Normalized richness') +
  scale_fill_viridis_d() +
  # guides(fill = guide_legend(nrow = 1, title = "Normalized diazotroph richness", title.position = "top")) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),

    # Plot background color
    panel.background = element_rect(fill = "darkgrey")) +

    # Axis labels
    xlab("") +
    ylab("") +

  # Panel labels
  labs(tag = "A")


# Save plot
ggsave(
  # plot = gg_diversity,
  filename = paste0(wd_out, "1A_EnsembleRichness_PacificCentered.png"),
  dpi = 300
)


### ===================================================================
### Figure 1, panel B: Latitudinal richness gradient
### ===================================================================

# Add latitudinal gradient, we plot each diversity gradient from each model in the ensemble
wd <- "/home/deriksson/Projects/Notion_DE/Code/8_Diversity/1_Output/Ensemble_mean/Total_dataset/lat_annually/pa/"

# Open
fnames <- list.files(wd)
fnames <- grep("normalized", fnames, value = TRUE)

# Load files into list
df_lat <- data.frame()
for(i in seq_along(fnames)){
  d <- data.frame( get(load(paste0(wd, fnames[i]))))[, -7]
  d$model <- as.factor(i)
  df_lat <- rbind(df_lat, d)
}

# Compute the mean latitudinal richness across the 18 mnodels
df <- aggregate(. ~ lat, data = df_lat, FUN = mean)

# Plot
gg_latGradient <- ggplot() +
  
  # Add the each ensemble model
  geom_line(data = df_lat, aes(x = lat, y = mean, color = model), show.legend = FALSE, size = 0.3, alpha = 0.6)+  
  # Add mean richness
  geom_line(data = df, aes(x = lat, y = mean), color = "black", show.legend = FALSE, size = 0.6) + 

  xlab("Latitude") +
  ylab("Normalized diazotroph richness") +
  theme(
    panel.background = element_rect(fill = "white"),
    axis.text.x = element_text(size = 32),
    axis.text.y = element_text(size = 32),
    axis.title.x = element_text(size = 32),
    axis.title.y = element_text(size = 32)
    ) + 
  xlim(-75, 75) +
  labs(tag = 'B')

# We use the coord_flip option to rotate the figure
gg_latGradient <- gg_latGradient + coord_flip()

# Save plot
ggsave(
  plot = gg_latGradient,
  filename = paste0(wd_out, "1B_RichnessLatitudinalGradient.png"),
  height = 7.51,
  width = 3,
  units = "cm",
  dpi = 300  
)


### ===================================================================
### 1.C Global ensemble Species Turnover
### ===================================================================

# Directories
wd_in <- "/home/deriksson/Projects/Notion_DE/Code/8_Diversity/6_Output/"

# Load data
r <- brick(paste0(wd_in, "TotalEnsemble_Betadiversity.grd"))

# Convert raster to spatRast object
dat <- rast(r$mean_jtu)
# Set reference system
crs(dat) <- 'epsg:4326'

# Calculate upper bound 70th percentile of coefficient of variation
df <- as.data.frame( (r$sd_jtu/r$mean_jtu * 100), xy = TRUE ) 
names(df)[3] <- c("coefficient_of_variation")
percentile <- quantile(df$coefficient_of_variation, probs = 0.7, na.rm = TRUE) # 0.33
df1 <- df[which(!is.na(df$coefficient_of_variation)), ]
df1 <- df1[which(df1$coefficient_of_variation > percentile  ), ]
hb  <- erode(hexbin(df1$x, df1$y, xbins = 50))
df70 <- as.data.frame(hcell2xy(hb))
# Convert stippling points to the Pacific-centered CRS
# stipple_points_sf <- st_as_sf(stipple_points, coords = c("lon", "lat"), crs = 4326)
stipple_points_sf <- st_as_sf(df70, coords = c("x", "y"), crs = 4326)
# stipple_points_sf <- st_transform(stipple_points_sf, crs = pacific_crs)
stipple_points_sf <- st_transform(stipple_points_sf, crs = target_crs)

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

  # Add stippling upper 70th percentile bound
  geom_sf(data = stipple_points_sf, color = "white", alpha = 0.5, size = 0.5) +

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
    labs(tag = 'C')

# Save plot
ggsave(
  plot = gg_diversity,
  filename = paste0(wd_out, "1C_EnsembleTurnover_pacificCentered.png"),
  dpi = 300
)


# Create legend separately
plot <- 
  ggplot(dat, aes(x, y, fill = mean_jtu)) +
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
fln <- paste0( wd_out, 'Legend_1C.png' )
ggsave(
  plot = legend,
  filename = fln,
  dpi = 300
)  



### ===================================================================
### Figure 1 panel D: Global ensemble Nestedness
### ===================================================================

# Convert to spatRast object
dat <- rast(r$mean_jne)
# Set reference system
crs(dat) <- 'epsg:4326'

# Calculate upper bound 70th percentile of coefficient of variation
df <- as.data.frame( (r$sd_jne/r$mean_jne * 100), xy = TRUE ) 
names(df)[3] <- c("coefficient_of_variation")
percentile <- quantile(df$coefficient_of_variation, probs = 0.7, na.rm = TRUE) # 0.33
df1 <- df[which(!is.na(df$coefficient_of_variation)), ]
df1 <- df1[which(df1$coefficient_of_variation > percentile  ), ]
hb  <- erode(hexbin(df1$x, df1$y, xbins = 50))
df70 <- as.data.frame(hcell2xy(hb))
# Convert stippling points to the Pacific-centered CRS
# stipple_points_sf <- st_as_sf(stipple_points, coords = c("lon", "lat"), crs = 4326)
stipple_points_sf <- st_as_sf(df70, coords = c("x", "y"), crs = 4326)
# stipple_points_sf <- st_transform(stipple_points_sf, crs = pacific_crs)
stipple_points_sf <- st_transform(stipple_points_sf, crs = target_crs)

# Plot using ggplot2
gg_diversity <- 
  
  # Plot world map
  ggplot(sf.ocean_pacific) + 
  geom_sf(aes(fill = name),
    fill = "white") +
  
  # Plot countour plot & and contour line settings
  geom_spatraster_contour_filled(data = dat, show.legend = TRUE) + 
  geom_spatraster_contour(
    data = dat,
    # color = "grey",
    linewidth = 0.1) +

  # Add st.1ippling upper 70th percentile bound
  geom_sf(data = stipple_points_sf, color = "white", alpha = 0.5, size = 0.5) +

  # Legend settings
  scale_fill_viridis_d(guide = guide_legend()) +
  guides(fill = guide_legend(nrow = 1, title = "Ensemble species nestedness", title.position = "top")) +
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
    labs(tag = 'D')

# Save plot
ggsave(
  plot = gg_diversity,
  filename = paste0(wd_out, "1D_EnsembleNestedness_pacificCentered.png"),
  dpi = 300
)


# Create legend separately
plot <- 
  ggplot(dat, aes(x, y, fill = mean_jne)) +
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
fln <- paste0( wd_out, 'Legend_1D.png' )
ggsave(
  plot = legend,
  filename = fln,
  dpi = 300
)  


###==============================================================
### END
###==============================================================