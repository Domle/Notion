### In this script we plot main figure 1 from the Manuscript "Nitrogen fixation increases with diazotroph richness in the global ocean"

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





#################################################################################################
## GLOBAL MAP - Ensemble Jaccard Dissimilarity (PA Based) - ANNUAL
#################################################################################################
library(hexbin)
## Total community -------------------------------------------------------------------------------
sort(names(r))
#
dat <- as.data.frame(r$Diversity_Jaccard_ensembleMean_paBased, xy = TRUE)
dat <- rast(dat)
# Set reference system
crs(dat) <- 'epsg:4326'

# Calculate upper bound 70th percentile of coefficient of variation

df <- cbind(
  as.data.frame(r$Diversity_Jaccard_ensembleMean_paBased, xy = TRUE),
  as.data.frame(r$Diversity_Jaccard_ensembleSD_paBased)
)
df$coefficient_of_variation <- df$Diversity_Jaccard_ensembleSD_paBased/df$Diversity_Jaccard_ensembleMean_paBased
percentile <- quantile(df$coefficient_of_variation, probs = 0.7, na.rm = TRUE) # 0.33
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
  geom_spatraster_contour_filled(data = dat, breaks = seq(0, 1, 0.05), show.legend = TRUE) + 
  geom_spatraster_contour(
    data = dat,
    color = "grey",
    linewidth = 0.1) +

  # Extend legend values to the break values above, froom 0 to 1
  scale_fill_viridis_d(drop = FALSE) +


  # Add stippling upper 70th percentile bound
  geom_sf(data = stipple_points_sf, color = "white", alpha = 0.5, size = 0.5) +

  # Legend settings
  # guides(fill = guide_legend(nrow = 1, title = "Ensemble Jaccard dissimilarity", title.position = "top")) +
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
fln <- paste0( "/home/deriksson/Projects/Notion_DE/Code/Tables/Visualizations/", "GlobalMap_Jaccard_totalDiaz_annualEnsembleMean.png" )
ggsave(
  filename = fln,
  plot = gg_diversity,
  # width = 
  # height = 
  dpi = 300
)





#################################################################################################
## GLOBAL MAP - Ensemble Species Turnover (PA Based)
#################################################################################################
sort(names(r))
## Total community
dat <- as.data.frame(r$Diversity_Turnover_ensembleMean_paBased, xy = TRUE)
dat <- rast(dat)
# Set reference system
crs(dat) <- 'epsg:4326'

# Calculate upper bound 70th percentile of coefficient of variation
df <- cbind(
  as.data.frame(r$Diversity_Turnover_ensembleMean_paBased, xy = TRUE),
  as.data.frame(r$Diversity_Turnover_ensembleSD_paBased)
)
df$coefficient_of_variation <- df$Diversity_Turnover_ensembleSD_paBased/df$Diversity_Turnover_ensembleMean_paBased
percentile <- quantile(df$coefficient_of_variation, probs = 0.7, na.rm = TRUE) 
df1 <- df[which(!is.na(df$coefficient_of_variation)), ]
df1 <- df1[which(df1$coefficient_of_variation > percentile  ), ]
hb  <- erode(hexbin(df1$x, df1$y, xbins = 50))
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
fln <- paste0( "/home/deriksson/Projects/Notion_DE/Code/Tables/Visualizations/", "GlobalMap_SpeciesTurnover_totalDiaz_annualEnsembleMean.png" )
ggsave(
  filename = fln,
  plot = gg_diversity,
  # width = 
  # height = 
  dpi = 300
)





#################################################################################################
## GLOBAL MAP - Ensemble Nestedness (PA Based)
#################################################################################################
sort(names(r))
## Total community
dat <- as.data.frame(r$Diversity_Nestedness_ensembleMean_paBased, xy = TRUE)
dat <- rast(dat)
# Set reference system
crs(dat) <- 'epsg:4326'

# Calculate upper bound 70th percentile of coefficient of variation
df <- cbind(
  as.data.frame(r$Diversity_Nestedness_ensembleMean_paBased, xy = TRUE),
  as.data.frame(r$Diversity_Nestedness_ensembleSD_paBased)
)
df$coefficient_of_variation <- df$Diversity_Nestedness_ensembleSD_paBased/df$Diversity_Nestedness_ensembleMean_paBased
percentile <- quantile(df$coefficient_of_variation, probs = 0.7, na.rm = TRUE) 
df1 <- df[which(!is.na(df$coefficient_of_variation)), ]
df1 <- df1[which(df1$coefficient_of_variation > percentile  ), ]
hb  <- erode(hexbin(df1$x, df1$y, xbins = 50))
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
fln <- paste0( "/home/deriksson/Projects/Notion_DE/Code/Tables/Visualizations/", "GlobalMap_Nestedness_totalDiaz_annualEnsembleMean.png" )
ggsave(
  filename = fln,
  plot = gg_diversity,
  # width = 
  # height = 
  dpi = 300
)