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



### Plot latitudinal SD
## Latitudinal ensemble spread annual richness WITH RECENT DATA --------------------------------------------------------------------
# total community --------------------------------------
library(matrixStats)
r <- brick( "/home/deriksson/Projects/Notion_DE/Code/Diversity/Diversities/pa/RichnessNormalized_annual_totalDiaz_Ensemble_and_EnsembleMembers.grd" )
# Convert to data frame
df <- as.data.frame( r, xy = TRUE )
df <- df[, -c(1, 3, 4)]
# Aggregate
df <- aggregate( .~y, data = df, FUN = mean )

# Wide to long
lat_ensemble_mean <- rowMeans( as.matrix(df[, -1]) )
lat_ensemble_sd <- matrixStats::rowSds( as.matrix(df[, -1]) )
df_ensemble <- data.frame(
  y = df$y,
  mean = lat_ensemble_mean,
  sd = lat_ensemble_sd
)

# Maximum mean north and south of the equator
t <- df_ensemble[which(df_ensemble$y > 0), ]
t[which(t$mean == max(t$mean)), ]
#
t <- df_ensemble[which(df_ensemble$y < 0), ]
t[which(t$mean == max(t$mean)), ]

# Richness drop at equator
t <- df_ensemble[which(df_ensemble$y < 28.5 & df_ensemble$y > -21.5), ]
t[which(t$mean == min(t$mean)), ]

names(df)

# Aggregate
test <- aggregate( .~y, data = df, mean )

# Calculate the sd for each row across all columns
df_sd <- data.frame(test$y)
# Assuming df_sd is your data frame
df_sd$sd <- as.vector(matrixStats::rowSds(as.matrix(test[, -1])))

# Plot
# Assuming df_sd is your data frame
gg_lineplot <- ggplot( data = df_ensemble ) +
  geom_path( aes(x = sd, y = y) ) +
  theme(axis.text.x = element_text(size = 32, angle = 90, vjust = 0.5),
    axis.text.y = element_text(size = 32),
    axis.title.x = element_text(size = 32),
    axis.title.y = element_text(size = 32)
    ) +
  xlab( "Ensemble standard deviation" ) + 
  ylab( "Latitude [°]" ) +
  scale_y_continuous(breaks = seq(min(df_ensemble$y), max(df_ensemble$y), by = 5))

# Save 
fln <- "/home/deriksson/Projects/Notion_DE/Code/Tables/Visualizations/"
fln <- paste0(fln, "Lineplot_Latitudinal_standardDeviation_totalDiaz.png")
#
ggsave(
  filename = fln,
  plot = gg_lineplot,
  width = 15,
  height = 25,
  bg = "white",
  dpi = 300
)

# cyano  --------------------------------------
library(matrixStats)
r <- brick( "/home/deriksson/Projects/Notion_DE/Code/Diversity/Diversities/pa/RichnessNormalized_annual_cyanoDiaz_Ensemble_and_EnsembleMembers.grd" )
# Convert to data frame
df <- as.data.frame( r, xy = TRUE )
df <- df[, -c(1, 3, 4)]
# Aggregate
df <- aggregate( .~y, data = df, FUN = mean )

# Wide to long
lat_ensemble_mean <- rowMeans( as.matrix(df[, -1]) )
lat_ensemble_sd <- matrixStats::rowSds( as.matrix(df[, -1]) )
df_ensemble <- data.frame(
  y = df$y,
  mean = lat_ensemble_mean,
  sd = lat_ensemble_sd
)

# Maximum mean north and south of the equator
t <- df_ensemble[which(df_ensemble$y > 0), ]
t[which(t$mean == max(t$mean)), ]
#
t <- df_ensemble[which(df_ensemble$y < 0), ]
t[which(t$mean == max(t$mean)), ]

# Richness drop at equator
t <- df_ensemble[which(df_ensemble$y < 28.5 & df_ensemble$y > -21.5), ]
t[which(t$mean == min(t$mean)), ]

names(df)

# Aggregate
test <- aggregate( .~y, data = df, mean )

# Calculate the sd for each row across all columns
df_sd <- data.frame(test$y)
# Assuming df_sd is your data frame
df_sd$sd <- as.vector(matrixStats::rowSds(as.matrix(test[, -1])))

# Plot
# Assuming df_sd is your data frame
gg_lineplot <- ggplot( data = df_ensemble ) +
  geom_path( aes(x = sd, y = y) ) +
  theme(axis.text.x = element_text(size = 32, angle = 90, vjust = 0.5),
    axis.text.y = element_text(size = 32),
    axis.title.x = element_text(size = 32),
    axis.title.y = element_text(size = 32)
    ) +
  xlab( "Ensemble standard deviation" ) + 
  ylab( "Latitude [°]" ) +
  scale_y_continuous(breaks = seq(min(df_ensemble$y), max(df_ensemble$y), by = 5))

# Save 
fln <- "/home/deriksson/Projects/Notion_DE/Code/Tables/Visualizations/"
fln <- paste0(fln, "Lineplot_Latitudinal_standardDeviation_cyanoDiaz.png")
#
ggsave(
  filename = fln,
  plot = gg_lineplot,
  width = 15,
  height = 25,
  bg = "white",
  dpi = 300
)


# non-cyano  --------------------------------------
library(matrixStats)
r <- brick( "/home/deriksson/Projects/Notion_DE/Code/Diversity/Diversities/pa/RichnessNormalized_annual_nonCyanoDiaz_Ensemble_and_EnsembleMembers.grd" )
# Convert to data frame
df <- as.data.frame( r, xy = TRUE )
df <- df[, -c(1, 3, 4)]
# Aggregate
df <- aggregate( .~y, data = df, FUN = mean )

# Wide to long
lat_ensemble_mean <- rowMeans( as.matrix(df[, -1]) )
lat_ensemble_sd <- matrixStats::rowSds( as.matrix(df[, -1]) )
df_ensemble <- data.frame(
  y = df$y,
  mean = lat_ensemble_mean,
  sd = lat_ensemble_sd
)

# Maximum mean north and south of the equator
t <- df_ensemble[which(df_ensemble$y > 0), ]
t[which(t$mean == max(t$mean)), ]
#
t <- df_ensemble[which(df_ensemble$y < 0), ]
t[which(t$mean == max(t$mean)), ]

# Richness drop at equator
t <- df_ensemble[which(df_ensemble$y < 28.5 & df_ensemble$y > -21.5), ]
t[which(t$mean == min(t$mean)), ]

names(df)

# Aggregate
test <- aggregate( .~y, data = df, mean )

# Calculate the sd for each row across all columns
df_sd <- data.frame(test$y)
# Assuming df_sd is your data frame
df_sd$sd <- as.vector(matrixStats::rowSds(as.matrix(test[, -1])))

# Plot
# Assuming df_sd is your data frame
gg_lineplot <- ggplot( data = df_ensemble ) +
  geom_path( aes(x = sd, y = y) ) +
  theme(axis.text.x = element_text(size = 32, angle = 90, vjust = 0.5),
    axis.text.y = element_text(size = 32),
    axis.title.x = element_text(size = 32),
    axis.title.y = element_text(size = 32)
    ) +
  xlab( "Ensemble standard deviation" ) + 
  ylab( "Latitude [°]" ) +
  scale_y_continuous(breaks = seq(min(df_ensemble$y), max(df_ensemble$y), by = 5))

# Save 
fln <- "/home/deriksson/Projects/Notion_DE/Code/Tables/Visualizations/"
fln <- paste0(fln, "Lineplot_Latitudinal_standardDeviation_nonCyanoDiaz.png")
#
ggsave(
  filename = fln,
  plot = gg_lineplot,
  width = 15,
  height = 25,
  bg = "white",
  dpi = 300
)