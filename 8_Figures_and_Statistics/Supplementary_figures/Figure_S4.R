

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
## BOXPLOTS and LINEPLOTS RICHNESS ACROSS ENSEMBLE MEMBERS
#################################################################################################

## This script plot all plots in the manuscript and additionally used anywhere else. 
dark_green_hex <- "#2ca25f" # Cyanobacterial diazotroph
purple_hex <- "#c51b8a" # Non-cyanobacterial diazotroph community
orange_hex <- "#f1a340" # Total diazotroph community

## Total community --------------------------------------------
r <- brick("/home/deriksson/Projects/Notion_DE/Code/Diversity/Diversities/pa/RichnessNormalized_annual_totalDiaz_Ensemble_and_EnsembleMembers.grd")
r <- r[[ 3:nlayers(r) ]]
names(r) <- paste0( "SDM", 1:90 )

# Convert to dataframe
df <- as.data.frame(r)
df <- na.omit(df)
# From wide to long foramt
df <- gather(df, SDM, richness, SDM1:SDM90, factor_key = TRUE)


# Plot boxplot
gg_boxplot <- ggplot( data = df ) +
  geom_boxplot( aes(x = SDM, y = richness), fill = orange_hex ) +
  theme(axis.text.x = element_text(size = 32, angle = 90, vjust = 0.5),
    axis.text.y = element_text(size = 32),
    axis.title.x = element_text(size = 32),
    axis.title.y = element_text(size = 32)
    ) +
  xlab( "Ensemble members" ) + 
  ylab( "Annual richness" ) 


# Save
wd_out <- "/home/deriksson/Projects/Notion_DE/Code/Tables/Visualizations/"
fln <- paste0(wd_out, "Boxplot_totalDiaz_medianRichness_across_ensembleMembers.png")
ggsave(
  filename = fln,
  plot = gg_boxplot,
  width = 40,
  height = 20,
  bg = "white",
  dpi = 300
)

## cyanos --------------------------------------------
r <- brick("/home/deriksson/Projects/Notion_DE/Code/Diversity/Diversities/pa/RichnessNormalized_annual_cyanoDiaz_Ensemble_and_EnsembleMembers.grd")
r <- r[[ 3:nlayers(r) ]]
names(r) <- paste0( "SDM", 1:90 )

# Convert to dataframe
df <- as.data.frame(r)
df <- na.omit(df)
# From wide to long foramt
df <- gather(df, SDM, richness, SDM1:SDM90, factor_key = TRUE)


# Plot boxplot
gg_boxplot <- ggplot( data = df ) +
  geom_boxplot( aes(x = SDM, y = richness), fill = dark_green_hex ) +
  theme(axis.text.x = element_text(size = 32, angle = 90, vjust = 0.5),
    axis.text.y = element_text(size = 32),
    axis.title.x = element_text(size = 32),
    axis.title.y = element_text(size = 32)
    ) +
  xlab( "Ensemble members" ) + 
  ylab( "Annual richness" ) 


# Save
wd_out <- "/home/deriksson/Projects/Notion_DE/Code/Tables/Visualizations/"
fln <- paste0(wd_out, "Boxplot_cyanoDiaz_medianRichness_across_ensembleMembers.png")
ggsave(
  filename = fln,
  plot = gg_boxplot,
  width = 40,
  height = 20,
  bg = "white",
  dpi = 300
)


## non-cyanos --------------------------------------------
r <- brick("/home/deriksson/Projects/Notion_DE/Code/Diversity/Diversities/pa/RichnessNormalized_annual_nonCyanoDiaz_Ensemble_and_EnsembleMembers.grd")
r <- r[[ 3:nlayers(r) ]]
names(r) <- paste0( "SDM", 1:90 )

# Convert to dataframe
df <- as.data.frame(r)
df <- na.omit(df)
# From wide to long foramt
df <- gather(df, SDM, richness, SDM1:SDM90, factor_key = TRUE)


# Plot boxplot
gg_boxplot <- ggplot( data = df ) +
  geom_boxplot( aes(x = SDM, y = richness), fill = purple_hex) +
  theme(axis.text.x = element_text(size = 32, angle = 90, vjust = 0.5),
    axis.text.y = element_text(size = 32),
    axis.title.x = element_text(size = 32),
    axis.title.y = element_text(size = 32)
    ) +
  xlab( "Ensemble members" ) + 
  ylab( "Yearly averaged diazotroph richness" ) 


# Save
wd_out <- "/home/deriksson/Projects/Notion_DE/Code/Tables/Visualizations/"
fln <- paste0(wd_out, "Boxplot_nonCyanoDiaz_medianRichness_across_ensembleMembers.png")
ggsave(
  filename = fln,
  plot = gg_boxplot,
  width = 40,
  height = 20,
  bg = "white",
  dpi = 300
)
