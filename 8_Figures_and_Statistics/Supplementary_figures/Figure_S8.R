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
## Total community --------------------------------------------
r <- brick("/home/deriksson/Projects/Notion_DE/Code/Ensembles/Diversity/Annual_BetaRatio_Ensemble_and_ensembleMembers.grd")
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
  ylab( "Beta ratio" ) 


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

# Stats
v_me <- aggregate( richness ~ SDM, data = df, FUN = median )
v_iqr <- aggregate( richness ~ SDM, data = df, FUN = IQR )
d <- data.frame(
  median = v_me[, 2],
  iqr = v_iqr[, 2]
)
# Get interquantile ranges
d[which( d$median == min(d$median) ), ]
d[which( d$median == max(d$median) ), ]

