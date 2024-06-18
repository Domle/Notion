
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



### ============================================================
### HEATMAP - TSS SCORES ACROSS TAXA AND ENSEMBLE MEMBERS
### ============================================================
# Directories
wd_in <- "/home/deriksson/Projects/Notion_DE/Code/6_SDM_fit/3_Output_750/"
wd_out <- "/home/deriksson/Projects/Notion_DE/Code/Tables/Visualizations/"


# Vector of filenames
filenames <- list.files( paste0(wd_in) )

# Open file
l <- list()
for(i in seq_along(filenames)){

    # Print progress
    print( i )

    # Open file
    d <- read.csv( paste0(wd_in, filenames[i]) )
    d$id <- as.factor(rep(1:5, (nrow(d)/5) ))
    d$id <- paste0( d$id, filenames[i] )

    # Add to list
    l[[i]] <- d

} # end of loop across filenames

# Merge list into dataframe
library(plyr)
df <- do.call( "rbind.fill", l )

# Reorder taxon and vars based on their original order in the dataframe
df$taxon <- factor(df$taxon, levels = unique(df$taxon))
df$id <- factor(df$id, levels = unique(df$id))




#####
## Cyanobacterial community ------------------------
vec_cyanos <- c(
  "Atelocyanobacterium",
  "Calothrix",
  "Calothrix confervicola",
  "Richelia",
  "Richelia intracellularis",
  "Trichodesmium",
  "Trichodesmium erythraeum",
  "Trichodesmium thiebautii",
  "UCYN.A",
  "UCYN.A1",
  "UCYN.A2",
  "UCYN.B",
  "UCYN.C"
)

# Subset for only cyanos
dd_cyano <- df[which( df$taxon %in% vec_cyanos ), ]

# Plot heatmap
gg_heatmap <- ggplot(dd_cyano, aes(x = id, y = taxon, fill = tss.full)) +
  geom_tile(color = "white") +
  viridis::scale_fill_viridis(na.value = "black") +  # Use the viridis color scale for fill
  # scale_fill_gradient(low = "blue", high = "yellow", na.value = "black") +  # Define color scale
  scale_x_discrete(labels = c("SDM1", "SDM2", "SDM3", "SDM4", "SDM5", "SDM6", "SDM7", "SDM8", "SDM9", "SDM10", "SDM11", 
                              "SDM12", "SDM13", "SDM14", "SDM15", "SDM16", "SDM17", "SDM18", "SDM19", "SDM20", "SDM21", 
                              "SDM22", "SDM23", "SDM24", "SDM25", "SDM26", "SDM27", "SDM28", "SDM29", "SDM30", "SDM31", 
                              "SDM32", "SDM33", "SDM34", "SDM35", "SDM36", "SDM37", "SDM38", "SDM39", "SDM40", "SDM41", 
                              "SDM42", "SDM43", "SDM44", "SDM45", "SDM46", "SDM47", "SDM48", "SDM49", "SDM50", "SDM51", 
                              "SDM52", "SDM53", "SDM54", "SDM55", "SDM56", "SDM57", "SDM58", "SDM59", "SDM60", "SDM61", 
                              "SDM62", "SDM63", "SDM64", "SDM65", "SDM66", "SDM67", "SDM68", "SDM69", "SDM70", "SDM71", 
                              "SDM72", "SDM73", "SDM74", "SDM75", "SDM76", "SDM77", "SDM78", "SDM79", "SDM80", "SDM81", 
                              "SDM82", "SDM83", "SDM84", "SDM85", "SDM86", "SDM87", "SDM88", "SDM89", "SDM90")) +
    theme(
    legend.key.size = unit(2, 'cm'),
    panel.background = element_rect(fill = "white"),
    axis.text.y = element_text(size = 36),   
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 18),  
    axis.title.y = element_text(size = 70),
    axis.title.x = element_text(size = 70),
    legend.title = element_text(size = 36, angle = 90, vjust = 0.5),  # Rotate the legend title
    legend.text = element_text(size = 36) 
    ) +
    guides(fill=guide_legend(title="TSS Score")) +
    ylab( "" ) +
    xlab( "Ensemble members" ) +
    labs(fill = "TSS-Score") +   # Change the legend name here
    guides(fill = guide_colorbar(title.position = "right", title.hjust = 0.5))  # Use guide_colorbar for continuous legend


## Save
ggsave( paste0(wd_out, "Heatmap_TSS_scores_all_models_cyanoDiaz.png"), 
    plot = gg_heatmap, 
    width = 30,
    height = 10,
    dpi = 300,
    bg = "white"
    )


## Non-cyanobacterial community ------------------------
vec_nonCyanos <- c(
  "Gamma.A",
  "HBD-01",
  "HBD-02",
  "HBD-03",
  "HBD-04",
  "HBD-05",
  "HBD-06",
  "HBD-07",
  "HBD-09"
)
# Subset for only cyanos
dd_nonCyano <- df[which( df$taxon %in% vec_nonCyanos ), ]

# Plot heatmap
gg_heatmap <- ggplot(dd_nonCyano, aes(x = id, y = taxon, fill = tss.full)) +
  geom_tile(color = "white") +
  viridis::scale_fill_viridis(na.value = "black") +  # Use the viridis color scale for fill
  # scale_fill_gradient(low = "blue", high = "yellow", na.value = "black") +  # Define color scale
  scale_x_discrete(labels = c("SDM1", "SDM2", "SDM3", "SDM4", "SDM5", "SDM6", "SDM7", "SDM8", "SDM9", "SDM10", "SDM11", 
                              "SDM12", "SDM13", "SDM14", "SDM15", "SDM16", "SDM17", "SDM18", "SDM19", "SDM20", "SDM21", 
                              "SDM22", "SDM23", "SDM24", "SDM25", "SDM26", "SDM27", "SDM28", "SDM29", "SDM30", "SDM31", 
                              "SDM32", "SDM33", "SDM34", "SDM35", "SDM36", "SDM37", "SDM38", "SDM39", "SDM40", "SDM41", 
                              "SDM42", "SDM43", "SDM44", "SDM45", "SDM46", "SDM47", "SDM48", "SDM49", "SDM50", "SDM51", 
                              "SDM52", "SDM53", "SDM54", "SDM55", "SDM56", "SDM57", "SDM58", "SDM59", "SDM60", "SDM61", 
                              "SDM62", "SDM63", "SDM64", "SDM65", "SDM66", "SDM67", "SDM68", "SDM69", "SDM70", "SDM71", 
                              "SDM72", "SDM73", "SDM74", "SDM75", "SDM76", "SDM77", "SDM78", "SDM79", "SDM80", "SDM81", 
                              "SDM82", "SDM83", "SDM84", "SDM85", "SDM86", "SDM87", "SDM88", "SDM89", "SDM90")) +
    theme(
    legend.key.size = unit(2, 'cm'),
    panel.background = element_rect(fill = "white"),
    axis.text.y = element_text(size = 36),   
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 18),  
    axis.title.y = element_text(size = 70),
    axis.title.x = element_text(size = 70),
    legend.title = element_text(size = 36, angle = 90, vjust = 0.5),  # Rotate the legend title
    legend.text = element_text(size = 36) 
    ) +
    guides(fill=guide_legend(title="TSS Score")) +
    ylab( "" ) +
    xlab( "Ensemble members" ) +
    labs(fill = "TSS-Score") +   # Change the legend name here
    guides(fill = guide_colorbar(title.position = "right", title.hjust = 0.5))  # Use guide_colorbar for continuous legend


## Save
ggsave( paste0(wd_out, "Heatmap_TSS_scores_all_models_nonCyanoDiaz.png"), 
    plot = gg_heatmap, 
    width = 30,
    height = 10,
    dpi = 300,
    bg = "white"
    )
