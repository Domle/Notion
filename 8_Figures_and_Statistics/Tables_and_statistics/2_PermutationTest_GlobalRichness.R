### In this script we perform a permutation test to see how the global richness patterns are influenced, when removing one species at a time.

### NOTE: 	When running the code, please always double check input and output directories, to ensure you know where exactly files get saved so they
###			    can be referred as the correct input files within each consecutive script.

# Author: 	
#     Dominic Eriksson
#			Environmental Physics Group, UP
# 		ETH Zurich
#			Switzerland


# Input files:
#       1.  

# Output files:
#       1. 

# 18th of June, 2024, deriksson@ethz.ch ----------------------------------------------



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


# Load functions
source("/home/deriksson/Projects/Notion_DE/Code/Functions/Dominic/subset_ocean_regions.R")


# Load Table
r <- brick("/home/deriksson/Projects/Notion_DE/Code/Tables/1_Output/FinalTable_AllMetrics.grd")
names(r)

df <- as.data.frame(r, xy = TRUE)

# Vectors
vec_taxa <- grep( "Species", names(df), value = TRUE )
vec_taxa <- grep( "hsi" , vec_taxa, value = TRUE)
vec_taxa <- grep( "Mean" , vec_taxa, value = TRUE)

# Subset dataframe
df <- df[, c( vec_taxa )]
# Remove NAs
df <- na.omit(df)

# Load required library
library(cluster)

# Function to calculate Silhouette index for given data
calculate_silhouette <- function(data) {
  dist_matrix <- dist(data)
  silhouette_index <- silhouette(cutree(hclust(dist_matrix), k = 2), dist_matrix)
  mean(silhouette_index[, 3])  # Mean Silhouette width
}

# Initialize variables
set.seed(123)  # For reproducibility
total_species <- 15
removal_numbers <- 1:14  # Adjust this range as needed
iterations <- 100  # Number of permutations

# Create a mock dataset (for example purposes)
species_data <- as.matrix(df)


# Initialize a list to store results
results <- list()

# Perform permutation test
for (num_remove in removal_numbers) {

  # Print progress
  print( paste0( "Remove ", num_remove, " taxa" ) )
  #
  silhouette_scores <- numeric(iterations)
  
  for (i in 1:iterations) {

    # Print progress
    print( paste0( "Iteration ", i, " out of ",  iterations) )

    remaining_indices <- sample(1:total_species, total_species - num_remove)
    remaining_data <- species_data[, remaining_indices]
    
    silhouette_scores[i] <- calculate_silhouette(remaining_data)
  }
  
  results[[paste0("Remove_", num_remove, "_species")]] <- mean(silhouette_scores)
}


# Convert results to a data frame for plotting
results_df <- data.frame(
  Num_Removed_Species = removal_numbers,
  Mean_Silhouette_Index = unlist(results)
)

# Plot the results
gg_silhoutte <- ggplot(results_df, aes(x = Num_Removed_Species, y = Mean_Silhouette_Index)) +
  geom_line() +
  geom_point() +
  labs(
    x = "Number of Removed Species",
    y = "Mean Silhouette Index",
    title = "Effect of Removing Species on Silhouette Index"
  )


# Save 
fln <- "/home/deriksson/Projects/Notion_DE/Code/Diversity/Visualizations/test.png"
ggsave(
  plot = gg_silhoutte,
  filename = fln
)