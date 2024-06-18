### MERGE SPECIES DISTRIBUTION MODEL (SDM) PROJECTIONS INTO MONTHLY SLOTS. Output will be one Taxon file that is a list of length 90 (for each ensemble member). 
### Each list element itself is a raster stack of length 12 (January - December, 1 - 12).

# This code works on lobally stored files ("../9_SDM_fit/Output_9_proj_NOTION") and products are stored locally ("10_SDM_prep/Output_11_merged_NOTION")
# No parallel computing involved

# Author: 	Dominic Eriksson
#			Environmental Physics Group, UP
# 			ETH Zurich
#			Switzerland
#
# 5th of September, 2022, deriksson@ethz.ch
#
# This script is based on the work of Damiano Righetti.


# Input files: 
#       1. An evaluation file, that contains the TSS score statistic for each taxon and ensemble member run.
#       2. The SDM output for each taxon and ensemble member either hsi or pa converted.

# Output files: 
#       1.  "Taxon_all_models.rds": Taxon files that are list objects of length 90 (for each ensemble member).
#           Each list element, is a raster stack of length 12 (January - December, 1-12)


### =========================================================================
### Initialize system
### =========================================================================
rm(list = ls())
lib_vec <- c("raster", "doParallel")
sapply(lib_vec, library, character.only = TRUE)

## Directories
input.eval <-  "/home/deriksson/Projects/Notion_DE/Code/6_SDM_fit/3_Output/Total_dataset/"
input.dir <- "/home/deriksson/Projects/Notion_DE/Code/6_SDM_fit/2_Output"
output.dir <- "/net/kryo/work/deriksson/Projects/Notion_DE/Code/7_SDM_prep/1_Output/"

# Define number of models (variable combinations) fitted per taxon in input (i.e. lists with n taxa * n model runs à 24 projections)
nr <- 5

## Tax table: We only include taxa that have been modelled successfully across all 90 ensemble members 
## (/home/deriksson/Projects/Notion_DE/Code/6_SDM_fit/3_Output/Statistics/Heatmap_TSS_scores_all_models.png).
# taxa <- read.csv("/home/deriksson/Projects/Notion_DE/Code/6_SDM_fit/3_Output/Total_dataset/glm_eval_pres_abs,tot_bg_overl.csv")
# taxa <- unique(taxa$taxon)
taxa <- c(
	"UCYN.C",
	"UCYN.B",
	"UCYN.A2",
	"UCYN.A1",
	"Trichodesmium thiebautii",
	"Trichodesmium erythraeum",
	"Richelia intracellularis",
	"Richelia",
	"HBD-06",
	"HBD-05",
	"HBD-04",
	"HBD-03",
	"HBD-02",
	"Gamma.A",
	"Calothrix"
)

# Projection type
proj <- c("hsi", "pa") 

# Vectors of folder
vec_folders <- c("Total_dataset", "MicroscopyBased_dataset", "SequenceBased_dataset")
f <- 1

# Vector of algorithms used (low, intermediate, high complexity)
mod.type <- c("gam", "glm", "rf")

# Get filenames
filenames <- list.files( paste0(input.dir, "/Total_dataset/") )
filenames <- grep( "RData" , filenames, value = TRUE)
#
filenames_eval <- list.files( paste0(input.eval) )

## Loop through
for(p in seq_along(proj)){

	for(t in seq_along(taxa)){

		# Create new list object
		list_final <- list()

		for(f in seq_along(filenames)){

			# Create raster
			raster_stack <- raster()

			l <- get(load( paste0(input.dir, "/Total_dataset/", filenames[f]) ))
			eval <- gsub( "RData", "csv", filenames[f] )
			eval <- gsub( "proj", "eval", eval )
			eval <- read.csv( paste0( input.eval, eval ) )

			if( taxa[t] %in% eval[which(!is.na(eval$tss.full)), ]$taxon ){

				# Name the list elements
				names(l) <- eval$taxon

				# Only keep tss score above 0.3
				l <- l[which(!is.na(eval$tss.full))]

				# Number of predictor sets modelled
				number_predictorSets <- seq_along(which(names(l) == taxa[t]))
				l <- l[which(names(l) == taxa[t])]

				if( proj[p] == "hsi" ){

					# List names
					list_names <- paste0( filenames[f], "_", number_predictorSets )
					# Subset list
					list_to_append <- lapply(l[number_predictorSets], function(x) x[[1:12]])
					names(list_to_append) <- list_names
					# Append objects into list
					list_final <- append( list_final, list_to_append )

				}

				if( proj[p] == "pa" ){
					
					# List names
					list_names <- paste0( filenames[f], "_", number_predictorSets )
					# Subset list
					list_to_append <- lapply(l[number_predictorSets], function(x) x[[13:24]])
					names(list_to_append) <- list_names
					# Append objects into list
					list_final <- append( list_final, list_to_append )
				}
			
			} # end of if statement

		} # end of loop across filenames

		# Save list object with rasters
		saveRDS( list_final, file = paste0(output.dir, proj[p], "/", taxa[t], "_all_models.rds" ) )

	} # enf of loop across taxa

} # end of loop across projection type, hsi or pa




# ## Sanity checks
# l <- list_final
# names(l)


# levelplot(l[[1]], l[[2]], l[[3]], l[[4]], l[[5]],
#           layout = c(1, 5), col.regions = terrain.colors(100),
#           main = c("Raster 1", "Raster 2", "Raster 3", "Raster 4", "Raster 5"))

# plot(l[[28]])

# test <- calc(stack(l), mean, na.rm = TRUE)



# ### Calculate final ensemble 
# # Load the raster package
# library(raster)

# # Assuming 'raster_list' is your list of raster objects

# # Function to calculate mean, standard deviation, and coefficient of variation
# calculate_stats <- function(l) {
  
#   # Extract the number of months from the first raster in the list
#   num_months <- nlayers(l[[1]])
  
#   # Initialize lists to store results for each month
#   mean_list <- vector("list", length = num_months)
#   sd_list <- vector("list", length = num_months)
#   cv_list <- vector("list", length = num_months)
  
#   # Iterate over each month
#   for (month in 1:num_months) {
    
#     # Extract raster layers corresponding to the current month across all 90 ensemble members
#     month_layers <- lapply(l, function(r) r[[month]])
    
#     # Calculate mean raster for the current month
#     mean_raster <- calc(stack(month_layers), mean, na.rm = TRUE)
#     mean_list[[month]] <- mean_raster
    
#     # Calculate standard deviation raster for the current month
#     sd_raster <- calc(stack(month_layers), sd, na.rm = TRUE)
#     sd_list[[month]] <- sd_raster
#   }
  
#   return(list(mean = mean_list, sd = sd_list))
# #   return(list(mean = mean_list, sd = sd_list, cv = cv_list))
# }

# # Call the function to calculate statistics for each month
# stats_by_month <- calculate_stats(l)
# test <- do.call("stack", stats_by_month$sd)
