### In this script we compute the normalized annual diazotroph richness for each ensemble member. Therefore the script has
### three sections, 1 = total diazotroph community, 2 = cyanobacterial and 3 = non-cyanobacterial diazotrophic community.Note to calculate
###	the annual richness from monthly SDM outputs, we intepret an annual presences, at being present at least once within twelve month.
### Actual richness is then the sum of taxa within one grid cell, normalized by the number of taxa modeled in total.
### Note we calculate the normalized richness, there fore we use proj = 1, as in the manuscript we show the normalized richness
###	from the presence absence converted data (pattern is conserved when using the hsi SDM output).

### NOTE: 	When running the code, please always double check input and output directories, to ensure you know where exactly files get saved so they
###			can be referred as the correct input files within each consecutive script.

# Author: 	Dominic Eriksson
#			Environmental Physics Group, UP
# 			ETH Zurich
#			Switzerland


# Input files:
#       1.  The taxon files that contain all 90 ensemble member SDM outputs, where each list element
#           is a raster stack from 1-12 (January - December)

# Output files:
#       1. 	Output files contain one diversity index of interest, for a specific group (total, cyanobacterial or non-cyanobacterial diazotrophs).
#			Files are saved as raster files (.grd) with 92 layers, that contain the ensemble mean, ensemble sd and all 90 ensemble members.

# 18th of June, 2024, deriksson@ethz.ch ----------------------------------------------


# Load functions
source("/home/deriksson/Projects/Notion_DE/Code/Functions/Righetti/common_functions.R") # Load functions to derive latitudinal zonal mean pattern
source("/home/deriksson/Projects/Notion_DE/Code/Functions/Dominic/call_libraries.R")


# Libraries
lib_vec <- c("fields", "raster", "maps", "sfsmisc", "doParallel", "purrr", "tidyverse", "betapart", "reshape2", "data.table")
call_libraries(lib_vec)

# # Directories
# wd_in <- "/home/deriksson/Projects/Notion_DE/Code/7_SDM_prep/1_Output_v2/"
# wd_out <- "/home/deriksson/Projects/Notion_DE/Code/Diversity/Diversities/"

# Projection type
proj.type <- c("hsi", "pa")


# Taxa included in richness
vector_taxa <- c(
    "Calothrix",
    "Cyanothece",
    "Gamma.A",
    "HBD-02",
    "HBD-03",
    "HBD-04",
    "HBD-05",
    "HBD-06",
    "Richelia intracellularis",
    "Richelia",
    "Trichodesmium erythraeum",
    "Trichodesmium thiebautii",
    "UCYN.A1",
    "UCYN.A2",
    "UCYN.B"
)

# Vector of filenames
p <- 2 # we use the presence absence converted SDM outputs
filenames <- list.files( paste0(wd_in, proj.type[p]) )

# Load data
# data <- readRDS( paste0(wd_in, proj.type[p], "/", filenames[f]) )


### ========================================================
### Total diazotrophs - Calculate annual richness normalized
### ========================================================
# Directories
wd_in <- "/home/deriksson/Projects/Notion_DE/Code/7_SDM_prep/1_Output/"
wd_out <- "/home/deriksson/Projects/Notion_DE/Code/6_Diversity/1_Output/"

## Total community --------------------------------------------------------
# Taxa included in richness
vector_taxa_all <- c(
    "Calothrix",
    "Cyanothece",
    "Gamma.A",
    "HBD-02",
    "HBD-03",
    "HBD-04",
    "HBD-05",
    "HBD-06",
    "Richelia intracellularis",
    "Richelia",
    "Trichodesmium erythraeum",
    "Trichodesmium thiebautii",
    "UCYN.A1",
    "UCYN.A2",
    "UCYN.B"
)

# Set presence absence based
p <- 2 # we use the presence absence converted outputs

# Filenames
filenames <- list.files( paste0(wd_in, proj.type[p]) )

# List that is used to calculate average
l_test_mean <- list()
# total diazotroph community
for(f in seq_along(filenames)){

	l <- readRDS( paste0(wd_in, proj.type[p], "/", filenames[f] ))

	# Create an empty list to store averaged rasters
	averaged_raster_list_mean <- vector("list", length(l))

	# Loop through each raster in the list
	for (i in seq_along(l)) {
		# Calculate average of layers
		averaged_raster_mean <- calc(l[[i]], mean, na.rm = TRUE)
		
		# Store averaged raster in the new list
		averaged_raster_list_mean[[i]] <- averaged_raster_mean
	}
	if(f == 1){
		l_test_mean <- averaged_raster_list_mean
	}else{
		l_test_mean <- Map(function(x, y) x + y, l_test_mean, averaged_raster_list_mean)
	}
}

# Normalize by number of species, n = 15
list_mean <- lapply(l_test_mean, function(r) {
  calc(r, function(x) x / length(filenames))
})

# Stack rasters
raster_stack <- do.call("stack", list_mean)
raster_stack_mean <- calc(raster_stack, mean, na.rm = TRUE)
names(raster_stack_mean) <- "mean"
raster_stack_sd <- calc(raster_stack, sd, na.rm = TRUE)
names(raster_stack_sd) <- "sd"
raster_stack <- stack( raster_stack_mean, raster_stack_sd, raster_stack )

# Save 
fln <- paste0( wd_out, proj.type[p], "/RichnessNormalized_annual_totalDiaz_Ensemble_and_EnsembleMembers.grd" )
writeRaster(
	raster_stack,
	fln,
	overwrite = TRUE
)

### =================================================================
### Cyanobacterial diazotrophs - Calculate annual richness normalized
### =================================================================
# Filenames
filenames_cyanos <- c(
	"Calothrix_all_models.rds",                              
	"Richelia intracellularis_all_models.rds",
	"Richelia_all_models.rds",                
	"Trichodesmium erythraeum_all_models.rds",
	"Trichodesmium thiebautii_all_models.rds",
	"UCYN.A1_all_models.rds",                 
	"UCYN.A2_all_models.rds",                 
	"UCYN.B_all_models.rds",                  
	"UCYN.C_all_models.rds")
# List that is used to calculate average
l_test_mean <- list()
for(f in seq_along(filenames_cyanos)){

	l <- readRDS( paste0(wd_in, proj.type[p], "/", filenames_cyanos[f] ))

	# Create an empty list to store averaged rasters
	averaged_raster_list_mean <- vector("list", length(l))

	# Loop through each raster in the list
	for (i in seq_along(l)) {
		# Calculate average of layers
		averaged_raster_mean <- calc(l[[i]], mean, na.rm = TRUE)
		
		# Store averaged raster in the new list
		averaged_raster_list_mean[[i]] <- averaged_raster_mean
	}
	if(f == 1){
		l_test_mean <- averaged_raster_list_mean
	}else{
		l_test_mean <- Map(function(x, y) x + y, l_test_mean, averaged_raster_list_mean)
	}
}

# Normalize by number of species
list_mean <- lapply(l_test_mean, function(r) {
  calc(r, function(x) x / length(filenames_cyanos))
})

# Stack rasters
raster_stack <- do.call("stack", list_mean)
raster_stack_mean <- calc(raster_stack, mean, na.rm = TRUE)
names(raster_stack_mean) <- "mean"
raster_stack_sd <- calc(raster_stack, sd, na.rm = TRUE)
names(raster_stack_sd) <- "sd"
raster_stack <- stack( raster_stack_mean, raster_stack_sd, raster_stack )

# Save 
fln <- paste0( wd_out, proj.type[p], "/RichnessNormalized_annual_cyanoDiaz_Ensemble_and_EnsembleMembers.grd" )
writeRaster(
	raster_stack,
	fln,
	overwrite = TRUE
)

### =================================================================
### Non-cyanobacterial diazotrophs - Calculate annual richness normalized
### =================================================================
# Filenames
filenames_nonCyanos <- c(
	"Gamma.A_all_models.rds",
	"HBD-02_all_models.rds",
	"HBD-03_all_models.rds",
	"HBD-04_all_models.rds",                  
	"HBD-05_all_models.rds",                  
	"HBD-06_all_models.rds"                  
	)
# List that is used to calculate average
l_test_mean <- list()
for(f in seq_along(filenames_nonCyanos)){

	l <- readRDS( paste0(wd_in, proj.type[p], "/", filenames_nonCyanos[f] ))

	# Create an empty list to store averaged rasters
	averaged_raster_list_mean <- vector("list", length(l))

	# Loop through each raster in the list
	for (i in seq_along(l)) {
		# Calculate average of layers
		averaged_raster_mean <- calc(l[[i]], mean, na.rm = TRUE)
		
		# Store averaged raster in the new list
		averaged_raster_list_mean[[i]] <- averaged_raster_mean
	}
	if(f == 1){
		l_test_mean <- averaged_raster_list_mean
	}else{
		l_test_mean <- Map(function(x, y) x + y, l_test_mean, averaged_raster_list_mean)
	}
}

# Compute average
list_mean <- lapply(l_test_mean, function(r) {
  calc(r, function(x) x / length(filenames_nonCyanos))
})

# Stack rasters
raster_stack <- do.call("stack", list_mean)
raster_stack_mean <- calc(raster_stack, mean, na.rm = TRUE)
names(raster_stack_mean) <- "mean"
raster_stack_sd <- calc(raster_stack, sd, na.rm = TRUE)
names(raster_stack_sd) <- "sd"
raster_stack <- stack( raster_stack_mean, raster_stack_sd, raster_stack )

# Save 
fln <- paste0( wd_out, proj.type[p], "/RichnessNormalized_annual_nonCyanoDiaz_Ensemble_and_EnsembleMembers.grd" )
writeRaster(
	raster_stack,
	fln,
	overwrite = TRUE
)

