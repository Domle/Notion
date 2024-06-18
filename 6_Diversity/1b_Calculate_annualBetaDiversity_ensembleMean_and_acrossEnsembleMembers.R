### In this script we compute the betadiversity (Jaccard dissimilarity index, species turnover and nestedness) for each ensemble member. To calculate this
### we use the package betapart. For each location, we calculate betadiversity components across all possible combinations. We then compute
### a global average for each grid cell, that represents the grid cell based mean state of the betadiversity component, when compared to the global ocean (spatial effect of
### community composition).
### To analyze the temporal effect of compunity, we also compute the diversity changes at one location across the twelve month of the year and compute an average at the end. This
### helps us to understand where temporal effect might be more important.
### Note we intepret an annual presences, at being present at least once within twelve month.

### NOTE: 	When running the code, please always double check input and output directories, to ensure you know where exactly files get saved so they
###			can be referred as the correct input files within each consecutive script.

### NOTE:   Because the computation and the global averaging of the betadiversity takes up to 24 hours, we created temporal files 
###         of each ensemble member separately. At the very end we compile the raster file that will contain ensemble mean and all ensemble members.

# Author: 	Dominic Eriksson
#			Environmental Physics Group, UP
# 			ETH Zurich
#			Switzerland


# Input files:
#       1.  The taxon files that contain all 90 ensemble member SDM outputs, where each list element
#           is a raster stack from 1-12 (January - December)

# Output files:
#       1.	The output files are raster stacks with 92 layers, the first 2 are the ensemble mean and standard deviation, followed
#			by 90 ensemble member based betadiversity estimates.

# 18th of June, 2024, deriksson@ethz.ch ----------------------------------------------


### ==========================================================================================
## Total diazotrophs - Calculate Global Annual Betadiversity (Jaccard, Turnover & Nestedness)
### ==========================================================================================
# Libraries
library(betapart)

# We work on presence absence data, therefore we choose p <- 2
proj.type <- c("hsi", "pa")
p <- 2

# Directories
wd_in <- "/home/deriksson/Projects/Notion_DE/Code/7_SDM_prep/1_Output_v2/"
wd_out <- "/home/deriksson/Projects/Notion_DE/Code/Diversity/Diversities/pa/Betadiversity/"

# Vector of filenames
filenames <- list.files( paste0(wd_in, proj.type[p]))

# Filenames
filenames <- c(
	"Calothrix_all_models.rds",
 	"Gamma.A_all_models.rds",              
 	"HBD-02_all_models.rds",                
 	"HBD-03_all_models.rds",                 
 	"HBD-04_all_models.rds",                 
 	"HBD-05_all_models.rds",                 
 	"HBD-06_all_models.rds",                 
 	"Richelia intracellularis_all_models.rds",
 	"Richelia_all_models.rds",      
	"Trichodesmium erythraeum_all_models.rds",
	"Trichodesmium thiebautii_all_models.rds",
	"UCYN.A1_all_models.rds", 
	"UCYN.A2_all_models.rds",                 
	"UCYN.B_all_models.rds",                  
	"UCYN.C_all_models.rds" 
)


# Create an empty list to store averaged rasters
l <- readRDS( paste0(wd_in, proj.type[p], "/", filenames[1] ))
averaged_raster_list_mean <- vector("list", length(l))
# total diazotroph community
for(f in seq_along(filenames)){

	# Print progress
	print(f)

	# Load data
	l <- readRDS( paste0(wd_in, proj.type[p], "/", filenames[f] ))

	# Loop through each raster in the list
	for (i in seq_along(l)) {
		
		if(f == 1){
			# Calculate average of layers
			averaged_raster_mean <- calc(l[[i]], mean, na.rm = TRUE)
			
			# Store averaged raster in the new list
			averaged_raster_list_mean[[i]] <- averaged_raster_mean
		}
		# Stack biogeography from other taxa
		if( f > 1 ){
			# Calculate average of layers/months
			averaged_raster_mean <- calc(l[[i]], mean, na.rm = TRUE)
			# Stack files
			averaged_raster_list_mean[[i]] <- stack( averaged_raster_list_mean[[i]], averaged_raster_mean )
		}

	} # end of loop across ensemble members
} # end loop across filenames/taxa files


# Convert to annual presence absence. NOTE: For an annual presence, it is enough if the species appeared once in a month,
# and is missing the other 11 months.
averaged_raster_list_mean <- lapply(averaged_raster_list_mean, function(x) {
	x[x > 0] <- 1 # We interpret an annual presence, with presence at least once within 12 month
	return(x)
})

l <- list()
# We loop through the 90 ensemble members
for(i in 1:90){

	# Print progress
	print(i)
	
	# Convert to dataframe
    df <- as.data.frame( averaged_raster_list_mean[[i]], xy = TRUE )

	# Save locations with NAs, beta.pair function only works without NAs, we rowbind the NAs again later
	na <- df[!complete.cases(df[, 3:ncol(df)]), ]
	na <- na[, c(1:2)]
	na$jac <- NA
	na$jtu <- NA
	na$jne <- NA

	# Remove NAs
	df <- df[which(!is.na(df[, 3])), ]

	# Calculate beta diversity using beta.pair function
	beta <- beta.pair(df[, 3:ncol(df)], index.family = "jaccard")

	## Extract Jacchard Index, species turnover, nestedness
	jac <- as.matrix(beta$beta.jac)
	jac <- reshape2::melt(jac)
	# Convert the dataframe to a data.table for efficient processing
	dt <- as.data.table(jac)
	# Specify the columns for which you want to calculate the mean
	# Replace 'column_name' with the actual column name
	columns_to_aggregate <- c("value")
	# Aggregate the mean for specified columns
	jac <- dt[, .(jac = mean(value, na.rm = TRUE)), by = Var1]

	jtu <- as.matrix(beta$beta.jtu)
	jtu <- reshape2::melt(jtu)
	# Convert the dataframe to a data.table for efficient processing
	dt <- as.data.table(jtu)
	# Calculate the mean for specified columns
	jtu <- dt[, .(jtu = mean(value, na.rm = TRUE)), by = Var1]

	jne <- as.matrix(beta$beta.jne)
	jne <- reshape2::melt(jne)
	# Convert the dataframe to a data.table for efficient processing
	dt <- as.data.table(jne)
	# Specify the columns for which you want to calculate the mean
	# Replace 'column_name' with the actual column name
	columns_to_aggregate <- c("value")
	# Calculate the mean for specified columns
	jne <- dt[, .(jne = mean(value, na.rm = TRUE)), by = Var1]

	jne_beta <- as.matrix(beta$beta.jne)
	jne_beta <- reshape2::melt(jne_beta)
	jac_beta <- as.matrix(beta$beta.jac)
	jac_beta <- reshape2::melt(jac_beta)
	betaRatio <- jne_beta[, 1:2] # first two columns are the aggregated id variables
	betaRatio$value <- jne_beta$value/jac_beta$value
	# Convert the dataframe to a data.table for efficient processing
	dt <- as.data.table(betaRatio)
	# Calculate the mean for specified columns
	betaRatio <- dt[, .(betaRatio = mean(value, na.rm = TRUE)), by = Var1]

	# Create dataframe
	df <- data.frame(
	x = df$x,
	y = df$y,
	jac = jac[, 2],
	jtu = jtu[, 2],
	jne = jne[, 2],
	betaRatio = betaRatio[, 2]
	)

	# Save in list and convert to raster
	r <- rasterFromXYZ( df )

    # Save 
	writeRaster(
		r,
		filename = paste0("/home/deriksson/Projects/Notion_DE/Code/Diversity/Diversities/pa/BetaDiversity_annual_updatedBetaRaio/Betadiversity_annual_totalDiaz_ensembleMember", i, ".grd"),
		overwrite = TRUE
	)
} # end loop across ensemble members


## Create one file containing ensemble mean and ensemble members
# Get filenames from the temporary files
filenames <- list.files("/home/deriksson/Projects/Notion_DE/Code/Diversity/Diversities/pa/BetaDiversity_annual_updatedBetaRaio/")
filenames <- grep("grd", filenames, value = TRUE)
filenames <- grep("total", filenames, value = TRUE)
# Load data
raster_stack <- lapply(paste0("/home/deriksson/Projects/Notion_DE/Code/Diversity/Diversities/pa/BetaDiversity_annual_updatedBetaRaio/", filenames), brick)

# Merge
l <- raster_stack
l_jac <- lapply(l, function(r) r[[1]])
l_tu <- lapply(l, function(r) r[[2]])
l_ne <- lapply(l, function(r) r[[3]])
l_betaRatio <- lapply(l, function(r) r[[4]])

## Stack and calculate statistics
# Define the WGS84 CRS
crs_wgs84 <- CRS("+proj=longlat +datum=WGS84")

## Calculate ensembles (mean and standard deviation)
# Jaccard dissimilarity
stack_jac <- do.call( "stack", l_jac )
stack_mean <- calc(stack_jac, mean, na.rm = TRUE)
stack_sd <- calc(stack_jac, sd, na.rm = TRUE)
rm(stack_jac)
stack_jac <- stack(stack_mean, stack_sd)
names(stack_jac) <- c("ensemble_mean", "ensemble_sd")
# Add individual ensemble members
stack_jac <- stack( stack_jac, do.call("stack", l_jac) )
# Assign the CRS to the raster
crs(stack_jac) <- crs_wgs84

# Species turnover
stack_tu <- do.call( "stack", l_tu )
stack_mean <- calc(stack_tu, mean, na.rm = TRUE)
stack_sd <- calc(stack_tu, sd, na.rm = TRUE)
rm(stack_tu)
stack_tu <- stack(stack_mean, stack_sd)
names(stack_tu) <- c("ensemble_mean", "ensemble_sd")
# Add individual ensemble members
stack_tu <- stack( stack_tu, do.call("stack", l_tu) )
# Assign the CRS to the raster
crs(stack_tu) <- crs_wgs84

# Nestedness
stack_ne <- do.call( "stack", l_ne )
stack_mean <- calc(stack_ne, mean, na.rm = TRUE)
stack_sd <- calc(stack_ne, sd, na.rm = TRUE)
rm(stack_ne)
stack_ne <- stack(stack_mean, stack_sd)
names(stack_ne) <- c("ensemble_mean", "ensemble_sd")
# Add individual ensemble members
stack_ne <- stack( stack_ne, do.call("stack", l_ne) )
# Assign the CRS to the raster
crs(stack_ne) <- crs_wgs84

# Beta ratio
stack_betaRatio <- do.call( "stack", l_betaRatio )
stack_mean <- calc(stack_betaRatio, mean, na.rm = TRUE)
stack_sd <- calc(stack_betaRatio, sd, na.rm = TRUE)
rm(stack_betaRatio)
stack_ne <- stack(stack_mean, stack_sd)
names(stack_betaRatio) <- c("ensemble_mean", "ensemble_sd")
# Add individual ensemble members
stack_ne <- stack( stack_betaRatio, do.call("stack", l_betaRatio) )
# Assign the CRS to the raster
crs(stack_betaRatio) <- crs_wgs84


# Save raster files
writeRaster(
	stack_jac,
	file = paste0( "/home/deriksson/Projects/Notion_DE/Code/Ensembles/Diversity/", "Annual_Jaccard_totalDiaz_Ensemble_and_ensembleMembers.grd" ),
	overwrite = TRUE
)
writeRaster(
	stack_tu,
	file = paste0( "/home/deriksson/Projects/Notion_DE/Code/Ensembles/Diversity/", "Annual_Turnover_totalDiaz_Ensemble_and_ensembleMembers.grd" ),
	overwrite = TRUE
)
writeRaster(
	stack_ne,
	file = paste0( "/home/deriksson/Projects/Notion_DE/Code/Ensembles/Diversity/", "Annual_Nestedness_totalDiaz_Ensemble_and_ensembleMembers.grd" ),
	overwrite = TRUE
)
writeRaster(
	stack_betaRatio,
	file = paste0( "/home/deriksson/Projects/Notion_DE/Code/Ensembles/Diversity/", "Annual_BetaRatio_totalDiaz_Ensemble_and_ensembleMembers.grd" ),
	overwrite = TRUE
)
