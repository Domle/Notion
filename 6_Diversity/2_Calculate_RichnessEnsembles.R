##  This script calculates ensemble means and sd of global marine diazotroph richness. Ensemble are created from diversity estimates 
##  from each month and ensemble member. We only calculate the ensembles on normalized richness.

# Input files:
#       1.  Files containing either annual or monthly richness estimates from each model.

# Output files:
#       1.  Files containing the ensemble richness either on an annual basis or monthly basis. We further calculated 
#           richness for the two groups cyanobacterial and non-cyanobacterial separately. 



# Strategy: We load the data and stack all relevant models into one raster stack and use the calc function from raster package
#           to calculate the mean across all stacked raster layers.

#           The sections include:
#               1.  Annual ensemble richness (across all 18 models) normalized for all dataset types, projection types and also the majority vote.
#               2.  Monthly ensemble richness (across all 18 models) normalized based on presence absence converted data, only for total dataset type and 
#                   also not done for majority vote option.
#               3.  Annual normalized richness for cyanobacterial and non-cyanobacterial diazotrophs only for presence absence converted projections and not done for majority
#                   vote option, only for total dataset option.

## Author: 	Dominic Eriksson
##			Environmental Physic Group
##			ETH, Zurich 
##			Switzerland 

# 5th of October 2022, deriksson@ethz.ch ------------------------------------------------------

### =========================================================================
### Preparations
### =========================================================================
# Clean system
rm(list = ls())

# Libraries
library(raster)

# Directories
wd_in <- "/net/kryo/work/deriksson/Projects/Notion_DE/Code/8_Diversity/1_Output/"
wd_out <- "/net/kryo/work/deriksson/Projects/Notion_DE/Code/8_Diversity/3_Output/"

# Load functions
source("/home/deriksson/Projects/Notion_DE/Code/Functions/Righetti/common_functions.R") # Load functions to derive latitudinal zonal mean pattern

## Load data
# Load ocean mask
msk <- subset(get(load("/home/deriksson/Projects/Notion_DE/Data/Predictor_variables/Environmental_data/VarSet07.2016.RData"))[[1]], "Mask") ## Define mask (continents)
na.layer <- raster()
na.layer[] <- NA

# Vecotrs
proj.type <- c("pa", "prb")
vec_folders <- c("Total_dataset", "MicroscopyBased_dataset", "SequenceBased_dataset")
mod.type <- c("glm", "gam", "rf")

# Define data sets (strategies)
nms.sets <- c(
"pres_abs,tot_bg_overl",
"pres_abs,gr_bg_overl",
"pres_abs,cr_bg_overl",
"pres_abs,gr_bg_nonov",
"pres_abs,cr_bg_nonov",
"pres_abs,tot_bg_nonov"
)

# Loop vector
prb_pa <- c("PA_majority_vote", "Ensemble_mean")

### =========================================================================
### 1. CALCULATE RICHNESS ANNUAL Ensemble (over all models 3x6) MEAN OF annual species richness
### =========================================================================
# Loop 
for(p in seq_along(prb_pa)){

	# Print progress
	print( paste0('Working on ', prb_pa[p]) )

    # Loop across dataset types
	for(fo in seq_along(vec_folders)){

        # Loop across projection types
		for(u in seq_along(proj.type)){

			# Directory
			if(p == 1){
				wd_div_global <- paste0(wd_in, prb_pa[p], "/", vec_folders[fo], "/div_annually/", proj.type[1] )
			}else{
				wd_div_global <- paste0(wd_in, prb_pa[p], "/", vec_folders[fo], "/div_annually/", proj.type[u] )}
			wd <- wd_div_global

			# file name vector
			vec_fn <- list.files(wd)
			vec_fn <- grep("grd", vec_fn, value = TRUE)
			if(p == 2){ 
				ind <- "_normalized" 
				vec_fn <- grep(ind, vec_fn, value = TRUE)}

			# Load data
			l <- list()
			for(f in 1:length(vec_fn)){
				# Print progress
				print(f)
				x <- raster( paste0(wd, "/", vec_fn[f]) )
				# Save data
				l[[f]] <- x
			}

			# Stack to one raster
			r <- do.call(stack, l)

			# Calculate ensemble mean and sd 
			r_ensemble_mean <- calc(r, fun = mean, na.rm = TRUE)
			names(r_ensemble_mean) <- "Total_ensemble_mean"
			r_ensemble_sd <- calc(r, fun = sd, na.rm = TRUE)
			names(r_ensemble_sd) <- "Total_ensemble_sd"
			r_ensemble_max <- calc(r, fun = max, na.rm = TRUE)
			names(r_ensemble_max) <- "Total_ensemble_max"
			r_ensemble_min <- calc(r, fun = min, na.rm = TRUE)
			names(r_ensemble_min) <- "Total_ensemble_min"
			# Calculate uncertaintiy to visualize which models keep agree with overall pattern --> divide standard deviation by model range
			ensemble_range <- r_ensemble_max - r_ensemble_min
			r_esnemble_sdRelativeToModelRange <- r_ensemble_sd/ensemble_range
			names(r_esnemble_sdRelativeToModelRange) <- "Standard_deviation_divided_by_modelRange"

			# Lon/lat gradients
			dflat <- by_lat(r_ensemble_mean) # Enable if needed
			# 
			dflon <- by_lon(r_ensemble_mean)

			# Merge and save
			r_final <- raster::stack(r_ensemble_mean, r_ensemble_sd, r_ensemble_max, r_ensemble_min, r_esnemble_sdRelativeToModelRange)
			#
			if(p == 1){
				fln <- paste0(wd_out, prb_pa[p], "/", vec_folders[fo], "/div_annually/", proj.type[1], "/Ensemble_richness_acrossAllModels.grd" )
			}else{
				fln <- paste0(wd_out, prb_pa[p], "/", vec_folders[fo], "/div_annually/", proj.type[u], "/Ensemble_richness_acrossAllModels.grd" )}
			writeRaster(r_final, filename = fln, format = "raster", overwrite = TRUE)
			#
			if(p == 1){
				fln <- paste0(wd_out, prb_pa[p], "/", vec_folders[fo], "/lon_annually/", proj.type[1], "/Longitudinal_Ensemble_richness_acrossAllModels.RData" )
			}else{
				fln <- paste0(wd_out, prb_pa[p], "/", vec_folders[fo], "/lon_annually/", proj.type[u], "/Longitudinal_Ensemble_richness_acrossAllModels.RData" )}
			save(dflon, file = fln)
			#
			if(p == 1){
				fln <- paste0(wd_out, prb_pa[p], "/", vec_folders[fo], "/lat_annually/", proj.type[1], "/Latitudinal_Ensemble_richness_acrossAllModels.RData" )
			}else{
				fln <- paste0(wd_out, prb_pa[p], "/", vec_folders[fo], "/lat_annually/", proj.type[u], "/Latitudinal_Ensemble_richness_acrossAllModels.RData" )}
			save(dflat, file = fln)
		} # Close loop across proj.type 
	} # Close loop acroos vec_folders
} # Close loop prb and pa majority vote




### ============================================================================
### 2. Richness ensemble - Monthly - Normalized
### ============================================================================
# Settings
p <- 2 # prb_pa
fo <- 1 # vec_folders
u <- 1 # proj.type


# Directory
wd_div_global <- paste0(wd_in, prb_pa[p], "/", vec_folders[fo], "/div_monthly/", proj.type[u], '/' )
wd <- wd_div_global


# Filenames
filenames <- list.files(wd)
filenames <- grep("normalized.grd", filenames, value = TRUE)

## Loop: 	Loop across month and then loop across filenames, opening each file, extracting the same month. Stack the same months
##			and use calc to get stats for each month. Stack all stats of each month as separate layer and save each month as grid file
##			separately.

# Functions: Calculate statistics
mean_na_rm <- function(x) {mean(x, na.rm = TRUE)}
#
sd_na_rm <- function(x) {sd(x, na.rm = TRUE)}
#
max_na_rm <- function(x) {max(x, na.rm = TRUE)}
#
min_na_rm <- function(x) {min(x, na.rm = TRUE)}


# Vectors for months
vec_month <- paste0('layer.', 1:12, sep = '')

# Open data
l_month <- list()
for(m in 1:12){

	# Print progress
	print(paste0('Month number: ', m))
	# Raster stack to save rasters
	r_month <- stack()
	for(f in seq_along(filenames)){

		# Print progress
		print(f)

		# Load data
		r <- stack(paste0(wd, filenames[f]))[[m]]
		# Stack data
		r_month <- stack(r_month, r)
	} # end of loop across filenames

	# Calculate stats
	mean_r <- calc(r_month, mean_na_rm)
	sd_r <- calc(r_month, sd_na_rm)
	max_r <- calc(r_month, max_na_rm)
	min_r <- calc(r_month, min_na_rm)

	# Stack stats of specific month in one raster
	raster_stack <- stack(mean_r, sd_r, max_r, min_r)
	names(raster_stack) <- c('ensemble_mean_richness', 'standard_deviation', 'maximum', 'minimum')

	# Save in list
	l_month[[m]] <- raster_stack
	rm(raster_stack)
} # end of loop across months


# Save monthly richness as list object
fln <- paste0(wd_out, prb_pa[p], "/", vec_folders[fo], "/div_monthly/", proj.type[u], "/EnsembleRichness_Monthly_AllModsels.rds")
saveRDS(l_month, file = fln)




### =========================================================================
### 3. CD's vs. NCD's
### CALCULATE RICHNESS ANNUAL Ensemble (over all models 3x6) MEAN OF annual species richness
### =========================================================================
# Settings
p <- 2 # prb_pa
fo <- 1 # vec_folders
u <- 1 # proj.type

# Directory
wd_div_global <- paste0(wd_in, prb_pa[p], "/", vec_folders[fo], "/div_annually/", proj.type[u] )
wd <- wd_div_global

# file name vector
vec_fn <- list.files(wd)
vec_fn <- grep("Cyanos", vec_fn, value = TRUE)
vec_fn <- grep('grd', vec_fn, value = TRUE)

# Load data
l_cd <- list()
l_ncd <- list()
for(f in 1:length(vec_fn)){
				
    # Print progress
	print(f)
	x_cd <- brick( paste0(wd, "/", vec_fn[f]) )$richness_cyanobacterial_diazotrophs
	x_ncd <- brick( paste0(wd, "/", vec_fn[f]) )$richness_nonCyanobacterial_diazotrophs
	# Save data
	l_cd[[f]] <- x_cd
	l_ncd[[f]] <- x_ncd
}

# Stack to one raster
r_cd <- do.call(stack, l_cd)
r_ncd <- do.call(stack, l_ncd)

# Calculate ensemble mean and sd 
r_ensemble_mean <- calc(r_cd, fun = mean, na.rm = TRUE)
names(r_ensemble_mean) <- "Total_ensemble_mean_richness_cyanos"
r_to_add <- calc(r_ncd, fun = mean, na.rm = TRUE)
r_ensemble_mean <- stack(r_ensemble_mean, r_to_add)
names(r_ensemble_mean)[2] <- "Total_ensemble_mean_richness_ncds"


# Merge and save
fln <- paste0(wd_out, prb_pa[p], "/", vec_folders[fo], "/div_annually/", proj.type[u], "/Ensemble_richness_acrossAllModels_CyanosVsNonCyanos.grd" )
writeRaster(r_ensemble_mean, filename = fln, format = "raster", overwrite = TRUE)


###==============================================================
### END
###==============================================================