### MERGE SPECIES DISTRIBUTION MODEL (SDM) PROJECTIONS INTO MONTHLY SLOTS

# This code works on locally stored files ("../9_SDM_fit/Output_9_proj_NOTION") and products are stored locally ("10_SDM_prep/Output_11_merged_NOTION")
# No parallel computing involved

# Input files: Lists of projection results and projection evaluation (TSS scores).

# Output files: List containing twelve elements (for each month), and each month element is another list
# 	containing raster objects for each modelled taxon (each taxon is modelled 5 times with varying predictor sets).

# Strategy: Model raw projection outputs consist of a stack of "number of taxa" * 5 ensemble member-models as elements
# Each of the  elements or models consists of 24 (2x12 months, 12 PRB, 12 PA) projections plotted on globally gridded (1°) and monthly climatological fields
# This code re-arranges the raw projection outputs by: 1. months and 2. taxa
# At this point no cross-validation based selection of projections applies, stacks are merely re-arranged

# Author: 	Dominic Eriksson
#			Environmental Physics Group, UP
# 			ETH Zurich
#			Switzerland
#
# 5th of September, 2022, deriksson@ethz.ch
#
# This script is based on the work of Damiano Righetti.


### =========================================================================
### Initialize system
### =========================================================================
rm(list = ls())
lib_vec <- c("raster", "doParallel")
sapply(lib_vec, library, character.only = T)

## Directories
# Input and output; NOTE: Implement also the other algorithms evaluations rf, glm and ensemble
input.eval <-  "/home/deriksson/Projects/Notion_DE/Code/6_SDM_fit/3_Output"
input.dir <- "/home/deriksson/Projects/Notion_DE/Code/6_SDM_fit/2_Output"
output.dir <- "/home/deriksson/Projects/Notion_DE/Code/7_SDM_prep/1_Output"


# Define number of models (variable combinations) fitted per taxon in input (i.e. lists with n taxa * n model runs à 24 projections)
nr <- 5

# Vectors of folder - Note: We re-run the analysis based on separation of the datasets to microscopy and sequence based datasets. Meaning dataset only 
# contains observation from either microscopy or sequence related methodologies. 
vec_folders <- c("Total_dataset", "MicroscopyBased_dataset", "SequenceBased_dataset")
f <- 1

# Vector of algorithms used (low, intermediate, high complexity)
mod.type <- c("gam", "glm", "rf")

## Define data sets (strategies)
# f = 1, Note: For f = 1 we don't need to create a separate nms.sets_eval
nms.sets <- list.files(paste0(input.dir, "/", vec_folders[f]))
nms.sets <- gsub("\\gam_proj_", "", nms.sets)
nms.sets <- gsub("\\glm_proj_", "", nms.sets)
nms.sets <- gsub("rf_proj_", "", nms.sets)
nms.sets <- gsub("\\.RData", "", nms.sets)
nms.sets <- unique(nms.sets)

# f = 2 or f = 3
nms.sets_eval <- list.files(paste0(input.eval, "/", vec_folders[f]))
nms.sets_eval <- gsub("\\gam_eval_", "", nms.sets_eval)
nms.sets_eval <- gsub("\\glm_eval_", "", nms.sets_eval)
nms.sets_eval <- gsub("rf_eval_", "", nms.sets_eval)
nms.sets_eval <- gsub("\\.csv", "", nms.sets_eval)
nms.sets_eval <- unique(nms.sets_eval)


### =========================================================================
###  PROBABILITY OF PRESENCE & PRESENCE-ABSENCE: REORDER PROJECTIONS
### =========================================================================

 # Loop across algorithms
for(v in 1:3){

	# Loop across data sets (strategies)
	for ( d in seq_along(nms.sets) ){

	# Get frames with predictive skill (each row corresponds to a taxon-specific and variable-specific model fit)
	if(f == 1){
		ls.eval <- list(); for (p in 1:length(nms.sets) ){ls.eval[[length(ls.eval) + 1]] <- read.csv(  paste0(input.eval, "/", vec_folders[f], "/", mod.type[v], "_eval_", nms.sets[p], ".csv")  , header = TRUE, sep = ",")}
	}else{
		ls.eval <- list(); for (p in 1:length(nms.sets_eval) ){ls.eval[[length(ls.eval) + 1]] <- read.csv(  paste0(input.eval, "/", vec_folders[f], "/", mod.type[v], "_eval_", nms.sets_eval[p], ".csv")  , header = TRUE, sep = ",")}
	}
	# Get data: Read in list with plankton PRB and PA projections (takes some seconds)
	ls.set <- get(load( paste0(input.dir, "/", vec_folders[f], "/", mod.type[v], "_proj_", nms.sets[d], ".RData")  )) 

	# Preparatory
	ens.PRB <-list()
	ens.PA <- list()

			# Loop across months
			for (k in 1:12){
				lisi.taxa.ens.PRB <- list()
				lisi.taxa.ens.PA <- list()
					# Loop across taxa: get five corresponding layers for each (or one in case of homogenuous approch)
					for (i in 1:length(unique(ls.eval[[d]]$taxon))){
						# PRB - Probability of presence: ensemble members; j = 4:0 is deduced from taxa(i)*5 projections to get the 5 ensemble model member projections
						for( j in (nr - 1):(nr - 5) ) {lisi.taxa.ens.PRB[length(lisi.taxa.ens.PRB) + 1] <- if(is.null(ls.set[[i * nr - j]])){list(NULL)}else{ls.set[[i * nr - j]][[k]]}}
						# PA - Presence absence: ensemble members; j = 4:0 is deduced from taxa(i)*5 projections to get the first 5 model member projections; for each member we add 12 to get PA of the corresponding month
						for( j in (nr - 1):(nr - 5) ) {lisi.taxa.ens.PA[length(lisi.taxa.ens.PA) + 1] <- if(is.null(ls.set[[i * nr - j]])){list(NULL)}else{ls.set[[i * nr - j]][[k + 12]]}}
					}
				# Arrange monthly lists in list
				ens.PRB[[k]] <- lisi.taxa.ens.PRB
				ens.PA[[k]] <- lisi.taxa.ens.PA
			}

	# Progress
	print(paste0(mod.type[v], ": set", d))

	# Save framed collection as brick (to avoid trouble with reloading the files)
	fln <- paste0(output.dir,"/", vec_folders[f], "/", mod.type[v] ,"_ens_prb_", nms.sets[d], ".grd")
	writeRaster(stack(unlist(ens.PRB)), filename = fln, format = "raster", overwrite = TRUE)

	# Save merged list-files PA (takes some seconds)
	fln <- paste0(output.dir, "/", vec_folders[f], "/", mod.type[v] , "_ens_pa_", nms.sets[d], ".grd")
	writeRaster(stack(unlist(ens.PA)), filename = fln, format = "raster", overwrite = TRUE)

	# Close loop across data set (strategies)
	}

# Close loop across algorithms
}

###==============================================================
### END
###==============================================================