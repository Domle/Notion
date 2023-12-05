### Generate ensemble across the predictor sets. We will follow two strategies. First strategy will create files that contain the mean,
### standard deviation, maximum and minimum of the habitat suitability values and the presence absence converted projections. As a second 
### strategy we will also introduce a majority vote transformation of the data. We do this to later analyze beta diversity where a binary 
### output is needed of either a species being present or absent. So we use a majority vote of assigning either a presence or absence, with 
### a presence needing at least a 50% agreement across models.

# Input files:
#   Input 1:    Merged ensemble model lists with 12 slots, list contains 12 months à n taxa * model members used per ensemble (i.e. 5)
#               There is  an input list for probability-of-presence layers and another for presence-absence layers, each.
#   Input 2:    Evaluation dataframes (note that 5 member models were fitted per species or taxon; these members are selected for further analysis, 
#               accordingly true skill statistic - TSS)
# Output files:
#   Output 1:   'Species occurrence or probability lists' or rasters (12 monthly slots), each containing the species' ensemble mean layers (raster projections) as a list in regard to the occurrence or probability of occurrence.
#               The ensembles are means of members (min or max projection of single members, across the members; this can be used to calculate the variable-based spread), or standard deviations of the members.


# Strategy: Member models with reasonable predictive skill are selected (based on TSS > 0.30) towards ensemble mean 
#           models (n members = 1 : 5) for each species.

# Author: Dominic Eriksson
#         Environmental Physics Group, UP
#         ETH, Zurich
#         Switzerland
#
# deriksson@ethz.ch, 15th of September 2022 --------------------------------------------------------------------------


### =========================================================================
### Initialize system
### =========================================================================
rm(list = ls())
require(raster); require(sfsmisc); require(doParallel)

# Directories
input.eval <- "/home/deriksson/Projects/Notion_DE/Code/6_SDM_fit/3_Output"
input.dir <- "/home/deriksson/Projects/Notion_DE/Code/7_SDM_prep/1_Output"
output.dir <- "/home/deriksson/Projects/Notion_DE/Code/7_SDM_prep/2a_Output/Ensemble_means"
output.dir_majorityVote <- "/home/deriksson/Projects/Notion_DE/Code/7_SDM_prep/2a_Output/PA_majority_vote"


# Define number of models (variable combinations) fitted per taxon in input
nr <- 5 # Number of ensemble members fitted per taxon
thres <- 0.3 # ** ADJUST IN CASE: TSS based selection threshold

# Vectors
vec_folders <- c("Total_dataset", "MicroscopyBased_dataset", "SequenceBased_dataset")

# Create NA.layer and zero.layer
na.layer <- raster()
na.layer[] <- NA
zero.layer <- na.layer
zero.layer [] <- 0

# Vector of algorithms used (low, intermediate, high complexity)
mod.type <- c("gam", "glm", "rf")

# Vector of projection types used
proj.type <- c("pa", "prb")

# Vector of months
months <- c(
	"January",
	"February",
	"March",
	"April",
	"May",
	"June",
	"July",
	"August",
	"September",
	"October",
	"November",
	"December"
)



### =========================================================================
### Prepare ensemble: later used to generate diversity estimate from stacked species layers - No backtransformation from pa data after averaging (majority vote)
### =========================================================================

# Loop acroos datasets (total, microscopy based, and sequence based)
for(f in seq_along(vec_folders)){

	# If condition to get correct file names
	if(f == 1){
		# Define data sets (strategies)
		nms.sets <- list.files(paste0(input.dir, "/", vec_folders[f]))
		nms.sets <- gsub("\\gam_ens_pa_", "", nms.sets)
		nms.sets <- gsub("\\glm_ens_pa_", "", nms.sets)
		nms.sets <- gsub("rf_ens_pa_", "", nms.sets)
		#
		nms.sets <- gsub("\\gam_ens_prb_", "", nms.sets)
		nms.sets <- gsub("\\glm_ens_prb_", "", nms.sets)
		nms.sets <- gsub("rf_ens_prb_", "", nms.sets)
		#
		nms.sets <- gsub("\\.grd", "", nms.sets)
		nms.sets <- gsub("\\.gri", "", nms.sets)
		nms.sets <- unique(nms.sets)
	}else{
		# f = 2 or f = 3
		nms.sets_eval <- list.files(paste0(input.eval, "/", vec_folders[f]))
		nms.sets_eval <- gsub("\\gam_eval_", "", nms.sets_eval)
		nms.sets_eval <- gsub("\\glm_eval_", "", nms.sets_eval)
		nms.sets_eval <- gsub("rf_eval_", "", nms.sets_eval)
		nms.sets_eval <- gsub("\\.csv", "", nms.sets_eval)
		nms.sets_eval <- unique(nms.sets_eval)
	}
	# Loop across projection types
	for (u in 1:2){
	# u = 1 # 1 = pa, 2 = prb

		# LOOP ACROSS ALGORITHMS
		for(v in 1:3){
		# v = 1

			# Loop across data sets (strategies)
			for (d in seq_along(nms.sets) ){
				# Get frames with predictive skill (each row corresponds to a taxon-specific and variable-specific model fit)
				if(f == 1){
					ls.eval <- list(); for (p in 1:length(nms.sets) ){ls.eval[[length(ls.eval) + 1]] <- read.csv(  paste0(input.eval, "/", vec_folders[f], "/", mod.type[v],"_eval_", nms.sets[p], ".csv"), header = TRUE, sep = ",")}
					}
				if(f == 2){
					ls.eval <- list(); for (p in 1:length(nms.sets_eval) ){ls.eval[[length(ls.eval) + 1]] <- read.csv(  paste0(input.eval, "/", vec_folders[f], "/", mod.type[v], "_eval_", nms.sets_eval[p], ".csv"), header = TRUE, sep = ",")}
				}
				if(f == 3){
					ls.eval <- list(); for (p in 1:length(nms.sets_eval) ){ls.eval[[length(ls.eval) + 1]] <- read.csv(  paste0(input.eval, "/", vec_folders[f], "/", mod.type[v], "_eval_", nms.sets_eval[p], ".csv"), header = TRUE, sep = ",")}
				}
				# Get data from brick (to avoid trouble with reloading the files) and allocate into list of 12 monthly slots with pre-gathered projections
				prj.stack <- brick(paste0(input.dir, "/", vec_folders[f], "/", mod.type[v], "_ens_", proj.type[u], "_", nms.sets[d], ".grd"))
				ls.dat.full <- list()   ;   nelem <- nlayers(prj.stack)/12   ;   from <- seq(1, nlayers(prj.stack), by = nelem)   ;   for(q in 1:12){   ls.dat.full[[q]] <- as.list(prj.stack[[ from[q]:(from[q] + nelem - 1)]])   }
				# Get data: read in list (12 monthly slots) with pre-merged projections (takes some 2 minutes)
				# Prepare lists to store monthly data
				lisi.months.mean <- list()
				lisi.months.max <- list()
				lisi.months.min <- list()
				lisi.months.sd <- list()

				# Loop across months
				for (k in 1:12){
					ls.dat <- ls.dat.full[[k]] # get monthly subset of data

					# PARALLEL COMPUTING ACROSS TAXA (cluster may not work due to memory constraints; try with 10 cores remotely, 4 locally for PA and 3 for PRB)
					n.cores <- 20
					cl = makeCluster(n.cores, outfile = "")
					registerDoParallel(cl)
					lisi.taxa <- foreach(i = 1 : length(unique(ls.eval[[d]]$taxon)) , .packages = c('raster', 'sfsmisc'))%dopar%{

						# Define the number of ensemble members that pass a selective skill-threshold
						len <- sum(   ls.eval[[d]][,"tss.xval4"]    [  c((i * nr - (nr - 1)):(i * nr))  ]  >= thres, na.rm = TRUE )

						if(len > 0){
							# Calculate mean projection of useful member models of taxon "i" in a monthly list, subselect members with TSS >= thres
							# From individual monthly lists, the species monthly member model projections (5 layers per species) are selected as follows:
							# i*5 - 4 (the first ensemble layer of the species) to i*5 (the last layer of the species); however from the eval sheet the corresponding positions are: i*6 - 5 (the first position) to i*6 -1 (the 5th position)
							# An additional point might be to transfom the derived ensemble means to pres if value >= 0.6, and abs if value < 0.6, yet the result is better (or at leas, smoother) when keeping the mean.
							if (len > 1){
								lisi.taxa.mean <- calc(  stack(     ls.dat    [c((i * nr - (nr - 1)):(i * nr)) [  which(ls.eval[[d]][,"tss.xval4"] [  c((i * nr - (nr - 1)):(i * nr)) ] >= thres)  ]]     ), mean, na.rm = TRUE)
								lisi.taxa.max <- calc(  stack(     ls.dat     [c((i * nr - (nr - 1)):(i * nr)) [  which(ls.eval[[d]][,"tss.xval4"] [  c((i * nr - (nr - 1)):(i * nr)) ] >= thres)  ]]     ), max, na.rm = TRUE)
								lisi.taxa.min <- calc(  stack(     ls.dat    [c((i * nr - (nr - 1)):(i * nr)) [  which(ls.eval[[d]][,"tss.xval4"] [  c((i * nr - (nr - 1)):(i * nr)) ] >= thres)  ]]      ), min, na.rm = TRUE)
								lisi.taxa.sd <- calc(  stack(   ls.dat    [c((i * nr - (nr - 1)):(i * nr)) [  which(ls.eval[[d]][,"tss.xval4"] [  c((i * nr - (nr - 1)):(i * nr)) ] >= thres)  ]]     ), sd, na.rm = TRUE)
							}else{ # if len = 1
								lisi.taxa.mean <- ls.dat [[  c((i * nr - (nr - 1)):(i * nr))  [   which(ls.eval[[d]][, "tss.xval4"]  [c((i * nr - (nr - 1)):(i * nr))]  >= thres)  ]  ]]
								lisi.taxa.min <- ls.dat [[  c((i * nr - (nr - 1)):(i * nr))  [   which(ls.eval[[d]][, "tss.xval4"]  [c((i * nr - (nr - 1)):(i * nr))]  >= thres)  ]  ]]
								lisi.taxa.max <- ls.dat [[  c((i * nr - (nr - 1)):(i * nr))  [   which(ls.eval[[d]][, "tss.xval4"]  [c((i * nr - (nr - 1)):(i * nr))]  >= thres)  ]  ]]
								lisi.taxa.sd <- zero.layer # standard deviation of single value is 0
							}
						}else{ # if len = 0
							lisi.taxa.mean <- na.layer
							lisi.taxa.min <- na.layer
							lisi.taxa.max <- na.layer
							lisi.taxa.sd <- na.layer
						}
						# Progress
						print(paste0(mod.type[v], " ", d , " tax", i, " month", k, " ", proj.type[u]))
						# Define output
						tax.merged <-list()
						tax.merged[[1]] <- lisi.taxa.mean
						tax.merged[[2]] <- lisi.taxa.max
						tax.merged[[3]] <- lisi.taxa.min
						tax.merged[[4]] <- lisi.taxa.sd
						tax.merged
						} # Close parallel computing across taxa
						stopCluster(cl)

						# Name the list elements according to the taxa
						names(lisi.taxa) <- unique(ls.eval[[d]]$taxon)

						# Store in list
						lisi.months.mean[[k]] <- lapply(lisi.taxa, '[[',1)
						lisi.months.max[[k]] <- lapply(lisi.taxa, '[[',2)
						lisi.months.min[[k]] <- lapply(lisi.taxa, '[[',3)
						lisi.months.sd[[k]] <- lapply(lisi.taxa, '[[',4)

						# Progress
						print(paste0(mod.type[v], " ", d, " ", proj.type[u], " month ", k," completed..."))
					} # Close loop across months

					# Name the months in the list
					names(lisi.months.mean) <- months
					names(lisi.months.max) <- months
					names(lisi.months.min) <- months
					names(lisi.months.sd) <- months

					# Save framed collection (as raster stack to avoid trouble with reloading the files)
					# Save ensemble mean
					fln <- paste0(output.dir, "/", vec_folders[f], "/", proj.type[u], "/", mod.type[v], "_", proj.type[u], "_mean_", nms.sets[d], ".grd")
					writeRaster(stack(unlist(lisi.months.mean)), filename = fln, format = "raster", overwrite = TRUE) 
					fln <- paste0(output.dir, "/", vec_folders[f], "/", proj.type[u], "/", mod.type[v] , "_", proj.type[u], "_mean_", nms.sets[d], ".RData") 
					save(lisi.months.mean, file = fln)
					# Save minimum and maximum
					fln <- paste0(output.dir,"/", vec_folders[f], "/", proj.type[u], "/", mod.type[v] ,"_",proj.type[u],"_min_", nms.sets[d], ".RData")
					save(lisi.months.min, file = fln)
					fln <- paste0(output.dir,"/", vec_folders[f], "/", proj.type[u], "/", mod.type[v] , "_", proj.type[u],"_max_", nms.sets[d], ".RData") 
					save(lisi.months.max, file = fln)
					# Save standard deviation
					fln <- paste0(output.dir,"/", vec_folders[f], "/", proj.type[u], "/", mod.type[v], "_", proj.type[u], "_sd_", nms.sets[d], ".grd")
					writeRaster(stack(unlist(lisi.months.sd)), filename = fln, format = "raster", overwrite = TRUE) 
					fln <- paste0(output.dir,"/", vec_folders[f], "/", proj.type[u], "/", mod.type[v] ,"_", proj.type[u], "_sd_", nms.sets[d], ".RData") 
					save(lisi.months.sd, file = fln)

				} # Close loop across data strategies
			} # Close loop across algorithms
	} # Close loop across projection type
} # Close loop across datasets (total, microscopy based, sequence based)


### =========================================================================
### Prepare ensemble: later used to generate diversity estimate from stacked species layers - With backtransformation from pa data after averaging (majority vote) 
### =========================================================================

# Directories
input.eval <- "/home/deriksson/Projects/Notion_DE/Code/6_SDM_fit/3_Output"
input.dir <- "/home/deriksson/Projects/Notion_DE/Code/7_SDM_prep/1_Output"
output.dir <- "/home/deriksson/Projects/Notion_DE/Code/7_SDM_prep/2a_Output/PA_majority_vote"


# LOOP across projection types
# Select presence absence data
u <- 1

for(f in seq_along(vec_folders)){

	if(f == 1){
		# Define data sets (strategies)
		nms.sets <- list.files(paste0(input.dir, "/", vec_folders[f]))
		nms.sets <- gsub("\\gam_ens_pa_", "", nms.sets)
		nms.sets <- gsub("\\glm_ens_pa_", "", nms.sets)
		nms.sets <- gsub("rf_ens_pa_", "", nms.sets)
		#
		nms.sets <- gsub("\\gam_ens_prb_", "", nms.sets)
		nms.sets <- gsub("\\glm_ens_prb_", "", nms.sets)
		nms.sets <- gsub("rf_ens_prb_", "", nms.sets)
		#
		nms.sets <- gsub("\\.grd", "", nms.sets)
		nms.sets <- gsub("\\.gri", "", nms.sets)
		nms.sets <- unique(nms.sets)
	}else{
		# f = 2 or f = 3
		nms.sets_eval <- list.files(paste0(input.eval, "/", vec_folders[f]))
		nms.sets_eval <- gsub("\\gam_eval_", "", nms.sets_eval)
		nms.sets_eval <- gsub("\\glm_eval_", "", nms.sets_eval)
		nms.sets_eval <- gsub("rf_eval_", "", nms.sets_eval)
		nms.sets_eval <- gsub("\\.csv", "", nms.sets_eval)
		nms.sets_eval <- unique(nms.sets_eval)
	}

	# LOOP ACROSS ALGORITHMS
	for(v in 1:3){

		# Loop across data sets (strategies)
		for (d in seq_along(nms.sets) ){

			# Get frames with predictive skill (each row corresponds to a taxon-specific and variable-specific model fit)
			if(f == 1){
				ls.eval <- list(); for (p in 1:length(nms.sets) ){ls.eval[[length(ls.eval) + 1]] <- read.csv(  paste0(input.eval, "/", vec_folders[f] ,"/", mod.type[v], "_eval_", nms.sets[p], ".csv")  , header = TRUE, sep = ",")}
			}else{
				ls.eval <- list(); for (p in 1:length(nms.sets_eval) ){ls.eval[[length(ls.eval) + 1]] <- read.csv(  paste0(input.eval, "/", vec_folders[f], "/", mod.type[v], "_eval_", nms.sets_eval[p], ".csv")  , header = TRUE, sep = ",")}
			}
			# Get data from brick (to avoid trouble with reloading the files) and allocate into list of 12 monthly slots with pre-gathered projections
			prj.stack <- brick(paste0(input.dir, "/", vec_folders[f], "/", mod.type[v], "_ens_", proj.type[u], "_", nms.sets[d], ".grd"))
			ls.dat.full <- list()   ;   nelem <- nlayers(prj.stack)/12   ;   from <- seq(1, nlayers(prj.stack), by = nelem)   ;   for(q in 1:12){   ls.dat.full[[q]] <- as.list(prj.stack[[ from[q]:(from[q] + nelem - 1)]])   }

			# Prepare lists to store monthly data
			lisi.months.mean <- list()

			# Loop across months
			for (k in 1:12){
				ls.dat <- ls.dat.full[[k]] # get monthly subset of data

				# PARALLEL COMPUTING ACROSS TAXA (cluster may not work due to memory constraints; try with 10 cores remotely, 4 locally for PA and 3 for PRB)
				n.cores <- 20
				cl = makeCluster(n.cores, outfile = "")
				registerDoParallel(cl)
				lisi.taxa <- foreach(i = 1 : length(unique(ls.eval[[d]]$taxon)) , .packages = c('raster', 'sfsmisc'))%dopar%{

					# Define the number of ensemble members that pass a selective skill-threshold
					len <- sum(   ls.eval[[d]][, "tss.xval4"]    [  c((i * nr - (nr - 1)):(i * nr))  ]  >= thres, na.rm = TRUE )

					if(len > 0){
						# Calculate mean projection of useful member models of taxon "i" in a monthly list, subselect members with TSS >= thres
						# From individual monthly lists, the species monthly member model projections (5 layers per species) are selected as follows:
						# i*5 - 4 (the first ensemble layer of the species) to i*5 (the last layer of the species); however from the eval sheet the corresponding positions are: i*6 - 5 (the first position) to i*6 -1 (the 5th position)
						# An additional point might be to transfom the derived ensemble means to pres if value >= 0.6, and abs if value < 0.6, yet the result is better (or at leas, smoother) when keeping the mean.
						if (len > 1){
							lisi.taxa.mean <- calc(  stack(     ls.dat    [c((i * nr - (nr - 1)):(i * nr)) [  which(ls.eval[[d]][, "tss.xval4"] [  c((i * nr - (nr - 1)):(i * nr)) ] >= thres)  ]]     ), mean, na.rm = TRUE)
							lisi.taxa.mean[lisi.taxa.mean >= 0.6] <- 1
							lisi.taxa.mean[lisi.taxa.mean < 0.6] <- 0
						}else{ # if len = 1
							lisi.taxa.mean <- ls.dat [[  c((i * nr - (nr - 1)):(i * nr))  [   which(ls.eval[[d]][,"tss.xval4"]  [c((i * nr - (nr - 1)):(i * nr))]  >= thres)  ]  ]]
						}
					}else{ # if len = 0
						lisi.taxa.mean <- na.layer
					}
					# Progress
					print(paste0(mod.type[v], " ", d , " tax", i, " month", k, " ", proj.type[u]))
					# Define output
					tax.merged <- list()
					tax.merged[[1]] <- lisi.taxa.mean

					tax.merged
				} # Close parallel computing across taxa
				stopCluster(cl)

				# Name the list elements according to the taxa
				names(lisi.taxa) <- unique(ls.eval[[d]]$taxon)

				# Store in list
				lisi.months.mean[[k]] <- lapply(lisi.taxa, '[[', 1)

				# Progress
				print(paste0(mod.type[v], " ", d, " ", proj.type[u], " month ", k, " completed..."))
			} # Close loop across months

			# Name the list elements
			names(lisi.months.mean) <- months


			# Save framed collection (as raster stack to avoid trouble with reloading the files)
			fln <- paste0(output.dir, "/", vec_folders[f], "/", mod.type[v], "_", proj.type[u], "_paBacktransformed_mean_", nms.sets[d], ".grd")		
			writeRaster(stack(unlist(lisi.months.mean)), filename = fln, format = "raster", overwrite = TRUE) 
			fln <- paste0(output.dir, "/", vec_folders[f], "/", mod.type[v], "_", proj.type[u], "_paBacktransformed_mean_", nms.sets[d], ".RData") 
			save(lisi.months.mean, file = fln)

		} # Close loop across data strategies
	} # Close loop across algorithms
} # Close loop across vec_folders (total, microscopy based, sequence based datasets)


###==============================================================
### END
###==============================================================