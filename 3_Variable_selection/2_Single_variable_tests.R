### UNIVARIATE MODEL TESTS FOR VARIABLE RANKING (SUBSEQUENT SELECTION) USING GENERAL LINEAR MODEL, GENERAL ADDITIVE MODEL AND RANDOM FOREST

# Input files:	1. Diazotroph dataset containing presences and pseudoabsences and matched up environmental variables.

# Output files: 1. CSV file containing the goodness of fit for each variable, algorithm and taxon.

# Strategy: Single glm, gam and RF-based variable test.



# Authors: 	Dominic Eriksson
#			Environmental Physics Group, UP
#			ETH Zurich
#			Switzerland


# Script developed by Damiano Righetti

# deriksson@ethz.ch, 4th of September 2023 --------------------------------------------------------


### =========================================================================
### Preparation
### =========================================================================
# Clear work space
rm(list = ls())

# Libraries
lib_vec <- c("doParallel", "doBy", "mgcv", "randomForest")
new_packages <- lib_vec[!(lib_vec %in% installed.packages()[, "Package"])] # check which package is not installed
if(length(new_packages) > 0 ) {install.packages(new_packages)} # install missing packages
sapply(lib_vec, library, character.only = TRUE) # load packages

# Directories
wd_var_selection <- "/net/gladiolus/backups/lunaria/work--2021_05_25_17.00.01--complete/rdamiano/legacy/Righetti_et_al_2019_sciadv/Data/8_Variable_Selection/"
wd_dat <- "/net/kryo/work/deriksson/Projects/Notion_DE/Code/4_Generate_absences/1_Output/"
output.dir <- "/net/kryo/work/deriksson/Projects/Notion_DE/Code/5_Variable_selection/2_Output/"


## Functions to test goodness of fit of single variable models
# gam uses statistics for goodness of fit obtained via summary(gam_full)$r.sq
source("/net/kryo/work/deriksson/Projects/Notion_DE/Code/Functions/Righetti/adj.D2.glm.fun", chdir = TRUE)
# rf uses oob mean out of bag error as a "goodness of fit" obtained via mean(rf_full$er[, "OOB"])


### =========================================================================
### Load lists with presence-pseudoabsence data, split by species (takes time)
### =========================================================================


# Enter loop here
vec_folders <- list.files(wd_dat)
for(v in seq_along(vec_folders)){

	# Print progress
	print(paste0("Dataset: ", vec_folders[v]))

	# Overview: There are 3 background selection strategies (tot, gr, cr) * 2 background overlap strategies (overlapping, nonoverlapping backgrounds, meaning species modeled included or excluded from bg) * 2 data type strategies (pres-abs, counts-abs) yielding 12 data scenarios
	# Define data types * background strategies
	nms.sets <- gsub("\\(T_MLD1).RData", "", list.files( paste0(wd_dat, vec_folders[v]) ))

	# Variables with multicollinearity less than 0.7 after using vifcor function of usdm package (VIF Analysis)
	vc.vars<-c("x","y","T","Sal","N","P","Si",
			"MLD1", "MLD2","PAR","Chl","Wind.CCMP","pCO2",
			"MLPAR1", "MLPAR2", "Nstar", "Sistar",
			"dT_dt", "dN_dt", "dP_dt", "dSi_dt","dMLD1_dt",
			"logMLD1", "logMLD2", "logChl", "logN", "logP", "logSi")
			# I remove NPP, as well reflected in Chl


	# Select one data version (not containing thinned data), to get the list of all taxa contained!
	# taxa.ref.list <- get(load(paste0(wd_dat, vec_folders[v],"/",   list.files( paste0(wd_dat, vec_folders[v]))   [6]))) # gridded, tot-bg, homogenised taxa, overlapping sites, T-MLD1 => Standard approach

	# Define the minimum of presences (obs) required per taxon
	min.obs <- 14

	# # Get an index for the selected taxa with minimum number of observarions or higher
	# index.useful <- which(as.numeric(unlist(sapply(taxa.ref.list, function(x){x1 <- NA; if(!is.null(x)){x1 <- nrow(x[which(x$obs == 1), ])}; return(x1)}))) >= min.obs)
	# length(index.useful)

	### =========================================================================
	###  Single variable test
	### =========================================================================

	# Loop across data strategies
	for( d in seq_along(nms.sets) ){

		# Print progress
		print(d)

		# Load data
		taxa.lists <- list()
		# Load data based on different background sampling strategies
		for (l in seq_along(nms.sets)){taxa.lists[[length(taxa.lists) + 1]] <- get(load(paste0(wd_dat, vec_folders[v], "/", nms.sets[l], "(T_MLD1).RData"  )))}

		# Get list with taxa
		taxa.ref.list <- taxa.lists[[d]]
		# Get an index for the selected taxa with minimum number of observarions or higher
		index.useful <- which(as.numeric(unlist(sapply(taxa.ref.list, function(x){x1 <- NA; if(!is.null(x)){x1 <- nrow(x[which(x$obs == 1), ])}; return(x1)}))) >= min.obs)
		length(index.useful)
		#
		taxa.ls <- taxa.lists[[d]]

		# Exclude elements below min.obs based on observation points in the standard taxa list
		taxa.list <- taxa.ls[index.useful] # contains only taxa with more than minimum observations
						
		# Parallel computing across taxa
		n.cores <- 10
		cl = makeCluster(n.cores, outfile = "")
		registerDoParallel(cl)
		list.cl <- foreach(i = c(1:length(taxa.list)), .packages = c("doBy", "mgcv", "randomForest"))%dopar%{

			# Single variable glm test:
			# Filter taxa with zero obs or less than the minimum abs, create NA's
			if( length(taxa.list[[i]]) == 0 | nrow(as.data.frame(taxa.list[[i]][which(taxa.list[[i]]$obs == 0),])) < min.obs ){
			
				# Create vectors to store results
				vec.glm <- rep(NA, length(vc.vars))
				vec.gam <- rep(NA, length(vc.vars))
				vec.rf <- rep(NA, length(vc.vars))
				#
				ranks.glm <- rep(NA, length(vc.vars))
				ranks.gam <- rep(NA, length(vc.vars))
				ranks.rf <- rep(NA, length(vc.vars))
				obs <- c(NA, NA)
			}else{

				# Filter taxa with obs < minimum obs --> Unneccassary step, is done right befor this step above
				if(  length(taxa.list[[i]]) > 0 & nrow(taxa.list[[i]][which(taxa.list[[i]]$obs == 1),]) < min.obs ) {
				vec.glm <- rep(NA, length(vc.vars))
				vec.gam <- rep(NA, length(vc.vars))
				vec.rf <- rep(NA, length(vc.vars))
				#
				ranks.glm <- rep(NA, length(vc.vars))
				ranks.gam <- rep(NA, length(vc.vars))
				ranks.rf <- rep(NA, length(vc.vars))
				obs <- c(NA, NA)
			} else {

				# Useful taxa: create vector to insert value for each variable k
				vec.glm <- rep(NA, length(vc.vars))
				vec.gam <- rep(NA, length(vc.vars))
				vec.rf <- rep(NA, length(vc.vars))

				# Loop across variables
				for( k in seq_along(vc.vars) ){

					# Get data: removing NA's with regard to the variable tested
					dat <- taxa.list[[i]][ which(is.na(taxa.list[[i]][, vc.vars[k]]) == FALSE), ]
					if( nrow(dat[which(dat$obs == 1),]) < min.obs ){ # we want at least min.obs points to fit a model
						vec.gam[k] <- NA
						vec.glm[k] <- NA
						vec.rf[k] <- NA
						next
					}
						
					# Re-adjust the weighting of the data: not entirely clear to me if weighting is really needed. As adjusted D-squared and adjusted R-squared take into account the sample size?
					dat$weights_dec <- 1
					weight_abs <- length(which(dat$obs == 1)) / length(which(dat$obs == 0)) # New ratio of presence count divided by absences counts, in case some presences or absences were lost due to variable constraints
					dat$weights_dec[dat$obs == 0] <- weight_abs
						
					
#-------------------------------------------------------------- FIT MODELS ---------------------------------------------------------------
					# fit glm and gam: weighted (note: the weighted models yield more comparable and similar results between gam and glm)
					mod.gam <- gam(formula(paste0("obs~s(", vc.vars[k], ",k=5)")), data = dat, na.action = na.omit, family = binomial, weights = weights_dec)
					mod.glm <- glm( dat$obs ~ dat[  ,  vc.vars[k] ] + I(  (dat[,  vc.vars[k] ])^2), data = dat, na.action = na.omit, family = binomial, weights = weights_dec)

					# fit random forest
					if( weight_abs > 1 ){mod.rf <- randomForest( eval(parse(text = paste0("as.factor(obs)~", vc.vars[k]))), data = dat, importance = TRUE, nodesize = 1, ntree = 300, sampsize = c( nrow(dat[which(dat $obs != 1),]) , nrow(dat[which(dat$obs == 1),])))
					}else{mod.rf <- randomForest( eval(parse(text = paste0("as.factor(obs)~", vc.vars[k]))), data = dat, importance = TRUE, nodesize = 1, ntree = 300, sampsize = c( nrow(dat[which(dat$obs == 1), ]) , nrow(dat[which(dat$obs == 1), ])))}
						# evaluate models (simply via: goodness of fit) - no cross validation is implemented here
						vec.gam[k] <- round(summary(mod.gam)$r.sq, 4)
						vec.glm[k] <- round(adj.D2.glm(mod.glm), 4)
						vec.rf[k] <- round(1 - mean(mod.rf$er[, "OOB"]), 4)
					} # end of loop across variables

				# Assemble the results - Note we use a minus here, because the rank function ranks from small to big, but our best ranking starts from high to low
				# Note - If many value equal zero and no deviance is explained by the predictors they all get the same rank example Atelocyanobacterium, cruise specific overlapping glm has nearly all zeros as Dsquared
				ranks.glm <- c(NA, NA, rank(-vec.glm[ 3:length(vec.glm)], ties.method = c("max")))
				ranks.gam <- c(NA, NA, rank(-vec.gam[ 3:length(vec.gam)], ties.method = c("max")))
				ranks.rf <- c(NA, NA, rank(-vec.rf[ 3:length(vec.rf)], ties.method = c("max")))
				obs <- rep( nrow( dat[which(dat$obs > 0), ] ), 2)
				} # end of else condition
			} # end of else condition

			## Merge test metrics (adjusted R-squared of model), ranking of test metric, obs, stratification-way, taxon name, group name and NA's row for better structure
			# GAM
			df.gam <- data.frame(
				rbind(vec.gam, ranks.gam, c(rep(NA, length(vec.gam)))),
				"obs" = c(obs, NA),
				"pa_strat" = c(rep(nms.sets[d], 2), "NA"), # Important: * Add the correct name of the data strategy *
				"value" = c("adj.Rsq", "rank", "NA"),
				"taxon" = c(rep(names(taxa.list)[i], 2), "NA"),
				"group" = c(rep(as.character(taxa.list[[i]]$group[1]), 2), "NA"))
			colnames(df.gam)[1:length(vc.vars)] <- vc.vars
			# Print progress
			print(paste(i, "gam"))

			# GLM
			df.glm <- data.frame(
				rbind(vec.glm, ranks.glm, c(rep(NA, length(vec.glm)))),
				"obs" = c(obs, NA),
				"pa_strat"= c(rep(nms.sets[d], 2), "NA"), # Important: * Add the correct name of the data strategy *
				"value" = c("adj.Dsq", "rank", "NA"),
				"taxon" = c(rep(names(taxa.list)[i], 2), "NA"),
				"group" = c(rep(as.character(taxa.list[[i]]$group[1]), 2), "NA"))
			colnames(df.glm)[1:length(vc.vars)] <- vc.vars
			# Print progress
			print(paste(i, "glm"))

			# RF
			df.rf <- data.frame(
				rbind(vec.rf, ranks.rf, c(rep(NA, length(vec.rf)))),
				"obs" = c(obs, NA),
				"pa_strat"= c(rep(nms.sets[d], 2), "NA"), # Important: * Add the correct name of the data strategy *
				"value" = c("1-OOB.error", "rank", "NA"),
				"taxon" = c(rep(names(taxa.list)[i], 2), "NA"),
				"group" = c(rep(as.character(taxa.list[[i]]$group[1]), 2), "NA"))
			colnames(df.rf)[1:length(vc.vars)] <- vc.vars 
			# Print progress
			print(paste(i, "rf"))

			# Define output
			list.mod <- list()
			list.mod[[1]] <- df.gam
			list.mod[[2]] <- df.glm
			list.mod[[3]] <- df.rf
			list.mod # defines output of the cluster run i that is fed into the list called "list.cl"

		} # Close parallel computing across taxa
		stopCluster(cl)

		#-----------------------------------------------------------------------------------------------------------------------------
		# Save species results: gam
		df.new <- as.data.frame(do.call("rbind", lapply(list.cl, "[[", 1)))
		df.new <- cbind(df.new[, 1:2], round(df.new[ , 3:(ncol(df.new) - 4)], 3), df.new[, (ncol(df.new) - 3):ncol(df.new)]) # round values (optional)
		df.new <- df.new[, c(vc.vars, "obs",  "taxon", "group", "pa_strat", "value")]
		# Save
		fln <- paste0(output.dir, vec_folders[v], "/Taxa_Gam_", nms.sets[d], "_minobs_", min.obs, ".csv")
		write.csv(df.new, file = fln, row.names = FALSE)

		# Save species results: glm
		df.new <- as.data.frame(do.call("rbind", lapply(list.cl, "[[", 2)))
		df.new <- cbind(df.new[, 1:2], round(df.new[ ,3:(ncol(df.new)-4)], 3), df.new[,(ncol(df.new)-3):ncol(df.new)]) # round values (optional)
		df.new <- df.new[, c(vc.vars, "obs",  "taxon", "group", "pa_strat", "value")]
		# Save
		fln <- paste0(output.dir, vec_folders[v], "/Taxa_Glm_", nms.sets[d], "_minobs_", min.obs, ".csv")
		write.csv(df.new, file = fln, row.names = FALSE)

		# Save species results: rf
		df.new <- as.data.frame(do.call("rbind", lapply(list.cl, "[[", 3)))
		df.new <- cbind(df.new[, 1:2], round(df.new[ ,3:(ncol(df.new)-4)], 3), df.new[,(ncol(df.new)-3):ncol(df.new)]) # round values (optional)
		df.new <- df.new[, c(vc.vars, "obs",  "taxon", "group", "pa_strat", "value")]
		# Save
		fln <- paste0(output.dir, vec_folders[v], "/Taxa_Rf_", nms.sets[d], "_minobs_", min.obs,".csv")
		write.csv(df.new, file = fln, row.names = FALSE)
	} # Close loop across data scenarios
} # close loop over vec_folders

###==============================================================
### END
###==============================================================