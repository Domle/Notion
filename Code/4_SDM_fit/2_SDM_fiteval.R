### SPECIES DISTRIBUTION MODEL (GLM, GAM, RF) - FIT AND EVALUATE

# This code works on locally stored files ("../4_Generate_Absences/Output_7" and "../8_Variable_Selection/Output_12_index_NOTION") and products are stored locally ("./Output_10_eval_NOTION")
# Approximate runtime on local machine with 6 cores in use is 20 minutes

# Author: 	Dominic Eriksson
#			Environmental Physics Group UP
#			ETH Zurich
# 			Switzerland
# 5th of September, 2022, deriksson@ethz.ch
# This script is based on a script by Daminao Righetti.

# Description: Load list with presence/pseudoabsences, with each list element containing the data-frame for an individual taxa
# Fit generalized additive models (GAM), generalized linear models (GLM), and Random Forests (RF), to species P/A data, and project the niche for each taxon*member-model on monthly-climatological global conditions across taxa*variable-members
# There are 5 background selection strategies * 2 data type strategies strategies, times three model algorithm, times 5 member-models (variable combinations) per taxon
# Test fit of model by deviance explained and adj. R-squared, test predictive skill of models by split-sample cross validation (10-fold, 8-fold, 6-fold, 4-fold, 2-fold)
#
# Input files:
# - Lists containing presences-pseudoabsences data for each taxon (n taxa = 30). Taxa with at least 14 observations have been processed, others deleted.
# - Raster objects containing environmental variables
# - Dataframe containing the ensembles of environmental predictor sets to be used in the index object
#
# Output files:
# - CSV files. Dataframe with results on evaluation for each model (multiple tss and auc (full model, and n-fold crossvalidated models), dev.expl. and adjusted R-squared (full model), obs used)

# NOTE regarding insufficient data and RF (random forest):
# The number of absences of a focal species must be at least the number of its presences (ideally even ten times the number of its presences, Barbet-Massin et al., 2012).
# The presences of a focal species may exceed its absences, due to insufficient absences cells in the target-group sampling step (see previous code on absence generation),
# In such cases, I apply a random subselection of presences to 10% below the number of absences.
# In case the presences of a focal species exceed the absesences by a factor of greater than 2, the taxon is left out from further modeling.
# Cases per strategy S:
# S1 = 0, S2 = 0, S3 = 1, S4 = 5, S5 = 5, S6 = 4, S7 = 1, S8 = 6, S9 = 8 (chaetoceros? min 18 abs...), S10 = 6  (random subselection of presences will be applied, per member model, fitted below)
# In case the presences of a focal species clearly exceed the absences (by more than a factor of 2), the species is not further modeled

### =========================================================================
### Initialize system
### =========================================================================
rm(list = ls())
lib_vec <- c("doParallel", "rJava", "dismo", "mgcv", "randomForest", "PresenceAbsence")
sapply(lib_vec, library, character.only = TRUE)

# Directories
setwd("/net/kryo/work/deriksson/Projects/Notion_DE/Code/6_SDM_fit") # Needs adjustment under different environment
#
# Input and output
input.dir <- "/home/deriksson/Projects/Notion_DE/Code/4_Generate_absences/1_Output" # Needs adjustment under different environment
output.dir <- "/home/deriksson/Projects/Notion_DE/Code/6_SDM_fit/3_Output"

# Vectors
vec_folders <- c("Total_dataset", "MicroscopyBased_dataset", "SequenceBased_dataset")
f <- 3

# Get functions for cross-validation (maintaining prevalence, i.e. "strat")
source(paste0("/home/deriksson/Projects/Notion_DE/Code/Functions/Righetti/cv.gam.strat.fun")) # gam cross validation function
source(paste0("/home/deriksson/Projects/Notion_DE/Code/Functions/Righetti/cv.glm.strat.lgit.fun")) # glm cross validation function
source("/home/deriksson/Projects/Notion_DE/Code/Functions/Righetti/ecospat.cv.rf.fun.sampsize.adj.R") # rf, adjusted (sampsize arg.), cross validation function, from ecospat package
# gam uses summary statistics for goodness of fit
source(paste0("/home/deriksson/Projects/Notion_DE/Code/Functions/Righetti/adj.D2.glm.fun")) # glm, goodness of fit function, copied from ecospat package
# rf uses oob mean out of bag error as a "goodness of fit"

### =========================================================================
### Get the data
### =========================================================================

# Strategy vector
nms.sets <- list.files( paste0(input.dir, "/", vec_folders[f]) )
nms.sets <- gsub("\\(T_MLD1).RData", "", nms.sets)

# Taxa in the list
ls.set.reference <- get(load(paste0(input.dir, "/", vec_folders[f], "/pres_abs,tot_bg_nonov(T_MLD1).RData")))
vec.no <- as.numeric(  unlist( lapply(ls.set.reference, function(x){y <- nrow(x[which(x$obs == 1), ]); return (y)})  )  )
length(vec.no[which(vec.no >= 24)]) 
length(vec.no[which(vec.no >= 14)])

# Get environmental data
env.stack <-  get (load(paste0("/net/kryo/work/deriksson/Projects/Notion_DE/Data/Predictor_variables/Environmental_data/VarSet07.2016.RData"))) # Stack with only new variables and (WOA13)

# Create NA.layer
na.layer <- raster()
na.layer[] <- NA

# Vector of algorithms used (low, intermediate, high complexity)
mod.type <- c("gam", "glm", "rf")

# Vector used for loading different index files
vec_index <- c("gam", "glm", "rf", "ensemble")

# Define number of split-samples (adjust column names in loop below)
vec.splits <- c("full", 4)

# Define threshold for species inclusion (note that taxa with < 15 observations have been previously replaced by null elements)
min.obs <- 24

### =========================================================================
###  Loop over all taxa and predefined variable combinations to fit GAMs
### =========================================================================
# Loop 
if(f == 1){index <- read.csv("/home/deriksson/Projects/Notion_DE/Code/5_Variable_selection/4_Output/Total_dataset/Species_minobs_15_index_5_members_a_4_predictors_(ensemble)_pres_abs,tot_bg_nonov.csv", header = TRUE, sep = ",")} # ensemble model approach (i.e., mean rank of gam, glm, rf single var tests)
if(f == 2){index <- read.csv("/home/deriksson/Projects/Notion_DE/Code/5_Variable_selection/4_Output/MicroscopyBased_dataset/Species_minobs_15_index_5_members_a_4_predictors_(ensemble)_microscopyBased_abs,tot_bg_non.csv", header = TRUE, sep = ",")}
if(f == 3){index <- read.csv("/home/deriksson/Projects/Notion_DE/Code/5_Variable_selection/4_Output/SequenceBased_dataset/Species_minobs_15_index_5_members_a_4_predictors_(ensemble)_sequenceBased_abs,tot_bg_non.csv", header = TRUE, sep = ",") }

# Loop across algorithms
for(v in seq_along(mod.type)){ # if error appears just run the loop below, with v manually adjusted, do not know the cause of the problem yet

	# Print progress
	print(paste0("Algorithm: ", mod.type[v]))

		# Loop across data strategies
		for (h in seq_along(nms.sets)){

		# Make the cluster setup faster, load dataset into position
		ls.sets <- list()
		d <- 1 # keep to position one, while adjusting h (to select dataset of correct data strategy)
	  	ls.sets[[d]] <- get(load(paste0(input.dir, "/", vec_folders[f], "/", nms.sets[h] , "(T_MLD1).RData"))) 

		# Parallel computing across taxa and model runs
		n.cores <- detectCores() - 10
		cl = makeCluster(n.cores, outfile = "")
		registerDoParallel(cl)
		#
		list.cl <- foreach(i = 1:nrow(index), .packages = c('raster', 'PresenceAbsence', ifelse(mod.type[v]=='gam', 'mgcv', ifelse(mod.type[v]=='rf', 'randomForest', 'sp'))))%dopar%{

		# i = 6 # CIII
		# i = 106 # Trichodesmium

			# Print progress
			print(i)
			# Observations per taxon
			if( length(ls.sets[[d]][[which(names(ls.sets[[d]]) == index[i, "taxon"])]]) > 0 ){
				ob <- length( na.omit(  ls.sets[[d]][[which(names(ls.sets[[d]]) == index[i, "taxon"])]][ which( ls.sets[[d]][[which(names(ls.sets[[d]]) == index[i, "taxon"])]]$obs == 1), c("obs")]  ) )

				# Flexibl choice of number of predictors, at least 7 degrees of freedom per predictor, specify formula (gam v == 1, glm v == 2, rf v == 3)
				if(ob < 21) {vrs <- c( as.character(index[  i,  "var1"]),  as.character(index[  i ,  "var2"]))} # bivariate model
				if(ob >= 21 && ob < 36) {vrs <- c( as.character(index[  i,  "var1"]),  as.character(index[  i ,  "var2"]), as.character(index[  i ,  "var3"]))} # three-variables per model
				if(ob >= 36) {vrs <- c( as.character(index[  i,  "var1"]),  as.character(index[  i ,  "var2"]), as.character(index[  i ,  "var3"]),  as.character(index[  i ,  "var4"]))} # four variables per model

				# --------------------------------------------------- Gam ----------------------------------------------------
				# k == 4 base dimensions (previously k == 5) to be a bit less responsive
				if (v == 1){if( length(vrs) == 2) {tt <- formula(paste0("obs~s(", vrs[1], ",k=4)+s(", vrs[2], ",k=4)"))}
				if( length(vrs) == 3) {tt <- formula(paste0("obs~s(", vrs[1], ",k=4)+s(", vrs[2], ",k=4)+s(", vrs[3], ",k=4)"))}
				if( length(vrs) == 4) {tt <- formula(paste0("obs~s(", vrs[1], ",k=4)+s(", vrs[2], ",k=4)+s(", vrs[3], ",k=4)+s(", vrs[4], ",k=4)"))}
				if( length(vrs) == 5) {tt <- formula(paste0("obs~s(", vrs[1], ",k=4)+s(", vrs[2], ",k=4)+s(", vrs[3], ",k=4)+s(", vrs[4], ",k=4)+s(", vrs[5], ",k=4)"))}
				if( length(vrs) == 6) {tt <- formula(paste0("obs~s(", vrs[1], ",k=4)+s(", vrs[2], ",k=4)+s(", vrs[3], ",k=4)+s(", vrs[4], ",k=4)+s(", vrs[5], ",k=4)+s(", vrs[6], ",k=4)"))}}
				# --------------------------------------------------- Glm vars ----------------------------------------------------
				if (v == 2){if( length(vrs) == 2) {tt <- paste0("obs~", vrs[1], "+I(", vrs[1], "^2)+", vrs[2], "+I(", vrs[2], "^2)")}
				if( length(vrs) == 3) {tt <- paste0("obs~", vrs[1], "+I(", vrs[1], "^2)+", vrs[2], "+I(", vrs[2], "^2)+", vrs[3], "+I(", vrs[3], "^2)")}
				if( length(vrs) == 4) {tt <- paste0("obs~", vrs[1], "+I(", vrs[1], "^2)+", vrs[2], "+I(", vrs[2], "^2)+", vrs[3], "+I(", vrs[3], "^2)+", vrs[4], "+I(", vrs[4], "^2)")}
				if( length(vrs) == 5) {tt <- paste0("obs~", vrs[1], "+I(", vrs[1], "^2)+", vrs[2], "+I(", vrs[2], "^2)+", vrs[3], "+I(", vrs[3], "^2)+", vrs[4], "+I(", vrs[4], "^2)+", vrs[5], "+I(", vrs[5], "^2)")}
				if( length(vrs) == 6) {tt <- paste0("obs~", vrs[1], "+I(", vrs[1], "^2)+", vrs[2], "+I(", vrs[2], "^2)+", vrs[3], "+I(", vrs[3], "^2)+", vrs[4], "+I(", vrs[4], "^2)+", vrs[5], "+I(", vrs[5], "^2)", vrs[6], "+I(", vrs[6], "^2)")}}
				# --------------------------------------------------- RF vars ----------------------------------------------------
				if (v == 3){if( length(vrs) == 2) {tt <- paste0("as.factor(obs)~", vrs[1], "+", vrs[2])}
				if( length(vrs) == 3) {tt <- paste0("as.factor(obs)~", vrs[1], "+", vrs[2], "+", vrs[3])}
				if( length(vrs) == 4) {tt <- paste0("as.factor(obs)~", vrs[1], "+", vrs[2], "+", vrs[3], "+", vrs[4])}
				if( length(vrs) == 5) {tt <- paste0("as.factor(obs)~", vrs[1], "+", vrs[2], "+", vrs[3], "+", vrs[4], "+", vrs[5])}
				if( length(vrs) == 6) {tt <- paste0("as.factor(obs)~", vrs[1], "+", vrs[2], "+", vrs[3], "+", vrs[4], "+", vrs[5], "+", vrs[6])}}
				# ------------------------------------------------------------------------------------------------------------
				}

			# Open else condition, address taxa with zero obs or too few obs regarding specified variables
			if( length(ls.sets[[d]][[which(names(ls.sets[[d]]) == index[i, "taxon"])]]) == 0  ) {
				vec.eval <- rep(NA, 15) # Taxa with zero obs, note: nrow does not work (returns logical, whilst length(...) does the job)
			}else{
				# Exclude taxa for which the presences exceed the absences by a factor of more than two
				if( length(which(ls.sets[[d]][[which(names(ls.sets[[d]]) == index[i, "taxon"])]]$obs == 1)) / length(which(ls.sets[[d]][[which(names(ls.sets[[d]]) == index[i, "taxon"])]]$obs == 0)) > 2 ) {
					# layers <- stack(replicate(24, na.layer))
					vec.eval <- rep(NA, 15)
				}else{
					if(  nrow(na.omit(ls.sets[[d]][[which(names(ls.sets[[d]]) == index[i, "taxon"])]][ which( ls.sets[[d]][[which(names(ls.sets[[d]]) == index[i, "taxon"])]]$obs == 1), c(vrs, "obs")])) < min.obs ){
						vec.eval <- rep(NA, 15)
					}else{

						# Prepare data
						dat <- ls.sets[[d]][[which(names(ls.sets[[d]]) == index[i, "taxon"])]]
						dat <- dat[ complete.cases(dat[, which(  names(dat) %in% vrs)]), ] # Glm (unlike gam) and RF have trouble with NA's in step-wise reduction (rows change), rows with NA's are removed (pertaining to variables of interest)
						rownames(dat) <- 1:nrow(dat) # Rownames need to be clean, else RF may have trouble with cv.test
						dat <- cbind(dat, weights_dec = 1) # Apparently gam has problems with the notation "weights" (a function) in cross-val, we call it weights_dec

						# For reproducibility, needs to be consistent with model evaluation code
						set.seed(99)

						# Make sure that there are at least 10% more absences than presences in the data: else subsample/thin the presences at random
						if( length(which(dat$obs == 1))/length(which(dat$obs == 0)) > 0.9 ) {
							dat <- dat[ c(sample(   which(dat$obs == 1) , round(length(which(dat$obs == 0))*0.9,0) ), which(dat$obs == 0)) , ]
							rownames(dat) <- 1:nrow(dat) # Rownames need to be clean, else RF may have trouble with cv.test
						}

						# Adjust weights after dropping points (product of weights times points shall be equal; usually presences = 1, absences = 0.1 - as ten times more absences had been generated)
						weight_abs <-length(which(dat$obs == 1))/length(which(dat$obs == 0)) # New ratio of presence count divided by absences counts, in case some presences or absences were lost due to variable constraints
						dat$weights_dec[dat$obs == 0] <- weight_abs

						# --------------------------------------------------- Gam ----------------------------------------------------
						if (mod.type[v] == "gam"){
							gam_full <- mgcv::gam(tt, family="binomial", data = dat, weights = weights_dec, select = FALSE)
						}
						# --------------------------------------------------- Glm ----------------------------------------------------
						if (mod.type[v] == "glm"){
							glm_prec <- glm( eval(parse(text = tt)),family = "binomial", data = dat, weights = weights_dec)
							glm_full <- step(glm_prec, direction = "both", trace = FALSE)
						} # Stepwise variable reduction (removing non-significant terms using AIC as selection criterion)
						# --------------------------------------------------- Rf (sampsize balanced) ----------------------------------------------------
						if(mod.type[v] == "rf"){
							if(weight_abs > 1){
								rf_full <- randomForest( eval(parse(text = tt)), data = dat, importance = TRUE, nodesize = 1, ntree = 4000, sampsize = c( nrow(dat[which(dat$obs != 1), ]), nrow(dat[which(dat$obs == 1), ]))) # ntree = 3000 is currently implemented in ecospat.cv.rf.adj function; however 4000 might be optimal
							}else{rf_full <- randomForest( eval(parse(text = tt)), data = dat, importance = TRUE, nodesize = 1, ntree = 4000, sampsize = c( nrow(dat[which(dat$obs == 1),]) , nrow(dat[which(dat$obs == 1), ])))}
						}
						# ------------------------------------------------------------------------------------------------------------
						# Evaluate model - calculate TSS and AUC by using k-fold cross validation (for full and leave-out data trained models)
						vec.auc <- vector('numeric')
						vec.tss <- vector('numeric')
						# --------------------------------------------------- Gam ----------------------------------------------------
						if (mod.type[v] == "gam"){
							for(k in 1:length(vec.splits)) { #loop: data of differential partitioning (train data on k-1 splits, predict on the left-out split, k-times, each split once predicted, test model skill)
								if (vec.splits[k] == "full"){
									gam.xval_test <- data.frame("id" = 1:length(gam_full$y), "observed" = gam_full$y, "predicted" = gam_full$fitted) #full model
								}else{
									gam_xval <- cv.gam.strat(dat, gam_full, as.numeric(vec.splits[k]), lvo.cv = FALSE) # cross validation gam
									gam.xval_test <- data.frame("id" = 1:length(gam_full$y),"observed" = gam_full$y, "predicted" = gam_xval[complete.cases(gam_xval), ]$predicted) # test frame
								}
								paa <- presence.absence.accuracy(gam.xval_test, threshold = optimal.thresholds(gam.xval_test)[3,2], st.dev = FALSE) # tss optimizing binarization threshold
								vec.auc[length(vec.auc) + 1] <- round(paa[, 7], 4) # auc
								vec.tss[length(vec.tss) + 1] <- round(paa[, 4] + paa[, 5] - 1, 4) # tss
							}
						vec.eval <- c(vec.auc, vec.tss)
						vec.eval[length(vec.eval) + 1] <- round(summary(gam_full)$dev.exp, 4) # dev expl. (full model)
						vec.eval[length(vec.eval) + 1] <- round(summary(gam_full)$r.sq, 4) # adj. R-squared (full model)
						vec.eval[length(vec.eval) + 1] <- length(gam_full$y[gam_full$y == 1])} # number of presence observations used for gam run
						
						# --------------------------------------------------- Glm ----------------------------------------------------
						if (mod.type[v] == "glm"){
						for(k in 1:length(vec.splits)){
							if (vec.splits[k] == "full"){
								glm.xval_test <- data.frame("id" = 1:length(glm_full$y),"observed" = glm_full$y,"predicted" = glm_full$fitted) # full model
							}else{
								glm_xval <- cv.glm.strat.lgit(dat, glm_full, as.numeric(vec.splits[k]), lvo.cv = FALSE) # cross validation glm
								glm.xval_test <- na.omit(data.frame("id" = 1:length(glm_full$y), "observed" = glm_full$y, "predicted" = glm_xval$predicted)) # New: Set to na.omit as NAs cause trouble
							}
							paa <- presence.absence.accuracy(glm.xval_test, threshold = optimal.thresholds(glm.xval_test)[3, 2], st.dev = FALSE)
							vec.auc[length(vec.auc) + 1] <- round(paa[, 7], 4)
							vec.tss[length(vec.tss) + 1] <- round(paa[, 4] + paa[, 5] - 1, 4)
						}
						vec.eval <- c(vec.auc, vec.tss)
						vec.eval[length(vec.eval) + 1] <- round(summary(glm_full)$aic, 4) # akaike information criterion (full model)
						vec.eval[length(vec.eval) + 1] <- round(adj.D2.glm(glm_full), 4) # adj. D-squared function (full model)
						vec.eval[length(vec.eval) + 1] <- length(glm_full$y[glm_full$y == 1])} #??Add: number of used observations for glm run
						
						# ---------------------------------------------------- Rf -----------------------------------------------------
						if (mod.type[v] == "rf"){
							for(k in 1:length(vec.splits)){
								if (vec.splits[k] == "full"){
									rf.xval_test <- data.frame("id" = 1:length(rf_full$y),"observed" = as.numeric(as.character(rf_full$y)), "predicted" = as.numeric(predict(rf_full, type = "prob")[, 2]))
								}else{
								rf.xval <- ecospat.cv.rf.adj(rf.obj = rf_full, data.cv = data.frame("obs" = dat[ , c("obs")], dat[ , which(names(dat) %in% vrs)]), cv.lim = 1, K = as.numeric(vec.splits[k]), jack.knife = FALSE) # cross validation RF, adjusted function from ecospat package
								rf.xval_test <- data.frame("id" = 1:length(rf_full$y),"observed" = as.numeric(as.character(rf_full$y)), "predicted" = rf.xval$predictions)
							}
								paa <- presence.absence.accuracy(rf.xval_test, threshold = optimal.thresholds(rf.xval_test)[3,2], st.dev = FALSE)
								vec.auc[length(vec.auc) + 1]<-round(paa[, 7], 4)
								vec.tss[length(vec.tss) + 1]<-round(paa[, 4] + paa[, 5] - 1, 4)
							}
						vec.eval <- c(vec.auc, vec.tss)
						vec.eval[length(vec.eval) + 1] <- round(mean(rf_full$er[, "OOB"]), 4) # there is NO deviance expl. (full model) for RF; i take the mean OOB (mean out of bag error) instead
						vec.eval[length(vec.eval) + 1] <- NA # There is no ajusted r-squared for RF.
						vec.eval[length(vec.eval) + 1] <- length(rf_full$y[rf_full$y == 1]) # number of presence observations used for rf run
						}
						# ------------------------------------------------------------------------------------------------------------

						# Close else conditions
	    			}
	    		}
	    	}

	   		# Progress (internal logbook)
	   		if (length(ls.sets[[d]][[which(names(ls.sets[[d]]) == index[i,"taxon"])]]) == 0) {
	   			print(paste0(i, " Set", h, " ", as.character(index[i, 1]), " taxon ", index[i, 2][[1]], " : ", mod.type[v], " NA vector produced"))
	   		}else{
	   			if (nrow(na.omit(ls.sets[[d]][[which(names(ls.sets[[d]]) == index[i, "taxon"])]][ which( ls.sets[[d]][[which(names(ls.sets[[d]]) == index[i, "taxon"])]]$obs == 1), c(vrs, "obs")])) < min.obs){
	   				print(paste0(i, " Set", h, " ", as.character(index[i, 1]), " taxon ", index[i, 2][[1]], " : ", mod.type[v], "  NA vector produced, too few obs, < ", min.obs))
	   			}else{
	   				print(paste0 (i, " Set", h, " ", as.character(index[i, 1]), " taxon ", index[i, 2][[1]]," : ", mod.type[v], " ", paste(vrs[!vrs %in% "NA"], collapse = "-"), " n.obs=",
	   				nrow(na.omit(ls.sets[[d]][[which(names(ls.sets[[d]]) == index[i,"taxon"])]][ which( ls.sets[[d]][[which(names(ls.sets[[d]]) == index[i, "taxon"])]]$obs == 1), c(vrs, "obs")]))
	   				," ..tss.4xval = ", vec.eval[4])) # fourth position for tss.4-fold.xval , in case we have c("full", 4) in the "vec.splits"
	   			}
			}

	   		# Define output
			vec.eval

		# Close parallel computing across taxa
		}

		# Stop cluster
		stopCluster(cl)

		# Merge output data
		df <- data.frame(cbind(
			# index[,c(2,3)][1:length(list.cl),], # taxon, group
			index[, c("taxon", "group")][1:length(list.cl),], # taxon, group
			do.call("rbind", list.cl)[1:length(list.cl), 1:(2 * (length(vec.splits)) + 3)], # results, curtailed to useful length of eval.vec (as the NA's usef as eval.vecs for unsuitable species are longer with length == 15 )
			paste0( apply(data.frame(index[, c("var1", "var2", "var3", "var4")][1:length(list.cl),]), 2, as.character)[ , 1],  "-" , apply(data.frame(index[,c("var1", "var2", "var3", "var4")][1:length(list.cl),]), 2, as.character)[ , 2], "-" , apply(data.frame(index[,c("var1", "var2", "var3", "var4")][1:length(list.cl),]), 2, as.character)[ , 3], "-" , apply(data.frame(index[, c(4:7)][1:length(list.cl),]), 2, as.character)[ , 4]   ))) # variables
		vars_few_obs <- paste0( apply(data.frame(index[,c("var1", "var2", "var3", "var4")][1:length(list.cl),]), 2, as.character)[ ,1],  "-" , apply(data.frame(index[,c("var1", "var2", "var3", "var4")][1:length(list.cl),]), 2, as.character)[ , 2], "-" , apply(data.frame(index[,c("var1", "var2", "var3", "var4")][1:length(list.cl),]), 2, as.character)[ , 3]  ) # variables for species with rather few observations
		vars_very_few_obs <- paste0( apply(data.frame(index[,c("var1", "var2", "var3", "var4")][1:length(list.cl),]), 2, as.character)[ ,1],  "-" , apply(data.frame(index[,c("var1", "var2", "var3", "var4")][1:length(list.cl),]), 2, as.character)[ , 2]) # variables for species with very few observations

		# Add column-names
		names(df) <- c("taxon", "group",
		"auc.full","auc.xval4",
		"tss.full","tss.xval4",
		"Dev.expl","adj.Rsq", "obs.used", "vars")

		# Correct variable entries
		df$vars <- as.character(df$vars)
		df[which(df$obs.used >= 21 & df$obs.used < 36) , "vars" ] <- vars_few_obs[which(df$obs.used >= 21 & df$obs.used < 36)]
		df[which(df$obs.used < 21) , "vars" ] <- vars_very_few_obs[which(df$obs.used < 21)]

		if(v==3) {names(df) <- c("taxon", "group",
		"auc.full","auc.xval4",
		"tss.full","tss.xval4",
		"avg.OOB.error","adj.Rsq", "obs.used", "vars")}

		# Save data
		fln <- paste0(output.dir,"/", vec_folders[f], "/", mod.type[v],"_eval_", nms.sets[h], ".csv")
		write.csv(df,file = fln, row.names = FALSE)


		} # Close loop over data strategies

} # Close loop across algorithms


