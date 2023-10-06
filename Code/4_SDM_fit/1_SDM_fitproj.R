### SPECIES DISTRIBUTION MODEL (GLM, GAM, RF) - FIT MODEL and PROJECT PROBABILITY OF PRESENCE AND PRESENCE-ABSENCE ON GLOBAL MAP


# Input files:	1. Dataframe with taxonomic information, longitude and latitude and presence + pseudoabsences encoded as 1 and 0's.
#				2. Grid file containing environmental variables
# 				3. Dataframe containing the ensemble information on environmental predictors used for one species, to create and index object.

# Output files: 1. 
#				2. R data files that contain visreg objects containing the fits of each environmental variable across all ensemble members (saved as individual files)
#               3. CSV file containing the merged information of each environmental predictor fit across all ensemble members


# Strategy: Load list with presence/pseudoabsences, with each list element containing the data-frame for an individual taxa
# Fit generalized additive models (GAM), generalized linear models (GLM), and Random Forests (RF), to species P/A data, and project the niche for each taxon*member-model on monthly-climatological global conditions across taxa*variable-members
# There are 5 background selection strategies * 2 data type strategies strategies, times three model algorithm, times 5 member-models (variable combinations) per taxon
# Input: Lists containing presences-pseudoabsences data for each taxon (n taxa = 49, of which 30 contain more than 14 observations used here as lower threshold for modeling)
# Output: Lists addressing each model type (times data strategy): list contains taxa*variable-combinations (n total = 900) as elements (12 prb-layers, 12 pa-layers)


# NOTE regarding insufficient data and RF (random forest):
# The number of absences of a focal species must be at least the number of its presences (ideally even ten times the number of its presences, Barbet-Massin et al., 2012).
# The presences of a focal species may exceed its absences, due to insufficient absences cells in the target-group sampling step (see previous code on absence generation),
# In such cases, I apply a random subselection of presences to 10% below the number of absences.
# In case the presences of a focal species exceed the absesences by a factor of greater than 2, the taxon is left out from further modeling.
# We will merge the individual files containing the fit of the response curve into one table.

# Script developed by Damiano Righetti

# deriksson@ethz.ch, 19 of September 2023 ------------------------------------------------------------------------------------


### =========================================================================
### Preparation
### =========================================================================
lib <- c("doParallel", "rJava", "dismo", "mgcv", "randomForest", "PresenceAbsence", "colorRamps", "RColorBrewer", "visreg")
sapply(lib, library, character.only = TRUE)

# Input and output
input.dir <- "/home/deriksson/Projects/Notion_DE/Code/4_Generate_absences/1_Output/Total_dataset" # Needs adjustment under different environment
output.dir <- "/home/deriksson/Projects/Notion_DE/Code/6_SDM_fit/2_Output/"
wd_out_responseCurve <- '/net/kryo/work/deriksson/Projects/Notion_DE/Code/6_SDM_fit/2_750_Output/Total_dataset/ResponseCurves/'

# Vectors
vec_folders <- c("Total_dataset", "MicroscopyBased_dataset", "SequenceBased_dataset")
f <- 1 # SWITCH
#
nms.sets <- list.files(input.dir)
nms.sets <- gsub("\\(T_MLD1).RData", "", nms.sets)

### =========================================================================
### Get the data
### =========================================================================
# Taxa in the list
ls.set.reference <- get(load(paste0(input.dir, "/pres_abs,tot_bg_overl(T_MLD1).RData")))
vec.no <- as.numeric(  unlist( lapply(ls.set.reference, function(x){y <- nrow(x[which(x$obs==1), ]); return (y)})  )  )
length(vec.no[which(vec.no >= 24)])
length(vec.no[which(vec.no >= 14)])

# Get environmental data
prj.stack <- brick("/net/kryo/work/deriksson/Projects/Notion_DE/Data/Predictor_variables/Environmental_data/VarSet07.2016.grd")
nelem <- nlayers(prj.stack)/12
nm <- gsub("\\.1", "", names(prj.stack[[1:nelem]]))
env.stack<-list() ; from<-seq(1,nlayers(prj.stack),by=nelem) ; for(q in 1:12){pre.env.stack <- stack(prj.stack[[ from[q]:(from[q]+nelem-1)]]); names(pre.env.stack) <- nm; env.stack[[q]] <- pre.env.stack}
names(env.stack) <- c("JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC")

# Create NA.layer
na.layer <- raster()
na.layer[] <- NA

# Prepare ensemble member models (variable choice)
if(f == 1){index <- read.csv("/home/deriksson/Projects/Notion_DE/Code/5_Variable_selection/4_Output/Total_dataset/Species_minobs_15_index_5_members_a_4_predictors_(ensemble)_pres_abs,tot_bg_nonov.csv", header = TRUE, sep = ",")} # ensemble model approach (i.e., mean rank of gam, glm, rf single var tests)
if(f == 2){index <- read.csv("/home/deriksson/Projects/Notion_DE/Code/5_Variable_selection/4_Output/MicroscopyBased_dataset/Species_minobs_15_index_5_members_a_4_predictors_(ensemble)_microscopyBased_abs,tot_bg_non.csv", header = TRUE, sep = ",")}
if(f == 3){index <- read.csv("/home/deriksson/Projects/Notion_DE/Code/5_Variable_selection/4_Output/SequenceBased_dataset/Species_minobs_15_index_5_members_a_4_predictors_(ensemble)_sequenceBased_abs,tot_bg_non.csv", header = TRUE, sep = ",") }
# Vector of algorithms used (low, intermediate, high complexity)
mod.type <- c("gam", "glm", "rf")

# Define threshold for species inclusion (note that taxa with < 14 observations have been previously replaced by null elements)
min.obs <- 24

### =========================================================================
###  Loop over all taxa and predefined variable combinations to fit models
### =========================================================================

# Loop across algorithms
for(v in 1:3){

	# Loop across data strategies
	for (h in seq_along(nms.sets)){

	# Prepare ensemble member models (variable choice)
	if(f == 1){index <- read.csv(paste0("/home/deriksson/Projects/Notion_DE/Code/5_Variable_selection/4_Output/Total_dataset/Species_minobs_15_index_5_members_a_4_predictors_(ensemble)_", nms.sets[h], ".csv"), header = TRUE, sep = ",")} # ensemble model approach (i.e., mean rank of gam, glm, rf single var tests)
	if(f == 2){index <- read.csv(paste0("/home/deriksson/Projects/Notion_DE/Code/5_Variable_selection/4_Output/MicroscopyBased_dataset/Species_minobs_15_index_5_members_a_4_predictors_(ensemble)_", nms.sets[h], ".csv"), header = TRUE, sep = ",")}
	if(f == 3){index <- read.csv(paste0("/home/deriksson/Projects/Notion_DE/Code/5_Variable_selection/4_Output/SequenceBased_dataset/Species_minobs_15_index_5_members_a_4_predictors_(ensemble)_", nms.sets[h], ".csv"), header = TRUE, sep = ",")}
	

	# Make the cluster setup faster, load dataset into position
	ls.sets <- list()
	d <- 1 # keep to position one, while adjusting h (to select dataset of correct data strategy)
	ls.sets[[1]] <- get(load(paste0(input.dir, "/", list.files(input.dir)[h])))


	# Parallel computing across taxa and model runs
	n.cores <- detectCores() - 10
	cl = makeCluster(n.cores, outfile = "")
	registerDoParallel(cl)
	list.cl <- foreach(i = 1:nrow(index), .packages = c('raster', 'PresenceAbsence', 'visreg', ifelse(mod.type[v] == 'gam', 'mgcv', ifelse (mod.type[v]=='rf', 'randomForest', 'sp'))))%dopar%{

		# Observations per taxon
		if( length(ls.sets[[d]] [[which(names(ls.sets[[d]]) == index[i,"taxon"])]]) > 0 ){
			ob <- length( na.omit(  ls.sets[[d]] [[which(names(ls.sets[[d]]) == index[i,"taxon"])]] [ which( ls.sets[[d]] [[which(names(ls.sets[[d]]) == index[i, "taxon"])]]$obs == 1), c("obs")]  ) )
			# Flexible choice of number of predictors, at least 7 degrees of freedom per predictor, specify formula (gam v == 1, glm v == 2, rf v == 3)
			if(ob < 21) {vrs <- c( as.character(index[  i ,  "var1"]),  as.character(index[  i ,  "var2"]))} # bivariate model
			if(ob >= 21 && ob < 36) {vrs <- c( as.character(index[  i ,  "var1"]),  as.character(index[  i ,  "var2"]), as.character(index[  i ,  "var3"]))} # three-variables per model
			if(ob >= 36) {vrs <- c( as.character(index[  i ,  "var1"]),  as.character(index[  i ,  "var2"]), as.character(index[  i ,  "var3"]),  as.character(index[  i ,  "var4"]))} # four variables per model

			# NOTE: Be carefull of the response curves you are fitting, do they make sense!!!!
			# --------------------------------------------------- Gam ----------------------------------------------------
			# k == 4 base dimensions (previously k == 5) to be a bit less responsive
			if (v == 1){if( length(vrs) == 2) {tt <- formula(paste0("obs~s(", vrs[1], ",k=4)+s(", vrs[2], ",k=4)"))}
			if( length(vrs) == 3) {tt <- formula(paste0("obs~s(", vrs[1], ",k=4)+s(", vrs[2], ",k=4)+s(", vrs[3], ",k=4)"))}
			if( length(vrs) == 4) {tt <- formula(paste0("obs~s(", vrs[1], ",k=4)+s(", vrs[2], ",k=4)+s(", vrs[3], ",k=4)+s(", vrs[4], ",k=4)"))}
			if( length(vrs) == 5) {tt <- formula(paste0("obs~s(", vrs[1], ",k=4)+s(", vrs[2], ",k=4)+s(", vrs[3], ",k=4)+s(", vrs[4], ",k=4)+s(", vrs[5], ",k=4)"))}
			if( length(vrs) == 6) {tt <- formula(paste0("obs~s(", vrs[1], ",k=4)+s(", vrs[2], ",k=4)+s(", vrs[3], ",k=4)+s(", vrs[4], ",k=4)+s(", vrs[5], ",k=4)+s(", vrs[6], ",k=4)"))}}
			# --------------------------------------------------- Glm ----------------------------------------------------
			if (v == 2){if( length(vrs) == 2) {tt <- paste0("obs~", vrs[1], "+I(", vrs[1], "^2)+", vrs[2], "+I(", vrs[2], "^2)")}
			if( length(vrs) == 3) {tt <- paste0("obs~", vrs[1], "+I(", vrs[1], "^2)+", vrs[2], "+I(", vrs[2], "^2)+", vrs[3], "+I(", vrs[3], "^2)")}
			if( length(vrs) == 4) {tt <- paste0("obs~", vrs[1], "+I(", vrs[1], "^2)+", vrs[2], "+I(", vrs[2], "^2)+", vrs[3], "+I(", vrs[3], "^2)+", vrs[4], "+I(", vrs[4],"^2)")}
			if( length(vrs) == 5) {tt <- paste0("obs~", vrs[1], "+I(", vrs[1], "^2)+", vrs[2], "+I(", vrs[2], "^2)+", vrs[3], "+I(", vrs[3], "^2)+", vrs[4], "+I(", vrs[4],"^2)+", vrs[5], "+I(", vrs[5], "^2)")}
			if( length(vrs) == 6) {tt <- paste0("obs~", vrs[1], "+I(", vrs[1], "^2)+", vrs[2], "+I(", vrs[2], "^2)+", vrs[3], "+I(", vrs[3], "^2)+", vrs[4], "+I(", vrs[4],"^2)+", vrs[5], "+I(", vrs[5], "^2)", vrs[6], "+I(", vrs[6], "^2)")}}
			# --------------------------------------------------- Rf ----------------------------------------------------
			if (v == 3){if( length(vrs) == 2) {tt <- paste0("as.factor(obs)~", vrs[1], "+", vrs[2])}
			if( length(vrs) == 3) { tt <- paste0("as.factor(obs)~", vrs[1], "+", vrs[2], "+", vrs[3])}
			if( length(vrs) == 4) {tt <- paste0("as.factor(obs)~", vrs[1], "+", vrs[2], "+", vrs[3], "+", vrs[4])}
			if( length(vrs) == 5) {tt <- paste0("as.factor(obs)~", vrs[1], "+", vrs[2], "+", vrs[3], "+", vrs[4], "+", vrs[5])}
			if( length(vrs) == 6) {tt <- paste0("as.factor(obs)~", vrs[1], "+", vrs[2], "+", vrs[3], "+", vrs[4], "+", vrs[5], "+", vrs[6])}}
			# ------------------------------------------------------------------------------------------------------------
		}

		# Open else condition, address taxa with zero obs or too few obs regarding specified variables
		if( length(ls.sets[[d]][[which(names(ls.sets[[d]]) == index[i, "taxon"])]]) == 0  ){layers <- stack(replicate(24, na.layer)) # Taxa with Null entries, note: nrow does not work (returns logical), whilst length(...) does the job
		}else{
			
			# Exclude taxa for which the presences exceed the absences by a factor of more than two
			if( length(which(ls.sets[[d]][[which(names(ls.sets[[d]]) == index[i, "taxon"])]]$obs == 1)) / length(which(ls.sets[[d]][[which(names(ls.sets[[d]]) == index[i, "taxon"])]]$obs == 0)) > 2 ) {layers <- stack(replicate(24, na.layer))
			}else{

				if(  nrow( na.omit(ls.sets[[d]][[which(names(ls.sets[[d]]) == index[i, "taxon"])]][ which( ls.sets[[d]][[which(names(ls.sets[[d]]) == index[i, "taxon"])]]$obs == 1), c(vrs, "obs")]) ) < min.obs   )  {layers <- stack(replicate(24, na.layer))
				}else{

					# Prepare data
					dat <- ls.sets[[d]][[which(names(ls.sets[[d]]) == index[i,"taxon"])]]
					dat <- dat[ complete.cases(dat[, which(  names(dat) %in% vrs)]), ] # Glm (unlike gam) has trouble with NA's in step-wise reduction (rows change), rows with NA's regarding selected variables are removed
					rownames(dat) <- 1:nrow(dat) # Rownames need to be clean, else RF may have trouble with cv.test
					dat <- cbind(dat, weights_dec = 1) # Gam has problems with the notation "weights" (a function) in cross-val, we call it weights_dec

					# For reproducibility, needs to be consistent with model evaluation code
					set.seed(99)

					# Make sure that there are at least 10% more absences than presences in the data: else subsample/thin the presences at random
					if( length(which(dat$obs == 1)) / length(which(dat$obs == 0)) > 0.9) {
					dat <- dat[ c(sample(   which(dat$obs == 1) , round(length(which(dat$obs == 0))*0.9,0) ), which(dat$obs == 0)) , ]
					rownames(dat) <- 1:nrow(dat) # Rownames need to be clean, else RF may have trouble with cv.test
					}

					# Adjust weights after dropping points (product of weights times points shall be equal; usually presences = 1, absences = 0.1 - as ten times more absences had been generated)
					weight_abs <-length(which(dat$obs == 1))/length(which(dat$obs == 0)) # New ratio of presence count divided by absences counts, in case some presences or absences were lost due to variable constraints
					dat$weights_dec[dat$obs == 0] <- weight_abs

					# --------------------------------------------------- Gam ----------------------------------------------------
					if (mod.type[v] == "gam"){
					mod_full <- mgcv::gam(tt, family = "binomial", data = dat, weights = weights_dec, select = FALSE)
					mod_visreg <- try(visreg(mod_full, plot = F, scale = 'response'), silent = TRUE)
					}
					# --------------------------------------------------- Glm ----------------------------------------------------
					if (mod.type[v] == "glm"){
					glm_prec <- glm( eval(parse(text = tt)),family = "binomial", data = dat, weights = weights_dec)
					mod_full <- step(glm_prec, direction = "both", trace = FALSE)
					mod_visreg <- try(visreg(mod_full, plot = F, scale = 'response'), silent = TRUE)
					} # Stepwise variable reduction (removing non-significant terms using AIC as selection criterion)
					# --------------------------------------------------- Rf (sampsize balanced) ----------------------------------------------------
					# ** Note may be corrupt in case not enough absences are there. Adjust in absences selection script!!! Something does not work with v = 3, and d = 3,4,5
					if(mod.type[v] == "rf"){
						if(weight_abs > 1){
							# mod_full <- randomForest( eval(parse(text = tt)), data = dat, importance = TRUE, nodesize = 1, ntree = 4000, sampsize = c( nrow(dat[which(dat$obs != 1),]) , nrow(dat[which(dat$obs == 1),])))
							mod_full <- randomForest( eval(parse(text = tt)), data = dat, importance = TRUE, nodesize = 10, ntree = 750, sampsize = c( nrow(dat[which(dat$obs != 1),]) , nrow(dat[which(dat$obs == 1),])))
							mod_visreg <- try(visreg(mod_full, plot = F), silent = TRUE)
						}else{
							mod_full <- randomForest( eval(parse(text = tt)), data = dat, importance = TRUE, nodesize = 10, ntree = 750, sampsize = c( nrow(dat[which(dat$obs == 1),]) , nrow(dat[which(dat$obs == 1), ])))} # ntree = 3000 is currently implemented in ecospat.cv.rf.adj function; however 4000 might be optimal
							mod_visreg <- try(visreg(mod_full, plot = F), silent = TRUE)
					}
					# ------------------------------------------------------------------------------------------------------------


					# Save response curve data
					fln_responseCurve <- paste0(wd_out_responseCurve, 'ResponseCurveDataList_', index[i, 'taxon'], "_", nms.sets[h], '_', mod.type[v], '_', index[i, 'var1'], '_', index[i, 'var2'], '_', index[i, 'var3'], '_', index[i, 'var4'], '.RData')
					saveRDS(mod_visreg, fln_responseCurve)
					rm(mod_visreg)

					# Project: probabiliy of presence - PRB
					lisi.proj <- list()
					for (k in 1:12){
						if(mod.type[v] == "gam" | mod.type[v] == "glm") {lisi.proj[[k]] <- predict(env.stack[[k]], mod_full, type = "response") }
						if(mod.type[v] == "rf") {lisi.proj[[k]] <- predict(env.stack[[k]], mod_full, type = "prob", index = 2)} # RF warrants type == "prob" for a probabiliscitc response, index = 2 selects prob of presence (index = 1, absence)
					}

					# Project: binary projection - PA (by defining True Skill Statistic - TSS - optimizing threshold to convert probability to presence-absence (note; random forest is not binarized by internal default threshold)
					# Are we using the default option threshold 0.5?
					if (mod.type[v] == "gam") {aa <- optimal.thresholds(data.frame("id" = 1:length(mod_full$y), "observed" = mod_full$y, "predicted" = mod_full$fitted))[3,2]}
					if (mod.type[v] == "glm") {aa <- optimal.thresholds(data.frame("id" = 1:length(mod_full$y), "observed" = mod_full$y, "predicted" = mod_full$fitted))[3,2]}
					if (mod.type[v] == "rf") {aa <- optimal.thresholds(data.frame("id" = 1:length(mod_full$y), "observed" = as.numeric(as.character(mod_full$y)), "predicted" = as.numeric(predict(mod_full, type = "prob")[, 2])))[3, 2]}
					for (k in 1:12){
						lisi.proj[[k + 12]] <- reclassify(lisi.proj[[k]],  matrix(c(0.0, aa, 0, aa, 1.0, 1), ncol = 3 , byrow = TRUE))
					}

					# Merge both types of projections in stack
					layers <- stack(lisi.proj)
				} # Close else conditions
			}
    	}

   		# Progress (internal logbook)
    	if (length(ls.sets[[d]][[which(names(ls.sets[[d]]) == index[i,"taxon"])]]) == 0){
    		print(paste0(i," Set", h, " ", as.character(index[i, 1]), " taxon ", index[i, 2][[1]], " : ", mod.type[v], " NA.layers produced"))
		}else{
    		print(paste0(i," Set", h, " ", as.character(index[i, 1]), " taxon ", index[i, 2][[1]], " : ", mod.type[v], " ",  paste(vrs[!vrs %in% "NA"], collapse = "-")))
    	}

    	# Define output
    	layers

    # Close parallel computing across taxa
    }
	# Stop cluster
	stopCluster(cl)

	# list.cl is a list of length 140. We modelled 28 taxa time 5 predictor sets, 28*5=140.
	# Each list element is a raster with 24 layers (2 time twelve month, one time for probability output; HSI, and another time 
	# for the presence absence converted data).
	d <- list.cl[[1]]

	# Save data
	fln <- paste0(output.dir_rf_750, vec_folders[f], "/", mod.type[v], "_proj_", nms.sets[h], ".RData")
	save(list.cl, file = fln)

	# Close loop over data strategies
	}

# Close loop across algorithms
}


### ===========================================================================
### Format response curve data into one big table 
###============================================================================

# Libraries
library(tidyr)
library(stringr)

# Get filenames
file_names <- list.files('/home/deriksson/Projects/Notion_DE/Code/6_SDM_fit/2_750_Output/Total_dataset/ResponseCurves')
file_names <- file_names[1:(length(file_names) - 1)]


# Loop
df_final <- data.frame() # YOU NEED TO ADD ALSO THE PREDICTOR SET IN THE ID TO NOT DOUBLE THE LINES
for(f in seq_along(file_names)){

	# Print progress
	print( paste0('Opening file number ', f, ', from ', length(file_names)) )
	
	# Open file
	l <- readRDS( paste0('/home/deriksson/Projects/Notion_DE/Code/6_SDM_fit/2_750_Output/Total_dataset/ResponseCurves/', file_names[f]) )
	
	# Number of variables included
	var_length <- length(l)
    # Skip to next if less than two predictors were included
    if(var_length <= 1) {
        cat('Skipping iteration')
        next
    }
    # Extract variable names - Sometimes the element to be extracted is called 'fit' sometimes 'visregFit'
    if( 'fit' %in% names(l[[1]]) ){
        var <- names( l[[1]]$fit )[1:var_length]
    }
    # We skip this --> Not working not fully understood either, why the parameter to extract changes although same function used    
    if( 'visregFit' %in% names(l[[1]]) ) {
        cat('Skipping iteration')
        next
    }

    # Name the list elements according to the predictor testes
    names(l) <- var

    ## Extract predictor fit
    # 4 Variables
    if(var_length == 4){
        
        ## Save info in data frame
        # --------------------------------- Variable 1 ----------------------------------
        df <- l[[var[1]]]$fit[var[1]]
        names(df) <- 'measurement'
        # Add variable name
        df$variable <- var[1]
        # Add fit
        df$visregFit <- l[[var[1]]]$fit$visregFit
        df$visregLwr <- l[[var[1]]]$fit$visregLwr
        df$visregUpr <- l[[var[1]]]$fit$visregUpr
        # --------------------------------- Variable 2 ----------------------------------
        df_inter <- l[[var[2]]]$fit[var[2]]
        names(df_inter) <- 'measurement'
        # Add variable name
        df_inter$variable <- var[2]
        # Add fit
        df_inter$visregFit <- l[[var[2]]]$fit$visregFit
        df_inter$visregLwr <- l[[var[2]]]$fit$visregLwr
        df_inter$visregUpr <- l[[var[2]]]$fit$visregUpr
        # Merge
        df <- rbind(df, df_inter)
        # --------------------------------- Variable 3 ----------------------------------
        df_inter <- l[[var[3]]]$fit[var[3]]
        names(df_inter) <- 'measurement'
        # Add variable name
        df_inter$variable <- var[3]
        # Add fit
        df_inter$visregFit <- l[[var[3]]]$fit$visregFit
        df_inter$visregLwr <- l[[var[3]]]$fit$visregLwr
        df_inter$visregUpr <- l[[var[3]]]$fit$visregUpr
        # Merge
        df <- rbind(df, df_inter)
        # --------------------------------- Variable 4 ----------------------------------
        df_inter <- l[[var[4]]]$fit[var[4]]
        names(df_inter) <- 'measurement'
        # Add variable name
        df_inter$variable <- var[4]
        # Add fit
        df_inter$visregFit <- l[[var[4]]]$fit$visregFit
        df_inter$visregLwr <- l[[var[4]]]$fit$visregLwr
        df_inter$visregUpr <- l[[var[4]]]$fit$visregUpr

        # Merge
        df <- rbind(df, df_inter)

        ## Add additional information
        id <- str_extract( file_names[f], "ResponseCurveDataList_(.*?)_rf")
        if( grepl('glm', file_names[f]) ){
            id <- str_extract( file_names[f], "ResponseCurveDataList_(.*?)_glm")
        }
        if( grepl('gam', file_names[f]) ){
            id <- str_extract( file_names[f], "ResponseCurveDataList_(.*?)_gam")
        }
        df$id <- gsub('ResponseCurveDataList_', '', id)

        # Merge into final data frame
        df_final <- rbind(df_final, df)

    } # end of if condition 4 variables
    
    # 3 Variables
    if(var_length == 3){
        
        ## Save info in data frame
        # --------------------------------- Variable 1 ----------------------------------
        df <- l[[var[1]]]$fit[var[1]]
        names(df) <- 'measurement'
        # Add variable name
        df$variable <- var[1]
        # Add fit
        df$visregFit <- l[[var[1]]]$fit$visregFit
        df$visregLwr <- l[[var[1]]]$fit$visregLwr
        df$visregUpr <- l[[var[1]]]$fit$visregUpr
        #--------------------------------- Variable 2 ----------------------------------
        df_inter <- l[[var[2]]]$fit[var[2]]
        names(df_inter) <- 'measurement'
        # Add variable name
        df_inter$variable <- var[2]
        # Add fit
        df_inter$visregFit <- l[[var[2]]]$fit$visregFit
        df_inter$visregLwr <- l[[var[2]]]$fit$visregLwr
        df_inter$visregUpr <- l[[var[2]]]$fit$visregUpr
        # Merge
        df <- rbind(df, df_inter)
        # --------------------------------- Variable 3 ----------------------------------
        df_inter <- l[[var[3]]]$fit[var[3]]
        names(df_inter) <- 'measurement'
        # Add variable name
        df_inter$variable <- var[3]
        # Add fit
        df_inter$visregFit <- l[[var[3]]]$fit$visregFit
        df_inter$visregLwr <- l[[var[3]]]$fit$visregLwr
        df_inter$visregUpr <- l[[var[3]]]$fit$visregUpr

        # Merge
        df <- rbind(df, df_inter)

        ## Add additional information
        id <- str_extract( file_names[f], "ResponseCurveDataList_(.*?)_rf")
        if( grepl('glm', file_names[f]) ){
            id <- str_extract( file_names[f], "ResponseCurveDataList_(.*?)_glm")
        }
        if( grepl('gam', file_names[f]) ){
            id <- str_extract( file_names[f], "ResponseCurveDataList_(.*?)_gam")
        }
        df$id <- gsub('ResponseCurveDataList_', '', id)

        # Merge into final data frame
        df_final <- rbind(df_final, df)

    } # end of if condition 3 variables
    
    
    # 2 Variables
    if(var_length == 2){
        
        ## Save info in data frame
        # --------------------------------- Variable 1 ----------------------------------
        df <- l[[var[1]]]$fit[var[1]]
        names(df) <- 'measurement'
        # Add variable name
        df$variable <- var[1]
        # Add fit
        df$visregFit <- l[[var[1]]]$fit$visregFit
        df$visregLwr <- l[[var[1]]]$fit$visregLwr
        df$visregUpr <- l[[var[1]]]$fit$visregUpr
        # --------------------------------- Variable 2 ----------------------------------
        df_inter <- l[[var[2]]]$fit[var[2]]
        names(df_inter) <- 'measurement'
        # Add variable name
        df_inter$variable <- var[2]
        # Add fit
        df_inter$visregFit <- l[[var[2]]]$fit$visregFit
        df_inter$visregLwr <- l[[var[2]]]$fit$visregLwr
        df_inter$visregUpr <- l[[var[2]]]$fit$visregUpr

        # Merge
        df <- rbind(df, df_inter)

        ## Add additional information
        id <- str_extract( file_names[f], "ResponseCurveDataList_(.*?)_rf")
        if( grepl('glm', file_names[f]) ){
            id <- str_extract( file_names[f], "ResponseCurveDataList_(.*?)_glm")
        }
        if( grepl('gam', file_names[f]) ){
            id <- str_extract( file_names[f], "ResponseCurveDataList_(.*?)_gam")
        }
        df$id <- gsub('ResponseCurveDataList_', '', id)

        # Merge into final data frame
        df_final <- rbind(df_final, df)
        
    } # end of if condition 2 variables
    
    
    
    # 1 Variables
    if(var_length == 1){
        
        ## Save info in data frame
        # --------------------------------- Variable 1 ----------------------------------
        df_inter <- l[[var[1]]]$fit[var[1]]
        names(df_inter) <- 'measurement'
        # Add variable name
        df_inter$variable <- var[1]
        # Add fit
        df_inter$visregFit <- l[[var[1]]]$fit$visregFit
        df_inter$visregLwr <- l[[var[1]]]$fit$visregLwr
        df_inter$visregUpr <- l[[var[1]]]$fit$visregUpr

        ## Add additional information
        id <- str_extract( file_names[f], "ResponseCurveDataList_(.*?)_rf")
        if( grepl('glm', file_names[f]) ){
            id <- str_extract( file_names[f], "ResponseCurveDataList_(.*?)_glm")
        }
        if( grepl('gam', file_names[f]) ){
            id <- str_extract( file_names[f], "ResponseCurveDataList_(.*?)_gam")
        }
        df_inter$id <- gsub('ResponseCurveDataList_', '', id)

        # Merge into final data frame
        df_final <- rbind(df_final, df_inter)        

    } # end of if condition 3 variables

} # end of loop


## Formatting
# Harmonize species names
df_final$id <- gsub('Scytonematopsis_pilosa', 'Scytonematopsis pilosa', df_final$id)
# Split columns
df_final <- separate_wider_delim(df_final, cols = id, delim = '_', names = c('taxon', 'str1', 'str2', 'str3', 'str4', 'algorithm'))

# Save table
fln <- paste0(wd_out_responseCurve, "/", "Table_ResponseCurveData_AcrossAllEnsembleMembers.csv")
write.csv(df_final, fln, row.names = FALSE)

###==============================================================
### END
###==============================================================