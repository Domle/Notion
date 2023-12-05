### CREATE INDEX DEFINING VARIABLE CHOICE ; BY RANDOMIZING VARIABLES WITHIN AN ENSEMBLE STRUCTURE - TAKING A BALANCED APPROACH OF VARIABLE CHOICE (TAXON BY TAXON)

# Input files:  1. Global spearman correlation table
#               2. Dataframe with merged results of the variable ranking of all taxa.

# Output files: 1. CSV file containing the predictor sets for each taxon.

# Strategy: The index defines five models à 4 different predictors, each, to build an SDM ensemble
#               a) avoiding correlations > |0.7| among variables within members
#               b) not selecting the same variable more than twice across all members (i.e. balancing potential variable contribution)
#               c) selecting variables from the top ten, and if needed, subsequent rank

# We will also use another strategy, always keeping the three most important predictors and only replacing the last one.

# Author: 	Dominic Eriksson
#			Environmental Physics Group, UP
#			ETH, Zurich
#			Switzerland

# Scritp developed by Damiano Righetti

# deriksson@ethz.ch, 19th of September 2023 --------------------------------------------------------------


### =========================================================================
### Preparation
### =========================================================================

# Directories
input.dir <- "/net/kryo/work/deriksson/Projects/Notion_DE/Code/5_Variable_selection/3_Output/"
output.dir <- "/net/kryo/work/deriksson/Projects/Notion_DE/Code/5_Variable_selection/4_Output/"

# Preparatory: data frame with the correlates for each variable
my.corr  <- read.csv("/home/deriksson/Projects/Notion_DE/Code/5_Variable_selection/1_Output/Corr_table_global_spearman.csv")
names(my.corr)
rownames(my.corr) <- names(my.corr[2:46])
my.corr$X <- NULL


# Vector for background selection strategies
vec_pseudoAbs <- c(

	# Microscopy based dataset
	"microscopyBased_abs,cr_bg_nonov",
	"microscopyBased_abs,cr_bg_overl",
	"microscopyBased_abs,gr_bg_nonov",
	"microscopyBased_abs,gr_bg_overl",
	"microscopyBased_abs,tot_bg_non",
	"microscopyBased_abs,tot_bg_overl",

	# Sequence based dataset
	"sequenceBased_abs,cr_bg_nonov",
	"sequenceBased_abs,cr_bg_overl",
	"sequenceBased_abs,gr_bg_nonov",
	"sequenceBased_abs,gr_bg_overl",
	"sequenceBased_abs,tot_bg_non",
	"sequenceBased_abs,tot_bg_overl",

	# Total dataset
	"pres_abs,cr_bg_nonov",
	"pres_abs,cr_bg_overl",
	"pres_abs,gr_bg_nonov",
	"pres_abs,gr_bg_overl",
	"pres_abs,tot_bg_non",
	"pres_abs,tot_bg_overl"
)

# Vector for sdms 
proj.type <- 'ensemble'
#
vec_folders <- c("MicroscopyBased_dataset", "SequenceBased_dataset", "Total_dataset")


### =========================================================================
### A. Spearman correlation based. Sites - based approach. Select variables for the first best model across species
### =========================================================================
for( v in seq_along(vec_folders)){

	print(v)

	# Get filenames
	# vec_fname <- list.files( paste0(input.dir, vec_folders[v]) )
	vec_fname <- list.files( paste0(input.dir, vec_folders[v]) )

	# Loop through
	for(d in seq_along(vec_fname)){

		# Print progress 
		print(paste0("Working on background selection strategy ", vec_fname[d], "."))

		# for(p in proj.type){

			# Load data
			dat.ran <- read.csv( paste0(input.dir, vec_folders[v], "/", vec_fname[d]) )

			# Set seed (for reproducibility)
			set.seed(67)

			# LOOP ACROSS TAXA
			lisi.vars <- list() # the variables selected per species
			lisi.frames <- list() # for the statistics
			for( s  in  1:length(na.omit(unique(dat.ran$taxon))) ){

				# LOOP ACROSS COMPOSED MEMBERS (VARIABLE COMBINATIONS): select five times 4 variables (without putting variables back)
                # The ensemble mean rank is rank is used and based on single variable tests using a glm, gam and RF
				dat.spp.raw <- dat.ran[           max(which(dat.ran$taxon == na.omit(unique(dat.ran$taxon))[s] )) + 1,   3: (ncol(dat.ran) - 5) ]
				#
				dat.spp <- c(  apply(   dat.spp.raw ,   2, function(x){as.numeric(x)  }   )  )
				df.var.raw <- data.frame("var" = names(dat.spp), "rank" = as.numeric(as.character(dat.spp)))
				df.var <- df.var.raw[order(df.var.raw$rank), ]
				rownames(df.var) <- 1:nrow(df.var)
				df.var$check.vector <- 0
				df.var$weight <- 1

				# Weighting of variables used for the initial selection (in case both log and non log versions of variables, and MLD1 and MLD2, is contained among the top ten ranked variables); thus temperature and physical variables count equal as nutrient variables and log versions of nutrients initially
				if(  sum(c("MLD1", "MLD2") %in% df.var[1:10, ]$var) == 2 ) {   df.var[ which(df.var$var == "MLD1") ,   "weight"] <- 0.5;  df.var[ which(df.var$var == "MLD2") ,   "weight"] <- 0.5         }
				if(  sum(c("MLD1", "logMLD1") %in% df.var[1:10, ]$var) == 2 ) {   df.var[ which(df.var$var == "MLD1") ,   "weight"] <- 0.5;  df.var[ which(df.var$var == "logMLD1") ,   "weight"] <- 0.5         }
				if(  sum(c("MLD2", "logMLD2") %in% df.var[1:10, ]$var) == 2 ) {   df.var[ which(df.var$var == "MLD2") ,   "weight"] <- 0.5;  df.var[ which(df.var$var == "logMLD2") ,   "weight"] <- 0.5         }
				if(  sum(c("MLD2", "logMLD1") %in% df.var[1:10, ]$var) == 2 ) {   df.var[ which(df.var$var == "MLD2") ,   "weight"] <- 0.5;  df.var[ which(df.var$var == "logMLD1") ,   "weight"] <- 0.5         }
				if(  sum(c("logMLD2", "logMLD1") %in% df.var[1:10, ]$var) == 2 ) {   df.var[ which(df.var$var == "logMLD2") ,   "weight"] <- 0.5;  df.var[ which(df.var$var == "logMLD1") ,   "weight"] <- 0.5         }
				if(  sum(c("MLPAR1", "MLPAR2") %in% df.var[1:10, ]$var) == 2 ) {   df.var[ which(df.var$var == "MLPAR1") ,   "weight"] <- 0.5;  df.var[ which(df.var$var == "MLPAR2") ,   "weight"] <- 0.5         }
				if(  sum(c("P", "logP") %in% df.var[1:10, ]$var) == 2 ) {   df.var[ which(df.var$var == "P") ,   "weight"] <- 0.5;  df.var[ which(df.var$var == "logP") ,   "weight"] <- 0.5         }
				if(  sum(c("N", "logN") %in% df.var[1:10, ]$var) == 2 ) {   df.var[ which(df.var$var == "N") ,   "weight"] <- 0.5;  df.var[ which(df.var$var == "logN") ,   "weight"] <- 0.5         }
				if(  sum(c("Si", "logSi") %in% df.var[1:10, ]$var) == 2 ) {   df.var[ which(df.var$var == "Si") ,   "weight"] <- 0.5;  df.var[ which(df.var$var == "logSi") ,   "weight"] <- 0.5         }
				if(  sum(c("Chl", "logChl") %in% df.var[1:10, ]$var) == 2 ) {   df.var[ which(df.var$var == "Chl") ,   "weight"] <- 0.5;  df.var[ which(df.var$var == "logChl") ,   "weight"] <- 0.5         }

				# Increase the weight of the top five predictors in order to consider them preferentially for initial selections of members
				#
				# df.var[ 1:5, "weight"] <- df.var [ 1:5, "weight"]*4 # Disabled for letter as differences are not huge

				lisi.members <- list()
				for(u in 1:5){ 	# We make 5 sets a 4 predictors

					# Prepare vector for member
					my.choice <- vector()

					# SELECT FIRST VARIABLE AT RANDOM (from the available ones); store in my.choice, with equal weights to log transformed plus normal nutrients (or mld, and chl) relative to other variables
					my.potential.vars <- df.var[ which(df.var[1:10,]$check.vector < 2), ] # ** CRITICAL DECISION, TOP 10
					my.choice[1] <- as.character(my.potential.vars[sample(nrow(my.potential.vars), 1, prob = my.potential.vars$weight  ), "var"]) # sample() takes a random sample 
					df.var[   which(my.choice[1] == df.var$var), "check.vector" ] <- df.var[   which(my.choice[1] == df.var$var) , "check.vector"  ] + 1 # Adjust the check.vector field

					# Loop across the remaining variables:
					while(   length(my.choice)  < 4  ){  # ** CRITICAL DECISION: NUMBER OR PREDICTORS

						# The remaining positions that potentially provide one predictor to the member are
						valid.positions <- df.var[which(df.var$check.vector < 2 & df.var$rank < 11), ] # ** CRITICAL DECISION, TOP 10

						# The already chosen variable(s) are correlated with
						df.corr.relevant <- my.corr[     which(rownames(my.corr) %in% my.choice),     ] # The useful part of the data frame
						columns.with.values.greater.zero.point.seven <- which(apply(df.corr.relevant, 2, max, na.rm = TRUE) > 0.7)
						variables.correlated <- names(columns.with.values.greater.zero.point.seven) # These variables are correlated with the already chosen variable

						# The remaining positions that potentially provide one predictor to the memer are hence
						if (   sum(variables.correlated %in% valid.positions$var) == 0 ) {valid.pos <- valid.positions } else { # Case: none of the variables correlated with already chosen variables is found in the potential variables
						valid.pos <- valid.positions[   -which(valid.positions$var %in% variables.correlated)   ,   ]} #

							# CASE A: there are positions left among the top *ten* ranked variables
							if( nrow(valid.pos) > 0 ) {
								my.choice[length(my.choice) + 1] <- as.character(valid.pos[sample(nrow(valid.pos), 1), "var"]) # select the variable
								df.var[   which(my.choice[length(my.choice)] == df.var$var), "check.vector" ] <- df.var[   which(my.choice[length(my.choice)] == df.var$var) , "check.vector"  ] + 1 # Adjust the check.vector field
								}

							# CASE B: no positions are left : loop across the lower ranked variables
							if( nrow(valid.pos) == 0 ) {
								# We start with the variable that is ranked 11 ** CRITICAL DECISION, TOP 10
								i = 11
								while(  length(my.choice) < 4  ){
									var.potential <- as.character(df.var[i, 1])
									print(i)

									if ( df.var[i, "check.vector"] > 1 ){ i = i + 1; next} # Each variable shall not be selected more than twice
										var.sel.pot <- c(my.choice, var.potential)  # test of correlation: the selected variables plus the new potential variable
										corr.subset <-  my.corr  [    which(   rownames(my.corr) %in% var.sel.pot   )    ,    which(   colnames(my.corr) %in% var.sel.pot   )  ]
										corr.subset [cbind(1:dim(corr.subset)[2], 1:dim(corr.subset)[2])]	<- NA # Replace diagonal elements by NA (but not all "1")
									if(   any ( corr.subset  > 0.7, na.rm = TRUE   ) == TRUE ) {i = i + 1; next} else {
										i = i + 1 
										my.choice[length(my.choice) + 1] <- var.potential
										df.var[   which(my.choice[length(my.choice)] == df.var$var), "check.vector" ] <- df.var[   which(my.choice[length(my.choice)] == df.var$var) , "check.vector"  ] + 1 # Adjust the check.vector field
									}
								}
							}
							# close across remaining variables
					}
					# Store member in list
					lisi.members[[u]] <- my.choice
				} # END LOOP ACROSS COMPOSED MEMBERS (VARIABLE COMBINATIONS) (u)

				# CREATE FRAME
				df.members <- data.frame("taxon" =  rep(as.character(na.omit(unique(dat.ran$taxon))[s]), 5) ,
				"group" =  rep( as.character(dat.ran[   max(which(dat.ran$taxon == na.omit(unique(dat.ran$taxon))[s] ))    ,   ]$group)    , 5) ,
				"var1" = unlist(lapply (lisi.members, '[', 1)),
				"var2" = unlist(lapply (lisi.members, '[', 2)),
				"var3" = unlist(lapply (lisi.members, '[', 3)),
				"var4" = unlist(lapply (lisi.members, '[', 4)))

				# END LOOP ACROSS TAXA
				lisi.vars[[length(lisi.vars)+1]] <- df.members
				lisi.frames[[length(lisi.frames)+1]] <- df.var
				lisi.frames$"Emiliania huxleyi"

				print(s)
			} # end loop over species

			# Add names to the frames
			names(lisi.frames) <- na.omit(unique(dat.ran$taxon))
			# df.frames <- do.call(rbind, lisi.frames)

			# Assemble variables index frame
			index <- cbind(do.call(rbind, lisi.vars), "id" = 1:nrow(do.call(rbind, lisi.vars)))

			# A few stats (ensemble mean approach) * NOT UP TO DATE *
			nrow(   index[ which(index$var1 == "T"), ]   ) +  nrow(   index[ which(index$var2 == "T"), ]   ) + nrow(   index[ which(index$var3 == "T"), ]   )  +  nrow(   index[ which(index$var4 == "T"), ]   ) # 949 (8 fold weighting top 5);         932 (4 fold weighting top 5);         902 (even weights)       -> glm: 932 (4 fold)      -> gam: 911    -> rf:  899 (4 fold)
			nrow(   index[ which(index$var1 == "P"), ]   ) +  nrow(   index[ which(index$var2 == "P"), ]   ) + nrow(   index[ which(index$var3 == "P"), ]   )  +  nrow(   index[ which(index$var4 == "P"), ]   ) # 533 (8 fold weighting top 5);         534 (4 fold weighting top 5);         518 (even weights)       -> glm: 557 (4 fold)        -> gam: 465     -> rf:  475 (4 fold)
			nrow(   index[ which(index$var1 == "logP"), ]   ) +  nrow(   index[ which(index$var2 == "logP"), ]   ) + nrow(   index[ which(index$var3 == "logP"), ]   )  +  nrow(   index[ which(index$var4 == "logP"), ]   ) # 333 (8 fold weighting top 5);         360 (4 fold weighting top 5);         364 (even weights)            -> glm: 445 (4 fold)       -> gam: 292    -> rf: 433 (4 fold)
			nrow(   index[ which(index$var1 == "N"), ]   ) +  nrow(   index[ which(index$var2 == "N"), ]   ) + nrow(   index[ which(index$var3 == "N"), ]   )  +  nrow(   index[ which(index$var4 == "N"), ]   ) # 300 (8 fold weighting top 5);         313 (4 fold weighting top 5);         321 (even weights)             -> glm: 291   (4 fold)   -> gam: 404     -> rf: 199 (4 fold)
			nrow(   index[ which(index$var1 == "logN"), ]   ) +  nrow(   index[ which(index$var2 == "logN"), ]   ) + nrow(   index[ which(index$var3 == "logN"), ]   )  +  nrow(   index[ which(index$var4 == "logN"), ]   ) # 215 (8 fold weighting top 5);         217 (4 fold weighting top 5);         227 (even weights)             -> glm: 255 (4 fold)     -> gam: 228     -> rf: 237 (4 fold)
			nrow(   index[ which(index$var1 == "Sal"), ]   ) +  nrow(   index[ which(index$var2 == "Sal"), ]   ) + nrow(   index[ which(index$var3 == "Sal"), ]   )  +  nrow(   index[ which(index$var4 == "Sal"), ]   ) #      761 (8 fold weighting top 5);         752 (4 fold weighting top 5);         751 (even weights)           -> glm: 825   (4 fold)    -> gam: 595     -> rf: 902 (4 fold)


			# Save (the filename is adjusted automatically)
			# Save frames that contain stats on variables selected as .RData
			nms_set <- gsub("Analysis_ranking_Species_minobs_15_", "", vec_fname[d])
			nms_set <- gsub(".csv", "", nms_set)
			# Save index of selected variables as .csv
			fln <- paste0(output.dir, v, "/Species_minobs_15_index_5_members_a_4_predictors_(", proj.type,")_", nms_set, ".csv")
			write.csv(index, file = fln)

		# }  # end of loop proj.type
	} # end of loop vec_fname
} # end of loop vec_folders



###==============================================================
### END
###==============================================================







### =========================================================================
### B. Sites - VIF based. Select variables for the first best model across species
### =========================================================================
for( v in seq_along(vec_folders)){

	# Get filenames
	vec_fname <- list.files( paste0(input.dir, vec_folders[v]) )

	# Loop through
	for(d in seq_along(vec_fname)){

		# Print progress 
		print(paste0("Working on background selection strategy ", vec_fname[d], "."))

		for(p in proj.type){

			# Load data
			dat.ran <- read.csv( paste0(input.dir, vec_folders[v], "/", vec_fname[d]) )

			# Set seed (for reproducibility)
			set.seed(67)

			# LOOP ACROSS TAXA
			lisi.vars <- list() # the variables selected per species
			lisi.frames <- list() # for the statistics
			# s <- 22 # Trichodesmium
			for( s  in  1:length(na.omit(unique(dat.ran$taxon))) ){

				# LOOP ACROSS COMPOSED MEMBERS (VARIABLE COMBINATIONS): select five times 4 variables (without putting variables back)
				#
				# Preparatory: create a chart from which to draw variables (ranking uniquely GLM, GAM, RF based; or the mean ranking of the variables based on GLM, GAM, RF, equally considered)
				if(p == "glm") {dat.spp.raw <- dat.ran[          which(dat.ran$taxon == na.omit(unique(dat.ran$taxon))[s] )[2],   3:18 ]} # The second position is site glm
				if(p == "gam") {dat.spp.raw <- dat.ran[          which(dat.ran$taxon == na.omit(unique(dat.ran$taxon))[s] )[4],   3:18 ]} # The fourth position is site gam
				if(p == "rf") {dat.spp.raw <- dat.ran[          which(dat.ran$taxon == na.omit(unique(dat.ran$taxon))[s] )[6],   3:18 ]} # The sixth position is site rf
				if(p == "ensemble") {dat.spp.raw <- dat.ran[           max(which(dat.ran$taxon == na.omit(unique(dat.ran$taxon))[s] ))+1,   3:18 ]} # The ensemble mean rank is rank seven (the "max plus one")
					#
				dat.spp <- c(  apply(   dat.spp.raw ,   2, function(x){as.numeric(x)  }   )  )
				df.var.raw <- data.frame("var" = names(dat.spp), "rank" = as.numeric(as.character(dat.spp)))
				df.var <- df.var.raw[order(df.var.raw$rank), ]
				rownames(df.var) <- 1:nrow(df.var)
				df.var$check.vector <- 0
				df.var$weight <- 1
				head(df.var)
				# 	var rank check.vector weight
				# 1  logN    1            0      1
				# 2     T    2            0      1
				# 3     N    3            0      1
				# 4  logP    4            0      1
				# 5     P    5            0      1
				# 6 dN_dt    6            0      1

				# Weighting of variables used for the initial selection (in case both log and non log versions of variables, and MLD1 and MLD2, is contained among the top ten ranked variables); thus temperature and physical variables count equal as nutrient variables and log versions of nutrients initially
				if(  sum(c("MLD1", "MLD2") %in% df.var[1:10, ]$var) == 2 ) {   df.var[ which(df.var$var == "MLD1") ,   "weight"] <- 0.5;  df.var[ which(df.var$var == "MLD2") ,   "weight"] <- 0.5         }
				if(  sum(c("MLD1", "logMLD1") %in% df.var[1:10, ]$var) == 2 ) {   df.var[ which(df.var$var == "MLD1") ,   "weight"] <- 0.5;  df.var[ which(df.var$var == "logMLD1") ,   "weight"] <- 0.5         }
				if(  sum(c("MLD2", "logMLD2") %in% df.var[1:10, ]$var) == 2 ) {   df.var[ which(df.var$var == "MLD2") ,   "weight"] <- 0.5;  df.var[ which(df.var$var == "logMLD2") ,   "weight"] <- 0.5         }
				if(  sum(c("MLD2", "logMLD1") %in% df.var[1:10, ]$var) == 2 ) {   df.var[ which(df.var$var == "MLD2") ,   "weight"] <- 0.5;  df.var[ which(df.var$var == "logMLD1") ,   "weight"] <- 0.5         }
				if(  sum(c("logMLD2", "logMLD1") %in% df.var[1:10, ]$var) == 2 ) {   df.var[ which(df.var$var == "logMLD2") ,   "weight"] <- 0.5;  df.var[ which(df.var$var == "logMLD1") ,   "weight"] <- 0.5         }
				if(  sum(c("MLPAR1", "MLPAR2") %in% df.var[1:10, ]$var) == 2 ) {   df.var[ which(df.var$var == "MLPAR1") ,   "weight"] <- 0.5;  df.var[ which(df.var$var == "MLPAR2") ,   "weight"] <- 0.5         }
				if(  sum(c("P", "logP") %in% df.var[1:10, ]$var) == 2 ) {   df.var[ which(df.var$var == "P") ,   "weight"] <- 0.5;  df.var[ which(df.var$var == "logP") ,   "weight"] <- 0.5         }
				if(  sum(c("N", "logN") %in% df.var[1:10, ]$var) == 2 ) {   df.var[ which(df.var$var == "N") ,   "weight"] <- 0.5;  df.var[ which(df.var$var == "logN") ,   "weight"] <- 0.5         }
				if(  sum(c("Si", "logSi") %in% df.var[1:10, ]$var) == 2 ) {   df.var[ which(df.var$var == "Si") ,   "weight"] <- 0.5;  df.var[ which(df.var$var == "logSi") ,   "weight"] <- 0.5         }
				if(  sum(c("Chl", "logChl") %in% df.var[1:10, ]$var) == 2 ) {   df.var[ which(df.var$var == "Chl") ,   "weight"] <- 0.5;  df.var[ which(df.var$var == "logChl") ,   "weight"] <- 0.5         }

				# Increase the weight of the top five predictors in order to consider them preferentially for initial selections of members
				#
				# df.var[ 1:5, "weight"] <- df.var [ 1:5, "weight"]*4 # Disabled for letter as differences are not huge

				lisi.members <- list()
				for(u in 1:5){ 	# We make 5 sets a 4 predictors

					# Prepare vector for member
					my.choice <- vector()

					# SELECT FIRST VARIABLE AT RANDOM (from the available ones); store in my.choice, with equal weights to log transformed plus normal nutrients (or mld, and chl) relative to other variables
					my.potential.vars <- df.var[ which(df.var[1:10,]$check.vector < 2), ] # ** CRITICAL DECISION, TOP 10
					my.choice[1] <- as.character(my.potential.vars[sample(nrow(my.potential.vars), 1, prob = my.potential.vars$weight  ), "var"]) # sample() takes a random sample 
					df.var[   which(my.choice[1] == df.var$var), "check.vector" ] <- df.var[   which(my.choice[1] == df.var$var) , "check.vector"  ] + 1 # Adjust the check.vector field

					# Loop across the remaining variables:
					while(   length(my.choice)  < 4  ){  # ** CRITICAL DECISION: NUMBER OR PREDICTORS

						# The remaining positions that potentially provide one predictor to the member are
						valid.positions <- df.var[which(df.var$check.vector < 2 & df.var$rank < 11), ] # ** CRITICAL DECISION, TOP 10

						# The already chosen variable(s) are correlated with
						df.corr.relevant <- my.corr[     which(rownames(my.corr) %in% my.choice),     ] # The useful part of the data frame
						columns.with.values.greater.zero.point.seven <- which(apply(df.corr.relevant, 2, max, na.rm = TRUE) > 0.7)
						variables.correlated <- names(columns.with.values.greater.zero.point.seven) # These variables are correlated with the already chosen variable

						# The remaining positions that potentially provide one predictor to the memer are hence
						if (   sum(variables.correlated %in% valid.positions$var) == 0 ) {valid.pos <- valid.positions } else { # Case: none of the variables correlated with already chosen variables is found in the potential variables
						valid.pos <- valid.positions[   -which(valid.positions$var %in% variables.correlated)   ,   ]} #

							# CASE A: there are positions left among the top *ten* ranked variables
							if( nrow(valid.pos) > 0 ) {
								my.choice[length(my.choice) + 1] <- as.character(valid.pos[sample(nrow(valid.pos), 1), "var"]) # select the variable
								df.var[   which(my.choice[length(my.choice)] == df.var$var), "check.vector" ] <- df.var[   which(my.choice[length(my.choice)] == df.var$var) , "check.vector"  ]+1 # Adjust the check.vector field
								}

							# CASE B: no positions are left : loop across the lower ranked variables
							if( nrow(valid.pos) == 0 ) {
								# We start with the variable that is ranked 11 ** CRITICAL DECISION, TOP 10
								i = 11
								while(  length(my.choice) < 4  ){
									var.potential <- as.character(df.var[i, 1])
									print(i)

									if ( df.var[i, "check.vector"] > 1 ){ i = i + 1; next} # Each variable shall not be selected more than twice
										var.sel.pot <- c(my.choice, var.potential)  # test of correlation: the selected variables plus the new potential variable
										corr.subset <-  my.corr  [    which(   rownames(my.corr) %in% var.sel.pot   )    ,    which(   colnames(my.corr) %in% var.sel.pot   )  ]
										corr.subset [cbind(1:dim(corr.subset)[2], 1:dim(corr.subset)[2])]	<- NA # Replace diagonal elements by NA (but not all "1")
									if(   any ( corr.subset  > 0.7, na.rm = TRUE   ) == TRUE ) {i = i + 1; next} else {
										i = i + 1; my.choice[length(my.choice) + 1] <- var.potential
										df.var[   which(my.choice[length(my.choice)] == df.var$var), "check.vector" ] <- df.var[   which(my.choice[length(my.choice)] == df.var$var) , "check.vector"  ] + 1 # Adjust the check.vector field
									}
								}
							}
							# close across remaining variables
					}
					# Store member in list
					lisi.members[[u]] <- my.choice

					# END LOOP ACROSS COMPOSED MEMBERS (VARIABLE COMBINATIONS)
					}

				# CREATE FRAME
				df.members <- data.frame("taxon" =  rep(as.character(na.omit(unique(dat.ran$taxon))[s]), 5) ,
				"group" =  rep( as.character(dat.ran[   max(which(dat.ran$taxon == na.omit(unique(dat.ran$taxon))[s] ))    ,   ]$group)    , 5) ,
				"var1" = unlist(lapply (lisi.members, '[', 1)),
				"var2" = unlist(lapply (lisi.members, '[', 2)),
				"var3" = unlist(lapply (lisi.members, '[', 3)),
				"var4" = unlist(lapply (lisi.members, '[', 4)))

				# END LOOP ACROSS TAXA
				lisi.vars[[length(lisi.vars)+1]] <- df.members
				lisi.frames[[length(lisi.frames)+1]] <- df.var
				lisi.frames$"Emiliania huxleyi"

				print(s)
				}

			# Add names to the frames
			names(lisi.frames) <- na.omit(unique(dat.ran$taxon))
			# df.frames <- do.call(rbind, lisi.frames)

			# Assemble variables index frame
			index <- cbind(do.call(rbind, lisi.vars), "id" = 1:nrow(do.call(rbind, lisi.vars)))

			# A few stats (ensemble mean approach) * NOT UP TO DATE *
			nrow(   index[ which(index$var1 == "T"), ]   ) +  nrow(   index[ which(index$var2 == "T"), ]   ) + nrow(   index[ which(index$var3 == "T"), ]   )  +  nrow(   index[ which(index$var4 == "T"), ]   ) # 949 (8 fold weighting top 5);         932 (4 fold weighting top 5);         902 (even weights)       -> glm: 932 (4 fold)      -> gam: 911    -> rf:  899 (4 fold)
			nrow(   index[ which(index$var1 == "P"), ]   ) +  nrow(   index[ which(index$var2 == "P"), ]   ) + nrow(   index[ which(index$var3 == "P"), ]   )  +  nrow(   index[ which(index$var4 == "P"), ]   ) # 533 (8 fold weighting top 5);         534 (4 fold weighting top 5);         518 (even weights)       -> glm: 557 (4 fold)        -> gam: 465     -> rf:  475 (4 fold)
			nrow(   index[ which(index$var1 == "logP"), ]   ) +  nrow(   index[ which(index$var2 == "logP"), ]   ) + nrow(   index[ which(index$var3 == "logP"), ]   )  +  nrow(   index[ which(index$var4 == "logP"), ]   ) # 333 (8 fold weighting top 5);         360 (4 fold weighting top 5);         364 (even weights)            -> glm: 445 (4 fold)       -> gam: 292    -> rf: 433 (4 fold)
			nrow(   index[ which(index$var1 == "N"), ]   ) +  nrow(   index[ which(index$var2 == "N"), ]   ) + nrow(   index[ which(index$var3 == "N"), ]   )  +  nrow(   index[ which(index$var4 == "N"), ]   ) # 300 (8 fold weighting top 5);         313 (4 fold weighting top 5);         321 (even weights)             -> glm: 291   (4 fold)   -> gam: 404     -> rf: 199 (4 fold)
			nrow(   index[ which(index$var1 == "logN"), ]   ) +  nrow(   index[ which(index$var2 == "logN"), ]   ) + nrow(   index[ which(index$var3 == "logN"), ]   )  +  nrow(   index[ which(index$var4 == "logN"), ]   ) # 215 (8 fold weighting top 5);         217 (4 fold weighting top 5);         227 (even weights)             -> glm: 255 (4 fold)     -> gam: 228     -> rf: 237 (4 fold)
			nrow(   index[ which(index$var1 == "Sal"), ]   ) +  nrow(   index[ which(index$var2 == "Sal"), ]   ) + nrow(   index[ which(index$var3 == "Sal"), ]   )  +  nrow(   index[ which(index$var4 == "Sal"), ]   ) #      761 (8 fold weighting top 5);         752 (4 fold weighting top 5);         751 (even weights)           -> glm: 825   (4 fold)    -> gam: 595     -> rf: 902 (4 fold)


			# Save (the filename is adjusted automatically)
			# Save frames that contain stats on variables selected as .RData
			nms_set <- gsub("Analysis_ranking_Species_minobs_15_", "", vec_fname[d])
			nms_set <- gsub(".csv", "", nms_set)
			fln <- paste0(output.dir, vec_folders[v], "/Species_minobs_15_frames_5_members_a_4_predictors_(", p, ")_", nms_set, ".RData")
			save(lisi.frames, file = fln)
			# Save index of selected variables as .csv
			fln <- paste0(output.dir, vec_folders[v], "/Species_minobs_15_index_5_members_a_4_predictors_(", p,")_", nms_set, ".csv")
			write.csv(index, file = fln)

		}  # end of loop proj.type
	} # end of loop vec_fname
} # end of loop vec_folders





### =========================================================================
### B. Sites - VIF based. use best three predictors and choose last one randomly.
### =========================================================================
for( v in seq_along(vec_folders)){

	# Get filenames
	vec_fname <- list.files( paste0(input.dir, vec_folders[v]) )

	# Loop through
	for(d in seq_along(vec_fname)){

		# Print progress 
		print(paste0("Working on background selection strategy ", vec_fname[d], "."))

		for(p in proj.type){

			# Load data
			dat.ran <- read.csv( paste0(input.dir, vec_folders[v], "/", vec_fname[d]) )

			# Set seed (for reproducibility)
			set.seed(67)

			# LOOP ACROSS TAXA
			lisi.vars <- list() # the variables selected per species
			lisi.frames <- list() # for the statistics
			# s <- 22 # Trichodesmium
			for( s  in  1:length(na.omit(unique(dat.ran$taxon))) ){

				# LOOP ACROSS COMPOSED MEMBERS (VARIABLE COMBINATIONS): select five times 4 variables (without putting variables back)
				#
				# Preparatory: create a chart from which to draw variables (ranking uniquely GLM, GAM, RF based; or the mean ranking of the variables based on GLM, GAM, RF, equally considered)
				if(p == "glm") {dat.spp.raw <- dat.ran[          which(dat.ran$taxon == na.omit(unique(dat.ran$taxon))[s] )[2],   3:18 ]} # The second position is site glm
				if(p == "gam") {dat.spp.raw <- dat.ran[          which(dat.ran$taxon == na.omit(unique(dat.ran$taxon))[s] )[4],   3:18 ]} # The fourth position is site gam
				if(p == "rf") {dat.spp.raw <- dat.ran[          which(dat.ran$taxon == na.omit(unique(dat.ran$taxon))[s] )[6],   3:18 ]} # The sixth position is site rf
				if(p == "ensemble") {dat.spp.raw <- dat.ran[           max(which(dat.ran$taxon == na.omit(unique(dat.ran$taxon))[s] ))+1,   3:18 ]} # The ensemble mean rank is rank seven (the "max plus one")
					#
				dat.spp <- c(  apply(   dat.spp.raw ,   2, function(x){as.numeric(x)  }   )  )
				df.var.raw <- data.frame("var" = names(dat.spp), "rank" = as.numeric(as.character(dat.spp)))
				df.var <- df.var.raw[order(df.var.raw$rank), ]
				rownames(df.var) <- 1:nrow(df.var)
				df.var$check.vector <- 0
				df.var$weight <- 1
				head(df.var)
				# 	var rank check.vector weight
				# 1  logN    1            0      1
				# 2     T    2            0      1
				# 3     N    3            0      1
				# 4  logP    4            0      1
				# 5     P    5            0      1
				# 6 dN_dt    6            0      1

				# Weighting of variables used for the initial selection (in case both log and non log versions of variables, and MLD1 and MLD2, is contained among the top ten ranked variables); thus temperature and physical variables count equal as nutrient variables and log versions of nutrients initially
				if(  sum(c("MLD1", "MLD2") %in% df.var[1:10, ]$var) == 2 ) {   df.var[ which(df.var$var == "MLD1") ,   "weight"] <- 0.5;  df.var[ which(df.var$var == "MLD2") ,   "weight"] <- 0.5         }
				if(  sum(c("MLD1", "logMLD1") %in% df.var[1:10, ]$var) == 2 ) {   df.var[ which(df.var$var == "MLD1") ,   "weight"] <- 0.5;  df.var[ which(df.var$var == "logMLD1") ,   "weight"] <- 0.5         }
				if(  sum(c("MLD2", "logMLD2") %in% df.var[1:10, ]$var) == 2 ) {   df.var[ which(df.var$var == "MLD2") ,   "weight"] <- 0.5;  df.var[ which(df.var$var == "logMLD2") ,   "weight"] <- 0.5         }
				if(  sum(c("MLD2", "logMLD1") %in% df.var[1:10, ]$var) == 2 ) {   df.var[ which(df.var$var == "MLD2") ,   "weight"] <- 0.5;  df.var[ which(df.var$var == "logMLD1") ,   "weight"] <- 0.5         }
				if(  sum(c("logMLD2", "logMLD1") %in% df.var[1:10, ]$var) == 2 ) {   df.var[ which(df.var$var == "logMLD2") ,   "weight"] <- 0.5;  df.var[ which(df.var$var == "logMLD1") ,   "weight"] <- 0.5         }
				if(  sum(c("MLPAR1", "MLPAR2") %in% df.var[1:10, ]$var) == 2 ) {   df.var[ which(df.var$var == "MLPAR1") ,   "weight"] <- 0.5;  df.var[ which(df.var$var == "MLPAR2") ,   "weight"] <- 0.5         }
				if(  sum(c("P", "logP") %in% df.var[1:10, ]$var) == 2 ) {   df.var[ which(df.var$var == "P") ,   "weight"] <- 0.5;  df.var[ which(df.var$var == "logP") ,   "weight"] <- 0.5         }
				if(  sum(c("N", "logN") %in% df.var[1:10, ]$var) == 2 ) {   df.var[ which(df.var$var == "N") ,   "weight"] <- 0.5;  df.var[ which(df.var$var == "logN") ,   "weight"] <- 0.5         }
				if(  sum(c("Si", "logSi") %in% df.var[1:10, ]$var) == 2 ) {   df.var[ which(df.var$var == "Si") ,   "weight"] <- 0.5;  df.var[ which(df.var$var == "logSi") ,   "weight"] <- 0.5         }
				if(  sum(c("Chl", "logChl") %in% df.var[1:10, ]$var) == 2 ) {   df.var[ which(df.var$var == "Chl") ,   "weight"] <- 0.5;  df.var[ which(df.var$var == "logChl") ,   "weight"] <- 0.5         }

				# Increase the weight of the top five predictors in order to consider them preferentially for initial selections of members
				#
				# df.var[ 1:5, "weight"] <- df.var [ 1:5, "weight"]*4 # Disabled for letter as differences are not huge

				lisi.members <- list()
				for(u in 1:5){ 	# We make 5 sets a 4 predictors

					# Prepare vector for member
					my.choice <- vector()

					# SELECT FIRST VARIABLE AT RANDOM (from the available ones); store in my.choice, with equal weights to log transformed plus normal nutrients (or mld, and chl) relative to other variables
					my.potential.vars <- df.var[ which(df.var[1:10,]$check.vector < 2), ] # ** CRITICAL DECISION, TOP 10
					my.choice[1] <- as.character(my.potential.vars[sample(nrow(my.potential.vars), 1, prob = my.potential.vars$weight  ), "var"]) # sample() takes a random sample 
					df.var[   which(my.choice[1] == df.var$var), "check.vector" ] <- df.var[   which(my.choice[1] == df.var$var) , "check.vector"  ] + 1 # Adjust the check.vector field

					# Loop across the remaining variables:
					while(   length(my.choice)  < 4  ){  # ** CRITICAL DECISION: NUMBER OR PREDICTORS

						# The remaining positions that potentially provide one predictor to the member are
						valid.positions <- df.var[which(df.var$check.vector < 2 & df.var$rank < 11), ] # ** CRITICAL DECISION, TOP 10

						# The already chosen variable(s) are correlated with
						df.corr.relevant <- my.corr[     which(rownames(my.corr) %in% my.choice),     ] # The useful part of the data frame
						columns.with.values.greater.zero.point.seven <- which(apply(df.corr.relevant, 2, max, na.rm = TRUE) > 0.7)
						variables.correlated <- names(columns.with.values.greater.zero.point.seven) # These variables are correlated with the already chosen variable

						# The remaining positions that potentially provide one predictor to the memer are hence
						if (   sum(variables.correlated %in% valid.positions$var) == 0 ) {valid.pos <- valid.positions } else { # Case: none of the variables correlated with already chosen variables is found in the potential variables
						valid.pos <- valid.positions[   -which(valid.positions$var %in% variables.correlated)   ,   ]} #

							# CASE A: there are positions left among the top *ten* ranked variables
							if( nrow(valid.pos) > 0 ) {
								my.choice[length(my.choice) + 1] <- as.character(valid.pos[sample(nrow(valid.pos), 1), "var"]) # select the variable
								df.var[   which(my.choice[length(my.choice)] == df.var$var), "check.vector" ] <- df.var[   which(my.choice[length(my.choice)] == df.var$var) , "check.vector"  ]+1 # Adjust the check.vector field
								}

							# CASE B: no positions are left : loop across the lower ranked variables
							if( nrow(valid.pos) == 0 ) {
								# We start with the variable that is ranked 11 ** CRITICAL DECISION, TOP 10
								i = 11
								while(  length(my.choice) < 4  ){
									var.potential <- as.character(df.var[i, 1])
									print(i)

									if ( df.var[i, "check.vector"] > 1 ){ i = i + 1; next} # Each variable shall not be selected more than twice
										var.sel.pot <- c(my.choice, var.potential)  # test of correlation: the selected variables plus the new potential variable
										corr.subset <-  my.corr  [    which(   rownames(my.corr) %in% var.sel.pot   )    ,    which(   colnames(my.corr) %in% var.sel.pot   )  ]
										corr.subset [cbind(1:dim(corr.subset)[2], 1:dim(corr.subset)[2])]	<- NA # Replace diagonal elements by NA (but not all "1")
									if(   any ( corr.subset  > 0.7, na.rm = TRUE   ) == TRUE ) {i = i + 1; next} else {
										i = i + 1; my.choice[length(my.choice) + 1] <- var.potential
										df.var[   which(my.choice[length(my.choice)] == df.var$var), "check.vector" ] <- df.var[   which(my.choice[length(my.choice)] == df.var$var) , "check.vector"  ] + 1 # Adjust the check.vector field
									}
								}
							}
							# close across remaining variables
					}
					# Store member in list
					lisi.members[[u]] <- my.choice

					# END LOOP ACROSS COMPOSED MEMBERS (VARIABLE COMBINATIONS)
					}

				# CREATE FRAME
				df.members <- data.frame("taxon" =  rep(as.character(na.omit(unique(dat.ran$taxon))[s]), 5) ,
				"group" =  rep( as.character(dat.ran[   max(which(dat.ran$taxon == na.omit(unique(dat.ran$taxon))[s] ))    ,   ]$group)    , 5) ,
				"var1" = unlist(lapply (lisi.members, '[', 1)),
				"var2" = unlist(lapply (lisi.members, '[', 2)),
				"var3" = unlist(lapply (lisi.members, '[', 3)),
				"var4" = unlist(lapply (lisi.members, '[', 4)))

				# END LOOP ACROSS TAXA
				lisi.vars[[length(lisi.vars)+1]] <- df.members
				lisi.frames[[length(lisi.frames)+1]] <- df.var
				lisi.frames$"Emiliania huxleyi"

				print(s)
				}

			# Add names to the frames
			names(lisi.frames) <- na.omit(unique(dat.ran$taxon))
			# df.frames <- do.call(rbind, lisi.frames)

			# Assemble variables index frame
			index <- cbind(do.call(rbind, lisi.vars), "id" = 1:nrow(do.call(rbind, lisi.vars)))

			# A few stats (ensemble mean approach) * NOT UP TO DATE *
			nrow(   index[ which(index$var1 == "T"), ]   ) +  nrow(   index[ which(index$var2 == "T"), ]   ) + nrow(   index[ which(index$var3 == "T"), ]   )  +  nrow(   index[ which(index$var4 == "T"), ]   ) # 949 (8 fold weighting top 5);         932 (4 fold weighting top 5);         902 (even weights)       -> glm: 932 (4 fold)      -> gam: 911    -> rf:  899 (4 fold)
			nrow(   index[ which(index$var1 == "P"), ]   ) +  nrow(   index[ which(index$var2 == "P"), ]   ) + nrow(   index[ which(index$var3 == "P"), ]   )  +  nrow(   index[ which(index$var4 == "P"), ]   ) # 533 (8 fold weighting top 5);         534 (4 fold weighting top 5);         518 (even weights)       -> glm: 557 (4 fold)        -> gam: 465     -> rf:  475 (4 fold)
			nrow(   index[ which(index$var1 == "logP"), ]   ) +  nrow(   index[ which(index$var2 == "logP"), ]   ) + nrow(   index[ which(index$var3 == "logP"), ]   )  +  nrow(   index[ which(index$var4 == "logP"), ]   ) # 333 (8 fold weighting top 5);         360 (4 fold weighting top 5);         364 (even weights)            -> glm: 445 (4 fold)       -> gam: 292    -> rf: 433 (4 fold)
			nrow(   index[ which(index$var1 == "N"), ]   ) +  nrow(   index[ which(index$var2 == "N"), ]   ) + nrow(   index[ which(index$var3 == "N"), ]   )  +  nrow(   index[ which(index$var4 == "N"), ]   ) # 300 (8 fold weighting top 5);         313 (4 fold weighting top 5);         321 (even weights)             -> glm: 291   (4 fold)   -> gam: 404     -> rf: 199 (4 fold)
			nrow(   index[ which(index$var1 == "logN"), ]   ) +  nrow(   index[ which(index$var2 == "logN"), ]   ) + nrow(   index[ which(index$var3 == "logN"), ]   )  +  nrow(   index[ which(index$var4 == "logN"), ]   ) # 215 (8 fold weighting top 5);         217 (4 fold weighting top 5);         227 (even weights)             -> glm: 255 (4 fold)     -> gam: 228     -> rf: 237 (4 fold)
			nrow(   index[ which(index$var1 == "Sal"), ]   ) +  nrow(   index[ which(index$var2 == "Sal"), ]   ) + nrow(   index[ which(index$var3 == "Sal"), ]   )  +  nrow(   index[ which(index$var4 == "Sal"), ]   ) #      761 (8 fold weighting top 5);         752 (4 fold weighting top 5);         751 (even weights)           -> glm: 825   (4 fold)    -> gam: 595     -> rf: 902 (4 fold)


			# Save (the filename is adjusted automatically)
			# Save frames that contain stats on variables selected as .RData
			nms_set <- gsub("Analysis_ranking_Species_minobs_15_", "", vec_fname[d])
			nms_set <- gsub(".csv", "", nms_set)
			fln <- paste0(output.dir, vec_folders[v], "/Species_minobs_15_frames_5_members_a_4_predictors_(", p, ")_", nms_set, ".RData")
			save(lisi.frames, file = fln)
			# Save index of selected variables as .csv
			fln <- paste0(output.dir, vec_folders[v], "/Species_minobs_15_index_5_members_a_4_predictors_(", p,")_", nms_set, ".csv")
			write.csv(index, file = fln)

		}  # end of loop proj.type
	} # end of loop vec_fname
} # end of loop vec_folders