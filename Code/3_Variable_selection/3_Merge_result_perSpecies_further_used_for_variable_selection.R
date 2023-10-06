### CREATE TABLES WITH CANDIDATE VARIABLE SKILL RANKINGS PER SPECIES, INCORPORATING DATA AND ALGORITHMS UNCERTAINTIES BROADLY AND SYSTEMATICALLY

# Input files:

# Output files:

# Strategy: Biologically flexible, objective predictor variable selection per taxon or species is based upon a ranking of variables at the species level, implemented in this code.
# Variable selection at species-specific detail is cricitally important as not all species do have the same predictors as their major driving niche-factors.
# Variables are ranked based on goodness of fit to plankton presence vs. pseudoabsence data using {single variable GLM, GAM, RF models} times {total bg, group bg, cr bg choice} times {overlapping, and non-overlapping background selection choice} = 18 variable rankings, integrating all uncertainties of choices (i.e., aspirational gold standard)

# The 'gold standard' is to create a final table that ranks the variables from best (position 1) to worst (last position) for each diazotroph species, embracing all uncertainties from {6 different background selection approaches} times {3 SDM algorithm choiceds}, i.e. 18 variable rankings per species. That is the
# GLM, GAM, and RF-based rankings (using D-squared, R-squared, or out-of-bag error, respectively) * 6 background strategies, being averaged to obtain the final variable ranking per species.


# NOTE: We will remove the CIII cluster since it is a combination of various unidentified diazotrophs and not representative for a niche of a specific taxon.


# Authors: 	Dominic Eriksson
#			Environmental Physics Group, UP
#			ETH Zurich
#			Switzerland

# Script developed by Damiano Righetti, damiano.righetti@ethz.ch

# deriksson@ethz.ch, 1st of December 2022 ------------------------------------------------------------------------



### =========================================================================
### Preparation
### =========================================================================
# Clear working space
rm(list = ls())

# Set locale
setwd("/net/kryo/work/deriksson/Projects/Notion_DE/Code/5_Variable_selection")

# Directories
wd_dat <- "/home/deriksson/Projects/Notion_DE/Code/5_Variable_selection/2_Output/"
output.dir <- "/net/kryo/work/deriksson/Projects/Notion_DE/Code/5_Variable_selection/3_Output/"

# Create vectors
vec_folders <- c("MicroscopyBased_dataset", "SequenceBased_dataset", "Total_dataset")

# Loop
for(v in seq_along(vec_folders)){

	# Print progress
	print(paste0("Datatype ", v))

	nms.sets <- grep("Taxa", list.files(paste0(wd_dat, vec_folders[v])), value = TRUE)
	nms.sets <- gsub("\\_minobs_14.csv", "", nms.sets)
	nms.sets <- gsub("Taxa_Gam_", "", nms.sets)
	nms.sets <- gsub("Taxa_Glm_", "", nms.sets)
	nms.sets <- gsub("Taxa_Rf_", "", nms.sets)
	nms.sets <- unique(nms.sets)
	# Environmental parameters of interest
	vc.vars <- c(
        "x","y","T","Sal","N","P","Si",
		"MLD1", "MLD2","PAR","Chl","Wind.CCMP","pCO2",
		"MLPAR1", "MLPAR2", "Nstar", "Sistar",
		"dT_dt", "dN_dt", "dP_dt", "dSi_dt","dMLD1_dt",
		"logMLD1", "logMLD2", "logChl", "logN", "logP", "logSi")
		# I remove NPP, as well reflected in Chl


	## Load data
	vec_fname <- grep("pres_abs", list.files(paste0(wd_dat, vec_folders[v])), ignore.case = TRUE, value = TRUE)
	l <- list()
	for (i in seq_along(vec_fname)){
		print(i)
		df <- read.csv(paste0(wd_dat, vec_folders[v], "/", vec_fname[i]))
		df$mod <- vec_fname[i]
		l[[i]] <- df
	}

	###=========================================================================
	### Total bg version: Create table - test metric1, rank, test metric 2, rank, test metric 3, rank, mean rank
	### =========================================================================
	# Loop
	for( n in seq_along(nms.sets) ){

		# Print progress
		print(n)

		# Load data
		glm.sit1 <- read.csv( paste0(wd_dat, vec_folders[v], "/Taxa_Glm_", nms.sets[n], "_minobs_14.csv") )
		gam.sit1 <- read.csv( paste0(wd_dat, vec_folders[v], "/Taxa_Gam_", nms.sets[n], "_minobs_14.csv") )
		rf.sit1 <- read.csv( paste0(wd_dat, vec_folders[v], "/Taxa_Rf_", nms.sets[n], "_minobs_14.csv") )

		# Get vector of taxa
		vec_taxon <- na.omit(unique(glm.sit1$taxon))

		# Loop
		lisi.species <- list()
		for ( i in seq_along(vec_taxon) ) {

			# Create species data frame on ranks
			df.glm <- glm.sit1[ which(glm.sit1$taxon == vec_taxon[i]), ]
			df.gam <- gam.sit1[ which(gam.sit1$taxon == vec_taxon[i]), ]
			df.rf <- rf.sit1[ which(rf.sit1$taxon == vec_taxon[i]), ]
			df.combined <- rbind(df.glm, df.gam, df.rf)

			# Essence of rank across the three algorithms
			vec2 <- c(apply(   df.combined [ c(2, 4, 6) , 3:(length(names(df.combined)) - 5)]  , 2, function(x){mean(as.numeric(x[which(x > 0)] ))   }   )  ) # mean across species' rank

			# I add 1/1000 of the gam rank to the mean rank (thus gam decides on the final ranking for equally ranked variables)
			vec2 <- vec2 + 0.001 * apply(   df.combined [ 4 , 3:(length(names(df.combined)) - 5)]  , 2, function(x){ as.numeric(x)   } )

			# Add essence to data frame
			df.combined[ nrow(df.combined) + 1, ] <- c(rep(NA, 2), rank(vec2[1:(length(vec2))], ties.method = c("min")), "NA", "NA", "NA", "mean_rank", "rank_mean" ) # rank of mean rank

			# Order total
			lisi.species[[i]] <- cbind(
				df.combined[, vc.vars],
				df.combined[, (ncol(df.combined) - 4) : ncol(df.combined)])
			# Progress
			print(paste(i))

		} # end of loop vec_taxon

		# Merge to frame
		df.fin <- do.call(rbind, lisi.species)

		# Save
		fln <- paste0(output.dir, vec_folders[v], "/Analysis_ranking_Species_minobs_15_", nms.sets[n], ".csv")
		fln <- write.csv(df.fin, file = fln, row.names = FALSE)

	}
} # end of loop vec_folders


###==============================================================
### END
###==============================================================