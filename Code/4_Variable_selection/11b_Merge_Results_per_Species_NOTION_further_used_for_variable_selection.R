### CREATE TABLE ON VARIABLE RANKING

# $Date: 2020-06-10
# Authors: Damiano Righetti, damiano.righetti@ethz.ch
#
# Description: Rank the variables based on GLM, GAM, RF models. Create a table for each phytoplankton species:
# 1. GLM sites - overlapping - total and group-specific background choice
# 2. GAM sites - overlapping - total and group-specific background choice 
# 3. RF sites - overlapping - total and group-specific background choice

# Notes:
# "Points approach" or "pts" (in filenames) refers to a selection approach of background data at 1° monthly resolution whereby different species are considered (this strategy was later discarded, as it introduces a richness signal into the pseudoabsences selection procedure). Example: In this approach, a monthly 1° cell containing 10 observation records from ten different species provides potentially ten pseudoabsences. 
# "Sites approach" or "sit" (in filenames) stands for selection of background data from "observational cells" -> briefly "sites" (i.e., pooling all species or taxa present within each 1° monthly resolution cell), as the basis to select pseodoabsences ("sit" in filenames). Example: In this approach, a monthly 1° cell containing 10 observation records from ten different species will just provide one potential pseudoabsence, as species are pooled. 

### =========================================================================
### Initialize system
### =========================================================================

# Set locale
setwd("/Users/rdamiano/Desktop/Data/8_Variable_Selection")
# kryo or clusters
# if (getwd() == "/home/rdamiano/") {setwd("/home/rdamiano/Data/8_Variable_Selection")}

# Input and output
output.dir <- "./Output_11"
             
### =========================================================================
### Load results
### =========================================================================

# glm
glm.sit <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_11/Taxa_Glm_1_Diazos_gridded_pres_abs_hom,tot_bg_overlapping(T_MLD1)_minobs_14.csv")
glm.sit$mod <- "glm,totbg"
glm.sit2 <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_11/Taxa_Glm_3_Diazos_gridded_pres_abs_hom,gr_bg_overlapping(T_MLD1)_minobs_14.csv")
glm.sit2$mod <- "glm,grbg"

# gam
gam.sit <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_11/Taxa_Gam_1_Diazos_gridded_pres_abs_hom,tot_bg_overlapping(T_MLD1)_minobs_14.csv")
gam.sit$mod <- "gam,totbg"
gam.sit2 <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_11/Taxa_Gam_3_Diazos_gridded_pres_abs_hom,gr_bg_overlapping(T_MLD1)_minobs_14.csv")
gam.sit2$mod <- "gam,grbg"

# rf
rf.sit <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_11/Taxa_Rf_1_Diazos_gridded_pres_abs_hom,tot_bg_overlapping(T_MLD1)_minobs_14.csv")
rf.sit$mod <- "rf,totbg"
rf.sit2 <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_11/Taxa_Rf_3_Diazos_gridded_pres_abs_hom,gr_bg_overlapping(T_MLD1)_minobs_14.csv")
rf.sit2$mod <- "rf,grbg"

### =========================================================================
### Total bg version: Create table - test metric1, rank, test metric 2, rank, test metric 3, rank, mean rank 
### =========================================================================

lisi.species <- list()
for ( i in 1:30 ) {
	
	# Create species data frame on ranks
	df.glm <- glm.sit[ which (glm.sit$taxon == na.omit(unique(glm.sit$taxon))[i] ), ]
	df.gam <- gam.sit[ which (gam.sit$taxon == na.omit(unique(gam.sit$taxon))[i] ), ]
	df.rf <- rf.sit[ which (rf.sit$taxon == na.omit(unique(rf.sit$taxon))[i] ), ]
	df.combined <- rbind(df.glm, df.gam, df.rf)

	# Essence of rank across the three algorithms
	vec2 <- c(apply(   df.combined [ c(2,4,6) , 3:(length(names(df.combined))-6)]  , 2, function(x){mean(as.numeric(x[which(x>0)] ))   }   )  ) # mean across species' rank
	
	# I add 1/1000 of the gam rank to the mean rank (thus gam decides on the final ranking for equally ranked variables)
	vec2 <- vec2 + 0.001*apply(   df.combined [ 4 , 3:(length(names(df.combined))-6)]  , 2, function(x){ as.numeric(x)   } )
	
	# Add essence to data frame	
	df.combined[ nrow(df.combined)+1, ] <- c(rep(NA,2), rank(vec2[1:(length(vec2))], ties.method=c("max")), "NA", "NA", "NA", "NA", "mean_rank", "rank_mean" ) # rank of mean rank
	
	# Order total
	lisi.species[[i]] <- cbind(df.combined[, c("x", "y")], df.combined[ , c("T", "P", "logP", "Sal", "pCO2" ,"Nstar", "dMLD1_dt", "N", "Wind.CCMP", "Si", "logChl", "logN", "logSi" ,"Sistar", "logMLD1", "logMLD2", "MLD1", "MLD2", 
	"Chl", "PAR", "MLPAR1", "MLPAR2" ,"dN_dt", "dT_dt", "dP_dt", "dSi_dt")], df.combined[, 29:34])
	
	# Progress
	print(paste(i))
	
	}

# Merge to frame
df.fin <- do.call(rbind, lisi.species)
df.fin[1:20, ]

# Save
fln <- paste0(output.dir,"/Analysis ranking/Species_minobs_14_sites_tot_bg.csv")
fln <- write.csv(df.fin,file=fln,row.names=F)

### =========================================================================
### Group bg version: Create table - test metric1, rank, test metric 2, rank, test metric 3, rank, mean rank 
### =========================================================================

lisi.species <- list()
for ( i in 1:30 ) {
	
	# Create species data frame on ranks
	df.glm <- glm.sit2[ which (glm.sit2$taxon == na.omit(unique(glm.sit2$taxon))[i] ), ]
	df.gam <- gam.sit2[ which (gam.sit2$taxon == na.omit(unique(gam.sit2$taxon))[i] ), ]
	df.rf <- rf.sit2[ which (rf.sit2$taxon == na.omit(unique(rf.sit2$taxon))[i] ), ]
	df.combined <- rbind(df.glm, df.gam, df.rf)

	# Essence of rank across the three algorithms
	vec2 <- c(apply(   df.combined [ c(2,4,6) , 3:(length(names(df.combined))-6)]  , 2, function(x){mean(as.numeric(x[which(x>0)] ))   }   )  ) # mean across species' rank
	
	# I add 1/1000 of the gam rank to the mean rank (thus gam decides on the final ranking for equally ranked variables)
	vec2 <- vec2 + 0.001*apply(   df.combined [ 4 , 3:(length(names(df.combined))-6)]  , 2, function(x){ as.numeric(x)   } )
	
	# Add essence to data frame	
	df.combined[ nrow(df.combined)+1, ] <- c(rep(NA,2), rank(vec2[1:(length(vec2))], ties.method=c("max")), "NA", "NA", "NA", "NA", "mean_rank", "rank_mean" ) # rank of mean rank
	
	# Order total
	lisi.species[[i]] <- cbind(df.combined[, c("x", "y")], df.combined[ , c("T", "P", "logP", "Sal", "pCO2" ,"Nstar", "dMLD1_dt", "N", "Wind.CCMP", "Si", "logChl", "logN", "logSi" ,"Sistar", "logMLD1", "logMLD2", "MLD1", "MLD2", 
	"Chl", "PAR", "MLPAR1", "MLPAR2" ,"dN_dt", "dT_dt", "dP_dt","dSi_dt")], df.combined[, 29:34])
	
	# Progress
	print(paste(i))
	
	}

# Merge to frame
df.fin <- do.call(rbind, lisi.species)
df.fin[1:20, ]

# Save
fln <- paste0(output.dir,"/Analysis ranking/Species_minobs_14_sites_gr_bg.csv")
fln <- write.csv(df.fin,file=fln,row.names=F)



