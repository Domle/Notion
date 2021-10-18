### CREATE TABLE ON VARIABLE RANKING

# $Date: 2020-10-23
#
# Authors: Damiano Righetti, rdamiano@ethz.ch
#
# Description: Rank the variables based on GLM, GAM, RF models. Create a table for each phytoplankton group:
# 1. GLM points - overlapping 
# 2. GAM points - overlapping
# 3. RF points - overlapping
# Same for sites - overlapping

# Notes:
# "Points approach" or "pts" (in filenames) refers to a selection approach of background data at 1° monthly resolution whereby different species are considered (this strategy was later discarded, as it introduces a richness signal into the pseudoabsences selection procedure). Example: In this approach, a monthly 1° cell containing 10 observation records from ten different species provides potentially ten pseudoabsences. 
# "Sites approach" or "sit" (in filenames) stands for selection of background data from "observational cells" -> briefly "sites" (i.e., pooling all species or taxa present within each 1° monthly resolution cell), as the basis to select pseodoabsences ("sit" in filenames). Example: In this approach, a monthly 1° cell containing 10 observation records from ten different species will just provide one potential pseudoabsence, as species are pooled. 
# The difference between total bg and group bg seems not evident. recheck sourcing of results.

### =========================================================================
### Initialize system
### =========================================================================

# Set locale
setwd("/Users/rdamiano/Desktop/Data/8_Variable_Selection")
# kryo or clusters
if (getwd() == "/home/rdamiano/") {setwd("/home/rdamiano/Data/8_Variable_Selection")}

# Input and output
input.dir <- "../7_Generate_Absences/Output_7"
output.dir <- "./Output_11"
 
# Define vector of groups
index <- read.csv(paste0("./Output_12_index_NOTION/Species_minobs_14_index_5_members_a_4_predictors_(ensemble)_tot_bg_NOTION.csv"), header =T, sep = ",") # member models that were fitted
as.character(unique(index$taxon))
vec.total <- c(1:30) # Positions 1:30 are all the taxa, ** ADJUST, IN CASE TAXA NUMBER CHANGES
vec.diazos <- c(1:5, 11:15, 20:30) # Positions 6 to 10 are Chaetoceros, 16 to 19 are Rhizosolenia, ** ADJUST, IN CASE TAXA NUMBER CHANGES
vec.hosts <- c(7:10, 17:19) # Positions 7 to 10 are Chaetoceros species, 17 to 19 are Rhizosolenia species, ** ADJUST, IN CASE TAXA NUMBER CHANGES

vec.gr <- c("total", "diazos", "hosts")
list.vecs <- list()
list.vecs[[1]] <- vec.total
list.vecs[[2]] <- vec.diazos
list.vecs[[3]] <- vec.hosts

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
### Total bg version and group bg version: create tables - full approach
### =========================================================================

# groups
for(i in c(1, 2, 3)) {
	
	# Get data of species/taxa 
	df.group <- rbind(	
	# tot bg
	glm.sit[ which( glm.sit$taxon%in%(as.character(unique(index$taxon))[list.vecs[[i]]]) ), ],
	gam.sit[which( gam.sit$taxon%in%(as.character(unique(index$taxon))[list.vecs[[i]]]) ), ],
	rf.sit[which( rf.sit$taxon%in%(as.character(unique(index$taxon))[list.vecs[[i]]]) ), ],	
	rep(NA, length(names(rf.sit2))),
	# gr bg
	glm.sit2[ which( glm.sit2$taxon%in%(as.character(unique(index$taxon))[list.vecs[[i]]]) ), ],
	gam.sit2[which( gam.sit2$taxon%in%(as.character(unique(index$taxon))[list.vecs[[i]]]) ), ],
	rf.sit2[which( rf.sit2$taxon%in%(as.character(unique(index$taxon))[list.vecs[[i]]]) ), ],	
	rep(NA, length(names(rf.sit2)))
	)
		
	# Essence rank total_bg, sit
	vec1basic <- c(apply(   df.group[ which(df.group$mod%in%c("glm,totbg", "gam,totbg", "rf,totbg") & df.group$value == "rank") , 3:(length(names(df.group))-6)]  , 2, function(x){mean(as.numeric(x[which(x>0)]), na.rm = T)}))  # mean across species' rank across 
	vec1 <- as.character(rank(vec1basic[1:(length(vec1basic))], ties.method=c("max")))

	# Essence value gam total bg, sit
	vec2 <- as.character(round(c(apply(   df.group[ which( df.group$mod%in%c("gam,totbg") & df.group$value%in%"adj.Rsq" ) , 3:(length(names(df.group))-6)]  , 2, function(x){mean(as.numeric(x[which(x>0)]), na.rm = T)})),3))  # mean across species' rank across 

	# Essence SD of value gam total bg, sit
	vec3 <- as.character(round(c(apply(   df.group[ which( df.group$mod%in%c("gam,totbg") & df.group$value%in%"adj.Rsq" ) , 3:(length(names(df.group))-6)]  , 2, function(x){sd(as.numeric(x[which(x>0)]), na.rm = T)})),3))  # mean across species' rank across 

	# Essence rank group_bg, sit
	vec4basic <- c(apply(   df.group[ which(df.group$mod%in%c("glm,grbg", "gam,grbg", "rf,grbg") & df.group$value == "rank") , 3:(length(names(df.group))-6)]  , 2, function(x){mean(as.numeric(x[which(x>0)]), na.rm = T)}))  # mean across species' rank across  
	vec4 <- rank(vec4basic[1:(length(vec4basic))], ties.method=c("max"))

	# Essence value gam group bg, sit
	vec5 <- as.character(round(c(apply(   df.group[ which(df.group$mod%in%c("gam,grbg") & df.group$value == "adj.Rsq") , 3:(length(names(df.group))-6)]  , 2, function(x){mean(as.numeric(x[which(x>0)]), na.rm = T)})),3))  # mean across species' rank across 

	# Essence SD of value gam total bg, sit
	vec6 <- as.character(round(c(apply(   df.group[ which(df.group$mod%in%c("gam,grbg") & df.group$value == "adj.Rsq") , 3:(length(names(df.group))-6)]  , 2, function(x){sd(as.numeric(x[which(x>0)]), na.rm = T)})),3))  # mean across species' rank across 

	# Rbind
	df2 <- data.frame(rbind(vec1, vec2, vec3, vec4, vec5, vec6)) 
	df2$x <- rep("NA", nrow(df2))
	df2$y <- rep("NA", nrow(df2))
	df2$obs <- rep("NA", nrow(df2))
	df2$taxon <- rep("NA", nrow(df2))
	df2$group <- rep("NA", nrow(df2))
	df2$pa_strat <-  rep("NA", nrow(df2))
	df2$value <-  c("rank", "adj.Rsq", "SD.adj.Rsq","rank", "adj.Rsq", "SD.adj.Rsq")
	df2$mod <-  c("tot_bg", "gam_only_tot_bg", "gam_only_tot_bg", "gr_bg", "gam_only_gr_bg", "gam_only_gr_bg")
	df.final <- rbind(df.group, df2[, names(df.group)])
	
	# Order
	df.final <- cbind(df.final[, c("x", "y")], df.final[ , c("T", "N", "dN_dt", "logN", "Sistar" ,"logChl", "Chl", "P", "MLPAR2", "logP", "PAR", "MLPAR1", "pCO2" ,"dMLD1_dt", "Wind.CCMP", "Sal", "dT_dt", "dP_dt", "logMLD1", "MLD1", "Nstar" ,"Si", "logMLD2", "MLD2", "logSi", "dSi_dt")], df.final[, 30:34])

	# Progress
	print(paste(i))
	
	# Save
	fln <- paste0(output.dir,"/Analysis ranking/",vec.gr[i],"_minobs_14.csv")
	fln <- write.csv(df.final,file=fln,row.names=F)
}

