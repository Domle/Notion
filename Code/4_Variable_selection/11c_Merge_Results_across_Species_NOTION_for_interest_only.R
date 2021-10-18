### CREATE TABLE ON VARIABLE RANKING

# $Date: 2017-12-30
#
# Authors: Damiano Righetti, rdamiano@ethz.ch
#
# Description: Rank the variables based on GLM, GAM, RF models. Create a table for each phytoplankton group:
# 1. GLM points - overlapping 
# 2. GAM points - overlapping
# 3. RF points - overlapping
# Same for sites - overlapping

# Additional table on full species to test sensitivity of results regarding: min number of observations, stratification (MLD SST vs. N MLD vs. N Sal vs. P Wind.CCMP), background full vs. group

# Notes:
# "Points approach" or "pts" (in filenames) refers to a selection approach of background data at 1° monthly resolution whereby different species are considered (this strategy was later discarded, as it introduces a richness signal into the pseudoabsences selection procedure). Example: In this approach, a monthly 1° cell containing 10 observation records from ten different species provides potentially ten pseudoabsences. 
# "Sites approach" or "sit" (in filenames) stands for selection of background data from "observational cells" -> briefly "sites" (i.e., pooling all species or taxa present within each 1° monthly resolution cell), as the basis to select pseodoabsences ("sit" in filenames). Example: In this approach, a monthly 1° cell containing 10 observation records from ten different species will just provide one potential pseudoabsence, as species are pooled. 

### =========================================================================
### Initialize system
### =========================================================================

# Set locale
setwd("/Users/rdamiano/Desktop/Data/8_Variable_Selection")
# kryo or clusters
if (getwd() == "/home/rdamiano/") {setwd("/home/rdamiano/Data/8_Variable_Selection")}

# Input and output
input.dir <- "../7_Generate_Absences/Output_7"
output.dir <- "./Output_9"
 
 # Define vector of groups
vec.gr <- c( "bacillariophyceae", "chlorophyta","chrysophyceae", "cryptophyta", "cyanobacteria", "dinoflagellata", "euglenoidea",
             "haptophyta","raphidophyceae", "total")
             
### =========================================================================
### Load results
### =========================================================================

# glm
glm.pts <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Glm_gridded_hom,gr_bg,pts_overlapping(T_MLD1)_minobs_24.csv")
glm.pts$mod <- "glm,pts"
glm.sit <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Glm_gridded_hom,gr_bg,sit_overlapping(T_MLD1)_minobs_24.csv")
glm.sit$mod <- "glm,sit"

# gam
gam.pts <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Gam_gridded_hom,gr_bg,pts_overlapping(T_MLD1)_minobs_24.csv")
gam.pts$mod <- "gam,pts"
gam.sit <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Gam_gridded_hom,gr_bg,sit_overlapping(T_MLD1)_minobs_24.csv")
gam.sit$mod <- "gam,sit"

# rf
rf.pts <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Rf_gridded_hom,gr_bg,pts_overlapping(T_MLD1)_minobs_24.csv")
rf.pts$mod <- "rf,pts"
rf.sit <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Rf_gridded_hom,gr_bg,sit_overlapping(T_MLD1)_minobs_24.csv")
rf.sit$mod <- "rf,sit"

### =========================================================================
### Create tables - full approach
### =========================================================================


# groups

for(i in c(1, 2,3,4,5,6,8,9, 10)){
	# Models x points / sites	
	df.group <- rbind(	
	glm.pts[which(glm.pts$group==vec.gr[[i]]), ]
	,
	gam.pts[which(gam.pts$group==vec.gr[[i]]), ],
	rf.pts[which(gam.pts$group==vec.gr[[i]]), ],
	rep(NA, length(names(rf.pts))),
	glm.sit[which(glm.sit$group==vec.gr[[i]]), ],
	gam.sit[which(gam.sit$group==vec.gr[[i]]), ],
	rf.sit[which(gam.sit$group==vec.gr[[i]]), ],
	rep(NA, length(names(rf.pts)))
	)

	# Essence rank pts
	vec2 <- c(apply(   df.group[ c(4,8,12) , 3:(length(names(df.group))-6)]  , 2, function(x){mean(as.numeric(x[which(x>0)]))}))  # mean across species' rank
	df.group[ nrow(df.group)+1, ]<- c(rep(NA,2), rank(vec2[1:(length(vec2))], ties.method=c("max")), "NA", "NA", "NA", "NA", "mean_rank", "points" )

	# Essence rank sit
	vec2 <- c(apply(   df.group[ c(17,21,25) , 3:(length(names(df.group))-6)]  , 2, function(x){mean(as.numeric(x[which(x>0)]))}))  # mean across species' rank
	df.group[ nrow(df.group)+1, ]<- c(rep(NA,2), rank(vec2[1:(length(vec2))], ties.method=c("max")), "NA", "NA", "NA", "NA", "mean_rank", "sites" )
	
	# Mean of both
	vec3 <- c(apply(   df.group[ c(27, 28) , 3:(length(names(df.group))-6)]  , 2, function(x){mean(as.numeric(x[which(x>0)]))}))
	df.group[ nrow(df.group)+1, ]<- c(rep(NA,2), rank(vec3[1:(length(vec2))], ties.method=c("max")), "NA", "NA", "NA", "NA", "mean_rank", "mean" )


	# Order bacillario
	if (i ==1) {df.final <- cbind(df.group[, c("x", "y")], df.group[ , c("T", "P", "logP", "Sal", "N", "Sistar", "logChl", "dMLD1_dt", "Nstar", "logN", 
	"pCO2", "Wind.CCMP", "Si", "Chl", "logMLD1", "logMLD1", "MLD1", "MLD2" ,"logSi", "dN_dt", "PAR", "MLPAR1", "MLPAR2", "dT_dt", "dP_dt"   )], df.group[, 29:34])}

	# Order total
	df.final <- cbind(df.group[, c("x", "y")], df.group[ , c("T", "P", "logP", "Sal", "pCO2" ,"Nstar", "dMLD1_dt", "N", "Wind.CCMP", "Si", "logChl", "logN", "logSi" ,"Sistar", "logMLD1", "logMLD2", "MLD1", "MLD2", 
	"Chl", "PAR", "MLPAR1", "MLPAR2" ,"dN_dt", "dT_dt", "dP_dt")], df.group[, 29:34])

	# Save
	fln <- paste0(output.dir,"/Analysis ranking/",vec.gr[i],"_minobs_24.csv")
	fln <- write.csv(df.final,file=fln,row.names=F)
}



# euglenoidea

tt1 <- glm.pts[which(glm.pts$group==vec.gr[[7]]), ]
tt1 <- c(rep("NA",2), rank(      tt1[ 1, 3:(length(tt1)-6)]   , ties.method=c("max"))  , "NA", "NA", "NA", "NA", "mean_rank", "points" )

tt2 <- gam.pts[which(gam.pts$group==vec.gr[[7]]), ]
tt2 <- c(rep("NA",2), rank(      tt2[ 1, 3:(length(tt2)-6)]   , ties.method=c("max"))  , "NA", "NA", "NA", "NA", "mean_rank", "points" )

tt3 <- rf.pts[which(rf.pts$group==vec.gr[[7]]), ]
tt3 <- c(rep("NA",2), rank(      tt3[ 1, 3:(length(tt3)-6)]   , ties.method=c("max"))  , "NA", "NA", "NA", "NA", "mean_rank", "points" )

tt4 <- glm.sit[which(glm.sit$group==vec.gr[[7]]), ]
tt4 <- c(rep("NA",2), rank(      tt4[ 1, 3:(length(tt4)-6)]   , ties.method=c("max"))  , "NA", "NA", "NA", "NA", "mean_rank", "sites" )

tt5 <- gam.sit[which(gam.sit$group==vec.gr[[7]]), ]
tt5 <- c(rep("NA",2), rank(      tt5[ 1, 3:(length(tt5)-6)]   , ties.method=c("max"))  , "NA", "NA", "NA", "NA", "mean_rank", "sites" )

tt6 <- rf.sit[which(rf.sit$group==vec.gr[[7]]), ]
tt6 <- c(rep("NA",2), rank(      tt6[ 1, 3:(length(tt6)-6)]   , ties.method=c("max"))  , "NA", "NA", "NA", "NA", "mean_rank", "sites" )

df.final <- rbind(tt1, tt2, tt3, tt4, tt5, tt6)
# Save
fln <- paste0(output.dir,"/Analysis ranking/",vec.gr[7],"_minobs_24.csv")
fln <- write.csv(df.final,file=fln,row.names=F)


### =========================================================================
### Create tables - sensitivity test
### =========================================================================

# Starndard: 15 obs

# glm
glm.pts <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Glm_gridded_hom,gr_bg,pts_overlapping(T_MLD1)_minobs_15.csv")
glm.pts$mod <- "glm,pts"
glm.sit <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Glm_gridded_hom,gr_bg,sit_overlapping(T_MLD1)_minobs_15.csv")
glm.sit$mod <- "glm,sit"

# gam
gam.pts <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Gam_gridded_hom,gr_bg,pts_overlapping(T_MLD1)_minobs_15.csv")
gam.pts$mod <- "gam,pts"
gam.sit <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Gam_gridded_hom,gr_bg,sit_overlapping(T_MLD1)_minobs_15.csv")
gam.sit$mod <- "gam,sit"

# rf
rf.pts <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Rf_gridded_hom,gr_bg,pts_overlapping(T_MLD1)_minobs_15.csv")
rf.pts$mod <- "rf,pts"
rf.sit <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Rf_gridded_hom,gr_bg,sit_overlapping(T_MLD1)_minobs_15.csv")
rf.sit$mod <- "rf,sit"


# Models x points / sites	
	i = 10
	df.group <- rbind(	
	glm.pts[which(glm.pts$group==vec.gr[[i]]), ],
	gam.pts[which(gam.pts$group==vec.gr[[i]]), ],
	rf.pts[which(gam.pts$group==vec.gr[[i]]), ],
	rep(NA, length(names(rf.pts))),
	glm.sit[which(glm.sit$group==vec.gr[[i]]), ],
	gam.sit[which(gam.sit$group==vec.gr[[i]]), ],
	rf.sit[which(gam.sit$group==vec.gr[[i]]), ],
	rep(NA, length(names(rf.pts)))
	)

	# Essence rank pts
	vec2 <- c(apply(   df.group[ c(4,8,12) , 3:(length(names(df.group))-6)]  , 2, function(x){mean(as.numeric(x[which(x>0)]))}))  # mean across species' rank
	df.group[ nrow(df.group)+1, ]<- c(rep(NA,2), rank(vec2[1:(length(vec2))], ties.method=c("max")), "NA", "NA", "NA", "NA", "mean_rank", "points_15obs" )

	# Essence rank sit
	vec2 <- c(apply(   df.group[ c(17,21,25) , 3:(length(names(df.group))-6)]  , 2, function(x){mean(as.numeric(x[which(x>0)]))}))  # mean across species' rank
	df.group[ nrow(df.group)+1, ]<- c(rep(NA,2), rank(vec2[1:(length(vec2))], ties.method=c("max")), "NA", "NA", "NA", "NA", "mean_rank", "sites_15obs" )
	
	# Order total
	df.final <- cbind(df.group[, c("x", "y")], df.group[ , c("T", "P", "logP", "Sal", "pCO2" ,"Nstar", "dMLD1_dt","N", "Wind.CCMP", "Si", "logChl", "logN", "logSi" ,"Sistar", "logMLD1", "MLD1",
	"Chl", "PAR", "MLPAR1", "dN_dt", "dT_dt", "dP_dt")], df.group[, 29:34])

	part1 <- df.final[ c(nrow(df.final):(nrow(df.final)-1) ) , ]

# ------------


# Starndard: 24 obs, but tot bg

# glm
glm.pts <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Glm_gridded_hom,tot_bg,pts_overlapping(T_MLD1)_minobs_24.csv")
glm.pts$mod <- "glm,pts"
glm.sit <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Glm_gridded_hom,tot_bg,sit_overlapping(T_MLD1)_minobs_24.csv")
glm.sit$mod <- "glm,sit"

# gam
gam.pts <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Gam_gridded_hom,tot_bg,pts_overlapping(T_MLD1)_minobs_24.csv")
gam.pts$mod <- "gam,pts"
gam.sit <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Gam_gridded_hom,tot_bg,sit_overlapping(T_MLD1)_minobs_24.csv")
gam.sit$mod <- "gam,sit"

# rf
rf.pts <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Rf_gridded_hom,tot_bg,pts_overlapping(T_MLD1)_minobs_24.csv")
rf.pts$mod <- "rf,pts"
rf.sit <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Rf_gridded_hom,tot_bg,sit_overlapping(T_MLD1)_minobs_24.csv")
rf.sit$mod <- "rf,sit"


# Models x points / sites	
	i = 10
	df.group <- rbind(	
	glm.pts[which(glm.pts$group==vec.gr[[i]]), ],
	gam.pts[which(gam.pts$group==vec.gr[[i]]), ],
	rf.pts[which(gam.pts$group==vec.gr[[i]]), ],
	rep(NA, length(names(rf.pts))),
	glm.sit[which(glm.sit$group==vec.gr[[i]]), ],
	gam.sit[which(gam.sit$group==vec.gr[[i]]), ],
	rf.sit[which(gam.sit$group==vec.gr[[i]]), ],
	rep(NA, length(names(rf.pts)))
	)

	# Essence rank pts
	vec2 <- c(apply(   df.group[ c(4,8,12) , 3:(length(names(df.group))-6)]  , 2, function(x){mean(as.numeric(x[which(x>0)]))}))  # mean across species' rank
	df.group[ nrow(df.group)+1, ]<- c(rep(NA,2), rank(vec2[1:(length(vec2))], ties.method=c("max")), "NA", "NA", "NA", "NA", "mean_rank", "points_24obs_totbg" )

	# Essence rank sit
	vec2 <- c(apply(   df.group[ c(17,21,25) , 3:(length(names(df.group))-6)]  , 2, function(x){mean(as.numeric(x[which(x>0)]))}))  # mean across species' rank
	df.group[ nrow(df.group)+1, ]<- c(rep(NA,2), rank(vec2[1:(length(vec2))], ties.method=c("max")), "NA", "NA", "NA", "NA", "mean_rank", "sites_24obs_totbg" )
	
	# Order total
	df.final <- cbind(df.group[, c("x", "y")], df.group[ , c("T", "P", "logP", "Sal", "pCO2" ,"Nstar", "dMLD1_dt","N", "Wind.CCMP", "Si", "logChl", "logN", "logSi" ,"Sistar", "logMLD1", "MLD1",
	"Chl", "PAR", "MLPAR1", "dN_dt", "dT_dt", "dP_dt")], df.group[, 29:34])

	part2 <- df.final[ c(nrow(df.final):(nrow(df.final)-1) ) , ]


# ------------


# Starndard: 50 obs

# glm
glm.pts <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Glm_gridded_hom,gr_bg,pts_overlapping(T_MLD1)_minobs_50.csv")
glm.pts$mod <- "glm,pts"
glm.sit <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Glm_gridded_hom,gr_bg,sit_overlapping(T_MLD1)_minobs_50.csv")
glm.sit$mod <- "glm,sit"

# gam
gam.pts <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Gam_gridded_hom,gr_bg,pts_overlapping(T_MLD1)_minobs_50.csv")
gam.pts$mod <- "gam,pts"
gam.sit <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Gam_gridded_hom,gr_bg,sit_overlapping(T_MLD1)_minobs_50.csv")
gam.sit$mod <- "gam,sit"

# rf
rf.pts <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Rf_gridded_hom,gr_bg,pts_overlapping(T_MLD1)_minobs_50.csv")
rf.pts$mod <- "rf,pts"
rf.sit <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Rf_gridded_hom,gr_bg,sit_overlapping(T_MLD1)_minobs_50.csv")
rf.sit$mod <- "rf,sit"


# Models x points / sites	
	i = 10
	df.group <- rbind(	
	glm.pts[which(glm.pts$group==vec.gr[[i]]), ],
	gam.pts[which(gam.pts$group==vec.gr[[i]]), ],
	rf.pts[which(gam.pts$group==vec.gr[[i]]), ],
	rep(NA, length(names(rf.pts))),
	glm.sit[which(glm.sit$group==vec.gr[[i]]), ],
	gam.sit[which(gam.sit$group==vec.gr[[i]]), ],
	rf.sit[which(gam.sit$group==vec.gr[[i]]), ],
	rep(NA, length(names(rf.pts)))
	)
	
	# Change to numeric
	class(df.group)
	df.group$y <- as.numeric(as.character(df.group$y))
	df.group[ , 1:(dim(df.group)[2]-6) ]  <- apply( df.group[ , 1:(dim(df.group)[2]-6) ] , 2, function(x) as.numeric(as.character(x)))

	# Essence rank pts
	vec2 <- c(apply(   df.group[ c(4,8,12) , 3:(length(names(df.group))-6)]  , 2, function(x){mean(as.numeric(x[which(x>0)]))}))  # mean across species' rank
	df.group[ nrow(df.group)+1, ]<- c(rep(NA,2), rank(vec2[1:(length(vec2))], ties.method=c("max")), "NA", "NA", "NA", "NA", "mean_rank", "points_50obs" )

	# Essence rank sit
	vec2 <- c(apply(   df.group[ c(17,21,25) , 3:(length(names(df.group))-6)]  , 2, function(x){mean(as.numeric(x[which(x>0)]))}))  # mean across species' rank
	df.group[ nrow(df.group)+1, ]<- c(rep(NA,2), rank(vec2[1:(length(vec2))], ties.method=c("max")), "NA", "NA", "NA", "NA", "mean_rank", "sites_50obs" )
	
	# Order total
	df.final <- cbind(df.group[, c("x", "y")], df.group[ , c("T", "P", "logP", "Sal", "pCO2" ,"Nstar", "dMLD1_dt","N", "Wind.CCMP", "Si", "logChl", "logN", "logSi" ,"Sistar", "logMLD1", "MLD1",
	"Chl", "PAR", "MLPAR1", "dN_dt", "dT_dt", "dP_dt")], df.group[, 29:34])

	part3 <- df.final[ c(nrow(df.final):(nrow(df.final)-1) ) , ]


# ------------


# N - MLD1:

# glm
glm.pts <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Glm_gridded_hom,gr_bg,pts_overlapping(N_MLD1)_minobs_24.csv")
glm.pts$mod <- "glm,pts"
glm.sit <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Glm_gridded_hom,gr_bg,sit_overlapping(N_MLD1)_minobs_24.csv")
glm.sit$mod <- "glm,sit"

# gam
gam.pts <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Gam_gridded_hom,gr_bg,pts_overlapping(N_MLD1)_minobs_24.csv")
gam.pts$mod <- "gam,pts"
gam.sit <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Gam_gridded_hom,gr_bg,sit_overlapping(N_MLD1)_minobs_24.csv")
gam.sit$mod <- "gam,sit"

# rf
rf.pts <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Rf_gridded_hom,gr_bg,pts_overlapping(N_MLD1)_minobs_24.csv")
rf.pts$mod <- "rf,pts"
rf.sit <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Rf_gridded_hom,gr_bg,sit_overlapping(N_MLD1)_minobs_24.csv")
rf.sit$mod <- "rf,sit"


# Models x points / sites	
	i = 10
	df.group <- rbind(	
	glm.pts[which(glm.pts$group==vec.gr[[i]]), ],
	gam.pts[which(gam.pts$group==vec.gr[[i]]), ],
	rf.pts[which(gam.pts$group==vec.gr[[i]]), ],
	rep(NA, length(names(rf.pts))),
	glm.sit[which(glm.sit$group==vec.gr[[i]]), ],
	gam.sit[which(gam.sit$group==vec.gr[[i]]), ],
	rf.sit[which(gam.sit$group==vec.gr[[i]]), ],
	rep(NA, length(names(rf.pts)))
	)

	# Change to numeric
	class(df.group)
	df.group$y <- as.numeric(as.character(df.group$y))
	df.group[ , 1:(dim(df.group)[2]-6) ]  <- apply( df.group[ , 1:(dim(df.group)[2]-6) ] , 2, function(x) as.numeric(as.character(x)))

	# Essence rank pts
	vec2 <- c(apply(   df.group[ c(4,8,12) , 3:(length(names(df.group))-6)]  , 2, function(x){mean(as.numeric(x[which(x>0)]))}))  # mean across species' rank
	df.group[ nrow(df.group)+1, ]<- c(rep(NA,2), rank(vec2[1:(length(vec2))], ties.method=c("max")), "NA", "NA", "NA", "NA", "mean_rank", "N_MLD_points" )

	# Essence rank sit
	vec2 <- c(apply(   df.group[ c(17,21,25) , 3:(length(names(df.group))-6)]  , 2, function(x){mean(as.numeric(x[which(x>0)]))}))  # mean across species' rank
	df.group[ nrow(df.group)+1, ]<- c(rep(NA,2), rank(vec2[1:(length(vec2))], ties.method=c("max")), "NA", "NA", "NA", "NA", "mean_rank", "N_MLD_sites" )
	
	# Order total
	df.final <- cbind(df.group[, c("x", "y")], df.group[ , c("T", "P", "logP", "Sal", "pCO2" ,"Nstar", "dMLD1_dt","N", "Wind.CCMP", "Si", "logChl", "logN", "logSi" ,"Sistar", "logMLD1", "MLD1",
	"Chl", "PAR", "MLPAR1", "dN_dt", "dT_dt", "dP_dt")], df.group[, 29:34])

	part4 <- df.final[ c(nrow(df.final):(nrow(df.final)-1) ) , ]


# ------------


# P - Wind:

# glm
glm.pts <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Glm_gridded_hom,gr_bg,pts_overlapping(P_Wind.CCMP)_minobs_24.csv")
glm.pts$mod <- "glm,pts"
glm.sit <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Glm_gridded_hom,gr_bg,sit_overlapping(P_Wind.CCMP)_minobs_24.csv")
glm.sit$mod <- "glm,sit"

# gam
gam.pts <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Gam_gridded_hom,gr_bg,pts_overlapping(P_Wind.CCMP)_minobs_24.csv")
gam.pts$mod <- "gam,pts"
gam.sit <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Gam_gridded_hom,gr_bg,sit_overlapping(P_Wind.CCMP)_minobs_24.csv")
gam.sit$mod <- "gam,sit"

# rf
rf.pts <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Rf_gridded_hom,gr_bg,pts_overlapping(P_Wind.CCMP)_minobs_24.csv")
rf.pts$mod <- "rf,pts"
rf.sit <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Rf_gridded_hom,gr_bg,sit_overlapping(P_Wind.CCMP)_minobs_24.csv")
rf.sit$mod <- "rf,sit"


# Models x points / sites	
	i = 10
	df.group <- rbind(	
	glm.pts[which(glm.pts$group==vec.gr[[i]]), ],
	gam.pts[which(gam.pts$group==vec.gr[[i]]), ],
	rf.pts[which(gam.pts$group==vec.gr[[i]]), ],
	rep(NA, length(names(rf.pts))),
	glm.sit[which(glm.sit$group==vec.gr[[i]]), ],
	gam.sit[which(gam.sit$group==vec.gr[[i]]), ],
	rf.sit[which(gam.sit$group==vec.gr[[i]]), ],
	rep(NA, length(names(rf.pts)))
	)

	# Essence rank pts
	vec2 <- c(apply(   df.group[ c(4,8,12) , 3:(length(names(df.group))-6)]  , 2, function(x){mean(as.numeric(x[which(x>0)]))}))  # mean across species' rank
	df.group[ nrow(df.group)+1, ]<- c(rep(NA,2), rank(vec2[1:(length(vec2))], ties.method=c("max")), "NA", "NA", "NA", "NA", "mean_rank", "P_Wind_points" )

	# Essence rank sit
	vec2 <- c(apply(   df.group[ c(17,21,25) , 3:(length(names(df.group))-6)]  , 2, function(x){mean(as.numeric(x[which(x>0)]))}))  # mean across species' rank
	df.group[ nrow(df.group)+1, ]<- c(rep(NA,2), rank(vec2[1:(length(vec2))], ties.method=c("max")), "NA", "NA", "NA", "NA", "mean_rank", "P_Wind_sites" )
	
	# Order total
	df.final <- cbind(df.group[, c("x", "y")], df.group[ , c("T", "P", "logP", "Sal", "pCO2" ,"Nstar", "dMLD1_dt","N", "Wind.CCMP", "Si", "logChl", "logN", "logSi" ,"Sistar", "logMLD1", "MLD1",
	"Chl", "PAR", "MLPAR1", "dN_dt", "dT_dt", "dP_dt")], df.group[, 29:34])

	part5 <- df.final[ c(nrow(df.final):(nrow(df.final)-1) ) , ]


# Starndard: 15 obs, but gridded to 300 km minimum distance

# glm
glm.pts <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Glm_mt_300_hom,gr_bg,pts_overlapping(T_MLD1)_minobs_15.csv")
glm.pts$mod <- "glm,pts"
glm.sit <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Glm_mt_300_hom,gr_bg,sit_overlapping(T_MLD1)_minobs_15.csv")
glm.sit$mod <- "glm,sit"

# gam
gam.pts <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Gam_mt_300_hom,gr_bg,pts_overlapping(T_MLD1)_minobs_15.csv")
gam.pts$mod <- "gam,pts"
gam.sit <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Gam_mt_300_hom,gr_bg,sit_overlapping(T_MLD1)_minobs_15.csv")
gam.sit$mod <- "gam,sit"

# rf
rf.pts <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Rf_mt_300_hom,gr_bg,pts_overlapping(T_MLD1)_minobs_15.csv")
rf.pts$mod <- "rf,pts"
rf.sit <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Rf_mt_300_hom,gr_bg,sit_overlapping(T_MLD1)_minobs_15.csv")
rf.sit$mod <- "rf,sit"


# Models x points / sites	
	i = 10
	df.group <- rbind(	
	glm.pts[which(glm.pts$group==vec.gr[[i]]), ],
	gam.pts[which(gam.pts$group==vec.gr[[i]]), ],
	rf.pts[which(gam.pts$group==vec.gr[[i]]), ],
	rep(NA, length(names(rf.pts))),
	glm.sit[which(glm.sit$group==vec.gr[[i]]), ],
	gam.sit[which(gam.sit$group==vec.gr[[i]]), ],
	rf.sit[which(gam.sit$group==vec.gr[[i]]), ],
	rep(NA, length(names(rf.pts)))
	)

	# Change to numeric
	class(df.group)
	df.group$y <- as.numeric(as.character(df.group$y))
	df.group[ , 1:(dim(df.group)[2]-6) ]  <- apply( df.group[ , 1:(dim(df.group)[2]-6) ] , 2, function(x) as.numeric(as.character(x)))

	# Essence rank pts
	vec2 <- c(apply(   df.group[ c(4,8,12) , 3:(length(names(df.group))-6)]  , 2, function(x){mean(as.numeric(x[which(x>0)]))}))  # mean across species' rank
	df.group[ nrow(df.group)+1, ]<- c(rep(NA,2), rank(vec2[1:(length(vec2))], ties.method=c("max")), "NA", "NA", "NA", "NA", "mean_rank", "15obs_300_points" )

	# Essence rank sit
	vec2 <- c(apply(   df.group[ c(17,21,25) , 3:(length(names(df.group))-6)]  , 2, function(x){mean(as.numeric(x[which(x>0)]))}))  # mean across species' rank
	df.group[ nrow(df.group)+1, ]<- c(rep(NA,2), rank(vec2[1:(length(vec2))], ties.method=c("max")), "NA", "NA", "NA", "NA", "mean_rank", "15obs_300_sites" )
	
	# Order total
	df.final <- cbind(df.group[, c("x", "y")], df.group[ , c("T", "P", "logP", "Sal", "pCO2" ,"Nstar", "dMLD1_dt","N", "Wind.CCMP", "Si", "logChl", "logN", "logSi" ,"Sistar", "logMLD1", "MLD1",
	"Chl", "PAR", "MLPAR1", "dN_dt", "dT_dt", "dP_dt")], df.group[, 29:34])

	part6 <- df.final[ c(nrow(df.final):(nrow(df.final)-1) ) , ]



# Starndard: 15 obs, but gridded to 600 km minimum distance

# glm
glm.pts <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Glm_mt_600_hom,gr_bg,pts_overlapping(T_MLD1)_minobs_15.csv")
glm.pts$mod <- "glm,pts"
glm.sit <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Glm_mt_600_hom,gr_bg,sit_overlapping(T_MLD1)_minobs_15.csv")
glm.sit$mod <- "glm,sit"

# gam
gam.pts <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Gam_mt_600_hom,gr_bg,pts_overlapping(T_MLD1)_minobs_15.csv")
gam.pts$mod <- "gam,pts"
gam.sit <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Gam_mt_600_hom,gr_bg,sit_overlapping(T_MLD1)_minobs_15.csv")
gam.sit$mod <- "gam,sit"

# rf
rf.pts <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Rf_mt_600_hom,gr_bg,pts_overlapping(T_MLD1)_minobs_15.csv")
rf.pts$mod <- "rf,pts"
rf.sit <- read.csv("~/Desktop/Data/8_Variable_Selection/Output_9/Groups_Rf_mt_600_hom,gr_bg,sit_overlapping(T_MLD1)_minobs_15.csv")
rf.sit$mod <- "rf,sit"


# Models x points / sites	
	i = 10
	df.group <- rbind(	
	glm.pts[which(glm.pts$group==vec.gr[[i]]), ],
	gam.pts[which(gam.pts$group==vec.gr[[i]]), ],
	rf.pts[which(gam.pts$group==vec.gr[[i]]), ],
	rep(NA, length(names(rf.pts))),
	glm.sit[which(glm.sit$group==vec.gr[[i]]), ],
	gam.sit[which(gam.sit$group==vec.gr[[i]]), ],
	rf.sit[which(gam.sit$group==vec.gr[[i]]), ],
	rep(NA, length(names(rf.pts)))
	)

	# Change to numeric
	class(df.group)
	df.group$y <- as.numeric(as.character(df.group$y))
	df.group[ , 1:(dim(df.group)[2]-6) ]  <- apply( df.group[ , 1:(dim(df.group)[2]-6) ] , 2, function(x) as.numeric(as.character(x)))

	# Essence rank pts
	vec2 <- c(apply(   df.group[ c(4,8,12) , 3:(length(names(df.group))-6)]  , 2, function(x){mean(as.numeric(x[which(x>0)]))}))  # mean across species' rank
	df.group[ nrow(df.group)+1, ]<- c(rep(NA,2), rank(vec2[1:(length(vec2))], ties.method=c("max")), "NA", "NA", "NA", "NA", "mean_rank", "15obs_600_points" )

	# Essence rank sit
	vec2 <- c(apply(   df.group[ c(17,21,25) , 3:(length(names(df.group))-6)]  , 2, function(x){mean(as.numeric(x[which(x>0)]))}))  # mean across species' rank
	df.group[ nrow(df.group)+1, ]<- c(rep(NA,2), rank(vec2[1:(length(vec2))], ties.method=c("max")), "NA", "NA", "NA", "NA", "mean_rank", "15obs_600_sites" )
	
	# Order total
	df.final <- cbind(df.group[, c("x", "y")], df.group[ , c("T", "P", "logP", "Sal", "pCO2" ,"Nstar", "dMLD1_dt","N", "Wind.CCMP", "Si", "logChl", "logN", "logSi" ,"Sistar", "logMLD1", "MLD1",
	"Chl", "PAR", "MLPAR1", "dN_dt", "dT_dt", "dP_dt")], df.group[, 29:34])

	part7 <- df.final[ c(nrow(df.final):(nrow(df.final)-1) ) , ]


### Save

ttt <- rbind(part1, part2, part3, part4, part5, part6, part7)

fln <- paste0(output.dir,"/The_sensitivity_test.csv")
write.csv(ttt, file = fln)

# --- end ---