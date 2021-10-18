### GENERATE PSEUDOABSENCES (TARGET-GROUP AND TARGET-CRUISE APPROACH, WITH NONOVERLAPPING BACKGROUND, SELECTION VIA ENVIRONMENTALLY STRATIFIED DESIGN)

# $Date: 24.09.2021
# Author: Dominic Eriksson, This script was developed by Damiano Righetti

# Description: Generate pseudoabsences using the target-group approach (Phillips et al., 2009) for each species with ≥ 14 observations.
# Absences are generated proportionally to the number of species presences within environmental strata of species in a defined TARGET GROUP (comprising also points of species <15 observations)
# The environmental stratification shall ensure that environmental space along two key variables is reasonably addressed in terms of balancing presences and pseudo-absences.
# The density surface of points/cells which pseudo-absences are selected from is defined as all SPECIES PRESENCES (monthly cells 1° lat x 1° lon) of a corresponding larger "target-group" (class, phylum or total plankton)
# Counting presences of multiple species at same location count as one "presence" or "sampling cell" (here termed "sites", preferred approach) as taking species-resolved presences ("points") induces some sort of richness effect)
# Distribution of pseudo-absences per taxon shall thereby exhibit similar bias as its presences (Phillips et al., 2009, Ecol. Appl.)
# Note: To allow for subsequent modeling, the number of presences of the focal species must not exceed the total number potential absences from the target-group.
# Note: Non-detection records are retained in the sampling-surface for background selection, but of course excluded from the focal species' presences

# Input 1: Merged data frames (prepared as list of taxa) containing gridded presences - unthinned or previously thinned on a spatial or monthly-spatial basis
# Input 2: Environmental variable stack containing data for each monthly (climatological) 1°-cell in the ocean
# Output: Data frames (list of taxa) containing presences and pseudoabsences for each taxon ≥ 15 obs (else Null element). I.e. taxa below 15 obs are replaced by NULL elements, keeping the original length or number of taxa
# Note: The standard approach for the Letter (Sci Adv, Righetti et al., 2019) was "group background", "sites overlapping", and "gam"

## Initialize the system
rm(list = ls())
require(raster); require(classInt); require(doParallel); require(dismo)
setwd("/net/kryo/work/deriksson/notion/")

#Operational
rasterOptions(todisk=TRUE,timer=T,progress="text") # Instead of saving data in memory, save on disc

### Get data, define strategies
# Input gridded plankton data (thinning strategies) (n = 2 elements)
lisi.dat.s <- list()
lisi.dat.s[[1]] <- read.csv("./3_Match_up/Output_2/Obs_gridded_diazo_pres_abs.csv")
lisi.dat.s[[2]] <- read.csv("./3_Match_up/Output_2/Obs_gridded_diazo_nifHcounts_abs.csv")
lisi.dat.s <- rep(lisi.dat.s, 2)
lapply(lisi.dat.s, nrow) # 5042; 3020

# Input ungridded plankton data (to construct target-cruise approach, including all methods - here sources - that are capable of detecting the focal taxon)
lisi.dat.raw <- list()
lisi.dat.raw[[1]] <- read.csv("./3_Match_up/Output_2/Obs_raw_diazo_pres_abs.csv")
lisi.dat.raw[[2]] <- read.csv("./3_Match_up/Output_2/Obs_raw_diazo_nifHcounts_abs.csv")
lisi.dat.raw <- rep(lisi.dat.raw, 2)
lapply(lisi.dat.raw, nrow) # 25048; 9644

# Concept: 2 background selection strategies * 2 data type strategies
nms.sets <- c(
"pres_abs,gr_bg_nonov",
"counts_abs,gr_bg_nonov",
"pres_abs,cr_bg_nonov",
"counts_abs,cr_bg_nonov"
)

# Specify a minimum observations number per taxon (skip taxa below, and replace by NULL element)
min.obs <- 14

# Test: What stratifying variable looses the least data?
dat.s <- lisi.dat.s[[1]]
length(which(is.na(dat.s$T))) # CHOOSE
length(which(is.na(dat.s$N))) # (choose)
length(which(is.na(dat.s$MLD1))) # CHOOSE
length(which(is.na(dat.s$Wind.CCMP))) # avoid
length(which(is.na(dat.s$PAR))) # CHOOSE
length(which(is.na(dat.s$MLPAR1))) # (choose)
# plot(  dat.s[which(is.na(dat.s$pCO2)), ]$x, dat.s[which(is.na(dat.s$pCO2)), ]$y ) # pCO2 excludes the Mediterranean and the Arctic sea

# Specify parameters to stratify the sampled environment: ** ADJUST IN CASE **
vec.strat <- c("T","MLD1") # corr = .... , SELECTED
# vec.strat <- c("N","MLD1")
# vec.strat <- c("N","Sal")
# vec.strat <- c ("P","Wind.CCMP")
# vec.strat <- c("T","PAR")

# Preparatory
random_num <- seq(0,1,0.01)

### Loop across data * background selection strategies: generate pseudo-absences proportional to presences within a maximum of 81 environmental strata

## Loop across data set or backgound selection strategies
for(d in 1:4){

	# Get gridded data and raw data (used to derive the cruise specific backgrounds for d == 3 and 4)
	dat.s <- lisi.dat.s[[d]]
	dat.raw <- lisi.dat.raw[[d]]

	# Prepare a dataframe of the focal taxa (for which background data shall be sampled)
	df.taxa.obs <- data.frame("taxon" = unique(as.character(lisi.dat.s[[1]]$scientificName)))
	dat.s$scientificName <- as.character(dat.s$scientificName)
	dat.splitted <- split(dat.s, f = dat.s$scientificName)
	# lapply(dat.splitted, nrow)
	df.taxa <- data.frame("taxon" = names(unlist(lapply(dat.splitted, nrow))), "obs" = as.numeric(unlist(lapply(dat.splitted, nrow))) )
	df.taxa <- df.taxa[order(df.taxa$taxon),]
	rownames(df.taxa) <- 1:nrow(df.taxa)
	df.taxa.obs$obs <- df.taxa[match(df.taxa.obs$taxon, df.taxa$taxon),][,2]
	df.taxa.obs # 49, 17, 49, 17 taxa with data

	# ---------------------------------------- PREPARE POTENTIAL BACKGROUND SURFACES FOR THE SPECIES: GROUP-SPECIFIC TARGET GROUP SITES (monthly cells), and TARGET-CRUISE BASED SITES (monthly cells) ------------------------------------
	# Prepare total target group background and species data: stratification, using total species gridded presences AND gridded non-detections ** ADJUST IN CASE **
	all_spec_id <- dat.s # all_spec_id_meaning the species level informatoin is still retained

	# Choice of two stratifying parameters defined above in vec.strat, remove NAs from background data regarding these two variables
	if (length(which(is.na(   all_spec_id[,which(names(all_spec_id)==vec.strat[1])]  )))!= 0) {all_spec_id <- all_spec_id[-which(is.na(    all_spec_id[,which(names(all_spec_id)==vec.strat[1])]    )),]} # remove NAs regarding Var 1
	if (length(which(is.na(   all_spec_id[,which(names(all_spec_id)==vec.strat[2])]  )))!= 0) {all_spec_id <- all_spec_id[-which(is.na(    all_spec_id[,which(names(all_spec_id)==vec.strat[2])]    )),]} # remove NAs regarding Var 2
	if (length(which(is.na(   all_spec_id[,which(names(all_spec_id)==vec.strat[1])]  ))) == 0 & length(which(is.na(   all_spec_id[,which(names(all_spec_id)==vec.strat[2])]   ))) == 0) {all_spec_id <- all_spec_id}

	# Specify the range of variable-values and strata into which the values fall (this will allow to drive sampling of absences proportional to overall presences in strata)
	x_envir <- all_spec_id[,which(names(all_spec_id)==vec.strat[1])]
	y_envir <- all_spec_id[,which(names(all_spec_id)==vec.strat[2])]
	# Split ranges into environmental strata (input Brun, take 3 variables?)
	breaks <- 9
	# Create matrix that divides range into 9 equal parts; with two variables we get a maximum of 81 strata
	x_breaks <- classIntervals(x_envir, breaks, style="equal")
	x_matrix <- cbind(x_breaks$brks[1:breaks],x_breaks$brks[2:(breaks + 1)],ID =1:breaks)
	y_breaks <- classIntervals(y_envir,breaks, style="equal")
	y_matrix <- cbind(y_breaks$brks[1:breaks],y_breaks$brks[2:(breaks + 1)],ID =1:breaks)
	# A. SPECIES x CELLS: define vector of length of total points of environmental variable
	x_reclass <- vector(); x_reclass <- c(1:length(x_envir))
	y_reclass <- vector(); y_reclass <- c(1:length(y_envir))
		# Allocate points from full data to one of the nine environmental strata per variable
 		for(i in 1:breaks){
  		x_reclass[which(x_envir >= x_matrix[i,1] & x_envir <= x_matrix[i,2] )] <- x_matrix[i,3]
  		y_reclass[which(y_envir >= y_matrix[i,1] & y_envir <= y_matrix[i,2] )] <- y_matrix[i,3]
  		}
	# Create ID denoting the stratum (unique combination of variables) into which each point falls in full data-frame
	all_species_ids_reclass <- data.frame(all_spec_id, x_rcls = x_reclass, y_rcls= y_reclass, xy_rcls=x_reclass+10*y_reclass)

	# Reproducibility
	set.seed(47)

	# ----------------------------------------------------------- Parallel computing across 24 species to generate pseudo-absences -----------------------------------------------------------
	n.cores <- detectCores()-4
	cl = makeCluster(n.cores, outfile="")
	registerDoParallel(cl)
	pres_abs_results <- foreach (s = 1:length(df.taxa.obs$taxon), .packages = c('dismo', 'classInt'))%dopar%{ # Note, packages must be provided specifically to the do parallel command

	# Example species
	# s = 24 # Gamma
	# s = 36 # Richelia

		# A. Total target-group background (tot bg), collapse the prepared background cells to monthly cells (exclude taxon information)
		if ( d%in%c(1,2) ) {
			# Now: EXCLUDE the presence sampling cells of the focal taxon from the target group sampling-surface, to create "non-overlapping" background. ** ADJUST IN CASE **
			collap_ids <- all_species_ids_reclass[ !duplicated(all_species_ids_reclass[ , c("x", "y", "month")]), ]
			dat_focal_taxon <- all_species_ids_reclass[which(as.character(all_species_ids_reclass$scientificName) == df.taxa.obs$taxon[s]), ] # get species data grom the target group
			pres_cells_focal_taxon <- dat_focal_taxon[which(dat_focal_taxon$occurrenceStatus == "PRESENT"), ]
			collapsed_ids <- collap_ids[ -which(  paste0(collap_ids$x,"_",collap_ids$y,"_",collap_ids$month)%in% paste0(pres_cells_focal_taxon$x,"_", pres_cells_focal_taxon$y,"_", pres_cells_focal_taxon$month) )  , ]
		}

		# B. Target cruise bg (cr bg), collapse background to monthly cells (exclude taxon level information)
			 if (nrow(all_species_ids_reclass[all_species_ids_reclass$scientificName == df.taxa.obs$taxon[s], ]) > 0){
		if ( d%in%c(3,4) ) {
			# Prepare target cruise background approach and species data: stratification, using total species gridded presences AND gridded non-detections ** ADJUST IN CASE **
			methods_underlying_focal_species <- unique(dat.raw[which(dat.raw$scientificName == as.character(df.taxa.obs$taxon[s])), ]$measurementMethod) # Specific for: "Direct biomass analysis", "Epifluorescence microscopy", "SSU_rRNA_detection", "Standard light microscopy", "Standard light or epifluorescence microscopy", "qPCR_nifH_detection"
			target_cruise_points <- dat.raw[ which(dat.raw$measurementMethod%in%methods_underlying_focal_species), ]
			all_cruise_id <- target_cruise_points[!duplicated(target_cruise_points[ , c("scientificName", "x", "y", "month")]),]
			dim(all_cruise_id) # 1 728
			# Choice of two stratifying parameters defined above in vec.strat, remove NAs from background data regarding these two variables
			if (length(which(is.na(   all_cruise_id[,which(names(all_cruise_id)==vec.strat[1])]  )))!= 0) {all_cruise_id <- all_cruise_id[-which(is.na(    all_cruise_id[,which(names(all_cruise_id)==vec.strat[1])]    )),]} # remove NAs regarding Var 1
			if (length(which(is.na(   all_cruise_id[,which(names(all_cruise_id)==vec.strat[2])]  )))!= 0) {all_cruise_id <- all_cruise_id[-which(is.na(    all_cruise_id[,which(names(all_cruise_id)==vec.strat[2])]    )),]} # remove NAs regarding Var 2
			if (length(which(is.na(   all_cruise_id[,which(names(all_cruise_id)==vec.strat[1])]  ))) == 0 & length(which(is.na(   all_cruise_id[,which(names(all_cruise_id)==vec.strat[2])]   ))) == 0) {all_cruise_id <- all_cruise_id}
			dim(all_cruise_id) # 1 721
			# Specify the range of variable-values and strata into which the values fall (this will allow to drive sampling of absences proportional to overall presences in strata)
			x_envir <- all_cruise_id[,which(names(all_cruise_id)==vec.strat[1])]
			y_envir <- all_cruise_id[,which(names(all_cruise_id)==vec.strat[2])]
			# Split ranges into environmental strata (input Brun, take 3 variables?)
			breaks <- 9
			# Create matrix that divides range into 9 equal parts; with two variables we get a maximum of 81 strata
			x_breaks <- classIntervals(x_envir, breaks, style="equal")
			x_matrix <- cbind(x_breaks$brks[1:breaks],x_breaks$brks[2:(breaks + 1)],ID =1:breaks)
			y_breaks <- classIntervals(y_envir,breaks, style="equal")
			y_matrix <- cbind(y_breaks$brks[1:breaks],y_breaks$brks[2:(breaks + 1)],ID =1:breaks)
			# A. SPECIES x CELLS: define vector of length of total points of environmental variable
			x_reclass <- vector(); x_reclass <- c(1:length(x_envir))
			y_reclass <- vector(); y_reclass <- c(1:length(y_envir))
			# Allocate points from full data to one of the nine environmental strata per variable
 				for(i in 1:breaks){
  				x_reclass[which(x_envir >= x_matrix[i,1] & x_envir <= x_matrix[i,2] )] <- x_matrix[i,3]
		  		y_reclass[which(y_envir >= y_matrix[i,1] & y_envir <= y_matrix[i,2] )] <- y_matrix[i,3]
  				}
  			# Create ID denoting the stratum (unique combination of variables) into which each point falls in full data-frame
  			all_cruise_ids_reclass <- data.frame(all_cruise_id, x_rcls = x_reclass, y_rcls= y_reclass, xy_rcls=x_reclass+10*y_reclass)
  			# Now: EXCLUDE focal taxon from target cruise sampling-surface, to create "non-overlapping" background. Thu, pseudoabsences are treated as "true absences". ** ADJUST IN CASE **
  			collap_ids <- all_cruise_ids_reclass[ !duplicated(all_cruise_ids_reclass[ , c("x", "y", "month")]), ] # Reduce background to x, y, month
			dat_focal_taxon <- all_cruise_ids_reclass[which(as.character(all_cruise_ids_reclass$scientificName) == df.taxa.obs$taxon[s]), ] # Get species data from the target cruises
			pres_cells_focal_taxon <- dat_focal_taxon[which(dat_focal_taxon$occurrenceStatus == "PRESENT"), ]
			collapsed_ids <- collap_ids[ -which(  paste0(collap_ids$x,"_",collap_ids$y,"_",collap_ids$month)%in% paste0(pres_cells_focal_taxon$x,"_", pres_cells_focal_taxon$y,"_", pres_cells_focal_taxon$month) )  , ]

  			# Add lacking columns as NAs, needed for merger with presences - data below
  			collapsed_ids$av.nifHcount <- NA
  			collapsed_ids$sd.nifHcounts <- NA
  			collapsed_ids$observations <- NA
  			# dim(collapsed_ids)
			# plot(collapsed_ids$x, collapsed_ids$y)
		}
			}

		# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		# Preparatory 2: Skip taxa below minimum number of presences (e.g. 15)
		ob <- nrow(all_species_ids_reclass[all_species_ids_reclass$scientificName == df.taxa.obs$taxon[s], ]) #
		# ob <- nrow(all_cruise_id[all_cruise_id$scientificName == df.taxa.obs$taxon[s], ]) # Same
		# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

		# Open If: observatios are sufficient
		if (ob<min.obs){
		print( paste0(s," ",as.character(all_species_ids_reclass[all_species_ids_reclass$scientificName == df.taxa.obs$taxon[s], c("scientificName")]) [1], ".. ", ob," yielded too few obs: skipped")); dat <- NULL

		}else{
  		print( paste0(s," ",as.character(all_species_ids_reclass[all_species_ids_reclass$scientificName == df.taxa.obs$taxon[s], c("scientificName")]) [1], ": ", ob," obs -------------- "))

  		# Preparatory 3: Obtain the presences of the target species (presence cells)
  		sp_table <- all_species_ids_reclass[all_species_ids_reclass$scientificName== df.taxa.obs$taxon[s], ]
  		sp_table <- sp_table[which(sp_table$occurrenceStatus == "PRESENT"), ] # ** ADJUST IN CASE**
  		dim(sp_table)

		# Preparatory 4: Extract frequencies by which points/sites of the target group fall into environmental strata. Then, derive the number of desired absences for the focal species per stratum
		# Here we include the focal in this estimate of general sampling effort
		xy_rcls_freq <- data.frame(table(collapsed_ids$xy_rcls)/length(collapsed_ids$x_rcls)) # Frequency of points or samples of total target-group-species-occurrences per stratum
		names(xy_rcls_freq)[1] <-"xy_rcls" # Give name to column
		xy_rcls_freq[,1] <- as.numeric(as.character(xy_rcls_freq[,1])) # Define reclassification frequency as numeric
  		table_abs1 <- data.frame(xy_rcls_freq, prop_abs = (dim(sp_table)[1]*10)*xy_rcls_freq$Freq) # Add background points to be sampled per stratum: 10 x more than presences, driven by general sampling effort
  		table_abs2 <- data.frame(table_abs1, prop_abs_0 = table_abs1$prop_abs - floor(table_abs1$prop_abs)) # To round desired absences to integer: adds column difference between smaller closest integer and desired number
  		table_abs3 <- data.frame(table_abs2, prob = sample(random_num, dim(table_abs1)[1], replace=T)) # Adds column with random number between 0 and 1 (with steps 0.01)
  		table_abs4 <- data.frame(table_abs3, absences =1) # Adds column with "1" s
  		table_abs4$absences[which(table_abs4$prop_abs_0 > table_abs4$prob)] <- ceiling(table_abs4$prop_abs[which(table_abs4$prop_abs_0 > table_abs4$prob)]) # Round absences up for random subset
  		table_abs4$absences[which(table_abs4$prop_abs_0 < table_abs4$prob)] <- floor(table_abs4$prop_abs[which(table_abs4$prop_abs_0 < table_abs4$prob)]) # Round absences down for random subset
  		absence_groups <- table_abs4[table_abs4$absences>0,] # Skip strata without presences

	  	# Preparatory 5: Select backround data, here including the points/sites of the focal species ('overlapping background') in case of the general target-group approach,
	  	# i.e. we treat pseudo-absences as a general background, not as confirmed absences
	  	# we also include the species in case of the more specific target-cruise approach,
	  	# i.e., here we treat pseudo-absences as a general background, not as confirmed absences, but other options may be tested as well (i.e., excluding the species)
  		absence_table <- collapsed_ids[ , colnames(sp_table)] # needed to select exactly same columns (from the cruise specific background), as else it wont match below

  		# 1. Select background points for target species in each stratum proportionally to the density of samples/points in the background surface (in a random way)
  		absences_per_group <- list()
  		for (i in 1:nrow(absence_groups)){
  			# 1. Select available absences within stratum in question
  			group_absence_table <- absence_table[absence_table$xy_rcls == absence_groups[i,"xy_rcls"],]
	  		# 2. Define number of absences that can be drawn maximally (if the number is smaller than desired, define maximum possible number of absences (limited by points per stratum):
  			absence_num <- ifelse(
  			absence_groups[i,"absences"] > dim(absence_table[absence_table$xy_rcls == absence_groups[i,"xy_rcls"],])[1] # Test: is the number of desired background points bigger than the available background points?
  			,
  			dim(absence_table[absence_table$xy_rcls == absence_groups[i,"xy_rcls"],])[1] # if TRUE the potential points are insufficient - however, save the number of available points as absence_num
  			,
  			absence_groups[i,"absences"] # ELSE: save the number of desired background points as absence_num
  			)
  			# 3. Randomly sample the background points from the table containing all possible absences for the stratum in question
  			sampled_group_absence_table <- group_absence_table[sample(1:dim(group_absence_table)[1], absence_num),]
  			absences_per_group[[i]] <- sampled_group_absence_table # we store the product for the stratum or group
  		}

  		# 2. Merge presences, denoted with obs = 1 and absences with obs  = 0 (i.e. rowbound absence strata)
  		sp_table_4_models <- rbind(data.frame(sp_table, obs=1), data.frame(do.call("rbind", absences_per_group), obs=0))

      df1 <- data.frame(sp_table, obs = 1)
      df2 <- data.frame(do.call("rbind", absences_per_group), obs = 0)

  		dat <- cbind(sp_table_4_models, weights = 1) # Create column with weights = 1; weights are associated with presences and absences for modelling
  		prec_abs_stats <- matrix(table(dat$obs)) # Number or absences versus number of presences (does not cause trouble if we do not get 10* more absences)
  		weights_abs <- prec_abs_stats[2,1]/prec_abs_stats[1,1] # Ratio of presence count divided by absences count
  		dat$weights[dat$obs == 0] <- weights_abs # For observation that are absences we add the weight
  		row.names(dat)<-seq(1:nrow(dat)) # # Add row ID (order is important for later cross-validation procedure/TSS calculation), important to have an order

	  	# Stats
  		# nrow(dat[which(dat$obs == 1), ]) # pres
  		# nrow(dat[which(dat$obs == 0), ]) # abs

  		# Plot example
  		# plot(dat[which(dat$obs == 0), ]$x, dat[which(dat$obs == 0), ]$y , col = "blue")
  		# points(dat[which(dat$obs == 1), ]$x, dat[which(dat$obs == 1), ]$y , col = "red")
		} # Close else condition

		# Define output of parallel computing: parallel stores table with presences and absences in the list of taxa
		dat

	} # Close parallel computing across taxa
	stopCluster(cl)

	# Statistics
	length(pres_abs_results)
	lapply(pres_abs_results, function(x) nrow(x[which(x$obs==1),])) # presences
	lapply(pres_abs_results, function(x) nrow(x[which(x$obs==0),])) # pseudoabsences
	as.numeric(unlist(lapply(pres_abs_results, function(x) nrow(x[which(x$obs==0),]))))/as.numeric(unlist(lapply(pres_abs_results, function(x) nrow(x[which(x$obs==1),])))) # absences to presences ratio, ideally ≥ ten (Barbet-Massin et al., 2010)

	# Add names to list of species
	names(pres_abs_results) <- as.character(df.taxa.obs$taxon)

	# Save list
	fln <- paste0("./4_Generate_absences/Output_1/",d+6," ",nms.sets[d],"(",vec.strat[1],"_",vec.strat[2],").RData")
	save(pres_abs_results,file=fln)

# Close loop across data strategies (different background and thinning)
}







##==============================================================================
## Now check the object you have created

## Check results
wd_DE <- "/net/kryo/work/deriksson/notion/4_Generate_absences/Output_1/"
wd_DR <- "/net/gladiolus/backups/lunaria/work--2021_05_25_17.00.01--complete/rdamiano/legacy/Righetti_et_al_notion/4_Generate_absences/Output_1/"

## Create data frame from my data
setwd(wd_DE)
dir()
# [3] "3 pres_abs,gr_bg_overl(T_MLD1).RData"
# [4] "4 counts_abs,gr_bg_overl(T_MLD1).RData"
# [5] "5 pres_abs,cr_bg_overl(T_MLD1).RData"

df_1_DE <- get(load("3 pres_abs,gr_bg_overl(T_MLD1).RData"))
df_2_DE <- get(load("4 counts_abs,gr_bg_overl(T_MLD1).RData"))
df_3_DE <- get(load("5 pres_abs,cr_bg_overl(T_MLD1).RData"))

setwd(wd_DR)
df_1_DR <- get(load("3 pres_abs,gr_bg_overl(T_MLD1).RData"))
df_2_DR <- get(load("4 counts_abs,gr_bg_overl(T_MLD1).RData"))
df_3_DR <- get(load("5 pres_abs,cr_bg_overl(T_MLD1).RData"))

# Precence data
setwd(wd_DE)
l1 <- lapply(df_1_DE, function(x) nrow(x[which(x$obs == 1), ])) # presences
l2 <- lapply(df_1_DE, function(x) nrow(x[which(x$obs == 0), ])) # absences
l3 <- lapply(df_1_DR, function(x) nrow(x[which(x$obs == 1), ])) # presences
l4 <- lapply(df_1_DR, function(x) nrow(x[which(x$obs == 0), ])) # absences
# l1; l2; l3; l4

df1 <- data.frame("presence_DE" =sort(unlist(l1)),
  "absences_DE" = sort(unlist(l2)))

df2 <- data.frame("presence_DR" =sort(unlist(l3)),
  "absences_DR" = sort(unlist(l4)))

# # which taxa are missing
# mis_names <- rownames(df1)[!rownames(df1)%in%rownames(df2)]
# # add missing taxa
# df2[mis_names[1], ] <- NA
# df2[mis_names[2], ] <- NA
# df2[mis_names[3], ] <- NA
# arrange(df2, rownames(df2))

# merge
df_presence_Comp_DE_DR <- data.frame(
  df1[ order(row.names(df1)), ],
  df2[ order(row.names(df2)), ]
)


##==============================================================================
# nifH counts --> DIFFERENCES IN FEW PRESENCES COMPARED TO DAMIANOS RESULT
l1 <- lapply(df_2_DE, function(x) nrow(x[which(x$obs == 1), ])) # presences
l2 <- lapply(df_2_DE, function(x) nrow(x[which(x$obs == 0), ])) # absences
l3 <- lapply(df_2_DR, function(x) nrow(x[which(x$obs == 1), ])) # presences
l4 <- lapply(df_2_DR, function(x) nrow(x[which(x$obs == 0), ])) # absences


## Check results
df1 <- data.frame("counts_DE" =sort(unlist(l1)),
  "absences_DE" = sort(unlist(l2)))

df2 <- data.frame("counts_DR" =sort(unlist(l3)),
  "absences_DR" = sort(unlist(l4)))

# # which taxa are missing
# mis_names <- rownames(df1)[!rownames(df1)%in%rownames(df2)]
# # add missing taxa
# df2[mis_names[1], ] <- NA
# df2[mis_names[2], ] <- NA
# df2[mis_names[3], ] <- NA
# arrange(df2, rownames(df2))

# merge
df_count_Comp_DE_DR <- data.frame(
  df1[ order(row.names(df1)), ],
  df2[ order(row.names(df2)), ]
)

# Perfect, same result as Damiano
###==========================================================================

# nifH counts
l1 <- lapply(df_3_DE, function(x) nrow(x[which(x$obs == 1), ])) # presences
l2 <- lapply(df_3_DE, function(x) nrow(x[which(x$obs == 0), ])) # absences
l3 <- lapply(df_3_DR, function(x) nrow(x[which(x$obs == 1), ])) # presences
l4 <- lapply(df_3_DR, function(x) nrow(x[which(x$obs == 0), ])) # absences


## Check results
df1 <- data.frame("counts_DE" =sort(unlist(l1)),
  "absences_DE" = sort(unlist(l2)))

df2 <- data.frame("counts_DR" =sort(unlist(l3)),
  "absences_DR" = sort(unlist(l4)))

# # which taxa are missing
# mis_names <- rownames(df1)[!rownames(df1)%in%rownames(df2)]
# # add missing taxa
# df2[mis_names[1], ] <- NA
# df2[mis_names[2], ] <- NA
# df2[mis_names[3], ] <- NA
# arrange(df2, rownames(df2))

# merge
df_count_Comp_DE_DR <- data.frame(
  df1[ order(row.names(df1)), ],
  df2[ order(row.names(df2)), ]
)
