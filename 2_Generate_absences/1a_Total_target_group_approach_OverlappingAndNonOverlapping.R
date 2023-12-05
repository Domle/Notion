### GENERATE PSEUDOABSENCES (TARGET-GROUP BACKGROUND, SELECTION VIA ENVIRONMENTALLY STRATIFIED DESIGN)

# Input files:  1. Diazotroph database (we use three different files here, one dataset the total dadtaset, and the other two 
#               separate between microscopy and sequence-based origins of the observations) (prepared as list of taxa) containing 
#               gridded presences - unthinned or previously thinned on a spatial or monthly-spatial basis
#               2. Environmental variable stack containing data for each monthly (climatological) 1°-cell in the ocean
#               3. The Phytobase dataset (Righetti et al. 2020)

# Output files: Data frames (list of taxa) containing presences and pseudoabsences for each taxon ≥ 15 obs (else Null element). I.e. taxa below 
#               15 obs are replaced by NULL elements, keeping the original length or number of taxa

# Note: The standard approach for the Letter (Sci Adv, Righetti et al., 2019) was "group background", "sites overlapping", and "gam"


# Strategy: Generate pseudoabsences using the target-group approach (Phillips et al., 2009) for each species with ≥ 14 observations.
#           Absences are generated proportionally to the number of species presences within environmental strata of species in a defined target group 
#           (comprising also points of species <15 observations). The environmental stratification shall ensure that environmental space along two key 
#           variables is reasonably addressed in terms of balancing presences and pseudo-absences. The density surface of points/cells which 
#           pseudo-absences are selected from is defined as all SPECIES PRESENCES (monthly cells 1° lat x 1° lon) of a corresponding larger "target-group" 
#           (class, phylum or total plankton). Counting presences of multiple species at same location count as one "presence" or "sampling cell" 
#           (here termed "sites", preferred approach) as taking species-resolved presences ("points") induces some sort of richness effect)
#           Distribution of pseudo-absences per taxon shall thereby exhibit similar bias as its presences (Phillips et al., 2009, Ecol. Appl.)
#           Backgrounds of the three largest groups are defined group-internally (i.e. taking presences of species within these groups as potential absences)
#           Backgrounds for smaller groups are defined by total occurrences excluding the diatoms (that have a strongly imbalanced sampling pattern in the NA)
#           We perform the procedure three times, one time using the complete data sets, and one time splitting the data
#           in two, one dataset only including microscopy-based observations and one dataset only containing 
#           sequence based data.


# Author:	Dominic Eriksson
#			Environmental Physics Group, UP
#			ETHZ, Zürich
#			Switzerland

# script developed by Damiano Righetti

# deriksson@ethz.ch, 4th of September 2023---------------------------------------------------------

### =========================================================================
### Preparatory steps
### =========================================================================
# Clear workspace
rm(list = ls())

# Packages
lib_vec <- c("raster", "classInt", "doParallel", "dismo")
sapply(lib_vec, library, character.only = TRUE)

# Directories
wd_dat <- "/net/kryo/work/deriksson/Projects/Notion_DE/Code/3_Match_up/2_Output/"
wd_dat_phytobase <- "/net/kryo/work/deriksson/Projects/Notion_DE/Data/Righetti_et_al_2020/Original_data_download/"
wd_out <- "/net/kryo/work/deriksson/Projects/Notion_DE/Code/4_Generate_absences/1_Output/"


### =========================================================================
### Load data
### =========================================================================
# Input gridded plankton data (thinning strategies) (n=2 elements). Three times based on dataset type
lisi.dat.s <- list()
lisi.dat.s[[1]] <- read.csv( paste0( wd_dat, "Obs_gridded_diazo_pres_abs.csv"))                
lisi.dat.s[[2]] <- read.csv( paste0( wd_dat, "Obs_gridded_diazo_microscopyBased_pres_abs.csv"))
lisi.dat.s[[3]] <- read.csv( paste0( wd_dat, "Obs_gridded_diazo_sequenceBased_pres_abs.csv"))
lapply(lisi.dat.s, nrow)

# Input ungridded plankton data from PHYTOBASE (to construct a generic target-group approach, including all methods - here sources - that are capable of detecting phytoplankton)
lisi.dat.raw <- list()
lisi.dat.raw[[1]] <- get(load('/home/deriksson/Projects/Notion_DE/Data/Righetti_et_al_2020/Copy_to_work_with_raw/PhytoBase.RData'))
lisi.dat.raw[[2]] <- get(load('/home/deriksson/Projects/Notion_DE/Data/Righetti_et_al_2020/Copy_to_work_with_raw/PhytoBase.RData'))
lisi.dat.raw[[3]] <- get(load('/home/deriksson/Projects/Notion_DE/Data/Righetti_et_al_2020/Copy_to_work_with_raw/PhytoBase.RData'))
lapply(lisi.dat.raw, nrow)

# Concept: 2 background selection strategies * 3 data type strategies (total observartions, microscopy based and sequence based annotations)
nms.sets <- c(
"pres_abs,tot_bg_overl",
"microscopyBased_abs,tot_bg_overl",
"sequenceBased_abs,tot_bg_overl"
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

# Specify parameters to stratify the sampled environment: ** ADJUST IN CASE **
vec.strat <- c("T","MLD1") 

# Preparatory
random_num <- seq(0, 1, 0.01)

### =========================================================================
### Overlapping: Loop across data * background selection strategies: generate pseudo-absences proportional to presences within a maximum of 81 environmental strata
### =========================================================================

# Loop across data set or background selection strategies
for(d in seq_along(nms.sets)){
  # Strategies
	bg.vec <- c("tot", "tot", "tot") # target-group approach (x 3 data-types)
        # d = 1 --> total observations
        # d = 2 --> observations originating from microscopy based annotations
        # d = 3 --> observations originating from sequence based annotations

	# Get gridded data and raw data
	dat.s <- lisi.dat.s[[d]]
	dat.r <- lisi.dat.raw[[d]]

	# Select and grid dat.r (PhytoBase) excluding records below MLD and using species level records only
	dat.raw <- dat.r[which(is.na(dat.r$depth)|dat.r$depth <= dat.r$MLD1), ]
	dat.raw <- dat.raw[which(dat.raw$taxonRank == "SPECIES" | dat.raw$taxonRank == "GENUS"),] # Note, occurrenceStatus are all "PRESENT", hence presence is not selected for in addition
	dat.raw$x <- dat.raw$decimalLongitude
	dat.raw$y <- dat.raw$decimalLatitude
	dat.raw$x<-round(dat.raw$x +0.5)-0.5
	dat.raw$y<-round(dat.raw$y +0.5)-0.5
	dat.raw <- dat.raw[!duplicated(dat.raw[ , c("scientificName", "x", "y", "month")]) ,  ]
	dat.raw <- dat.raw[ -which(dat.raw$class=="Dinophyceae" | dat.raw$class=="Bacillariophyceae"), ] # ** CRITICAL DECISION: Exclude classes (diato, dino) which do not fit diazo size class and sampling design)

	# Prepare a dataframe of the focal taxa (for which background data shall be sampled)
	df.taxa.obs <- data.frame("taxon" = unique(as.character(lisi.dat.s[[d]]$scientificName)))
	dat.s$scientificName <- as.character(dat.s$scientificName)
	dat.splitted <- split(dat.s, f = dat.s$scientificName)
	#
	df.taxa <- data.frame("taxon" = names(unlist(lapply(dat.splitted, nrow))), "obs" = as.numeric(unlist(lapply(dat.splitted, nrow))) )
	df.taxa <- df.taxa[order(df.taxa$taxon),]
	rownames(df.taxa) <- 1:nrow(df.taxa)
	df.taxa.obs$obs <- df.taxa[match(df.taxa.obs$taxon, df.taxa$taxon),][,2]
	df.taxa.obs

	# ---------------------------------------- PREPARE POTENTIAL BACKGROUND SURFACES FOR THE SPECIES: TOTAL TARGET GROUP SITES (monthly cells) ------------------------------------
	# Prepare total target group background and species data: stratification, using total species gridded presences AND gridded non-detections ** ADJUST IN CASE **
	all_spec_id <- rbind( dat.s[, colnames(dat.s)[ which(colnames(dat.s)%in%colnames(dat.raw) )]] , dat.raw[, colnames(dat.s)[ which(colnames(dat.s)%in%colnames(dat.raw) )]] ) # species informatoin retained, add phytobase
	all_spec_id <- all_spec_id[!duplicated(all_spec_id[     ,   c("scientificName", "x", "y", "month")]), ] # Remove duplicates with respect to useful columns (excluding the column phylum and class)

	# Choice of two stratifying parameters defined above in vec.strat, remove NAs from background data regarding these two variables
	if (length(which(is.na(   all_spec_id[,which(names(all_spec_id)==vec.strat[1])]  )))!= 0) {all_spec_id <- all_spec_id[-which(is.na(    all_spec_id[,which(names(all_spec_id)==vec.strat[1])]    )),]} # remove NAs regarding Var 1
	if (length(which(is.na(   all_spec_id[,which(names(all_spec_id)==vec.strat[2])]  )))!= 0) {all_spec_id <- all_spec_id[-which(is.na(    all_spec_id[,which(names(all_spec_id)==vec.strat[2])]    )),]} # remove NAs regarding Var 2
	if (length(which(is.na(   all_spec_id[,which(names(all_spec_id)==vec.strat[1])]  ))) == 0 & length(which(is.na(   all_spec_id[,which(names(all_spec_id)==vec.strat[2])]   ))) == 0) {all_spec_id <- all_spec_id}

	# Specify the range of variable-values and strata into which the values fall (this will allow to drive sampling of absences proportional to overall presences in strata)
	x_envir <- all_spec_id[,which(names(all_spec_id)==vec.strat[1])]
	y_envir <- all_spec_id[,which(names(all_spec_id)==vec.strat[2])]
	# Split ranges into environmental strata
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

	# ----------------------------------------------------------- Parallel computing across species to generate pseudo-absences -----------------------------------------------------------
	n.cores <- detectCores() - 4
	cl = makeCluster(n.cores, outfile = "")
	registerDoParallel(cl)
	pres_abs_results <- foreach (s = 1:length(df.taxa.obs$taxon), .packages = c('dismo', 'classInt'))%dopar%{

		# A. Total target-group background (tot bg), collapse the prepared background cells to monthly cells (exclude taxon information)
		collapsed_ids <- all_species_ids_reclass[ !duplicated(all_species_ids_reclass[ , c("x", "y", "month")]), ]

		# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		# Preparatory 2: Skip taxa below minimum number of presences (e.g. 15)
		de <- all_species_ids_reclass[all_species_ids_reclass$scientificName == df.taxa.obs$taxon[s], ]
		ob <- nrow(de[which(de$occurrenceStatus == "PRESENT"), ] )
		# ob <- nrow(all_cruise_id[all_cruise_id$scientificName == df.taxa.obs$taxon[s], ]) # Same
		# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

		# Open If: observatios are sufficient
		if (ob < min.obs){
		    print( paste0(s," ",as.character(all_species_ids_reclass[all_species_ids_reclass$scientificName == df.taxa.obs$taxon[s], c("scientificName")]) [1], ".. ", ob," yielded too few obs: skipped")); dat <- NULL

		}else{
  		    print( paste0(s," ",as.character(all_species_ids_reclass[all_species_ids_reclass$scientificName == df.taxa.obs$taxon[s], c("scientificName")]) [1], ": ", ob," obs -------------- "))

            # Preparatory 3: Obtain the presences of the target species (presence cells)
            sp_table <- all_species_ids_reclass[all_species_ids_reclass$scientificName == df.taxa.obs$taxon[s], ]
            sp_table <- sp_table[which(sp_table$occurrenceStatus == "PRESENT"), ]
            #
            if( nrow(sp_table) == 0 ) next

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
            dim(absence_table)

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
            } # end of loop absence_groups

            # 2. Merge presences, denoted with obs = 1 and absences with obs  = 0 (i.e. rowbound absence strata)
            sp_table_4_models <- rbind(data.frame(sp_table, obs=1), data.frame(do.call("rbind", absences_per_group), obs=0))
            dat <- cbind(sp_table_4_models, weights = 1) # Create column with weights = 1; weights are associated with presences and absences for modelling
            prec_abs_stats <- matrix(table(dat$obs)) # Number or absences versus number of presences (does not cause trouble if we do not get 10* more absences)
            weights_abs <- prec_abs_stats[2,1]/prec_abs_stats[1,1] # Ratio of presence count divided by absences count
            dat$weights[dat$obs == 0] <- weights_abs # For observation that are absences we add the weight
            row.names(dat)<-seq(1:nrow(dat)) # # Add row ID (order is important for later cross-validation procedure/TSS calculation), important to have an order

		} # Close else condition

		# Define output of parallel computing: parallel stores table with presences and absences in the list of taxa
		dat

	} # Close parallel computing across taxa
	stopCluster(cl)

	# Statistics
	length(pres_abs_results)
	lapply(pres_abs_results, function(x) nrow(x[which(x$obs==1),])) # presences
	lapply(pres_abs_results, function(x) nrow(x[which(x$obs==0),])) # pseudoabsences
	as.numeric(unlist(lapply(pres_abs_results, function(x) nrow(x[which(x$obs==1),])))) / as.numeric(unlist(lapply(pres_abs_results, function(x) nrow(x[which(x$obs==0),]))))# ratio

	# Add names to list of species
	names(pres_abs_results) <- as.character(df.taxa.obs$taxon)

	# Save list
	if(d == 1){ fln <- paste0(wd_out, "Total_dataset/", nms.sets[d], "(", vec.strat[1], "_", vec.strat[2], ").RData") }
	if(d == 2){ fln <- paste0(wd_out, "MicroscopyBased_dataset/", nms.sets[d], "(", vec.strat[1], "_", vec.strat[2], ").RData") }
	if(d == 3){ fln <- paste0(wd_out, "SequenceBased_dataset/", nms.sets[d], "(", vec.strat[1], "_", vec.strat[2], ").RData") }
	save(pres_abs_results, file = fln)
} # Close loop across data strategies (different background and thinning)




### =========================================================================
### Non-overlapping: Loop across data*background selection strategies: generate pseudo-absences proportional to presences within a maximum of 81 environmental strata
### =========================================================================

# Define background selection strategies, overlapping vs. non-overlapping
nms.sets <- c(
"pres_abs,tot_bg_nonov",
"microscopyBased_abs,tot_bg_nonov",
"sequenceBased_abs,tot_bg_nonov"
)

# Loop across data set or background selection strategies
for(d in seq_along(nms.sets)){

	# Strategies
	bg.vec <- c("tot", "tot", "tot") # target-group approach (x 2 data-types)
        # d = 1 --> total observations
        # d = 2 --> observations originating from microscopy based annotations
        # d = 3 --> observations originating from sequence based annotations

	# Get gridded data and raw data
	dat.s <- lisi.dat.s[[d]]
	dat.r <- lisi.dat.raw[[d]]

	# Select and grid dat.r (PhytoBase)
	dat.raw <- dat.r[which(is.na(dat.r$depth) | dat.r$depth <= dat.r$MLD1), ] # excluding depth records below MLD and using species level records only
	dat.raw <- dat.raw[which(dat.raw$taxonRank == "SPECIES" | dat.raw$taxonRank == "GENUS"),] # Note, occurrenceStatus are all "PRESENT", hence presence is not selected for in addition
	dat.raw$x <- dat.raw$decimalLongitude
	dat.raw$y <- dat.raw$decimalLatitude
	dat.raw$x <- round(dat.raw$x + 0.5) - 0.5
	dat.raw$y <- round(dat.raw$y + 0.5) - 0.5
	dat.raw <- dat.raw[!duplicated(dat.raw[ , c("scientificName", "x", "y", "month")]) ,  ]
	dim(dat.raw) # 269'682
	dat.raw <- dat.raw[ -which(dat.raw$class=="Dinophyceae" | dat.raw$class =="Bacillariophyceae"), ] # ** CRITICAL DECISION: Exclude classes (diato, dino) which do not fit diazo size class and sampling design), corrected from "Cyanophyceae" to "Dinophyceae" in first place, 2022-01-14
	dim(dat.raw) # Now 14'398 monthly absence cells; previously (excluding Cyanophyceae and not Dinophyceae) there were 123 653 potential monthly absence cells

	# Prepare a dataframe of the focal taxa (for which background data shall be sampled)
	df.taxa.obs <- data.frame("taxon" = unique(as.character(lisi.dat.s[[d]]$scientificName)))
	dat.s$scientificName <- as.character(dat.s$scientificName)
	dat.splitted <- split(dat.s, f = dat.s$scientificName)
	#
	df.taxa <- data.frame("taxon" = names(unlist(lapply(dat.splitted, nrow))), "obs" = as.numeric(unlist(lapply(dat.splitted, nrow))) )
	df.taxa <- df.taxa[order(df.taxa$taxon),]
	rownames(df.taxa) <- 1:nrow(df.taxa)
	df.taxa.obs$obs <- df.taxa[match(df.taxa.obs$taxon, df.taxa$taxon),][,2]
	df.taxa.obs 

	# ---------------------------------------- PREPARE POTENTIAL BACKGROUND SURFACES FOR THE SPECIES: TOTAL TARGET GROUP SITES (monthly cells) ------------------------------------
	# Prepare total target group background and species data: stratification, using total species gridded presences AND gridded non-detections ** ADJUST IN CASE **
	all_spec_id <- rbind( dat.s[, colnames(dat.s)[ which(colnames(dat.s)%in%colnames(dat.raw) )]] , dat.raw[, colnames(dat.s)[ which(colnames(dat.s)%in%colnames(dat.raw) )]] ) # species informatoin retained, add phytobase
	all_spec_id <- all_spec_id[!duplicated(all_spec_id[     ,   c("scientificName", "x", "y", "month")]), ] # Remove duplicates with respect to useful columns (excluding the column phylum and class)


	# Choice of two stratifying parameters defined above in vec.strat, remove NAs from background data regarding these two variables
	if (length(which(is.na(   all_spec_id[,which(names(all_spec_id)==vec.strat[1])]  )))!= 0) {all_spec_id <- all_spec_id[-which(is.na(    all_spec_id[,which(names(all_spec_id)==vec.strat[1])]    )),]} # remove NAs regarding Var 1
	if (length(which(is.na(   all_spec_id[,which(names(all_spec_id)==vec.strat[2])]  )))!= 0) {all_spec_id <- all_spec_id[-which(is.na(    all_spec_id[,which(names(all_spec_id)==vec.strat[2])]    )),]} # remove NAs regarding Var 2
	if (length(which(is.na(   all_spec_id[,which(names(all_spec_id)==vec.strat[1])]  ))) == 0 & length(which(is.na(   all_spec_id[,which(names(all_spec_id)==vec.strat[2])]   ))) == 0) {all_spec_id <- all_spec_id}

	# Specify the range of variable-values and strata into which the values fall (this will allow to drive sampling of absences proportional to overall presences in strata)
	x_envir <- all_spec_id[,which(names(all_spec_id)==vec.strat[1])]
	y_envir <- all_spec_id[,which(names(all_spec_id)==vec.strat[2])]
	# Split ranges into environmental strata
	breaks <- 9
	# Create matrix that divides range into 9 equal parts; with two variables we get a maximum of 81 strata
	x_breaks <- classIntervals(x_envir, breaks, style = "equal")
	x_matrix <- cbind(x_breaks$brks[1:breaks],x_breaks$brks[2:(breaks + 1)],ID = 1:breaks)
	y_breaks <- classIntervals(y_envir,breaks, style = "equal")
	y_matrix <- cbind(y_breaks$brks[1:breaks],y_breaks$brks[2:(breaks + 1)],ID = 1:breaks)
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
	n.cores <- detectCores() - 4
	cl = makeCluster(n.cores, outfile = "")
	registerDoParallel(cl)
	pres_abs_results <- foreach (s = 1:length(df.taxa.obs$taxon), .packages = c('dismo', 'classInt'))%dopar%{ # Note, packages must be provided specifically to the do parallel command


		# A. Total target-group background (tot bg), collapse the prepared background cells to monthly cells (exclude taxon information)
		# EXCLUDE the presence sampling cells of the focal taxon from the target group sampling-surface, to create "non-overlapping" background. ** ADJUST IN CASE **
		collap_ids <- all_species_ids_reclass[ !duplicated(all_species_ids_reclass[ , c("x", "y", "month")]), ]
		dat_focal_taxon <- all_species_ids_reclass[which(as.character(all_species_ids_reclass$scientificName) == df.taxa.obs$taxon[s]), ] # get species data grom the target group
		pres_cells_focal_taxon <- dat_focal_taxon[which(dat_focal_taxon$occurrenceStatus == "PRESENT"), ]
		collapsed_ids <- collap_ids[ -which(  paste0(collap_ids$x,"_",collap_ids$y,"_",collap_ids$month)%in% paste0(pres_cells_focal_taxon$x,"_", pres_cells_focal_taxon$y,"_", pres_cells_focal_taxon$month) )  , ]

		# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		# Preparatory 2: Skip taxa below minimum number of presences (e.g. 15)
		# ob <- nrow(all_species_ids_reclass[all_species_ids_reclass$scientificName == df.taxa.obs$taxon[s], ]) #
		de <- all_species_ids_reclass[all_species_ids_reclass$scientificName == df.taxa.obs$taxon[s], ] #
		ob <- nrow( de[which(de$occurrenceStatus == "PRESENT"), ] )
		# ob <- nrow(all_cruise_id[all_cruise_id$scientificName == df.taxa.obs$taxon[s], ]) # Same
		# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

		# Open If: observatios are sufficient
		if (ob < min.obs){
		    print( paste0(s," ",as.character(all_species_ids_reclass[all_species_ids_reclass$scientificName == df.taxa.obs$taxon[s], c("scientificName")]) [1], ".. ", ob," yielded too few obs: skipped")); dat <- NULL

		}else{
  		    print( paste0(s," ",as.character(all_species_ids_reclass[all_species_ids_reclass$scientificName == df.taxa.obs$taxon[s], c("scientificName")]) [1], ": ", ob," obs -------------- "))

            # Preparatory 3: Obtain the presences of the target species (presence cells)
            sp_table <- all_species_ids_reclass[all_species_ids_reclass$scientificName== df.taxa.obs$taxon[s], ]
            sp_table <- sp_table[which(sp_table$occurrenceStatus == "PRESENT"), ] # ** ADJUST IN CASE**

            # Preparatory 4: Extract frequencies by which points/sites of the target group fall into environmental strata. Then, derive the number of desired absences for the focal species per stratum
            # Here we include the focal in this estimate of general sampling effort
            xy_rcls_freq <- data.frame(table(collapsed_ids$xy_rcls)/length(collapsed_ids$x_rcls)) # Frequency of points or samples of total target-group-species-occurrences per stratum
            names(xy_rcls_freq)[1] <-"xy_rcls" # Give name to column
            xy_rcls_freq[,1] <- as.numeric(as.character(xy_rcls_freq[, 1])) # Define reclassification frequency as numeric
            table_abs1 <- data.frame(xy_rcls_freq, prop_abs = (dim(sp_table)[1]*10)*xy_rcls_freq$Freq) # Add background points to be sampled per stratum: 10 x more than presences, driven by general sampling effort
            table_abs2 <- data.frame(table_abs1, prop_abs_0 = table_abs1$prop_abs - floor(table_abs1$prop_abs)) # To round desired absences to integer: adds column difference between smaller closest integer and desired number
            table_abs3 <- data.frame(table_abs2, prob = sample(random_num, dim(table_abs1)[1], replace=T)) # Adds column with random number between 0 and 1 (with steps 0.01)
            table_abs4 <- data.frame(table_abs3, absences =1) # Adds column with "1" s
            table_abs4$absences[which(table_abs4$prop_abs_0 > table_abs4$prob)] <- ceiling(table_abs4$prop_abs[which(table_abs4$prop_abs_0 > table_abs4$prob)]) # Round absences up for random subset
            table_abs4$absences[which(table_abs4$prop_abs_0 < table_abs4$prob)] <- floor(table_abs4$prop_abs[which(table_abs4$prop_abs_0 < table_abs4$prob)]) # Round absences down for random subset
            absence_groups <- table_abs4[table_abs4$absences > 0,] # Skip strata without presences

            # Preparatory 5: Select backround data, here including the points/sites of the focal species ('overlapping background') in case of the general target-group approach,
            # i.e. we treat pseudo-absences as a general background, not as confirmed absences
            # we also include the species in case of the more specific target-cruise approach,
            # i.e., here we treat pseudo-absences as a general background, not as confirmed absences, but other options may be tested as well (i.e., excluding the species)
            absence_table <- collapsed_ids[ , colnames(sp_table)] # needed to select exactly same columns (from the cruise specific background), as else it wont match below

            # 1. Select background points for target species in each stratum proportionally to the density of samples/points in the background surface (in a random way)
            absences_per_group <- list()
            for (i in 1:nrow(absence_groups)){
                # 1. Select available absences within stratum in question
                group_absence_table <- absence_table[absence_table$xy_rcls == absence_groups[i, "xy_rcls"],]
                # 2. Define number of absences that can be drawn maximally (if the number is smaller than desired, define maximum possible number of absences (limited by points per stratum):
                absence_num <- ifelse(
                absence_groups[i,"absences"] > dim(absence_table[absence_table$xy_rcls == absence_groups[i,"xy_rcls"],])[1] # Test: is the number of desired background points bigger than the available background points?
                ,
                dim(absence_table[absence_table$xy_rcls == absence_groups[i,"xy_rcls"],])[1] # if TRUE the potential points are insufficient - however, save the number of available points as absence_num
                ,
                absence_groups[i,"absences"] # ELSE: save the number of desired background points as absence_num
                )
                # 3. Randomly sample the background points from the table containing all possible absences for the stratum in question
                sampled_group_absence_table <- group_absence_table[sample(1:dim(group_absence_table)[1], absence_num), ]
                absences_per_group[[i]] <- sampled_group_absence_table # we store the product for the stratum or group
            } # close for loop absence_groups

            # 2. Merge presences, denoted with obs = 1 and absences with obs  = 0 (i.e. rowbound absence strata)
            sp_table_4_models <- rbind(data.frame(sp_table, obs=1), data.frame(do.call("rbind", absences_per_group), obs=0))
            dat <- cbind(sp_table_4_models, weights = 1) # Create column with weights = 1; weights are associated with presences and absences for modelling
            prec_abs_stats <- matrix(table(dat$obs)) # Number or absences versus number of presences (does not cause trouble if we do not get 10* more absences)
            weights_abs <- prec_abs_stats[2,1]/prec_abs_stats[1,1] # Ratio of presence count divided by absences count
            dat$weights[dat$obs == 0] <- weights_abs # For observation that are absences we add the weight
            row.names(dat)<-seq(1:nrow(dat)) # # Add row ID (order is important for later cross-validation procedure/TSS calculation), important to have an order
		} # Close else condition

		# Define output of parallel computing: parallel stores table with presences and absences in the list of taxa
		dat

	} # Close parallel computing across taxa
	stopCluster(cl)

	# Statistics
	length(pres_abs_results)
	lapply(pres_abs_results, function(x) nrow(x[which(x$obs==1),])) # presences
	lapply(pres_abs_results, function(x) nrow(x[which(x$obs==0),])) # pseudoabsences
	as.numeric(unlist(lapply(pres_abs_results, function(x) nrow(x[which(x$obs==1),])))) / as.numeric(unlist(lapply(pres_abs_results, function(x) nrow(x[which(x$obs==0),])))) # ratio

	# Add names to list of species
	names(pres_abs_results) <- as.character(df.taxa.obs$taxon)

	# Save list
	if(d == 1){ fln <- paste0(wd_out, "Total_dataset/", nms.sets[d], "(", vec.strat[1], "_", vec.strat[2], ").RData") }
	if(d == 2){ fln <- paste0(wd_out, "MicroscopyBased_dataset/", nms.sets[d], "(", vec.strat[1], "_", vec.strat[2], ").RData") }
	if(d == 3){ fln <- paste0(wd_out, "SequenceBased_dataset/", nms.sets[d], "(", vec.strat[1], "_", vec.strat[2], ").RData") }
	save(pres_abs_results, file = fln)
} # Close loop across data strategies

###==============================================================
### END
###==============================================================