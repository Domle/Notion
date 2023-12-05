##  In this script we calculate the jaccard dissimilarity index for each model. We further save
##	the two components that make up the jaccard dissimilarity, species turnover and nestedness. We will perform
## 	analysis only on the total dataset, and not on the splits of dataset based on methodology. We will 
##	calculate jaccard, species turnover and nestedness also for the cyanobacterial and non-cyanobacterial fraction of 
##	diazotroph species included.

## Input files:
##      1. Majority vote backtransformed data: To calculate the betadiversity we need presence and absence
##			inputs. Since we already created an ensemble across predictor sets, we use a dataset that is based on 
##			a majority vote, where a presence has been assigned if more than 50% of the models agreed on a presence.

## Output files:
##		2. CSV files containing longitude, latitude, jaccard index, species turnover and nestedness.


## Strategy:    We compute betadiversity from one grid cell with all the other grid cells available and then we compute an 
##			    compute an average betadiversity for that grid cell. We use the beta.pair package to calculate jaccard dissimilarity.

## Author:  Dominic Eriksson
##          Environmental Physics Group, UP
##          ETH Zurich
##          Switzerland

## deriksson@ethz.ch; 3rd of May 2023 ---------------------------------------------------------------------------------


### =========================================================================
### Initialize system
### =========================================================================
rm(list = ls())

# Load functions
source("/home/deriksson/Projects/Notion_DE/Code/Functions/Righetti/common_functions.R") # Load functions to derive latitudinal zonal mean pattern
source("/home/deriksson/Projects/Notion_DE/Code/Functions/Dominic/call_libraries.R")

# Libraries
lib_vec <- c("fields", "raster", "maps", "sfsmisc", "doParallel", "purrr", "tidyverse", "betapart", "reshape2", "data.table")
call_libraries(lib_vec)

# Directories
input.dir_backtransformed <- "/net/kryo/work/deriksson/Projects/Notion_DE/Code/7_SDM_prep/2a_Output/PA_majority_vote/"
output.dir_backtransformed <- "/net/kryo/work/deriksson/Projects/Notion_DE/Code/8_Diversity/5_Output/"

## Load data
# Load ocean mask
msk <- subset(get(load("/home/deriksson/Projects/Notion_DE/Data/Predictor_variables/Environmental_data/VarSet07.2016.RData"))[[1]], "Mask") ## Define mask (continents)
na.layer <- raster()
na.layer[] <- NA

# Number of ensemble member models per taxon (n = 5)
nr <- 5
nr.models <- 5
id.months <- list(); id.months[[1]]<-c(12,1,2); id.months[[2]]<-c(3,4,5); id.months[[3]]<-c(6,7,8); id.months[[4]]<-c(9,10,11) ## Definition of the seasons
vec.seas <- c("winter (boreal)", "spring (boreal)", "summer (boreal)", "autumn (boreal)")

# Model threshold
thres <- 0.3

## Preparatory: load a list with the evaluation metrics for model predictive skill (just to have a list of the species)
## Names of algorithms used (low, intermediate, high complexity)
mod.type <- c("gam", "glm", "rf")
v <- 1

## Names of projection types
proj.type <- c("pa", "prb")
u <- 1

## Folders
vec_folders <- c("Total_dataset", "MicroscopyBased_dataset", "SequenceBased_dataset")
fo <- 1

# Define data sets (strategies)
nms.sets <- c(
	"pres_abs,cr_bg_nonov",
	"pres_abs,cr_bg_overl",
	"pres_abs,gr_bg_nonov",
	"pres_abs,gr_bg_overl",
	"pres_abs,tot_bg_nonov",
	"pres_abs,tot_bg_overl"
)


# Load tax table
tax_table <- read.csv("/home/deriksson/Projects/Notion_DE/Data/commons/Taxon_table.csv", sep = ";")


### =========================================================================
### 1. CALCULATE Betadiversity monthly mean
### =========================================================================
# Set proj.type to 'pa'
u <- 1

# Loop across sdm algorithms
for(v in 1:3){

	# Loop across data sets (strategies)
	for (d in 1:length(nms.sets) ){

        # Read in data
		ls.dat.mean <- get(load(  paste0(input.dir_backtransformed, vec_folders[fo], "/pa/", mod.type[v],"_", proj.type[u], "_paBacktransformed_mean_", nms.sets[d], ".RData") ))

		# Loop across months
		for(k in 1:12){

			# What taxa should be included to do the richness estimate --> We look for only species + gamma.A (phylotype)
			# We need the evaluation file
			if(vec_folders[fo] == "Total_dataset"){
				eval <- read.csv(paste0("/home/deriksson/Projects/Notion_DE/Code/6_SDM_fit/3_Output/", vec_folders[fo], "/", mod.type[v], "_eval_", nms.sets[d], ".csv" ))
				eval <- eval[which(!is.na(eval$"tss.xval4")), ]
				eval <- eval[which(eval$"tss.xval4" > 0.3), ]
				positions <- unique(eval$taxon)}
			if(vec_folders[fo] == "MicroscopyBased_dataset"){
				eval <- read.csv(paste0("/home/deriksson/Projects/Notion_DE/Code/6_SDM_fit/3_Output/", vec_folders[fo], "/", mod.type[v], "_eval_microscopyBased_", gsub("pres_", "", nms.sets[d]), ".csv" ))
				eval <- eval[which(!is.na(eval$"tss.xval4")), ]
				eval <- eval[which(eval$"tss.xval4" > 0.3), ]
				positions <- unique(eval$taxon)}
			if(vec_folders[fo] == "SequenceBased_dataset"){
				eval <- read.csv(paste0("/home/deriksson/Projects/Notion_DE/Code/6_SDM_fit/3_Output/", vec_folders[fo], "/", mod.type[v], "_eval_sequenceBased_", gsub("pres_", "", nms.sets[d]), ".csv" ))
				eval <- eval[which(!is.na(eval$"tss.xval4")), ]
				eval <- eval[which(eval$"tss.xval4" > 0.3), ]
				positions <- unique(eval$taxon)}
					
			## We only keep species level info, we need to check if Gamma A is included in order to add it
			# Add UCYN-A1 and UCYN-A2 if UCYN-A has not been able to be modelled
			if( length(positions[which(positions == "UCYN.A")]) == 0 & length(positions[which(positions == "UCYN.A1")]) == 1 & length(positions[which(positions == "UCYN.A2")]) == 1){		
				#  Add phylotype Gamma-A
				if( length(grep("Gamma.A", positions, value = TRUE)) > 0 ){
					positions <- unique(tax_table[which(tax_table$scientificName_harmonized %in% positions & tax_table$taxonRank_of_scientificName == "SPECIES" | tax_table$scientificName_harmonized == "Gamma.A" ), "scientificName_harmonized"])
				}else{
					positions <- unique(tax_table[which(tax_table$scientificName_harmonized %in% positions & tax_table$taxonRank_of_scientificName == "SPECIES"), "scientificName_harmonized"])}
				# We need to remove Scytonematopsis_pilosa, this is a coastal benthic diazotroph
				if( length(grep("Scytonematopsis_pilosa", positions, value = TRUE)) > 0){
					positions <- positions[-grep("Scytonematopsis_pilosa", positions)]}
				# Add UCYN.A1 and UCYN.A2 if UCYN-A was not modelled successfully
				positions <- c(positions, "UCYN.A1", "UCYN.A2")
			}

			# If UCYN.A was modelled successfully don't include the ecotypes UCYN.A1 and UCYN.A2
			if( length(positions[which(positions == "UCYN.A")]) == 1){		
				#  Add phylotype Gamma-A
				if( length(grep("Gamma.A", positions, value = TRUE)) > 0 ){
					positions <- unique(tax_table[which(tax_table$scientificName_harmonized %in% positions & tax_table$taxonRank_of_scientificName == "SPECIES" | tax_table$scientificName_harmonized == "Gamma.A" ), "scientificName_harmonized"])
				}else{
					positions <- unique(tax_table[which(tax_table$scientificName_harmonized %in% positions & tax_table$taxonRank_of_scientificName == "SPECIES"), "scientificName_harmonized"])}
				# We need to remove Scytonematopsis_pilosa, this is a coastal benthic diazotroph
				if( length(grep("Scytonematopsis_pilosa", positions, value = TRUE)) > 0){
					positions <- positions[-grep("Scytonematopsis_pilosa", positions)]}}

			## Calculate nestedness and turnover
			# Convert raster to dataframe
			raster_stack <- stack(ls.dat.mean[[k]][positions])
			df <- as.data.frame(raster_stack, xy = TRUE)

			# Save locations with NAs, beta.pair function only works without NAs, we rowbind the NAs again later
			na <- df[!complete.cases(df[, 3:ncol(df)]), ]
			na <- na[, c(1:2)]
			na$jac <- NA
			na$jtu <- NA
			na$jne <- NA

			# Remove NAs
			df <- df[which(!is.na(df[, 3])), ]

			# Calculate beta diversity using beta.pair function
			print(paste0("Calculating betadiversity for month ", k))
			beta <- beta.pair(df[, 3:ncol(df)], index.family = "jaccard")

			## Extract Jacchard Index, species turnover, nestedness
			jac <- as.matrix(beta$beta.jac)
			jac <- reshape2::melt(jac)
			# Convert the dataframe to a data.table for efficient processing
			dt <- as.data.table(jac)
			# Specify the columns for which you want to calculate the mean
			# Replace 'column_name' with the actual column name
			columns_to_aggregate <- c("value")
			# Aggregate the mean for specified columns
			jac <- dt[, .(jac = mean(value)), by = Var1]

			jtu <- as.matrix(beta$beta.jtu)
			jtu <- reshape2::melt(jtu)
			# Convert the dataframe to a data.table for efficient processing
			dt <- as.data.table(jtu)
			# Calculate the mean for specified columns
			jtu <- dt[, .(jtu = mean(value, na.rm = TRUE)), by = Var1]

			jne <- as.matrix(beta$beta.jne)
			jne <- reshape2::melt(jne)
			# Convert the dataframe to a data.table for efficient processing
			dt <- as.data.table(jne)
			# Specify the columns for which you want to calculate the mean
			# Replace 'column_name' with the actual column name
			columns_to_aggregate <- c("value")
			# Calculate the mean for specified columns
			jne <- dt[, .(jne = mean(value, na.rm = TRUE)), by = Var1]

			# Create dataframe
			df <- data.frame(
				x = df$x,
				y = df$y,
				jac = jac[, 2],
				jtu = jtu[, 2],
				jne = jne[, 2]
			)

			df_final <- rbind(df, na)
				
			# Save layers (richness map), df (lat), df (lon)
			print("saving file")
			fln <- paste0(output.dir_backtransformed, "Diazotrophs_Betadiversity_", mod.type[v], "_", proj.type[u], "_", nms.sets[d], "_month_", k, ".csv")
			write.csv(df_final, fln, row.names = FALSE)
				
			# Remove objects
			rm(beta)
			rm(jac)
			rm(jtu)
			rm(jne)

		} # close loop across months
	} # Close loop across dataset strategies
} # Close loop across algorithms



### =========================================================================
### 2. CALCULATE Betadiversity cyanos vs. NCDs
### =========================================================================
# Set proj.type to 'pa'
u <- 1


# Loop across sdm algorithms
for(v in 1:3){

	# Loop across data sets (strategies)
	for (d in 1:length(nms.sets) ){

        # Read in data
		ls.dat.mean <- get(load(  paste0(input.dir_backtransformed, vec_folders[fo], "/", proj.type[u],'/', mod.type[v],"_", proj.type[u], "_paBacktransformed_mean_", nms.sets[d], ".RData") ))

		# Loop across months
		for(k in 1:12){

			# What taxa should be included to do the richness estimate --> We look for only species + gamma.A (phylotype)
			# We need the evaluation file
			if(vec_folders[fo] == "Total_dataset"){
				eval <- read.csv(paste0("/home/deriksson/Projects/Notion_DE/Code/6_SDM_fit/3_Output/", vec_folders[fo], "/", mod.type[v], "_eval_", nms.sets[d], ".csv" ))
				eval <- eval[which(!is.na(eval$"tss.xval4")), ]
				eval <- eval[which(eval$"tss.xval4" > 0.3), ]
				positions <- unique(eval$taxon)}
			if(vec_folders[fo] == "MicroscopyBased_dataset"){
				eval <- read.csv(paste0("/home/deriksson/Projects/Notion_DE/Code/6_SDM_fit/3_Output/", vec_folders[fo], "/", mod.type[v], "_eval_microscopyBased_", gsub("pres_", "", nms.sets[d]), ".csv" ))
				eval <- eval[which(!is.na(eval$"tss.xval4")), ]
				eval <- eval[which(eval$"tss.xval4" > 0.3), ]
				positions <- unique(eval$taxon)}
			if(vec_folders[fo] == "SequenceBased_dataset"){
				eval <- read.csv(paste0("/home/deriksson/Projects/Notion_DE/Code/6_SDM_fit/3_Output/", vec_folders[fo], "/", mod.type[v], "_eval_sequenceBased_", gsub("pres_", "", nms.sets[d]), ".csv" ))
				eval <- eval[which(!is.na(eval$"tss.xval4")), ]
				eval <- eval[which(eval$"tss.xval4" > 0.3), ]
				positions <- unique(eval$taxon)}						
						
					
			## We only keep species level info, we need to check if Gamma A is included in order to add it
			# Add UCYN-A1 and UCYN-A2 if UCYN-A has not been able to be modelled
			if( length(positions[which(positions == "UCYN.A")]) == 0 & length(positions[which(positions == "UCYN.A1")]) == 1 & length(positions[which(positions == "UCYN.A2")]) == 1){		
				#  Add phylotype Gamma-A
				if( length(grep("Gamma.A", positions, value = TRUE)) > 0 ){
					positions <- unique(tax_table[which(tax_table$scientificName_harmonized %in% positions & tax_table$taxonRank_of_scientificName == "SPECIES" | tax_table$scientificName_harmonized == "Gamma.A" ), "scientificName_harmonized"])
				}else{
					positions <- unique(tax_table[which(tax_table$scientificName_harmonized %in% positions & tax_table$taxonRank_of_scientificName == "SPECIES"), "scientificName_harmonized"])}
				# We need to remove Scytonematopsis_pilosa, this is a coastal benthic diazotroph
				if( length(grep("Scytonematopsis_pilosa", positions, value = TRUE)) > 0){
					positions <- positions[-grep("Scytonematopsis_pilosa", positions)]}
				# Add UCYN.A1 and UCYN.A2 if UCYN-A was not modelled successfully
				positions <- c(positions, "UCYN.A1", "UCYN.A2")
			}

			# Separate CD's and NCD's
			positions_cd <- tax_table[which(tax_table$Group == 'Cyano'), 'scientificName_harmonized']
			positions_cd <- positions_cd[which(positions_cd %in% positions)]
			#
			positions_ncd <- tax_table[which(tax_table$Group == 'NCD'), 'scientificName_harmonized']
			positions_ncd <- positions_ncd[which(positions_ncd %in% positions)]

			# If UCYN.A was modelled successfully don't include the ecotypes UCYN.A1 and UCYN.A2
			if( length(positions_cd[which(positions_cd == "UCYN.A")]) == 1){		
				#  Add phylotype Gamma-A
				if( length(grep("Gamma.A", positions_cd, value = TRUE)) > 0 ){
					positions_cd <- unique(tax_table[which(tax_table$scientificName_harmonized %in% positions_cd & tax_table$taxonRank_of_scientificName == "SPECIES" | tax_table$scientificName_harmonized == "Gamma.A" ), "scientificName_harmonized"])
				}else{
					positions_cd <- unique(tax_table[which(tax_table$scientificName_harmonized %in% positions_cd & tax_table$taxonRank_of_scientificName == "SPECIES"), "scientificName_harmonized"])}
				# We need to remove Scytonematopsis_pilosa, this is a coastal benthic diazotroph
				if( length(grep("Scytonematopsis_pilosa", positions, value = TRUE)) > 0){
					positions_cd <- positions_cd[-grep("Scytonematopsis_pilosa", positions_cd)]}
			}


			## Calculate nestedness and turnover
			# Convert raster to dataframe
			raster_stack_cd <- stack(ls.dat.mean[[k]][positions_cd])
			raster_stack_ncd <- stack(ls.dat.mean[[k]][positions_ncd])
			df_cd <- as.data.frame(raster_stack_cd, xy = TRUE)
			df_ncd <- as.data.frame(raster_stack_ncd, xy = TRUE)

			# Save locations with NAs, beta.pair function only works without NAs, we rowbind the NAs again later
			na_cd <- df_cd[!complete.cases(df_cd[, 3:ncol(df_cd)]), ]
			na_cd <- na_cd[, c(1:2)]
			na_cd$jac_cd <- NA
			na_cd$jtu_cd <- NA
			na_cd$jne_cd <- NA
			#
			na_ncd <- df_ncd[!complete.cases(df_ncd[, 3:ncol(df_ncd)]), ]
			na_ncd <- na_ncd[, c(1:2)]
			na_ncd$jac_ncd <- NA
			na_ncd$jtu_ncd <- NA
			na_ncd$jne_ncd <- NA
					
			# Remove NAs
			df_cd <- df_cd[which(!is.na(df_cd[, 3])), ]
			df_ncd <- df_ncd[which(!is.na(df_ncd[, 3])), ]

			# Calculate beta diversity using beta.pair function
			print(paste0("Calculating betadiversity for month ", k))
			beta_cd <- beta.pair(df_cd[, 3:ncol(df_cd)], index.family = "jaccard")
			beta_ncd <- beta.pair(df_ncd[, 3:ncol(df_ncd)], index.family = "jaccard")

			## Extract Jacchard Index, species turnover, nestedness
			jac_cd <- as.matrix(beta_cd$beta.jac)
			jac_cd <- reshape2::melt(jac_cd)
			# Convert the dataframe to a data.table for efficient processing
			dt <- as.data.table(jac_cd)
			# Specify the columns for which you want to calculate the mean
			# Replace 'column_name' with the actual column name
			columns_to_aggregate <- c("value")
			# Aggregate the mean for specified columns
			jac_cd <- dt[, .(jac_cd = mean(value)), by = Var1]
			#
			jac_ncd <- as.matrix(beta_ncd$beta.jac)
			jac_ncd <- reshape2::melt(jac_ncd)
			# Convert the dataframe to a data.table for efficient processing
			dt <- as.data.table(jac_ncd)
			# Specify the columns for which you want to calculate the mean
			# Replace 'column_name' with the actual column name
			columns_to_aggregate <- c("value")
			# Aggregate the mean for specified columns
			jac_ncd <- dt[, .(jac_ncd = mean(value)), by = Var1]

			jtu_cd <- as.matrix(beta_cd$beta.jtu)
			jtu_cd <- reshape2::melt(jtu_cd)
			# Convert the dataframe to a data.table for efficient processing
			dt <- as.data.table(jtu_cd)
			# Calculate the mean for specified columns
			jtu_cd <- dt[, .(jtu_cd = mean(value, na.rm = TRUE)), by = Var1]
			#
			jtu_ncd <- as.matrix(beta_ncd$beta.jtu)
			jtu_ncd <- reshape2::melt(jtu_ncd)
			# Convert the dataframe to a data.table for efficient processing
			dt <- as.data.table(jtu_ncd)
			# Calculate the mean for specified columns
			jtu_ncd <- dt[, .(jtu_ncd = mean(value, na.rm = TRUE)), by = Var1]

			jne_cd <- as.matrix(beta_cd$beta.jne)
			jne_cd <- reshape2::melt(jne_cd)
			# Convert the dataframe to a data.table for efficient processing
			dt <- as.data.table(jne_cd)
			# Specify the columns for which you want to calculate the mean
			# Replace 'column_name' with the actual column name
			columns_to_aggregate <- c("value")
			# Calculate the mean for specified columns
			jne_cd <- dt[, .(jne_cd = mean(value, na.rm = TRUE)), by = Var1]
			#
			jne_ncd <- as.matrix(beta_ncd$beta.jne)
			jne_ncd <- reshape2::melt(jne_ncd)
			# Convert the dataframe to a data.table for efficient processing
			dt <- as.data.table(jne_ncd)
			# Specify the columns for which you want to calculate the mean
			# Replace 'column_name' with the actual column name
			columns_to_aggregate <- c("value")
			# Calculate the mean for specified columns
			jne_ncd <- dt[, .(jne_ncd = mean(value, na.rm = TRUE)), by = Var1]

			# Create dataframe
			df_cd <- data.frame(
				x = df_cd$x,
				y = df_cd$y,
				jac_cd = jac_cd[, 2],
				jtu_cd = jtu_cd[, 2],
				jne_cd = jne_cd[, 2]
			)
			#
			df_cd <- rbind(df_cd, na_cd)
			df_cd$id <- paste0( df_cd$x, '_', df_cd$y )

			df_ncd <- data.frame(
				x = df_ncd$x,
				y = df_ncd$y,
				jac_ncd = jac_ncd[, 2],
				jtu_ncd = jtu_ncd[, 2],
				jne_ncd = jne_ncd[, 2]
			)
			#
			df_ncd <- rbind(df_ncd, na_ncd)
			df_ncd$id <- paste0( df_ncd$x, '_', df_ncd$y )

			# Column bind dataframes
			df_final <- merge(df_cd, df_ncd, by = "id", all = TRUE)
				
			# Save layers (richness map), df (lat), df (lon)
			print("saving file")
			fln <- paste0(output.dir_backtransformed, "Diazotrophs_CyanosVsNCDs_Betadiversity_", mod.type[v], "_", proj.type[u], "_", nms.sets[d], "_month_", k, ".csv")
			write.csv(df_final, fln, row.names = FALSE)
				
			# Remove objects
			rm(beta_cd)
			rm(beta_ncd)
			rm(jac_cd)
			rm(jac_ncd)
			rm(jtu_cd)
			rm(jtu_ncd)
			rm(jne_cd)
			rm(jne_ncd)
				

		} # close loop across months
	} # Close loop across dataset strategies
} # Close loop across algorithms


###==============================================================
### END
###==============================================================







