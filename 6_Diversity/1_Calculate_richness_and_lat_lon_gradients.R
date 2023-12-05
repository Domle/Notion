### CREATE GLOBAL RASTER LAYERS FOR SPECIES RICHNESS AND DERIVE LAT, LON ZONAL MEANS


# - Input files: 
#           1. Lists with monthly projections (each list contains 12 slots or months with species as layers).
#           Note: Member models have been TSS selected in previous code (to construct the ensemble means)

# - Output files:
#		    1. Lists containing (12 monthly) richness projections in 4 slots (mean, minima, maxima, and sd)
# 		    2. Derived latitudinal and longitudinal one-degree zonal means: a. monthly, b. seasonal means, c. annual means. And also the object has length 3, 1= mean, 2= sd, 3 = absolute error
#

# Strategy: Calculate richness (sum of species' projections) and addition: sums of  SD, and an error propagation across the ensembles member models.
            # a. monthly
            # b. seasonal mean of months
            # c. annual mean of months

            # We will compute each richness estimate unnormalized and normalized (divided by the number of species successfully modeled). We do this for the 
            # HSI, presence absence converted data and for the majority vote option. We futher separated and modeled data depending on the methodology used. 
            # There fore we will retrieve richness estimates based on the whole dataset as input, only microscopy based data and sequence based data, which can
            # be later compared to see if richness estimates differe depending on methodology used. We further split the groups cyanobacteria and non-cyanobacteria to 
            # also calculate richness within the two groups. 

            # Sections refer to

            # 1. Monthly non-normalized richness
            # 2. Annual non-normalized richness
            # 3. Monthly normalized richness
            # 4. Annual normalized richness
            # 5. Annual majority vote richness normalized
            # 6. Monthly majority vote richness normalized
            # 7. Annual cyanobacterial and non-cyanobacterial richness normalized



# Author: 	Dominic Eriksson
#			Environmental Physics Group
# 			ETH, Zurich
#			Switzerland


# deriksson@ethz.ch, 15th of September 2022 -------------------------------------------------------------

### =========================================================================
### Initialize system
### =========================================================================
rm(list = ls())

# Load functions
source("/home/deriksson/Projects/Notion_DE/Code/Functions/Righetti/common_functions.R") # Load functions to derive latitudinal zonal mean pattern
source("/home/deriksson/Projects/Notion_DE/Code/Functions/Dominic/call_libraries.R")

# Libraries
lib_vec <- c("fields", "raster", "maps", "sfsmisc", "doParallel", "purrr", "tidyverse")
call_libraries(lib_vec)

# Directories
input.dir <- "/net/kryo/work/deriksson/Projects/Notion_DE/Code/7_SDM_prep/2a_Output/Ensemble_mean/"
input.dir_backtransformed <- "/net/kryo/work/deriksson/Projects/Notion_DE/Code/7_SDM_prep/2a_Output/PA_majority_vote/" # refers to majority vote
output.dir <- "/home/deriksson/Projects/Notion_DE/Code/8_Diversity/1_Output/Ensemble_mean/"
output.dir_backtransformed <- "/home/deriksson/Projects/Notion_DE/Code/8_Diversity/1_Output/PA_majority_vote/" # refers to majority vote

# Load ocean mask
msk <- subset (get(load("/home/deriksson/Projects/Notion_DE/Data/Predictor_variables/Environmental_data/VarSet07.2016.RData"))[[1]], "Mask") ## Define mask (continents)
na.layer <- raster()
na.layer[] <- NA

# Number of ensemble member models per taxon (n = 5)
nr <- 5
nr.models <- 5
id.months <- list(); id.months[[1]] <- c(12,1,2); id.months[[2]] <- c(3,4,5); id.months[[3]] <- c(6,7,8); id.months[[4]] <- c(9,10,11) ## Definition of the seasons
vec.seas <- c("winter (boreal)", "spring (boreal)", "summer (boreal)", "autumn (boreal)")

# Model threshold
thres <- 0.3

## Preparatory: load a list with the evaluation metrics for model predictive skill (just to have a list of the species)
## Names of algorithms used (low, intermediate, high complexity)
mod.type <- c("gam", "glm", "rf")
# v <- 1

## Names of projection types
proj.type <- c("pa", "prb")
# u <- 1

## Folders
vec_folders <- c("Total_dataset", "MicroscopyBased_dataset", "SequenceBased_dataset")
# fo <- 1

# Define data sets (strategies)
nms.sets <- c(
	"pres_abs,cr_bg_nonov",
	"pres_abs,cr_bg_overl",
	"pres_abs,gr_bg_nonov",
	"pres_abs,gr_bg_overl",
	"pres_abs,tot_bg_nonov",
	"pres_abs,tot_bg_overl"
)
# p <- 1

# Load tax table
tax_table <- read.csv("/home/deriksson/Projects/Notion_DE/Data/commons/Taxon_table.csv", sep = ";")



### =========================================================================
### 1. CALCULATE RICHNESS MONTHLY - Non-normalized by taxa
### =========================================================================
# Loop across folders
for(fo in seq_along(vec_folders)){

	for(u in seq_along(proj.type)){

		# Loop across sdm algorithms; v = 1 corresponds to gam; v = 2 glm; v = 3 rf
		for(v in seq_along(mod.type)){

			# Loop over data strategies
			for ( d in 1:length(nms.sets) ){
				
				# Read in data
				ls.dat.mean <- get(load(  paste0(input.dir, vec_folders[fo], "/", proj.type[u], "/", mod.type[v],"_", proj.type[u], "_mean_", nms.sets[d], ".RData") ))

				# Parallel computing across months
				n.cores <- detectCores()/2
				cl = makeCluster(n.cores, outfile = "")
				registerDoParallel(cl)
				list.cl <- foreach(k = 1:12 , .packages = c('raster', 'sfsmisc', 'tidyverse'))%dopar%{

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

					# Let's save the information on how many taxa are included in each richness map
					fln <- paste0(output.dir, vec_folders[fo], "/div_monthly/", proj.type[u], "/", "Table_RichnessTaxaIncluded_", mod.type[v], "_", proj.type[u], "_", nms.sets[d], "_", k, ".csv")
					df_positions <- data.frame(taxa_included = positions)
					write.csv(df_positions, fln, row.names = FALSE) 

					# Calculate richness
					lisi.rasters.mean <- reclassify(calc(stack(ls.dat.mean[[k]][positions]), sum, na.rm = TRUE), rcl = cbind(0, NA))

					# Derive df (lat)
					lisi.dflat.mean <- by_lat(lisi.rasters.mean) # Enable if needed

					# Derive df (lon)
					lisi.dflon.mean <- by_lon(lisi.rasters.mean) # Enable if needed

					# Define output list
					data.stored <- list(lisi.rasters.mean, lisi.dflat.mean, lisi.dflon.mean) # adjust in case

					# Close parallel computing across months
				}
				# Stop parallel computing
				stopCluster(cl)

				# Merge lists
				evthing.raster <- list()
				evthing.raster[[1]] <- lapply(list.cl,'[[',1) # adjust in case

				evthing.dflat <- list()
				evthing.dflat[[1]] <- lapply(list.cl,'[[',2) # adjust in case

				evthing.dflon <- list()
				evthing.dflon[[1]] <- lapply(list.cl,'[[',3) # adjust in case

				# Save layers (richness map), df (lat), df (lon)
				fln <- paste0(output.dir, vec_folders[fo], "/div_monthly/", proj.type[u], "/", "Diazotrophs_speciesRichness_", mod.type[v], "_", proj.type[u], "_", nms.sets[d], "_mon.grd")
				writeRaster(stack(unlist(evthing.raster)), filename = fln, format = "raster", overwrite = TRUE)
				fln <- paste0(output.dir, vec_folders[fo], "/lat_monthly/", proj.type[u], "/", "Diazotrophs_speciesRichness_", mod.type[v], "_", proj.type[u], "_", nms.sets[d], "_mon.RData")
				save(evthing.dflat, file = fln)
				fln <- paste0(output.dir, vec_folders[fo], "/lon_monthly/", proj.type[u], "/", "Diazotrophs_speciesRichness_", mod.type[v], "_", proj.type[u], "_", nms.sets[d], "_mon.RData")
				save(evthing.dflon, file = fln)

				# Progress
				print(paste0("Model ", v, " strategy ", d, " monthly"))	
			} # Close loop across dataset strategies
		} # Close loop across algorithms
	} # end loop across projection type
} # end loop across folders (total dataset, microscopy based and sequence based dataset)


### =========================================================================
### 2. CALCULATE RICHNESS ANNUAL MEAN - Non normalized
### =========================================================================
# Loop
for(fo in seq_along(vec_folders)){

	for(u in seq_along(proj.type)){


		# Loop across sdm algorithms
		for(v in 1:3){

			# Loop across data sets (strategies)
			for (d in 1:length(nms.sets) ){

				# Get frames with predictive skill (each row corresponds to a taxon-specific and variable-specific model fit)
				ls.dat.mean <- get(load(  paste0(input.dir, vec_folders[fo], "/", proj.type[u], "/", mod.type[v],"_", proj.type[u], "_mean_", nms.sets[d], ".RData"))  )

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
							positions <- positions[-grep("Scytonematopsis_pilosa", positions)]}
					}


				# Let's save the information on how many taxa are included in each richness map
				fln <- paste0(output.dir, vec_folders[fo], "/div_monthly/", proj.type[u], "/", "Table_RichnessTaxaIncluded_", mod.type[v], "_", proj.type[u], "_", nms.sets[d], "_annual.csv")
				df_positions <- data.frame(taxa_included = positions)
				write.csv(df_positions, fln, row.names = FALSE) 

				# Prepare diversity rasters and lat frames
				lisi.taxa.mean <- list ()
				for (i in positions){
					lisi.taxa.mean[[length(lisi.taxa.mean) + 1]] <- calc(stack(unlist(lapply(ls.dat.mean, '[', i))), mean, na.rm = TRUE)
					print(paste(i))
				} # close loop across positions

				# Calculate richness layer (sum of the taxa)
				raster.mean <- calc(stack(lisi.taxa.mean), sum, na.rm = TRUE)

				# Merge lists
				evthing.raster <- list()
				evthing.raster[[1]] <- reclassify(raster.mean, rcl = cbind(0, NA))

				evthing.dflat <- list()
				evthing.dflat[[1]] <- by_lat(evthing.raster[[1]])

				evthing.dflon <- list()
				evthing.dflon[[1]] <- by_lon(evthing.raster[[1]])

				# Save layers (richness map), frames (lat), frames (lon)
				fln <- paste0(output.dir, vec_folders[fo], "/div_annually/", proj.type[u], "/", "Diazotrophs_speciesRichness_", mod.type[v], "_", proj.type[u], "_", nms.sets[d], "_annual.grd")
				writeRaster(stack(unlist(evthing.raster)), filename = fln, format = "raster", overwrite = T)
				fln <- paste0(output.dir, vec_folders[fo], "/lat_annually/", proj.type[u], "/", "Diazotrophs_speciesRichness_", mod.type[v], "_", proj.type[u], "_", nms.sets[d], "_annual.RData")
				save(evthing.dflat, file = fln)
				fln <- paste0(output.dir, vec_folders[fo], "/lon_annually/", proj.type[u], "/", "Diazotrophs_speciesRichness_", mod.type[v], "_", proj.type[u], "_", nms.sets[d], "_annual.RData")
				save(evthing.dflon, file = fln)

				# Progress
				print(paste0(" model ", v, " strategy ", d, " annually"))				
			} # Close loop across dataset strategies
		} # Close loop across algorithms
	} # Close loop across proj.type
} # Close loop across vec_folders	



### =========================================================================
### 3. CALCULATE RICHNESS MONTHLY MEAN - Normalized
### =========================================================================
# Loop
for(fo in seq_along(vec_folders)){

	for(u in seq_along(proj.type)){

		for(v in seq_along(mod.type)){

			# Parallel computing across data strategies
			for ( d in 1:length(nms.sets) ){

				# Get frames with predictive skill (each row corresponds to a taxon-specific and variable-specific model fit)
				ls.dat.mean <- get(load(  paste0(input.dir, vec_folders[fo], "/", proj.type[u], "/", mod.type[v],"_", proj.type[u], "_mean_", nms.sets[d], ".RData"))  )

				# Parallel computing across months
				n.cores <- 10
				cl = makeCluster(n.cores, outfile = "")
				registerDoParallel(cl)
				list.cl <- foreach(k = 1:12 , .packages = c('raster', 'sfsmisc', 'tidyverse'))%dopar%{

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
							positions <- positions[-grep("Scytonematopsis_pilosa", positions)]}
					}

					# Calculate richness layer of month k
					lisi.rasters.mean <- reclassify(calc(stack(ls.dat.mean[[k]][positions]), sum, na.rm = TRUE), rcl = cbind(0, NA))
					# Normalize richness (divide by number of species included)
					lisi.rasters.mean_taxa_corrected <- lisi.rasters.mean/length(positions)

					# Derive df (lat)
					lisi.dflat.mean_corrected <- by_lat(lisi.rasters.mean_taxa_corrected)

					# Derive df (lon)
					lisi.dflon.mean_corrected <- by_lon(lisi.rasters.mean_taxa_corrected)

					# Define output list
					data.stored <- list(lisi.rasters.mean_taxa_corrected,	lisi.dflat.mean_corrected, lisi.dflon.mean_corrected) # adjust in case

				} # Close parallel computing across months
				stopCluster(cl)

				# Merge lists
				evthing.raster <- list()
				evthing.raster[[1]] <- lapply(list.cl,'[[',1) # adjust in case

				evthing.dflat <- list()
				evthing.dflat[[1]] <- lapply(list.cl,'[[',2)	 # adjust in case

				evthing.dflon <- list()
				evthing.dflon[[1]] <- lapply(list.cl,'[[',3) # adjust in case

				# Save layers (richness map), frames (lat), frames (lon)
				fln <- paste0(output.dir, vec_folders[fo], "/div_monthly/", proj.type[u], "/", "Diazotrophs_speciesRichness_", mod.type[v], "_", proj.type[u], "_", nms.sets[d], "_mon_normalized.grd")
				writeRaster(stack(unlist(evthing.raster)), filename = fln, format = "raster", overwrite = T)
				fln <- paste0(output.dir, vec_folders[fo], "/lat_monthly/", proj.type[u], "/", "Diazotrophs_speciesRichness_", mod.type[v], "_", proj.type[u], "_", nms.sets[d], "_mon_normalized.RData")
				save(evthing.dflat, file = fln)
				fln <- paste0(output.dir, vec_folders[fo], "/lon_monthly/", proj.type[u], "/", "Diazotrophs_speciesRichness_", mod.type[v], "_", proj.type[u], "_", nms.sets[d], "_mon_normalized.RData")
				save(evthing.dflon, file = fln)
				save(evthing.dflon, file = fln)

				# Progress
				print(paste0( "model ", v, " strategy ", d, " monthly"))
			} # Close loop across dataset strategies
		} # Close loop across algorithms
	} # Close loop over proj.type	
} # Close loop over vec_folders

### =========================================================================
### 4. CALCULATE RICHNESS ANNUAL MEAN - Normalized
### =========================================================================
# Loop
for(fo in seq_along(vec_folders)){
	
	for(u in seq_along(proj.type)){

		for(v in seq_along(mod.type)){

			# Loop across data sets (strategies)
			for (d in seq_along(nms.sets) ){

				# Read in list with 12 stacks: data from brick (to avoid trouble with reloading the files) and allocate into list of 12 monthly slots with pre-gathered projections
				ls.dat.mean <- get(load(  paste0(input.dir, vec_folders[fo], "/", proj.type[u], "/", mod.type[v],"_", proj.type[u], "_mean_", nms.sets[d], ".RData"))  )

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
							positions <- positions[-grep("Scytonematopsis_pilosa", positions)]}
					}

				# Prepare diversity rasters and lat frames
				lisi.taxa.mean <- list ()
				for (i in positions){
					lisi.taxa.mean[[length(lisi.taxa.mean) + 1]] <- calc(stack(unlist(lapply(ls.dat.mean, '[', i))), mean, na.rm = TRUE)
					print(paste(i))
				}

				# Calculate richness layer (sum of the taxa)
				raster.mean <- calc(stack(lisi.taxa.mean), sum, na.rm = TRUE) # adjust in case
				# Normalize by the number of species
				raster.mean_taxa_corrected <- raster.mean/length(positions)

				# Merge lists
				evthing.raster <- list()
				evthing.raster[[1]] <- reclassify(raster.mean_taxa_corrected, rcl = cbind(0, NA))

				evthing.dflat <- list()
				evthing.dflat[[1]] <- by_lat(evthing.raster[[1]])

				evthing.dflon <- list()
				evthing.dflon[[1]] <- by_lon(evthing.raster[[1]])


				# Save layers (richness map), frames (lat), frames (lon)
				fln <- paste0(output.dir, vec_folders[fo], "/div_annually/", proj.type[u], "/", "Diazotrophs_speciesRichness_", mod.type[v], "_", proj.type[u], "_", nms.sets[d], "_annual_normalized.grd")
				writeRaster(stack(unlist(evthing.raster)), filename = fln, format = "raster", overwrite = TRUE)
				fln <- paste0(output.dir, vec_folders[fo], "/lat_annually/", proj.type[u], "/", "Diazotrophs_speciesRichness_", mod.type[v], "_", proj.type[u], "_", nms.sets[d], "_annual_normalized.RData")
				save(evthing.dflat, file = fln)
				fln <- paste0(output.dir, vec_folders[fo], "/lon_annually/", proj.type[u], "/", "Diazotrophs_speciesRichness_", mod.type[v], "_", proj.type[u], "_", nms.sets[d], "_annual_normalized.RData")
				save(evthing.dflon, file = fln)

				# Progress
				print(paste0(" model ", v, " strategy ", d, " annually"))
			} # Close loop across dataset strategies
		} # Close loop across algorithms
	} # Close loop across proj.type
} # Close loop across vec_folders


### =========================================================================
### 5. CALCULATE RICHNESS ANNUAL MEAN Majority vote - NORMALISE THE DIVERSITY FOR EACH TAXON RELATIVE TO TAXON COMPLETENESS
### =========================================================================
# Set projection type
u <- 1

# Loop
for(fo in seq_along(vec_folders)){

	for(v in seq_along(mod.type)){

		# Loop across data sets (strategies)
		for (d in seq_along(nms.sets) ){

			# Read in list with 12 stacks: data from brick (to avoid trouble with reloading the files) and allocate into list of 12 monthly slots with pre-gathered projections
			ls.dat.mean <- get(load(  paste0(input.dir_backtransformed, vec_folders[fo], "/", mod.type[v],"_", proj.type[1], "_paBacktransformed_mean_", nms.sets[d], ".RData"))  )

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
							positions <- positions[-grep("Scytonematopsis_pilosa", positions)]}
					}

			# Let's save the information on how many taxa are included in each richness map
			fln <- paste0(output.dir_backtransformed, vec_folders[fo], "/div_monthly/", proj.type[1], "/", "Table_RichnessTaxaIncluded_", mod.type[v], "_", proj.type[u], "_", nms.sets[d], "_annual.csv")
			df_positions <- data.frame(taxa_included = positions)
			write.csv(df_positions, fln, row.names = FALSE) 

			# Prepare diversity rasters and lat frames
			lisi.taxa.mean <- list ()
			for (i in positions){
				lisi.taxa.mean[[length(lisi.taxa.mean) + 1]] <- calc(stack(unlist(lapply(ls.dat.mean, '[', i))), mean, na.rm = TRUE)
				print(paste(i))
			}


			# Calculate annual mean richness
			raster.mean <- calc(stack(unlist(lisi.taxa.mean)), sum, na.rm = TRUE) # adjust in case
			# Normalize by the number of species
			raster.mean_taxa_corrected <- raster.mean/length(positions)

			# Merge lists
			evthing.raster <- list()
			evthing.raster[[1]] <- reclassify(raster.mean_taxa_corrected,rcl = cbind(0, NA))

			evthing.dflat <- list()
			evthing.dflat[[1]] <- by_lat(evthing.raster[[1]])

			evthing.dflon <- list()
			evthing.dflon[[1]] <- by_lon(evthing.raster[[1]])


			# Save layers (richness map), frames (lat), frames (lon)
			fln <- paste0(output.dir_backtransformed, vec_folders[fo], "/", "div_annually", "/", proj.type[1], "/", "Diazotrophs_speciesRichness_", mod.type[v], "_", nms.sets[d], "_paBacktransformed_taxa_corrected_ann.grd")
			writeRaster(stack(unlist(evthing.raster)), filename = fln, format = "raster", overwrite = TRUE)
			fln <- paste0(output.dir_backtransformed, vec_folders[fo], "/", "lat_annually", "/", proj.type[1], "/", "Diazotrophs_speciesRichness_", mod.type[v], "_", nms.sets[d], "_paBacktransformed_taxa_corrected_ann.RData")
			save(evthing.dflat, file = fln)
			fln <- paste0(output.dir_backtransformed, vec_folders[fo], "/", "lon_annually", "/", proj.type[1], "/", "Diazotrophs_speciesRichness_", mod.type[v], "_", nms.sets[d], "_paBacktransformed_taxa_corrected_ann.RData")
			save(evthing.dflon, file = fln)

			# Progress
			print(paste0(" model ", v, " strategy ", d, " annually"))
				} # Close loop across dataset strategies
	} # Close loop across algorithms
} # Close loop over vec_folders


### =========================================================================
### 6. CALCULATE RICHNESS MONTHLY - Majority vote - NORMALISE THE DIVERSITY FOR EACH TAXON RELATIVE TO TAXON COMPLETENESS
### =========================================================================
# Select projection type
u <- 1

# Loop
for(fo in seq_along(vec_folders)){

	for(v in seq_along(mod.type)){

		# Parallel computing across data strategies
		for ( d in 1:length(nms.sets) ){

			# Get frames with predictive skill (each row corresponds to a taxon-specific and variable-specific model fit)
			ls.dat.mean <- get(load(  paste0(input.dir_backtransformed, vec_folders[fo], "/", mod.type[v],"_", proj.type[u], "_paBacktransformed_mean_", nms.sets[d], ".RData"))  )

			# Parallel computing across months
			n.cores <- 5
			cl = makeCluster(n.cores, outfile="")
			registerDoParallel(cl)
			list.cl <- foreach(k = 1:12 , .packages = c('raster', 'sfsmisc', 'tidyverse'))%dopar%{
					
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
							positions <- positions[-grep("Scytonematopsis_pilosa", positions)]}
					}

				# Let's save the information on how many taxa are included in each richness map
				fln <- paste0(output.dir_backtransformed, vec_folders[fo], "/div_monthly/", proj.type[1], "/", "Table_RichnessTaxaIncluded_", mod.type[v], "_", proj.type[u], "_", nms.sets[d], "_", k, ".csv")
				df_positions <- data.frame(taxa_included = positions)
				write.csv(df_positions, fln, row.names = FALSE) 				

				# Calculate richness layer of month k
				lisi.rasters.mean <- reclassify(calc(stack(ls.dat.mean[[k]][positions]), sum, na.rm = TRUE), rcl = cbind(0, NA))
				# Normalize richness (divide by number of species included)
				lisi.rasters.mean_taxa_corrected <- lisi.rasters.mean/length(positions)

				# Derive df (lat)
				lisi.dflat.mean_corrected <- by_lat(lisi.rasters.mean_taxa_corrected)

				# Derive df (lon)
				lisi.dflon.mean_corrected <- by_lon(lisi.rasters.mean_taxa_corrected)

				# Define output list
				data.stored <- list(lisi.rasters.mean_taxa_corrected,	lisi.dflat.mean_corrected, lisi.dflon.mean_corrected) # adjust in case

			} # Close parallel computing across months
			stopCluster(cl)

			# Merge lists
			evthing.raster <- list()
			evthing.raster[[1]] <- lapply(list.cl,'[[',1) # adjust in case

			evthing.dflat <- list()
			evthing.dflat[[1]] <- lapply(list.cl,'[[',2)	 # adjust in case

			evthing.dflon <- list()
			evthing.dflon[[1]] <- lapply(list.cl,'[[',3) # adjust in case

			# Save layers (richness map), frames (lat), frames (lon)
			fln <- paste0(output.dir_backtransformed, vec_folders[fo], "/", "div_monthly", "/", proj.type[u], "/", "Diazotrophs_speciesRichness_", mod.type[v], "_", nms.sets[d], "_paBacktransformed_taxa_corrected_mon.grd")
			writeRaster(stack(unlist(evthing.raster)), filename = fln, format = "raster", overwrite = TRUE)
			fln <- paste0(output.dir_backtransformed, vec_folders[fo], "/", "lat_monthly", "/", proj.type[u], "/", "Diazotrophs_speciesRichness_", mod.type[v], "_", nms.sets[d], "_paBacktransformed_taxa_corrected_mon.RData")
			save(evthing.dflat, file = fln)
			fln <- paste0(output.dir_backtransformed, vec_folders[fo], "/", "lon_monthly", "/", proj.type[u], "/", "Diazotrophs_speciesRichness_", mod.type[v], "_", nms.sets[d], "_paBacktransformed_taxa_corrected_mon.RData")
			save(evthing.dflon, file = fln)

			# Progress
			print(paste0( "model ", v, " strategy ", d, " monthly"))

			
		} # Close loop across dataset strategies
	} # Close loop across algorithms
} # Close loop acroos vec_folders


### =========================================================================
### 7. NCD's vs. CD's - CALCULATE RICHNESS ANNUAL MEAN - Normalized
### =========================================================================
# Loop
for(fo in seq_along(vec_folders)){
	
	for(u in seq_along(proj.type)){

		for(v in seq_along(mod.type)){

			# Loop across data sets (strategies)
			for (d in seq_along(nms.sets) ){

				# Read in list with 12 stacks: data from brick (to avoid trouble with reloading the files) and allocate into list of 12 monthly slots with pre-gathered projections
				ls.dat.mean <- get(load(  paste0(input.dir, vec_folders[fo], "/", proj.type[u], "/", mod.type[v],"_", proj.type[u], "_mean_", nms.sets[d], ".RData") ))

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

				# Prepare diversity rasters and lat frames
				lisi.taxa.mean_cd <- list()
				for (i in positions_cd){
					lisi.taxa.mean_cd[[length(lisi.taxa.mean_cd) + 1]] <- calc(stack(unlist(lapply(ls.dat.mean, '[', i))), mean, na.rm = TRUE)
					print(paste(i))
				}
				#
				lisi.taxa.mean_ncd <- list ()
				for (i in positions_ncd){
					lisi.taxa.mean_ncd[[length(lisi.taxa.mean_ncd) + 1]] <- calc(stack(unlist(lapply(ls.dat.mean, '[', i))), mean, na.rm = TRUE)
					print(paste(i))
				}

				# Calculate richness layer (sum of the taxa)
				raster.mean_cd <- calc(stack(lisi.taxa.mean_cd), sum, na.rm = TRUE) # adjust in case
				# Normalize by the number of species
				raster.mean_taxa_corrected_cd <- raster.mean_cd/length(positions_cd)
				#
				raster.mean_ncd <- calc(stack(lisi.taxa.mean_ncd), sum, na.rm = TRUE) # adjust in case
				# Normalize by the number of species
				raster.mean_taxa_corrected_ncd <- raster.mean_ncd/length(positions_ncd)

				# Merge lists
				evthing.raster <- list()
				evthing.raster[[1]] <- reclassify(raster.mean_taxa_corrected_cd, rcl = cbind(0, NA))
				names(evthing.raster[[1]]) <- 'richness_cyanobacterial_diazotrophs'
				evthing.raster[[2]] <- reclassify(raster.mean_taxa_corrected_ncd, rcl = cbind(0, NA))
				names(evthing.raster[[2]]) <- 'richness_nonCyanobacterial_diazotrophs'


				# Save layers (richness map), frames (lat), frames (lon)
				fln <- paste0(output.dir, vec_folders[fo], "/div_annually/", proj.type[u], "/", "Diazotrophs_speciesRichness_CyanosVsNonCyanos", mod.type[v], "_", proj.type[u], "_", nms.sets[d], "_annual_normalized.grd")
				writeRaster(stack(unlist(evthing.raster)), filename = fln, format = "raster", overwrite = TRUE)


				# Progress
				print(paste0(" model ", v, " strategy ", d, " annually"))
			} # Close loop across dataset strategies
		} # Close loop across algorithms
	} # Close loop across proj.type
} # Close loop across vec_folders




###==============================================================
### END
###==============================================================