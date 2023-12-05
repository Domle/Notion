# In this script we compute individual ensemble biogeographies of each taxa on a monthly and also
# annual temporal resolutions. NOTE: TAKE CARE OF MAJORITY VOTE ALSO

# Input files:
#       1. R-files containing lists of all month and each taxon as raster objects.
#       2. Evaluation table with TSS scores of each taxon and run.
#       3. A table with taxonomic information on all diazotrophs used in the study

# Output files:
#       1.  Lists containing raster elements that store at least mean and standard deviation from
#           from each taxon either on annual or monthly temporal resolutions.


# Strategy: We use the evalution table that contains the TSS score for each taxon and
#           model run, to subset the list elements, and re-stack the object so that the calc-function
#           from the raster object can be used to compute annual means. For the monthly means we 
#           we will extract each taxon and loop through each month separately saving the rasters 
#           of each model as in one stack and apply the calc function to that stack. We will save a list of 
#           length twelve containing rasters for one species and mean, standard deviation, minimum and maximum stats.
#           For the majority vote option we want to end up with either a species being present or absent when using an 
#           ensemble. Each location will be assigned a presence or absence based on a threshold of more than half of the 
#           ensemble members agreement.

# Author: 	Dominic Eriksso
#			Environmental Physics Group, UP
#			ETH, Zurich
#			Switzerland

# deriksson@ethz.ch 16th of June, 2023 ----------------------------------------------------------



### ==============================================================================
### Preparation
### ==============================================================================
rm(list = ls())

# Libraries
library(raster)

# Load functions
source("/home/deriksson/Projects/Notion_DE/Code/Functions/Righetti/common_functions.R") # Load functions to derive latitudinal zonal mean pattern
source("/home/deriksson/Projects/Notion_DE/Code/Functions/Dominic/call_libraries.R")


## Directories
wd_in <- "/net/kryo/work/deriksson/Projects/Notion_DE/Code/7_SDM_prep/2a_Output/"
wd_out <- "/net/kryo/work/deriksson/Projects/Notion_DE/Code/7_SDM_prep/3_Output/"


# Vector for looping
ensemble_pa <- c("Ensemble_mean", "PA_majority_vote")
p <- 1
#
vec_folders <- c("Total_dataset", "SequenceBased_dataset", "MicroscopyBased_dataset")
f <- 1
#
proj.type <- c("prb", "pa")

# Load tax table
tax_table <- read.csv('/home/deriksson/Projects/Notion_DE/Data/commons/Taxon_table.csv', sep = ';')


### ===============================================================
### Annual ensemble biogeographies
###================================================================

## Loop 
for(u in seq_along(proj.type)){ #  not used yet

    # Get filenames to be included in the ensmebles
    filenames <- list.files( paste0( wd_in, ensemble_pa[1], "/", vec_folders[f], "/", proj.type[u]) )
    filenames <- grep(".RData", filenames, value = TRUE)
    filenames <- grep("mean", filenames, value = TRUE)

    # Loop over filenames
    l_fi <- list()
    l_pos <- list()
    for(fi in seq_along(filenames)){    # NOTE: Not all species have been modelled successfully by all 18 ensemble members. The mean function accounts for that
                                        # and only calculates the mean across those ensemble members where species/taxon has been modeled successfully.
                
        # Load biogeographies
        l_rich <- get(load( paste0( wd_in, ensemble_pa[1], "/", vec_folders[f], "/", proj.type[u], "/", filenames[fi]) ))
                
        # Load which taxa pass quality control
        # Format filename to load evaluation table
        if(proj.type[u] == "prb"){
            filenames_eval <- gsub("prb_mean_", "eval_", filenames[fi])
        }else{
            filenames_eval <- gsub("pa_mean_", "eval_", filenames[fi])
        }
        filenames_eval <- gsub(".RData", ".csv", filenames_eval)
        # Load evaluation statistics 
        eval <- read.csv( paste0("/home/deriksson/Projects/Notion_DE/Code/6_SDM_fit/3_Output/", vec_folders[f],  "/", filenames_eval) )
        eval <- eval[which(!is.na(eval$"tss.xval4")), ]
        eval <- eval[which(eval$"tss.xval4" > 0.3), ]
        positions <- unique(eval$taxon)
       
        # We create a list of raster stacks containing one species and all month as layers
        lisi.taxa.mean <- list()
        # Loop through each taxon and compute annual average
        for (i in positions){
        
            # Sanity check
            stack(unlist( lapply(l_rich, '[', i)  )) 

            # Extract taxon and month and stack
            lisi.taxa.mean[[length(lisi.taxa.mean) + 1]] <- stack(unlist( lapply(l_rich, '[', i)  ))

            # Print progress
            print(paste(i))

        } # close loop across positions

        # Name list elements, as those are the taxa
        names(lisi.taxa.mean) <- positions

        # Print progress
        print(paste0("Finishing file: ", filenames[fi]))

        # Store result
        l_fi[[fi]] <- lisi.taxa.mean

    } # close loop across filenames

    # Now we need to create a vector with all taxa present in each model
    vec_taxa <- as.character()
    for(i in seq_along(l_fi) ){
        n <- names(l_fi[[i]])
        vec_taxa <- c(vec_taxa, n)
    }
    # Collapse names
    positions <- unique(vec_taxa)

    # Calculate annual ensemble by averaging across
    lisi_ensemble_mean <- list ()
    lisi_ensemble_sd <- list()
    for (i in positions){
                
        # Sanity check: 12 (months) * number of successful ensemble members (SDMs)
        nlayers(stack(unlist(lapply(l_fi, '[', i))))
                
        # Compute annual ensemble mean --------------------------------------------------------------
        lisi_ensemble_mean[[length(lisi_ensemble_mean) + 1]] <- calc(stack(unlist(lapply(l_fi, '[', i))), mean, na.rm = TRUE)
        names(lisi_ensemble_mean)[length(lisi_ensemble_mean)] <- i
        # Ensemble SD -------------------------------------------------------------------------------
        lisi_ensemble_sd[[length(lisi_ensemble_sd) + 1]] <- calc(stack(unlist(lapply(l_fi, '[', i))), sd, na.rm = TRUE)
        names(lisi_ensemble_sd)[length(lisi_ensemble_sd)] <- i
        # Print progress
        print(paste(i))
    } # close loop across positions

    # Save object
    fname_mean <- paste0( wd_out, ensemble_pa[1], "/", vec_folders[f], "/", proj.type[u], "/Annual_ensembles/Annual_HSI_Ensemble_mean_AllModels_IndividualTaxa.R") 
    fname_sd <- paste0( wd_out, ensemble_pa[1], "/", vec_folders[f], "/", proj.type[u], "/Annual_ensembles/Annual_HSI_Ensemble_sd_AllModels_IndividualTaxa.R") 
    save(lisi_ensemble_mean, file = fname_mean)
    save(lisi_ensemble_sd, file = fname_sd)

} # end loop across proj.type


### ===============================================================
### We compute ensemble for each month separately
###================================================================
# Set proj.type index
u <- 1

# Get filenames to be included in the ensmebles
filenames <- list.files( paste0( wd_in, ensemble_pa[1], "/", vec_folders[f], "/", proj.type[u]) )
filenames <- grep(".RData", filenames, value = TRUE)
filenames <- grep("mean", filenames, value = TRUE)


# In total we have 27 species that potentially were able to be modeled. Retrieve a vector containing all taxa
# available. 
species <- get(load( paste0( wd_in, ensemble_pa[1], "/", vec_folders[f], "/", proj.type[u], "/", filenames[1]) ))
species <- names(species[[1]])


for(s in seq_along(species)){
       
    # We save each month as a list object
    stack_final <- list()
    for(m in 1:12){


        # Create stacks to save results
        stack_monthly <- stack()
        stack_stats <- stack()
        # Strategy: We open each filename while keeping the month index (m) constant. Then compile a new raster stack containing all files from one specific month 
        # where we compute the average, sd, min and max.
        for(fi in seq_along(filenames)){

            # Open list of rasters
            l <- get(load( paste0( wd_in, ensemble_pa[1], "/", vec_folders[f], "/", proj.type[u], "/", filenames[fi]) ))
            # Extract data
            sp_name <- names( l[[m]] [s] )
            r <- l[[m]][[ sp_name ]]

            # Print progress
            print(paste0( 'Working on species ', sp_name, ', month ', m ))

            # Stack data
            stack_stats <- addLayer(stack_stats, r)
        } # close loop across filenames

        # Calculate monthly ensemble stats
        stack_monthly <- addLayer(stack_monthly, calc(stack_stats, fun = mean, na.rm = TRUE) )
        stack_monthly <- addLayer(stack_monthly, calc(stack_stats, fun = sd, na.rm = TRUE) )
        stack_monthly <- addLayer(stack_monthly, calc(stack_stats, fun = max, na.rm = TRUE) )
        stack_monthly <- addLayer(stack_monthly, calc(stack_stats, fun = min, na.rm = TRUE) )
        # Name layers
        names(stack_monthly) <- c('mean', 'sd', 'max', 'min')

        # Save in list
        stack_final[[m]] <- stack_monthly

        # Remove 
        rm(stack_monthly)
        rm(stack_stats)
    } # close loop across month (m)

    # Save 
    saveRDS(
        object = stack_final, 
        file = paste0(wd_out, ensemble_pa[1], "/", vec_folders[f], "/", proj.type[u],'/Monthly_ensembles/', sp_name, '_ensembleStatsMonthly.R' )
        )
} # close loop across species




### ===============================================================
### Annual ensemble biogeographies - Majority vote
###================================================================
# Set proj.type index to pa (U = 2)
u <- 2

# Get filenames to be included in the ensmebles
filenames <- list.files( paste0( wd_in, ensemble_pa[2], "/", vec_folders[f], "/", proj.type[u]) )
filenames <- grep(".RData", filenames, value = TRUE)
filenames <- grep("mean", filenames, value = TRUE)

# Loop over filenames
l_fi <- list()
l_pos <- list()
    for(fi in seq_along(filenames)){    # NOTE: Not all species have been modelled successfully by all 18 ensemble members. The mean function accounts for that
                                        # and only calculates the mean across those ensemble members where species/taxon has been modeled successfully.
                
    # Load biogeographies
    l_rich <- get(load( paste0( wd_in, ensemble_pa[2], "/", vec_folders[f], "/", proj.type[u], "/", filenames[fi]) ))
                
    # Load which taxa pass quality control
    # Format filename to load evaluation table
    filenames_eval <- gsub("pa_paBacktransformed_mean_", "eval_", filenames[fi])
    filenames_eval <- gsub(".RData", ".csv", filenames_eval)
        
    # Load evaluation statistics 
    eval <- read.csv( paste0("/home/deriksson/Projects/Notion_DE/Code/6_SDM_fit/3_Output/", vec_folders[f],  "/", filenames_eval) )
    eval <- eval[which(!is.na(eval$"tss.xval4")), ]
    eval <- eval[which(eval$"tss.xval4" > 0.3), ]
    positions <- unique(eval$taxon)
       
    # We create a list of raster stacks containing one species and all month as layers
    lisi.taxa.mean <- list()
    # Loop through each taxon and compute annual average
    for (i in positions){
        
        # Sanity check
        stack(unlist( lapply(l_rich, '[', i)  )) 

        # Extract taxon and month and stack
        lisi.taxa.mean[[length(lisi.taxa.mean) + 1]] <- stack(unlist( lapply(l_rich, '[', i)  ))

        # Print progress
        print(paste(i))

    } # close loop across positions

    # Name list elements, as those are the taxa
    names(lisi.taxa.mean) <- positions

    # Print progress
    print(paste0("Finishing file: ", filenames[fi]))

    # Store result
    l_fi[[fi]] <- lisi.taxa.mean

} # close loop across filenames

# Now we need to create a vector with all taxa present in each model
vec_taxa <- as.character()
for(i in seq_along(l_fi) ){
    n <- names(l_fi[[i]])
    vec_taxa <- c(vec_taxa, n)
}
# Collapse names
positions <- unique(vec_taxa)

# Calculate annual ensemble by averaging across
lisi_ensemble_mean <- list ()
for (i in positions){
                
    # Sanity check: 12 (months) * number of successful ensemble members (SDMs)
    nlayers(stack(unlist(lapply(l_fi, '[', i))))
                
    ## Compute annual ensemble mean --------------------------------------------------------------
    # NOTE: This ensemble contains all twelve month and all 18 models. We will assign a presence to the annual distribution estimate if
    # the average is equal or above 0.5
    lisi_ensemble_mean[[length(lisi_ensemble_mean) + 1]] <- calc(stack(unlist(lapply(l_fi, '[', i))), mean, na.rm = TRUE)
    # Reclassify to presence or absence 
    lisi_ensemble_mean[[length(lisi_ensemble_mean)]] <- reclassify( lisi_ensemble_mean[[length(lisi_ensemble_mean)]], cbind(0.5, Inf, 1))
    lisi_ensemble_mean[[length(lisi_ensemble_mean)]] <- reclassify( lisi_ensemble_mean[[length(lisi_ensemble_mean)]], cbind(-Inf, 0.5, 0))
    names(lisi_ensemble_mean)[length(lisi_ensemble_mean)] <- i

    # Print progress
    print(paste(i))
} # close loop across positions

# Save object
fname_mean <- paste0( wd_out, ensemble_pa[2], "/", vec_folders[f], "/", proj.type[u], "/Annual_ensembles/Annual_MajorityVote_Ensemble_mean_AllModels_IndividualTaxa.R") 
save(lisi_ensemble_mean, file = fname_mean)



### ===============================================================
### We compute ensemble for each month separately - Majority vote
###================================================================
# Set proj.type index to pa (U = 2)
u <- 2

# Get filenames to be included in the ensmebles
filenames <- list.files( paste0( wd_in, ensemble_pa[2], "/", vec_folders[f], "/", proj.type[u]) )
filenames <- grep(".RData", filenames, value = TRUE)


# In total we have 27 species that potentially were able to be modeled. Retrieve a vector containing all taxa
# available. 
species <- get(load( paste0( wd_in, ensemble_pa[2], "/", vec_folders[f], "/", proj.type[u], "/", filenames[1]) ))
species <- names(species[[1]])


for(s in seq_along(species)){
       
    # We save each month as a list object
    stack_final <- list()
    for(m in 1:12){


        # Create stacks to save results
        stack_monthly <- stack()
        stack_stats <- stack()
        # Strategy: We open each filename while keeping the month index (m) constant. Then compile a new raster stack containing all files from one specific month 
        # where we compute the average, sd, min and max.
        for(fi in seq_along(filenames)){

            # Open list of rasters
            l <- get(load( paste0( wd_in, ensemble_pa[2], "/", vec_folders[f], "/", proj.type[u], "/", filenames[fi]) ))
            # Extract data
            sp_name <- names( l[[m]] [s] )
            r <- l[[m]][[ sp_name ]]

            # Print progress
            print(paste0( 'Working on species ', sp_name, ', month ', m ))

            # Stack data
            stack_stats <- addLayer(stack_stats, r)
        } # close loop across filenames

        # Calculate monthly ensemble mean
        r <- calc(stack_stats, fun = mean, na.rm = TRUE)
        # Reclassify values, everything above 0.5 becomes a presence, below and absence
        r <- reclassify(r, cbind(0.5, Inf, 1))
        r <- reclassify(r, cbind(-Inf, 0.5, 0))
        # Add monthly reclassified ensemble (majority vote)
        stack_monthly <- addLayer(stack_monthly,  r)

        # Name layers
        names(stack_monthly) <- c('majority_vote')

        # Save in list
        stack_final[[m]] <- stack_monthly

        # Remove 
        rm(stack_monthly)
        rm(stack_stats)
    } # close loop across month (m)

    # Save 
    saveRDS(
        object = stack_final, 
        file = paste0(wd_out, ensemble_pa[2], "/", vec_folders[f], "/", proj.type[u],'/Monthly_ensembles/', sp_name, '_MajorityVote_Monthly.R' )
        )
} # close loop across species



###==============================================================
### END
###==============================================================