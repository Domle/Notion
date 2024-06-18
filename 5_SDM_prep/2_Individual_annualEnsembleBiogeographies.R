# In this script we compute individual ensemble biogeographies of each taxa on a monthly and also
# annual temporal resolutions. NOTE: TAKE CARE OF MAJORITY VOTE ALSO

# Input files: 
#       1. 

# Output files:  
#       2.  We create raster stacks of length 92, with ensemble mean, standard deviation and the output of the 90 ensemble members.


# Author: 	Dominic Eriksso
#			Environmental Physics Group, UP
#			ETH, Zurich
#			Switzerland

## deriksson@ethz.ch; 18th of June

### =======================================================================
### Calculate ensemble mean, sd annual
### =======================================================================
## We calculate the annual global distribution for each model and average at the end to show
## uncertainties

# Directory
wd <- "/home/deriksson/Projects/Notion_DE/Code/7_SDM_prep/1_Output_v2/"
wd_out <- "/home/deriksson/Projects/Notion_DE/Code/Ensembles/Annual_individual_ensemble_biogeographies/"

# Projection vector
proj <- c("hsi", "pa")

for(p in seq_along(proj)){

  # Get filenames
  filenames <- list.files( paste0(wd, proj[p], "/" ) )
  # Loop through filenames
  for(f in seq_along(filenames)){

    # Print progress
    print( filenames[f])

    # Load data
    l <- readRDS( paste0(wd, proj[p], "/", filenames[f]) )
    if(length(l) > 0){
      r <- raster()
      for(ll in seq_along(l)){
        rr <- l[[ll]] 
        rr <- calc(rr, mean, na.rm = TRUE)
        r <- stack(r, rr)
      } # end of loop across list length

      # Compute annual ensemble
      raster_mean <- calc(r, mean, na.rm = TRUE)
      raster_sd <- calc(r, sd, na.rm = TRUE)
      names(r) <- (paste0("SDM_", 1:nlayers(r)))

      # Save grid
      raster_stack <- stack(raster_mean, raster_sd, r)
      names(raster_stack)[1:2] <- c("Annual ensemble mean", "Annual ensemble standard deviation")
      fln <- paste0( wd_out, proj[p], "/", gsub( "_all_models.rds", "", filenames[f] ), "_annual_ensemble_mean_sd_individualEnsembleMembers.grd")
      writeRaster(
        raster_stack, 
        fln,
        overwrite = TRUE
      )
    } # end of if statement, if list of length zero
  } # end of loop across filenames
} # end of loop across projection types
