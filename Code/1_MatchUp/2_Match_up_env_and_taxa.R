## In this script we match our diazotroph observations with the environmental variables.

## Input files: 1. CSV files containing gridded and ungridded/raw presence data with observations
##              as rows and columns following Darwin Core standards (https://www.gbif.org/darwin-core).
##              2. Environmental data saved as a grid file.

## Output files:    Two CSV files each containing diazotroph presences either gridded or ungridded/raw 
##                  and matched environmental parameters.

## Strategy:    We will use a loop to perform steps on the gridded and raw datasets
##                  
##              1. Load the data sets, environmental parameters and diazotroph observations.
##              2. Conversion of diazotroph dataframe to spatial object.
##              3. Match up presences with environmental parameters using a monthly resolution and
##              subset the dataset for oceanic conditions (remove coastal observations, depth < 200 meters
##              and salinity less than 20).
##              4. Re-order the dataframe based on number of observations (highest to low) and save as
##              CSV file containing observations as rows and columns as Darwin Core standards 
##              (https://www.gbif.org/darwin-core).

## Author:  Dominic Eriksson
##          Environmental Physics Group, UP
##          ETH Zurich
##          Switzerland

# dominic.eriksson@usys.ethz.ch, 17th of November 2021 ------------------------------



##--------------------------------------------------------------------------
## Preparatory steps
##--------------------------------------------------------------------------

# Clear workspace
rm(list = ls())

# Load libraries
lib_vec <- c("raster", "sfsmisc", "grid", "spam", "fields", "base", "doBy", "ncdf4", "grid")
if(length(new_packages)) {install.packages(new.packages)} # install missing packages
sapply(lib_vec, library, character.only = TRUE) # load packages

# Directories
wd_dat      <- "/net/kryo/work/deriksson/Projects/Notion_DE/Code/3_Match_up/1_Output/"
wd_pred_var <- "/net/kryo/work/deriksson/Projects/Notion_DE/Data/Predictor_variables/Environmental_data/"
wd_out      <- "/net/kryo/work/deriksson/Projects/Notion_DE/Code/1_MatchUp/2_Output/"


##--------------------------------------------------------------------------
## 1. Load data
##--------------------------------------------------------------------------

# Get environmental data
prj.stack <- brick(paste0(wd_pred_var, "VarSet07.2016.grd"))
dim(prj.stack) # 180 360 576 --> monthly resolution 1°x 1°
#
nelem <- nlayers(prj.stack)/12 # 48
nm <- gsub("\\.1", "", names(prj.stack[[1:nelem]])) # correct the variable names
env.stack <- list()
from <- seq(1, nlayers(prj.stack), by = nelem)
for(q in 1:12){
  pre.env.stack <- stack(prj.stack[[ from[q]: (from[q] + nelem - 1)]])
  names(pre.env.stack) <- nm;
  env.stack[[q]] <- pre.env.stack # we created a list. 12 raster objects (from January to December) containing the env. variables
}
names(env.stack) <- c("JAN", "FEB", "MAR", "ARP", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC")
class(env.stack) # list
class(env.stack[[1]]) # raster

## Load data on total presences of diazotrophs
## Grided versions
df_gridded <- read.csv(paste0(wd_dat, "diazos_pres_abs_gridded.csv"))
# Raw versions
df_raw <- read.csv(paste0(wd_dat, "diazos_pres_abs_raw.csv"))

## Add environmental variables using a for loop based on dataset (gridded vs. raw)
# Preparatory
vec.dat <- c(

    # Gridded versions
    "pres_abs",
    
    # Raw versions
    "pres_abs")

### LOOP across data-type versions (presence data versus count data using consistent methodology)
for(s in seq_along(vec.dat)){
  
    ## Get data
    if(s == 1){dat <- df_gridded}else{dat <- df_raw}


    ## 2. Convert dataframe to spatial object -------------------------------------------------
    # Formatting for spatial objects
    llCRS <- CRS("+proj=longlat +ellps=WGS84") # add projection type
    crd.index <- which(colnames(dat)%in%c(c("x", "y")))
    
    ## Convert data points to SpatialDataFrame
    spdf <- SpatialPointsDataFrame(dat[, crd.index], dat[, -crd.index], proj4string = llCRS)
    
    ### 3. Monthly match up and apply open oceanic conditions ---------------------------------
    ## Loop across months: mask the observations and prepare list to extract points, column-bind with mask, keep rows (points) with mask equal 1 (ocean condition)
    mups <- list()
    for(k in 1:12){
        
        # 1. Define monthly subset, for monthly match-up with mask
        if(sum(spdf$month == k, na.rm = TRUE) > 0){ss <- subset(spdf, spdf$month == k)} else {next}
        
        # 2. Extract data points overlapping with "mask"
        extr <- raster::extract(x = env.stack[[k]], y = ss, method = "simple") # returns the values of a raster object (msk) for cells in which points fall
        extr <- as.data.frame(extr) # the mask values are converted to data frame
        ss.new <- cbind(as.data.frame(ss), extr)[which(extr$Mask!=0), ] # CORE OPERATION: SELECT ROWS/DATA for which mask == 1
        
        # 3. Store monthly specific ss.new in list (for later row-bind)
        mups[[length(mups) + 1]] <- ss.new
        print(paste(k))
    } # end of loop across months

    # Merge match-ups
    tot <- do.call("rbind", mups)
    
    # Usefule gridded 1° by 1° monthly data records for modelling
    dim(tot); class(tot) # Note that any information on the measurement method is lost at this point, and coastal data are excluded

    ## 4. Reorder and save dataframe ------------------------------------------------------------
    # Overview of taxa
    dat.spl <- split(tot, f = tot$scientificName)
    df.taxa.obs <- data.frame("taxon" = names(unlist(lapply(dat.spl, nrow))), "obs" = as.numeric(unlist(lapply(dat.spl, nrow))))
    df.taxa.obs <- df.taxa.obs[order(df.taxa.obs$obs, decreasing = T), ]
    rownames(df.taxa.obs) <- 1:nrow(df.taxa.obs)
    df.taxa.obs <- df.taxa.obs[which(df.taxa.obs$obs > 0), ]
    df.taxa.obs <- df.taxa.obs[order(df.taxa.obs$taxon), ]
    rownames(df.taxa.obs) <- 1:nrow(df.taxa.obs)
    df.taxa.obs

    ## Order data
    # Define taxa as character
    tot$scientificName <- as.character(tot$scientificName)

    # Re-order the data along taxa
    lisi.taxa <- list()
    # Loop across taxa: we keep taxa == NA such as picoeukaryotes for the moment
    for(i in 1:nrow(df.taxa.obs)){
        
        # 1. get data of taxon
        lisi.taxa[[length(lisi.taxa) + 1]] <- tot[which(tot$scientificName == df.taxa.obs$taxon[i]), ]
        names(lisi.taxa[[1]])
        
        # Re-merge groups
        tot.new <- do.call("rbind", lisi.taxa)
        rownames(tot.new) <- 1:nrow(tot.new)

        # Add cell ID
        tot.new$month <- as.numeric(as.character(tot.new$month)) # change months from integer to numeric
        
        # Preparatory: as.character(problem: function paste, collapse = "_" below introduces some spaces into the IDs if the components x, y and motnh are numeric; so I try it with as character here)
        tot.new$month <- as.character(tot.new$month)
        tot.new$x <- as.character(tot.new$x)
        tot.new$y <- as.character(tot.new$y)
        tot.new$ID.cell <- apply(tot.new, 1, FUN = function(x){x1 <- paste(x[c("x", "y", "month")], sep = "", collapse = "_"); return(x1)}) # Create ID of monthly sampling cell
        tot.new$month <- as.numeric(tot.new$month) # change back to numeric
        tot.new$x <- as.numeric(tot.new$x)
        tot.new$y <- as.numeric(tot.new$y)

    } # close loop across taxa

    # Save as CSV file
    if(s == 1) {fln <- paste0(wd_out, "Obs_gridded_diazo_",vec.dat[s],".csv")}
    if(s == 2) {fln <- paste0(wd_out, "Obs_raw_diazo_",vec.dat[s],".csv")}

    write.csv(tot.new, file = fln, row.names = FALSE)
} # Close loop over data strategies

###==============================================================
### END
###==============================================================