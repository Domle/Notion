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



### =========================================================================
### Preparatory steps
### =========================================================================

# Clear workspace
rm(list = ls())

# Packages
lib_vec <- c("raster", "sp")

# Install and load packages
lapply(lib_vec, function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
})

# Directories
wd_dat <- "/home/deriksson/Projects/Notion_DE/Code/GitHub_PublicRepo/1_MatchUp/1_Output/"
wd_pred_var <- "/net/kryo/work/deriksson/Projects/Notion_DE/Data/Predictor_variables/Environmental_data/"
wd_out <- "/home/deriksson/Projects/Notion_DE/Code/GitHub_PublicRepo/1_MatchUp/2_Output/"

# Create the output directory if it doesn't exist
if (!dir.exists(wd_out)) {
  dir.create(wd_out, recursive = TRUE)
}

### =========================================================================
### 1. Load data
### =========================================================================

# Get environmental data
prj.stack <- raster::brick(paste0(wd_pred_var, "VarSet07.2016.grd"))
dim(prj.stack) # Check dimensions: 180 360 576 --> monthly resolution 1°x 1°

# Correct the variable names
nelem <- raster::nlayers(prj.stack) / 12
nm <- gsub("\\.1", "", names(prj.stack[[1:nelem]]))

# Create a list of 12 raster objects (one for each month) containing the environmental variables
env.stack <- list()
from <- seq(1, raster::nlayers(prj.stack), by = nelem)
for (q in 1:12) {
  pre.env.stack <- raster::stack(prj.stack[[from[q]:(from[q] + nelem - 1)]])
  names(pre.env.stack) <- nm
  env.stack[[q]] <- pre.env.stack
}
names(env.stack) <- c("JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC")

# Check classes
class(env.stack) # list
class(env.stack[[1]]) # raster

## Load data on total presences of diazotrophs
# Total methods
df_gridded <- read.csv(paste0(wd_dat, "diazos_pres_abs_gridded.csv"))
df_raw <- read.csv(paste0(wd_dat, "diazos_pres_abs_raw.csv"))
# Microscopy based observations
df_gridded_microscopy <- read.csv(paste0(wd_dat, "diazos_microscopyBased_pres_abs_gridded.csv"))
df_raw_microscopy <- read.csv(paste0(wd_dat, "diazos_microscopyBased_pres_abs_raw.csv"))
# Sequence based observations
df_gridded_sequence <- read.csv(paste0(wd_dat, "diazos_sequenceBased_pres_abs_gridded.csv"))
df_raw_sequence <- read.csv(paste0(wd_dat, "diazos_sequenceBased_pres_abs_raw.csv"))

# Vector for data types, because we create one gridded and one raw version
vec.dat <- c("pres_abs", "pres_abs")

### =========================================================================
### 2. Add environmental variables using a loop based on dataset (gridded vs. raw)
### =========================================================================

# Total methods ------------------------------------------------------------------------
for (s in seq_along(vec.dat)) {
  
  # Get data
  dat <- if (s == 1) df_gridded else df_raw
  
  # Convert dataframe to spatial object
  llCRS <- sp::CRS("+proj=longlat +ellps=WGS84")
  crd.index <- which(colnames(dat) %in% c("x", "y"))
  spdf <- sp::SpatialPointsDataFrame(dat[, crd.index], dat[, -crd.index], proj4string = llCRS)
  
  # Monthly match up and apply open oceanic conditions
  mups <- list()
  for (k in 1:12) {
    if (sum(spdf$month == k, na.rm = TRUE) > 0) {
      ss <- subset(spdf, spdf$month == k)
    } else {
      next
    }
    
    extr <- raster::extract(x = env.stack[[k]], y = ss, method = "simple")
    extr <- as.data.frame(extr)
    ss.new <- cbind(as.data.frame(ss), extr)[which(extr$Mask != 0), ]
    mups[[length(mups) + 1]] <- ss.new
    print(paste(k))
  }
  
  # Merge match-ups
  tot <- do.call("rbind", mups)
  dim(tot)
  class(tot)
  
  # Reorder and save dataframe
  dat.spl <- split(tot, f = tot$scientificName)
  df.taxa.obs <- data.frame("taxon" = names(unlist(lapply(dat.spl, nrow))), "obs" = as.numeric(unlist(lapply(dat.spl, nrow))))
  df.taxa.obs <- df.taxa.obs[order(df.taxa.obs$obs, decreasing = TRUE), ]
  rownames(df.taxa.obs) <- 1:nrow(df.taxa.obs)
  df.taxa.obs <- df.taxa.obs[df.taxa.obs$obs > 0, ]
  df.taxa.obs <- df.taxa.obs[order(df.taxa.obs$taxon), ]
  rownames(df.taxa.obs) <- 1:nrow(df.taxa.obs)
  df.taxa.obs
  
  tot$scientificName <- as.character(tot$scientificName)
  
  lisi.taxa <- list()
  for (i in 1:nrow(df.taxa.obs)) {
    lisi.taxa[[length(lisi.taxa) + 1]] <- tot[tot$scientificName == df.taxa.obs$taxon[i], ]
    tot.new <- do.call("rbind", lisi.taxa)
    rownames(tot.new) <- 1:nrow(tot.new)
    
    tot.new$month <- as.numeric(as.character(tot.new$month))
    tot.new$month <- as.character(tot.new$month)
    tot.new$x <- as.character(tot.new$x)
    tot.new$y <- as.character(tot.new$y)
    tot.new$ID.cell <- apply(tot.new, 1, FUN = function(x) {
      paste(x[c("x", "y", "month")], collapse = "_")
    })
    tot.new$month <- as.numeric(tot.new$month)
    tot.new$x <- as.numeric(tot.new$x)
    tot.new$y <- as.numeric(tot.new$y)
  }
  
  # Save as CSV file
  fln <- if (s == 1) paste0(wd_out, "Obs_gridded_diazo_", vec.dat[s], ".csv") else paste0(wd_out, "Obs_raw_diazo_", vec.dat[s], ".csv")
  write.csv(tot.new, file = fln, row.names = FALSE)
}


# Microscopy based methods ------------------------------------------------------------------------
for (s in seq_along(vec.dat)) {
  
  # Get data
  dat <- if (s == 1) df_gridded_microscopy else df_raw_microscopy
  
  # Convert dataframe to spatial object
  llCRS <- sp::CRS("+proj=longlat +ellps=WGS84")
  crd.index <- which(colnames(dat) %in% c("x", "y"))
  spdf <- sp::SpatialPointsDataFrame(dat[, crd.index], dat[, -crd.index], proj4string = llCRS)
  
  # Monthly match up and apply open oceanic conditions
  mups <- list()
  for (k in 1:12) {
    if (sum(spdf$month == k, na.rm = TRUE) > 0) {
      ss <- subset(spdf, spdf$month == k)
    } else {
      next
    }
    
    extr <- raster::extract(x = env.stack[[k]], y = ss, method = "simple")
    extr <- as.data.frame(extr)
    ss.new <- cbind(as.data.frame(ss), extr)[which(extr$Mask != 0), ]
    mups[[length(mups) + 1]] <- ss.new
    print(paste(k))
  }
  
  # Merge match-ups
  tot <- do.call("rbind", mups)
  dim(tot)
  class(tot)
  
  # Reorder and save dataframe
  dat.spl <- split(tot, f = tot$scientificName)
  df.taxa.obs <- data.frame("taxon" = names(unlist(lapply(dat.spl, nrow))), "obs" = as.numeric(unlist(lapply(dat.spl, nrow))))
  df.taxa.obs <- df.taxa.obs[order(df.taxa.obs$obs, decreasing = TRUE), ]
  rownames(df.taxa.obs) <- 1:nrow(df.taxa.obs)
  df.taxa.obs <- df.taxa.obs[df.taxa.obs$obs > 0, ]
  df.taxa.obs <- df.taxa.obs[order(df.taxa.obs$taxon), ]
  rownames(df.taxa.obs) <- 1:nrow(df.taxa.obs)
  df.taxa.obs
  
  tot$scientificName <- as.character(tot$scientificName)
  
  lisi.taxa <- list()
  for (i in 1:nrow(df.taxa.obs)) {
    lisi.taxa[[length(lisi.taxa) + 1]] <- tot[tot$scientificName == df.taxa.obs$taxon[i], ]
    tot.new <- do.call("rbind", lisi.taxa)
    rownames(tot.new) <- 1:nrow(tot.new)
    
    tot.new$month <- as.numeric(as.character(tot.new$month))
    tot.new$month <- as.character(tot.new$month)
    tot.new$x <- as.character(tot.new$x)
    tot.new$y <- as.character(tot.new$y)
    tot.new$ID.cell <- apply(tot.new, 1, FUN = function(x) {
      paste(x[c("x", "y", "month")], collapse = "_")
    })
    tot.new$month <- as.numeric(tot.new$month)
    tot.new$x <- as.numeric(tot.new$x)
    tot.new$y <- as.numeric(tot.new$y)
  }
  
  # Save as CSV file
  fln <- if (s == 1) paste0(wd_out, "Obs_gridded_diazo_microscopyBased_", vec.dat[s], ".csv") else paste0(wd_out, "Obs_raw_diazo_microscopyBased_", vec.dat[s], ".csv")
  write.csv(tot.new, file = fln, row.names = FALSE)
}


# Sequence based methods ------------------------------------------------------------------------
for (s in seq_along(vec.dat)) {
  
  # Get data
  dat <- if (s == 1) df_gridded_sequence else df_raw_sequence
  
  # Convert dataframe to spatial object
  llCRS <- sp::CRS("+proj=longlat +ellps=WGS84")
  crd.index <- which(colnames(dat) %in% c("x", "y"))
  spdf <- sp::SpatialPointsDataFrame(dat[, crd.index], dat[, -crd.index], proj4string = llCRS)
  
  # Monthly match up and apply open oceanic conditions
  mups <- list()
  for (k in 1:12) {
    if (sum(spdf$month == k, na.rm = TRUE) > 0) {
      ss <- subset(spdf, spdf$month == k)
    } else {
      next
    }
    
    extr <- raster::extract(x = env.stack[[k]], y = ss, method = "simple")
    extr <- as.data.frame(extr)
    ss.new <- cbind(as.data.frame(ss), extr)[which(extr$Mask != 0), ]
    mups[[length(mups) + 1]] <- ss.new
    print(paste(k))
  }
  
  # Merge match-ups
  tot <- do.call("rbind", mups)
  dim(tot)
  class(tot)
  
  # Reorder and save dataframe
  dat.spl <- split(tot, f = tot$scientificName)
  df.taxa.obs <- data.frame("taxon" = names(unlist(lapply(dat.spl, nrow))), "obs" = as.numeric(unlist(lapply(dat.spl, nrow))))
  df.taxa.obs <- df.taxa.obs[order(df.taxa.obs$obs, decreasing = TRUE), ]
  rownames(df.taxa.obs) <- 1:nrow(df.taxa.obs)
  df.taxa.obs <- df.taxa.obs[df.taxa.obs$obs > 0, ]
  df.taxa.obs <- df.taxa.obs[order(df.taxa.obs$taxon), ]
  rownames(df.taxa.obs) <- 1:nrow(df.taxa.obs)
  df.taxa.obs
  
  tot$scientificName <- as.character(tot$scientificName)
  
  lisi.taxa <- list()
  for (i in 1:nrow(df.taxa.obs)) {
    lisi.taxa[[length(lisi.taxa) + 1]] <- tot[tot$scientificName == df.taxa.obs$taxon[i], ]
    tot.new <- do.call("rbind", lisi.taxa)
    rownames(tot.new) <- 1:nrow(tot.new)
    
    tot.new$month <- as.numeric(as.character(tot.new$month))
    tot.new$month <- as.character(tot.new$month)
    tot.new$x <- as.character(tot.new$x)
    tot.new$y <- as.character(tot.new$y)
    tot.new$ID.cell <- apply(tot.new, 1, FUN = function(x) {
      paste(x[c("x", "y", "month")], collapse = "_")
    })
    tot.new$month <- as.numeric(tot.new$month)
    tot.new$x <- as.numeric(tot.new$x)
    tot.new$y <- as.numeric(tot.new$y)
  }
  
  # Save as CSV file
  fln <- if (s == 1) paste0(wd_out, "Obs_gridded_diazo_sequenceBased_", vec.dat[s], ".csv") else paste0(wd_out, "Obs_raw_diazo_sequenceBased_", vec.dat[s], ".csv")
  write.csv(tot.new, file = fln, row.names = FALSE)
}


###==============================================================
### END
###==============================================================
