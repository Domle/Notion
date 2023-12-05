## Load the compiled dataset from diazotrophs and create a raw ungridded version and a gridded one.
## We do this procedure two times, one only for presence/absence data and one time for volumetric measurements.
## To create the raw data set that will be saved, we will adjust column names so they fit for later analysis and
## round coordinates so they fit the grid required which is -179.5 to 179.5 longitude and -89.5 to 89.5 latitude.
## The gridded datasets are created by aggregating (function aggregate) using longitude, latitude and month. Discarding
## year as a factor aggregates cells (1° x 1° resolution) of different year but same month to one cell.


# Input files:  CSV file containing observations as rows and columns following Darwin Core standards 
#               (https://www.gbif.org/darwin-core) (e.g.: depth, month, taxonRank, phylum, ...).


# Output files: Two CSV files containing raw observations and gridded data with observations as 
#               rows and columns following Darwin Core standards (https://www.gbif.org/darwin-core) 
#               (e.g.: depth, month, taxonRank, phylum, ...)

# Strategy: The gridded versions are based on the removal of duplicates when working with presence absence data 
#           using !duplicated(tt1abs[ , c("scientificName", "x", "y", "year", "month").
#           1. Load data
#           2. Regrid (target grid: 1 by 1 degree resolution, -179.5, 179.5, -89.5, 89.5)
#           3. Remove duplicates
#           4. Save as csv file.

## Dominic Eriksson
## Environmental Physics Group
## ETH Zurich
## Switzerland

## dominic.eriksson@usys.ethz.ch, 25th of October 2021

### ================================================================
### Preparatory steps
### ================================================================
# Clear workspace
rm(list = ls())

## Packages
lib_vec <- c("raster", "sfsmisc", "grid", "spam", "fields", "base", "doBy", "ncdf4", "dplyr")
#
new_packages <- lib_vec[!(lib_vec %in% installed.packages()[, "Package"])] # check which package is not installed
if(length(new_packages)) {install.packages(new.packages)} # install missing packages
#
sapply(lib_vec, library, character.only = TRUE) # load packages

# Directories
wd_dat <- "/net/kryo/work/deriksson/Projects/Notion_DE/Code/2_Merge_Data/2_Output/"
wd_out <- "/net/kryo/work/deriksson/Projects/Notion_DE/Code/1_MatchUp/1_Output/"

# Help functions for aggregation to 1° monthly cell
vec.extr_mean <- function(x) {if(all(is.na(x)) == FALSE) {y <- mean(x, na.rm = T)} else {y <- NA}; return(y)} # extract vector for aggregated cells (mean)
vec.extr_sd <- function(x) {if(all(is.na(x)) == FALSE) {y <- sd(x, na.rm = T)} else {y <- NA}; return (y)} # extract vector for aggregated cells (standard deviation)
vec.extr_max <- function(x) {if(all(is.na(x)) == FALSE) {y <- max(x, na.rm = T)} else {y <- NA}; return (y)} # extract vector for aggregated cells (maximum)
vec.extr_min <- function(x) {if(all(is.na(x)) == FALSE) {y <- min(x, na.rm = T)} else {y <- NA}; return (y)} # extract vector for aggregated cells (minimum)

###
### 1. Load data
###

# Get data
tt <- read.csv(paste0(wd_dat, "Diazo_hierarchic_taxa.csv"))
dim(tt)
names(tt)
# Exploration of data types
as.character(unique(tt$occurrenceStatus))

### ======================================================================
### 2. Re-grid to common spatial resolution & 3. Remove duplicates
### ======================================================================
## Presences
# Raw
tt1pres <- tt[which(tt$occurrenceStatus%in%c("PRESENT")), ] # subset each row of dataframe for present/rare taxa
# Harmonize coordinates
tt1pres$x <- tt1pres$decimalLongitude
tt1pres$y <- tt1pres$decimalLatitude
tt1pres$taxon <- tt1pres$scientificName
#
tt1pres$x <- round(tt1pres$x + 0.5) - 0.5
tt1pres$y <- round(tt1pres$y + 0.5) - 0.5

# Check duplicates - We work with monthly climatologies, to month is our highest resolution
sum(duplicated(tt1pres[, c("scientificName", "x", "y", "year", "month")]))

# Remove duplicates in regard to scientificName, x, y, year and month
tt1presgridded <- tt1pres[!duplicated(tt1pres[, c("scientificName", "x", "y", "year", "month")]), 
    c(
    "phylum", "class", "taxon", "scientificName",
    "x", "y",
    "year", "month",
    "cellQuantity", "cellQuantityType",
    "trichomeQuantity", "trichomeQuantityType",
    "nifQuantity", "nifQuantityType",
    "mOTUcounts",
    "DNA_sequence_reads",
    "occurrenceStatus", "measurementMethod", "basisOfRecord")]

## Non-detections
# raw
tt1abs <-tt[which(tt$occurrenceStatus%in%c("ABSENT")), ]
tt1abs$occurrenceStatus <- "NOT_DETECTED"
# Harmonize coordinates
tt1abs$x <- tt1abs$decimalLongitude
tt1abs$y <- tt1abs$decimalLatitude
tt1abs$taxon <- tt1abs$scientificName
#
tt1abs$x <- round(tt1abs$x +0.5) - 0.5
tt1abs$y <- round(tt1abs$y +0.5) - 0.5

# Gridded - We work with monthly climatologies
sum(duplicated(tt1abs[ , c("scientificName", "x", "y", "year", "month")]))

tt1absgridded <- tt1abs[!duplicated(tt1abs[ , c("scientificName", "x", "y", "year", "month")]) ,  
  c(
      "phylum", "class", "taxon", "scientificName",
      "x", "y",
      "year", "month",
      "cellQuantity", "cellQuantityType",
      "trichomeQuantity", "trichomeQuantityType",
      "nifQuantity", "nifQuantityType",
      "mOTUcounts",
      "DNA_sequence_reads",
      "occurrenceStatus", "measurementMethod", "basisOfRecord")]

### ======================================================================
### 4. Save as CSV file
### ======================================================================
tt1raw <- rbind(tt1pres, tt1abs)
fln <- paste0(wd_out, "diazos_pres_abs_raw.csv")
write.csv(tt1raw, file = fln, row.names = FALSE)
tt1gri <- rbind(tt1presgridded, tt1absgridded)
fln <- paste0(wd_out, "diazos_pres_abs_gridded.csv")
write.csv(tt1gri, fln, row.names = FALSE)

### --------------------------------------------------------------------------
### END
### --------------------------------------------------------------------------
