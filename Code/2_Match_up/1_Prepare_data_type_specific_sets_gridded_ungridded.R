## Load the compiled dataset from diazotrophs and create a raw ungridded version and a gridded one.
## We do this procedure two times, one only for presence/absence data and one time for volumetric measurements.
## To create the raw data set that will be saved, we will adjust column names so they fit for later analysis and
## round coordinates so they fit the grid required which is -179.5 to 179.5 longitude and -89.5 to 89.5 latitude.
## The gridded datasets are created by aggregating (function aggregate) using longitude, latitude and month. Discarding
## year as a factor aggregates cells (1° x 1° resolution) of different year but same month to one cell, therefore gridding the
## data.

## This script has been developed by Damiano Righetti.

## Dominic Eriksson
## Environmental Physics Group
## ETH Zurich
## Switzerland

## dominic.eriksson@usys.ethz.ch   


# Preparatory steps
rm(list = ls())
library(raster); library(sfsmisc); library(grid); library(spam); library(fields); library(base); library(doBy); library(ncdf4); library(grid)

# Help functions for aggregation to 1° monthly cell
vec.extr_mean <- function(x){if(all(is.na(x)) == F) {y <- mean(x, na.rm = T)} else {y <- NA}; return(y)} # extract vector for aggregated cells (mean)
vec.extr_sd <- function(x){if(all(is.na(x)) == F) {y <- sd(x, na.rm = T)} else {y <- NA}; return (y)} # extract vector for aggregated cells (standard deviation)

## Create datasets (gridded, ungridded) for two different data-types (presence vs. quantitative-volumetric)
# Get data
tt <- read.csv("/net/kryo/work/deriksson/notion/1_Merge_data/Output_2/Diazo_hierarchic_taxa.csv")
dim(tt) # 30528    55
names(tt)
# Exploration of data types
as.character(unique(tt$occurrenceStatus))

### ======================================================================
### 1. Create dataset (ungridded, gridded) for presence - absence
### ======================================================================

# Presences
# raw
tt1pres <- tt[which(tt$occurrenceStatus%in%c("PRESENT", "Rare (p < 1%)")), ] # subset each row of dataframe for present/rare taxa
dim(tt1pres) # 29602 55
tt1pres$occurrenceStatus <- "PRESENT" # *** Critical decision ***
tt1pres$x <- tt1pres$decimalLongitude
tt1pres$y <- tt1pres$decimalLatitude
tt1pres$taxon <- tt1pres$scientificName
is.numeric(tt1pres$month)
unique(tt1pres$month) # {0 < m < 14} CAREFULL there are more than 12 month
tt1pres$x <- round(tt1pres$x +0.5) - 0.5
tt1pres$y <- round(tt1pres$y +0.5) - 0.5
# gridded, this means to remove the duplicates
sum(duplicated(tt1pres[, c("taxon", "scientificName", "x", "y", "year", "month")])) 22591
# Now damiano removes all duplicate rows and only keeps the column c("phylum", "class", "taxon", "scientificName", "x" ...
tt1presgridded <- tt1pres[!duplicated(tt1pres[, c("taxon", "scientificName", "x", "y", "year", "month")]), c("phylum", "class", "taxon", "scientificName", "x", "y", "year", "month", "nifQuantity", "nifQuantityType", "occurrenceStatus")]
dim(tt1presgridded) # 7011 11

# Non-detections
# raw
tt1abs <-tt[which(tt$occurrenceStatus%in%c("Desconocido", "NOT_DETECTED", "ABSENT")), ] # 926 55
tt1abs$occurrenceStatus <- "NOT_DETECTED" # *** Critical decision ***
tt1abs$x <- tt1abs$decimalLongitude
tt1abs$y <- tt1abs$decimalLatitude
tt1abs$taxon <- tt1abs$scientificName
is.numeric(tt1abs$month)
unique(tt1abs$month) # {0 < m < 13} month are in correct range
tt1abs$x <- round(tt1abs$x +0.5) -0.5
tt1abs$y <- round(tt1abs$y +0.5) -0.5
# gridded
sum(duplicated(tt1abs[ , c("taxon", "scientificName", "x", "y", "year", "month")])) # 761
tt1absgridded <- tt1abs[!duplicated(tt1abs[ , c("taxon", "scientificName", "x", "y", "year", "month")]) ,  c("phylum", "class", "taxon", "scientificName", "x", "y", "year", "month", "nifQuantity", "nifQuantityType", "occurrenceStatus")]
dim(tt1absgridded) # 165 11

# Save raw and gridded
tt1raw <- rbind(tt1pres, tt1abs)
fln <- "/net/kryo/work/deriksson/notion/3_Match_up/Output_1/diazos_pres_abs_raw.csv"
write.csv(tt1raw, file = fln, row.names = F)
tt1gri <- rbind(tt1presgridded, tt1absgridded)
fln <- "/net/kryo/work/deriksson/notion/3_Match_up/Output_1/diazos_pres_abs_gridded.csv"
write.csv(tt1gri, fln, row.names = F)

# Explore
dim(tt1raw) # 30528
dim(tt1gri) # 7176
unique(as.character(tt1raw$scientificName)) # 57
unique(as.character(tt1gri$scientificName)) # 57

### ======================================================================
### 1. Create dataset (ungridded, gridded) for nifH counts - absence
### ======================================================================

# Volumetric nifH-abundances (based on qPCR)
# raw
tt2pres <- tt[which(tt$nifQuantity > 0 & tt$nifQuantityType == "nifH_gene_copies_x10^6_per_m3"), ]
dim(tt2pres) # 9926
tt2pres$x <- tt2pres$decimalLongitude
tt2pres$y <- tt2pres$decimalLatitude
tt2pres$taxon <- tt2pres$scientificName
is.numeric(tt2pres$month)
unique(tt2pres$month) # {0 < m < 13} correct range of months
tt2pres$x<-round(tt2pres$x +0.5)-0.5
tt2pres$y<-round(tt2pres$y +0.5)-0.5
# gridded: here, during the gridding, information needs to be aggregated; specifically, we want the mean of the abundance counts per monthly-cell (per taxon), hence preparatory: split data into taxa
dat.split <- split(tt2pres, f = as.character(tt2pres$taxon)) # 18 Taxa
lisi.dat <- list()
for(i in 1:length(dat.split)){
  # get taxon-wise data
  raw <- dat.split[[i]]
  # 1-4. Get coordinates and taxon: we discard year as a factor
  tax <- aggregate(scientificName ~ x + y + month, data = raw, na.action = na.pass, function(x){paste(x[1])})
  phylum <- aggregate(phylum ~ x + y + month, data = raw, na.action = na.pass, function(x){paste(x[1])})[4]
  class <- aggregate(class ~ x + y + month, data = raw, na.action = na.pass, function(x){paste(x[1])})[4]
  # 5-6. Get nifH gene counts in an aggregate manner (mean, sd)
  vec.nifHcounts.mean <- aggregate(nifQuantity ~ x + y + month, data = raw, na.action = na.pass, FUN = vec.extr_mean)[4]
  vec.nifHcounts.sd <- aggregate(nifQuantity ~ x + y + month, data = raw, na.action = na.pass, FUN = vec.extr_sd)[4]
  # 7 Number of observation (observation, human ovservation) records per grid cell
  vec.obs.no <- aggregate(nifQuantity ~ x + y + month, data = raw, na.action = na.pass, FUN = function(x){y <- length(x)})[4]
  # occurrenceStatus
  vec.occStatus <- rep("PRESENT", nrow(tax)) # no pull argument needed
  # Create gridded data frame: containing grid-cell aggregated gactors of interest, use pull argument because vec.nifHcounts.mean is still a dataframe and it will not rename the column to "av.nifHcounts", pull argument from dplyr package converts dataframe columns to vector sequences with no column names attached anymore
  gridded1 <- data.frame(tax[, 1:4], phylum, class,
    "av.nifHcount" = pull(vec.nifHcounts.mean, nifQuantity),
    "sd.nifHcounts" = pull(vec.nifHcounts.sd, nifQuantity),
    "observations" = pull(vec.obs.no, nifQuantity),
    "occurrenceStatus" = vec.occStatus)
  # store
  lisi.dat[[i]] <- gridded1
  print(paste(i)) # progress
}
tt2presgridded <- do.call("rbind", lisi.dat)
dim(tt2presgridded) # 1864 10
names(tt2presgridded)
head(tt2presgridded)
tail(tt2presgridded)

# Non-detections
# raw
tt2abs <- tt1abs
# gridded: during the gridding, information needs to be aggregated; specifically, we want the mean of the abundance counts per monthly per monthly-cell (per taxon), hence preparatory: split data into taxa
dat.split <- split(tt2abs, f = as.character(tt2abs$taxon))
length(dat.split) # 8 taxa show absences
lisi.dat <- list()
for(i in 1: length(dat.split)){
  # get taxon-wise data
  raw <- dat.split[[i]]
  # 1-4. Get coordinates and taxon: we discard year as a factor
  tax <- aggregate(scientificName ~ x + y + month, data = raw, na.action = na.pass, function(x){paste(x[1])})
  phylum <- aggregate(phylum ~ x + y + month, data = raw, na.action = na.pass, function(x){paste(x[1])})[4]
  class <- aggregate(class ~ x + y + month, data = raw, na.action = na.pass, function(x){paste(x[1])})[4]
  # 5-6. Get nifH gene counts in an aggregate manner (mean, sd)
  vec.nifHcounts.mean <- aggregate(nifQuantity ~ x + y + month, data = raw, na.action = na.pass, FUN = vec.extr_mean)[4]
  vec.nifHcounts.sd <- aggregate(nifQuantity ~ x + y + month, data = raw, na.action = na.pass, FUN = vec.extr_sd)[4]
  # 7 Number of observation (observation, human ovservation) records per grid cell
  vec.obs.no <- aggregate(nifQuantity ~ x + y + month, data = raw, na.action = na.pass, FUN = function(x){y <- length(x)})[4]
  # occurrenceStatus
  vec.occStatus <- rep("PRESENT", nrow(tax)) # no pull argument needed
  # Create gridded data frame: containing grid-cell aggregated gactors of interest, use pull argument because vec.nifHcounts.mean is still a dataframe and it will not rename the column to "av.nifHcounts", pull argument from dplyr package converts dataframe columns to vector sequences with no column names attached anymore
  gridded1 <- data.frame(tax[, 1:4], phylum, class,
    "av.nifHcount" = pull(vec.nifHcounts.mean, nifQuantity),
    "sd.nifHcounts" = pull(vec.nifHcounts.sd, nifQuantity),
    "observations" = pull(vec.obs.no, nifQuantity),
    "occurrenceStatus" = vec.occStatus)
  # store
  lisi.dat[[i]] <- gridded1
  print(paste(i)) # progress
}
tt2absgridded <- do.call("rbind", lisi.dat)
dim(tt2absgridded) # 1740 10
names(tt2absgridded)
head(tt1absgridded)

# Save raw and gridded
tt2raw <- rbind(tt2pres, tt2abs)
fln <- "/net/kryo/work/deriksson/notion/3_Match_up/Output_1/diazos_nifHcounts_abs_raw.csv"
write.csv(tt2raw, file = fln, row.names = F)
tt2gri <- rbind(tt2presgridded, tt2absgridded)
fln <- "/net/kryo/work/deriksson/notion/3_Match_up/Output_1/diazos_nifHcounts_abs_gridded.csv"
write.csv(tt2gri, fln, row.names = F)
