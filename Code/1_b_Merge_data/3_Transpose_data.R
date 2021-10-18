### PREPARE RAW DATA FRAME (PRESENCE-ONYLY vs. ASSUMED ABSENCE MATRIX, OR ABUNDANCE MATRIX) WITH SAMPLE-ID AS ROWS AND SPECIES AS COLUMNS
### MARINE DIAZOTRIOPHS

# Date: 2020-06-09
# Author: Damiano Righetti
# Environmental Physics Group, Institute of Biogeochemistry
# and Pollutant Dynamics, ETH Zurich, CH

### =========================================================================
### Initialize system
### =========================================================================
rm(list = ls())
library(sfsmisc); library(reshape2)
setwd("C:/Users/Dominic Eriksson/Desktop/ETH Internship/Data_notion") # User one
setwd("~/Desktop/Data_notion") # User two

### =========================================================================
### Load and prepare masked (raw) data (PhytoBase) with x, y, year, month, day and depth as factors
### =========================================================================

### Preparatory: Load data
dat <- read.csv("~/Desktop/Data_notion/1_Merge_data/Output_2/Diazo_hierarchic_taxa.csv")
dat <- dat[which(dat$occurrenceStatus=="PRESENT"), ]
dim(dat) # 29 594  x  55 

# Create a pseudo-sample ID: same time and place of sampling
dat$ID.sample <- paste0(as.character(dat$decimalLongitude), "_", as.character(dat$decimalLatitude), "_", as.character(dat$year), "_", as.character(dat$month),"_",as.character(dat$day),"_", as.character(dat$depth))

# Number of taxa
unique(dat$scientificName) # 57

# Check number of unique samples
length(unique(dat$ID.sample)) # 9 773

### =========================================================================
### DIAZOS PA DATASET: FULL DEPTHS OR NA DEPTHS (I.E., RECORDS NOT LIMITED TO MLD)
### =========================================================================

	# Rename data
	tot <- dat
	
	## Transpose data (takes some 2 minutes)
	frame <-as.data.frame.matrix(table(tot$ID.sample, tot$scientificName))
	dim(frame) # 9 773 samples x 57 species or taxa
	
	## Test if class of columns as numeric changes the values
	class(frame$"Chaetoceros") # integer
	frame$"Chaetoceros"
	as.numeric(as.character(frame$"Chaetoceros")) # ok, looks the same
	
	# Preparatory 1: IDs of samples, merged later as first column with data frame
	df.1 <- as.data.frame(rownames(frame)) # looks reasonable
	names(df.1) <- c("ID.sample")
	df.1$ID.sample <- as.character(df.1$ID.sample)
	df.1$running_nr <- 1:nrow(df.1)
		
	# Preparatory 2: add information on source dataset
	set_of_data_collection <- tot[ , c("ID.sample", "set")] # 
	df.1$set <- set_of_data_collection$set[match(df.1$ID.sample, set_of_data_collection$ID.sample)] # Match information from the raw data set with IDs of the melted dataset
	
	# Implement: Merge tot.df with frame
	frame.final <- cbind(df.1, frame)
	rownames(frame.final) <- 1:nrow(frame.final)
	frame.final$ID.sample <- as.character(frame.final$ID.sample)
	frame.final[1:10, ]
	
	# Explore the new data frame
	dim(frame.final) # 13'008 samples x 62 species
	sum(duplicated(frame.final$running_nr)) # 0
	frame.final$running_nr <- NULL
	frame.final[1:10, ]
	
	# Get all values as numeric: change duplicated entries per pseudo-sample (e.g. same record obtained from various sources) to one entry
	frame.final[ , c(3:ncol(frame.final)) ] <- apply(frame.final[ , c(3:ncol(frame.final)) ] , 2, function(x) sign(as.numeric(as.character(x)))) # takes 1 min
	
	# Split ID into x, y, year, month, day, depth again
	class(frame.final$ID) # Ok
	datsplits<- strsplit(frame.final$ID.sample,split="_")
	co.dat <-as.data.frame(do.call(rbind,datsplits))
	names(co.dat)<-c("decimalLongitude", "decimalLatitude", "year","month","day", "depth")
	co.dat$decimalLongitude <- as.numeric(as.character(co.dat$decimalLongitude))
	co.dat$decimalLatitude <- as.numeric(as.character(co.dat$decimalLatitude))
	co.dat$year <- as.numeric(as.character(co.dat$year))
	co.dat$month <- as.numeric(as.character(co.dat$month))
	co.dat$day <- as.numeric(as.character(co.dat$day))
	co.dat$depth <- as.numeric(as.character(co.dat$depth))
	summary(co.dat) # *** COMMENT: Some strange longitude values still there.. and some strange year, month, day values persist.. and deep depths.. these are all subject to further cleaning
	
	# Merge with frame.final
	tt.new <- cbind(co.dat, frame.final)
	tt.new[1:10, ]
	dim(tt.new) # 9 774   x  65   (species or taxa starting from column nine)
	
	# Save
	fln <-"./1_Merge_data/Output_3/Diazo_transposed.csv"
	write.csv(tt.new, file = fln, row.names = F) # be careful not to overwrite
	
	
### =========================================================================
### DIAZOS ABUNDANCE DATASET: FULL DEPTHS OR NA DEPTHS (I.E., RECORDS NOT LIMITED TO MLD) - KEEPING ABUNDANCES INSIDE THE MATRIX - TO BE IMPLEMENTED
### =========================================================================
	
	# From here: yet to be implemented.....
	# select data with abundance information

	# tot <- tt[which(tt$organismQuantity >=1 & tt$organismQuantityType=="number_of_cells_per_L"), ]
	# dim(tot) # 33'490 

	## Test set: Long to wide format data using the "dcast funtion"
	# dat <- tot[200:220,]
	# md1 <- dcast(dat[ , c("ID.sample", "organismQuantity", "taxon")], ID.sample ~ taxon, value.var = "organismQuantity", mean, na.rm=T)	
	# md1 <- dcast(dat[ , c("ID.sample", "cellsPerLitre", "taxon")], ID.sample ~ taxon, value.var = "cellsPerLitre", mean, na.rm=T)
	
	## Full set:
	# frame <- dcast(tot[ , c("ID.sample", "taxon", "organismQuantity")], ID.sample ~ taxon, value.var = "organismQuantity", mean, na.rm=T)
	# dim(frame) # 8'517 samples
	
	## Add  "x", "y", "year", "month", "day", "depth" back to data
	# date.data.frame <- data.frame(do.call("rbind", strsplit(as.character (frame$ID.sample),"_" ,fixed =T)  ))
	# names(date.data.frame) <- c("x","y","year","month","day","depth")
	# date.data.frame$x <- as.numeric(as.character(date.data.frame$x))
	# date.data.frame$y <- as.numeric(as.character(date.data.frame$y))
	# date.data.frame$month <- as.numeric(as.character(date.data.frame$month))
	# date.data.frame$day <- as.numeric(as.character(date.data.frame$day))
	# date.data.frame$depth <- as.numeric(as.character(date.data.frame$depth))
	# frame <- cbind(date.data.frame, frame)

	## Test if class of columns as numeric changes the values
	# class(frame$"Emiliania huxleyi") # numeric. OK
	# class(frame$x) # numeric. OK
	
	# Explore the new data frame
	# frame[1:10, 1:30]
	# dim(frame) # 8517 samples x 329 species and 7 positional variables
	
	# Save: (note: data are to large to save as .csv)
	# fln <- paste0("/..._ABUNDANCES.csv")
	# write.csv(frame, file = fln, row.names = F) # be careful not to overwrite

