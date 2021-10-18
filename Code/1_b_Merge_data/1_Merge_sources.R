### MERGER: TANG AND CASSAR (2019), LUO ET AL (2012, ESSD), GRADOVILLE ET AL (2020)
### PHYTOBASE (2020), GBIF (MAY 2020), OBIS (MAY 2020), TARA (MgnifY OTU BASED), MARTINEZ-PEREZ ET AL (2016)

## Note: no cleaning with regards to missing coordinate, negative depth, unrealistic date or non-oceanic taxa has been performed up to this stage
# Date: 2020-06-02
# Author: Dominic Eriksson
# Environmental Physics Group, ETH Zurich, CH

# This script is maintained via github sdgds

### ===================================================================================
### Preparatory steps
### ===================================================================================
rm(list = ls())
library(doBy);library(plyr);library(doParallel);library(foreach);library(tidyverse)

setwd("C:/Users/Dominic Eriksson/Desktop/ETH Internship/Data_notion/1_Merge_data") # user one
setwd("~/Desktop/Data_notion/1_Merge_data") # user two

# get list of focal taxa
tax.table<-read.csv("../1_Download_data/commons/Taxon_table.csv")

# get desired columns
source("../1_Download_data/commons/desired_columns.R", chdir = TRUE)

# *** Comment: overview of columns desired, can be deleted again
# [1] "decimalLongitude"                       "decimalLatitude"                        "year"
# [4] "month"                                  "day"                                    "depth"
# [7] "depthAccuracy"                          "depthIntegral"                          "scientificName"
# [10] "taxonRank"                              "phylum"                                 "class",     "family",    # *** Note: family added, also to the commons sheet!
# [13] "genus"                                  "species"                                "lifeForm"
# [16] "associatedTaxa"                         "occurrenceStatus"                       "trichomeStatus"
# [19] "trichomeQuantity"                       "trichomeQuantityType"                   "cellStatus"
# [22] "cellQuantity"                           "cellQuantityType"                       "nifStatus"
# [25] "nifQuantity"                            "nifQuantityType"                        "hitStatus"
# [28] "hitQuantity"                            "hitQuantityType"                        "biomass_conversion_factor"
# [31] "biomass_conversion_factor_QuantityType" "biomassQuantity"                        "biomassQuantityType"
# [34] "basisOfRecord"                          "measurementMethod"                      "references"
# [37] "occurrenceID"                           "sampleID"                               "stationID"
# [40] "cruiseID"                               "Sea_surface_temperature_C"              "Temperature_C"
# [43] "Sea_surface_Salinity"                   "Salinity"                               "Sea_surface_nitrate_micromolar"
# [46] "Nitrate_micromolar"                     "Sea_surface_phosphate_micromolar"       "Phosphate_micromolar"
# [49] "Sea_surface_Fe_nanomolar"               "Fe_nanomolar"                           "ChlorophyllQuantity"
# [52] "ChlorophyllQuantityType"                "flag"

# In addition, we add a column called "set". Indicating the data source / dataset downloaded

### ===================================================================================
### Get data sets adjust columns
### ===================================================================================

## 1. Luo et al.
ll<-read.csv("../1_Download_data/Luo_et_al_2012/Output_2/Merged_luo.csv",header = TRUE)
desired.cols[desired.cols%in%colnames(ll)==FALSE] # print the columns missing
dim(ll)# 9696 41

# *** Some comments and priorites for further cleaning:
# A key important remaining task is probably to complete and verify the occurrenceStatus (as presences are useful for modeling)
# A second priority is to clean/check the day entries and/or month entries (these are mistaken with year, at current, in many instances). I do a first / quick cleaning here!
# A third remaining task may be to check what happened with the coordinates stemming from OBIS. There are some 2526 NA entries....
# A fourth task is to check that volumentric entries are comparable (i.e. scaled to the target volume, e..g. counts per m^3 , as defined in the commons sheet)
# and perhaps the basis of record (second priority), and to retain biomass from tang and cassar (now already implemented),
# and to retain the environmental data from tang and cassar (now I already implemented this)
# Now implemented: skript follows the order of the desired columns of the commons sheet, per dataset, in a systematic manner. This avoids mistakes.
# Take special care with unit conversions: copies != copies*10^6
# Next important step (partially addressed and solved below): retain "month" in obis data. date_start and date_end are corrupted in the raw file already. eventDate looks better.
# Another interesting next step is to retain possible / or key / host taxa. This is implemented now below for the TARA oceans dataset, not the others. We may have to add a column: isDiazotroph (YES, NO...) via commons sheet.
# Own insight: replacing and changing vlaues for year, month, day should be done as characer. only in the end change to numeric. else the operations may not work. (I did this for Tara otu based data, but not yet clean-checked for other sets)
# Naming conventions: assign Cyanothece to UCYN.C in gbif or obis dataset. (following the scientificName entries in our taxon sheet)
# Finally, to keep track for ourselves I suggest to add a colum termed "setOrigin" - denotes if a record stems from the dataset of Luo, Tang and Cassar, Gradoville etc. or their combinations

# a. Get longitude, latitude, year, month, day, depth
ll$decimalLongitude<-ll$LONGITUDE
ll$decimalLatitude<-ll$LATITUDE
# --- inset: get year, month, day ---
# *** comment: year. the year needs to be extracted from the ll$DATE..yyyy.mm.dd. column -> there are different formats in the raw data, I have carefully standardized those for tang and cassar going through the excel file
# but... here in luo, it seems that there are two different formats (ideally to be cleaned beforehands and satandardized) i.e., YEAR-month-day versus month/day/YEAR...  I encode all things here (but careful check missing!) towards: YEAR-month-day
ll$DATE..yyyy.mm.dd. <- as.character(ll$DATE..yyyy.mm.dd.) # adjust class of date-column to make it workable
dates.with.slash <- gsub("/", "-", x = ll$DATE..yyyy.mm.dd. [  grep("/", x = ll$DATE..yyyy.mm.dd.)]  ) # adjusting all dates with separator "-"
df.thereof <- data.frame(do.call("rbind", strsplit(dates.with.slash, "-", fixed = T)))
reassembled.dates <- paste0(df.thereof[,3], "-",df.thereof[,2], "-",df.thereof[,1]) # assemble to desired format
ll[ grep("/", x = ll$DATE..yyyy.mm.dd.) , "DATE..yyyy.mm.dd."] <- reassembled.dates # replace the dates-to-be adjusted in ll
df.date <- data.frame(do.call("rbind", strsplit(ll$DATE..yyyy.mm.dd., "-", fixed = T))) # day, month, year
colnames(df.date) <- c("year", "month", "day")
# --- end inset ---
ll$year <- df.date$year # *** comment: added
ll$month <- df.date$month # *** comment: added
ll$day <- df.date$day # *** comment: added
class(ll$year) # *** comment: inserted from below
class(ll$month) # *** comment: inserted from below
class(ll$day) # *** comment: inserted from below
ll$year <- as.numeric(as.character(ll$year)) # *** comment: inserted from below
ll$month <- as.numeric(as.character(ll$month)) # *** comment: inserted from below
ll$day <- as.numeric(as.character(ll$day)) # *** comment: inserted from below
ll$depth<-as.numeric(as.character(ll$DEPTH..m.)) # *** comment: I define as numeric
ll$depthAccuracy<-NA
ll$depthIntegral<-as.numeric(as.character(ll$INTEGRAL.DEPTH..m.)) # *** comment: I define as numeric

# b. Get lifeForm from specific dataset (perhaps, rathern than from the generic Taxon_table file) and associated taxa
ll$lifeForm<-"NA" # information will be added after big merger . *** Comment: I think we preferrably add this info at the dataset level - if data area available, as sometimes the same species live as symbionts, sometimes not? (However: This is a second priority information!)
ll$associatedTaxa <- "NA" # Is there some information on the associated taxa in Luo et al? I think there might be some information on the diatom hosts (second priority task!)

# c. Get occurrence status
ll$occurrenceStatus <- "NA" # *** Comment: will be filled in below. Column shall ideally indicate "PRESENT", "ABSENT", "NOT_DETECTED" - to be defined in commons sheet, and/or according to information available

# d. Get trichome-based abundances into target format, target column "trichomeQuantity" and "trichomeQuantityType"
ll$Total.Trichomes..x10.6.m.3.<-as.numeric(as.character(ll$Total.Trichomes..x10.6.m.3.))
ll$Total.Trichomes..x10.6.m.2.<-as.numeric(as.character(ll$Total.Trichomes..x10.6.m.2.))
ll$trichomeStatus<-"NA"
ll$trichomeQuantity<-NA
ll$trichomeQuantityType<-"NA"
ll[which(ll$Total.Trichomes..x10.6.m.3.> 0), "trichomeQuantity"]<-ll[which(ll$Total.Trichomes..x10.6.m.3.>0),"Total.Trichomes..x10.6.m.3."]
ll[which(ll$Total.Trichomes..x10.6.m.3. > 0), "trichomeQuantityType"] <- "Total_trichomes_x10^6_per_m3" # *** Comment: In the end, we probably want to use a clean version e.g. "Trichomes_x_10^6_per_m^3". I am implementing this here.
ll[which(ll$Total.Trichomes..x10.6.m.2.> 0), "trichomeQuantity"]<-ll[which(ll$Total.Trichomes..x10.6.m.2.>0),"Total.Trichomes..x10.6.m.2."]
ll[which(ll$Total.Trichomes..x10.6.m.2. > 0), "trichomeQuantityType"] <- "Total_trichomes_x10^6_per_m2" # *** Comment: cleaned
ll[which(ll$Total.Trichomes..x10.6.m.3. > 0), "trichomeStatus"]<-"PRESENT"
ll[which(ll$Total.Trichomes..x10.6.m.3. > 0), "occurrenceStatus"]<-"PRESENT"
ll[which(ll$Total.Trichomes..x10.6.m.2. > 0), "trichomeStatus"]<-"PRESENT"
ll[which(ll$Total.Trichomes..x10.6.m.2. > 0), "occurrenceStatus"]<-"PRESENT"

# e. Get cell-based abundances into target format, target column "cellQuantity" and "cellQuantityType"
ll$Cells..x10.6.m.3.<-as.numeric(as.character(ll$Cells..x10.6.m.3.))
ll$Cells..x10.6.m.2.<-as.numeric(as.character(ll$Cells..x10.6.m.2.))
ll$cellStatus<-"NA"
ll$cellQuantity<-NA
ll$cellQuantityType<-"NA"
ll[which(ll$Cells..x10.6.m.3.> 0), "cellQuantity"]<-ll[which(ll$Cells..x10.6.m.3.>0),"Cells..x10.6.m.3."]
ll[which(ll$Cells..x10.6.m.3. > 0), "cellQuantityType"] <- "Cells_x10^6_per_m3"
ll[which(ll$Cells..x10.6.m.2. > 0), "cellQuantity"] <- ll[which(ll$Cells..x10.6.m.2. > 0), "Cells..x10.6.m.2."]
ll[which(ll$Cells..x10.6.m.2. > 0), "cellQuantityType"] <- "Cells_x10^6_per_m2"
ll[which(ll$Cells..x10.6.m.3. > 0), "cellStatus"]<- "PRESENT"
ll[which(ll$Cells..x10.6.m.3. > 0), "occurrenceStatus"]<- "PRESENT"
ll[which(ll$Cells..x10.6.m.2. > 0), "cellStatus"]<- "PRESENT"
ll[which(ll$Cells..x10.6.m.2. > 0), "occurrenceStatus"] <- "PRESENT"

# f. Get gene-based abundances into target format, target column "nifQuantity" and "nifQuantityType"
ll$nifH.Gene..x10.6.copies.m.3.<-as.numeric(as.character(ll$nifH.Gene..x10.6.copies.m.3.))
ll$nifH.Gene..x10.6.copies.m.2.<-as.numeric(as.character(ll$nifH.Gene..x10.6.copies.m.2.))
ll$nifStatus<-"NA"
ll$nifQuantity<-NA
ll$nifQuantityType<-"NA"
ll[which(ll$nifH.Gene..x10.6.copies.m.3.> 0), "nifQuantity"]<-ll[which(ll$nifH.Gene..x10.6.copies.m.3.>0),"nifH.Gene..x10.6.copies.m.3."]
ll[which(ll$nifH.Gene..x10.6.copies.m.3. > 0), "nifQuantityType"] <- "nifH_gene_copies_x10^6_per_m3" # *** Comment: cleaned
ll[which(ll$nifH.Gene..x10.6.copies.m.2. > 0), "nifQuantity"] <- ll[which(ll$nifH.Gene..x10.6.copies.m.2. > 0), "nifH.Gene..x10.6.copies.m.2."]
ll[which(ll$nifH.Gene..x10.6.copies.m.2. > 0), "nifQuantityType"] <- "nifH_gene_copies_x10^6_per_m2" # *** Comment: cleaned
ll[which(ll$nifH.Gene..x10.6.copies.m.3. > 0), "nifStatus"]<-"PRESENT"
ll[which(ll$nifH.Gene..x10.6.copies.m.3. > 0), "occurrenceStatus"] <- "PRESENT"
ll[which(ll$nifH.Gene..x10.6.copies.m.2. > 0), "nifStatus"]<-"PRESENT"
ll[which(ll$nifH.Gene..x10.6.copies.m.2. > 0), "occurrenceStatus"] <- "PRESENT"

# g. Get hit-based abundances into target format, target column "nifQuantity" and "nifQuantityType"
ll$hitStatus<-"NA"
ll$hitQuantity<-NA
ll$hitQuantityType<-"NA"

# h. Get biomass conversion factor, and biomass conversion factor quantity type
ll$Biomass.Conversion.factor..mg.C.10.6.cells.<-as.numeric(as.character(ll$Biomass.Conversion.factor..mg.C.10.6.cells.))
ll$Biomass.Conversion.factor..mg.C.10.6.nifH.copies.<-as.numeric(as.character(ll$Biomass.Conversion.factor..mg.C.10.6.nifH.copies.))
ll$Biomass.Conversion.factor..mg.C.10.6.trichomes.<-as.numeric(as.character(ll$Biomass.Conversion.factor..mg.C.10.6.trichomes.))
ll$biomass_conversion_factor<-NA
ll$biomass_conversion_factor_QuantityType<-"NA"
ll[which(ll$Biomass.Conversion.factor..mg.C.10.6.cells. > 0), "biomass_conversion_factor"]<-ll[which(ll$Biomass.Conversion.factor..mg.C.10.6.cells. > 0), "Biomass.Conversion.factor..mg.C.10.6.cells."]
ll[which(ll$Biomass.Conversion.factor..mg.C.10.6.cells. > 0),"biomass_conversion_factor_QuantityType"]<- "mgC_per_10^6_cells" # *** Comment: cleaned
ll[which(ll$Biomass.Conversion.factor..mg.C.10.6.nifH.copies. > 0), "biomass_conversion_factor"]<-ll[which(ll$Biomass.Conversion.factor..mg.C.10.6.nifH.copies. > 0), "Biomass.Conversion.factor..mg.C.10.6.nifH.copies."]
ll[which(ll$Biomass.Conversion.factor..mg.C.10.6.nifH.copies. > 0),"biomass_conversion_factor_QuantityType"]<-"mmg.C./10^6.nifH.copies." ## COMMENT: Why mmg.. ? Please check/adjust.
ll[which(ll$Biomass.Conversion.factor..mg.C.10.6.trichomes. > 0),"biomass_conversion_factor"]<-ll[which(ll$Biomass.Conversion.factor..mg.C.10.6.trichomes. > 0),"Biomass.Conversion.factor..mg.C.10.6.trichomes."]
ll[which(ll$Biomass.Conversion.factor..mg.C.10.6.trichomes. > 0), "biomass_conversion_factor_QuantityType"]<-"mgC_per_10^6_trichomes" # *** Comment: cleaned

# i. Get biomass quantity and quantity type to differentiate integrated measurements biomass_mgC_per_m2
ll$Biomass..mg.C.m.3.<-as.numeric(as.character(ll$Biomass..mg.C.m.3.))
ll$Biomass..mg.C.m.2.<-as.numeric(as.character(ll$Biomass..mg.C.m.2.))
ll$biomassQuantity<-NA
ll$biomassQuantityType<-"NA"
ll[which(ll$Biomass..mg.C.m.3. > 0), "biomassQuantity"]<-ll[which(ll$Biomass..mg.C.m.3.>0),"Biomass..mg.C.m.3."]
ll[which(ll$Biomass..mg.C.m.3. > 0), "biomassQuantityType"] <- "mgC_per_m3" # *** Comment: cleaned
ll[which(ll$Biomass..mg.C.m.2. > 0), "biomassQuantity"] <- ll[which(ll$Biomass..mg.C.m.2. > 0), "Biomass..mg.C.m.2."]
ll[which(ll$Biomass..mg.C.m.2. > 0), "biomassQuantityType"] <- "mgC_per_m2" # *** Comment: cleaned
ll[which(ll$biomassQuantity > 0), "occurrenceStatus"] <- "PRESENT" # Note: there are entries in luo et al that only contain information about tbe biomass --> information on cell counts, colonies etc is missing. I will interpret a biomass entry with no further information as a PRESENCE of the organism. *** Comment: excellent. Thus, we get more presences!

# j. Put basis of record, measurement method, and references in place
ll$basisOfRecord<-NA # *** Comment: see commons sheet. Possible to define? Standard light microscopy may correspond to OBSERVATION, qpcr as well
ll$measurementMethod<-ll$METHODS...................Sampling.Analysis
ll$references<-ll$SOURCE..Data # *** Comment: There are a few "NAs" -> could these be completed / filled in / in the raw sheet or not? Please check (second priority). I did this "gap-filling" in tang and cassar!

# k. Get IDs, add flag column
ll$occurrenceID <- paste0(as.character(ll$scientificName),"_",as.character(ll$family),"_", as.character(ll$genus),"_",as.character(ll$decimalLongitude),"_",   # *** Comment: I add genus. Does it make sense? NOTE: "family" is excluded below, as it is not in the target columns!!!! I hence add family to the target columns -> correct?
as.character(ll$decimalLatitude),"_",as.character(ll$year),"_",as.character(ll$month),"_",as.character(ll$day),"_",as.character(ll$depth))
ll$sampleID<-paste0(as.character(ll$decimalLongitude),"_",as.character(ll$decimalLatitude),"_",as.character(ll$year),"_",as.character(ll$month),"_",as.character(ll$day),"_",as.character(ll$depth))
ll$stationID<-NA
ll$cruiseID<-NA
ll$flag<-NA

# l. Get environmental data (if available)
ll$Chlorophyll...mg.m.3.<-as.numeric(as.character(ll$Chlorophyll...mg.m.3.))
ll$Chlorophyll...mg.m.2.<-as.numeric(as.character(ll$Chlorophyll...mg.m.2.))
ll$ChlorophyllQuantity<-NA
ll$ChlorophyllQuantityType<-"NA"
ll[which(ll$Chlorophyll...mg.m.3. > 0), "ChlorophyllQuantity"]<-ll[which(ll$Chlorophyll...mg.m.3. > 0), "Chlorophyll...mg.m.3."]
ll[which(ll$Chlorophyll...mg.m.3. > 0), "ChlorophyllQuantityType"]<-"mg_per_m3" # *** Comment: cleaned
ll[which(ll$Chlorophyll...mg.m.2. > 0), "ChlorophyllQuantity"]<-ll[which(ll$Chlorophyll...mg.m.2. > 0), "Chlorophyll...mg.m.2."]
ll[which(ll$Chlorophyll...mg.m.2. > 0), "ChlorophyllQuantityType"]<-"mg_per_m2" # *** Comment: cleaned
ll$Sea_surface_temperature_C<-ll$Surface.Temperature.._C.
ll$Temperature_C<-ll$Temperature.._C.
ll$Sea_surface_Salinity<-ll$Surface.Salinity..PSU.
ll$Salinity<-ll$Salinity..PSU.
ll$Sea_surface_nitrate_micromolar<-ll$Surface.Nitrate..micromolar.
ll$Nitrate_micromolar<-ll$Nitrate..micromolar.
ll$Sea_surface_phosphate_micromolar<-ll$Surface.Phosphate..micromolar.
ll$Phosphate_micromolar<-ll$Phosphate..micromolar.
ll$Sea_surface_Fe_nanomolar<-ll$Surface.Fe..nanomolar.
ll$Fe_nanomolar<-ll$Fe..nanomolar.

## m. Check if all of the desired columns were addressed, and order dataframe at same time
ll <- ll[, desired.cols] # No complaints, fine.

## n. Check if/where occurrence status has been assigned
unique(ll$occurrenceStatus) # *** Comment: There are some NA's, let us assign these, as well. This is a higher priority task, I would say.
nrow(ll[which(ll$occurrenceStatus=="NA"),]) # 3480 NA's
nrow(ll[which(ll$occurrenceStatus=="PRESENT"),]) # 6216 "PRESENT"
head(ll[which(ll$occurrenceStatus=="NA"),]) # *** Comment: I take a look at the NA cases: there is no biomass entry, no other entry. Are these "non detections" or "presences"....? Can we verfiy with the original luo et al dataset?
unique(ll[which(ll$occurrenceStatus=="NA"),"measurementMethod"]) # *** Comment: Different methods were involved in the potential non-detection measurements.

## o. add set source
ll$set <- "Luo"

# PS: check availability of date and if date makes sense
dim(ll[is.na(ll$month), ]) # 0, fine
unique(ll$day) # Ok, looks fine
unique(ll$month) # Ok, looks fine
unique(ll$year) # Ok, looks fine

##----------------------------------------------------------------------------------------------------------------

## 2. Tang and Cassar 2019
tt <- read.csv("../1_Download_data/Tang_and_Cassar_2019/Output_2/Tang_merged_raw.csv")
desired.cols[desired.cols%in%colnames(tt)==FALSE] # print the columns missing
dim(tt)# 10058 41

# a. Get longitude, latitude, year, month, day, depth
tt$decimalLongitude<- tt$"LONGITUDE"
tt$decimalLatitude <- tt$"LATITUDE"
df.date <- data.frame(do.call("rbind", strsplit(as.character (tt$DATE..yyyy.mm.dd.),"/" ,fixed =T)  ))
names(df.date) <- c("month", "day", "year")
tt$year <- as.numeric(as.character(df.date$year))
tt$month <- as.character(df.date$month) # *** Comment: these are characters, not numeric, from source
tt[ which(tt$month == "Jan") , "month" ] <- 1
tt[ which(tt$month == "Feb") , "month" ] <- 2
tt[ which(tt$month == "Mar") , "month" ] <- 3
tt[ which(tt$month == "Apr") , "month" ] <- 4
tt[ which(tt$month == "May") , "month" ] <- 5
tt[ which(tt$month == "Jun") , "month" ] <- 6
tt[ which(tt$month == "Jul") , "month" ] <- 7
tt[ which(tt$month == "Aug") , "month" ] <- 8
tt[ which(tt$month == "Sep") , "month" ] <- 9
tt[ which(tt$month == "Oct") , "month" ] <- 10
tt[ which(tt$month == "Nov") , "month" ] <- 11
tt[ which(tt$month == "Dec") , "month" ] <- 12
tt[ which(tt$month == "Nov-Dec") , "month" ] <- 11 # ***Comment: Let's try to check with original cruise (e.g. from where to where the cruise went over what time, then assign the correct month... - second priority task)
tt[ which(tt$month == "Jun-Jul") , "month" ] <- 6 # *** (same comment)
tt[ which(tt$month == "Oct-Nov") , "month"] <- 10 # *** (same comment)
tt[ which(tt$month == "Dec-Jan") , "month" ] <- 12 # *** (same comment)
tt$month <- as.numeric(as.character(tt$month)) # *** Comment: inserted
tt$day <- as.numeric(as.character(df.date$day))
tt$depth <- as.numeric(as.character(tt$"DEPTH..m."))
tt$depthAccuracy <- NA # *** comment: added
tt$depthIntegral <- as.numeric(as.character(tt$"Integral.DEPTH..m."))

# b. Get lifeForm from specific dataset (perhaps, rathern than from the generic Taxon_table file) and associated taxa (is already available)
tt$lifeForm<-"NA" # *** Comment: information may be derived from column: associatedTaxa. Second priority task
tt$family <- "NA" # *** COMMENT, family is needed  for the occurrenceID. Please add family to the raw sources / retain family in the raw sources

# c. Get occurrence status
tt$occurrenceStatus <- "NA" # *** Comment: Will be filled in below. Column shall ideally indicate "PRESENT", "ABSENT", "NON-DETECTED" or others.. to be defined in commons sheet.

# d. Get trichome-based abundances into target format, target column "trichomeQuantity" and "trichomeQuantityType"
tt$trichomeStatus<-"NA"
tt$trichomeQuantity<-NA
tt$trichomeQuantityType<-"NA"

# e. Get cell-based abundances into target format, target column "cellQuantity" and "cellQuantityType"
tt$cellStatus<-"NA"
tt$cellQuantity<-NA
tt$cellQuantityType<-"NA"

# f. Get gene-based abundances into target format, target column "nifQuantity" and "nifQuantityType"
tt$nifH_Gene_x10.6_copies_m.3 <- as.numeric(as.character(tt$nifH_Gene_x10.6_copies_m.3))
tt$nifH_Gene_x10.6_copies_m.2 <- as.numeric(as.character(tt$nifH_Gene_x10.6_copies_m.2))
tt$nifStatus <- "NA" # Character vector -> note there are many types of data that need to be carefully considered here in a next version of this script: N.d. - not detected - nd - etc...  e.g.,
tt$nifQuantity <- NA
tt$nifQuantityType <- "NA"
tt[which(tt$nifH_Gene_x10.6_copies_m.3 > 0), "nifQuantity"] <- tt[which(tt$nifH_Gene_x10.6_copies_m.3 > 0), "nifH_Gene_x10.6_copies_m.3"]
tt[which(tt$nifH_Gene_x10.6_copies_m.3 > 0), "nifQuantityType"] <- "nifH_gene_copies_x10^6_per_m3" # *** Comment: this was wrong, I change from "nifH_gene_copies_per_m3" to "nifH_gene_copies_x10^6_per_m3" in accordance with Luo.
tt[which(tt$nifH_Gene_x10.6_copies_m.2 > 0), "nifQuantity"] <- tt[which(tt$nifH_Gene_x10.6_copies_m.2 > 0), "nifH_Gene_x10.6_copies_m.2"]
tt[which(tt$nifH_Gene_x10.6_copies_m.2 > 0), "nifQuantityType"] <- "nifH_gene_copies_x10^6_per_m2" # *** Comment: I correct from "nifH_gene_copies_per_m2" to "nifH_gene_copies_x10^6_per_m2" in accordance with Luo.
tt[which(tt$nifH_Gene_x10.6_copies_m.3 == "n.d."), "nifStatus"] <- "NOT_DETECTED"
tt[which(tt$nifH_Gene_x10.6_copies_m.3 == "n.d."), "occurrenceStatus"] <- "NOT_DETECTED" # *** Commment: I added this column. (But is this correct?)
tt[which(tt$nifH_Gene_x10.6_copies_m.3 > 0 ), "nifStatus"] <- "PRESENT"
tt[which(tt$nifH_Gene_x10.6_copies_m.3 > 0 ), "occurrenceStatus"] <- "PRESENT"
tt[which(tt$nifH_Gene_x10.6_copies_m.2 > 0 ), "nifStatus"] <- "PRESENT"
tt[which(tt$nifH_Gene_x10.6_copies_m.2 > 0 ), "occurrenceStatus"] <- "PRESENT"

# g. Get hit-based abundances into target format, target column "nifQuantity" and "nifQuantityType"
tt$hitStatus<-"NA"
tt$hitQuantity<-NA
tt$hitQuantityType<-"NA"

# h. Get biomass conversion factor, and biomass conversion factor quantity type
tt$biomass_conversion_factor<-as.numeric(as.character(tt$Biomass_conversion_factor_mg_C.10.6_nifH_copies)) # *** Comment: added
tt$biomass_conversion_factor_QuantityType<-"NA" #  *** Comment: added - there is no such information in tang and cassar, right? Perhaps this column is only needed in the context of trichomes. Perhaps briefly check with commons file!

# i. Get biomass quantity and quantity type to differentiate integrated measurements biomass_mgC_per_m2
tt$Biomass_mgC_m.3<-as.numeric(as.character(tt$Biomass_mgC_m.3)) # *** Comment: added
tt$Biomass_mg.C_m.2<-as.numeric(as.character(tt$Biomass_mgC_m.2)) # *** Comment: added
tt$biomassQuantity<-NA # *** Comment: added
tt$biomassQuantityType<-"NA" # *** Comment: added
tt[which(tt$Biomass_mgC_m.3 > 0), "biomassQuantity"]<-tt[which(tt$Biomass_mgC_m.3>0),"Biomass_mgC_m.3"] # *** Comment: added
tt[which(tt$Biomass_mgC_m.3 > 0), "biomassQuantityType"] <- "mgC_per_m3" # *** Comment: added, matching the expression of Luo
tt[which(tt$Biomass_mgC_m.2 > 0), "biomassQuantity"] <- tt[which(tt$Biomass_mgC_m.2 > 0), "Biomass_mgC_m.2"] # *** Comment: added
tt[which(tt$Biomass..mg.C.m.2. > 0), "biomassQuantityType"] <- "mgC_per_m2" # *** Comment: added, matching the expression of Luo
tt[which(tt$biomassQuantity > 0), "occurrenceStatus"] <- "PRESENT" # *** Comment: added

# j. Put basis of record, measurement method, and references in place
tt$basisOfRecord<-NA # *** Comment: Standard light microscopy may correspond to "OBSERVATION", qpcr as well
tt$references <- tt$SOURCE..Data
tt$measurementMethod <- "qPCR" # A check of methodological consistency across sources within Tang and Cassar may be very valuable... (second order priority!)

# k. Get IDs, add flag column
tt$occurrenceID <- paste0(as.character(tt$scientificName),"_",as.character(tt$family),"_", as.character(tt$genus),"_",as.character(tt$decimalLongitude),"_",   # *** Comment: I add genus. Does it make sense?
as.character(tt$decimalLatitude),"_",as.character(tt$year),"_",as.character(tt$month),"_",as.character(tt$day),"_",as.character(tt$depth))
tt$sampleID<-paste0(as.character(tt$decimalLongitude),"_",as.character(tt$decimalLatitude),"_",as.character(tt$year),"_",as.character(tt$month),"_",as.character(tt$day),"_",as.character(tt$depth))
tt$stationID<-NA
tt$cruiseID<-NA
tt$flag<-NA

# l. Get environmental data (if available)
tt$Chlorophyll...mg.m.3.<-as.numeric(as.character(tt$Chlorophyll...mg.m.3.))
tt$Chlorophyll...mg.m.2.<-as.numeric(as.character(tt$Chlorophyll...mg.m.2.))
tt$ChlorophyllQuantity<-NA #  *** Note: This and similar lines could be avoided in the script, by direcly using: ll$ChlorophyllQuantity<-as.numeric(as.character(tt$Chlorophyll...mg.m.3.))
tt$ChlorophyllQuantityType<-"NA"
tt[which(tt$Chlorophyll...mg.m.3. > 0), "ChlorophyllQuantity"]<-tt[which(tt$Chlorophyll...mg.m.3. > 0), "Chlorophyll...mg.m.3."]
tt[which(tt$Chlorophyll...mg.m.3. > 0), "ChlorophyllQuantityType"]<-"mg_m3" # *** Comment: cleaned
tt[which(tt$Chlorophyll...mg.m.2. > 0), "ChlorophyllQuantity"]<-tt[which(tt$Chlorophyll...mg.m.2. > 0), "Chlorophyll...mg.m.2."]
tt[which(tt$Chlorophyll...mg.m.2. > 0), "ChlorophyllQuantityType"]<-"mg_per_m2" # *** Comment: cleaned
tt$Sea_surface_temperature_C<-tt$Surface.Temperature.._C.
tt$Temperature_C<-tt$Temperature.._C.
tt$Sea_surface_Salinity<-tt$Surface.Salinity..PSU.
tt$Salinity<-tt$Salinity..PSU.
tt$Sea_surface_nitrate_micromolar<-as.numeric(as.character(tt$Surface.Nitrate..micromolar.))
tt$Nitrate_micromolar<-tt$Nitrate..micromolar.
tt$Sea_surface_phosphate_micromolar<-tt$Surface.Phosphate..micromolar.
tt$Phosphate_micromolar<-tt$Phosphate..micromolar.
tt$Sea_surface_Fe_nanomolar<-tt$Surface.Fe..nanomolar.
tt$Fe_nanomolar<-tt$Fe..nanomolar.

## m. Check if all of the desired columns were addressed, and order dataframe at same time
tt <- tt[, desired.cols] # No complaints, fine.

## n. Check if/where occurrence status has been assigned
unique(tt$occurrenceStatus) # *** Comment: There are some NA's, let us assign these, as well. This is a high priority task, I would say.
nrow(tt[which(tt$occurrenceStatus=="NA"),]) # 4133 NA's
nrow(tt[which(tt$occurrenceStatus=="PRESENT"),]) # 5925 "PRESENT"
head(tt[which(tt$occurrenceStatus=="NA"),]) # *** Comment: I take a look at the NA cases: there is no biomass entry, no other entry. Are these "non detections" or "presences"....? Please verfiy
unique(tt[which(tt$occurrenceStatus=="NA"),"measurementMethod"]) # *** Comment: all qpcr

## o. add set source
tt$set <- "Tan"

# PS: check availability of date
dim(tt[is.na(tt$month), ]) # 0, fine
unique(tt$day) # Ok, looks fine, one NA, though
unique(tt$month) # Ok, looks fine
unique(tt$year) # Ok, looks fine

##-----------------------------------------------------------------------------------------------

## 3. Gradoville et al. 2020
gg <-  read.csv("../1_Download_data/Gradoville_Farnelid_2020/Output_2/Gradoville_merged.csv")
desired.cols[desired.cols%in%colnames(gg)==FALSE] # print the columns missing

# a. Get longitude, latitude, year, month, day, depth
gg$decimalLongitude <- -gg$"Longitude..W." # Longitude West becomes negative longitude
gg$decimalLatitude <- gg$"Latitude..N."
gg$year <- NA
gg[ which(gg$Year == "April-May/19-03/2016") , "year" ] <- 2016
gg[ which(gg$Year == "April-May/19-03/2016") , "month" ] <- 4 # ***Comment: The current choice is not based on facts. Please check with Mary R, how the original cruise went (e.g. from where to where the cruise went over what time, then specify the month)
gg[ which(gg$Year == "April-May/19-03/2016") , "day" ] <- NA
gg[ which(gg$Year == "27-13/May-June/2017") , "year" ] <- 2017
gg[ which(gg$Year == "27-13/May-June/2017") , "month" ] <- 5 # # ***Comment: The current choice is not based on facts. Please check with Mary R, how the original cruise went (e.g. from where to where the cruise went over what time, then specify the month)
gg[ which(gg$Year == "27-13/May-June/2017") , "day" ] <- NA
gg$depth <- NA # *** Comment: Please check in the paper or with Mary R Gracoville, if these were sea surface samples, collected at what depth (?)
gg$depthAccuracy<-NA # *** Comment: dadded
gg$depthIntegral <- NA
# get class correct
gg$year <- as.numeric(as.character(gg$year))
gg$month <- as.numeric(as.character(gg$month))
gg$day <- as.numeric(as.character(gg$day))
gg$depth <- as.numeric(as.character(gg$depth))

# b. Get lifeForm and associated taxa, and genus plus species
gg$lifeForm<-"NA" # *** Comment: added. No information, I guess (However: This is a second priority information!)
gg$associatedTaxa <- "NA" # *** Comment: added. Is there some information? (second priority task!)
gg$species <- "Atelocyanobacterium thalassa" # *** Comment: added. Useful for our purposes
gg$genus <- "Atelocyanobacterium" # *** Comment: added. Useful for our purposes
tt$family <- "Aphanothecaceae" # *** COMMENT, family is needed  for the occurrenceID. Please add family to the raw sources / retain family in the raw sources

# c. Get occurrence status
gg$occurrenceStatus <- "NA" # *** Comment: will be filled in below. Column shall ideally indicate "PRESENT", "ABSENT", "NOT_DETECTED", or similar

# d. Get trichome-based abundances into target format, target column "trichomeQuantity" and "trichomeQuantityType"
gg$trichomeStatus<-NA # *** Comment: added
gg$trichomeQuantity<-NA # *** Comment: added
gg$trichomeQuantityType<-"NA" # *** Comment: added

# e. Get cell-based abundances into target format, target column "cellQuantity" and "cellQuantityType"
gg$cellStatus<-"NA" # *** Comment: added
gg$cellQuantity<-NA #  *** Comment: added
gg$cellQuantityType<-"NA" #  *** Comment: added

# f. Get gene-based abundances into target format, target column "nifQuantity" and "nifQuantityType"
help.vector <- as.character(gg$nifH_Gene_copies_L.1) # *** Comment: Inserted, to retain the information on non-detection (this info lost when changing to numeric)
gg$nifH_Gene_copies_L.1 <- as.numeric(as.character(gg$nifH_Gene_copies_L.1)) # *** ATTENTION: as.numeric(as.character(....)) is needed here, not just as character.. all values had changed
gg$nifStatus<-"NA"
gg$nifQuantity<-NA
gg$nifQuantityType<-"NA"
gg[which(help.vector == "UD"), "nifStatus"] <- "NOT_DETECTED" # *** COMMENT: Now we can use the help.vector here to add the information on non-detection!
gg[which(help.vector == "DNQ"), "nifStatus"] <- "DETECTED_BUT_NOT_QUANTIFIED" # *** COMMENT: evaluate if it is worth to use and add to option in the commons file -> could we change this entry simply to PRESENT ?
gg[which(help.vector == "ND"), "nifStatus"] <- "NOT_DETERMINED" # *** COMMENT: Does this mean taxon present, or absent, or present, but not measured (same as just above)? Please check. If it means nothing, just exclude...
gg[which(gg$nifH_Gene_copies_L.1 > 0), "nifStatus"] <- "PRESENT"
gg[which(gg$nifH_Gene_copies_L.1 > 0), "occurrenceStatus"] <- "PRESENT"
gg$nifH_Gene_x10.6_copies_m.3 <- (as.numeric(as.character(gg$nifH_Gene_copies_L.1)) * 1000)/(10^6) # Multiply by 1000 (to go from Litres to m^3) - *** COMMENT: ATTENTION - on the left you say copies x 10^6. You thus also have to divide by 10^6
gg[which(gg$nifH_Gene_x10.6_copies_m.3 > 0), "nifQuantity"] <- gg[which(gg$nifH_Gene_x10.6_copies_m.3 > 0), "nifH_Gene_x10.6_copies_m.3"] # *** NOTE: to keep things short, you could just use:  gg$nifQuantity <- gg$"nifH_Gene_x10.6_copies_m.3"
gg[which(gg$nifH_Gene_x10.6_copies_m.3 > 0), "nifQuantityType"] <- "nifH_gene_copies_x10^6_per_m3" # *** COMMENT: I correct from "nifH_gene_copies_per_m3" to "nifH_gene_copies_x10^6_per_m3"

# g. Get hit-based abundances into target format, target column "nifQuantity" and "nifQuantityType"
gg$hitStatus<-"NA"
gg$hitQuantity<-NA
gg$hitQuantityType<-"NA"

# h. Get biomass conversion factor, and biomass conversion factor quantity type
gg$biomass_conversion_factor<-NA
gg$biomass_conversion_factor_QuantityType<-"NA"

# i. Get biomass quantity and quantity type to differentiate integrated measurements biomass_mgC_per_m2
gg$biomassQuantity<-NA #  *** Comment: added
gg$biomassQuantityType<-"NA" #  *** Comment: added

# j. Put basis of record, measurement method, and references in place
gg$basisOfRecord<-NA # *** Comment: see commons sheet. Possible to define? qPCR may correspond to OBSERVATION?
gg$measurementMethod<- paste0(gg$measurementMethod," using ",gg$Notes)  # *** Comment: I add the information on the size fraction here! Makes sense?
# gg$references # perfekt. already there

# k. Get IDs, add flag column
gg$occurrenceID <- paste0(as.character(gg$scientificName),"_",as.character(gg$family),"_", as.character(gg$genus),"_",as.character(gg$decimalLongitude),"_",   # *** Comment: I add genus. Does it make sense?
as.character(gg$decimalLatitude),"_",as.character(gg$year),"_",as.character(gg$month),"_",as.character(gg$day),"_",as.character(gg$depth))
gg$sampleID<-paste0(as.character(gg$decimalLongitude),"_",as.character(gg$decimalLatitude),"_",as.character(gg$year),"_",as.character(gg$month),"_",as.character(gg$day),"_",as.character(gg$depth))
gg$stationID <- "NA" # Comment: Not available...?
gg$cruiseID <- "NA" # Comment: Not available...?
gg$flag <- "NA"

# l. Get environmental data (if available)
# Not available, I guess (?)
gg$"Sea_surface_temperature_C"  <- NA
gg$"Temperature_C"  <- NA
gg$"Sea_surface_Salinity" <- NA
gg$"Salinity" <- NA
gg$"Sea_surface_nitrate_micromolar" <- NA
gg$"Nitrate_micromolar" <- NA
gg$"Sea_surface_phosphate_micromolar" <- NA
gg$"Phosphate_micromolar" <- NA
gg$"Sea_surface_Fe_nanomolar" <- NA
gg$"Fe_nanomolar" <- NA
gg$"ChlorophyllQuantity" <- NA
gg$"ChlorophyllQuantityType" <- NA

# m. Check if all of the desired columns were addressed, and order dataframe at same time
gg <- gg[, desired.cols] # No complaints, fine.

## n. Check if/where occurrence status has been assigned
unique(gg$occurrenceStatus) # *** Comment: There are some NA's, let us assign these, as well. This is a high priority task, I would say.
nrow(gg[which(gg$occurrenceStatus=="NA"),]) # 615 NA's
nrow(gg[which(gg$occurrenceStatus=="PRESENT"),]) # 73 "PRESENT"
head(gg[which(gg$occurrenceStatus=="NA"),]) # *** Comment: I take a look at the NA cases. Please verfiy with original dataset
unique(gg[which(gg$occurrenceStatus=="NA"),"measurementMethod"]) # *** Comment: qPCR based

## o. add set source
gg$set <- "Gra"

# PS: check availability of date
dim(gg[is.na(gg$month), ]) # 0, fine
unique(gg$day) #  NA
unique(gg$month) # 4 or 5, needs to be specified or verified if possible. As unclear in raw data.
unique(gg$year) # Ok, looks fine
##---------------------------------------------------------------------------------------------

# 4. GBIF
dat.gbif<-read.csv("../1_Download_data/GBIF_2020/Copies_to_work_with/Collect_raw_data(gbif)_version2020_DE_copy_raw.csv")
dim(dat.gbif)# 11963 x 175
dat.gbif$year <- as.numeric(as.character(dat.gbif$year)) # *** comment: added
dat.gbif$month <- as.numeric(as.character(dat.gbif$month)) # *** comment: added
dat.gbif$day <- as.numeric(as.character(dat.gbif$day)) # *** comment: added
dat.gbif$depth <- as.numeric(as.character(dat.gbif$depth)) # *** comment: added
dat.gbif$occurrenceID<-paste0(as.character(dat.gbif$scientificName),"_",as.character(dat.gbif$family),"_",as.character(dat.gbif$genus),"_",as.character(dat.gbif$decimalLongitude),"_",
as.character(dat.gbif$decimalLatitude),"_",as.character(dat.gbif$year),"_",as.character(dat.gbif$month),"_",as.character(dat.gbif$day),"_",as.character(dat.gbif$depth))
dat.gbif$sampleID<-paste0(as.character(dat.gbif$decimalLongitude),"_",as.character(dat.gbif$decimalLatitude),"_",as.character(dat.gbif$year),"_",as.character(dat.gbif$month),"_",as.character(dat.gbif$day),"_",as.character(dat.gbif$depth))
dat.gbif$basisOfRecord <- as.character(dat.gbif$basisOfRecord) # *** Comment: needs to be explicitly defined as character, else the information is lost in the big merger below, as other sets are defined as factor
dat.gbif$set <- "Gbi"
# *** Comment: please retain key information on dataset or cruise ID, cell or individual counts, and perhaps other potentially useful data-fields (e.g., quickly check in PhytoBase paper, table 1)

# PS: check availability of date
dim(dat.gbif[is.na(dat.gbif$month), ]) # 547, some..
unique(dat.gbif$day) # Ok, looks fine, one NA, though
unique(dat.gbif$month) # Ok, looks fine, one NA, tough
unique(dat.gbif$year) # Ok, looks fine
##------------------------------------------------------------------------------------------------

# 5. OBIS
dat.obis<-read.csv("../1_Download_data/OBIS_2020/Copies_to_work_with/Collect_raw_data(obis)_version2020_DE_copy_raw.csv")
dim(dat.obis)# 5616 147

# --- inset: get year, month, day ---
# A. Get the sampling date via the "date_start" field: does currently not work as raw data are corrupted regarding this data field
# B. Get the sampling date via the "eventDate" field, this is a quick correction that helps for most (but perhaps a more careful check would be valid)
dat.obis$eventDate <- as.character(dat.obis$eventDate) # adjust class of date-column to make it workable. the format is DAY-MONTH-YEAR
dates.with.slash <- gsub("/", "-", x = dat.obis$eventDate [  grep("/", x = dat.obis$eventDate)]  ) # adjusting all dates with separator "-"
# adjust special cases and/or exceptions manually that do not have the correct day-month-year order (unlike the majority of the date entries)
dates.with.slash[which(dates.with.slash == "2012-06-06-2012-06-09" )] <- "NA-06-2012"
dates.with.slash[which(dates.with.slash == "1992-06-18-1992-07-07" )] <- "NA-NA-1992"
dates.with.slash[which(dates.with.slash == "1996-04-06-1998-06-20" )] <- "NA-NA-1996"
dates.with.slash[which(dates.with.slash == "1996-04-06-1998-06-20" )] <- "NA-NA-1996"
dates.with.slash[which(dates.with.slash == "1985-01-1985-05" )] <- "NA-NA-1985"
dates.with.slash[which(dates.with.slash == "1985-01-1985-05" )] <- "NA-NA-1985"
dates.with.slash[which(dates.with.slash == "2012-06-06-2012-06-09" )] <- "NA-06-2012"
dates.with.slash[which(dates.with.slash == "1968-1972-1968" )] <- "NA-NA-NA"
dates.with.slash # Looks good! No detailed check was done, though.
df.thereof <- data.frame(do.call("rbind", strsplit(dates.with.slash, "-", fixed = T)))
reassembled.dates <- paste0(   lapply(strsplit(as.character(df.thereof[,3])," ", fixed = T), '[[', 1), "-",df.thereof[,2], "-",df.thereof[,1]) # assemble to desired format
dat.obis[ grep("/", x = dat.obis$eventDate) , "eventDate"] <- reassembled.dates # replace the dates-to-be adjusted in dat.obis
# continue with the cleaning
df.date <- data.frame(do.call("rbind", strsplit(dat.obis$eventDate, "-", fixed = T))) # day, month, year
reassembled.dates.final <- paste0(  df.date[,1], "-" ,df.date[,2], "-" , lapply(strsplit(as.character(df.date[,3]),"T", fixed = T), '[[', 1)  )
# adjust special cases and/or exceptions: note that this was a quick check, there are certainly some more mistakes in the data
reassembled.dates.final[which(reassembled.dates.final == "2008-05-2008" )] <- "2008-05-NA"
df.dat <- data.frame(do.call("rbind", strsplit(reassembled.dates.final, "-", fixed = T)))
colnames(df.dat) <- c("day", "month", "year")
# --- end inset ---

# fill in the date attributes (year, month, day) from the just established event date data.frame - where missing
dat.obis[which(is.na(dat.obis$year)), "year"]<-df.dat[which(is.na(dat.obis$year)),1] # *** comment: added
dat.obis[which(is.na(dat.obis$month)), "month"]<-df.dat[which(is.na(dat.obis$month)),2] # *** comment: added
dat.obis[which(is.na(dat.obis$day)), "day"]<-df.dat[which(is.na(dat.obis$day)),3] # *** comment: added
dat.obis$year <- as.numeric(as.character(dat.obis$year)) # *** comment: added
dat.obis$month <- as.numeric(as.character(dat.obis$month)) # *** comment: added
dat.obis$day <- as.numeric(as.character(dat.obis$day)) # *** comment: added
dat.obis$depth <- as.numeric(as.character(dat.obis$depth)) # *** comment: added
dat.obis$occurrenceID<-paste0(as.character(dat.obis$scientificName),"_",as.character(dat.obis$family),"_",as.character(dat.obis$genus),"_",as.character(dat.obis$decimalLongitude),"_",
as.character(dat.obis$decimalLatitude),"_",as.character(dat.obis$year),"_",as.character(dat.obis$month),"_",as.character(dat.obis$day),"_",as.character(dat.obis$depth))
dat.obis$sampleID<-paste0(as.character(dat.obis$decimalLongitude),"_",as.character(dat.obis$decimalLatitude),"_",as.character(dat.obis$year),"_",as.character(dat.obis$month),"_",as.character(dat.obis$day),"_",as.character(dat.obis$depth))
dat.obis$basisOfRecord <- as.character(dat.obis$basisOfRecord) # *** Comment: needs to be explicitly defined as character, else the information is lost in the big merger below, as other sets are defined as factor
dat.obis$set <- "Obi"
# *** Comment: perhaps retain key information on dataset or cruise ID, cell or individual counts, and perhaps other potentially useful data-fields (e.g., quickly check in PhytoBase paper, table 1)

# PS: check availability of date
dim(dat.obis[is.na(dat.obis$month), ]) # 0, perfect (note: this was 2776 before!!! without inferring the date from the eventDate field. So, this is excellent progress)
unique(dat.obis$day) # Ok, looks fine, one NA, though, yet some day numbers are too large!
dat.obis[which(dat.obis$day > 31), "day"] <- NA # Adjustment
unique(dat.obis$month) # Ok, looks fine
unique(dat.obis$year) # Comment: Ok, this does not look good... some years are out of range *** further claning/checking needed***
##-------------------------------------------------------------------------------------------------

# 6.Phytobase
dat.phyto<-read.csv("../1_Download_data/Righetti_et_al_2020/Output_2/PhytoBase_subset.csv")
dim(dat.phyto)# 3091 41
dat.phyto$year <- as.numeric(as.character(dat.phyto$year)) # *** comment: added
dat.phyto$month <- as.numeric(as.character(dat.phyto$month)) # *** comment: added
dat.phyto$day <- as.numeric(as.character(dat.phyto$day)) # *** comment: added
dat.phyto$depth <- as.numeric(as.character(dat.phyto$depth)) # *** comment: added
dat.phyto$references<-dat.phyto$originDatabase_maredat
dat.phyto$cellStatus<-"NA"
dat.phyto[which(dat.phyto$cellQuantity > 0), "cellStatus" ]<-"PRESENT"
dat.phyto$cellQuantity<-dat.phyto$organismQuantity
dat.phyto$cellQuantityType<-dat.phyto$organismQuantityType
dat.phyto$occurrenceID<-paste0(as.character(dat.phyto$scientificName),"_",as.character(dat.phyto$family),"_",as.character(dat.phyto$genus),"_",as.character(dat.phyto$decimalLongitude),"_",
as.character(dat.phyto$decimalLatitude),"_",as.character(dat.phyto$year),"_",as.character(dat.phyto$month),"_",as.character(dat.phyto$day),"_",as.character(dat.phyto$depth))
dat.phyto$sampleID<-paste0(as.character(dat.phyto$decimalLongitude),"_",as.character(dat.phyto$decimalLatitude),"_",as.character(dat.phyto$year),"_",as.character(dat.phyto$month),"_",as.character(dat.phyto$day),"_",as.character(dat.phyto$depth))
dat.phyto$basisOfRecord <- as.character(dat.phyto$basisOfRecord) # *** Comment: needs to be explicitly defined as character, else the information is lost in the big merger below, as other sets are defined as factor
dat.phyto$set <- "Phy"
# *** Comment: please retain key information on dataset or cruise ID, cell or individual counts, and perhaps other potentially useful data-fields (e.g., quickly check in PhytoBase paper, table 1)

# PS: check availability of date
dim(dat.phyto[is.na(dat.phyto$month), ]) # 0, fine
unique(dat.phyto$day) # Ok, looks fine
unique(dat.phyto$month) # Ok, looks fine, one NA, tough
unique(dat.phyto$year) # Ok, looks fine
##-------------------------------------------------------------------------------------------------

# 7. Tara (OTU based, obtained via Mgnify)
dat.tar<-read.csv("../1_Download_data/MGnify/Tara Ocean/Output_4/Tara_OTUbased_merged.csv")
dim(dat.tar)# 3669 x 95, among which only 45 records belong to diazotrophs, and the rest to potential host taxa

# *** COMMENT: This is a quick and preliminary file load, if possible check through the columns systematically, to get what can be used... as for luo
dat.tar$decimalLongitude <- dat.tar$"geographic.location..longitude."
dat.tar$decimalLatitude <- dat.tar$"geographic.location..latitude."
# A. Get the sampling date via the "event date time start" field
df.time.start <- data.frame(do.call("rbind", strsplit(as.character (dat.tar$event.date.time.start),"-" ,fixed =T)  )) # *** Comment: please check, no complete time information available? there are many NAs
# df.time.end <- data.frame(do.call("rbind", strsplit(as.character (dat.tar$event.date.time.end),"-" ,fixed =T)  )) # *** Comment: no better coverage with respect to date
dat.tar$year <- as.character(df.time.start[,1]) # *** Comment: note, even though year should be numeric in the end, operations with replacing vlaues are much more savely and better performed as character
dat.tar$month <- as.character(df.time.start[,2])
dat.tar$day <- as.character(data.frame(do.call("rbind", strsplit(as.character (df.time.start[,3]),"T" ,fixed =T)  ))[,1])
# B. Get the sampling date via the "collection date" field
df.date <- data.frame(do.call("rbind", strsplit(as.character (dat.tar$collection.date),"-" ,fixed =T)  )) # *** Comment: please check, no complete time information available? there are many NAs
dat.tar[which(is.na(dat.tar$year)), "year"] <- as.character(df.date[which(is.na(dat.tar$year)),1])
dat.tar[which(is.na(dat.tar$month)), "month"]<-df.date[which(is.na(dat.tar$month)),2]
dat.tar[which(is.na(dat.tar$day)), "day"]<-df.date[which(is.na(dat.tar$day)),3]
# adjust class
dat.tar$year <- as.numeric(as.character(dat.tar$year))
dat.tar$month <- as.numeric(as.character(dat.tar$month))
dat.tar$day <- as.numeric(as.character(dat.tar$day))
# add measurement method: *** comment: ideally done before in the set composition of the raw dataset
dat.tar$scientificName <- dat.tar$species.name
unique(dat.tar$scientificName) # fine
dat.tar$occurrenceID<-paste0(as.character(dat.tar$scientificName),"_",as.character(dat.tar$family),"_",as.character(dat.tar$genus),"_",as.character(dat.tar$decimalLongitude),"_",
as.character(dat.tar$decimalLatitude),"_",as.character(dat.tar$year),"_",as.character(dat.tar$month),"_",as.character(dat.tar$day),"_",as.character(dat.tar$depth))
dat.tar$sampleID<-paste0(as.character(dat.tar$decimalLongitude),"_",as.character(dat.tar$decimalLatitude),"_",as.character(dat.tar$year),"_",as.character(dat.tar$month),"_",as.character(dat.tar$day),"_",as.character(dat.tar$depth))
# *** COMMENT: There is also a column called SampleID in the dataset, we should retain this TARA SAMPLE ID, probably (e.g, adding a column, originalSampleID ? -> add to commons sheet)
dat.tar$occurrenceStatus <- "PRESENT" # *** COMMENT: This is an assumption. Please verify! Are all the record confirmed presences...? Do they include quantitative statements..?
dat.tar$measurementMethod <-  as.character(dat.tar$"sample.collection.device") # *** Comment: Necessary to change to character! Else with factors the below replacement does not work. if possible add reference or details on method (otu sequencing based on v...) I think there were two different analyses methods used?
head(dat.tar[is.na(dat.tar$measurementMethod), ]) # *** Comment: to check
dat.tar[ which(is.na(dat.tar$measurementMethod)), "measurementMethod" ] <- "SSU_rRNA_detection" # *** COMMENT: IF lacking, I insert a method here, manually as "SSU_rRNA_detection"
head(dat.tar[is.na(dat.tar$measurementMethod), ]) # *** Comment: None, fine!
dat.tar$set <- "TarOTU"

# PS: check availability of date
dim(dat.tar[is.na(dat.tar$month), ]) # 92, only some.. (this information on the month of collection, in theory, should be available...)
unique(dat.tar$day) # Ok, looks fine, but some day numbers are too large! *** Further cleaning/checks needed ***. Something seems wrong..
dat.tar[which(dat.tar$day > 31), "day"] <- NA # Adjustment: may be wrong to do, here, quick implementation. Check original day entries.
unique(dat.tar$month) # Ok, looks fine, one NA, tough, which should be possible to obtain!
unique(dat.tar$year) # Ok, looks fine, but one NA, should be possible to obtain!
##-------------------------------------------------------------------------------------------------

# 8. Martinez-Perez (cell count and nifH based) # *** COMMENT: NEW
mm <-read.csv("../1_Download_data/Martinez_et_al_2016/Output_1/Martinez_reformat_raw.csv")
dim(mm)# 504 x 95

# a. Get longitude, latitude, year, month, day, depth
mm$decimalLongitude<-mm$"LONGITUDE"
mm$decimalLatitude<-mm$"LATITUDE"
class(mm$year)
class(mm$month)
class(mm$day)
class(mm$depth)
mm$year <- as.numeric(as.character(mm$year))
mm$month <- as.numeric(as.character(mm$month))
mm$day <- as.numeric(as.character(mm$day))
mm$depth<-as.numeric(as.character(mm$depth))
mm$depthAccuracy<-NA
mm$depthIntegral<-NA

mm$family <- "NA" # *** Comment: added here, ideally implemented/completed in raw file!

# b. Get lifeForm from specific dataset (perhaps, rathern than from the generic Taxon_table file) and associated taxa
mm$lifeForm<-"NA" # Is there some information?
mm$associatedTaxa <- "NA" # Is there some information?

# c. Get occurrence status
# Already set in raw file

# d. Get trichome-based abundances into target format, target column "trichomeQuantity" and "trichomeQuantityType"
mm$trichomeStatus<-"NA"
mm$trichomeQuantity<-NA
mm$trichomeQuantityType<-"NA"

# e. Get cell-based abundances into target format, target column "cellQuantity" and "cellQuantityType"
unique(mm$cellQuantityType) # native volume = per litre
mm$cellQuantity <- (1000*(as.numeric(as.character(mm$cellQuantity))))/(10^6) # scale up from cells per litre to cells per m^3, by multiplying with 1000,  and divide by 1'000'000 to get to target format
mm$cellQuantityType <- as.character(mm$cellQuantityType)
mm[which(mm$cellQuantityType == "cells_l"), "cellQuantityType"]<- "Cells_x10^6_per_m3" # Adjust the unit (cellQuantityType)

# f. Get gene-based abundances into target format, target column "nifQuantity" and "nifQuantityType"
unique(mm$nifQuantityType) # native volume = per litre
mm$nifQuantity <- (1000*(as.numeric(as.character(mm$nifQuantity))))/(10^6) # scale up from cells per litre to cells per m^3, by multiplying with 1000,  and divide by 1'000'000 to get to target format
mm$nifQuantityType <- as.character(mm$nifQuantityType)
mm[which(mm$nifQuantityType == "nifH_copies_l"), "nifQuantityType"]<- "nifH_gene_copies_x10^6_per_m3" # Adjust the unit (cellQuantityType)

# g. Get hit-based abundances into target format, target column "nifQuantity" and "nifQuantityType"
mm$hitStatus<-"NA"
mm$hitQuantity<-NA
mm$hitQuantityType<-"NA"

# h. Get biomass conversion factor, and biomass conversion factor quantity type
mm$biomass_conversion_factor<-NA
mm$biomass_conversion_factor_QuantityType<-"NA"

# i. Get biomass quantity and quantity type to differentiate integrated measurements biomass_mgC_per_m2
mm$biomassQuantity<-NA
mm$biomassQuantityType<-"NA"

# j. Put basis of record, measurement method, and references in place
mm$basisOfRecord<-NA # *** Comment: see commons sheet. Possible to define? Epifluorescence microscopy may correspond to OBSERVATION, qpcr as well
# mm$references # already in raw file

# k. Get IDs, add flag column
mm$occurrenceID <- paste0(as.character(mm$scientificName),"_",as.character(mm$family),"_", as.character(mm$genus),"_",as.character(mm$decimalLongitude),"_",   # *** Comment: I add genus. Does it make sense? NOTE: "family" is excluded below, as it is not in the target columns!!!! I hence add family to the target columns -> correct?
as.character(mm$decimalLatitude),"_",as.character(mm$year),"_",as.character(mm$month),"_",as.character(mm$day),"_",as.character(mm$depth))
mm$sampleID<-paste0(as.character(mm$decimalLongitude),"_",as.character(mm$decimalLatitude),"_",as.character(mm$year),"_",as.character(mm$month),"_",as.character(mm$day),"_",as.character(mm$depth))
mm$stationID<-NA # Perhaps can be assigned to specific stations?
mm$cruiseID<-"M96" # Perhaps can be specified
mm$flag<-NA

# l. Get environmental data (if available): *** implement this carefully *** currently set to NA
mm$ChlorophyllQuantity<-NA
mm$ChlorophyllQuantityType<-"NA"
mm$Sea_surface_temperature_C<-NA
mm$Temperature_C<-NA
mm$Sea_surface_Salinity<-NA
mm$Salinity<-NA
### Include Environmental data on "PO43.", "NO3.",  "NO2." , and   "NH4." are there. I quickly include the data here - but the units need a more careful check if this is correct. and if it is sea surface or other measurement...
mm$Sea_surface_nitrate_micromolar<- as.numeric(as.character(mm$"NO3."))/1000 # Assumption nanomolar needs to be divided by 1000 to correspond to micromolar
mm$Nitrate_micromolar<-NA
mm$Sea_surface_phosphate_micromolar<-  as.numeric(as.character(mm$"PO43."))/1000 # Assumption nanomolar needs to be divided by 1000 to correspond to micromolar
mm$Phosphate_micromolar<-NA
### Other environmental data
mm$Sea_surface_Fe_nanomolar<-NA
mm$Fe_nanomolar<-NA

## m. Check if all of the desired columns were addressed, and order dataframe at same time
mm <- mm[, desired.cols] # No complaints, fine.

## n. Check if/where occurrence status has been assigned
unique(mm$occurrenceStatus) # looks fine
unique(mm[which(mm$occurrenceStatus=="NA"),"measurementMethod"]) # looks fine

## o. add set source
mm$set <- "Mar"

# PS: check availability of date
dim(mm[is.na(mm$month), ]) # 0, fine
unique(mm$day) # NA
unique(mm$month) # Looks fine
unique(mm$year) # Looks fine
# plot(mm$decimalLatitude, mm$decimalLongitude, xlim = c(-180,180), ylim = c(-90,90)) # An atlantic meridional cruise on the Southern Hemisphere. This is valuable!

### ===================================================================================
### 3. Merge the source datasets, add lifeForm information, remove NAs & duplicates, save file
### ===================================================================================

dat <- rbind.fill( ll , tt , gg , dat.gbif , dat.obis , dat.phyto , dat.tar , mm)
dat <- dat [,c(desired.cols, "set")] # *** NOTE: I add "setOrigin" to the desired columns to keep trace on our sources

# Add column lifeForm: *** COMMENT: I would probably not generalize, but just keep relevant information from specific sources (second order priority). Else, this is no information gain relative to the taxon table
dat[which(dat$scientificName == "Calothrix"), "lifeForm"] <- as.character(tax.table[which(tax.table$scientificName_used_in_dataset == "Calothrix"), "Singular_or_colony_forming_reference" ])
dat[which(dat$scientificName == "Richelia"), "lifeForm"] <- as.character(tax.table[which(tax.table$scientificName_used_in_dataset == "Richelia"), "Singular_or_colony_forming_reference" ])
dat[which(dat$scientificName == "Richelia intracellularis"), "lifeForm"] <- as.character(tax.table[which(tax.table$scientificName_used_in_dataset == "Richelia intracellularis"), "Singular_or_colony_forming_reference" ])
dat[which(dat$scientificName == "Trichodesmium"), "lifeForm"] <- as.character(tax.table[which(tax.table$scientificName_used_in_dataset == "Trichodesmium"), "Singular_or_colony_forming_reference" ])
dat[which(dat$scientificName == "Trichodesmium contortum "), "lifeForm"] <- as.character(tax.table[which(tax.table$scientificName_used_in_dataset == "Trichodesmium contortum "), "Singular_or_colony_forming_reference" ])
dat[which(dat$scientificName == "Trichodesmium erythraeum"), "lifeForm"] <- as.character(tax.table[which(tax.table$scientificName_used_in_dataset == "Trichodesmium erythraeum"), "Singular_or_colony_forming_reference" ])
dat[which(dat$scientificName == "Trichodesmium lacustre"), "lifeForm"] <- as.character(tax.table[which(tax.table$scientificName_used_in_dataset == "Trichodesmium lacustre"), "Singular_or_colony_forming_reference" ])
dat[which(dat$scientificName == "Trichodesmium thiebautii"), "lifeForm"] <- as.character(tax.table[which(tax.table$scientificName_used_in_dataset == "Trichodesmium thiebautii"), "Singular_or_colony_forming_reference" ])
dat[which(dat$scientificName == "Atelocyanobacterium"), "lifeForm"] <- as.character(tax.table[which(tax.table$scientificName_used_in_dataset == "Atelocyanobacterium"), "Singular_or_colony_forming_reference" ])
dat[which(dat$scientificName == "UCYN.A"), "lifeForm"] <- as.character(tax.table[which(tax.table$scientificName_used_in_dataset == "UCYN.A"), "Singular_or_colony_forming_reference" ])
dat[which(dat$scientificName == "UCYN.A1"), "lifeForm"] <- as.character(tax.table[which(tax.table$scientificName_used_in_dataset == "UCYN.A1"), "Singular_or_colony_forming_reference" ])
dat[which(dat$scientificName == "UCYN.A2"), "lifeForm"] <- as.character(tax.table[which(tax.table$scientificName_used_in_dataset == "Trichodesmium lacustre"), "Singular_or_colony_forming_reference" ])
dat[which(dat$scientificName == "UCYN.B"), "lifeForm"] <- as.character(tax.table[which(tax.table$scientificName_used_in_dataset == "UCYN.B"), "Singular_or_colony_forming_reference" ])
dat[which(dat$scientificName == "UCYN.C"), "lifeForm"] <- as.character(tax.table[which(tax.table$scientificName_used_in_dataset == "UCYN.C"), "Singular_or_colony_forming_reference" ])
dat[which(dat$scientificName == "Gamma.A"), "lifeForm"] <- as.character(tax.table[which(tax.table$scientificName_used_in_dataset == "Gamma.A"), "Singular_or_colony_forming_reference" ])
dat[which(dat$scientificName == "Azotobacter salinestris"), "lifeForm"] <- as.character(tax.table[which(tax.table$scientificName_used_in_dataset == "Azotobacter salinestris"), "Singular_or_colony_forming_reference" ])

# Remove NA's using the occurrenceStatus information: *** Comment -> harmonize/clean the occurrence status information?
unique(dat$occurrenceStatus)
#[1] "PRESENT"                   "NA"                        NA
#[4] "Presente"                  "present"                   "Present"
#[7] "Desconocido"               "Pr�sent"                   "Present (1% <= p < 5%)"
#[10] "Rare (p < 1%)"             "Common (5% <= p < 10%)"    "Abundant (10% <= p < 20%)"
#[13] "Dominant (20% <= p)"       "absent"
dat[which(dat$occurrenceStatus=="Presente"), "occurrenceStatus"] <- "PRESENT"
dat[which(dat$occurrenceStatus=="present"), "occurrenceStatus"] <- "PRESENT"
dat[which(dat$occurrenceStatus=="Present"), "occurrenceStatus"] <- "PRESENT"
dat[which(dat$occurrenceStatus=="Pr�sent"), "occurrenceStatus"] <- "PRESENT"
dat[which(dat$occurrenceStatus=="Pr\351sent"), "occurrenceStatus"] <- "PRESENT"
dat[which(dat$occurrenceStatus=="absent"), "occurrenceStatus"] <- "ABSENT"

# Let us move such information into the notes/flag column for the moment
dat[which(dat$occurrenceStatus=="Present (1% <= p < 5%)"), "occurrenceStatus"] <- "PRESENT"
dat[which(dat$occurrenceStatus=="Present (1% <= p < 5%)"), "flag"] <- "Present (1% <= p < 5%)"
dat[which(dat$occurrenceStatus=="Common (5% <= p < 10%)"), "occurrenceStatus"] <- "PRESENT"
dat[which(dat$occurrenceStatus=="Common (5% <= p < 10%)"), "flag"] <- "Common (5% <= p < 10%)"
dat[which(dat$occurrenceStatus=="Abundant (10% <= p < 20%)"), "occurrenceStatus"] <- "PRESENT"
dat[which(dat$occurrenceStatus=="Abundant (10% <= p < 20%)"), "flag"] <- "Abundant (10% <= p < 20%)"
dat[which(dat$occurrenceStatus=="Dominant (20% <= p)"), "occurrenceStatus"] <- "PRESENT"
dat[which(dat$occurrenceStatus=="Dominant (20% <= p)"), "flag"] <- "Dominant (20% <= p)"
unique(dat$occurrenceStatus)
# "PRESENT"       "NA"            NA              "Desconocido"   "Rare (p < 1%)" "ABSENT"

# Convert "NA" to NA since the first is not recognized as a Na entry: *** Comment: I would leave as is.. or is this a problem?
dat[dat=="NA"] = NA
dat.final<- dat[  -which( is.na(dat$occurrenceStatus) ),]
dim(dat.final)# 26308 x 55

# Save removed/extracted(?) NAs in an additional dataframe: *** Comment -> we will gain substantially more data by specifiying this information ->  see above
na.extracted.subset<- dat[  which( is.na(dat$occurrenceStatus) ),]
dim(na.extracted.subset)# 18977 x 55

# Change class of columns
dat.final$scientificName <- as.character(dat.final$scientificName)

# Adjust certain scientific names: *** Comment: IDEALLY implement already at level datasource (I guess these author names come from gbif and obis only)
unique(dat.final$scientificName)
dat.final[which(dat.final$scientificName=="Richelia intracellularis J.Schmidt"), "scientificName"] <- "Richelia"
dat.final[which(dat.final$scientificName=="Richelia J.Schmidt, 1901"), "scientificName"] <- "Richelia"
dat.final[which(dat.final$scientificName=="Trichodesmium erythraeum Ehrenb. ex Gomont"), "scientificName"] <- "Trichodesmium"
dat.final[which(dat.final$scientificName=="Trichodesmium thiebautii Gomont ex Gomont"), "scientificName"] <- "Trichodesmium thiebautii"
dat.final[which(dat.final$scientificName=="Calothrix C.Agardh ex Bornet & Flahault, 1886"), "scientificName"] <- "Calothrix"
dat.final[which(dat.final$scientificName=="Calothrix scopulorum C.Agardh ex Bornet & Flahault"), "scientificName"] <- "Calothrix scopulorum"
dat.final[which(dat.final$scientificName=="Calothrix scopulorum (Weber & Mohr) C.Agardh, 1824"), "scientificName"] <- "Calothrix scopulorum"
dat.final[which(dat.final$scientificName=="Calothrix confervicola C.A.Agardh Ex Born. & Flah"), "scientificName"] <- "Calothrix confervicola"
dat.final[which(dat.final$scientificName=="Calothrix fusca Born. & Flah."), "scientificName"] <- "Calothrix fusca"
dat.final[which(dat.final$scientificName=="Calothrix stagnalis Gomont, 1895"), "scientificName"] <- "Calothrix stagnalis"
dat.final[which(dat.final$scientificName=="Calothrix brevissima G.S.West"), "scientificName"] <- "Calothrix brevissima"
dat.final[which(dat.final$scientificName=="Calothrix braunii Bornet & Flahault" ), "scientificName"] <- "Calothrix braunii"
dat.final[which(dat.final$scientificName=="Calothrix nodulosa Setchell & N.L.Gardner, 1924"), "scientificName"] <- "Calothrix nodulosa"
dat.final[which(dat.final$scientificName=="Calothrix crustacea Thuret"), "scientificName"] <- "Calothrix crustacea"
dat.final[which(dat.final$scientificName=="Calothrix parietina Thur. ex Bornet & Flahault"), "scientificName"] <- "Calothrix parietina"
dat.final[which(dat.final$scientificName== "Calothrix parasitica (Chauv.) Thur. ex Bornet & Flahault"), "scientificName"] <-  "Calothrix parasitica"
dat.final[which(dat.final$scientificName=="Calothrix confervicola (Roth) C.Agardh"), "scientificName"] <- "Calothrix confervicola"
dat.final[which(dat.final$scientificName=="Calothrix pulvinata Ag."), "scientificName"] <- "Calothrix pulvinata"
dat.final[which(dat.final$scientificName== "Calothrix robusta Setch. & N.L.Gardner" ), "scientificName"] <- "Calothrix robusta"
dat.final[which(dat.final$scientificName=="Calothrix pulvinata C.Agardh ex Bornet & Flahault"), "scientificName"] <- "Calothrix pulvinata"
dat.final[which(dat.final$scientificName== "Trichodesmium Ehrenberg ex Gomont, 1892"), "scientificName"] <- "Trichodesmium"
dat.final[which(dat.final$scientificName=="Trichodesmium hildebrandtii Gomont"), "scientificName"] <- "Trichodesmium hildebrandtii"
dat.final[which(dat.final$scientificName=="Cyanothece majus (Schr\366t.) Kom\341rek"), "scientificName"] <- "Cyanothece majus"
dat.final[which(dat.final$scientificName=="Cyanothece Kom\341rek, 1976"), "scientificName"] <- "Cyanothece"
dat.final[which(dat.final$scientificName=="Heterocyst.Richelia.Calothris"), "scientificName"] <- "Heterocyst.Richelia.Calothrix"  # *** Comment an "x" instead of an "s". Please correct in the raw source (i.e., Tang and Cassar)
unique(dat.final$scientificName)

## *** Comment: Check if all records have geographic coordinates (the cleaning of records without coordinates is now implemented here, but it should be implemented in a next-step cleaning file)
unique(dat.final$decimalLongitude) # Note: There are some NAs - what happended... Let us check in the raw source, i.e., ... Tara OTU. Did something go wrong with matching up of flagships with sample information on lat lon?
dim(dat.final[is.na(dat.final$decimalLongitude), ]) # There are 23 NA-entries
unique(dat.final$decimalLatitude) # Note: There are some NAs - what happended... Let us check in the raw source, i.e., ... Tara OTU. Did something go wrong with matching up of flagships with sample information on lat lon?
dim(dat.final[is.na(dat.final$decimalLatitude), ]) # There are 23 NA-entries
dat.final <- dat.final[!is.na(dat.final$decimalLatitude), ] # Excluding all records that have NA-coordinate entries; note, hoewer, that this cleaning step may be implemented together with all other cleaning steps in a next script.

## *** Comment: Check if all records have a month of collection indication (the cleaning of records without coordinates is now implemented here, but it should be implemented in a next-step cleaning file)
unique(dat.final$month) # Note: There are some NAs -
dim(dat.final[is.na(dat.final$month), ]) # There are 116 NA entries...
dat.final[is.na(dat.final$month), "set"] # These records stem largely from GBIF and taraOTU. Note that dates were not checked carefully yet!!!
dat.final <- dat.final[!is.na(dat.final$month), ] # Excluding all records that have lacking month entries; note, hoewer, that this cleaning step may be implemented together with all other cleaning steps in a next script

## Harmonize the measurement method, if it makes sense *** Comment: inserted -> ideally check carefully in each source dataset beforehands
unique(dat.final$measurementMethod) # Note there are NA's! These need to be filled in.
# [1] "Standard Light Microscopy"                           "Epifluorescence Microscopy "                         "Epiflourescence Microscopy "                         "Microscope, transmitted light or epifluorescence"
# [5] "Direct Biomass Analysis"                             "qPCR"                                                "qPCR using Size fraction filtered for <3 micrometer" "qPCR using Size fraction filtered for >3 micrometer"
# [9] NA                                                    "PUMP (High Volume Peristaltic Pump)"
# I adjust typos and homogenise the method entries, where useful
dat.final[which(dat.final$measurementMethod%in%c("Epifluorescence Microscopy ", "Epiflourescence Microscopy ", "Epifluorescence Microscopy")), "measurementMethod"] <- "Epifluorescence microscopy"

dat.final[which(dat.final$measurementMethod == "Microscope, transmitted light or epifluorescence"), "measurementMethod"] <-  "Standard light or epifluorescence microscopy"
dat.final[which(dat.final$measurementMethod%in%c("qPCR",   "qPCR using Size fraction filtered for <3 micrometer", "qPCR using Size fraction filtered for >3 micrometer")), "measurementMethod"] <-  "qPCR_nifH_detection" # Add information on size class in another column, eg flag or notes
dat.final[which(dat.final$measurementMethod == "PUMP (High Volume Peristaltic Pump)" & dat.final$set=="TarOTU"), "measurementMethod"] <- "SSU_rRNA_detection" # For Tara Ocean
dat.final[which(dat.final$measurementMethod == "Standard Light Microscopy"), "measurementMethod"] <-  "Standard light microscopy"
dat.final[which(dat.final$measurementMethod == "Direct Biomass Analysis"), "measurementMethod"] <-  "Direct biomass analysis"
# Check again
unique(dat.final$measurementMethod) # *** COMMENT: Note there are NA's. If possible, let us fill in these
# check the origin of the NAs
unique(dat.final[is.na(dat.final$measurementMethod), "set"]) # GBIF, OBIS, PhytoBase
# Check basis of record for these
dat.final[which(dat.final$set=="Gbi"), "basisOfRecord"]
dat.final[which(dat.final$set=="Obi"), "basisOfRecord"]
dat.final[which(dat.final$set=="Phy"), "basisOfRecord"]
# probably the method for most of these is traditional light microscopy, and some via flow cytometry, but we do not have any specific evidence at the moment

## Check for and remove duplicates
## *** Comment: when removing duplicates, first a backbone of occurrenceIDs needs to be established, together with information on the set, method, or any other information of interest. Then, after removing duplicates, such information can be re-assigned to the remaining records,
## *** For example: one and a same record may stem from two sources (tang and luo), then the information on the "set" of the remaining point (e.g., luo) needs to change from "Luo" to "Tan_Luo" (see skripts of ESSD dataset, for specific example how this is done, or just ask Damiano directly)
dim(dat.final[duplicated(dat.final$occurrenceID),]) # 4215 x 55, be carefull --> taxonomic resolution when including family in the occurrenceID --> loss of data? # *** Comment: I added the family, for the moment being
dat.final<-dat.final[!duplicated(dat.final$occurrenceID),]
dim(dat.final) # 21 954 x 55 (state: 09 June 2020), 21 622 x 55 (state: 03 June 2020)

# Save file
file.name="./Output_1/Diazo.csv"
write.csv(dat.final,file = file.name,row.names = F)

#  ------------- PS: Take a quick look at the key statistics  -------------
#  ------------------------------------------------------------------------------

# Points per taxon
dat.spl<- split(dat.final, f = dat.final$scientificName)
df.tax.obs <- data.frame("taxon" = names(unlist(lapply(dat.spl, nrow))), "obs" = as.numeric(unlist(lapply(dat.spl, nrow))) )
df.tax.obs <- df.tax.obs[order(df.tax.obs$obs, decreasing = T),]
rownames(df.tax.obs) <- 1:nrow(df.tax.obs)
df.tax.obs

# Number of records
nrow(dat.final) # 21 954 (before: 18 075)
nrow(dat.final[which(dat.final$occurrenceStatus=="PRESENT"),]) # 21 527 (before: 17 879) not bad

# Number of taxa
length(unique(dat.final$scientificName)) # 59 taxa
length(unique(dat.final$genus)) # 10 genera
