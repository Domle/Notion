## This script formats the dataobservations from Juan Karlusich et al. 2021 on diazotrophs retrieved from Tara Ocean Samples.
## We will format all columns to a universal formating according to Darwin Standards (Righetti et al. 2019).

## Dominic Eriksson
## Environmental Physics Group
## ETH Zurich
## Switzerland

## dominic.eriksson@usys.ethz.ch,   xx-xx-2021

#################
## Preparatory ##
#################
# Set directory
setwd("C:/Users/domin/Documents/PhD/Notion/1_Downloaded_data/TaraOceans/Juan_Karlusich_nature_communications_2021/Confocal_microscopy_UVP5/S-BSST529/Copies_to_work_with/")
# load master data frame with all the correct column names according to Darwin Standard
df_master <- get(load("C:/Users/domin/Documents/PhD/Notion/commons/df_master.RData"))
col_names <- names(df_master)

################
## eHCFM Data ##
################
# load data Environmental High-Content Fluorescent Microscopy
dat <- read.table("N2fixers_eHCFM_20-180.tsv", sep = "\t", header = TRUE)
names(dat) # we use columns 1:15
dim(dat) # 961 569
# merge with master_df
library("plyr")
df <- rbind.fill(dat, df_master)

####################
## Format columns ##
####################
# location
df$decimalLongitude <- df$object_lon
df$decimalLatitude <- df$object_lat

# date
df$object_date <- as.character(df$object_date)
df$object_date <- as.Date(df$object_date, format = "%Y%m%d")
df$year <- as.numeric(format(df$object_date, format="%Y"))
df$month <- as.numeric(format(df$object_date, format = "%m"))
df$day <- as.numeric(format(df$object_date, format = "%d"))

# depth better: fill up the rest of the columns once with automatic function inserting NAs
df$depthIntegral <- paste(df$object_depth_min, df$object_depth_max, sep = "_")
# remove old depth columns

# Taxonomy
unique(df$object_annotation_category)
# [1] "multiple<Chaetoceros inter. Calothrix"
# [2] "Chaetoceros inter. Calothrix"
# [3] "Climacodium inter. Crocosphaera"
# [4] "Richelia"
# [5] "Hemiaulus inter. Richelia"
# [6] "multiple<Hemiaulus inter. Richelia"
# [7] "Rhizosolenia inter. Richelia"
# [8] "Trichodesmium"
# [9] "multiple<Trichodesmium"
# [10] "multiple<Rhizosolenia inter. Richelia"
tax_table <- read.csv("C:/Users/domin/Documents/PhD/Notion/commons/Taxon_table.csv")

# Calothrix, here R sometimes inserts a list that causes problems when saving, so check and unlist if neccassary
df$scientificName[which(df$object_annotation_category == "multiple<Chaetoceros inter. Calothrix")] <- "Calothrix"
class(df$scientificName) # character
df$taxonRank[which(df$object_annotation_category == "Calothrix")] <- tax_table[which(tax_table$scientificName_used_in_dataset == "Calothrix"), ]["taxonRank_of_scientificName"]
class(df$taxonRank) # list, unlist
df$taxonRank <- unlist(df$taxonRank)
df$phylum[which(df$object_annotation_category == "Calothrix")] <- tax_table[which(tax_table$scientificName == "Calothrix"), ]["Phylum"]
class(df$phylum)
df$phylum <- unlist(df$phylum)
df$class[which(df$object_annotation_category == "Calothrix")] <- tax_table[which(tax_table$scientificName == "Calothrix"), ]["Class"]
class(df$class)
df$class <- unlist(df$class)
df$family[which(df$object_annotation_category == "Calothrix")] <- unlist(tax_table[which(tax_table$scientificName == "Calothrix"), ]["Family"])
class(df$family)
df$family <- unlist(df$family)
df$genus[which(df$object_annotation_category == "Calothrix")] <- tax_table[which(tax_table$scientificName == "Calothrix"), ]["Genus"]
class(df$genus)
df$genus <- unlist(df$genus)
df$species[which(df$object_annotation_category == "Calothrix")] <- tax_table[which(tax_table$scientificName == "Calothrix"), ]["Species"]
class(df$species)
df$species <- unlist(df$species)

df$lifeForm[which(df$object_annotation_category == "Calothrix")] <- paste0(tax_table[which(tax_table$scientificName == "Calothrix"), ]["Singular_or_colony_forming_reference"], "_" ,tax_table[which(tax_table$scientificName == "Calothrix"), ]["In_symbiosis_or_freeliving_reference"])
df$associatedTaxa[which(df$object_annotation_category == "Calothrix")] <- "Chaetoceros"

# [2] "Chaetoceros inter. Calothrix"
# [3] "Climacodium inter. Crocosphaera"
# [4] "Richelia"
# [5] "Hemiaulus inter. Richelia"
# [6] "multiple<Hemiaulus inter. Richelia"
# [7] "Rhizosolenia inter. Richelia"
# [8] "Trichodesmium"
# [9] "multiple<Trichodesmium"
# [10] "multiple<Rhizosolenia inter. Richelia"
df$scientificName[which(df$object_annotation_category == "Chaetoceros inter. Calothrix")] <- "Calothrix"
df$taxonRank[which(df$object_annotation_category == "Calothrix")] <- tax_table[which(tax_table$scientificName_used_in_dataset == "Calothrix"), ]["taxonRank_of_scientificName"]
df$taxonRank <- unlist(df$taxonRank)
df$phylum[which(df$object_annotation_category == "Calothrix")] <- tax_table[which(tax_table$scientificName == "Calothrix"), ]["Phylum"]
df$phylum <- unlist(df$phylum)
df$class[which(df$object_annotation_category == "Calothrix")] <- tax_table[which(tax_table$scientificName == "Calothrix"), ]["Class"]
df$class <- unlist(df$class)
df$family[which(df$object_annotation_category == "Calothrix")] <- tax_table[which(tax_table$scientificName == "Calothrix"), ]["Family"]
df$family <- unlist(df$family)
df$genus[which(df$object_annotation_category == "Calothrix")] <- tax_table[which(tax_table$scientificName == "Calothrix"), ]["Genus"]
df$species[which(df$object_annotation_category == "Calothrix")] <- tax_table[which(tax_table$scientificName == "Calothrix"), ]["Species"]


df$lifeForm[which(df$object_annotation_category == "Calothrix")] <- paste0(tax_table[which(tax_table$scientificName == "Calothrix"), ]["Singular_or_colony_forming_reference"], "_" ,tax_table[which(tax_table$scientificName == "Calothrix"), ]["In_symbiosis_or_freeliving_reference"])
df$associatedTaxa[which(df$object_annotation_category == "Calothrix")] <- "Chaetoceros"

# [3] "Climacodium inter. Crocosphaera"
# [4] "Richelia"
# [5] "Hemiaulus inter. Richelia"
# [6] "multiple<Hemiaulus inter. Richelia"
# [7] "Rhizosolenia inter. Richelia"
# [8] "Trichodesmium"
# [9] "multiple<Trichodesmium"
# [10] "multiple<Rhizosolenia inter. Richelia"
df$scientificName[which(df$object_annotation_category == "Climacodium inter. Crocosphaera")] <- "Climacodium"
df$taxonRank[which(df$object_annotation_category == "Climacodium")] <- tax_table[which(tax_table$scientificName_used_in_dataset == "Climacodium"), ]["taxonRank_of_scientificName"]
df$phylum[which(df$object_annotation_category == "Climacodium")] <- tax_table[which(tax_table$scientificName == "Climacodium"), ]["Phylum"]
df$class[which(df$object_annotation_category == "Climacodium")] <- tax_table[which(tax_table$scientificName == "Climacodium"), ]["Class"]
df$family[which(df$object_annotation_category == "Climacodium")] <- tax_table[which(tax_table$scientificName == "Climacodium"), ]["Family"]
df$genus[which(df$object_annotation_category == "Climacodium")] <- tax_table[which(tax_table$scientificName == "Climacodium"), ]["Genus"]
df$species[which(df$object_annotation_category == "Climacodium")] <- tax_table[which(tax_table$scientificName == "Climacodium"), ]["Species"]
df$lifeForm[which(df$object_annotation_category == "Climacodium")] <- paste0(tax_table[which(tax_table$scientificName == "Climacodium"), ]["Singular_or_colony_forming_reference"], "_" ,tax_table[which(tax_table$scientificName == "Climacodium"), ]["In_symbiosis_or_freeliving_reference"])
df$associatedTaxa[which(df$object_annotation_category == "Climacodium")] <- "Calothrix"

# [4] "Richelia"
# [5] "Hemiaulus inter. Richelia"
# [6] "multiple<Hemiaulus inter. Richelia"
# [7] "Rhizosolenia inter. Richelia"
# [8] "Trichodesmium"
# [9] "multiple<Trichodesmium"
# [10] "multiple<Rhizosolenia inter. Richelia"
df$scientificName[which(df$object_annotation_category == "Richelia")] <- "Richelia"
df$taxonRank[which(df$object_annotation_category == "Richelia")] <- tax_table[which(tax_table$scientificName_used_in_dataset == "Richelia"), ]["taxonRank_of_scientificName"]
df$phylum[which(df$object_annotation_category == "Richelia")] <- tax_table[which(tax_table$scientificName == "Richelia"), ]["Phylum"]
df$class[which(df$object_annotation_category == "Richelia")] <- tax_table[which(tax_table$scientificName == "Richelia"), ]["Class"]
df$family[which(df$object_annotation_category == "Richelia")] <- tax_table[which(tax_table$scientificName == "Richelia"), ]["Family"]
df$genus[which(df$object_annotation_category == "Richelia")] <- tax_table[which(tax_table$scientificName == "Richelia"), ]["Genus"]
df$species[which(df$object_annotation_category == "Richelia")] <- tax_table[which(tax_table$scientificName == "Richelia"), ]["Species"]
df$lifeForm[which(df$object_annotation_category == "Richelia")] <- paste0(tax_table[which(tax_table$scientificName == "Richelia"), ]["Singular_or_colony_forming_reference"], "_" ,tax_table[which(tax_table$scientificName == "Richelia"), ]["In_symbiosis_or_freeliving_reference"])
df$associatedTaxa[which(df$object_annotation_category == "Richelia")] <- "Rhizosolenia, Hemiaulus"

# [5] "Hemiaulus inter. Richelia"
# [6] "multiple<Hemiaulus inter. Richelia"
# [7] "Rhizosolenia inter. Richelia"
# [8] "Trichodesmium"
# [9] "multiple<Trichodesmium"
df$scientificName[which(df$object_annotation_category == "Hemiaulus inter. Richelia")] <- "Hemiaulus"
df$taxonRank[which(df$object_annotation_category == "Hemiaulus")] <- tax_table[which(tax_table$scientificName_used_in_dataset == "Hemiaulus"), ]["taxonRank_of_scientificName"]
df$phylum[which(df$object_annotation_category == "Hemiaulus")] <- tax_table[which(tax_table$scientificName == "Hemiaulus"), ]["Phylum"]
df$class[which(df$object_annotation_category == "Hemiaulus")] <- tax_table[which(tax_table$scientificName == "Hemiaulus"), ]["Class"]
df$family[which(df$object_annotation_category == "Hemiaulus")] <- tax_table[which(tax_table$scientificName == "Hemiaulus"), ]["Family"]
df$genus[which(df$object_annotation_category == "Hemiaulus")] <- tax_table[which(tax_table$scientificName == "Hemiaulus"), ]["Genus"]
df$species[which(df$object_annotation_category == "Hemiaulus")] <- tax_table[which(tax_table$scientificName == "Hemiaulus"), ]["Species"]
df$lifeForm[which(df$object_annotation_category == "Hemiaulus")] <- paste0(tax_table[which(tax_table$scientificName == "Hemiaulus"), ]["Singular_or_colony_forming_reference"], "_" ,tax_table[which(tax_table$scientificName == "Hemiaulus"), ]["In_symbiosis_or_freeliving_reference"])
df$associatedTaxa[which(df$object_annotation_category == "Hemiaulus")] <- "Richelia"

# [6] "multiple<Hemiaulus inter. Richelia"
# [7] "Rhizosolenia inter. Richelia"
# [8] "Trichodesmium"
# [9] "multiple<Trichodesmium"
# [10] "multiple<Rhizosolenia inter. Richelia"
df$scientificName[which(df$object_annotation_category == "multiple<Hemiaulus inter. Richelia")] <- "Hemiaulus"
df$taxonRank[which(df$object_annotation_category == "Hemiaulus")] <- tax_table[which(tax_table$scientificName_used_in_dataset == "Hemiaulus"), ]["taxonRank_of_scientificName"]
df$phylum[which(df$object_annotation_category == "Hemiaulus")] <- tax_table[which(tax_table$scientificName == "Hemiaulus"), ]["Phylum"]
df$class[which(df$object_annotation_category == "Hemiaulus")] <- tax_table[which(tax_table$scientificName == "Hemiaulus"), ]["Class"]
df$family[which(df$object_annotation_category == "Hemiaulus")] <- tax_table[which(tax_table$scientificName == "Hemiaulus"), ]["Family"]
df$genus[which(df$object_annotation_category == "Hemiaulus")] <- tax_table[which(tax_table$scientificName == "Hemiaulus"), ]["Genus"]
df$species[which(df$object_annotation_category == "Hemiaulus")] <- tax_table[which(tax_table$scientificName == "Hemiaulus"), ]["Species"]
df$lifeForm[which(df$object_annotation_category == "Hemiaulus")] <- paste0(tax_table[which(tax_table$scientificName == "Hemiaulus"), ]["Singular_or_colony_forming_reference"], "_" ,tax_table[which(tax_table$scientificName == "Hemiaulus"), ]["In_symbiosis_or_freeliving_reference"])
df$associatedTaxa[which(df$object_annotation_category == "Hemiaulus")] <- "Richelia"

# [7] "Rhizosolenia inter. Richelia"
# [8] "Trichodesmium"
# [9] "multiple<Trichodesmium"
# [10] "multiple<Rhizosolenia inter. Richelia"
df$scientificName[which(df$object_annotation_category == "Rhizosolenia inter. Richelia")] <- "Rhizosolenia"
df$taxonRank[which(df$object_annotation_category == "Rhizosolenia")] <- tax_table[which(tax_table$scientificName_used_in_dataset == "Rhizosolenia"), ]["taxonRank_of_scientificName"]
df$phylum[which(df$object_annotation_category == "Rhizosolenia")] <- tax_table[which(tax_table$scientificName == "Rhizosolenia"), ]["Phylum"]
df$class[which(df$object_annotation_category == "Rhizosolenia")] <- tax_table[which(tax_table$scientificName == "Rhizosolenia"), ]["Class"]
df$family[which(df$object_annotation_category == "Rhizosolenia")] <- tax_table[which(tax_table$scientificName == "Rhizosolenia"), ]["Family"]
df$genus[which(df$object_annotation_category == "Rhizosolenia")] <- tax_table[which(tax_table$scientificName == "Rhizosolenia"), ]["Genus"]
df$species[which(df$object_annotation_category == "Rhizosolenia")] <- tax_table[which(tax_table$scientificName == "Rhizosolenia"), ]["Species"]
df$lifeForm[which(df$object_annotation_category == "Rhizosolenia")] <- paste0(tax_table[which(tax_table$scientificName == "Rhizosolenia"), ]["Singular_or_colony_forming_reference"], "_" ,tax_table[which(tax_table$scientificName == "Rhizosolenia"), ]["In_symbiosis_or_freeliving_reference"])
df$associatedTaxa[which(df$object_annotation_category == "Rhizosolenia")] <- "Richelia"

# [10] "multiple<Rhizosolenia inter. Richelia"
df$scientificName[which(df$object_annotation_category == "multiple<Rhizosolenia inter. Richelia")] <- "Rhizosolenia"
df$taxonRank[which(df$object_annotation_category == "Rhizosolenia")] <- tax_table[which(tax_table$scientificName_used_in_dataset == "Rhizosolenia"), ]["taxonRank_of_scientificName"]
df$phylum[which(df$object_annotation_category == "Rhizosolenia")] <- tax_table[which(tax_table$scientificName == "Rhizosolenia"), ]["Phylum"]
df$class[which(df$object_annotation_category == "Rhizosolenia")] <- tax_table[which(tax_table$scientificName == "Rhizosolenia"), ]["Class"]
df$family[which(df$object_annotation_category == "Rhizosolenia")] <- tax_table[which(tax_table$scientificName == "Rhizosolenia"), ]["Family"]
df$genus[which(df$object_annotation_category == "Rhizosolenia")] <- tax_table[which(tax_table$scientificName == "Rhizosolenia"), ]["Genus"]
df$species[which(df$object_annotation_category == "Rhizosolenia")] <- tax_table[which(tax_table$scientificName == "Rhizosolenia"), ]["Species"]
df$lifeForm[which(df$object_annotation_category == "Rhizosolenia")] <- paste0(tax_table[which(tax_table$scientificName == "Rhizosolenia"), ]["Singular_or_colony_forming_reference"], "_" ,tax_table[which(tax_table$scientificName == "Rhizosolenia"), ]["In_symbiosis_or_freeliving_reference"])
df$associatedTaxa[which(df$object_annotation_category == "Rhizosolenia")] <- "Richelia"

# [8] "Trichodesmium"
# [9] "multiple<Trichodesmium"
df$scientificName[which(df$object_annotation_category == "Trichodesmium")] <- "Trichodesmium"
df$taxonRank[which(df$object_annotation_category == "Trichodesmium")] <- tax_table[which(tax_table$scientificName_used_in_dataset == "Trichodesmium"), ]["taxonRank_of_scientificName"]
df$phylum[which(df$object_annotation_category == "Trichodesmium")] <- tax_table[which(tax_table$scientificName == "Trichodesmium"), ]["Phylum"]
df$class[which(df$object_annotation_category == "Trichodesmium")] <- tax_table[which(tax_table$scientificName == "Trichodesmium"), ]["Class"]
df$family[which(df$object_annotation_category == "Trichodesmium")] <- tax_table[which(tax_table$scientificName == "Trichodesmium"), ]["Family"]
df$genus[which(df$object_annotation_category == "Trichodesmium")] <- tax_table[which(tax_table$scientificName == "Trichodesmium"), ]["Genus"]
df$species[which(df$object_annotation_category == "Trichodesmium")] <- tax_table[which(tax_table$scientificName == "Trichodesmium"), ]["Species"]
df$lifeForm[which(df$object_annotation_category == "Trichodesmium")] <- "free_filaments"
df$associatedTaxa[which(df$object_annotation_category == "Trichodesmium")] <- NA

# [9] "multiple<Trichodesmium"
df$scientificName[which(df$object_annotation_category == "multiple<Trichodesmium")] <- "Trichodesmium"
df$taxonRank[which(df$object_annotation_category == "Trichodesmium")] <- tax_table[which(tax_table$scientificName_used_in_dataset == "Trichodesmium"), ]["taxonRank_of_scientificName"]
df$phylum[which(df$scientificName == "multiple<Trichodesmium")] <- tax_table[which(tax_table$scientificName == "Trichodesmium"), ]["Phylum"]
df$class[which(df$scientificName == "Trichodesmium")] <- tax_table[which(tax_table$scientificName == "Trichodesmium"), ]["Class"]
df$family[which(df$scientificName == "Trichodesmium")] <- tax_table[which(tax_table$scientificName == "Trichodesmium"), ]["Family"]
df$genus[which(df$scientificName == "Trichodesmium")] <- tax_table[which(tax_table$scientificName == "Trichodesmium"), ]["Genus"]
df$species[which(df$scientificName == "Trichodesmium")] <- tax_table[which(tax_table$scientificName == "Trichodesmium"), ]["Species"]
df$lifeForm[which(df$scientificName == "Trichodesmium")] <- "free_filaments"
df$associatedTaxa[which(df$scientificName == "Trichodesmium")] <- NA

# categorical or quantitative field - measurements
# df$occurrence_of_cells <- "PRESENT" # not sure here if needed, or what damiano was aiming for
df$occurrenceStatus <- "PRESENT"

#
df$trichomeStatus[which(df$scientificName == "Trichodesmium")] <- "PRESENT"
df$trichomeQuantity[which(df$scientificName == "Trichodesmium")] <- "1-40"
df$trichomeQuantityType[which(df$scientificName == "Trichodesmium")] <- "L"

#
# df$biomass_conversion_factor <- NA
# df$biomass_conversion_factor_QuantityType <- NA
# df$biomassQuantity <- NA
# df$biomassQuantityType <- NA ## Juan Karlusich calculated biovolume but data is only shown in boxplot graph,

# Method
df$basisOfRecord <- "MachineObservation"
df$measurementMethod <- "Environmental_High_Content_Fluorescence_Microscopy"
df$references <- "Juan_Karlusich_nature_2021"

df$occurrenceID <- paste0(df$scientificName, "_", df$family, "_", df$genus, "_", df$decimalLongitude, "_", df$decimalLatitude, "_", df$year, "_", df$month, "_", df$day, "_", df$day)
df$sampleID <- df$object_htm_sample_barcode
df$stationID <- df$sample_station
df$cruiseID <- "Tara_Ocean"


## Quality control
unique(df$scientificName)
names(df)
summary(df[, c(16:20)])

# there are some lists as columns that cause problems when saving, so unlist those columns, find the columns with: sapply(df_final, class)
# Calothrix, here R sometimes inserts a list that causes problems when saving, so check and unlist if neccassary
class(df$taxonRank) # list, unlist
df$taxonRank <- unlist(df$taxonRank)
class(df$phylum)
df$phylum <- unlist(df$phylum)
class(df$class)
df$class <- unlist(df$class)
class(df$family)
df$family <- unlist(df$family)
class(df$genus)
df$genus <- unlist(df$genus)
class(df$species)
df$species <- unlist(df$species)

# Select columns of intrest
df_final <- df[, unlist(strsplit(names(df_master), " "))]
names(df_final) # 54 columns

# save as csv
filepath <- getwd()
fileName <- "/juan_karlusich_eHCFM_formatted.csv"
write.csv(df_final, paste0(filepath, fileName), row.names = FALSE)


################
## UVP 5 Data ##
################
# load data
dat <- read.table("Trichodesmium_UVP5.tsv", sep = "\t", header = TRUE)
names(dat) # we use columns 1:15
# merge with master_df
library("plyr")
df <- rbind.fill(dat, df_master)
names(df)

####################
## Format columns ##
####################
# Location
df$decimalLongitude <- df$object_lon
df$decimalLatitude <- df$object_lat

# date
df$object_date <- as.character(df$object_date)
df$object_date <- as.Date(df$object_date, format = "%Y%m%d")
df$year <- as.numeric(format(df$object_date, format = "%Y"))
df$month <- as.numeric(format(df$object_date, format = "%m"))
df$day <- as.numeric(format(df$object_date, format = "%d"))


# Taxonomy
unique(df$object_annotation_category)
# [1] "puff" "tuft"

# here R inserts a list type into the column that causes problems when we want to save the data later, unlist each column and check before
df$scientificName <- "Trichodesmium"
class(df$scientificName)
df$taxonRank[which(df$scientificName == "Trichodesmium")] <- tax_table[which(tax_table$scientificName_used_in_dataset == "Trichodesmium"), ]["taxonRank_of_scientificName"]
class(df$taxonRank)
df$taxonRank <- unlist(df$taxonRank)
df$phylum[which(df$scientificName == "Trichodesmium")] <- tax_table[which(tax_table$scientificName == "Trichodesmium"), ]["Phylum"]
class(df$phylum)
df$phylum <- unlist(df$phylum)
df$class[which(df$scientificName == "Trichodesmium")] <- tax_table[which(tax_table$scientificName == "Trichodesmium"), ]["Class"]
class(df$class)
df$class <- unlist(df$class)
df$family[which(df$scientificName == "Trichodesmium")] <- tax_table[which(tax_table$scientificName == "Trichodesmium"), ]["Family"]
class(df$family)
df$family <- unlist(df$family)
df$genus[which(df$scientificName == "Trichodesmium")] <- tax_table[which(tax_table$scientificName == "Trichodesmium"), ]["Genus"]
class(df$genus)
df$genus <- unlist(df$genus)
df$species[which(df$scientificName == "Trichodesmium")] <- tax_table[which(tax_table$scientificName == "Trichodesmium"), ]["Species"]
class(df$species)
df$species <- unlist(df$species)
# keep information on the different colony types
df$lifeForm[which(df$object_annotation_category == "puff")] <- "colony_puffColonies"
df$lifeForm[which(df$object_annotation_category == "tuft")] <- "colony_tuftColonies"
df$associatedTaxa[which(df$scientificName == "Trichodesmium")] <- NA

df$occurrenceStatus <- "PRESENT"
df$cellStatus <- "PRESENT"
df$cellQuantityType <- "colonies"

df$basisOfRecord <- "MachineObservation"
df$measurementMethod <- "Underwater_Vision_Profiler_5"
df$references <- "Juan_Karlusich_nature_2021"

df$occurrenceID <- paste0(df$scientificName, "_", df$family, "_", df$genus, "_", df$decimalLongitude, "_", df$decimalLatitude, "_", df$year, "_", df$month, "_", df$day, "_", df$day)
df$sampleID <- df$sample_id
df$stationID <- df$sample_stationid
df$cruiseID <- df$sample_cruise

# select columns of intrest
df_final <- df[, unlist(names(df_master))]
names(df_final)
dim(df_final) # 210 54

## Quality control
unique(df$scientificName)
names(df)
summary(df[, c(150:154)])

# save as csv
filepath <- getwd()
fileName <- "/juan_karlusich_uvp5_formatted.csv"
write.csv(df_final, paste0(filepath, fileName), row.names = FALSE)






## plot
library("raster")

nc <- raster("C:/Users/domin/Downloads/wod_apb_2020.nc")


# importing and overlaying species data
pinus_edulis <- read.table("tabular/species/p_edulis.txt", header = TRUE, sep = ",")
# check class of the object
class(pinus_edulis)
# convert into spaitial object by assigning coordinates
coordinates(pinus_edulis) <- c("lon", "lat")
coordinates(df_final) <- c("decimalLongitude", "decimalLatitude")
# access the coordinates
coordinates(pinus_edulis)
# assign the projection type
projection(pinus_edulis) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# get species observation directly from the web using dismo package
library("jsonlite")
library("dismo")

pinus_edulis <- gbif("pinus", "edulis", download = TRUE, geo = TRUE,
  sp = TRUE, removeZeros = TRUE)
# rename variable to download date
names(pinus_edulis)[1] <- "dwnld.date"
# save as csv file
write.table(data.frame(pinus_edulis@coords, pinus_edulis@data),
  "tabular/species/p_edulis.txt", sep = ",", row.names = FALSE)

# extract stacked climate data and species observation with bilinear interpolation method
pts.clim <- extract(world.stk, pinus_edulis, method = "bilinear")
pin_edu.clim <- data.frame(cbind(coordinates(pinus_edulis), pts.clim, pinus_edulis@data))
coordinates(pin_edu.clim) <- c("lon", "lat")
# plot
map.ext <- extent(-120, -100, 30, 44)
map.ext <- extent(0, 360, -90, 90)
plot(hillsh_na, col = grey(0:100/100), legend = FALSE, axes = FALSE,
  ext = map.ext)
plot(elev_na, col = cols, add = TRUE, ext = map.ext)
plot(pinus_edulis, pch = 16, cex = .5, add = TRUE)
plot(df_final, pch = 16, cex = .5, add = TRUE)
