### ORGANIZE TAXA INTO PARENT TAXA (WHERE NEEDED) 
### THUS, STRAINS OR SPECIES ARE CONTAINED WITHIN THEIR PARENTTAXON, AND THE TOTAL TAXA ARE LISTED UNDER "SCIENTIFIC NAME"
### THE GOAL IS HENCE A CROSS-HIERARCHIC DATABASE, i.e., taxa of different resolution build this dataset, while the taxonomic rank of each taxon is indicated by the field "taxonRank"

# $Date: 2020-06-09
# Author: Damiano Righetti, damiano.righetti@env.ethz.ch
# Environmental Physics Group, ETH Zurich, CH
# Note: An aggregation step has not been thorougly implemented so far.
# Aggregation step: Merge species and/or strains into larger groups useful for modelling, and/or reflecting taxonomic hierarchy

### ===================================================================================
### Preparatory steps
### ===================================================================================
rm(list = ls())
library(doBy);library(plyr);library(doParallel);library(foreach);library(tidyverse)
setwd("C:/Users/Dominic Eriksson/Desktop/ETH Internship/Data_notion/1_Merge_data") # user one
setwd("~/Desktop/Data_notion/1_Merge_data") # user two

# Get data
dat <- read.csv("./Output_1/Diazo.csv") # Merged data
tax.table <- read.csv("../1_Download_data/commons/Taxon_table.csv") # Taxonomic table 

# Explore species and their number of points
dat.spl <- split(dat, f = dat$scientificName) 
df.taxa.obs <- data.frame("taxon" = names(unlist(lapply(dat.spl, nrow))), "obs" = as.numeric(unlist(lapply(dat.spl, nrow))) ) 
df.taxa.obs <- df.taxa.obs[order(df.taxa.obs$obs, decreasing = T),]
rownames(df.taxa.obs) <- 1:nrow(df.taxa.obs)
df.taxa.obs <- df.taxa.obs[which(df.taxa.obs$obs > 0),]
df.taxa.obs # 59 taxa

### =========================================================================
### Corrections of taxa
### =========================================================================

## Assign Cyanothece spp. and Cyanothece to UCYN.C 
dat[which(dat$scientificName == "Cyanothece"), "scientificName"] <- "UCYN.C"
dat[which(dat$scientificName == "Cyanothece majus"), "scientificName"] <- "UCYN.C" # Questionable if correct, needs to be checked
dat[which(dat$scientificName == "UCYN.C"), "taxonRank"]  <- tax.table[which(tax.table$scientificName_used_in_dataset == "UCYN.C"), "taxonRank_of_scientificName"] 

### =========================================================================
### Aggregate taxa
### =========================================================================

## 1. ADD Gamma.A, Gamma_1, Gamma_4, and Gamma.P to Gamma
subset.gamma <- dat[which(dat$scientificName%in%c("Gamma.A", "Gamma_1", "Gamma_2", "Gamma_4", "Gamma.P")), ]
subset.gamma[ , "scientificName"] <- "Gamma"
subset.gamma$taxonRank <- tax.table[which(tax.table$scientificName_used_in_dataset == "Gamma"), "taxonRank_of_scientificName"]
## *** Note there is a BUG in line 11738, a record sourced from Gradoville -> a taxon called "gamma atelocyanobacterium", which makes no sense. Needs to be adjusted in gradoville raw merger file. 
##
## In addition: ADD Gamma_1, Gamma_4 to Gamma.A
subset.gamma.A <- dat[which(dat$scientificName%in%c("Gamma_1", "Gamma_2", "Gamma_4")), ]
subset.gamma.A[ , "scientificName"] <- "Gamma.A"
subset.gamma.A$taxonRank <- tax.table[which(tax.table$scientificName_used_in_dataset == "Gamma.A"), "taxonRank_of_scientificName"]

## 2. ADD UCYN-variants and strains to UCYN  
subset.ucyn <-  dat[which(dat$scientificName%in%c("UCYN", "UCYN.A", "UCYN.A1", "UCYN.A2", "UCYN.A2.A3", "UCYN.B", "UCYN.C", "Cyanothece", "Cyanothece majus")), ]
subset.ucyn[ , "scientificName"] <- "UCYN"
subset.ucyn$taxonRank <- tax.table[which(tax.table$scientificName_used_in_dataset == "UCYN"), "taxonRank_of_scientificName"]
##
## In addition: ADD UCYN.A1, UCYN.A2, UCYN.A2.A3 to UCYN.A
subset.ucyn.A <-  dat[which(dat$scientificName%in%c( "UCYN.A1", "UCYN.A2", "UCYN.A2.A3")), ]
subset.ucyn.A[ , "scientificName"] <- "UCYN.A"
subset.ucyn.A$taxonRank <- tax.table[which(tax.table$scientificName_used_in_dataset == "UCYN.A"), "taxonRank_of_scientificName"]

## 4. ADD Richelia intracellularis to Richelia
subset.ric <- dat[which(dat$scientificName%in%c("Richelia intracellularis")), ]
subset.ric[ , "scientificName"] <- "Richelia"
subset.ric$taxonRank <- tax.table[which(tax.table$scientificName_used_in_dataset == "Richelia"), "taxonRank_of_scientificName"]

## 5. ADD Calothrix confervicola, Calothrix crustacea and further species to Calothrix 
subset.cal <- dat[which(dat$scientificName%in%c("Calothrix confervicola", "Calothrix crustacea", "Calothrix scopulorum", "Calothrix aeruginea", "Calothrix pilosa", "Calothrix confervicola var. purpurea", "Calothrix fuscoviolacea", "Calothrix nodulosa", "Calothrix parietina")), ]
subset.cal[ , "scientificName"] <- "Calothrix"
subset.cal$taxonRank <- tax.table[which(tax.table$scientificName_used_in_dataset == "Calothrix"), "taxonRank_of_scientificName"]

## 6. ADD Trichodesmium species to Trichodesmium
subset.tri <- dat[which(dat$scientificName%in%c("Trichodesmium hildebrandtii", "Trichodesmium contortum", "Trichodesmium thiebautii", "Trichodesmium erythraeum")), ]
subset.tri[ , "scientificName"] <- "Trichodesmium"
subset.tri$taxonRank <- tax.table[which(tax.table$scientificName_used_in_dataset == "Trichodesmium"), "taxonRank_of_scientificName"]

## 7. ADD Azotobacter armeniacus to Azotobacter
subset.azo <- dat[which(dat$scientificName%in%c("Azotobacter armeniacus")), ]
subset.azo[ , "scientificName"] <- "Azotobacter"
subset.azo$taxonRank <- tax.table[which(tax.table$scientificName_used_in_dataset == "Azotobacter"), "taxonRank_of_scientificName"]

## ----------- Diatom host species (heterocyst forming nitrogen fixers, such as Richelia and Calothrix are associated with certain diatoms, listed here below) -----------
## 8. ADD Chaetoceros species to Chaetoceros
subset.cha <- dat[which(dat$scientificName%in%c("Chaetoceros affinis", "Chaetoceros radicans", "Chaetoceros didymus", "Chaetoceros debilis", "Chaetoceros muellerii", "Chaetoceros rostratus", "Chaetoceros calcitrans", "Chaetoceros brevis", " Chaetoceros elegans")), ]
subset.cha[ , "scientificName"] <- "Chaetoceros"
subset.cha$taxonRank <- tax.table[which(tax.table$scientificName_used_in_dataset == "Chaetoceros"), "taxonRank_of_scientificName"]

## 9. ADD Rhizosolenia species to Rhizosolenia
subset.rhi <- dat[which(dat$scientificName%in%c("Rhizosolenia shrubsolei", "Rhizosolenia fallax", "Rhizosolenia setigera")), ]
subset.rhi[ , "scientificName"] <- "Rhizosolenia"
subset.rhi$taxonRank <- tax.table[which(tax.table$scientificName_used_in_dataset == "Rhizosolenia"), "taxonRank_of_scientificName"]

### =========================================================================
### Merge data
### =========================================================================

# Merge
dat.new <- rbind(dat, subset.gamma, subset.gamma.A, subset.ucyn, subset.ucyn.A, subset.ric, subset.cal, subset.tri, subset.azo, subset.cha, subset.rhi)

# Re-explore species and their number of points
dat.spl <- split(dat.new, f = dat.new$scientificName) 
df.taxa.obs <- data.frame("taxon" = names(unlist(lapply(dat.spl, nrow))), "obs" = as.numeric(unlist(lapply(dat.spl, nrow))) ) 
df.taxa.obs <- df.taxa.obs[order(df.taxa.obs$obs, decreasing = T),]
df.taxa.obs <- df.taxa.obs[which(df.taxa.obs$obs > 0),]
df.taxa.obs <- df.taxa.obs[order(df.taxa.obs$taxon), ]
rownames(df.taxa.obs) <- 1:nrow(df.taxa.obs)
df.taxa.obs # 57 taxa

# Sort alphabetically (taxon-wise)
lisi.taxa <- list()
for(i in 1:nrow(df.taxa.obs)){
lisi.taxa[[i]] <- dat.new[which(dat.new$scientificName == df.taxa.obs[ i , "taxon"]), ]
}
dat.hierarchical <- do.call(rbind, lisi.taxa)
dim(dat.hierarchical) # 30 528  x  55
rownames(dat.hierarchical) <- 1:nrow(dat.hierarchical)

# Save 
fln <- "./Output_2/Diazo_hierarchic_taxa.csv"
write.csv(dat.hierarchical, file = fln, row.names = F)






