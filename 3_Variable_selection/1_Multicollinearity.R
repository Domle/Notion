## 8. CORRELATION ANALYSIS AMONG VARIABLES IN THE OCEAN.


# Input files: Environmental parameters saved as a .grd file.

# Output files: Table containing the pearson and spearman rank correlation coefficients.

# Description:  We will check the correlation between environmental parameters using pearson and spearman rank correlation coefficients and save the 
#               results in a table.

# Author:   Dominic Eriksson
#           Environmental Physics Group, UP
#           ETH Zurich
#           Switzerland

# deriksson@ethz.ch, May 3rd 2023 -------------------------------------------------------------



### =========================================================================
### Preparation
### =========================================================================

# Libraries
lib_vec <- c("psych", "raster", "doParallel", "doBy", "plyr", "doParallel", "foreach", "tidyverse", "usdm")
new_packages <- lib_vec[!(lib_vec %in% installed.packages()[, "Package"])] # check which package is not installed
if(length(new_packages)) {install.packages(new.packages)} # install missing packages
sapply(lib_vec, library, character.only = TRUE) # load packages

# Directories
wd_dat <- "/home/deriksson/Projects/Notion_DE/Code/3_Match_up/2_Output/"
output.dir <- "/net/kryo/work/deriksson/Projects/Notion_DE/Code/5_Variable_selection/5_Output/"

# Functions
source("/net/kryo/work/deriksson/Projects/Notion_DE/Code/Functions/Righetti/correl.group.fun")
source("/net/kryo/work/deriksson/Projects/Notion_DE/Code/Functions/Righetti/pairs.panel.col.fun") # contains a modified (colored) pairs.panel function (pairs.panel.col)

### =========================================================================
###  Load and format the data
### =========================================================================

## Preparatory
# Load data: *** select gridded or thinned data points from the matchup step --> also tried the ungridded dataset, leads to same results
list.files(wd_dat)
tot <- read.csv(paste0(wd_dat, "Obs_gridded_diazo_pres_abs.csv"))

# Definition of predictors considered
vc.vars <- c(
    # Physical 
        "MLD1", "MLD2", "PAR", "Wind.CCMP", "T", "deltaT", "sdT", "dT_dt", 
    # Chemical
        "Sal", "N", "P", "Si","pCO2",
    # Biological 
        "Chl", "NPP", 
        "MLPAR1", "MLPAR2", "logMLD1", "logMLD2", "deltaWind.CCMP","monsdWind.CCMP", "cvWind.CCMP","deltaMLD1", "sdMLD1", "cvMLD1", "deflMLD1",  "dMLD1_dt", # Derived physical (inserted for variability min stack; cv.Wind.CCMP), sd.MLD1, cv.MLD1
        "logN","logP","logSi","deltaN", "deltaP", "deltaSi","Nstar","Pstar","Sistar", "deltaNstar", "deltaSistar", "dN_dt", "dP_dt", "dSi_dt", # Derived chemical
        "deltapCO2", "logChl","deltaChl", "dNPP_dt")

# => Comment: This variable stack can be expanded by further predictors suited to the organisms analyzed. E.g. dissolved oxygen (02) as a driver of heterotrophic bacteria.
# => Comment: Several variables, such as sdXY, cvXY, or deltaXY are excluded later in the modeling, since they have not monthly (but annual) data resolution

# Definition of reduced predictors considered for pairs panels colored plot
vc.vars.reduced <- c(
    "T","Sal","N","P","Si","MLD1","PAR","Chl","Wind.CCMP","pCO2",
    "MLPAR1", "Nstar", "Sistar",
    "dT_dt", "dN_dt", "dMLD1_dt",
    "logMLD1",  "logChl", "logN", "logP", "logSi")
# => Comment: This variable stack can be expanded by further predictors suited to the organisms analyzed. E.g. dissolved oxygen (02) as a driver of heterotrophic bacteria.
# => Comment: These and additional variables are used for monthly plankton distribution modeling


 # Get environmental predictor data status Sci Adv Article 2019 (Righetti et al, 10.1126/sciadv.aau6253)
prj.stack <- brick("/net/kryo/work/deriksson/Projects/Notion_DE/Data/Predictor_variables/Environmental_data/VarSet07.2016.grd")
nelem <- nlayers(prj.stack)/12
nm <- gsub("\\.1", "", names(prj.stack[[1:nelem]]))
env.stack<-list() ; from<-seq(1,nlayers(prj.stack),by=nelem) ; for(q in 1:12){pre.env.stack <- stack(prj.stack[[ from[q]:(from[q]+nelem-1)]]); names(pre.env.stack) <- nm; env.stack[[q]] <- pre.env.stack}
names(env.stack) <- c("JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC")



### =========================================================================
###  Calculate correlation coefficients and save table
### =========================================================================

# Create data frame denoting all cells containing variable information for global ocean, takes 10 sec
lisi.df <- list ()
for (k in 1:12 ) { 
    lisi.df [[k]] <- as.data.frame(rasterToPoints(env.stack [[k]] ))}
    df.tot <- do.call("rbind", lisi.df)
# Inspect
head(df.tot)

## Pearson
# full ocean
df <- correl.group (na.omit(df.tot[ , vc.vars ]), PlotFileName=paste0(output.dir, "Plot_to_delete"), cor.Method = "pearson")$AbsoluteValueCorrelations
fln <- paste0(output.dir, "Corr_table_global_pearson.csv")
write.csv(df, file = fln)

## Spearman
# full ocean
df <- correl.group (na.omit(df.tot[ , vc.vars ]), PlotFileName=paste0(output.dir, "Plot_to_delete"), cor.Method = "spearman")$AbsoluteValueCorrelations
fln <- paste0(output.dir, "Corr_table_global_spearman.csv") # * CURRENT CHOICE * for subsequent variable exclusion due to high correlations within individual models
write.csv(df, file = fln)

###==============================================================
### END
###==============================================================