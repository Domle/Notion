


## Preparatory steps
rm(list = ls())
library(raster); library(sfsmisc); library(grid); library(spam); library(fields); library(base); library(doBy); library(ncdf4); library(grid)
setwd("/net/kryo/work/deriksson/notion/")

# Get environmental data
prj.stack <- brick(paste0(".", "/2_Environmental_data/VarSet07.2016.grd"))
dim(prj.stack) # 180 360 576
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

## Load data on total presences of diazotroph phytoplankton
# Get gridded data
list.datasets <- list()
list.datasets[[1]] <- read.csv("./3_Match_up/Output_1/diazos_pres_abs_gridded.csv")
list.datasets[[2]] <- read.csv("./3_Match_up/Output_1/diazos_nifHcounts_abs_gridded.csv")
list.datasets[[3]] <- read.csv("./3_Match_up/Output_1/diazos_pres_abs_raw.csv")
list.datasets[[4]] <- read.csv("./3_Match_up/Output_1/diazos_nifHcounts_abs_raw.csv")
lapply(list.datasets, nrow)

## Add environmental variables based on criterion (gridded x, gridded y, month)
# Preparatory
vec.dat <- c("pres_abs", "nifHcounts_abs", "pres_abs", "nifHcounts_abs") # names of dataset versions

### LOOP across data-type versions (presence data versus count data using consistent methodology)
for(s in 1:length(vec.dat)){
  ## Add environmental information
  dat <- list.datasets[[s]]
  llCRS <- CRS("+proj=longlat +ellps=WGS84")
  crd.index <- which(colnames(dat)%in%c(c("x", "y")))
    ## Convert data points to SpatialDataFrame
    spdf <- SpatialPointsDataFrame(dat[, crd.index], dat[, -crd.index], proj4string = llCRS) # coords not shown, only as data.frame: spdf.text <- as.data.frame(spdf)
    ## Loop across months: mask the observations and prepare list to extract points, column-bind with mask, keep rows (points) with mask equal 1 (ocean condition)
    mups <- list()
      for(k in 1:12){
        # 1. Define monthly subset, for monthly match-up with mask
        if(sum(spdf$month == k)> 0){ss <-subset(spdf, spdf$month == k)} else {next}
        # 2. Extract data points overlapping with "mask"
        extr <- extract(x = env.stack[[k]], y = ss, method = "simple") # returns the values of a raster object (msk) for cells in which points fall
        extr <- as.data.frame(extr) # the mask values are converted to data frame
        ss.new <- cbind(as.data.frame(ss), extr)[which(extr$Mask!=0), ] # CORE OPERATION: SELECT ROWS/DATA for which mask == 1
        # 3. Store monthly specific ss.new in list (for later row-bind)
        mups[[length(mups) + 1]] <- ss.new
        print(paste(k))
      }
    # Merge match-ups
    tot <- do.call("rbind", mups)
    # Usefule gridded 1° by 1° monthly data records for modelling
    dim(tot); class(tot) # 5042   59; data frame

    # Note that any information on the measurement method is lost at this point, and coastal data are excluded
    # plot(tot$x, tot$y, xlim = c(-180, 180), ylim = c(-90, 90))

    # overview of taxa
    dat.spl <- split(tot, f = tot$scientificName)
    df.taxa.obs <- data.frame("taxon" = names(unlist(lapply(dat.spl, nrow))), "obs" = as.numeric(unlist(lapply(dat.spl, nrow))))
    df.taxa.obs <- df.taxa.obs[order(df.taxa.obs$obs, decreasing = T), ]
    rownames(df.taxa.obs) <- 1:nrow(df.taxa.obs)
    df.taxa.obs <- df.taxa.obs[which(df.taxa.obs$obs > 0), ]
    df.taxa.obs <- df.taxa.obs[order(df.taxa.obs$taxon), ]
    rownames(df.taxa.obs) <- 1:nrow(df.taxa.obs)
    df.taxa.obs # 49 taxa

# Here Damiano made a choice on a minimum on presences that are required to model the species, minimum not know at this stage
    # Dataset I: -------------------------------------------------------------------------------------------------------------
    # The potentially modelable diazotroph taxa (in bracket, the number of gridded monthly presences or non-detections) are:
    # 1. Azotobacter (18)
    # 2. CIII (18)
    # 3. Calothrix (109)
    # 4. Calothrix confervicola (21)
    # 5. Gamma (112)
    # 6. Gamma.A (86)
    # 7. Gamma.P (26)
    # 8. Gamma_1 (16)
    # 9. Heterocyst.Richelia.Calothrix (95)
    # 10. Richelia (508)
    # 11. Richelia intracellularis (43)
    # 12. Trichodesmium (1882)
    # 13. Trichodesmium erythraeum (216)
    # 14. Trichodesmium thiebautii (112)
    # 15. UCYN (315)
    # 16. UCYN.A (286)
    # 17. UCYN.A1 (99)
    # 18. UCYN.A2 (51)
    # 19. UCYN.B (205)
    # 20. UCYN.C (107)
    #
    # The following taxa are independent, chosing the most detailed resolution possible, except there were only one species per genus available (case: Richelia)
    # 1. Azotobacter (18)
    # 2. CIII (18)
    # 3. Calothrix confervicola (21)
    # 4. Gamma.A (112)
    # 5. Gamma.P (26)
    # 6. Richelia (508)
    # 7. Trichodesmium erythraeum (216)
    # 8. Trichodesmium thiebautii (112)
    # 9. UCYN.A1 (99)
    # 10. UCYN.A2 (51)
    # 11. UCYN.B (205)
    # 12. UCYN.C (107)
    #
    # Potentially modelable host taxa (in bracket, the number of gridded monthly presences or non-detections) are :
    # 1. Chaetoceros (143)
    # 2. Chaetoceros affinis (124)
    # 3. Chaetoceros debilis (28)
    # 4. Chaetoceros didymus (42)
    # 5. Rhizosolenia (112)
    # 6. Rhizosolenia fallax (43)
    # 7. Rhizosolenia setigera (19)
    # 8. Rhizosolenia shrubsolei (51)

    # Dataset II: -------------------------------------------------------------------------------------------------------------
    # The potentially modelable diazotroph taxa (in bracket, the number of gridded monthly presences or non-detections) are:
    # 1. CIII (18)
    # 2. Gamma (112)
    # 3. Gamma.A (86)
    # 4. Gamma.P (26)
    # 5. Gamma_1 (16)
    # 6. Richelia (164)
    # 7. Trichodesmium (251)
    # 8. UCYN (305)
    # 9. UCYN.A (275)
    # 10. UCYN.A1 (70)
    # 11. UCYN.A2 (42)
    # 12. UCYN.B (203)
    # 13. UCYN.C (104)
    #
    # The following taxa are independent, chosing the most detailed resolution possible, except there were only one species per genus available (case: Richelia)
    # 1. CIII (18)
    # 2. Gamma.A (86)
    # 3. Gamma.P (26)
    # 4. Richelia (164)
    # 5. Trichodesmium (251)
    # 6. UCYN.A1 (70)
    # 7. UCYN.A2 (42)
    # 8. UCYN.B (203)
    # 9. UCYN.C (104)
    #
    # Potentially modelable host taxa (in bracket, the number of gridded monthly presences or non-detections) are : none

  ## Order data
  # Define taxa as character
  tot$scientificName <- as.character(tot$scientificName)

  # Re-order the data along taxa
  lisi.taxa <- list()
    # Loop across taxa: we keep taxa == NA such as picoeukaryotes for the moment
    for(i in 1:nrow(df.taxa.obs)){
        # 1. get data of taxon
        lisi.taxa[[length(lisi.taxa) + 1]] <- tot[which(tot$scientificName == df.taxa.obs$taxon[i]), ]
        # Re-merge groups
        # lapply(lisi.taxa, head)
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

        # Save list of obs with env.data (without year as factor)
        if(s%in%c(1, 2)){fln <- paste("./3_Match_up/Output_2/Obs_gridded_diazo_",vec.dat[s],".csv", sep="")}
        if(s%in%c(3,4)){fln <- paste("./3_Match_up/Output_2/Obs_raw_diazo_",vec.dat[s],".csv", sep="")}

        write.csv(tot.new, file = fln, row.names = F)
        print(paste(s)) # progress
    } # close loop across data-type versions
  }






### ================================= PS: outdated - but not out of interest ====================================



### Notes: Plot of environmental variables correlation analysis new env stack versus maredat env stack
# Explore correlation of Env.Stack implemented from WOA etc. versus Env.Stack from maredat --> CAN'T FIND THAT MAREDAT DATA
dat <- read.csv(".//3_Match_Up/Output/Total_Obs.csv", header = T, sep = ",") # total observation data (sp, genus, families, infrakingdoms), matched with env data

# Available variables: names(dat)
# "Bathymetry", "MLD.deBoyerMontegut.", "PAR.SeaWiFS.", "Temperature.1", "Salinity.1", "Nitrate.1", "Phosphate", "Silicate" (origin: new dta on environmental variables)
# "Nitrate.2", "Salinity.2", "Temperature.2", "log_MLD", "Pstar", "Sistar", "MLPAR" (origin: environmental data used with maredat by M. Vogt for species distribution modeling)

# Annually Integrated Climatology Plot
output.dir.var <- getwd()
png(filename = paste0(output.dir.var, "/", "Vars_Corr_tot.png"), width = 5300, height = 4000, units = "px", pointsize = 80)
par(mfrow = c(2, 3))
plot(dat$Temperature, dat$Temperature.2, cex = 0.4)
plot(dat$Salinity.1, dat$Salinity.2, cex = 0.4)
plot(dat$Nitrate.1, dat$Nitrate.2, cex = 0.4)
plot(dat$log_MLD, dat$MLD, cex = 0.4)
plot(dat$PAR.SeaWiFS, dat$MLPAR, cex = 0.4)
plot(dat$Pstar, dat$Sistar, cex = 0.4)
dev.off

# Annually Integrated Climatology Plot
png (filename = paste0(output.dir.var,"/","Vars_Corr_tot.png"), width = 5300, height = 4000, units = "px", pointsize = 80)
par (mfrow=c(2,3))
plot (dat$Temperature.1, dat$Temperature.2, cex=0.4)
plot (dat$Salinity.1, dat$Salinity.2, cex = 0.4)
plot (dat$Nitrate.1, dat$Nitrate.2, cex = 0.4)
plot (dat$log_MLD, dat$MLD, cex = 0.4)
plot (dat$PAR.SeaWiFS, dat$MLPAR, cex = 0.4)
plot (dat$Pstar, dat$Sistar, cex = 0.4)
dev.off()

# Monthly Integrated Climatology Plots
x <- c("1.JAN","2.FEB","3.MAR","4.APR","5.MAY","6.JUN","7.JUL","8.AUG","9.SEP","10.OCT","11.NOV","12.DEC")
for (k in 1:12) {
png (filename = paste0(output.dir.var,"/","Vars_Corr","_",x[k],".png"), width = 5300, height = 4000, units = "px", pointsize = 80)
par (mfrow=c(2,3))
dat.m <- subset(dat,dat$Month==k)
plot (dat.m$Temperature.1, dat.m$Temperature.2, cex=0.4)
plot (dat.m$Salinity.1, dat.m$Salinity.2, cex = 0.4)
plot (dat.m$Nitrate.1, dat.m$Nitrate.2, cex = 0.4)
plot (dat.m$log_MLD, dat.m$MLD, cex = 0.4)
plot (dat.m$PAR.SeaWiFS, dat.m$MLPAR, cex = 0.4)
plot (dat.m$Pstar, dat.m$Sistar, cex = 0.4)
dev.off()
}
