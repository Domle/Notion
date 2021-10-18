### PLOT DIAZOTROPH RECORDS IN SPACE AND TIME (following Righetti et al., 2020; Phytobase: A global synthesis of open ocean phytoplankton occurrences, ESSD)

# Date: 2020-06-09
# Author: Damiano Righetti, rdamiano@ethz.ch
# Environmental Physics Group, ETH Zurich, CH

### =========================================================================
### Initialize system
### =========================================================================

## Clear home
rm(list = ls())
setwd("D:/Research/NOTION/") # User one

## Libraries
library(rgdal); library(fields); library(raster); library(maps);
library(colorRamps); library(RColorBrewer); library(sfsmisc); library(matlab);
library(rasterVis); library(sp); library(ggplot2); library(RColorBrewer);
library(plyr); library(reshape2); library(akima); library(rworldmap) # this package has way better world shapefiles
library(spatialEco); library(raster); library(ncdf4); library(maps); library(ncdf4); library(grid); library(spam); library(fields); library(base); library(doBy)

# Load functions from "commons"
source('~/Desktop/Data/COMMONS/common_functions.R', chdir = TRUE)


###  ===============================================================
### FIGURE 1A | PRESENCE RECORDS IN SPACE
###  ===============================================================

# Preparatory: Load raster layer of ocen temperature, climatological monthly data Locarini et al., 2013
prj.stack <- brick(paste0(".","/2_Environmental_data/VarSet07.2016.grd"))
nelem <- nlayers(prj.stack)/12
nm <- gsub("\\.1", "", names(prj.stack[[1:nelem]]))
env.stack<-list() ; from<-seq(1,nlayers(prj.stack),by=nelem) ; for(q in 1:12){pre.env.stack <- stack(prj.stack[[ from[q]:(from[q]+nelem-1)]]); names(pre.env.stack) <- nm; env.stack[[q]] <- pre.env.stack}
names(env.stack) <- c("JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC")
layer.var <- calc( stack(lapply(env.stack, '[[', which(names(env.stack[[1]]) == "T" ) ) ), mean, na.rm = T)
layer.var[env.stack[[1]]$Mask== 0] <- NA
temp.levels <- layer.var
temp.levels[temp.levels>20] <- 30
temp.levels[temp.levels>=10 & temp.levels<20  ] <- 20
temp.levels[temp.levels<10 ] <- 10

### Preparatory: Load data
dat <- read.csv("~/Desktop/Data_notion/1_Merge_data/Output_2/Diazo_hierarchic_taxa.csv")
dat <- dat[which(dat$occurrenceStatus=="PRESENT"), ]
head(dat)

# Stats on amount of data in North Atlantic (discussion point)
dim(dat[which(dat$decimalLongitude> -90   & dat$decimalLongitude < 0  & dat$decimalLatitude > 0  ), ])[1] / dim(dat)[1]
dt <- dat[which(dat$decimalLongitude> -90   & dat$decimalLongitude < 0  & dat$decimalLatitude > 0  ), ]
# plot(dt$decimalLongitude, dt$decimalLatitude, xlim = c(-180,180), ylim = c(-90,90))

# Peruvian coast cluster
dim(dat[which(dat$decimalLongitude> -90  & dat$decimalLongitude < -60 & dat$decimalLatitude > -30 & dat$decimalLatitude < 0  ), ])[1] / dim(dat)[1]
dt <- dat[which(dat$decimalLongitude> -90  & dat$decimalLongitude < -60  & dat$decimalLatitude > -30 & dat$decimalLatitude < 0  ), ]
# plot(dt$decimalLongitude, dt$decimalLatitude, xlim = c(-180,180), ylim = c(-90,90))

# Make it workable
dat$x <- dat$decimalLongitude
dat$y <- dat$decimalLatitude
dat$taxon <- dat$scientificName

# Varia
cea <- 1.2 # axis lables type size
sz <- 0.3 # axes in plot
# c.col <- "grey86" # continents color
c.col <- "#D6EFED" # continents color
c.col <- "grey63" # continents color
s.col <- "white" # s line color (?)
col.pts <- "black" # color of points
col.pts <- "navy" # color of points
col.pts <- "black" # color of points
point <- 1
point <- 16 # point symbol
col.line <- "firebrick1"
col.line <- "grey65" # color of help-lines
#
sz <- 0.3 # axes numbers in plot
nr.size <- 0.3 # axes numbers in lat grad plot
# c.col <- "grey86"
c.col <- "#D6EFED"
s.col <- "white"
col.pts <- "navy"
col.pts <- "black"
pt.sym <- 16
lwd.box <- 15
dist <- 0.5
#
# Lables
cea <- 1
lab.size <- 0.93
tck.dist <- -0.006 # originally -0.007
dist <- 0.22 # distance of numbers (axis lables) from axis
#
# Points
pt.size <- 0.45 # size of points
lwd.size <- 4 # width of lines of points
pt.sym <- 1
pt.col <- "black"
#
# Lines and continents
lwd.box <- 10
c.col <- "black"
#
###  Preparatory: Colors und definitions
col.hot <- "lemonchiffon3"
col.cold <- colorRampPalette(c("#B7C8DE",              colorRampPalette(c("black", "royalblue"))(5)[2])                        ) (20)[2]

### PLOT ###
#
png (filename =  paste0("./1_Visualize_data/Output_3/ESSD_style_map_sources_in_color.png"), width = 6100, height = 4100, units = "px", pointsize = 80) # original width 5600, original height 3600
# 0. Empty plot
plot(dat$decimalLongitude, dat$decimalLatitude,  pch = pt.sym, xlab="",  ylab="",  ylim = c(-77,84), xlim = c (-167, 167), col = pt.col, axes = F,  xaxt = "n", yaxt = "n", cex = pt.size, type ="n")
#
#
## 1. Optional: add some background shading/mask on the ocean regions: tropical, temperate and polar: the last one of the three is standard
# plot (temp.levels, add = T, legend = F, col = c(alpha("#003366", fact) , alpha("#CCCCCC", fact), alpha("#FFCC00", fact))) #2
# plot (temp.levels, add = T, legend = F, col = c( alpha(col.cold, 0.22) , alpha("snow", 1), alpha(col.hot, 0.3) )) #2
plot (temp.levels, add = T, legend = F, col = c(alpha("#003366", 0.093), alpha("grey97", 1), alpha(col.hot, 0.335) )) # 0.1, 1, 0.4 originally
#
#
## 2. Grid lines
axis (1, at = seq(-180, 180 ,30), labels = F, cex.axis = 0.4, tck = 1, lwd = 0.8, lty = 1, mgp = c(1, dist, 0), col = "grey30") # y-axis grid-lines
axis (2, at=seq(-90,90,15), labels = F, cex.axis = 0.4, tck = 1, lwd = 0.8, lty = 1, mgp = c(1, dist, 0), col = "grey30") # x-axis grid-lines
# axis (2, at=seq(-90,90,10), labels = F, cex.axis = 0.4, tck = 1, lwd = 0.8, lty = 1, mgp = c(1, dist, 0), col = "grey30") # x-axis grid-lines, originally I took this
#
#
## 3. Map
map ("world", interior=TRUE, fill = T, lwd = 0.5, boundary=TRUE, col="black", col.fill = "black",add = T)
# shape <- readOGR(dsn = "./shapefiles_world/ne_110m_coastline/ne_110m_coastline.shp")
# shape # it seems to have the same projection..
# plot(shape, add = T, lwd = 3)
#
#
## 4. Visualization of data: microscopic records = circle (pch 1), nifH records = filled triangle (pch 25), OTU records = cube (pch 22)
## Luo
tot <- dat
dt <- tot[which(tot$set=="Luo"),] # unique(dt$measurementMethod)
points(dt[which(dt$measurementMethod%in%c("Standard light microscopy", "Epifluorescence microscopy", "Standard light or epifluorescence microscopy", "Direct biomass analysis")),]$x,
dt[which(dt$measurementMethod%in%c("Standard light microscopy", "Epifluorescence microscopy", "Standard light or epifluorescence microscopy", "Direct biomass analysis")),]$y, col = "orange", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 1)
points(dt[which(dt$measurementMethod%in%c("qPCR_nifH_detection")),]$x, dt[which(dt$measurementMethod%in%c("qPCR_nifH_detection")),]$y, col = "orange", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 25)
#
## Tan
dt <- tot[which(tot$set=="Tan"),] # unique(dt$measurementMethod)
points(dt[which(dt$measurementMethod%in%c("qPCR_nifH_detection")),]$x, dt[which(dt$measurementMethod%in%c("qPCR_nifH_detection")),]$y, col = "darkblue", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 25)
#
## Phy
dt <-  tot[which(tot$set=="Phy"),]  # unique(dt$measurementMethod)
points(dt$x, dt$y, col = "yellow2", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 1) # Assuming microscopy basaed (information on method, NA)
#
## Gra
dt <- tot[which(tot$set=="Gra"),] # unique(dt$measurementMethod)
points(dt[which(dt$measurementMethod%in%c("qPCR_nifH_detection")),]$x, dt[which(dt$measurementMethod%in%c("qPCR_nifH_detection")),]$y, col = "firebrick2", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 25)
#
## Gbi
dt <-  tot[which(tot$set=="Gbi"),] # unique(dt$measurementMethod)
points(dt$x, dt$y, col = "seagreen2", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 1) # Assuming microscopy basaed (information on method, NA)
#
## Obi
dt <-  tot[which(tot$set=="Obi"),] # unique(dt$measurementMethod)
points(dt$x, dt$y, col = "pink", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 1) # Assuming microscopy basaed (information on method, NA)
#
## Tara OTU based
dt <-  tot[which(tot$set=="TarOTU"),] # unique(dt$measurementMethod)
points(dt[which(dt$measurementMethod%in%c("SSU_rRNA_detection")),]$x, dt[which(dt$measurementMethod%in%c("SSU_rRNA_detection")),]$y, col = "slateblue1", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 22)
#
## Martinez-Perez (M96) based
dt <-  tot[which(tot$set=="Mar"),] # unique(dt$measurementMethod)
points(dt$measurementMethod%in%c("Epifluorescence microscopy")),]$x, dt[which(dt$measurementMethod%in%c("Epifluorescence microscopy")),]$y,
col = "firebrick3", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 1)
points(
jitter(x = dt[which(dt$measurementMethod%in%c("qPCR_nifH_detection")),]$x, factor = 1, amount = 1.2),
jitter(x = dt[which(dt$measurementMethod%in%c("qPCR_nifH_detection")),]$y, factor = 1, amount = 0.6),
col = "firebrick3", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 25)
#

legend(-160, 91.7, c("Luo", "Tang", "PhytoBase", "Gradoville", "Gbif", "Obis", "Tara", "M96"), lty=c(1, 1, 1, 1, 1, 1, 1, 1,1), col=c("orange", "darkblue","yellow2", "firebrick2", "seagreen2", "pink", "slateblue1", "firebrick3" ), bty = "n", lwd=32, cex=lab.size, horiz =T, text.width=22)
legend(-160, - 72, c(
paste0("microscopic or NA (",nrow(tot[which(tot$measurementMethod%in%c("Standard Light Microscopy", "Epifluorescence Microscopy", "Standard light or epifluorescence microscopy", "Direct Biomass Analysis", "Assumed microscopic")),]),")"),
paste0("nifH records (",nrow(tot[which(tot$measurementMethod%in%c("qPCR_nifH_detection")),]),")"),
paste0("SSU rRNA records (",nrow(tot[which(tot$measurementMethod%in%c("SSU_rRNA_detection")),]),")        Total: ",nrow(tot))
), lty=c(NA, NA, NA), pch = c(1, 25, 22), col=c("firebrick1","firebrick1","firebrick1"), text.col = c("firebrick1", "firebrick1", "firebrick1"), bty = "n", lwd=lwd.size+0.7, cex=lab.size, horiz =T, text.width=60)
#
## 5. Lables
axis (1,at=seq(-180,180,30), cex.axis = lab.size, lwd= lwd.box,  tck = tck.dist, mgp = c(1,dist, 0), las = 1)
# axis (3, at=seq(-180,180,30), labels = F, cex.axis = lab.size, lwd=7,  tck = tck.dist, mgp = c(1,dist, 0), las = 1)
axis (2, at = c(-75, -60, -45, -30, -15, 0, 15, 30, 45, 60, 75, 90), cex.axis = lab.size, lwd= lwd.box,  tck = tck.dist, mgp = c(1,dist, 0), las = 1)
# axis (4, at = c(-80, -60, -40, -20, 0, 20, 40, 60, 80), labels = F, cex.axis = lab.size, lwd=7,  tck = tck.dist, mgp = c(1,dist, 0), las = 1)
box(lwd = lwd.box)
#
mtext(expression("Longitude ("*degree*")"), side = 1,  line = 1.2, cex = lab.size)
mtext(expression("Latitude ("*degree*")"), side=2, line=1.2,cex=lab.size)
#
dev.off()







###  ===============================================================
### FIG. 1B, C | PRESENCE RECORDS IN SPACE-TIME
###  ===============================================================

# Lables
cea <- 1
lab.size <- 1
lab.dist <- 0.45
tck.dist <- -0.007
xlab.dist <- 2.5
ylab.dist <- 2.5
dist <- 0.25 # distance of numbers (axis lables) from axis
# Points
pt.sym <- 1
pt.col <- "black"
# Lines
lwd.box <- 10
# Continents
c.col <- "grey70"

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------
## *** FIGURE 1 B | SAMPLE LEVEL RICHNESS IN COLOR

## Get transposed data
dat <- read.csv("./1_Merge_data/Output_3/Diazo_transposed.csv")
head(dat)
# Make it workable
dat$x <- dat$decimalLongitude
dat$y <- dat$decimalLatitude
dat$taxon <- dat$scientificName

## Get richness
dat$day <- as.numeric(dat$day)
unique(dat$day)
## Calculate sample level richness
dat$spp.no <- rowSums(sign(dat[,9:dim(dat)[2]]),na.rm=T)
max(dat$spp.no) # 11
## Order data along diversity
dat.new <- dat[order(dat$spp.no, decreasing = F),] # Top line will be plotted first; higher richness last
## Exclude richness outliers
# thr <- 0.9999 # Excluded data above 99.99th percentile (currently above 147.8512)
# upper.div.limit.to.be.considered <- quantile(dat.new$spp.no, thr, na.rm = T)
# dat.new <- dat.new[ which(dat.new$spp.no <= upper.div.limit.to.be.considered),]
max(dat.new$spp.no)
hist(dat.new$spp.no, breaks = 10)

## Choose a color palette
# col.vec <- rev(colorRampPalette(brewer.pal(11,'RdYlBu')[1:11])(20)) [ as.numeric(cut(dat.rest$spp.no, breaks = 20))] # blue - neutral - red
# col.vec <- rev(colorRampPalette(brewer.pal(11,'OrRd')[1:11])(20)) [ as.numeric(cut(dat.rest$spp.no, breaks = 20))] # blue - neutral - red
# col.vec <- rev(colorRampPalette(brewer.pal(9,'YlOrRd')[1:11])(20)) [ as.numeric(cut(dat.rest$spp.no, breaks = 20))] # blue - neutral - red
# col.vec <- rev(ygobb(20)) [ as.numeric(cut(dat.new$spp.no, breaks = 20))] # dark blues to yellows
col.vec <- rev(colorRampPalette(c("#FFFFAA", "coral","royalblue4", "black") )(12)[1:11]) [ as.numeric(cut(dat.new$spp.no, breaks = 11))] # dark blues to yellows

# Assign color to richness/samples
dat.new$colvec <- col.vec

# Order by spp.no
dat.new <- dat.new[order(dat.new$spp.no, decreasing = F),]
summary(dat.new$spp.no)

## Graphical specifications
cea <- 1
lab.size <- 0.95
lab.dist <- 0.35
tck.dist <- -0.007
xlab.dist <- 2.5
ylab.dist <- 2.5
pt.size <- 0.43
pt.sym <- 1
lwd.box <- 10

################# Plot time versus latitude: diversity in color
## *** FIGURE 1 B | SAMPLE LEVEL RICHNESS IN COLOR
#
png(filename = paste0("./1_Visualize_data/Output_3/Plot_1B_revised.png"), width = 3098, height = 3333, units = "px", pointsize = 80) #orig dim 2500x2690
#
#
par(mar=c(3,3,3,1))
#
# Test of color concept
plot(dat.new$y, dat.new$month+dat.new$day*1/31, pch = pt.sym, xlab = "", ylab = "", ylim = c(1, 13.5), xlim=c(-90,90), col = dat.new$colvec, axes = F,  xaxt = "n", yaxt = "n", cex = pt.size, lwd = 9.2, type = "n")
# Add grid lines
axis (2, at = c(1:13), labels = F, cex.axis = 0.4, tck = 1, lwd = 0.5, lty = 1, mgp = c(1, lab.dist, 0), col = "grey30") # y-axis grid-lines
axis (1, at=seq(-180,180,30), labels = F, cex.axis = 0.4, tck = 1, lwd = 0.8, lty = 1, mgp = c(1, lab.dist, 0), col = "grey30") # x-axis grid-lines
# Add dots
points(dat.new$y, dat.new$month+dat.new$day*1/31, pch = pt.sym, xlab = "", ylab = "", ylim = c(1, 13.5), xlim=c(-90,90), col = dat.new$colvec, axes = F,  xaxt = "n", yaxt = "n", cex = pt.size, lwd = 9.2)
axis (1,at=seq(-90,90,15), labels = c("-90",  "-75", "-60",  "-45", "-30", "-15","0","15",  "30", "45", "60", "75", "90"), tck = tck.dist, cex.axis = lab.size, lwd = lwd.box, mgp = c(1, lab.dist, 0), las = 1)# x-axis
axis (2, at=seq(1, 12, 1), labels = c("Jan", "Feb", "Mar","Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"), cex.axis = lab.size, lwd=lwd.box,  tck = tck.dist, mgp = c(1, lab.dist, 0), las = 2) # y-axis
# axis (3,at=seq(-90,90,7.5), labels = F, tck = -tck.dist, cex.axis = nr.size, lwd = lwd.box, mgp = c(1, dist, 0), las = 1)# x-axis
axis (4, at=seq(1, 12, 1), labels =F, cex.axis = 0.4, lwd=lwd.box,  tck = tck.dist, mgp = c(1, lab.dist, 0), las = 2) # y-axis
#
box(lwd = lwd.box)
mtext(expression("Latitude ("*degree*")"), side=1, line=1.4,cex=lab.size)
# Add color brick in plot
colorbar.plot (-58.5, 13.92, strip = c(1:13), strip.width = 0.021, strip.length = 0.38, adj.x = 0.5, adj.y = 1, col= rev(colorRampPalette(c("#FFFFAA", "coral","royalblue4", "black") )(12)[1:11]) , horizontal=T)
# Add values
legend(-55.4, 13.78, paste0(1,"                                                  ","11"), xjust = 0.54, bty = 'n', cex = lab.size-0.1)

# Add header
legend(-60.5, 13.78, "Within-sample richness", xjust = 0.54, bty = 'n', cex = lab.size-0.1)	# Close png
#
# -------------------- SUBPLOT CHARACTER
# text(-114, 13.85, c("b"), cex = cea+0.2, xpd = NA, font = 2) # par(legend( char.dist+1.6, 35, "a", bty = 'n', cex = cea+0.1), font = 2, xpd=F)
dev.off()
#

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## *** SINGLE COLOR HISTOGRAM SUPPLEMENT TO SUPPLEMENTARY FIGURE 1B - Histogram using the integral of months
png(filename = paste0("./1_Visualize_data/Output_3/1B_hist.png"), width = 2940, height = 1500, units = "px", pointsize = 80) # original dimensions 2500 x 2690
opar <- par(lwd=12)
hist(dat.new$y, freq = F, axes = F,  main = "", ylab = "", xlab = "", ylim = c(0,0.06), xlim = c(-90, 90), breaks = 180, col = "black") # cex.axis = 0.4, mgp = c(1,0.05,0), tck = -0.05, lwd = 6)
abline(v=median(dat.new$y, na.rm=T), col ="yellow2", lwd = 10)
axis (1, labels = T, at = c(-90, -75, -60, -45, -30, -15, 0, 15, 30, 45, 60, 75, 90), tck = -tck.dist*2.5, cex.axis = lab.size, mgp = c(1,lab.dist, 0), las = 1, lwd = lwd.box) # x-axis
axis (2, labels = T,  at = c(0.0, 0.02, 0.04, 0.06), tck = tck.dist*3, cex.axis = lab.size-0.1, mgp = c(1,0.25, 0), las = 1, lwd = lwd.box) # y-axis
# axis (2, at =  c(0.0, 0.02, 0.04, 0.06, 0.08), labels = F, cex.axis = 0.4, tck = 1, lwd = 0.8, lty = 1, mgp = c(1, dist, 0), col = "grey30") # y-axis grid-lines
# axis (1, at=seq(-90,90,15), labels = F, cex.axis = 0.4, tck = 1, lwd = 0.8, lty = 1, mgp = c(1, dist, 0), col = "grey30") # x-axis grid-lines
mtext("Fraction of samples",side=2,line=1.65,cex= lab.size-0.15)
# mtext("Latitude",side=1,line=0.15,cex=0.55)
# legend(28.2, 0.002, expression(""*degree*C*""), bty = "n",  cex= inside.tp.sz, xpd =NA)
# legend(-11, 22, c("per", "cent"), bty = "n",  cex= 0.65, xpd =NA)
dev.off()

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## *** TWO-COLORED HISTOGRAM SUPPLEMENT TO SUPPLEMENTARY FIGURE 1B - reviewer's choice - Histogram using the integral of months
# opar <- par(lwd=2)
# Red (Year 2015 and 2017)
# hist(dat.new$y, freq = F, axes = F, border = "grey60",  main = "", ylab = "", xlab = "", ylim = c(0,0.06), xlim = c(-90, 90), breaks = 180, col = "red2") # cex.axis = 0.4, mgp = c(1,0.05,0), tck = -0.05, lwd = 6)
# axis (1, labels = F, at = c(-90, -75, -60, -45, -30, -15, 0, 15, 30, 45, 60, 75, 90), tck = -tck.dist*2.5, cex.axis = lab.size, mgp = c(1,dist, 0), las = 1, lwd = lwd.box) # x-axis
# axis (2, labels = T,  at = c(0.0, 0.02, 0.04, 0.06), tck = tck.dist*3, cex.axis = lab.size-0.1, mgp = c(1,0.25, 0), las = 1, lwd = lwd.box) # y-axis
# axis (2, at =  c(0.0, 0.02, 0.04, 0.06, 0.08), labels = F, cex.axis = 0.4, tck = 1, lwd = 0.8, lty = 1, mgp = c(1, dist, 0), col = "grey30") # y-axis grid-lines
# axis (1, at=seq(-90,90,15), labels = F, cex.axis = 0.4, tck = 1, lwd = 0.8, lty = 1, mgp = c(1, dist, 0), col = "grey30") # x-axis grid-lines
# mtext("Samples portion",side=2,line=1.91,cex= lab.size)
# mtext("Latitude",side=1,line=0.15,cex=0.55)
# legend(28.2, 0.002, expression(""*degree*C*""), bty = "n",  cex= inside.tp.sz, xpd =NA)
# legend(-11, 22, c("per", "cent"), bty = "n",  cex= 0.65, xpd =NA)
# Blue on top (Year 2015 only)
# hist(dat.new[which(dat.new$yearOfDataAccess == "2015" | dat.new$yearOfDataAccess == "2015_2017"  ), ]$y, freq = F, axes = F, border = "black",  main = "", ylab = "", xlab = "", ylim = c(0,0.06), xlim = c(-90, 90), breaks = 180, col = "darkblue", add = T) #
# dev.off()

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## *** FIGURE 1 C | SAMPLE LEVEL RICHNESS IN COLOR
#
png(filename = paste0("./1_Visualize_data/Output_3/Plot_1C_revised.png"), width = 3098, height = 3333, units = "px", pointsize = 80) #orig dim 2500x2690
#
#
par(mar=c(3,3,3,1))
#
#
plot(dat.new$x, dat.new$month+dat.new$day*1/31, pch = 1, xlab = "", ylab = "", ylim = c(1, 13.5), xlim=c(-180,180), col = dat.new$colvec, axes = F,  xaxt = "n", yaxt = "n", cex = pt.size, lwd = 9.2, type = "n")
# Add grid lines
axis (2, at = c(1:13), labels = F, cex.axis = 0.4, tck = 1, lwd = 0.5, lty = 1, mgp = c(1, lab.dist, 0), col = "grey30") # y-axis grid-lines
axis (1, at=seq(-180,180,30), labels = F, cex.axis = 0.4, tck = 1, lwd = 0.8, lty = 1, mgp = c(1, lab.dist, 0), col = "grey30") # x-axis grid-lines
# Add dots
points(dat.new$x, dat.new$month+dat.new$day*1/31, pch = 1, xlab = "", ylab = "", ylim = c(1, 13.5), xlim=c(-180,180), col = dat.new$colvec, axes = F,  xaxt = "n", yaxt = "n", cex = pt.size, lwd = 9.2)
axis (1,at=seq(-180,180,30), labels = NULL, tck = tck.dist, cex.axis = lab.size, lwd = lwd.box, mgp = c(1, lab.dist, 0), las = 1)# x-axis
axis (2, at=seq(1, 12, 1), labels = F, cex.axis = lab.size, lwd=lwd.box,  tck = tck.dist, mgp = c(1, lab.dist, 0), las = 2) # y-axis
# axis (3,at=seq(-180,180,30), labels = F, tck = -0.02, cex.axis = nr.size, lwd = lwd.box, mgp = c(1, dist, 0), las = 1)# x-axis
# axis (4, at=seq(1, 13, 1), labels =F, cex.axis = nr.size, lwd=lwd.box,  tck = -0.02, mgp = c(1, dist, 0), las = 2) # y-axis
box(lwd = lwd.box)
mtext(expression("Longitude ("*degree*")"), side=1, line=1.4, cex=lab.size)
# mtext("Time (month of the year)", side = 2, line = 1.8, cex = nr.size)
# Add color brick in plot
colorbar.plot (-117, 13.92, strip = c(min(dat.new$spp.no, na.rm=T):max(dat.new$spp.no, na.rm=T)), strip.width = 0.021, strip.length = 0.38, adj.x = 0.5, adj.y = 1,col= rev(colorRampPalette(c("#FFFFAA","coral", "royalblue4", "black") )(12)[1:11]) , horizontal=T)
# Add values
legend(-110.7, 13.78, paste0(1,"                                                  ","11"), xjust = 0.54, bty = 'n', cex = lab.size-0.1)
# Add header
legend(-119, 13.78, "Within-sample richness", xjust = 0.54, bty = 'n', cex = lab.size-0.1)	# Close png
#
# -------------------- SUBPLOT CHARACTER
# text(-203, 13.85, c(""), cex = cea+0.2, xpd = NA, font = 2) # par(legend( char.dist+1.6, 35, "a", bty = 'n', cex = cea+0.1), font = 2, xpd=F)
dev.off()

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## *** SINGLE COLOR HISTOGRAM SUPPLEMENT TO SUPPLEMENTARY FIGURE 1C - Histogram using the integral of months
png(filename = paste0("./1_Visualize_data/Output_3/1C_hist.png"), width = 2940, height = 1500, units = "px", pointsize = 80) # original dimensions 2500 x 2690
opar <- par(lwd=12)
# Red
hist(dat.new$x, freq = F, axes = F,  main = "", ylab = "", xlab = "", ylim = c(0,0.06), xlim = c(-180, 180), breaks = 180, col = "black") # cex.axis = 0.4, mgp = c(1,0.05,0), tck = -0.05, lwd = 6)
abline(v=median(dat.new$x, na.rm=T), col ="yellow2", lwd = 10)
axis (1, labels = T, at = c(-180, -150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150, 180), tck = -tck.dist*2.5, cex.axis = lab.size, mgp = c(1,dist, 0), las = 1, lwd = lwd.box) # x-axis
axis (2, labels = T,  at = c(0.0, 0.02, 0.04, 0.06), tck = tck.dist*3, cex.axis = lab.size-0.1, mgp = c(1,0.25, 0), las = 1, lwd = lwd.box) # y-axis
# axis (2, at =  c(0.0, 0.02, 0.04, 0.06, 0.08), labels = F, cex.axis = 0.4, tck = 1, lwd = 0.8, lty = 1, mgp = c(1, dist, 0), col = "grey30") # y-axis grid-lines
# axis (1, at=seq(-180,180,30), labels = F, cex.axis = 0.4, tck = 1, lwd = 0.8, lty = 1, mgp = c(1, dist, 0), col = "grey30") # x-axis grid-lines
mtext("Fraction of samples",side=2,line=1.8,cex= inside.tp.sz)
# mtext("Longitude",side=1,line=0.15,cex=0.55)
# legend(28.2, 0.002, expression(""*degree*C*""), bty = "n",  cex= inside.tp.sz, xpd =NA)
# legend(-11, 22, c("per", "cent"), bty = "n",  cex= 0.65, xpd =NA)
# Blue on top
#
dev.off()

## *** TWO-COLORED HISTOGRAM SUPPLEMENT TO SUPPLEMENTARY FIGURE 1B - Histogram using the integral of months
# opar <- par(lwd=2)
# Red
# hist(dat.new$x, freq = F, axes = F, border = "grey60",  main = "", ylab = "", xlab = "", ylim = c(0,0.06), xlim = c(-180, 180), breaks = 180, col = "red2") # cex.axis = 0.4, mgp = c(1,0.05,0), tck = -0.05, lwd = 6)
# axis (1, labels = F, at = c(-180, -150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150, 180), tck = -tck.dist, cex.axis = lab.size, mgp = c(1,dist, 0), las = 1, lwd = lwd.box) # x-axis
# axis (2, labels = T,  at = c(0.0, 0.02, 0.04, 0.06), tck = tck.dist*3, cex.axis = lab.size-0.1, mgp = c(1,0.25, 0), las = 1, lwd = lwd.box) # y-axis
# axis (2, at =  c(0.0, 0.02, 0.04, 0.06, 0.08), labels = F, cex.axis = 0.4, tck = 1, lwd = 0.8, lty = 1, mgp = c(1, dist, 0), col = "grey30") # y-axis grid-lines
# axis (1, at=seq(-180,180,30), labels = F, cex.axis = 0.4, tck = 1, lwd = 0.8, lty = 1, mgp = c(1, dist, 0), col = "grey30") # x-axis grid-lines
# mtext("Presences portion",side=2,line=1.8,cex= inside.tp.sz)
# mtext("Latitude",side=1,line=0.15,cex=0.55)
# legend(28.2, 0.002, expression(""*degree*C*""), bty = "n",  cex= inside.tp.sz, xpd =NA)
# legend(-11, 22, c("per", "cent"), bty = "n",  cex= 0.65, xpd =NA)
# Blue on top
# hist(dat.new[which(dat.new$yearOfDataAccess == "2015" | dat.new$yearOfDataAccess == "2015_2017"  ), ]$x, freq = F, axes = F, border = "black",  main = "", ylab = "", xlab = "", ylim = c(0,0.06), xlim = c(-180, 180), breaks = 180, col = "darkblue", add = T) # cex.axis = 0.4, mgp = c(1,0.05,0), tck = -0.05, lwd = 6)
#
# dev.off()

# Fig. 1D Additional Histogram: YEAR OF SAMPLING -----------------------------------------------------------------------------------------------------------------------------------------------
## *** UPDATED VERSION - SINGLE COLOR HISTOGRAM SUPPLEMENT TO SUPPLEMENTARY FIGURE 1C - reviewer's choice - Histogram using the integral of months
png(filename = paste0("./1_Visualize_data/Output_3/hist_year.png"), width = 2900, height = 1400, units = "px", pointsize = 80)
opar <- par(lwd=12)
# Histogram
hist(as.numeric(dat.new[which(dat.new$year>1800),]$year), freq = F, axes = F, main = "", ylab = "", xlab = "", ylim = c(0,0.15), xlim = c(1800, 2017), breaks = 500, col = "black") # cex.axis = 0.4, mgp = c(1,0.05,0), tck = -0.05, lwd = 6)
abline(v=median(dat.new$year, na.rm=T), col ="yellow2", lwd = 10)
axis (1, labels = T, at = c(1800, 1850, 1900, 1950, 2000, 2017), tck = -tck.dist*2.5, cex.axis = lab.size-0.15, mgp = c(1,0.155, 0), las = 1, lwd = lwd.box) # x-axis
axis (2, labels = T,  at = c(0.0, 0.05, 0.1, 0.15), tck = tck.dist*3, cex.axis = lab.size-0.1, mgp = c(1,0.2, 0), las = 1, lwd = lwd.box) # y-axis
# axis (2, at =  c(0.0, 0.05, 0.1, 0.15), labels = F, cex.axis = 0.4, tck = 1, lwd = 0.8, lty = 1, mgp = c(1, dist, 0), col = "grey30") # y-axis grid-lines
# axis (1, at=seq(-180,180,30), labels = F, cex.axis = 0.4, tck = 1, lwd = 0.8, lty = 1, mgp = c(1, dist, 0), col = "grey30") # x-axis grid-lines
mtext(" Fraction of samples",side=2,line=1.65,cex= lab.size-0.15)
mtext("Year",side=1,line=0.86,cex=lab.size-0.19)
# legend(28.2, 0.002, expression(""*degree*C*""), bty = "n",  cex= inside.tp.sz, xpd =NA)
# legend(-11, 22, c("per", "cent"), bty = "n",  cex= 0.65, xpd =NA)
# Blue on top
#
dev.off()

# Fig. 1E Additional Histogram: DEPTH OF SAMPLING -----------------------------------------------------------------------------------------------------------------------------------------------
## *** UPDATED VERSION - SINGLE COLOR HISTOGRAM SUPPLEMENT TO SUPPLEMENTARY FIGURE 1C - reviewer's choice - Histogram using the integral of months
png(filename = paste0("./1_Visualize_data/Output_3/hist_depth.png"), width = 3000, height = 1500, units = "px", pointsize = 80)
opar <- par(lwd=12)
# Histogram
hist(as.numeric(dat.new[which(dat.new$depth>0),]$depth), freq = F, axes = F, main = "", ylab = "", xlab = "", ylim = c(0,0.75), xlim = c(0, 400), breaks = 6000, col = "black") # cex.axis = 0.4, mgp = c(1,0.05,0), tck = -0.05, lwd = 6)
abline(v=median(dat.new$depth, na.rm=T), col ="yellow2", lwd = 10)
axis (1, labels = T, at = c(0, 50, 100, 150, 200, 250, 300, 350, 400), tck = -tck.dist*2.5, cex.axis = lab.size-0.07, mgp = c(1,0.22, 0), las = 1, lwd = lwd.box) # x-axis
axis (2, labels = T,  at = c(0.0,  0.25, 0.5, 0.75), tck = tck.dist*3, cex.axis = lab.size-0.07, mgp = c(1,0.25, 0), las = 1, lwd = lwd.box) # y-axis
# axis (2, at =  c(0.0, 0.02, 0.04, 0.06, 0.08), labels = F, cex.axis = 0.4, tck = 1, lwd = 0.8, lty = 1, mgp = c(1, dist, 0), col = "grey30") # y-axis grid-lines
# axis (1, at=seq(-180,180,30), labels = F, cex.axis = 0.4, tck = 1, lwd = 0.8, lty = 1, mgp = c(1, dist, 0), col = "grey30") # x-axis grid-lines
# mtext("Fraction of samples",side=2,line=1.65,cex= lab.size-0.15)
mtext("Depth (m)",side=1,line=1.02,cex=lab.size-0.12)
# legend(28.2, 0.002, expression(""*degree*C*""), bty = "n",  cex= inside.tp.sz, xpd =NA)
# legend(-11, 22, c("per", "cent"), bty = "n",  cex= 0.65, xpd =NA)
# Blue on top
#
dev.off()



### ========================================================
### SUPPLEMENTARY FIGURE | GLOBAL DISTRIBUTION OF SAMPLED RICHNESS
### ========================================================

## Varia
cea <- 1.2 # axis lables type size
sz <- 0.3 # axes in plot
# c.col <- "grey86" # continents color
c.col <- "#D6EFED" # continents color
c.col <- "grey63" # continents color
s.col <- "white" # s line color (?)
col.pts <- "black" # color of points
col.pts <- "navy" # color of points
col.pts <- "black" # color of points
point <- 1
point <- 16 # point symbol
col.line <- "firebrick1"
col.line <- "grey65" # color of help-lines
##
sz <- 0.3 # axes numbers in plot
nr.size <- 0.3 # axes numbers in lat grad plot
# c.col <- "grey86"
c.col <- "#D6EFED"
s.col <- "white"
col.pts <- "navy"
col.pts <- "black"
pt.sym <- 16
lwd.box <- 15
dist <- 0.5
##
cea <- 1
lab.size <- 0.93
tck.dist <- -0.006 # originally -0.007
dist <- 0.22 # distance of numbers (axis lables) from axis
##
pt.size <- 0.7 # size of points
lwd.size <- 6 # width of lines of points
pt.sym <- 1
pt.col <- "black"
##
lwd.box <- 10
c.col <- "black"
##
col.hot <- "lemonchiffon3"
col.cold <- colorRampPalette(c("#B7C8DE",              colorRampPalette(c("black", "royalblue"))(5)[2])                        ) (20)[2]

### PLOT
##
png (filename = paste0("./1_Visualize_data/Output_3/Sampled_RICHNESS_in_color.png"), width = 6100, height = 4100, units = "px", pointsize = 80) # original width 5600, original height 3600
## 0. Empty plot
plot(dat.new$x, dat$y,  pch = pt.sym, xlab="",  ylab="",  ylim = c(-77,84), xlim = c (-167, 167), col = pt.col, axes = F,  xaxt = "n", yaxt = "n", cex = pt.size, type ="n")
## 1. Optional: add some background shading/mask on the ocean regions: tropical, temperate and polar: the last one of the three is standard
# plot (temp.levels, add = T, legend = F, col = c(alpha("#003366", fact) , alpha("#CCCCCC", fact), alpha("#FFCC00", fact))) #2
# plot (temp.levels, add = T, legend = F, col = c( alpha(col.cold, 0.22) , alpha("snow", 1), alpha(col.hot, 0.3) )) #2
plot (temp.levels, add = T, legend = F, col = c(alpha("#003366", 0.093), alpha("grey97", 1), alpha(col.hot, 0.335) )) # 0.1, 1, 0.4 originally
## 2. Grid lines
axis (1, at = seq(-180, 180 ,30), labels = F, cex.axis = 0.4, tck = 1, lwd = 0.8, lty = 1, mgp = c(1, dist, 0), col = "grey30") # y-axis grid-lines
axis (2, at=seq(-90,90,15), labels = F, cex.axis = 0.4, tck = 1, lwd = 0.8, lty = 1, mgp = c(1, dist, 0), col = "grey30") # x-axis grid-lines
# axis (2, at=seq(-90,90,10), labels = F, cex.axis = 0.4, tck = 1, lwd = 0.8, lty = 1, mgp = c(1, dist, 0), col = "grey30") # x-axis grid-lines, originally I took this
## 3. Map
map ("world", interior=TRUE, fill = T, lwd = 0.5, boundary=TRUE, col="black", col.fill = "black",add = T)
# shape <- readOGR(dsn = "./shapefiles_world/ne_110m_coastline/ne_110m_coastline.shp")
# shape # it seems to have the same projection..
# plot(shape, add = T, lwd = 3)
## 4. Points: ordered by richness
points(dat.new$x, dat.new$y, pch = 1, xlab = "", ylab = "", col = dat.new$colvec, axes = F,  xaxt = "n", yaxt = "n", cex = pt.size, lwd = 9.2, add = T)
## 5. Lables
axis (1,at=seq(-180,180,30), cex.axis = lab.size, lwd= lwd.box,  tck = tck.dist, mgp = c(1,dist, 0), las = 1)
## axis (3, at=seq(-180,180,30), labels = F, cex.axis = lab.size, lwd=7,  tck = tck.dist, mgp = c(1,dist, 0), las = 1)
axis (2, at = c(-75, -60, -45, -30, -15, 0, 15, 30, 45, 60, 75, 90), cex.axis = lab.size, lwd= lwd.box,  tck = tck.dist, mgp = c(1,dist, 0), las = 1)
# axis (4, at = c(-80, -60, -40, -20, 0, 20, 40, 60, 80), labels = F, cex.axis = lab.size, lwd=7,  tck = tck.dist, mgp = c(1,dist, 0), las = 1)
box(lwd = lwd.box)
##
mtext(expression("Longitude ("*degree*")"), side = 1,  line = 1.2, cex = lab.size)
mtext(expression("Latitude ("*degree*")"), side=2, line=1.2,cex=lab.size)
##
## Add color brick in plot
colorbar.plot (-137, 90, strip = c(min(dat.new$spp.no, na.rm=T):max(dat.new$spp.no, na.rm=T)), strip.width = 0.021, strip.length = 0.38, adj.x = 0.5, adj.y = 1,col= rev(colorRampPalette(c("#FFFFAA","coral", "royalblue4", "black") )(23)[1:20]) , horizontal=T)
## Add values
legend(-136, 89, paste0(1,"                                                         ","11"), xjust = 0.54, bty = 'n', cex = lab.size-0.1)
## Add header
legend(-140.7, 89.2, "Within-sample richness", xjust = 0.54, bty = 'n', cex = lab.size-0.1)	# Close png
##
dev.off()
