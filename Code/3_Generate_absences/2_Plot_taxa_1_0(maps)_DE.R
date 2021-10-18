## Put the presences and pseudoabsences records on the map

# Date 27.09.2021

screen -r 45199.pts-4.kryo

# Author: Dominic Eriksson, This script was developed by Damiano Righetti

### 1. Preparatory
rm(list = ls())
library(maps); library(ncdf4); library(RColorBrewer); library(colorRamps); library(raster); library(ncdf4)
setwd("/net/kryo/work/deriksson/notion")

# Define data types * background strategies
nms.sets <- c(
"pres_abs,tot_bg_overl",
"counts_abs,tot_bg_overl",
"pres_abs,gr_bg_overl",
"counts_abs,gr_bg_overl",
"pres_abs,cr_bg_overl",
"counts_abs,cr_bg_overl",
"pres_abs,gr_bg_nonov",
"counts_abs,gr_bg_nonov",
"pres_abs,cr_bg_nonov",
"counts_abs,cr_bg_nonov"
)

# Reference list of tqxa, which are to be plotted (baseline scenario total background, presence absence dataset)
dat.baseline <- get(load("./4_Generate_absences/Output_1/1 pres_abs,tot_bg_overl(T_MLD1).RData"))
df.tax.obs <- data.frame("taxon" = names(dat.baseline), "obs" = unlist(lapply(dat.baseline, function(x) {if(is.null(x)) {x1 <- 0} else {x1 <- nrow(x)}; return(x1)})))
rownames(df.tax.obs) <- 1:nrow(df.tax.obs)
df.tax.obs

## 2. Graphical specifications

# Graphical specifications
cea <- 1.2 # axis labels type size
sz <- 0.3 # axes in plot
col.pts <- "black" # color of points
point <- 16 # point symbol
col.line <- "grey65" # color of help-lines
#
sz <- 0.3 # yxes numbers in plot
nr.size <- 0.3 # axes numbers in lat grad plot
s.col <- "white"
col.pts <- "black"
pt.sym <- 16
lwd.box <- 15
dist <- 0.5
#
cea <- 1 # labels
lab.size <- 0.93
tck.dist <- -0.006 # originally -0.007
dist <- 0.22 # distance of number (axis lables) from axis
#
pt.size <- 0.8 # size of the points
lwd.size <- 12 # width of lines of points
pt.sym <- 1
pt.col <- "black"

lwd.box <- 10 # lines
c.col <- "#D6EFED" # Continents
c.col <- "black" # Continents
col.sh <- colorRampPalette(c("#B7C8DE", colorRampPalette(c("black", "royalblue"))(5)[2]))

### 3. Visualize taxon-wise geographic data distribution (put on map)
# Loop across data types and/or background strategies
for(k in 1:19){

  if(k == 1){tt <- get(load("./4_Generate_absences/Output_1/1 pres_abs,tot_bg_overl(T_MLD1).RData"))} # Presence-absence, tot bg, overlapping bg
  if(k == 2){tt <- get(load("./4_Generate_absences/Output_1/2 counts_abs,tot_bg_overl(T_MLD1).RData"))} # Counts-absence, tot bg, overlapping bg
  if(k == 3){tt <- get(load("./4_Generate_absences/Output_1/3 pres_abs,gr_bg_overl(T_MLD1).RData"))} # Presence-absence, tot bg, overlapping bg
  if(k == 4){tt <- get(load("./4_Generate_absences/Output_1/4 counts_abs,gr_bg_overl(T_MLD1).RData"))} # Counts-absence, tot bg, overlapping bg
  if(k == 5){tt <- get(load("./4_Generate_absences/Output_1/5 pres_abs,cr_bg_overl(T_MLD1).RData"))} # Presence-absence, tot bg, overlapping bg
  if(k == 6){tt <- get(load("./4_Generate_absences/Output_1/6 counts_abs,cr_bg_overl(T_MLD1).RData"))} # Counts-absence, cr bg, overlapping bg
  if(k == 7){tt <- get(load("./4_Generate_absences/Output_1/7 pres_abs_,gr_bg_nonov(T_MLD1).RData"))} # Presence-absence, gr bg, nonoverlapping bg
  if(k == 8){tt <- get(load("./4_Generate_absences/Output_1/8 counts_abs,gr_bg_nonov(T_MLD1).RData"))} # Counts-absence, gr bg, nonoverlapping bg
  if(k == 9){tt <- get(load("./4_Generate_absences/Output_1/9 pres_abs,cr_bg_nonov(T_MLD1).RData"))} # Presence-absence, cr bg, nonoverlapping bg
  if(k == 10){tt <- get(load("./4_Generate_absences/Output_1/10 counts_abs,cr_bg_nonov(T_MLD1).RData"))} # Counts-absence, cr bg, nonoverlapping bg

# Overview of taxa contained and to be plotted in strategy
  df.taxa.withdata <- data.frame("taxon" = names(unlist(lapply(tt, nrow))), "obs" = as.numeric(unlist(lapply(tt, nrow))))
  df.taxa.withdata

## Loop across taxa
for(i in 1:nrow(df.tax.obs)){

  ## Get taxon specific data
  dat.taxon <- tt[[i]]
  if(is.null(dat.taxon)){next}
s
  # Open the figure (.png)
  if(k == 1) { png (filename = paste0("./4_Generate_absences/Output_2_visualization/1_PA_tot_bg_overlapping/",as.character(df.tax.obs$taxon[i]),".png"), width = 6100, height = 4100, units = "px", pointsize = 80) }
  if(k == 2) { png (filename = paste0("./4_Generate_absences/Output_2_visualization/2_CA_tot_bg_overlapping/",as.character(df.tax.obs$taxon[i]),".png"), width = 6100, height = 4100, units = "px", pointsize = 80) }
  if(k == 3) { png (filename = paste0("./4_Generate_absences/Output_2_visualization/3_PA_gr_bg_overlapping/",as.character(df.tax.obs$taxon[i]),".png"), width = 6100, height = 4100, units = "px", pointsize = 80) }
  if(k == 4) { png (filename = paste0("./4_Generate_absences/Output_2_visualization/4_CA_gr_bg_overlapping/",as.character(df.tax.obs$taxon[i]),".png"), width = 6100, height = 4100, units = "px", pointsize = 80) }
  if(k == 5) { png (filename = paste0("./4_Generate_absences/Output_2_visualization/5_PA_cr_bg_overlapping/",as.character(df.tax.obs$taxon[i]),".png"), width = 6100, height = 4100, units = "px", pointsize = 80) }
  if(k == 6) { png (filename = paste0("./4_Generate_absences/Output_2_visualization/6_CA_cr_bg_overlapping/",as.character(df.tax.obs$taxon[i]),".png"), width = 6100, height = 4100, units = "px", pointsize = 80) }
  if(k == 7) { png (filename = paste0("./4_Generate_absences/Output_2_visualization/7_PA_gr_bg_nonoverlapping/",as.character(df.tax.obs$taxon[i]),".png"), width = 6100, height = 4100, units = "px", pointsize = 80) }
  if(k == 8) { png (filename = paste0("./4_Generate_absences/Output_2_visualization/8_CA_gr_bg_nonoverlapping/",as.character(df.tax.obs$taxon[i]),".png"), width = 6100, height = 4100, units = "px", pointsize = 80) }
  if(k == 9) { png (filename = paste0("./4_Generate_absences/Output_2_visualization/9_PA_cr_bg_nonoverlapping/",as.character(df.tax.obs$taxon[i]),".png"), width = 6100, height = 4100, units = "px", pointsize = 80) }
  if(k == 10) { png (filename = paste0("./4_Generate_absences/Output_2_visualization/10_CA_cr_bg_nonoverlapping/",as.character(df.tax.obs$taxon[i]),".png"), width = 6100, height = 4100, units = "px", pointsize = 80) }

  ## 1. Empty plot
  plot(dat.taxon[which(dat.taxon$obs ==1),]$x, dat.taxon[which(dat.taxon$obs ==1),]$y,  pch = pt.sym, xlab="",  ylab="",  ylim = c(-77,84), xlim = c (-167, 167), col = pt.col, axes = F,  xaxt = "n", yaxt = "n", cex = pt.size, type ="n")
  ## 2. Grid lines
  axis (1, at = seq(-180, 180 ,30), labels = F, cex.axis = 0.4, tck = 1, lwd = 0.8, lty = 1, mgp = c(1, dist, 0), col = "grey30") # y-axis grid-lines
  axis (2, at=seq(-90,90,15), labels = F, cex.axis = 0.4, tck = 1, lwd = 0.8, lty = 1, mgp = c(1, dist, 0), col = "grey30") # x-axis grid-lines
  ## 3. Map
  map ("world", interior=TRUE, fill = T, lwd = 0.00005, boundary=TRUE, col="grey96", col.fill = "grey96",add = T)
  ## 4. Points: Total
   col.sh <- colorRampPalette(c("#B7C8DE",  colorRampPalette(c("black", "royalblue"))(5)[2])                        ) (20)[3]
  points(dat.taxon[which(dat.taxon$obs ==0),]$x, dat.taxon[which(dat.taxon$obs ==0),]$y, col = "lightblue", cex = pt.size+0.1, lwd = lwd.size-6.4, pch = 1) # absences
  points(dat.taxon[which(dat.taxon$obs ==1),]$x, dat.taxon[which(dat.taxon$obs ==1),]$y, col = "red2", cex = pt.size+0.28, lwd = lwd.size+0.7, pch = 16) # presences
  ## Legend
  legend(-160, 91,
  	c(as.character(df.tax.obs$taxon[i]) , paste0("No. presences = ",nrow(dat.taxon[which(dat.taxon$obs ==1), ])), paste0("No. absences = ",nrow(dat.taxon[which(dat.taxon$obs == 0), ])), paste0("Ratio (abs/pres) = ",round(    nrow(dat.taxon[which(dat.taxon$obs == 0), ]    ) / nrow(dat.taxon[which(dat.taxon$obs ==1), ]), 1) )  ),
  	lty=c(NA, NA, NA),
  	col=c("orange", "darkblue", "firebrick2", "yellow", "green", "pink", "slateblue" ) # "coral",  "lightskyblue", "darkgreen"
  	# col=c("grey5", "darkblue", "red2", "greenyellow" ) # "coral",  "lightskyblue", "darkgreen"
  	, bty = "n", lwd=20, cex= 0.93, horiz =T)
    ## Axes
    axis (1,at=seq(-180,180,30), cex.axis = lab.size, lwd= lwd.box,  tck = tck.dist, mgp = c(1,dist, 0), las = 1)
    axis (2, at = c(-75, -60, -45, -30, -15, 0, 15, 30, 45, 60, 75, 90), cex.axis = lab.size, lwd= lwd.box,  tck = tck.dist, mgp = c(1,dist, 0), las = 1)
    box(lwd = lwd.box)
    mtext(expression("Longitude ("*degree*")"), side = 1,  line = 1.2, cex = lab.size)
    mtext(expression("Latitude ("*degree*")"), side=2, line=1.2,cex=lab.size)

    ## Close figure
    dev.off()
    print(paste(k))  # Progress

  }# Close loop across taxa

  } # close loop across data strategies
