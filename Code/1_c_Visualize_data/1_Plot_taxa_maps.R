

wd_data <- "D:/Research/NOTION/Righetti_et_al_notion/1_Merge_data/Output_1/"
wd_data <- "D:/Research/NOTION/Righetti_et_al_notion/1_Merge_data/Output_1/"

setwd(wd_data)

df <- read.csv("Diazo.csv")
dim(df) # 21954 55
names(df)

df_taxa <- read.csv("D:/Research/NOTION/Righetti_et_al_notion/1_Merge_data/Output_2/Diazo_hierarchic_taxa.csv")
dim(df_taxa) # 30528 55
names(df_taxa)

###  PUT THE MERGED RECORDS  ON THE MAP ###

# Date: 2020-06-02
# Author: Dominic Eriksson

### ===================================================================================
### Preparatory
### ===================================================================================
rm(list = ls())
library(maps); library(ncdf4); library(RColorBrewer); library(colorRamps); library(raster); library(ncdf4)
setwd("D:/Research/NOTION/") # User one
# setwd("~/Desktop/Data_notion") # User two

### ===================================================================================
### Get data
### ===================================================================================

# Import
tt <- read.csv("./Righetti_et_al_notion/1_Merge_data/Output_2/Diazo_hierarchic_taxa.csv")

# Make dataset workable for plotting by renaming columns
tt$x <- tt$decimalLongitude
tt$y <- tt$decimalLatitude
tt$taxon <- tt$scientificName

# Adjust class of columns
tt$x <- as.numeric(tt$x)
tt$y <- as.numeric(tt$y)
tt$month <- as.numeric(as.character(tt$month))

# Grid the data
tt$x<-round(tt$x +0.5)-0.5
tt$y<-round(tt$y +0.5)-0.5

# Remove duplicates
sum(duplicated(tt[ , c("scientificName", "decimalLongitude", "decimalLatitude", "family", "genus", "year", "month", "day", "depth")])) # 5919
tt.new <- tt[!duplicated(tt[ , c("scientificName", "decimalLongitude", "decimalLatitude", "family", "genus", "year", "month", "day", "depth")]) , ]
dim(tt.new) # 24609 58

# Take a glance
head(tt.new)
dim(tt.new) # 22955 x 58
length(unique(tt$taxon)) # 57

# Environmental data, is a list that contains rasters with the env. variables of intrest from Jan to Dec
env.stack <- get(load("D:/Research/NOTION/Righetti_et_al_notion/2_Environmental_data/VarSet07.2016.RData"))
# class(env.stack[[1]])
# r <- env.stack[[1]]
# dim(r)
# names(r)

# Apply an open ocean mask
dat <- tt
llCRS <- CRS("+proj=longlat +ellps=WGS84")
crd.index <- which(colnames(dat)%in%c(c("x", "y")))
			## Convert data points to SpatialPointsDataFrame
			spdf <- SpatialPointsDataFrame(dat[,crd.index] , dat[,-crd.index], proj4string = llCRS) # coords not shown, only as data.frame: spdf.test <- as.data.frame(spdf)
			## Loop across months: mask the observations and prepare list to extract points, column-bind with mask, keep rows (points) with mask equal 1 (ocean condition)
			mups<-list()
 		   	for (k in 1:12){
    				# 1. Define monthly subset, for monthly match-up with mask
    				if(sum(spdf$month==k)>0){ss<-subset(spdf,spdf$month==k)} else {next}
    				# 2. Extract data points overlapping with "mask"
    				extr<-extract(x=env.stack[[k]],y=ss,method="simple")# returns the values of a raster object (msk) for cells in which points fall
    				extr<-as.data.frame(extr)#??the mask values are converted to data frame
    				ss.new <-cbind(as.data.frame(ss),extr)[which(extr$Mask!=0),] # CORE OPERATION: SELECT ROWS / DATA for which mask == 1
    				# 3. Store monthly specific ss.new in list (for later row-bind)
    				mups[[length(mups)+1]]<-ss.new
    				print(paste(k))
 			}

# Merge match-ups
tot <-do.call("rbind",mups)

# Useful gridded 1° by 1°, monthly data records for modelling
dim(tot) # 24815   106


### ===================================================================================
### 1. Visualize total geographic data distribution (put on map)
### ===================================================================================

# Graphical specifications
cea <- 1.2 # axis lables type size
sz <- 0.3 # axes in plot
col.pts <- "black" # color of points
point <- 16 # point symbol
col.line <- "grey65" # color of help-lines
#
sz <- 0.3 # axes numbers in plot
nr.size <- 0.3 # axes numbers in lat grad plot
s.col <- "white"
col.pts <- "black"
pt.sym <- 16
lwd.box <- 15
dist <- 0.5
#
cea <- 1 # lables
lab.size <- 0.93
tck.dist <- -0.006 # originally -0.007
dist <- 0.22 # distance of numbers (axis lables) from axis
#
pt.size <- 0.8 # size of points
lwd.size <- 12 # width of lines of points
pt.sym <- 1
pt.col <- "black"
#
lwd.box <- 10 # Lines
c.col <- "#D6EFED" # Continents
c.col <- "black" # Continents
col.sh <- colorRampPalette(c("#B7C8DE",  colorRampPalette(c("black", "royalblue"))(5)[2])                        ) (20)[3]

## PLOT
png (filename = paste0("./Righetti_et_al_notion/1_Visualize_data/Output_1DE/1 A Total_data_",nrow(tot),"_records.png"), width = 6100, height = 4100, units = "px", pointsize = 80) # original width 5600, original height 3600

# Empty plot
plot(tot$x, tot$y,  pch = pt.sym, xlab="",  ylab="",  ylim = c(-77,84), xlim = c (-167, 167), col = pt.col, axes = F,  xaxt = "n", yaxt = "n", cex = pt.size, type ="n")
# Grid lines
axis (1, at = seq(-180, 180 ,30), labels = F, cex.axis = 0.4, tck = 1, lwd = 0.8, lty = 1, mgp = c(1, dist, 0), col = "grey30") # y-axis grid-lines
axis (2, at=seq(-90,90,15), labels = F, cex.axis = 0.4, tck = 1, lwd = 0.8, lty = 1, mgp = c(1, dist, 0), col = "grey30") # x-axis grid-lines
map ("world", interior=TRUE, fill = T, lwd = 0.5, boundary=TRUE, col="black", col.fill = "black",add = T)

# 1. Luo
dt <- tot[which(tot$set=="Luo"),]
points(dt$x, dt$y, col = "orange", cex = pt.size, lwd = lwd.size+0.7) # originally greenyellow
# 2. Tang and Cassar
dt <- tot[which(tot$set=="Tan"),]
points(dt$x, dt$y, col = "darkblue", cex = pt.size, lwd = lwd.size+0.7) # originally greenyellow
# 3. Gradoville
dt <- tot[which(tot$set=="Gra"),]
points(dt$x, dt$y, col = "firebrick2", cex = pt.size, lwd = lwd.size+0.7) # originally greenyellow
# 4. Phy
dt <-  tot[which(tot$set=="Phy"),]
points(dt$x, dt$y, col = "yellow", cex = pt.size, lwd = lwd.size+0.7)
# 5. GBIF
dt <-  tot[which(tot$set=="Gbi"),]
points(dt$x, dt$y, col = "green", cex = pt.size, lwd = lwd.size+0.7)
# 6. OBIS
dt <-  tot[which(tot$set=="Obi"),]
points(dt$x, dt$y, col = "pink", cex = pt.size, lwd = lwd.size+0.7)
# 7. Tara OTU based
dt <-  tot[which(tot$set=="TarOTU"),]  # *** COMMENT: Attention, something is wrong with tara ocean records, either the month is missing, or something else...
points(dt$x, dt$y, col = "slateblue", cex = pt.size+0.2, lwd = lwd.size+0.7, pch = 16)
# 8. Mar
dt <-  tot[which(tot$set=="Mar"),]  # *** COMMENT: Attention, something is wrong with tara ocean records, either the month is missing, or something else...
points(dt$x, dt$y, col = "firebrick3", cex = pt.size+0.2, lwd = lwd.size+0.7, pch = 16)

## Legend
legend(-160, 91,
	c("Luo", "Tang and Cassar", "Gradoville", "PhytoBase", "Gbif", "Obis", "Tara OTU", "M96"),
	lty=c(1, 1, 1, 1, 1, 1, 1, 1),
	col=c("orange", "darkblue", "firebrick2", "yellow", "green", "pink", "slateblue", "firebrick3" ) # "coral",  "lightskyblue", "darkgreen"
	# col=c("grey5", "darkblue", "red2", "greenyellow" ) # "coral",  "lightskyblue", "darkgreen"
	, bty = "n", lwd=20, cex= 0.75, horiz =T)

## Axes
axis (1,at=seq(-180,180,30), cex.axis = lab.size, lwd= lwd.box,  tck = tck.dist, mgp = c(1,dist, 0), las = 1)
axis (2, at = c(-75, -60, -45, -30, -15, 0, 15, 30, 45, 60, 75, 90), cex.axis = lab.size, lwd= lwd.box,  tck = tck.dist, mgp = c(1,dist, 0), las = 1)
box(lwd = lwd.box)
mtext(expression("Longitude ("*degree*")"), side = 1,  line = 1.2, cex = lab.size)
mtext(expression("Latitude ("*degree*")"), side=2, line=1.2, cex=lab.size)

dev.off()



### ===================================================================================
### 1. Visualize taxon-wise  geographic data distribution (put on map)
### ===================================================================================

dat.spl<- split(tot, f = tot$scientificName)
df.tax.obs <- data.frame("taxon" = names(unlist(lapply(dat.spl, nrow))), "obs" = as.numeric(unlist(lapply(dat.spl, nrow))) )
df.tax.obs <- df.tax.obs[order(df.tax.obs$obs, decreasing = T),]
rownames(df.tax.obs) <- 1:nrow(df.tax.obs)
df.tax.obs

### B. Plot taxa one by one ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
for ( i in 1:length(unique(tot$scientificName)) ){
  #
  png (filename = paste0("./1_Visualize_data/Output_1/", rownames(df.tax.obs[which(df.tax.obs$taxon==unique(tot$scientificName)[i]),])," ",unique(tot$scientificName)[i],".png"), width = 6100, height = 4100, units = "px", pointsize = 80)
  #
  ## 1. Empty plot
  dat.taxon <- tot[which(tot$scientificName == unique(tot$scientificName)[i]), ]
  plot(dat.taxon$decimalLongitude, dat.taxon$decimalLatitude,  pch = pt.sym, xlab="",  ylab="",  ylim = c(-77,84), xlim = c (-167, 167), col = pt.col, axes = F,  xaxt = "n", yaxt = "n", cex = pt.size, type ="n")
  ## 2. Grid lines
  axis (1, at = seq(-180, 180 ,30), labels = F, cex.axis = 0.4, tck = 1, lwd = 0.8, lty = 1, mgp = c(1, dist, 0), col = "grey30") # y-axis grid-lines
  axis (2, at=seq(-90,90,15), labels = F, cex.axis = 0.4, tck = 1, lwd = 0.8, lty = 1, mgp = c(1, dist, 0), col = "grey30") # x-axis grid-lines
  ## 3. Map
  map ("world", interior=TRUE, fill = T, lwd = 0.5, boundary=TRUE, col="black", col.fill = "black",add = T)
  ## 4. Points: Total
   col.sh <- colorRampPalette(c("#B7C8DE",  colorRampPalette(c("black", "royalblue"))(5)[2])                        ) (20)[3]
  points(dat.taxon$decimalLongitude, dat.taxon$decimalLatitude, col = "red", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 16) # originally greenyellow

 ## Legend
legend(-160, 91,
	c(as.character(unique(tot$scientificName)[i]) , paste0("No. records = ",nrow(dat.taxon) ) , paste0("No. presences = ",nrow(dat.taxon[which(dat.taxon$occurrenceStatus == "PRESENT"), ])) ),
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

  dev.off()

  print(paste(i))  # Progress
}
