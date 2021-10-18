###  PROJECT DIAZOTROPH OCCURRENCES OF THE (PRELIMINARY) MERGED DATASET ON GLOBAL GEOGRAPHIC AND ENVIRONMENTAL SPACE ###

# Date: 2020-06-02
# Author: Damiano Righetti
# Environmental Physics Group, ETH Zurich, CH

# Note: Currently we have ~12'100 gridded presences available for modelling, spanning 57 taxa or their variants (genera, species, strains), of which 31 have at least 16 presences
# These include also 6 host taxa, so far (Chaetoceros, and Rhizosolenia) and genera or nifH - detected species (e.g. UCYN.A) along with strains or variants (e.g. UCYN.A1)

### ===================================================================================
### Preparatory
### ===================================================================================
rm(list = ls())
library(maps); library(RColorBrewer); library(colorRamps); library(raster); library(ncdf4); library(cluster)
setwd("D:/Research/NOTION/") # User one
setwd("~/Desktop/Data_notion") # User two

### ===================================================================================
### Get diazotroph data
### ===================================================================================

# Presence-select and grid the data (for the plots)
tt <- read.csv("./Righetti_et_al_notion/1_Merge_data/Output_2/Diazo_hierarchic_taxa.csv")
dim(tt)# 30 528 x 55 (before, with no cleaned obis dates: 18 075 x 55)
tt <- tt[which(  tt$occurrenceStatus%in%c("PRESENT", "Rare (p < 1%)")  ), ] # Critical choice: presences (including very rare presences)
tt$x <- tt$decimalLongitude
tt$y <- tt$decimalLatitude
tt$taxon <- tt$scientificName
tt$x <- as.numeric(tt$x)
tt$y <- as.numeric(tt$y)
tt$month <- as.numeric(as.character(tt$month))
tt$x<-round(tt$x +0.5)-0.5
tt$y<-round(tt$y +0.5)-0.5

# Remove duplicates, not useful for modeling, at the aggregate level of monthly 1° latitude x monthly 1° longitude:
sum(duplicated(tt[ , c("scientificName", "x", "y", "year", "month", "day", "depth", "measurementMethod", "set")])) # 7284, here we want to include/ consider the measurementMethod and the sources, in the map plots!!
tt.new <- tt[!duplicated(tt[ , c("scientificName", "x", "y", "year", "month", "day", "depth", "measurementMethod", "set")]), ]

# Stats: total gridded (1° by 1°) monthly presences
dim(tt.new) # 22 318 x 58

# Stats: total taxa
length(unique(tt$taxon)) # 57

### ===================================================================================
### Add environmental data (without gridding the data)
### ===================================================================================

prj.stack <- brick(paste0(".","/Righetti_et_al_notion/2_Environmental_data/VarSet07.2016.grd"))
nelem <- nlayers(prj.stack)/12
nm <- gsub("\\.1", "", names(prj.stack[[1:nelem]]))
env.stack<-list() ; from<-seq(1,nlayers(prj.stack),by=nelem) ; for(q in 1:12){pre.env.stack <- stack(prj.stack[[ from[q]:(from[q]+nelem-1)]]); names(pre.env.stack) <- nm; env.stack[[q]] <- pre.env.stack}
names(env.stack) <- c("JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC")
dat <- tt.new
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

# Gridded 1° by 1°, monthly data records useful for modelling
dim(tot) # 18 155 x 106 (contains duplicates across different sources)

# Prepare measurement method
tot$measurementMethod <- as.character(tot$measurementMethod)
tot[is.na(tot$measurementMethod), "measurementMethod"] <- "Assumed microscopic" # *** COMMENT: ASSUMPTION made for the visualization below

# Prepare mask
msk <- subset (env.stack [[1]], "Mask")

# Split by taxa
dat.spl <- split(tot, f = tot$scientificName)
df.taxa.obs <- data.frame("taxon" = names(unlist(lapply(dat.spl, nrow))), "obs" = as.numeric(unlist(lapply(dat.spl, nrow))) )
df.taxa.obs <- df.taxa.obs[order(df.taxa.obs$taxon),]
df.tax.obs <- df.taxa.obs[which(df.taxa.obs$obs >= 9),]
rownames(df.tax.obs) <- 1:nrow(df.tax.obs)
df.tax.obs # 35 taxa with at least 9 gridded monthly-presences

### ===================================================================================
### 1. Visualize data of total diazotrophs geographically and ecologically
### ===================================================================================

cea <- 1.2 # axis lables type size
sz <- 0.3 # axes in plot
col.pts <- "black" # color of points
point <- 16 # point symbol
col.line <- "grey65" # color of help-lines
sz <- 0.3 # axes numbers in plot
nr.size <- 0.3 # axes numbers in lat grad plot
s.col <- "white"
col.pts <- "black"
pt.sym <- 16
lwd.box <- 15
dist <- 0.5
cea <- 1 # lables
lab.size <- 0.98
tck.dist <- -0.006 # originally -0.007
dist <- 0.22 # distance of numbers (axis lables) from axis
pt.size <- 0.8 # size of points
lwd.size <- 12 # width of lines of points
pt.sym <- 1
pt.col <- "black"
lwd.box <- 6 # Lines
c.col <- "#D6EFED" # Continents
c.col <- "black" # Continents
col.sh <- colorRampPalette(c("#B7C8DE",  colorRampPalette(c("black", "royalblue"))(5)[2])                        ) (20)[3]


### PLOT ###

png (filename = paste0("./Righetti_et_al_notion/1_Visualize_data/Output_2DE/1 Total (n=",nrow( tot[!duplicated(tot[ , c("scientificName", "x", "y", "year", "month", "day", "depth", "measurementMethod")]), ] ),").png"), width = 6100, height = 6100, units = "px", pointsize = 80) # original width 5600, original height 3600

# -------------------------geography: frame one-------------------------------------
par(mfrow=c(2,1),mar=c(0.7,0.7,1.05,0.7),oma=c(2.5,2.5,2,3.2),cex=0.9)
# Continents
plot(tot$x, tot$y,  pch = pt.sym, xlab="",  ylab="",  ylim = c(-77,84), xlim = c (-167, 167), col = pt.col, axes = F,  xaxt = "n", yaxt = "n", cex = pt.size, type ="n")
axis (1, at = seq(-180, 180 ,30), labels = F, cex.axis = 0.4, tck = 1, lwd = 0.8, lty = 1, mgp = c(1, dist, 0), col = "grey30") # y-axis grid-lines
axis (2, at=seq(-90,90,15), labels = F, cex.axis = 0.4, tck = 1, lwd = 0.8, lty = 1, mgp = c(1, dist, 0), col = "grey30") # x-axis grid-lines
plot(msk, add =T, col = c("grey76"), legend=F)
axis (1, at = seq(-180, 180 ,30), labels = F, cex.axis = 0.4, tck = 1, lwd = 1.2, lty = 1, mgp = c(1, dist, 0), col = "grey40") # y-axis grid-lines
axis (2, at=seq(-90,90,15), labels = F, cex.axis = 0.4, tck = 1, lwd = 1.2, lty = 1, mgp = c(1, dist, 0), col = "grey40") # x-axis grid-lines
map ("world", interior=T, fill = T, lwd = 0.00001, boundary=TRUE, col="snow1", col.fill = "snow1",add = T)
#
## Visualization of data: microscopic records = circle (pch 1), nifH records = filled triangle (pch 25), OTU records = cube (pch 22)
## Luo
dt <- tot[which(tot$set=="Luo"),] # unique(dt$measurementMethod)
points(
jitter(x = dt[which(dt$measurementMethod%in%c("Standard light microscopy", "Epifluorescence microscopy", "Standard light or epifluorescence microscopy", "Direct biomass analysis")),]$x, factor = 1, amount = 1.2),
jitter(x = dt[which(dt$measurementMethod%in%c("Standard light microscopy", "Epifluorescence microscopy", "Standard light or epifluorescence microscopy", "Direct biomass analysis")),]$y, factor = 1, amount = 0.6),
col = "orange", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 1)
points(
jitter(x = dt[which(dt$measurementMethod%in%c("qPCR_nifH_detection")),]$x, factor = 1, amount = 1.2),
jitter(x = dt[which(dt$measurementMethod%in%c("qPCR_nifH_detection")),]$y, factor = 1, amount = 0.6),
col = "orange", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 25)
#
## Tan
dt <- tot[which(tot$set=="Tan"),] # unique(dt$measurementMethod)
points(
jitter(x = dt[which(dt$measurementMethod%in%c("qPCR_nifH_detection")),]$x, factor = 1, amount = 1.2),
jitter(x = dt[which(dt$measurementMethod%in%c("qPCR_nifH_detection")),]$y, factor = 1, amount = 0.6),
col = "darkblue", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 25)
#
## Phy
dt <-  tot[which(tot$set=="Phy"),]  # unique(dt$measurementMethod)
points(
jitter(x = dt$x, factor = 1, amount = 1.2),
jitter(x = dt$y, factor = 1, amount = 0.6),
col = "yellow2", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 1) # Assuming microscopy basaed (information on method, NA)
#
## Gra
dt <- tot[which(tot$set=="Gra"),] # unique(dt$measurementMethod)
points(
jitter(x = dt[which(dt$measurementMethod%in%c("qPCR_nifH_detection")),]$x, factor = 1, amount = 1.2),
jitter (x = dt[which(dt$measurementMethod%in%c("qPCR_nifH_detection")),]$y, factor = 1, amount = 0.6),
col = "firebrick2", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 25)
#
## Gbi
dt <-  tot[which(tot$set=="Gbi"),] # unique(dt$measurementMethod)
points(
jitter(x = dt$x, factor = 1, amount = 1.2),
jitter(x = dt$y, factor = 1, amount = 0.6),
col = "seagreen2", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 1) # Assuming microscopy basaed (information on method, NA)
#
## Obi
dt <-  tot[which(tot$set=="Obi"),] # unique(dt$measurementMethod)
points(
jitter(x = dt$x, factor = 1, amount = 1.2),
jitter(x = dt$y, factor = 1, amount = 0.6),
col = "pink", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 1) # Assuming microscopy basaed (information on method, NA)
#
## Tara OTU based
dt <-  tot[which(tot$set=="TarOTU"),] # unique(dt$measurementMethod)
points(
jitter(x = dt[which(dt$measurementMethod%in%c("SSU_rRNA_detection")),]$x, factor = 1, amount = 1.2),
jitter (x = dt[which(dt$measurementMethod%in%c("SSU_rRNA_detection")),]$y, factor = 1, amount = 0.6),
col = "slateblue1", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 22)
#
## Martinez-Perez (M96) based
dt <-  tot[which(tot$set=="Mar"),] # unique(dt$measurementMethod)
points(
jitter(x = dt[which(dt$measurementMethod%in%c("Epifluorescence microscopy")),]$x, factor = 1, amount = 1.2),
jitter(x = dt[which(dt$measurementMethod%in%c("Epifluorescence microscopy")),]$y, factor = 1, amount = 0.6),
col = "firebrick3", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 1)
points(
jitter(x = dt[which(dt$measurementMethod%in%c("qPCR_nifH_detection")),]$x, factor = 1, amount = 1.2),
jitter(x = dt[which(dt$measurementMethod%in%c("qPCR_nifH_detection")),]$y, factor = 1, amount = 0.6),
col = "firebrick3", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 25)
#
## Legend
legend(-160, 91.7, c("Luo", "Tang", "PhytoBase", "Gradoville", "Gbif", "Obis", "Tara", "M96"), lty=c(1, 1, 1, 1, 1, 1, 1, 1, 1), col=c("orange", "darkblue","yellow2", "firebrick2", "seagreen2", "pink", "slateblue1", "firebrick3" ), bty = "n", lwd=32, cex=lab.size, horiz =T, text.width=22)
tot.clsets <- tot[!duplicated(tot[ , c("scientificName", "x", "y", "year", "month", "day", "depth", "measurementMethod")]), ] # set information is removed (i.e., the origin of the records), divergent methods are retained
legend(-160, - 72, c(
paste0("microscopic or NA (",nrow(tot.clsets[which(tot.clsets $measurementMethod%in%c("Standard light microscopy", "Epifluorescence microscopy", "Standard light or epifluorescence microscopy", "Direct biomass analysis", "Assumed microscopic")),]),")"),
paste0("nifH records (",nrow(tot.clsets[which(tot.clsets $measurementMethod%in%c("qPCR_nifH_detection")),]),")"),
paste0("SSU rRNA records (",nrow(tot.clsets[which(tot.clsets$measurementMethod%in%c("SSU_rRNA_detection")),]),
")       Total: ",nrow(tot.clsets), "    At monthly 1deg.: ", nrow(tot[!duplicated(tot[ , c("scientificName", "x", "y", "month")]), ] )    )
), lty=c(NA, NA, NA), pch = c(1, 25, 22), col=c(1,1,1), bty = "n", lwd=lwd.size+0.7, cex=lab.size, horiz =T, text.width=58)
axis (1,at=seq(-180,180,30), cex.axis = lab.size, lwd= lwd.box,  tck = tck.dist, mgp = c(1,dist, 0), las = 1)
axis (2, at = c(-75, -60, -45, -30, -15, 0, 15, 30, 45, 60, 75, 90), cex.axis = lab.size, lwd= lwd.box,  tck = tck.dist, mgp = c(1,dist, 0), las = 1); box(lwd = lwd.box)
mtext(expression("Longitude ("*degree*")"), side = 1,  line = 1.2, cex = lab.size)
mtext(expression("Latitude ("*degree*")"), side=2, line=1.2, cex=lab.size)
title(main="The realized superniche", cex = lab.size-0.2, line = 0.3)


# -------------------------environment temp log(mld) space: frame two-------------------------------------
	# Prepare data temperature and mld
	var1 <- vector (  )
	var2 <- vector (  )
	for (i in 1:length(env.stack)) {
	param1 <- subset (env.stack[[i]], "MLD1")
	param1[msk==0] <- NA
	var1 <- append (var1, values (param1))
	param2 <- subset (env.stack [[i]], "T")
	param2[msk==0] <- NA
	var2 <- append (var2, values (param2))}
	#
	# Prepare shading environmental background
	var1.new <- var1 [-which (is.na(var1)  |  is.na(var2)) ]
	var2.new <- var2 [-which (is.na(var1)  |  is.na(var2)) ]
	env.space <- SpatialPoints (coords = cbind(log10(var1.new), var2.new))
	ranges <- raster(xmn = round (min(log10(var1.new)), digits = 3) ,
	xmx = round (max(log10(var1.new)), digits = 3),
	ymn = round (min(var2.new), digits = 2),
	ymx = round (max(var2.new), digits = 2),
	res = c(0.05, 0.5))
	env.space.new <- rasterize ( x = env.space , y= ranges, fun = "count" )
	ranges2 <- focal( env.space.new , w=matrix(c(1,1,1,1,1,1,2,2,2,1,1,2,3,2,1,1,2,2,2,1,1,1,1,1,1),5,5), fun=mean)
# Plot environmental background (as kernel)
par(fig=c(0 , 0.49 , 0 , 0.49), new=TRUE)
smoothScatter(log10(var1),var2,colramp=colorRampPalette(c("white","grey30"),50), bandwidth=c(0.015,0.25), nrpoints=0, nbin=400, axes = F, xlab = "log (MLD) [m]", ylab = "Temperature", xlim = c(1,3), ylim = c(-2,34))
# contour (ranges2, xlim = c(1,400), levels = c(5, 20, 50, 100, 200, 300, 400, 600, 800, 1000),  xlab = "MLD.deBoyerMontegut", ylab = "Temperature", drawlabels = F, add = T, lwd = 4) # does not work
# Plot ecological super niche
dat <- tot
df <- data.frame(x = log10(dat$MLD1), y = dat$T, measurementMethod = dat$measurementMethod, set = dat$set, scientificName = dat$scientificName)
	df[ ,1][is.infinite(df[,1])] = NA # remove -Inf
	df[ ,2][is.infinite(df[,2])] = NA # remove -Inf
	df <- df[complete.cases(df),] # 17343 x 5
	#
	#
	## Luo
	dt <- df[which(df$set=="Luo"),] # unique(dt$measurementMethod)
	da <- dt[which(dt$measurementMethod%in%c("Standard light microscopy", "Epifluorescence microscopy", "Standard light or epifluorescence microscopy", "Direct biomass analysis")),]
	x <- jitter ( x = da$x, factor = 1, amount = 0.01)
	y <- jitter ( x = da$y, factor = 1, amount = 0.2)
	points(x, y, cex = 0.7, col = "orange", pch = 1, xlab = "log(MLD)", ylab = "Temperature", xlim = c(1,2.5), ylim = c(-2, 34), lwd = lwd.size+0.7)
	da <- dt[which(dt$measurementMethod%in%c("qPCR_nifH_detection")),]
	x <- jitter ( x = da$x, factor = 1, amount = 0.005)
	y <- jitter ( x = da$y, factor = 1, amount = 0.12)
	points(x, y, cex = 0.7, col = "orange", pch = 25, xlab = "log(MLD)", ylab = "Temperature", xlim = c(1,2.5), ylim = c(-2, 34), lwd = lwd.size+0.7)
	#
	## Tan
	dt <- df[which(df$set=="Tan"),] # unique(dt$measurementMethod)
	da <- dt[which(dt$measurementMethod%in%c("qPCR_nifH_detection")),]
	x <- jitter ( x = da$x, factor = 1, amount = 0.01)
	y <- jitter ( x = da$y, factor = 1, amount = 0.2)
	points(x, y, cex = 0.7, col = "darkblue", pch = 25, xlab = "log(MLD)", ylab = "Temperature", xlim = c(1,2.5), ylim = c(-2, 34), lwd = lwd.size+0.7)
	#
	## Phy
	dt <-  df[which(df$set=="Phy"),]  # unique(dt$measurementMethod)
	x <- jitter ( x = dt$x, factor = 1, amount = 0.01)
	y <- jitter ( x = dt$y, factor = 1, amount = 0.2)
	points(dt$x, dt$y, col = "yellow2", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 1) # Assuming microscopy basaed (information on method, NA)
	#
	## Gra
	dt <- df[which(df$set=="Gra"),] # unique(dt$measurementMethod)
	x <- jitter ( x = dt$x, factor = 1, amount = 0.01)
	y <- jitter ( x = dt$y, factor = 1, amount = 0.2)
	points(dt[which(dt$measurementMethod%in%c("qPCR_nifH_detection")),]$x, dt[which(dt$measurementMethod%in%c("qPCR_nifH_detection")),]$y, col = "firebrick2", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 25)
	#
	## Gbi
	dt <-  df[which(df$set=="Gbi"),] # unique(dt$measurementMethod)
	x <- jitter ( x = dt$x, factor = 1, amount = 0.01)
	y <- jitter ( x = dt$y, factor = 1, amount = 0.2)
	points(dt$x, dt$y, col = "seagreen2", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 1) # Assuming microscopy basaed (information on method, NA)
	#
	## Obi
	dt <-  df[which(df$set=="Obi"),] # unique(dt$measurementMethod)
	x <- jitter ( x = dt$x, factor = 1, amount = 0.01)
	y <- jitter ( x = dt$y, factor = 1, amount = 0.2)
	points(dt$x, dt$y, col = "pink", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 1) # Assuming microscopy basaed (information on method, NA)
	#
	## Tara OTU based
	dt <-  df[which(df$set=="TarOTU"),] # unique(dt$measurementMethod)
	x <- jitter ( x = dt$x, factor = 1, amount = 0.01)
	y <- jitter ( x = dt$y, factor = 1, amount = 0.2)
	points(dt[which(dt$measurementMethod%in%c("SSU_rRNA_detection")),]$x, dt[which(dt$measurementMethod%in%c("SSU_rRNA_detection")),]$y, col = "slateblue1", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 22)
	#
	## M96 (Mar)
	dt <- df[which(df$set=="Mar"),] # unique(dt$measurementMethod)
	da <- dt[which(dt$measurementMethod%in%c("Epifluorescence microscopy")),]
	x <- jitter ( x = da$x, factor = 1, amount = 0.01)
	y <- jitter ( x = da$y, factor = 1, amount = 0.2)
	points(x, y, cex = 0.7, col = "firebrick3", pch = 1, xlab = "log(MLD)", ylab = "Temperature", xlim = c(1,2.5), ylim = c(-2, 34), lwd = lwd.size+0.7)
	da <- dt[which(dt$measurementMethod%in%c("qPCR_nifH_detection")),]
	x <- jitter ( x = da$x, factor = 1, amount = 0.005)
	y <- jitter ( x = da$y, factor = 1, amount = 0.12)
	points(x, y, cex = 0.7, col = "firebrick3", pch = 25, xlab = "log(MLD)", ylab = "Temperature", xlim = c(1,2.5), ylim = c(-2, 34), lwd = lwd.size+0.7)
	#
	## Grand total -----------------------------------------
	dat <- tot.clsets # set information is removed (i.e., the origin of the records), divergent methods are retained
	df <- data.frame(x = log10(dat$MLD1), y = dat$T, measurementMethod = dat$measurementMethod, set = dat$set, scientificName = dat$scientificName)
		df[ ,1][is.infinite(df[,1])] = NA # remove -Inf
		df[ ,2][is.infinite(df[,2])] = NA # remove -Inf
		df <- df[complete.cases(df),] # 14685 x 5
	x <- jitter ( x = df$x, factor = 1, amount = 0.0)
	y <- jitter ( x = df$y, factor = 1, amount = 0.0)
	# points(x, y, cex = 0.7, col = "grey46", pch = 19, xlab = "log(MLD)", ylab = "Temperature", xlim = c(1,2.5), ylim = c(-2, 34))
	segments (  x0=mean(df$x)-sd(df$x), x1=mean(df$x)+sd(df$x) , y0=mean(df$y), y1=mean(df$y), lwd = 34, col = "grey20" )
	segments (  y0=mean(df$y)-sd(df$y), y1=mean(df$y)+sd(df$y) , x0=mean(df$x), x1=mean(df$x), lwd = 34, col = "grey20" )
	legend("topright", c("Full ocean (frequency)"), lty=c(NA), pch=c(22), col=c("grey80" ), bty = "n", lwd=42, cex= lab.size, horiz =T)
	points(pam(df, 1, metric = "euclidean")$medoids, cex = 2.5, col = "grey20", pch = 1, lwd = 48)
	axis (1,at=seq(0, 3, 0.5), cex.axis = lab.size, lwd= lwd.box,  tck = tck.dist, mgp = c(1,dist, 0), las = 1)
	axis (2, at=seq(-5, 35, 5), cex.axis = lab.size, lwd= lwd.box,  tck = tck.dist, mgp = c(1,dist, 0), las = 1); box(lwd = lwd.box)
	mtext(expression("log(mixed-layer depth [m])"), side = 1,  line = 1.4, cex = lab.size)
	mtext(expression("Temperature ("*degree*"C)"), side=2, line=1.2, cex=lab.size)


# -------------------------environment log(no3-) inverse(mlpar) space: frame two-------------------------------------
	# Prepare data temperature and mld
	var1 <- vector (  )
	var2 <- vector (  )
	for (i in 1:length(env.stack)) {
	param1 <- subset (env.stack [[i]], "MLPAR1")
	param1[msk==0] <- NA
	var1 <- append (var1, values (param1))
	param2 <- subset (env.stack [[i]], "N")
	param2[msk==0] <- NA
	var2 <- append (var2, values (param2))
	}
	# Prepare shading environmental background
	var1.new <- var1 [-which (is.na(var1)  |  is.na(var2)) ]
	var2.new <- var2 [-which (is.na(var1)  |  is.na(var2)) ]
	env.space <- SpatialPoints (coords = cbind(log10(var1.new), var2.new))
	ranges <- raster(xmn = round (min(log10(var1.new)), digits = 3) ,
	xmx = round (max(log10(var1.new)), digits = 3),
	ymn = round (min(var2.new), digits = 2),
	ymx = round (max(var2.new), digits = 2),
	res = c(0.05, 0.5))
	env.space.new <- rasterize ( x = env.space , y= ranges, fun = "count" )
	ranges2 <- focal( env.space.new , w=matrix(c(1,1,1,1,1,1,2,2,2,1,1,2,3,2,1,1,2,2,2,1,1,1,1,1,1),5,5), fun=mean)
# Plot environmental background (as kernel)
par(fig=c(0.51 , 1 , 0 , 0.49), new=TRUE)
smoothScatter(var1,log10(var2),colramp=colorRampPalette(c("white","grey30"),50),bandwidth=c(0.25,0.06),nrpoints=0,nbin=300, axes=F, xlab = expression("MLPAR [E"*~m^-3 * d^-1~"]"), ylab = expression("log (NO3-) ["*~mu*mol/L~"]"), ylim = c(-4.5,1.6), xlim = c(50, 0))
# Plot ecological super niche
dat <- tot
df <- data.frame(x = dat$MLPAR1, y = log10(dat$N), measurementMethod = dat$measurementMethod, set = dat$set, scientificName = dat$scientificName)
	df[ ,1][is.infinite(df[,1])] = NA # remove -Inf
	df[ ,2][is.infinite(df[,2])] = NA # remove -Inf
	df <- df[complete.cases(df),] # 16701 x 5
	## Luo -----------------------------------------
	dt <- df[which(df$set=="Luo"),] # unique(dt$measurementMethod)
	da <- dt[which(dt$measurementMethod%in%c("Standard light microscopy", "Epifluorescence microscopy", "Standard light or epifluorescence microscopy", "Direct biomass analysis")),]
	x <- jitter ( x = da$x, factor = 1, amount = 0.6)
	y <- jitter ( x = da$y, factor = 1, amount = 0.06)
	points(x, y, cex = 0.7, col = "orange", pch = 1, xlab = "log(MLD)", ylab = "Temperature", xlim = c(1,2.5), ylim = c(-2, 34), lwd = lwd.size+0.7)
	da <- dt[which(dt$measurementMethod%in%c("qPCR_nifH_detection")),]
	x <- jitter ( x = da$x, factor = 1, amount = 0.6)
	y <- jitter ( x = da$y, factor = 1, amount = 0.06)
	points(x, y, cex = 0.7, col = "orange", pch = 25, xlab = "log(MLD)", ylab = "Temperature", xlim = c(1,2.5), ylim = c(-2, 34), lwd = lwd.size+0.7)
	#
	## Tan
	dt <- df[which(df$set=="Tan"),] # unique(dt$measurementMethod)
	da <- dt[which(dt$measurementMethod%in%c("qPCR_nifH_detection")),]
	x <- jitter ( x = da$x, factor = 1, amount = 0.6)
	y <- jitter ( x = da$y, factor = 1, amount = 0.06)
	points(x, y, cex = 0.7, col = "darkblue", pch = 25, xlab = "log(MLD)", ylab = "Temperature", xlim = c(1,2.5), ylim = c(-2, 34), lwd = lwd.size+0.7)
	#
	## Phy
	dt <-  df[which(df$set=="Phy"),]  # unique(dt$measurementMethod)
	x <- jitter ( x = dt$x, factor = 1, amount = 0.01)
	y <- jitter ( x = dt$y, factor = 1, amount = 0.2)
	points(dt$x, dt$y, col = "yellow2", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 1) # Assuming microscopy basaed (information on method, NA)
	#
	## Gra
	dt <- df[which(df$set=="Gra"),] # unique(dt$measurementMethod)
	x <- jitter ( x = dt$x, factor = 1, amount = 0.01)
	y <- jitter ( x = dt$y, factor = 1, amount = 0.2)
	points(dt[which(dt$measurementMethod%in%c("qPCR_nifH_detection")),]$x, dt[which(dt$measurementMethod%in%c("qPCR_nifH_detection")),]$y, col = "firebrick2", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 25)
	#
	## Gbi
	dt <-  df[which(df$set=="Gbi"),] # unique(dt$measurementMethod)
	x <- jitter ( x = dt$x, factor = 1, amount = 0.01)
	y <- jitter ( x = dt$y, factor = 1, amount = 0.2)
	points(dt$x, dt$y, col = "seagreen2", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 1) # Assuming microscopy basaed (information on method, NA)
	#
	## Obi
	dt <-  df[which(df$set=="Obi"),] # unique(dt$measurementMethod)
	x <- jitter ( x = dt$x, factor = 1, amount = 0.01)
	y <- jitter ( x = dt$y, factor = 1, amount = 0.2)
	points(dt$x, dt$y, col = "pink", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 1) # Assuming microscopy basaed (information on method, NA)
	#
	## Tara OTU based
	dt <-  df[which(df$set=="TarOTU"),] # unique(dt$measurementMethod)
	x <- jitter ( x = dt$x, factor = 1, amount = 0.01)
	y <- jitter ( x = dt$y, factor = 1, amount = 0.2)
	points(dt[which(dt$measurementMethod%in%c("SSU_rRNA_detection")),]$x, dt[which(dt$measurementMethod%in%c("SSU_rRNA_detection")),]$y, col = "slateblue1", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 22)
	#
	## M96 (Mar)
	dt <- df[which(df$set=="Mar"),] # unique(dt$measurementMethod)
	da <- dt[which(dt$measurementMethod%in%c("Epifluorescence microscopy")),]
	x <- jitter ( x = da$x, factor = 1, amount = 0.01)
	y <- jitter ( x = da$y, factor = 1, amount = 0.2)
	points(x, y, cex = 0.7, col = "firebrick3", pch = 1, xlab = "log(MLD)", ylab = "Temperature", xlim = c(1,2.5), ylim = c(-2, 34), lwd = lwd.size+0.7)
	da <- dt[which(dt$measurementMethod%in%c("qPCR_nifH_detection")),]
	x <- jitter ( x = da$x, factor = 1, amount = 0.005)
	y <- jitter ( x = da$y, factor = 1, amount = 0.12)
	points(x, y, cex = 0.7, col = "firebrick3", pch = 25, xlab = "log(MLD)", ylab = "Temperature", xlim = c(1,2.5), ylim = c(-2, 34), lwd = lwd.size+0.7)
	#
	## Grand total -----------------------------------------
	dat <- tot.clsets # set information is removed (i.e., the origin of the records), divergent methods are retained
	df <- data.frame(x = dat$MLPAR1, y = log10(dat$N), measurementMethod = dat$measurementMethod, set = dat$set, scientificName = dat$scientificName)
		df[ ,1][is.infinite(df[,1])] = NA # remove -Inf
		df[ ,2][is.infinite(df[,2])] = NA # remove -Inf
		df <- df[complete.cases(df),] # 14685 x 5
	x <- jitter ( x = df$x, factor = 1, amount = 0.0)
	y <- jitter ( x = df$y, factor = 1, amount = 0.0)
	# points (x, y, cex = 0.7, col = "grey46", pch = 19,  xlab = expression("MLPAR [E"*~m^-3 * d^-1~"]"), ylab = expression("log (NO3-) ["*~mu*mol/L~"]"), ylim = c(-4.5,1.6), xlim = c(50, 0))  # original points
	segments (  x0=mean(df$x)-sd(df$x), x1=mean(df$x)+sd(df$x) , y0=mean(df$y), y1=mean(df$y), lwd = 34, col = "grey20" )
	segments (  y0=mean(df$y)-sd(df$y), y1=mean(df$y)+sd(df$y) , x0=mean(df$x), x1=mean(df$x), lwd = 34, col = "grey20" )
	# legend("topright", legend= "Realized superniche", bty="n", text.font = 2, cex = lab.size)
	legend("topleft", c("Full ocean (frequency)"), lty=c(NA), pch=c(22), col=c("grey80" ), bty = "n", lwd=42, cex= lab.size, horiz =T)
	points(pam(df, 1, metric = "euclidean")$medoids, cex = 2.5, col = "grey20", pch = 1, lwd = 48)
axis (1,at=seq(50, 0, -10), cex.axis = lab.size, lwd= lwd.box,  tck = tck.dist, mgp = c(1,dist, 0), las = 1)
axis (2, at=seq(-4.3, 1.6, 1), cex.axis = lab.size, lwd= lwd.box,  tck = tck.dist, mgp = c(1,dist, 0), las = 1); box(lwd = lwd.box)
mtext(expression("MLPAR (E"*~m^-3 * d^-1~")"), side = 1,  line = 1.4, cex = lab.size)
mtext(expression("log (NO3-) ["*~mu*mol/L~"]"), side=2, line=1.4, cex=lab.size)
#
#
#
dev.off()



### ===================================================================================
### 2. Visualize data of diazotroph taxa geographically and ecologically
### ===================================================================================


### B. Plot genera, species or strains one by one ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
for ( i in 1:length(df.tax.obs$taxon) ){

png (filename = paste0("./1_Visualize_data/Output_2/",i,"_",df.tax.obs$taxon[i]," (n=",nrow(tot[which(tot$taxon==df.tax.obs$taxon[i]),][!duplicated(tot[which(tot$taxon==df.tax.obs$taxon[i]),][ , c("scientificName", "x", "y", "year", "month", "day", "depth", "measurementMethod")]), ]),").png"), width = 6100, height = 6100, units = "px", pointsize = 80)

# -------------------------geography: frame one-------------------------------------
q <- i
par(mfrow=c(2,1),mar=c(0.7,0.7,1.05,0.7),oma=c(2.5,2.5,2,3.2),cex=0.9)
# Select data of taxon
dat.tax <- tot[which(tot$taxon==df.tax.obs$taxon[i]),]
# Continents
plot(dat.tax$x, dat.tax$y,  pch = pt.sym, xlab="",  ylab="",  ylim = c(-77,84), xlim = c (-167, 167), col = pt.col, axes = F,  xaxt = "n", yaxt = "n", cex = pt.size, type ="n")
axis (1, at = seq(-180, 180 ,30), labels = F, cex.axis = 0.4, tck = 1, lwd = 0.8, lty = 1, mgp = c(1, dist, 0), col = "grey30") # y-axis grid-lines
axis (2, at=seq(-90,90,15), labels = F, cex.axis = 0.4, tck = 1, lwd = 0.8, lty = 1, mgp = c(1, dist, 0), col = "grey30") # x-axis grid-lines
plot(msk, add =T, col = c("grey76"), legend=F)
axis (1, at = seq(-180, 180 ,30), labels = F, cex.axis = 0.4, tck = 1, lwd = 1.2, lty = 1, mgp = c(1, dist, 0), col = "grey40") # y-axis grid-lines
axis (2, at=seq(-90,90,15), labels = F, cex.axis = 0.4, tck = 1, lwd = 1.2, lty = 1, mgp = c(1, dist, 0), col = "grey40") # x-axis grid-lines
map ("world", interior=T, fill = T, lwd = 0.00001, boundary=TRUE, col="snow1", col.fill = "snow1",add = T)
#
## Visualization of data: microscopic records = circle (pch 1), nifH records = filled triangle (pch 25), OTU records = cube (pch 22)
## Luo
dt <- dat.tax[which(dat.tax$set=="Luo"),] # unique(dt$measurementMethod)
points(
jitter(x = dt[which(dt$measurementMethod%in%c("Standard light microscopy", "Epifluorescence microscopy", "Standard light or epifluorescence microscopy", "Direct biomass analysis")),]$x, factor = 1, amount = 1.2),
jitter(x = dt[which(dt$measurementMethod%in%c("Standard light microscopy", "Epifluorescence microscopy", "Standard light or epifluorescence microscopy", "Direct biomass analysis")),]$y, factor = 1, amount = 0.6),
col = "orange", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 1)
points(
jitter(x = dt[which(dt$measurementMethod%in%c("qPCR_nifH_detection")),]$x, factor = 1, amount = 1.2),
jitter(x = dt[which(dt$measurementMethod%in%c("qPCR_nifH_detection")),]$y, factor = 1, amount = 0.6),
col = "orange", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 25)
#
## Tan
dt <- dat.tax[which(dat.tax$set=="Tan"),] # unique(dt$measurementMethod)
points(
jitter(x = dt[which(dt$measurementMethod%in%c("qPCR_nifH_detection")),]$x, factor = 1, amount = 1.2),
jitter(x = dt[which(dt$measurementMethod%in%c("qPCR_nifH_detection")),]$y, factor = 1, amount = 0.6),
col = "darkblue", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 25)
#
## Phy
dt <-  dat.tax[which(dat.tax$set=="Phy"),]  # unique(dt$measurementMethod)
points(
jitter(x = dt$x, factor = 1, amount = 1.2),
jitter(x = dt$y, factor = 1, amount = 0.6),
col = "yellow2", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 1) # Assuming microscopy basaed (information on method, NA)
#
## Gra
dt <- dat.tax[which(dat.tax$set=="Gra"),] # unique(dt$measurementMethod)
points(
jitter(x = dt[which(dt$measurementMethod%in%c("qPCR_nifH_detection")),]$x, factor = 1, amount = 1.2),
jitter (x = dt[which(dt$measurementMethod%in%c("qPCR_nifH_detection")),]$y, factor = 1, amount = 0.6),
col = "firebrick2", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 25)
#
## Gbi
dt <-  dat.tax[which(dat.tax$set=="Gbi"),] # unique(dt$measurementMethod)
points(
jitter(x = dt$x, factor = 1, amount = 1.2),
jitter(x = dt$y, factor = 1, amount = 0.6),
col = "seagreen2", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 1) # Assuming microscopy basaed (information on method, NA)
#
## Obi
dt <-  dat.tax[which(dat.tax$set=="Obi"),] # unique(dt$measurementMethod)
points(
jitter(x = dt$x, factor = 1, amount = 1.2),
jitter(x = dt$y, factor = 1, amount = 0.6),
col = "pink", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 1) # Assuming microscopy basaed (information on method, NA)
#
## Tara OTU based
dt <-  dat.tax[which(dat.tax$set=="TarOTU"),] # unique(dt$measurementMethod)
points(
jitter(x = dt[which(dt$measurementMethod%in%c("SSU_rRNA_detection")),]$x, factor = 1, amount = 1.2),
jitter (x = dt[which(dt$measurementMethod%in%c("SSU_rRNA_detection")),]$y, factor = 1, amount = 0.6),
col = "slateblue1", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 22)
#
## Martinez-Perez (M96) based
dt <-  dat.tax[which(dat.tax$set=="Mar"),] # unique(dt$measurementMethod)
points(
jitter(x = dt[which(dt$measurementMethod%in%c("Epifluorescence microscopy")),]$x, factor = 1, amount = 1.2),
jitter(x = dt[which(dt$measurementMethod%in%c("Epifluorescence microscopy")),]$y, factor = 1, amount = 0.6),
col = "firebrick3", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 1)
points(
jitter(x = dt[which(dt$measurementMethod%in%c("qPCR_nifH_detection")),]$x, factor = 1, amount = 1.2),
jitter(x = dt[which(dt$measurementMethod%in%c("qPCR_nifH_detection")),]$y, factor = 1, amount = 0.6),
col = "firebrick3", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 25)
#
## Legend
legend(-160, 91.7, c("Luo", "Tang", "PhytoBase", "Gradoville", "Gbif", "Obis", "Tara", "M96"), lty=c(1, 1, 1, 1, 1, 1, 1, 1, 1), col=c("orange", "darkblue","yellow2", "firebrick2", "seagreen2", "pink", "slateblue1", "firebrick3" ), bty = "n", lwd=32, cex=lab.size, horiz =T, text.width=22)
tot.clsets <- dat.tax[!duplicated(dat.tax[ , c("scientificName", "x", "y", "year", "month", "day", "depth", "measurementMethod")]), ] # set information is removed (i.e., the origin of the records), divergent methods are retained
legend(-160, - 72, c(
paste0("microscopic or NA (",nrow(tot.clsets[which(tot.clsets $measurementMethod%in%c("Standard light microscopy", "Epifluorescence microscopy", "Standard light or epifluorescence microscopy", "Direct biomass analysis", "Assumed microscopic")),]),")"),
paste0("nifH records (",nrow(tot.clsets[which(tot.clsets $measurementMethod%in%c("qPCR_nifH_detection")),]),")"),
paste0("SSU rRNA records (",nrow(tot.clsets[which(tot.clsets$measurementMethod%in%c("SSU_rRNA_detection")),]),
")         Total: ",nrow(tot.clsets), "        At monthly 1deg.: ", nrow(dat.tax[!duplicated(dat.tax[ , c("scientificName", "x", "y", "month")]), ] )    )
), lty=c(NA, NA, NA), pch = c(1, 25, 22), col=c(1,1,1), bty = "n", lwd=lwd.size+0.7, cex=lab.size, horiz =T, text.width=58)
axis (1,at=seq(-180,180,30), cex.axis = lab.size, lwd= lwd.box,  tck = tck.dist, mgp = c(1,dist, 0), las = 1)
axis (2, at = c(-75, -60, -45, -30, -15, 0, 15, 30, 45, 60, 75, 90), cex.axis = lab.size, lwd= lwd.box,  tck = tck.dist, mgp = c(1,dist, 0), las = 1); box(lwd = lwd.box)
mtext(expression("Longitude ("*degree*")"), side = 1,  line = 1.2, cex = lab.size)
mtext(expression("Latitude ("*degree*")"), side=2, line=1.2, cex=lab.size)
title(main=paste0("Realized niche of ",df.tax.obs$taxon[i]), cex = lab.size-0.2, line = 0.3)

# -------------------------environment temp log(mld) space: frame two-------------------------------------
	# Prepare data temperature and mld
	var1 <- vector (  )
	var2 <- vector (  )
	for (i in 1:length(env.stack)) {
	param1 <- subset (env.stack[[i]], "MLD1")
	param1[msk==0] <- NA
	var1 <- append (var1, values (param1))
	param2 <- subset (env.stack [[i]], "T")
	param2[msk==0] <- NA
	var2 <- append (var2, values (param2))}
	#
	# Prepare shading environmental background
	var1.new <- var1 [-which (is.na(var1)  |  is.na(var2)) ]
	var2.new <- var2 [-which (is.na(var1)  |  is.na(var2)) ]
	env.space <- SpatialPoints (coords = cbind(log10(var1.new), var2.new))
	ranges <- raster(xmn = round (min(log10(var1.new)), digits = 3) ,
	xmx = round (max(log10(var1.new)), digits = 3),
	ymn = round (min(var2.new), digits = 2),
	ymx = round (max(var2.new), digits = 2),
	res = c(0.05, 0.5))
	env.space.new <- rasterize ( x = env.space , y= ranges, fun = "count" )
	ranges2 <- focal( env.space.new , w=matrix(c(1,1,1,1,1,1,2,2,2,1,1,2,3,2,1,1,2,2,2,1,1,1,1,1,1),5,5), fun=mean)
# Plot environmental background (as kernel)
par(fig=c(0 , 0.49 , 0 , 0.49), new=TRUE)
smoothScatter(log10(var1),var2,colramp=colorRampPalette(c("white","grey30"),50), bandwidth=c(0.015,0.25), nrpoints=0, nbin=400, axes = F, xlab = "log (MLD) [m]", ylab = "Temperature", xlim = c(1,3), ylim = c(-2,34))
# contour (ranges2, xlim = c(1,400), levels = c(5, 20, 50, 100, 200, 300, 400, 600, 800, 1000),  xlab = "MLD.deBoyerMontegut", ylab = "Temperature", drawlabels = F, add = T, lwd = 4) # does not work
# Plot ecological super niche
dat <- dat.tax
df <- data.frame(x = log10(dat$MLD1), y = dat$T, measurementMethod = dat$measurementMethod, set = dat$set, scientificName = dat$scientificName)
	df[ ,1][is.infinite(df[,1])] = NA # remove -Inf
	df[ ,2][is.infinite(df[,2])] = NA # remove -Inf
	df <- df[complete.cases(df),]
	#
	#
	## Luo
	dt <- df[which(df$set=="Luo"),] # unique(dt$measurementMethod)
	da <- dt[which(dt$measurementMethod%in%c("Standard light microscopy", "Epifluorescence microscopy", "Standard light or epifluorescence microscopy", "Direct biomass analysis")),]
	x <- jitter ( x = da$x, factor = 1, amount = 0.01)
	y <- jitter ( x = da$y, factor = 1, amount = 0.2)
	points(x, y, cex = 0.7, col = "orange", pch = 1, xlab = "log(MLD)", ylab = "Temperature", xlim = c(1,2.5), ylim = c(-2, 34), lwd = lwd.size+0.7)
	da <- dt[which(dt$measurementMethod%in%c("qPCR_nifH_detection")),]
	x <- jitter ( x = da$x, factor = 1, amount = 0.005)
	y <- jitter ( x = da$y, factor = 1, amount = 0.12)
	points(x, y, cex = 0.7, col = "orange", pch = 25, xlab = "log(MLD)", ylab = "Temperature", xlim = c(1,2.5), ylim = c(-2, 34), lwd = lwd.size+0.7)
	#
	## Tan
	dt <- df[which(df$set=="Tan"),] # unique(dt$measurementMethod)
	da <- dt[which(dt$measurementMethod%in%c("qPCR_nifH_detection")),]
	x <- jitter ( x = da$x, factor = 1, amount = 0.01)
	y <- jitter ( x = da$y, factor = 1, amount = 0.2)
	points(x, y, cex = 0.7, col = "darkblue", pch = 25, xlab = "log(MLD)", ylab = "Temperature", xlim = c(1,2.5), ylim = c(-2, 34), lwd = lwd.size+0.7)
	#
	## Phy
	dt <-  df[which(df$set=="Phy"),]  # unique(dt$measurementMethod)
	x <- jitter ( x = dt$x, factor = 1, amount = 0.01)
	y <- jitter ( x = dt$y, factor = 1, amount = 0.2)
	points(dt$x, dt$y, col = "yellow2", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 1) # Assuming microscopy basaed (information on method, NA)
	#
	## Gra
	dt <- df[which(df$set=="Gra"),] # unique(dt$measurementMethod)
	x <- jitter ( x = dt$x, factor = 1, amount = 0.01)
	y <- jitter ( x = dt$y, factor = 1, amount = 0.2)
	points(dt[which(dt$measurementMethod%in%c("qPCR_nifH_detection")),]$x, dt[which(dt$measurementMethod%in%c("qPCR_nifH_detection")),]$y, col = "firebrick2", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 25)
	#
	## Gbi
	dt <-  df[which(df$set=="Gbi"),] # unique(dt$measurementMethod)
	x <- jitter ( x = dt$x, factor = 1, amount = 0.01)
	y <- jitter ( x = dt$y, factor = 1, amount = 0.2)
	points(dt$x, dt$y, col = "seagreen2", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 1) # Assuming microscopy basaed (information on method, NA)
	#
	## Obi
	dt <-  df[which(df$set=="Obi"),] # unique(dt$measurementMethod)
	x <- jitter ( x = dt$x, factor = 1, amount = 0.01)
	y <- jitter ( x = dt$y, factor = 1, amount = 0.2)
	points(dt$x, dt$y, col = "pink", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 1) # Assuming microscopy basaed (information on method, NA)
	#
	## Tara OTU based
	dt <-  df[which(df$set=="TarOTU"),] # unique(dt$measurementMethod)
	x <- jitter ( x = dt$x, factor = 1, amount = 0.01)
	y <- jitter ( x = dt$y, factor = 1, amount = 0.2)
	points(dt[which(dt$measurementMethod%in%c("SSU_rRNA_detection")),]$x, dt[which(dt$measurementMethod%in%c("SSU_rRNA_detection")),]$y, col = "slateblue1", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 22)
	#
	## M96 (Mar)
	dt <- df[which(df$set=="Mar"),] # unique(dt$measurementMethod)
	da <- dt[which(dt$measurementMethod%in%c("Epifluorescence microscopy")),]
	x <- jitter ( x = da$x, factor = 1, amount = 0.01)
	y <- jitter ( x = da$y, factor = 1, amount = 0.2)
	points(x, y, cex = 0.7, col = "firebrick3", pch = 1, xlab = "log(MLD)", ylab = "Temperature", xlim = c(1,2.5), ylim = c(-2, 34), lwd = lwd.size+0.7)
	da <- dt[which(dt$measurementMethod%in%c("qPCR_nifH_detection")),]
	x <- jitter ( x = da$x, factor = 1, amount = 0.005)
	y <- jitter ( x = da$y, factor = 1, amount = 0.12)
	points(x, y, cex = 0.7, col = "firebrick3", pch = 25, xlab = "log(MLD)", ylab = "Temperature", xlim = c(1,2.5), ylim = c(-2, 34), lwd = lwd.size+0.7)
	#
	## Grand total -----------------------------------------
	dat <- tot.clsets # set information is removed (i.e., the origin of the records), divergent methods are retained
	df <- data.frame(x = log10(dat$MLD1), y = dat$T, measurementMethod = dat$measurementMethod, set = dat$set, scientificName = dat$scientificName)
		df[ ,1][is.infinite(df[,1])] = NA # remove -Inf
		df[ ,2][is.infinite(df[,2])] = NA # remove -Inf
		df <- df[complete.cases(df),]
	x <- jitter ( x = df$x, factor = 1, amount = 0.0)
	y <- jitter ( x = df$y, factor = 1, amount = 0.0)
	# points(x, y, cex = 0.7, col = "grey46", pch = 19, xlab = "log(MLD)", ylab = "Temperature", xlim = c(1,2.5), ylim = c(-2, 34))
	points(x = 1.614071, y = 26.5604, cex = 2.5, col = "black", pch = 1, lwd = 44) # superniche (as reference)
	segments (  x0=mean(df$x)-sd(df$x), x1=mean(df$x)+sd(df$x) , y0=mean(df$y), y1=mean(df$y), lwd = 34, col = "khaki1" )
	segments (  y0=mean(df$y)-sd(df$y), y1=mean(df$y)+sd(df$y) , x0=mean(df$x), x1=mean(df$x), lwd = 34, col = "khaki1" )
	# legend("topright", legend= "Realized superniche", bty="n", text.font = 2, cex = lab.size)
	legend("topright", c("Full ocean (frequency)", "Super niche (reference)", "Focal taxon"), lty=c(NA, NA, NA), pch=c(22, 22, 22), col=c("grey80", "black", "khaki2"), bty = "n", lwd=42, cex= lab.size, horiz = F)
	# if(nrow(dat.tax)>1){points(pam(df, 1, metric = "euclidean")$medoids, cex = 2.5, col = "khaki3", pch = 1, lwd = 50)} # focal taxon
	if(nrow(dat.tax)>1){points(pam(df, 1, metric = "euclidean")$medoids, cex = 2.5, col = "khaki1", pch = 1, lwd = 48)} # focal taxon, if combined with khaki3 choose size 30	axis (1,at=seq(0, 3, 0.5), cex.axis = lab.size, lwd= lwd.box,  tck = tck.dist, mgp = c(1,dist, 0), las = 1)
	axis (2, at=seq(-5, 35, 5), cex.axis = lab.size, lwd= lwd.box,  tck = tck.dist, mgp = c(1,dist, 0), las = 1); box(lwd = lwd.box)
	mtext(expression("log(mixed-layer depth [m])"), side = 1,  line = 1.4, cex = lab.size)
	mtext(expression("Temperature ("*degree*"C)"), side=2, line=1.2, cex=lab.size)


# -------------------------environment log(no3-) inverse(mlpar) space: frame two-------------------------------------
	# Prepare data temperature and mld
	var1 <- vector (  )
	var2 <- vector (  )
	for (i in 1:length(env.stack)) {
	param1 <- subset (env.stack [[i]], "MLPAR1")
	param1[msk==0] <- NA
	var1 <- append (var1, values (param1))
	param2 <- subset (env.stack [[i]], "N")
	param2[msk==0] <- NA
	var2 <- append (var2, values (param2))
	}
	# Prepare shading environmental background
	var1.new <- var1 [-which (is.na(var1)  |  is.na(var2)) ]
	var2.new <- var2 [-which (is.na(var1)  |  is.na(var2)) ]
	env.space <- SpatialPoints (coords = cbind(log10(var1.new), var2.new))
	ranges <- raster(xmn = round (min(log10(var1.new)), digits = 3) ,
	xmx = round (max(log10(var1.new)), digits = 3),
	ymn = round (min(var2.new), digits = 2),
	ymx = round (max(var2.new), digits = 2),
	res = c(0.05, 0.5))
	env.space.new <- rasterize ( x = env.space , y= ranges, fun = "count" )
	ranges2 <- focal( env.space.new , w=matrix(c(1,1,1,1,1,1,2,2,2,1,1,2,3,2,1,1,2,2,2,1,1,1,1,1,1),5,5), fun=mean)
# Plot environmental background (as kernel)
par(fig=c(0.51 , 1 , 0 , 0.49), new=TRUE)
smoothScatter(var1,log10(var2),colramp=colorRampPalette(c("white","grey30"),50),bandwidth=c(0.25,0.06),nrpoints=0,nbin=300, axes=F, xlab = expression("MLPAR [E"*~m^-3 * d^-1~"]"), ylab = expression("log (NO3-) ["*~mu*mol/L~"]"), ylim = c(-4.5,1.6), xlim = c(50, 0))
# Plot ecological super niche
dat <- dat.tax
df <- data.frame(x = dat$MLPAR1, y = log10(dat$N), measurementMethod = dat$measurementMethod, set = dat$set)
	df[ ,1][is.infinite(df[,1])] = NA # remove -Inf
	df[ ,2][is.infinite(df[,2])] = NA # remove -Inf
	df <- df[complete.cases(df),] # 13808 x 4

	## Luo -----------------------------------------
	dt <- df[which(df$set=="Luo"),] # unique(dt$measurementMethod)
	da <- dt[which(dt$measurementMethod%in%c("Standard light microscopy", "Epifluorescence microscopy", "Standard light or epifluorescence microscopy", "Direct biomass analysis")),]
	x <- jitter ( x = da$x, factor = 1, amount = 0.6)
	y <- jitter ( x = da$y, factor = 1, amount = 0.06)
	points(x, y, cex = 0.7, col = "orange", pch = 1, xlab = "log(MLD)", ylab = "Temperature", xlim = c(1,2.5), ylim = c(-2, 34), lwd = lwd.size+0.7)
	da <- dt[which(dt$measurementMethod%in%c("qPCR_nifH_detection")),]
	x <- jitter ( x = da$x, factor = 1, amount = 0.6)
	y <- jitter ( x = da$y, factor = 1, amount = 0.06)
	points(x, y, cex = 0.7, col = "orange", pch = 25, xlab = "log(MLD)", ylab = "Temperature", xlim = c(1,2.5), ylim = c(-2, 34), lwd = lwd.size+0.7)
	#
	## Tan
	dt <- df[which(df$set=="Tan"),] # unique(dt$measurementMethod)
	da <- dt[which(dt$measurementMethod%in%c("qPCR_nifH_detection")),]
	x <- jitter ( x = da$x, factor = 1, amount = 0.6)
	y <- jitter ( x = da$y, factor = 1, amount = 0.06)
	points(x, y, cex = 0.7, col = "darkblue", pch = 25, xlab = "log(MLD)", ylab = "Temperature", xlim = c(1,2.5), ylim = c(-2, 34), lwd = lwd.size+0.7)
	#
	## Phy
	dt <-  df[which(df$set=="Phy"),]  # unique(dt$measurementMethod)
	x <- jitter ( x = dt$x, factor = 1, amount = 0.01)
	y <- jitter ( x = dt$y, factor = 1, amount = 0.2)
	points(dt$x, dt$y, col = "yellow2", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 1) # Assuming microscopy basaed (information on method, NA)
	#
	## Gra
	dt <- df[which(df$set=="Gra"),] # unique(dt$measurementMethod)
	x <- jitter ( x = dt$x, factor = 1, amount = 0.01)
	y <- jitter ( x = dt$y, factor = 1, amount = 0.2)
	points(dt[which(dt$measurementMethod%in%c("qPCR_nifH_detection")),]$x, dt[which(dt$measurementMethod%in%c("qPCR_nifH_detection")),]$y, col = "firebrick2", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 25)
	#
	## Gbi
	dt <-  df[which(df$set=="Gbi"),] # unique(dt$measurementMethod)
	x <- jitter ( x = dt$x, factor = 1, amount = 0.01)
	y <- jitter ( x = dt$y, factor = 1, amount = 0.2)
	points(dt$x, dt$y, col = "seagreen2", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 1) # Assuming microscopy basaed (information on method, NA)
	#
	## Obi
	dt <-  df[which(df$set=="Obi"),] # unique(dt$measurementMethod)
	x <- jitter ( x = dt$x, factor = 1, amount = 0.01)
	y <- jitter ( x = dt$y, factor = 1, amount = 0.2)
	points(dt$x, dt$y, col = "pink", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 1) # Assuming microscopy basaed (information on method, NA)
	#
	## Tara OTU based
	dt <-  df[which(df$set=="TarOTU"),] # unique(dt$measurementMethod)
	x <- jitter ( x = dt$x, factor = 1, amount = 0.01)
	y <- jitter ( x = dt$y, factor = 1, amount = 0.2)
	points(dt[which(dt$measurementMethod%in%c("SSU_rRNA_detection")),]$x, dt[which(dt$measurementMethod%in%c("SSU_rRNA_detection")),]$y, col = "slateblue1", cex = pt.size+0.1, lwd = lwd.size+0.7, pch = 22)
	#
	## M96 (Mar)
	dt <- df[which(df$set=="Mar"),] # unique(dt$measurementMethod)
	da <- dt[which(dt$measurementMethod%in%c("Epifluorescence microscopy")),]
	x <- jitter ( x = da$x, factor = 1, amount = 0.01)
	y <- jitter ( x = da$y, factor = 1, amount = 0.2)
	points(x, y, cex = 0.7, col = "firebrick3", pch = 1, xlab = "log(MLD)", ylab = "Temperature", xlim = c(1,2.5), ylim = c(-2, 34), lwd = lwd.size+0.7)
	da <- dt[which(dt$measurementMethod%in%c("qPCR_nifH_detection")),]
	x <- jitter ( x = da$x, factor = 1, amount = 0.005)
	y <- jitter ( x = da$y, factor = 1, amount = 0.12)
	points(x, y, cex = 0.7, col = "firebrick3", pch = 25, xlab = "log(MLD)", ylab = "Temperature", xlim = c(1,2.5), ylim = c(-2, 34), lwd = lwd.size+0.7)
	#
	## Grand total -----------------------------------------
	dat <- tot.clsets # set information is removed (i.e., the origin of the records), divergent methods are retained
	df <- data.frame(x = dat$MLPAR1, y = log10(dat$N), measurementMethod = dat$measurementMethod, set = dat$set)
		df[ ,1][is.infinite(df[,1])] = NA # remove -Inf
		df[ ,2][is.infinite(df[,2])] = NA # remove -Inf
		df <- df[complete.cases(df),]
	x <- jitter ( x = df$x, factor = 1, amount = 0.0)
	y <- jitter ( x = df$y, factor = 1, amount = 0.0)
	# points (x, y, cex = 0.7, col = "grey46", pch = 19,  xlab = expression("MLPAR [E"*~m^-3 * d^-1~"]"), ylab = expression("log (NO3-) ["*~mu*mol/L~"]"), ylim = c(-4.5,1.6), xlim = c(50, 0))  # original points
	points(x= 20.4871, y=-0.7515354, cex = 2.5, col = "black", pch = 1, lwd = 44) # superniche (as reference)
	segments (  x0=mean(df$x)-sd(df$x), x1=mean(df$x)+sd(df$x) , y0=mean(df$y), y1=mean(df$y), lwd = 34, col = "khaki1" )
	segments (  y0=mean(df$y)-sd(df$y), y1=mean(df$y)+sd(df$y) , x0=mean(df$x), x1=mean(df$x), lwd = 34, col = "khaki1" )
	# legend("topright", legend= "Realized superniche", bty="n", text.font = 2, cex = lab.size)
	legend("topleft", c("Full ocean (frequency)", "Super niche (reference)", "Focal taxon"), lty=c(NA, NA, NA), pch=c(22, 22, 22), col=c("grey80", "black", "khaki2"), bty = "n", lwd=42, cex= lab.size, horiz = F)
	# if(nrow(dat.tax)>1){points(pam(df, 1, metric = "euclidean")$medoids, cex = 2.5, col = "khaki3", pch = 1, lwd = 50)} # focal taxon
	if(nrow(dat.tax)>1){points(pam(df, 1, metric = "euclidean")$medoids, cex = 2.5, col = "khaki1", pch = 1, lwd = 48)} # focal taxon, if combined with khaki3 choose size 30
axis (1,at=seq(50, 0, -10), cex.axis = lab.size, lwd= lwd.box,  tck = tck.dist, mgp = c(1,dist, 0), las = 1)
axis (2, at=seq(-4.3, 1.6, 1), cex.axis = lab.size, lwd= lwd.box,  tck = tck.dist, mgp = c(1,dist, 0), las = 1); box(lwd = lwd.box)
mtext(expression("MLPAR (E"*~m^-3 * d^-1~")"), side = 1,  line = 1.4, cex = lab.size)
mtext(expression("log (NO3-) ["*~mu*mol/L~"]"), side=2, line=1.4, cex=lab.size)
#
#
#
dev.off()

# Close loop across taxa
print(paste(q))
}
