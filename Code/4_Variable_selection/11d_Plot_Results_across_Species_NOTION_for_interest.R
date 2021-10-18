### PLOT VARIABLE RANKING (MEAN GAM RESULTS PLUS MINUS SD)

# $Date: 2020-10-23
# Authors: Damiano Righetti, rdamiano@ethz.ch

# Description: Rank the variables based on GLM, GAM, RF models. Create a table for each phytoplankton group:
# 1. GLM points - overlapping 
# 2. GAM points - overlapping
# 3. RF points - overlapping
# Same for sites - overlapping

# Notes:
# "Points approach" or "pts" (in filenames) refers to a selection approach of background data at 1° monthly resolution whereby different species are considered (this strategy was later discarded, as it introduces a richness signal into the pseudoabsences selection procedure). Example: In this approach, a monthly 1° cell containing 10 observation records from ten different species provides potentially ten pseudoabsences. 
# "Sites approach" or "sit" (in filenames) stands for selection of background data from "observational cells" -> briefly "sites" (i.e., pooling all species or taxa present within each 1° monthly resolution cell), as the basis to select pseodoabsences ("sit" in filenames). Example: In this approach, a monthly 1° cell containing 10 observation records from ten different species will just provide one potential pseudoabsence, as species are pooled. 
# The difference between total bg and group bg seems not evident. recheck sourcing of results.

### =========================================================================
### Initialize system
### =========================================================================

# Set locale
setwd("/Users/rdamiano/Desktop/Data/8_Variable_Selection")
# kryo or clusters
if (getwd() == "/home/rdamiano/") {setwd("/home/rdamiano/Data/8_Variable_Selection")}

# Input and output
input.dir <- "./Output_11"
 
### =========================================================================
### Load results
### =========================================================================

# gam, total bg, sites
gam.sit.diazos <- read.csv(paste0(input.dir,"/Analysis ranking/diazos_minobs_14.csv"))
gam.sit.diazos$pCO2 <- as.character(gam.sit.diazos$pCO2)
#
gam.sit.hosts <- read.csv(paste0(input.dir,"/Analysis ranking/hosts_minobs_14.csv"))
gam.sit.hosts$pCO2 <- as.character(gam.sit.hosts$pCO2)
#
gam.sit.total <- read.csv(paste0(input.dir,"/Analysis ranking/total_minobs_14.csv"))
gam.sit.total$pCO2 <- as.character(gam.sit.total$pCO2)

### =========================================================================
### BAR PLOT | METRIC OF SKILL (GOODNESS OF FIT) FOR DIAZOTROPH TAXA
### =========================================================================

### Specifications
line.width <- 8
sz <- 0.5

## Plot for quick exploration
vec.names <- names(gam.sit.diazos  [  which(gam.sit.diazos$value == "adj.Rsq" & gam.sit.diazos$mod == "gam_only_tot_bg"),  3:28  ]   ) 
vec.mean.vals <-  as.numeric(gam.sit.diazos  [  which(gam.sit.diazos$value == "adj.Rsq" & gam.sit.diazos$mod == "gam_only_tot_bg"),  3:28  ]   ) 
vec.sd.vals <-  as.numeric(gam.sit.diazos  [  which(gam.sit.diazos$value == "SD.adj.Rsq" & gam.sit.diazos$mod == "gam_only_tot_bg"),  3:28  ]   ) 
df.essence <- data.frame("var" = vec.names, "mean.adj.Rsq" = vec.mean.vals, "SD.adj.Rsq" = vec.sd.vals, "upper" = vec.mean.vals + vec.sd.vals, "lower" = vec.mean.vals - vec.sd.vals)
df.essence <- df.essence[order(df.essence$mean.adj.Rsq, decreasing = T), ]
df.essence.diazos <- df.essence # for ordering below
# barplot(df.essence$mean.adj.Rsq)

## Plot for detailed analysis
png(filename = paste0("./Output_11/Figures single predictor skill/Predictors_mean_adj.Rsqared_diazos.png"), width = 2700, height = 2700, units = "px", pointsize = 80)
#
opar <- par(lwd = line.width)
barplot(df.essence$mean.adj.Rsq, beside = TRUE, cex.axis = sz, cex.lab = sz, ylim = c(0, 0.8), xlim = c(0.45, 32), density = c(200, 200), 
lwd = 2, angle = c(30, 0), space = rep(0.2, length(df.essence$mean.adj.Rsq)), las = 2, names.arg = df.essence$var, axes = F, cex.names = sz)
## Add standard error of the mean (s.e.m) or confidence interval (ci)
for(  i in 1:length(df.essence$mean.adj.Rsq)  ){
arrows(x0 = 0.5+(i-1)+i*0.2, y0 = c(df.essence$mean.adj.Rsq[i]), x1 = 0.5+(i-1)+i*0.2, y1 = df.essence$upper[i],  length = 0, lwd = line.width+2) # upper CI of variable 1
arrows(x0 = 0.5+(i-1)+i*0.2, y0 = c(df.essence$mean.adj.Rsq[i]), x1 = 0.5+(i-1)+i*0.2, y1 = df.essence$lower[i],  length = 0, lwd = line.width+2) # upper CI of variable 1
}
## Add axis
axis(2, labels =T, cex.axis = sz + 0.1, las = 2, tck = 0.015, lwd = line.width, mgp = c(1,0.4, 0))
mtext( expression("Mean skill of gam to predict diazotrophic taxa (adjusted "*italic(R)^2*")"), side = 2, line = 2, cex = sz + 0.1) # Proportion of variance associated with one or more main effects (effect size)
#
dev.off()

### =========================================================================
### BAR PLOT | METRIC OF SKILL (GOODNESS OF FIT) FOR HOST TAXA
### =========================================================================

### Specifications
line.width <- 8
sz <- 0.5

## Plot for quick exploration
vec.names <- names(gam.sit.hosts  [  which(gam.sit.hosts$value == "adj.Rsq" & gam.sit.hosts$mod == "gam_only_tot_bg"),  3:28  ]   ) 
vec.mean.vals <-  as.numeric(gam.sit.hosts  [  which(gam.sit.hosts$value == "adj.Rsq" & gam.sit.hosts$mod == "gam_only_tot_bg"),  3:28  ]   ) 
vec.sd.vals <-  as.numeric(gam.sit.hosts  [  which(gam.sit.hosts$value == "SD.adj.Rsq" & gam.sit.hosts$mod == "gam_only_tot_bg"),  3:28  ]   ) 
df.essence <- data.frame("var" = vec.names, "mean.adj.Rsq" = vec.mean.vals, "SD.adj.Rsq" = vec.sd.vals, "upper" = vec.mean.vals + vec.sd.vals, "lower" = vec.mean.vals - vec.sd.vals)
df.essence <- df.essence[ match(df.essence.diazos$var, df.essence$var) , ]
# barplot(df.essence$mean.adj.Rsq)

## Plot for detailed analysis
png(filename = paste0("./Output_11/Figures single predictor skill/Predictors_mean_adj.Rsqared_hosts.png"), width = 2700, height = 2700, units = "px", pointsize = 80)
#
opar <- par(lwd = line.width)
barplot(df.essence$mean.adj.Rsq, beside = TRUE, cex.axis = sz, cex.lab = sz, ylim = c(0, 0.8), xlim = c(0.45, 32), density = c(200, 200), 
lwd = 2, angle = c(30, 0), space = rep(0.2, length(df.essence$mean.adj.Rsq)), las = 2, names.arg = df.essence$var, axes = F, cex.names = sz)
## Add standard error of the mean (s.e.m) or confidence interval (ci)
for(  i in 1:length(df.essence$mean.adj.Rsq)  ){
arrows(x0 = 0.5+(i-1)+i*0.2, y0 = c(df.essence$mean.adj.Rsq[i]), x1 = 0.5+(i-1)+i*0.2, y1 = df.essence$upper[i],  length = 0, lwd = line.width+2) # upper CI of variable 1
arrows(x0 = 0.5+(i-1)+i*0.2, y0 = c(df.essence$mean.adj.Rsq[i]), x1 = 0.5+(i-1)+i*0.2, y1 = df.essence$lower[i],  length = 0, lwd = line.width+2) # upper CI of variable 1
}
## Add axis
axis(2, labels =T, cex.axis = sz + 0.1, las = 2, tck = 0.015, lwd = line.width, mgp = c(1,0.4, 0))
mtext( expression("Mean skill of gam to predict host taxa (adjusted "*italic(R)^2*")"), side = 2, line = 2, cex = sz + 0.1) # Proportion of variance associated with one or more main effects (effect size)
#
dev.off()

