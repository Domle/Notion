

### In this script we plot main figure 1 from the Manuscript "Nitrogen fixation increases with diazotroph richness in the global ocean"

### NOTE: 	When running the code, please always double check input and output directories, to ensure you know where exactly files get saved so they
###			    can be referred as the correct input files within each consecutive script.

# Author: 	
#     Dominic Eriksson
#			Environmental Physics Group, UP
# 		ETH Zurich
#			Switzerland


# Input files:
#       1. shapefile goas

# Output files:
#       1. Raster stack as .grd file, later used for final plotting

# derikssong@ethz.ch; 18th of June 2024



###======================================================================
### Plot heatmap
###======================================================================
# Libraries
library(ggplot2)

# Directories
wd_dat <- "/home/deriksson/Projects/Notion_DE/Code/5_Variable_selection/3_Output/"
wd_out <- "/net/kryo/work/deriksson/Projects/Notion_DE/Code/5_Variable_selection/5_Output/"

# Vectors directories
vec_folders <- c("Total_dataset", "MicroscopyBased_dataset", "SequenceBased_dataset")
vc_vars <- c("T","Sal","N","P","Si",
			"MLD1", "MLD2","PAR","Chl","Wind.CCMP","pCO2",
			"MLPAR1", "MLPAR2", "Nstar", "Sistar",
			"dT_dt", "dN_dt", "dP_dt", "dSi_dt","dMLD1_dt",
			"logMLD1", "logMLD2", "logChl", "logN", "logP", "logSi")
			# I remove NPP, as well reflected in Chl
#
vec_nms <- c(
  "Total_overlapping",
  "Total_non_overlapping",
  "Group_overlapping",
  "Group_non_overlapping",
  "Cruise_overlapping",
  "Cruise_non_overlapping"
)

## Load data
# Tax table - Needed for assigning groups
tax_table <- read.csv("/home/deriksson/Projects/Notion_DE/Data/commons/Taxon_table.csv", sep = ";")
# set v for folder choice
v = 1
list.files(wd_dat)
#
file_names <- list.files(paste0(wd_dat, vec_folders[v], "/"))
# file_names <- grep("Total", file_names, value = TRUE)
l <- list()
# for(i in seq_along(file_names)){
for(i in 5:6){
  print(i)
  l[[length(l) + 1]] <- read.csv(paste0(wd_dat, vec_folders[v], "/", file_names[i]))
}
# Merge
df <- do.call("rbind", l)

x <- df[which(df$value == "rank"), ]
x <- df[, c(vc_vars, "taxon")]
x <- aggregate(. ~ taxon, data = x, mean)
x[, vc_vars] <- round(x[, vc_vars], digits = 0)

# We need to harmonize the variable names so they are all the same in the manuscript
vc_vars_harm <- c("SST","S","N","P","Si",
			"MLD1", "MLD2","PAR","Chl","SSW","pCO2",
			"MLPAR1", "MLPAR2", "Nstar", "Sistar",
			"dT/dt", "dN/dt", "dP/dt", "dSi/dt","dMLD1/dt",
			"logMLD1", "logMLD2", "logChl", "logN", "logP", "logSi")
#
names(x)[2:ncol(x)] <- vc_vars_harm



# Wide to long format 
library(tidyr)

df_long <- gather(x, EnvPar, Ranking, SST:logSi, factor_key = TRUE)

# Plot heatmap
heatmap_PredRanking <- ggplot(df_long, aes(x = taxon, y = EnvPar, fill = Ranking)) + 
  geom_tile(color = 'white', lwd = 1.5, linetype = 1) +
  coord_fixed() + 
  theme(axis.text.x = element_text(angle = 90)) + 
  ylab('Environmental Parameter') + 
  xlab('Taxon') + 
  ggtitle('Environmental Predictor Ranking') +
  scale_fill_gradient(limits = c(0, 13), low = "orange", 
    high = "blue",
    breaks = c(1, 3, 5, 7, 9, 11, 13)
    ) 


# Save plot
ggsave(
  filename = "/home/deriksson/Projects/Notion_DE/Code/5_Variable_selection/5_Output/Environmental_predictor_ranking/Heatmap2D_PredRanking.png",
  plot = heatmap_PredRanking,
  dpi = 300
)



### =========================================================
### DO A Canonical Correspondence Analysis (CCA) on the Annual HSI
### =========================================================
# Load Table
r <- brick("/home/deriksson/Projects/Notion_DE/Code/Tables/1_Output/FinalTable_AllMetrics.grd")

sort(names(r))
vec_env <- grep("EnvPar", names(r), value = TRUE)
vec_env <- vec_env[-c(1:2)]
#
vec_spec <- grep("Species", names(r), value = TRUE)
vec_spec <- grep("hsi", vec_spec, value = TRUE)
vec_spec <- grep("Mean", vec_spec, value = TRUE)

# Convert to dataframe
df <- as.data.frame(r)
df <- df[, c( vec_spec, vec_env )]
df <- na.omit(df)
# Only use environmental parameters used to model the taxa
d <- read.csv("/home/deriksson/Projects/Notion_DE/Code/6_SDM_fit/1_Output/Dataframe_ResponseCurves.csv")
vec_env <- unique(d$predictor)
vec_env <- paste0( "EnvPar_", vec_env )

df_env <- df[, vec_env]
df_env <- df_env[, !(names(df_env) %in% c("EnvPar_MLD2"))]
df_env <- df_env[, !(names(df_env) %in% c("EnvPar_dMLD1_dt"))]
df_env <- df_env[, !(names(df_env) %in% c("EnvPar_MLPAR2"))]
df_env <- df_env[, !(names(df_env) %in% c("EnvPar_dSi_dt"))]
df_env <- df_env[, !(names(df_env) %in% c("EnvPar_dN_dt"))]
df_env <- df_env[, !(names(df_env) %in% c("EnvPar_dP_dt"))]
df_env <- df_env[, !(names(df_env) %in% c("EnvPar_dT_dt"))]
df_env <- df_env[, !(names(df_env) %in% c(grep("log", names(df_env), value = TRUE)))]
names(df_env) <- gsub("EnvPar_", "", names(df_env))

df_spec <- df[, vec_spec]
# short names for species
names(df_spec) <- c(
  "Calotrhix",
  "Gamma-A",
  "HBD-02",
  "HBD-03",
  "HBD-04",
  "HBD-05",
  "HBD-06",
  "R. intracellularis",
  "Richelia",
  "T. erythraeum",
  "T. thiebautii",
  "UCYN-A1",
  "UCYN-A2",
  "UCYN-B",
  "UCYN-C"
)

# CCA
spec.cca <- cca( df_spec ~ ., df_env )
# summary(spec.cca)
# RsquareAdj(spec.cca)

# Define width and height of the plot in inches
width_inches <- 10  # Adjust as needed
height_inches <- 8  # Adjust as needed

# Open a PNG graphics device with DPI = 300
png("/home/deriksson/Projects/Notion_DE/Code/10_Niche_analysis/Ordination/NMDS.png", 
    width = width_inches, 
    height = height_inches, 
    units = "in", 
    res = 300)

# Plot the spec.cca object
plot(spec.cca, scaling = 2, display = c("sp", "cn"))

# Close the graphics device
dev.off()

