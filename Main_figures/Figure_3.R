## This script computes Figure 3 from Dominic Eriksson et al.(submitted). Correlation between
##  diazotroph richness and nitrogen fixation rates.


#   Input files:
#           1. Nitrogen fixation rates from Wang et al. (2019) and Shaozaki et al. (2023)
#           2. Annual diazotroph richness. 


#   Output files:
#           2. PNG files for each panel.


# Strategy: 


# Figure Caption:   Relationship between diazotroph species richness and biological nitrogen fixation (BNF). 
#                   Shown are correlations between nitrogen fixation rates and global annual diazotroph richness. 
#                   A) Nitrogen fixation rates originate from (44). The black line indicates a 2nd order polynomial 
#                   fit. Further indicated are Spearman’s rank correlation coefficient and the statistical p-value. 
#                   B) Correlation between grid cells that contain in situ nitrogen fixation measurements compiled by 
#                   Shao et al. (2023) and projected annual diazotroph richness from an ensemble of species distribution models. 
#                   Spearman correlation coefficient is given in the top left corner with the associated p-value. Grey shading 
#                   indicates the 0.95 confidence interval.

# Author:   Dominic Eriksson
#		    Environmental Physics Group, UP
#			ETH, Zurich
#			Switzerland

# 9th of October 2023, deriksson@ethz.ch --------------------------------------------------------------------------


### =========================================================================
### Preparations
### =========================================================================
# Clear workspace
rm(list = ls())

# Libraries
library(ggplot2)
library(ggpmisc)

## Directories
# Inputs
input.dir.nfix_wang <- paste0("/home/deriksson/Projects/Notion_DE/Data/n2fix_Wang_et_al_2019/1_Output/") # nitrogen fixation data (Wang 2019)
# Outputs
wd_out <- '/home/deriksson/Projects/Notion_DE/Figures/3_Output/'



### ==============================================================
### Figure 3 panel A: Correlation between annual ensemble richness and annual biological nitrogen fixation from Wang et al. 2019
### ==============================================================
### Load data
## Nitrogen fixation rate wang
nfix_wang <-  raster(paste0(input.dir.nfix_wang, "depth_integrated_n2_fix_top_2_layers.grd"))

# Load richness (esemble over all SDM and one backround selection)
df <- brick(paste0('/net/kryo/work/deriksson/Projects/Notion_DE/Code/8_Diversity/3_Output/PA_majority_vote/Total_dataset/div_annually/pa/Ensemble_richness_acrossAllModels.grd'))
df <- as.data.frame(df, xy = TRUE)
names(df)[3] <- "richness"
df$nfix_wang <- as.data.frame(nfix_wang)[, 1]
# Remove negative values
df <- df[which(df$richness > 0), ]
df <- df[which(df$nfix_wang > 0), ]

# Fit polynomial
p_wang <- ggplot(data = df, aes(x = richness, y = nfix_wang)) +
    geom_point(color = "grey", alpha = 0.5) +
    stat_poly_line(method = 'lm', formula = y ~ poly(x, 2), color = 'black') +
    stat_cor(method = "spearman", label.x = 0.1, label.y = 400, size = 8) +
    stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 8) +
    # stat_cor(method = "spearman", label.x = 0.1, label.y = 300, size = 22) 

    theme(
        panel.background = element_rect(fill = "white"),
        axis.text.x = element_text(size = 32),
        axis.text.y = element_text(size = 32),
        axis.title.x = element_text(size = 32),
        axis.title.y = element_text(size = 32)) +

    ylab(expression("N2 fix (mmol N m-2 yr-1)")) + 
    theme(legend.position = "none", text = element_text(size = 5)) + 
    labs(tag = "A", size = 32)


# Save
ggsave(
    plot = p_wang,
    filename = paste0(wd_out, '3A_Correlation_WangVsRichness.png'),
    dpi = 300
)


### ==============================================================
### 3.B Plot nitrogen fixation measurements in situ with richness - Luo 2023
### ==============================================================

# Luo nfix 2014
df_luo2014 <- read.csv('/home/deriksson/Projects/Notion_DE/Data/n2_fix_Luo_et_al_2014/Copy_to_work_with/dat_N2_fixation.csv', sep = ';')
df_luo2014 <- df_luo2014[, c(4:7, 12)]
df_luo2014[, 5] <- as.numeric(df_luo2014[, 5])
df_luo2014$x <- round(df_luo2014$LONGITUDE + 0.5) - 0.5
df_luo2014$y <- round(df_luo2014$LATITUDE + 0.5) - 0.5
df_luo2014$id <- paste0( df_luo2014$x, '_', df_luo2014$y )
names(df_luo2014)[5] <- 'Total.N2.Fixation..μmol.N.m.2.d.1.'

# Load Luo 2023 nfix
df_luo2023 <- read.csv('/home/deriksson/Projects/Notion_DE/Data/Zhibo_Shao_2023/Copy_to_work_with/NitrogenFixation_integrated_above24h.csv', row.names = NULL, sep = ';')
df_luo2023 <- df_luo2023[, c(4:7, 13)]
df_luo2023[, 5] <- as.numeric(df_luo2023[, 5])
df_luo2023$x <- round(df_luo2023$LONGITUDE + 0.5) - 0.5
df_luo2023$y <- round(df_luo2023$LATITUDE + 0.5) - 0.5
df_luo2023$id <- paste0( df_luo2023$x, '_', df_luo2023$y )


# Load richness
r <- brick('/home/deriksson/Projects/Notion_DE/Code/8_Diversity/3_Output/Final_ensemble_used/Ensemble_AcrossEveryModelandMonth.grd')
r_cd_ncd <- brick('/home/deriksson/Projects/Notion_DE/Code/8_Diversity/3_Output/Ensemble_mean/Total_dataset/div_annually/pa/Ensemble_richness_acrossAllModels_CyanosVsNonCyanos.grd')
rr <- stack(r, r_cd_ncd)
#
df_rich <- as.data.frame(rr, xy = TRUE)
df_rich$id <- paste0( df_rich$x, '_', df_rich$y )


# Rbind luo dataframes
df_luo <- rbind(df_luo2014, df_luo2023)
names(df_luo) <- c('Date', 'latitude', 'longitude', 'depth', 'nfix', 'x', 'y', 'id')

# Aggregate
df_luo <- aggregate(df_luo$nfix, by = list(df_luo$id), mean)
names(df_luo) <- c('id', 'nfix')
df_luo$log_nfix <- log10(df_luo$nfix)
# Merge back together
df <- merge(df_luo, df_rich, by = 'id')


# Read shapefile with geometries of ocean basins
sf.ocean <- st_read("/net/kryo/work/deriksson/Projects/Notion_DE/Data/Ocean_regions_shapefile/GOaS_v1_20211214/goas_v01.shp")

# Convert to from micromol per square meter per day to  mmol N m–2 yr–
df$nfix_annual <- df$nfix * 0.001 * 365
df$log_nfix_annual <- log10(df$nfix_annual)
df <- df[which(df$nfix != 0), ]

# Total richness 
gg <- ggplot(data = df, aes(x = mean_richness, y = log_nfix_annual)) + geom_point() +
    stat_poly_line(method = 'lm', formula = y ~ poly(x), color = 'black') +
    stat_poly_eq(use_label(c("eq", "adj.R2", "p")), size = 10) +
    stat_cor(method = "spearman", label.x = 0.1, label.y = 2.8, size = 8) +
    ylab(expression("log. N2 fix (mmol N m-2 yr-1)")) + 
    xlab('Normalized diazotroph richness') +
      # Some specific settings
  theme(
    # Plot background color
    panel.background = element_rect(fill = "white"),
    axis.text.x = element_text(size = 32),
    axis.text.y = element_text(size = 32),
    axis.title.x = element_text(size = 32),
    axis.title.y = element_text(size = 32))

# Save plot
fln <- paste0( wd_out, 'Correlation_inSituNfix_SDMrichness.png' )
ggsave(
    plot = gg,
    filename = fln,
    dpi = 300
)


###==============================================================
### END
###==============================================================