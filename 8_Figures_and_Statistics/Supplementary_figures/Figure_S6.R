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



### ===========================================================================
### SCATTERPLOT In situ N2 measurement correlation with richness  
### ===========================================================================  
## Now we use nitrogen fixation rate observations - Let's see if we find a relationship between measured 
## marine nitrogen fixation rates: the nfix layer in the raster is micromol m-2/m-3 yr-1
# Load Table
r <- brick("/hrome/deriksson/Projects/Notion_DE/Code/Tables/1_Output/FinalTable_AllMetrics.grd")


#
df <- as.data.frame(r, xy = TRUE)
df <- data.frame(
  x = df$x, 
  y = df$y, 
  richnessTotal = df$Diversity_Richness_ensembleMean_paBased, 
  richnessCyanos = df$Diversity_Richness_ensembleMean_paBased_cyanos, 
  richnessNCD = df$Diversity_Richness_ensembleMean_paBased_ncd, 
  nfix = df$nfix)

# # Total community ------------------------------
df <- na.omit(df)
df <- df[which(df$nfix > 0.0001), ]
df$nfix <- df$nfix * 365/1000
df$log <- log10(df$nfix)


# Fit linear regression model
fit <- lm(log ~ richnessTotal, data = df)

# Extract R-squared value
r_squared <- summary(fit)$r.squared

# Extract coefficients of the linear regression model
intercept <- coef(fit)[1]
slope <- coef(fit)[2]

# Create formula text
formula_text <- paste("y =", round(intercept, 2), "+", round(slope, 2), "* x")

# Perform correlation significance test
cor_test_result <- cor.test(df$richnessTotal, df$log, method = "spearman", exact = FALSE, continuity = TRUE)


# Create scatterplot with linear regression line
gg_scatter <- ggplot(data = df, aes(x = richnessTotal, y = log)) +
  geom_point() +  # Add scatterplot points
  geom_smooth(method = 'lm', se = FALSE) +  # Add linear regression line without confidence interval
  labs(x = "X-axis Label", y = "Y-axis Label") +  # Add axis labels
  annotate("text", x = max(df$richnessTotal), y = max(df$log), 
           label = formula_text, hjust = 1, vjust = 1, size = 22) +
  annotate("text", x = max(df$richnessTotal), y = min(df$log), 
           label = paste("Spearman's rho =", round(cor_test_result$estimate, 2)), 
           hjust = 1, vjust = 0, size = 18) +  # Add Spearman rank correlation coefficient to the plot with increased size
  annotate("text", x = max(df$richnessTotal), y = max(df$log) - 0.5, 
           label = paste("p-value =", signif(cor_test_result$p.value, 3)), 
           hjust = 1, vjust = 1, size = 18) +
  theme(
    axis.title = element_text(size = 70),  # Increase font size of axis labels
    axis.text = element_text(size = 50),  # Increase font size of tick labels
    legend.text = element_text(size = 32)  # Increase font size of legend text
  ) +
    # Axis labels
  ylab( expression(paste( ~log[10]~~N[2]~"fixation flux (mmol N"~ m^-2*y^-1*")")) ) + 
  xlab("Annual richness") 

## Save
fln <- "/home/deriksson/Projects/Notion_DE/Code/Tables/Visualizations/"
fln <- paste0( fln, "Scatterplot_correlation_totalDiaz_inSituN2Fixation_and_annualRichness.png" )
#
ggsave(
  filename = fln,
  plot = gg_scatter,
  bg = "white",
  width = 15,
  height = 25,
  dpi = 300
)


# # cyanos ------------------------------
# Fit linear regression model
fit <- lm(log ~ richnessCyanos, data = df)

# Extract R-squared value
r_squared <- summary(fit)$r.squared

# Extract coefficients of the linear regression model
intercept <- coef(fit)[1]
slope <- coef(fit)[2]

# Create formula text
formula_text <- paste("y =", round(intercept, 2), "+", round(slope, 2), "* x")

# Perform correlation significance test
cor_test_result <- cor.test(df$richnessCyanos, df$log, method = "spearman", exact = FALSE, continuity = TRUE)

# Create scatterplot with linear regression line
gg_scatter <- ggplot(data = df, aes(x = richnessCyanos, y = log)) +
  geom_point() +  # Add scatterplot points
  geom_smooth(method = 'lm', se = FALSE) +  # Add linear regression line without confidence interval
  labs(x = "X-axis Label", y = "Y-axis Label") +  # Add axis labels
  annotate("text", x = max(df$richnessCyanos), y = max(df$log), 
           label = formula_text, hjust = 1, vjust = 1, size = 22) +
  annotate("text", x = max(df$richnessCyanos), y = min(df$log), 
           label = paste("Spearman's rho =", round(cor_test_result$estimate, 2)), 
           hjust = 1, vjust = 0, size = 18) +  # Add Spearman rank correlation coefficient to the plot with increased size
  annotate("text", x = max(df$richnessCyanos), y = max(df$log) - 0.5, 
           label = paste("p-value =", signif(cor_test_result$p.value, 3)), 
           hjust = 1, vjust = 1, size = 18) + # add p value
  theme(
    axis.title = element_text(size = 70),  # Increase font size of axis labels
    axis.text = element_text(size = 50),  # Increase font size of tick labels
    legend.text = element_text(size = 50)  # Increase font size of legend text
  ) +
    # Axis labels
  ylab( expression(paste( ~log[10]~~N[2]~"fixation flux (mmol N"~ m^-2*y^-1*")")) ) + 
  xlab("Annual richness") 

## Save
fln <- "/home/deriksson/Projects/Notion_DE/Code/Tables/Visualizations/"
fln <- paste0( fln, "Scatterplot_correlation_cyanoDiaz_inSituN2Fixation_and_annualRichness.png" )
#
ggsave(
  filename = fln,
  plot = gg_scatter,
  bg = "white",
  width = 15,
  height = 25,
  dpi = 300
)



# # non-cyanos ------------------------------
# Fit linear regression model
fit <- lm(log ~ richnessNCD, data = df)

# Extract R-squared value
r_squared <- summary(fit)$r.squared

# Extract coefficients of the linear regression model
intercept <- coef(fit)[1]
slope <- coef(fit)[2]

# Create formula text
formula_text <- paste("y =", round(intercept, 2), "+", round(slope, 2), "* x")

# Perform correlation significance test
cor_test_result <- cor.test(df$richnessNCD, df$log, method = "spearman", exact = FALSE, continuity = TRUE)



# Create scatterplot with linear regression line
gg_scatter <- ggplot(data = df, aes(x = richnessNCD, y = log)) +
  geom_point() +  # Add scatterplot points
  geom_smooth(method = 'lm', se = FALSE) +  # Add linear regression line without confidence interval
  labs(x = "X-axis Label", y = "Y-axis Label") +  # Add axis labels
  annotate("text", x = max(df$richnessNCD), y = max(df$log), 
           label = formula_text, hjust = 1, vjust = 1, size = 22) +
  annotate("text", x = max(df$richnessNCD), y = min(df$log), 
           label = paste("Spearman's rho =", round(cor_test_result$estimate, 2)), 
           hjust = 1, vjust = 0, size = 18) +  # Add Spearman rank correlation coefficient to the plot with increased size
  annotate("text", x = max(df$richnessNCD), y = max(df$log) - 0.5, 
           label = paste("p-value =", signif(cor_test_result$p.value, 3)), 
           hjust = 1, vjust = 1, size = 18) + # add p value
  theme(
    axis.title = element_text(size = 70),  # Increase font size of axis labels
    axis.text = element_text(size = 50),  # Increase font size of tick labels
    legend.text = element_text(size = 50)  # Increase font size of legend text
  ) +
    # Axis labels
  ylab( expression(paste( ~log[10]~~N[2]~"fixation flux (mmol N"~ m^-2*y^-1*")")) ) + 
  xlab("Annual richness") 

## Save
fln <- "/home/deriksson/Projects/Notion_DE/Code/Tables/Visualizations/"
fln <- paste0( fln, "Scatterplot_correlation_nonCyanoDiaz_inSituN2Fixation_and_annualRichness.png" )
#
ggsave(
  filename = fln,
  plot = gg_scatter,
  bg = "white",
  width = 15,
  height = 25,
  dpi = 300
)


