### In this script we plot main figure 2 from the Manuscript "Nitrogen fixation increases with diazotroph richness in the global ocean"

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


# Libraries
library(raster)
library(gridExtra)


# Load Table
r <- brick("/home/deriksson/Projects/Notion_DE/Code/Tables/1_Output/FinalTable_AllMetrics.grd")

## This script plot all plots in the manuscript and additionally used anywhere else. 
dark_green_hex <- "#2ca25f" # Cyanobacterial diazotroph
dark_blue_hex <- "#c51b8a" # Non-cyanobacterial diazotroph community
orange_hex <- "#f1a340" # Total diazotroph community


## BOXPLOTS RICHNESS vs. N2-Fixation
# Boxplot richness total diaz. community and nitrogen fixation rate Wang et al. 2019
sort(names(r))
stack <- stack( 
  r$Diversity_Richness_ensembleMean_paBased, 
  r$Diversity_Richness_ensembleMean_paBased_ncd,
  r$Diversity_Richness_ensembleMean_paBased_cyanos,
  r$nfix_wang)
d <- as.data.frame(stack)
df <- d


# Fit linear regression model -------------------------
fit <- lm(nfix_wang ~ Diversity_Richness_ensembleMean_paBased, data = df)

# Extract R-squared value
r_squared <- summary(fit)$r.squared

# Extract coefficients of the linear regression model
intercept <- coef(fit)[1]
slope <- coef(fit)[2]

# Create formula text
formula_text <- paste("y =", round(intercept, 2), "+", round(slope, 2), "* x")

# Fit a linear regression model
model <- lm(nfix_wang ~ Diversity_Richness_ensembleMean_paBased, data = df)

# Extract the slope coefficient
slope <- coef(model)[2]  # 2nd coefficient corresponds to the slope


# Perform correlation significance test
cor_test_result <- cor.test(df$Diversity_Richness_ensembleMean_paBased, df$nfix_wang, method = "spearman", exact = FALSE, continuity = TRUE)


# Create bins
# df <- mutate(df, bin = cut_width(Diversity_Richness_ensembleMean_paBased, width = 0.05, boundary = 0) )


# Define the desired bin names
bin_names <- paste0("Bin ", 1:14)

# Convert the bin column to a factor with the desired levels
df <- df %>%
  mutate(bin = cut_width(Diversity_Richness_ensembleMean_paBased, width = 0.05, boundary = 0),
         bin_label = factor(bin, labels = bin_names))



# Remove NAs
df <- df[which(!is.na(df$Diversity_Richness_ensembleMean_paBased)), ]
# test <- aggregate(df, by = list(df$bin), FUN = "mean")
# test$T <- round(test$T, 0)
# df$mean_temp <- NA
# for(i in 1:nrow(test)){
#   df[which(df$bin == test[i, "Group.1"]), "mean_temp"] <- test[i, "T"]
# }


# 
sort(unique(df$bin))

# Boxplots: Total diazotroph community
gg_richness_total <- ggplot(data = df, aes(x = bin_label, y = nfix_wang)) +
  geom_boxplot(fill = orange_hex) +
  theme(
    panel.background = element_rect(fill = "white"),
    axis.text.x = element_text(
        margin = margin(t = 0, r = 0, b = 20, l = 0), # add some space between axis labels and axis ticks
        angle = 90, vjust = 0.5, 
        size = 40),
    axis.text.y = element_text(size = 70),
    # axis.text.x = element_text(size = 100),
    axis.title.x = element_text(size = 70),  # Adjust margin as needed
    axis.title.y = element_text(size = 70),
    legend.position = "none") +
    labs(
    x = expression(italic("Richness")),
    y = expression(italic(paste(N[2], " fixation flux (", "mmol N m"^{-2}, " y"^{-1}, ")")))) +
  annotate("text", x = as.numeric(sort(unique(df$bin_label)))[10], y = 400, 
           label = paste0("Spearman's rho = ", round(cor_test_result$estimate, 2)), 
           hjust = 1, vjust = 0, size = 12) + # Add Spearman rank correlation coefficient to the plot with increased size
  annotate("text", x =  as.numeric(sort(unique(df$bin_label)))[10], y = 380, 
           label = paste0("p-value < 0.001"), 
           hjust = 1, vjust = 0, size = 12) + # Add Spearman rank correlation coefficient to the plot with increased size
  annotate("text", x = as.numeric(sort(unique(df$bin_label)))[10], y = 360, 
           label = formula_text, 
           hjust = 1, vjust = 0, size = 12)





## BOXPLOTS RICHNESS vs. N2-Fixation
# Boxplot richness total diaz. community and nitrogen fixation rate Wang et al. 2019
sort(names(r))
stack <- stack( 
  r$Diversity_Richness_ensembleMean_paBased, 
  r$Diversity_Richness_ensembleMean_paBased_ncd,
  r$Diversity_Richness_ensembleMean_paBased_cyanos,
  r$nfix_wang)
d <- as.data.frame(stack)
df <- d


# Fit linear regression model -------------------------
fit <- lm(nfix_wang ~ Diversity_Richness_ensembleMean_paBased, data = df)

# Extract R-squared value
r_squared <- summary(fit)$r.squared

# Extract coefficients of the linear regression model
intercept <- coef(fit)[1]
slope <- coef(fit)[2]

# Create formula text
formula_text <- paste("y =", round(intercept, 2), "+", round(slope, 2), "* x")

# Fit a linear regression model
model <- lm(nfix_wang ~ Diversity_Richness_ensembleMean_paBased, data = df)

# Extract the slope coefficient
slope <- coef(model)[2]  # 2nd coefficient corresponds to the slope


# Perform correlation significance test
cor_test_result <- cor.test(df$Diversity_Richness_ensembleMean_paBased, df$nfix_wang, method = "spearman", exact = FALSE, continuity = TRUE)


# Create bins
# df <- mutate(df, bin = cut_width(Diversity_Richness_ensembleMean_paBased, width = 0.05, boundary = 0) )


# Define the desired bin names
bin_names <- paste0("Bin ", 1:14)

# Convert the bin column to a factor with the desired levels
df <- df %>%
  mutate(bin = cut_width(Diversity_Richness_ensembleMean_paBased, width = 0.05, boundary = 0),
         bin_label = factor(bin, labels = bin_names))



# Remove NAs
df <- df[which(!is.na(df$Diversity_Richness_ensembleMean_paBased)), ]
# test <- aggregate(df, by = list(df$bin), FUN = "mean")
# test$T <- round(test$T, 0)
# df$mean_temp <- NA
# for(i in 1:nrow(test)){
#   df[which(df$bin == test[i, "Group.1"]), "mean_temp"] <- test[i, "T"]
# }


# 
sort(unique(df$bin))

# Boxplots: Total diazotroph community
gg_richness_total <- ggplot(data = df, aes(x = bin_label, y = nfix_wang)) +
  geom_boxplot(fill = orange_hex) +
  theme(
    panel.background = element_rect(fill = "white"),
    axis.text.x = element_text(
        margin = margin(t = 0, r = 0, b = 20, l = 0), # add some space between axis labels and axis ticks
        angle = 90, vjust = 0.5, 
        size = 40),
    axis.text.y = element_text(size = 70),
    # axis.text.x = element_text(size = 100),
    axis.title.x = element_text(size = 70),  # Adjust margin as needed
    axis.title.y = element_text(size = 70),
    legend.position = "none") +
    labs(
    x = expression(italic("Richness")),
    y = expression(italic(paste(N[2], " fixation flux (", "mmol N m"^{-2}, " y"^{-1}, ")")))) +
  annotate("text", x = as.numeric(sort(unique(df$bin_label)))[10], y = 400, 
           label = paste0("Spearman's rho = ", round(cor_test_result$estimate, 2)), 
           hjust = 1, vjust = 0, size = 12) + # Add Spearman rank correlation coefficient to the plot with increased size
  annotate("text", x =  as.numeric(sort(unique(df$bin_label)))[10], y = 380, 
           label = paste0("p-value < 0.001"), 
           hjust = 1, vjust = 0, size = 12) + # Add Spearman rank correlation coefficient to the plot with increased size
  annotate("text", x = as.numeric(sort(unique(df$bin_label)))[10], y = 360, 
           label = formula_text, 
           hjust = 1, vjust = 0, size = 12)

# Boxplots: Cyanobacterial diazotroph community -------------------------
df <- d


# Define the desired bin names
bin_names <- paste0("Bin ", 1:15)

# Fit linear regression model -------------------------
fit <- lm(nfix_wang ~ Diversity_Richness_ensembleMean_paBased_cyanos, data = df)

# Extract R-squared value
r_squared <- summary(fit)$r.squared

# Extract coefficients of the linear regression model
intercept <- coef(fit)[1]
slope <- coef(fit)[2]

# Create formula text
formula_text <- paste("y =", round(intercept, 2), "+", round(slope, 2), "* x")

# Perform correlation significance test
cor_test_result <- cor.test(df$Diversity_Richness_ensembleMean_paBased_cyanos, df$nfix_wang, method = "spearman", exact = FALSE, continuity = TRUE)



# Convert the bin column to a factor with the desired levels
df <- df %>%
  mutate(bin = cut_width(Diversity_Richness_ensembleMean_paBased_cyanos, width = 0.05, boundary = 0),
         bin_label = factor(bin, labels = bin_names))

# Remove NAs
df <- df[which(!is.na(df$Diversity_Richness_ensembleMean_paBased_cyanos)), ]

# 
sort(unique(df$bin))

# Plot
gg_richness_cyano <- ggplot(data = df, aes(x = bin_label, y = nfix_wang)) +
  geom_boxplot(fill = dark_green_hex) +
  theme(
    panel.background = element_rect(fill = "white"),
    axis.text.x = element_text(
        margin = margin(t = 0, r = 0, b = 20, l = 0), # add some space between axis labels and axis ticks
        angle = 90, vjust = 0.5, 
        size = 40),
    axis.text.y = element_blank(),  # Remove y-axis tick labels
    axis.title.y = element_blank(), 
    axis.title.x = element_text(size = 70),  # Adjust margin as needed
    legend.position = "none") +
    labs(
    x = expression(italic("Richness"))) +
  annotate("text", x = as.numeric(sort(unique(df$bin_label)))[10], y = 400, 
           label = paste0("Spearman's rho = ", round(cor_test_result$estimate, 2)), 
           hjust = 1, vjust = 0, size = 12) + # Add Spearman rank correlation coefficient to the plot with increased size
  annotate("text", x = as.numeric(sort(unique(df$bin_label)))[10], y = 380, 
           label = paste0("p-value < 0.001"), 
           hjust = 1, vjust = 0, size = 12) + # Add Spearman rank correlation coefficient to the plot with increased size
  annotate("text", x = as.numeric(sort(unique(df$bin_label)))[10], y = 360, 
           label = formula_text, 
           hjust = 1, vjust = 0, size = 12)



# Spearman rank correlation coefficient
df_median <- aggregate(nfix_wang ~ bin, data = df, FUN = median)
# Calculate Spearman correlation
correlation <- cor.test(as.numeric(df_median$bin), df_median$nfix_wang, method = "spearman")
# Print the correlation coefficient and p-value
print(correlation)


# Boxplots: Non-cyanobacterial diazotroph community
df <- d

# Define the desired bin names
bin_names <- paste0("Bin ", 1:14)

# Fit linear regression model -------------------------
fit <- lm(nfix_wang ~ Diversity_Richness_ensembleMean_paBased_ncd, data = df)

# Extract R-squared value
r_squared <- summary(fit)$r.squared

# Extract coefficients of the linear regression model
intercept <- coef(fit)[1]
slope <- coef(fit)[2]

# Create formula text
formula_text <- paste("y =", round(intercept, 2), "+", round(slope, 2), "* x")

# Perform correlation significance test
cor_test_result <- cor.test(df$Diversity_Richness_ensembleMean_paBased_ncd, df$nfix_wang, method = "spearman", exact = FALSE, continuity = TRUE)


# Convert the bin column to a factor with the desired levels
df <- df %>%
  mutate(bin = cut_width(Diversity_Richness_ensembleMean_paBased_ncd, width = 0.05, boundary = 0),
         bin_label = factor(bin, labels = bin_names))

# Remove NAs
df <- df[which(!is.na(df$Diversity_Richness_ensembleMean_paBased_ncd)), ]

#
sort(unique(df$bin))

# Plot
gg_richness_ncd <- ggplot(data = df, aes(x = bin_label, y = nfix_wang)) +
  geom_boxplot(fill = dark_blue_hex) +
  theme(
    panel.background = element_rect(fill = "white"),
    axis.text.x = element_text(
        margin = margin(t = 0, r = 0, b = 20, l = 0), # add some space between axis labels and axis ticks
        angle = 90, vjust = 0.5, 
        size = 40),
    axis.text.y = element_blank(),
    # axis.text.x = element_text(size = 100),
    axis.title.x = element_text(size = 70),  # Adjust margin as needed
    axis.title.y = element_blank(),
    legend.position = "none") +
    labs(
    x = expression(italic("Richness")),
    y = expression(italic(paste(N[2], " fixation flux (", "mmol N m"^{-2}, " y"^{-1}, ")")))) +
  annotate("text", x = as.numeric(sort(unique(df$bin_label)))[10], y = 400, 
           label = paste0("Spearman's rho = ", round(cor_test_result$estimate, 2)), 
           hjust = 1, vjust = 0, size = 12) + # Add Spearman rank correlation coefficient to the plot with increased size
  annotate("text", x = as.numeric(sort(unique(df$bin_label)))[10], y = 380, 
           label = paste0("p-value < 0.001"), 
           hjust = 1, vjust = 0, size = 12) + # Add Spearman rank correlation coefficient to the plot with increased size
  annotate("text", x = as.numeric(sort(unique(df$bin_label)))[10], y = 360, 
           label = formula_text, 
           hjust = 1, vjust = 0, size = 12)



# Gridd arrange
gg_arranged <- gridExtra::grid.arrange(
  gg_richness_total,
  gg_richness_cyano,
  gg_richness_ncd,
  ncol = 3,
  nrow = 1
)


# Save plot
fln <- paste0( "/home/deriksson/Projects/Notion_DE/Code/Tables/Visualizations/Boxplot_BEF_nitrogenFixation_Richness_allGroups.png" )
ggsave(
  filename = fln,
  plot = gg_arranged,
  width = 35,
  height = 20,
  dpi = 300
)
