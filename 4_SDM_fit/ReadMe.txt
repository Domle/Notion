ReadMe file for section 4 (4_SDM_fit)

This folder is designed for species distribution modeling (SDM) using 
a variety of statistical and machine learning models, including Generalized 
Linear Models (GLM), Generalized Additive Models (GAM), and Random Forests 
(RF). The aim is to model species distributions based on environmental 
data and evaluate the model’s performance using various evaluation metrics 
and cross-validation strategies. The pipeline includes steps for preparing 
data, fitting models, testing predictive accuracy, and visualizing the 
results.

1. Species Distribution Modeling (GLM, GAM, RF)
To model the distribution of species using environmental variables 
and presence/absence (P/A) data with the help of three modeling 
algorithms: GLM, GAM, and RF. The species distribution models are 
evaluated through various metrics including deviance explained, 
adjusted R-squared, and predictive accuracy using cross-validation.

Data Preparation:
The pipeline begins by loading a list of species presence/absence data 
for each taxon.

Model Fitting:
It then fits three types of models (GLM, GAM, and RF) to the species 
data, using the environmental data as predictors using the predictors 
created in the folder before. These models are projected onto 
monthly climatological global conditions for each species and each 
ensemble member on a 1° longitude by 1° latitude grid.

Model Evaluation:
The performance of each model is evaluated using metrics such as 
deviance explained, adjusted R-squared, and TSS (True Skill Statistic). A 
cross-validation approach (n-fold) is employed to test the predictive 
skill of the models, with various sample splits.

2 & 3. Species Distribution Model Evaluation (Boxplot Visualization)
To visualize the TSS (True Skill Statistic) scores for each species 
across different predictor ensembles using boxplots. This provides 
a comparative overview of model performance for different taxonomic 
groups.

Authors:
Dominic Eriksson
Environmental Physics Group, UP, ETH Zurich, Switzerland
Email: deriksson@ethz.ch, 21st of December 2024
