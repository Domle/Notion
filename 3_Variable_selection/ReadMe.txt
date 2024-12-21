ReadMe file for section 3 (3_Variable_selection)

This folder focuses on the analysis of environmental data, primarily for understanding the correlations between various environmental parameters and their relationship with diazotroph distributions using statistical and machine learning models. The scripts contain several stages including correlation analysis, univariate model tests for variable ranking, species-specific variable selection, and the creation of indices for predictors used for each diazotroph in the final species distribution model.


1. Correlation Analysis Among Environmental Variables
To assess the correlation between environmental parameters using Pearson and Spearman rank correlation coefficients, and save the results in a table. Use the Spearman correlation table to ensure that selected variables have a correlation less than 0.7, avoiding highly correlated variables within the ensemble.

2. Univariate Model Tests for Variable Ranking
Rank environmental variables for each species using three different algorithms: General Linear Model (GLM), Generalized Additive Model (GAM), and Random Forest (RF). Rank the variables based on model performance metrics (e.g., D-squared, R-squared, or out-of-bag error).

3. Species-Specific Variable Ranking and Selection
To create tables that rank environmental variables based on their relevance to specific species, taking into account the uncertainties arising from different background selection strategies and algorithms.
We use the average of the rankings to obtain a final variable ranking per species.


Authors:
Dominic Eriksson
Environmental Physics Group, UP, ETH Zurich, Switzerland
Email: deriksson@ethz.ch
21st of December 2024

