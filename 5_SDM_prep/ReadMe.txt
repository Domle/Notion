ReadMe file for section 5 (5_SDM_prep)

This section of the pipeline focuses on processing Species Distribution Model (SDM) projections into monthly and annual temporal resolutions for each taxon, and then merging these projections to create ensemble-based biogeographies. The resulting outputs provide a comprehensive view of species distributions, with both spatial and temporal components.

1. Merging SDM Projections into Monthly Slots
Merging SDM projections for each taxon and ensemble member into monthly raster stacks, representing the species' predicted distributions across the year (January to December).An evaluation file containing the TSS score statistic for each taxon and ensemble member.
SDM output for each taxon and ensemble member (either HSI or P/A converted projections).


A file titled Taxon_all_models.rds: A list object for each taxon, with a length of 90 (representing each ensemble member). Each list element is a raster stack of length 12, corresponding to the months of the year (January to December).
Steps:

Data Preparation:
The SDM outputs for each taxon and ensemble member are loaded into the script. These outputs represent the projections for each species across different time periods and model configurations.
Monthly Projections:
The SDM projections are merged into monthly slots (1–12) to create raster stacks for each ensemble member. This enables a detailed analysis of the temporal patterns in species distributions.

2. Ensemble Biogeographies: Monthly and Annual Resolutions
Objective:
To calculate the global annual distribution of each species using ensemble-based SDM projections and compute the ensemble mean and standard deviation across all ensemble members for both monthly and annual time scales. This provides a quantitative assessment of the uncertainty and variability in the species' predicted distribution. The mean and standard deviation are calculated across the 90 ensemble members, providing insight into the consensus of the projections and the associated uncertainty.

Authors:
Dominic Eriksson
Environmental Physics Group, UP, ETH Zurich, Switzerland
Email: deriksson@ethz.ch, 21st of December 2024
