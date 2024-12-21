ReadMe file for section 6 (6_Diversity)

This section of the pipeline focuses on calculating biodiversity metrics: normalized diazotroph richness and betadiversity (Jaccard dissimilarity index, species turnover, and nestedness). The analysis is conducted for each ensemble member and averaged at the end across all ensemble members.

1. Normalized Annual Diazotroph Richness
We interpret the annual presence of a diazotroph as being present at least once within the twelve months. The richness is then normalized by the total number of taxa modeled (n = 15 for total diazotroph community).

Data Interpretation:
Each taxon’s monthly SDM outputs are interpreted as annual presences based on whether the species is present at least once within the twelve months.

Richness Calculation:
The richness is calculated by summing the taxa within each grid cell, normalized by the total number of taxa modeled.

2. Betadiversity Calculation (Jaccard Dissimilarity Index)
To calculate the beta diversity for each ensemble member using the Jaccard dissimilarity index, species turnover, and nestedness. The beta diversity components are computed for each location, and a global average is then calculated for each grid cell, reflecting the spatial variability of each grid cell to the global ocean. Temporal changes in community composition are also analyzed at each location across the twelve months to see the magnitude of community changes in each grid cell across the twelve months.

Betadiversity Calculation:
The betapart package is used to calculate the Jaccard dissimilarity index, species turnover, and nestedness for each location across all ensemble members. This involves evaluating community composition in each grid cell and calculating the spatial variability across the entire ocean.

Authors:
Dominic Eriksson
Environmental Physics Group, UP, ETH Zurich, Switzerland
Email: deriksson@ethz.ch, 21st of December 2024
