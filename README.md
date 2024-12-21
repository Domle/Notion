ReadMe for the manuscript:

Global Marine Diazotroph Diversity and Nitrogen Fixation Rates

This repository contains the scripts to reproduce the results from Eriksson et al. (in submission, Preprint: https://doi.org/10.32942/X2Z323). 


Overview

This project investigates the global distribution and diversity of marine diazotrophs (organisms that fix atmospheric nitrogen), exploring how their diversity correlates with nitrogen fixation rates. By integrating microscopic and sequence-based data (e.g.: qPCR, metagenomics) with machine learning techniques, the study provides a global-scale analysis of diazotroph diversity across 15 species including heterotrophic non-cyanobacterial diazotrophs as well as cyanobacterial ones. The project also links biodiversity patterns to nitrogen fixation rates, shedding light on the role of microbial diversity in global biogeochemical cycles.

The work was carried out using ensemble-based species distribution models (SDMs) that integrate data from various sources, such as microscopy, qPCR, and metagenomics, to estimate global diazotroph richness. Beta diversity indices (species turnover and nestedness) were computed to assess spatial patterns in community composition and how nestedness and species turnover structure the global diazotroph richness gradient. The project aims to emphasize the importance of microbial biodiversity in ecosystem function.

Key Questions Addressed:

What are the global patterns of pelagic diazotroph richness, and what major environmental drivers influence these patterns?
How do these patterns covary with marine biological nitrogen fixation (BNF)?
Do cyanobacterial and non-cyanobacterial diazotrophs exhibit distinct covariance patterns with annual BNF, and which group is more strongly linked to BNF?

Project Structure

This repository is organized into several sections, each focusing on different aspects of the analysis pipeline. Each section contains a ReadMe with an overview of what is done across the script. Additionally, each script contains a more detailed opening paragraph on what exactly is done within the script and additional commenting is used to increase clarity of each analysis step. 

The main results of this project are stored in raster files (.grd) containing beta diversity and richness indices for each grid cell. These results are further aggregated and averaged to show global patterns and temporal dynamics of diazotrophs in the oceans. 

This project utilizes R and several R libraries for data manipulation, analysis, and visualization.

The diazotroph dataset can be accessed via https://doi.org/10.3929/ethz-b-000635803. 
Environmental predictor sources are mentioned in the Material and Methods part. The environmental predictor set, if not downloaded by yourself, can be requested via deriksson@ethz.ch.

Dominic Eriksson
Environmental Physics Group, UP
ETH Zurich, Switzerland
Email: deriksson@ethz.ch, 21st of December 2024
