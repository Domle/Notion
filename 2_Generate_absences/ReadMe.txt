ReadMe file for section 2 (2_Generate_absences)

Generate Pseudoabsences

This folder contains scripts for generating pseudoabsences using various 
target-group approaches, ensuring environmentally stratified designs for 
modeling diazotroph presences and pseudoabsences. 
These methods balance environmental bias and ensure robust data 
preparation for subsequent analyses.

General Strategy
Use the target-group approach (Phillips et al., 2009) to generate 
pseudoabsences for taxa with ≥ 14 observations.
Absences are proportionally distributed relative to presences within 
environmental strata of species from a defined target group. 

Total Target Group:
Includes taxonomic records from Phytobase (Righetti et al., 2019) 
within the surface ocean mixed-layer.
Excludes diatoms and dinoflagellates to reduce bias caused by their 
larger size and sampling imbalance.

Group-Specific Target Group:
Uses all locations from the compiled diazotroph database in this study.
Provides a more focused background by considering only diazotroph 
observations.

Cruise-Specific Target Group:
Generates background data from specific cruises using the same sampling methods.
Ensures consistency in methodology and environmental context for background selection.


Dominic Eriksson
Environmental Physics Group, UP
ETH Zurich, Switzerland
Email: deriksson@ethz.ch; 21st of December 2024
