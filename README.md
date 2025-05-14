# README for the Manuscript  
**Global Marine Diazotroph Diversity and Nitrogen Fixation Rates**

This repository contains all scripts required to reproduce the results from:  
**Eriksson et al.** (in revision)  
**Preprint:** [https://doi.org/10.32942/X2Z323](https://doi.org/10.32942/X2Z323)

---

## 📘 Overview

This project investigates the global distribution and diversity of marine diazotrophs—organisms that fix atmospheric nitrogen—and how this diversity correlates with nitrogen fixation rates across the global ocean.

By integrating microscopy and sequence-based data (e.g., qPCR, metagenomics) with ensemble machine learning techniques, this study assesses the global-scale diversity of 15 key marine diazotroph species, including both heterotrophic non-cyanobacterial and cyanobacterial taxa. Key outputs include beta diversity indices, species richness estimates, and spatial/temporal models linking diversity to marine nitrogen fixation (BNF).

This analysis emphasizes the ecological significance of microbial biodiversity in shaping biogeochemical processes.

### Key Research Questions
- What are the global patterns of pelagic diazotroph richness, and what major environmental drivers influence these patterns?
- How do diazotroph richness and composition relate to spatial patterns of biological nitrogen fixation (BNF)?
- Are cyanobacterial and non-cyanobacterial diazotrophs differentially associated with global BNF?

---

## 📂 Project Structure

The repository is organized into subfolders representing key parts of the analysis pipeline. Each R script begins with a structured header describing the purpose, input/output files, and implementation strategy. In addition:

- All scripts are formatted following good coding practices.
- In-line comments improve clarity and reproducibility.
- Outputs are typically `.csv`, `.rds`, or raster files (`.grd`) showing richness or beta diversity.

**Main results** are stored in raster files representing:
- Richness indices (e.g., total, group-specific)
- Beta diversity (species turnover and nestedness)

These files are later used to generate global distribution maps and correlation plots with BNF.

---

## 🖥️ 1. System Requirements

### Software Dependencies
- **R** version ≥ 4.2.2 (2022-10-31)

- **Required R Packages**
  - `classInt` ≥ 0.4.10
  - `doParallel` ≥ 1.0.17
  - `parallel` ≥ 4.2.2
  - `raster` ≥ 3.6.26
  - `dismo` ≥ 1.3.14
  - `ggplot2` ≥ 3.5.1
  - `dplyr` ≥ 1.1.2
  - `gridExtra` ≥ 2.3
  - `sf` ≥ 1.0.14
  - `rnaturalearth` ≥ 1.0.1
  - `psych` ≥ 2.2.9
  - `tidyverse` ≥ 1.3.2
  - `doBy` ≥ 4.6.22
  - `mgcv` ≥ 1.8.41
  - `randomForest` ≥ 4.7.1.1
  - `tidyr` ≥ 1.3.1
  - `rJava` ≥ 1.0.6
  - `PresenceAbsence` ≥ 1.1.11
  - `colorRamps` ≥ 2.3.4
  - `RColorBrewer` ≥ 1.1.3
  - `visreg` ≥ 2.7.0
  - `patchwork` ≥ 1.3.0
  - `betapart` ≥ 1.6
  - `ggpubr` ≥ 0.6.0
  - `terra` ≥ 1.7.78
  - `tidyterra` ≥ 0.4.0
  - `hexbin` ≥ 1.28.2
  - `utils` ≥ 4.2.2

### Operating Systems Tested
macOS Sequoia 15.5

⚙️ 2. **Installation Guide**

## Installation Steps
1. Install **R** from [https://cran.r-project.org](https://cran.r-project.org)
2. Install required packages:
   ```r
   install.packages(c("tidyverse", "raster", "data.table", "caret", "sf", "sp", "terra", "ggplot2", "dplyr"))
