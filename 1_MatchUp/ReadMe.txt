ReadMe file for section 1 (1_MatchUp)

Script 1:
Purpose:
This script processes a compiled dataset of diazotroph observations, 
creating both raw ungridded and gridded versions. It performs this procedure three times: once with the total dataset and twice by separating observations based on taxonomic assignment (microscopy-based or sequence-based identification).

Steps:

- Data Loading and Adjustment:
    Reads the dataset and adjusts column names to ensure compatibility for later analysis.
    Rounds coordinates to fit a specific grid resolution (-179.5 to 179.5 longitude and -89.5 to 89.5 latitude).

- Creating Raw Dataset:
    Creates a raw dataset by adjusting and preparing the observations for subsequent analysis.

- Creating Gridded Datasets:
    Aggregates data using longitude, latitude, and month to create gridded datasets.
    Discards year as a factor to aggregate cells at a 1° x 1° resolution, consolidating different years but the same month into one cell.

Script 2:
Purpose:
This script matches diazotroph observations with environmental variables and prepares outputs in both gridded and ungridded/raw formats.

Steps:

- Data Preparation:
    Loads environmental parameters and diazotroph observations from CSV files formatted according to Darwin Core standards.

- Conversion and Matching:
    Converts the diazotroph dataframe into a spatial object.
    Matches diazotroph presences with environmental parameters at a monthly resolution.
    Subsets the dataset to focus on oceanic conditions by removing coastal observations, those with depths < 200 meters, and salinity less than 20.

- Data Processing:
    Re-orders the dataset based on the number of observations (from highest to lowest).
    Saves the processed dataset as CSV files following Darwin Core standards, ready for further analysis or publication.

Script 3:
Purpose: Visualize outputs from previous data processing steps and validate data integrity. Specifically, for Script 2, visualize diazotroph presences in temperature-nitrate space.

Steps:


- 1_Output:
    Directory Setup: Define the input directory for data files.
    Read Files: Load CSV files containing processed data.
    Plotting: Loop through each file, convert data to spatial objects, and plot diazotroph presences and absences on a pacific-centered world map. Each plot represents different datasets or stages of processing from Script 1.

- 2_Output:
    Directory Setup: Define the input directory for another set of output data files.
    Read Files: Load CSV files containing more processed data.
    Plotting: Loop through each file, extract relevant columns, and plot temperature-nitrate space for diazotroph presences. Each plot represents different datasets or stages of processing from Script 2.



Author:
Dominic Eriksson
Environmental Physics Group, UP
ETH Zurich, Switzerland
Email: deriksson@ethz.ch
Date: 19th of June, 2024 ------------------------------------------------------------------------
