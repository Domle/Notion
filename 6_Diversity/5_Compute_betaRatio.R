## This script calculates the beta ratio based on the nestedness (beta-ratio = nestedness/jaccard)
## We calculate for each enemble member and month and save it.

## Input files:
##      1.  CSV files containing the betadiversity components and location.

## Output files:
##      2. CSV files containing the calculated beta-ratio.

## Strategy: We open each file, calcluate beta-ration (nestedness/jaccard) and save the files again.


## Author:  Dominic Eriksson
##          Environmental Physics Group, UP
##          ETH Zürich, Zurich
##          Switzerland

## deriksson@ethz.ch; 9th of August 2023 ---------------------------------------------------------------------

# Clear workspace
rm(list = ls())

# Libraries
library(raster)

# Directories
wd_in <- "/home/deriksson/Projects/Notion_DE/Code/8_Diversity/5_Output/"
wd_out <- "/net/kryo/work/deriksson/Projects/Notion_DE/Code/8_Diversity/8_Output/"


# Vectors
filenames <- list.files(wd_in)

# Calculate betaRatio using loop
for(f in seq_along(filenames)){

    # Print progress
    print( paste0("File number ", f, ". Out of ", length(filenames), ".") )

    # Open data
    df <- read.csv( paste0(wd_in, filenames[f]) )


    # Calculate beta ratio (nestedness/jaccard) --> value higher than 0.5 indicates region dominated by nestedness and below 0.5 species turnover is dominating
    df$betaRatio_jne <- df$jne/df$jac

    # Save
    fln <- paste0(wd_out, gsub("Betadiversity", "BetaRatio_jne", filenames[f]))
    write.csv(df, fln, row.names = FALSE) 
} # close loop


###==============================================================
### END
###==============================================================