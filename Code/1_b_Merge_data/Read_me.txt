This folder contains contains the script for the merging of 
the data sets that are stored in the folder "1_Downloeaded_Data".

The R file "1_Merged_Data" contains the scripts to merge the
following data sets:
- Gradoville et al - 2020 - Latitudinal constraints on the abundance and activity of the cyanobacterium UCYN-A and other marine diazotrophs in the North Pacific
- Luo et al. - 2012 - Database of diazotrophs in global ocean: Abundance, biomass and nitrogen fixation rates
- Righetti et al - 2020 - Phytobase, a global synthesis of open-ocean phytoplankton occurrences
- Weiyi Tang and Nicolas Cassar - 2019 - Data-Driven Modeling of the Distribution of Diazotrophs in the Global Ocean
- GBIF (https://www.gbif.org/) - Data accesses on on 20 May, 2020
- OBIS (https://obis.org/) - Data accessed on 20 May, 2020 
- Tara oceans (via Mgnify) - OTU based data

Note that at this stage no stringent cleaning of dates and/or quality control has been implemented. 

The parameters retained from the data sets are:

 [1] "decimalLongitude"                       "decimalLatitude"                       
 [3] "year"                                   "month"                                 
 [5] "day"                                    "depth"                                 
 [7] "depthAccuracy"                          "depthIntegral"                         
 [9] "scientificName"                         "taxonRank"                             
[11] "phylum"                                 "class"                                 
[13] "genus"                                  "species"                               
[15] "lifeForm"                               "associatedTaxa"                        
[17] "occurrenceStatus"                       "trichomeStatus"                        
[19] "trichomeQuantity"                       "trichomeQuantityType"                  
[21] "cellStatus"                             "cellQuantity"                          
[23] "cellQuantityType"                       "nifStatus"                             
[25] "nifQuantity"                            "nifQuantityType"                       
[27] "hitStatus"                              "hitQuantity"                           
[29] "hitQuantityType"                        "biomass_conversion_factor"             
[31] "biomass_conversion_factor_QuantityType" "biomassQuantity"                       
[33] "biomassQuantityType"                    "basisOfRecord"                         
[35] "measurementMethod"                      "references"                            
[37] "occurrenceID"                           "sampleID"                              
[39] "stationID"                              "cruiseID"                              
[41] "Sea_surface_temperature_C"              "Temperature_C"                         
[43] "Sea_surface_Salinity"                   "Salinity"                              
[45] "Sea_surface_nitrate_micromolar"         "Nitrate_micromolar"                    
[47] "Sea_surface_phosphate_micromolar"       "Phosphate_micromolar"                  
[49] "Sea_surface_Fe_nanomolar"               "Fe_nanomolar"                          
[51] "ChlorophyllQuantity"                    "ChlorophyllQuantityType"               
[53] "flag"				      "set"

Datasets were merged. To remove duplicated records we use an occurrenceID with the information:
- occurrenceID: "scientificName"_"family"_"genus"_"decimalLongitude"_"decimalLatitude"_"year"_"month"_"day"_"depth"

--- Dominic Eriksson and Damiano Righetti, 03 June 2020 ---
