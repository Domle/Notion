### In this script we plot main figure 1 from the Manuscript "Nitrogen fixation increases with diazotroph richness in the global ocean"

### NOTE: 	When running the code, please always double check input and output directories, to ensure you know where exactly files get saved so they
###			    can be referred as the correct input files within each consecutive script.

# Author: 	
#     Dominic Eriksson
#			Environmental Physics Group, UP
# 		ETH Zurich
#			Switzerland


# Input files:
#       1. shapefile goas

# Output files:
#       1. Raster stack as .grd file, later used for final plotting

# derikssong@ethz.ch; 18th of June 2024




#################################################################################################
## LINEPLOTS: RESPONSE CURVES
#################################################################################################
d <- read.csv("/home/deriksson/Projects/Notion_DE/Code/6_SDM_fit/1_Output/Dataframe_ResponseCurves.csv")
unique(d$taxon)

# Subset for taxa of interest
## Choose between cyanos and NCDs
# Subset for species of interest
vector_taxaNames <- c(
  "Calothrix",  
  "Cyanothece",  
  "Richelia",
  "Richelia intracellularis",
  "Trichodesmium erythraeum",
  "Trichodesmium thiebautii",
  "UCYN.A1",
  "UCYN.A2",
  "UCYN.B",
  "Gamma.A",
  "HBD-02",
  "HBD-03",
  "HBD-04",
  "HBD-05",
  "HBD-06")
#
dd <- d[d$taxon %in% vector_taxaNames, ]

# Which predictors have been used by all of the models
pred <- unique(dd$predictor)
for(i in seq_along(pred)){
  nr <- dd[which(dd$predictor == pred[i]), "taxon"]
  nr <- length(unique(nr))
  print(nr)
}

# Subset for predictor
dd <- d[which(d$predictor == "T"), ]
dd <- d[which(d$algorithm == "gam"), ]
# dd <- dd[which(dd$taxon == "Trichodesmium"), ]

# Round
# test <- dd

# test1 <- test[which(test$ensembleID == "pres_abs,cr_bg_nonov"), ]
# sort(unique(test1$Predictor_values))
# test2 <- test[which(test$ensembleID == "pres_abs,gr_bg_nonov"), ] 
# sort( unique(test2$Predictor_values) )

# Round predictor values
dd$Predictor_values <- round( dd$Predictor_values, 2 )


vector_taxaNames <- c("HBD-02", "HBD-03", "HBD-04", "HBD-05", "HBD-06")
ddd <- dd[dd$taxon %in% vector_taxaNames, "pres_abs,cr_bg_nonov"]

# Aggregate on predictor values
ddd <- aggregate(dd$visregFit, by = list(dd$Predictor_values), FUN = mean )
names(ddd) <- c("Predictor_values", "visregFit")

# Plot
ggplot( data = ddd, aes(x = Predictor_values, y = visregFit, color = taxon) ) +
  geom_smooth()



## Choose between cyanos and NCDs
vector_taxaNames <- c(  
  "Calothrix",  
  "Cyanothece",  
  "Richelia",
  "Richelia intracellularis",
  "TRichodesmium eryhtraeum",
  "Trichodesmium thiebautii",
  "UCYN.A1",
  "UCYN.A2",
  "UCYN.B")

vector_taxaNames <- c(  
  "Gamma.A",
  "HBD-02",
  "HBD-03",
  "HBD-04",
  "HBD-05",
  "HBD-06")

dd <- d[d$taxon %in% vector_taxaNames, ]
dd <- dd[which(dd$algorithm == "gam"), ]
# Parameter vector
# vector_parameters <- unique(dd$predictor)
vector_parameters <- c("T", "P", "logP", "N", "logN", "Chl")
vector_parameters <- c("T", "P", "N", "Chl")
vector_xlab_labels <- c("Temperature [° Celsius]", "Phosphate [micromolar]", 
  "Nitrate [mircomolar]", "Chlorophyll a [microgram/liter]")
# 
vector_parameters <- c("T", "dN_dt", "Nstar", "dP_dt")
vector_xlab_labels <- c("Temperature [° Celsius]", "Temporal trend nitrate [micromolar]", 
  "Excess concetration of nitrate [micromolar]", "Temporal trend phosphate [micromolar]")
# Plot tags
# vector_tags <- c("A", "B", "C", "D", "E", "F")
# vector_tags <- c("G", "H", "I", "J", "K", "L")

#
list_plots <- list() 
for( p in  seq_along(vector_parameters)){

  # Subset for taxa and parameters
  ddd <- dd[which(dd$predictor == vector_parameters[p]), ]

  # Vector 
  vector_taxaNames <- unique(ddd$taxon)
  #
  l <- list()
  for(v in seq_along(vector_taxaNames)){

    # Subset for each taxa
    dddd <- ddd[which(ddd$taxon == vector_taxaNames[v]), ]

    # # Aggregate on predictor values
    # ddd$Predictor_values <- round(dddd$Predictor_values, 0)
    dddd <- aggregate(dddd$visregFit, by = list(dddd$Predictor_values), FUN = mean )
    names(dddd) <- c("Predictor_values", "visregFit")
    dddd$Taxon <- vector_taxaNames[v]

    # Add to list
    l[[v]] <- dddd
  } # end of loop across taxa

  # Merge list to dataframe
  df <- do.call("rbind", l)

  # Plot ensemble response curves
  # Plot
 list_plots[[p]] <- ggplot( data = df, aes(x = Predictor_values, y = visregFit, color = Taxon) ) +
    geom_smooth(se = FALSE) + ggtitle( "" ) +
    ylab("Species presence") + xlab( vector_xlab_labels[p] ) +
    theme(axis.text.x = element_text(size = 32),
          axis.text.y = element_text(size = 32),
          axis.title.x = element_text(size = 32),
          axis.title.y = element_text(size = 32),
          legend.title = element_text(size = 32),      # Increase legend title font size
          legend.text = element_text(size = 28), legend.position = "top") 

} # end of loop across parameters

# Arrange the plots
arranged_plot <- do.call(ggpubr::ggarrange, c(list_plots, nrow = 3, ncol = 2))


# Save plot
fln <- "/home/deriksson/Projects/Notion_DE/Code/Tables/Visualizations/Lineplots_ResponseCurves_nonCyano.png"
ggsave( 
  plot = arranged_plot,
  filename = fln,
  dpi = 300,
  width = 20,
  height = 20
)

