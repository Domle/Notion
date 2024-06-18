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
## BINNED BOXPLOT CORRELATION NITROGEN FIXATION AND RICHNESS - ACROSS DIFFERENT BIN SIZES
#################################################################################################
# Colors used for the different gourps

# Load Table
r <- brick("/home/deriksson/Projects/Notion_DE/Code/Tables/1_Output/FinalTable_AllMetrics.grd")



## BOXPLOTS RICHNESS vs. N2-Fixation
# Boxplot richness total diaz. community and nitrogen fixation rate Wang et al. 2019
sort(names(r))
stack <- stack( 
  r$Diversity_Richness_ensembleMean_paBased, 
  r$Diversity_Richness_ensembleMean_paBased_ncd,
  r$Diversity_Richness_ensembleMean_paBased_cyanos,
  r$nfix_wang)
d <- as.data.frame(stack)

binSizes <- seq(0.05, 0.8, by = 0.05)[1:4] # We only need the first four, upper ones binn in one or bins remain two
list_plots <- list()
for(b in seq_along(binSizes)){

  df <- d
  # Create bins
  df <- mutate(df, bin = cut_width(Diversity_Richness_ensembleMean_paBased, width = binSizes[b], boundary = 0) )

  # Remove NAs
  df <- df[which(!is.na(df$Diversity_Richness_ensembleMean_paBased)), ]
  # test <- aggregate(df, by = list(df$bin), FUN = "mean")
  # test$T <- round(test$T, 0)
  # df$mean_temp <- NA
  # for(i in 1:nrow(test)){
  #   df[which(df$bin == test[i, "Group.1"]), "mean_temp"] <- test[i, "T"]
  # }


  # Boxplots: Total diazotroph community
  if(b == 1){
  gg_richness_total <- ggplot(data = df, aes(x = bin, y = nfix_wang)) +
    geom_boxplot(fill = orange_hex) +
    theme(
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(
          margin = margin(t = 0, r = 0, b = 20, l = 0), # add some space between axis labels and axis ticks
          angle = 45, vjust = 0.5, 
          size = 18),
      axis.text.y = element_text(size = 18),
      axis.title.y = element_text(size = 50),
      axis.title.x = element_text(size = 50),
      legend.position = "none") +
      labs(
      x = expression(italic("Richness")),
      y = expression(italic(paste(N[2], " fixation flux (", "mmol N m"^{-2}, " y"^{-1}, ")")))
    )
  }else{
    gg_richness_total <- ggplot(data = df, aes(x = bin, y = nfix_wang)) +
    geom_boxplot(fill = orange_hex) +
    theme(
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(
          margin = margin(t = 0, r = 0, b = 20, l = 0), # add some space between axis labels and axis ticks
          angle = 45, vjust = 0.5, 
          size = 18),
      axis.text.y = element_text(size = 18),
      axis.title.y = element_text(size = 50),
      axis.title.x = element_text(size = 50),
      legend.position = "none") +
      labs(
      x = expression(italic("Richness")),
      y = ""
    )
  }

  # Save plot in list
  list_plots[[b]] <- gg_richness_total

}

multi_plot <- cowplot::plot_grid(
  plotlist = list_plots,
  ncol = 4
)

# Save 
fln <- "/home/deriksson/Projects/Notion_DE/Code/Tables/Visualizations/Boxplot_BEF_nitrogenFixation_Richness_AcrossBinSizes_totalDiaz.png"
ggsave(
  filename = fln,
  plot = multi_plot,
  height = 15,
  width = 25,
  dpi = 300
)


### Cyanos -------------------------------------------------------
list_plots <- list()
for(b in seq_along(binSizes)){

  df <- d
  # Create bins
  df <- mutate(df, bin = cut_width(Diversity_Richness_ensembleMean_paBased_cyanos, width = binSizes[b], boundary = 0) )

  # Remove NAs
  df <- df[which(!is.na(df$Diversity_Richness_ensembleMean_paBased_cyanos)), ]


  # Boxplots: Total diazotroph community
  if(b == 1){
  gg_richness_total <- ggplot(data = df, aes(x = bin, y = nfix_wang)) +
    geom_boxplot(fill = dark_green_hex) +
    theme(
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(
          margin = margin(t = 0, r = 0, b = 20, l = 0), # add some space between axis labels and axis ticks
          angle = 45, vjust = 0.5, 
          size = 18),
      axis.text.y = element_text(size = 18),
      axis.title.y = element_text(size = 50),
      axis.title.x = element_text(size = 50),
      legend.position = "none") +
      labs(
      x = expression(italic("Richness")),
      y = expression(italic(paste(N[2], " fixation flux (", "mmol N m"^{-2}, " y"^{-1}, ")")))
    )
  }else{
    gg_richness_total <- ggplot(data = df, aes(x = bin, y = nfix_wang)) +
    geom_boxplot(fill = dark_green_hex) +
    theme(
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(
          margin = margin(t = 0, r = 0, b = 20, l = 0), # add some space between axis labels and axis ticks
          angle = 45, vjust = 0.5, 
          size = 18),
      axis.text.y = element_text(size = 18),
      axis.title.y = element_text(size = 50),
      axis.title.x = element_text(size = 50),
      legend.position = "none") +
      labs(
      x = expression(italic("Richness")),
      y = ""
    )
  }

  # Save plot in list
  list_plots[[b]] <- gg_richness_total

}

multi_plot <- cowplot::plot_grid(
  plotlist = list_plots,
  ncol = 4
)

# Save 
fln <- "/home/deriksson/Projects/Notion_DE/Code/Tables/Visualizations/Boxplot_BEF_nitrogenFixation_Richness_AcrossBinSizes_cyanoDiaz.png"
ggsave(
  filename = fln,
  plot = multi_plot,
  height = 15,
  width = 25,
  dpi = 300
)

### Non-cyanos -------------------------------------------------------
list_plots <- list()
for(b in seq_along(binSizes)){

  df <- d
  # Create bins
  df <- mutate(df, bin = cut_width(Diversity_Richness_ensembleMean_paBased_ncd, width = binSizes[b], boundary = 0) )

  # Remove NAs
  df <- df[which(!is.na(df$Diversity_Richness_ensembleMean_paBased_ncd)), ]


  # Boxplots: Total diazotroph community
  if(b == 1){
  gg_richness_total <- ggplot(data = df, aes(x = bin, y = nfix_wang)) +
    geom_boxplot(fill = dark_blue_hex) +
    theme(
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(
          margin = margin(t = 0, r = 0, b = 20, l = 0), # add some space between axis labels and axis ticks
          angle = 45, vjust = 0.5, 
          size = 18),
      axis.text.y = element_text(size = 18),
      axis.title.y = element_text(size = 50),
      axis.title.x = element_text(size = 50),
      legend.position = "none") +
      labs(
      x = expression(italic("Richness")),
      y = expression(italic(paste(N[2], " fixation flux (", "mmol N m"^{-2}, " y"^{-1}, ")")))
    )
  }else{
    gg_richness_total <- ggplot(data = df, aes(x = bin, y = nfix_wang)) +
    geom_boxplot(fill = dark_blue_hex) +
    theme(
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(
          margin = margin(t = 0, r = 0, b = 20, l = 0), # add some space between axis labels and axis ticks
          angle = 45, vjust = 0.5, 
          size = 18),
      axis.text.y = element_text(size = 18),
      axis.title.y = element_text(size = 50),
      axis.title.x = element_text(size = 50),
      legend.position = "none") +
      labs(
      x = expression(italic("Richness")),
      y = ""
    )
  }

  # Save plot in list
  list_plots[[b]] <- gg_richness_total

}

multi_plot <- cowplot::plot_grid(
  plotlist = list_plots,
  ncol = 4
)

# Save 
fln <- "/home/deriksson/Projects/Notion_DE/Code/Tables/Visualizations/Boxplot_BEF_nitrogenFixation_Richness_AcrossBinSizes_nonCyanoDiaz.png"
ggsave(
  filename = fln,
  plot = multi_plot,
  height = 15,
  width = 25,
  dpi = 300
)
