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



### ==============================================================
### BETADIVERSITY - One location across 12 months averaged
### ==============================================================
# 
library(raster)





## Required data & data formatting
# Read shapefile with geometries of ocean basins
sf.ocean <- st_read("/net/kryo/work/deriksson/Projects/Notion_DE/Data/Ocean_regions_shapefile/GOaS_v1_20211214/goas_v01.shp")
# Pacific centered world map: Source: https://stackoverflow.com/questions/56146735/visual-bug-when-changing-robinson-projections-central-meridian-with-ggplot2
worldMap <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_make_valid()
# Set projection pacific centered
lon_pacific <- '200' 
target_crs <- st_crs( paste0("+proj=eqc +x_0=0 +y_0=0 +lat_0=0 +lon_0=", lon_pacific) )

# define a long & slim polygon that overlaps the meridian line & set its CRS to match
# that of world
# Centered in lon 200
offset <- 180 - as.numeric(lon_pacific)


polygon <- st_polygon(x = list(rbind(
  c(-0.0001 - offset, 90),
  c(0 - offset, 90),
  c(0 - offset, -90),
  c(-0.0001 - offset, -90),
  c(-0.0001 - offset, 90)
))) %>%
  st_sfc() %>%
  st_set_crs(4326)


# modify world dataset to remove overlapping portions with world's polygons
world2 <- worldMap %>% st_difference(polygon)
#> Warning: attribute variables are assumed to be spatially constant throughout all
#> geometries
# Transform
sf.ocean_pacific <- world2 %>% st_transform(crs = target_crs)



filenames <- list.files("/home/deriksson/Projects/Notion_DE/Code/Diversity/Diversities/pa/BetaDiversity_acrossOneGridCell_andAllTwelveMonth/")
filenames <- grep(".grd", filenames, value = TRUE)
r <- brick(paste0("/home/deriksson/Projects/Notion_DE/Code/Diversity/Diversities/pa/BetaDiversity_acrossOneGridCell_andAllTwelveMonth/", filenames[1]))

for(i in seq_along(filenames)){

  # Print progress
  print(i)

  if(i == 1){
    r <- brick(paste0("/home/deriksson/Projects/Notion_DE/Code/Diversity/Diversities/pa/BetaDiversity_acrossOneGridCell_andAllTwelveMonth/", filenames[i]))
    df <- as.data.frame(r, xy = TRUE)
  }else{
    r <- brick(paste0("/home/deriksson/Projects/Notion_DE/Code/Diversity/Diversities/pa/BetaDiversity_acrossOneGridCell_andAllTwelveMonth/", filenames[i]))
    df <- rbind(df, as.data.frame(r, xy = TRUE) )
  }
}

# Merge
test <- df
test$id <- paste0( test$x, "_", test$y )

# Aggregate
test_mean <- aggregate( .~id, data = test, FUN = mean )
test_mean <- test_mean[, -1]
#
test_sd <- aggregate( .~id, data = test, FUN = sd )
test_sd <- test_sd[, -1]
test_sd$x <- test_mean$x
test_sd$y <- test_mean$y


# Convert to raster
dat <- test_mean
dat <- dat[, c("x", "y", "jac")]
library(terra)
dat <- terra::rast(dat)
# Set reference system
crs(dat) <- 'epsg:4326'

# Calculate upper bound 70th percentile of coefficient of variation

# df <- cbind(
#   as.data.frame(r$Diversity_Richness_ensembleMean_paBased, xy = TRUE),
#   as.data.frame(r$Diversity_Richness_ensembleSD_paBased)
# )
# df$coefficient_of_variation <- df$Diversity_Richness_ensembleSD_paBased/df$Diversity_Richness_ensembleMean_paBased
# percentile <- quantile(df$coefficient_of_variation, probs = 0.7, na.rm = TRUE) # 0.33
# df1 <- df[which(!is.na(df$coefficient_of_variation)), ]
# df1 <- df1[which(df1$coefficient_of_variation > percentile  ), ]
# hb  <- erode(hexbin(df1$x, df1$y, xbins = 50))
# df70 <- as.data.frame(hcell2xy(hb))
# # Convert stippling points to the Pacific-centered CRS
# stipple_points_sf <- st_as_sf(df70, coords = c("x", "y"), crs = 4326)
# stipple_points_sf <- st_transform(stipple_points_sf, crs = target_crs)

# Plot using ggplot2
library(ggplot2)
gg_diversity <- 
  
  # Plot world map
  ggplot(sf.ocean_pacific) + 
  geom_sf(fill = "grey") +
  
  # Plot countour plot & and contour line settings
  geom_spatraster_contour_filled(data = dat, breaks = seq(0, 1, 0.05), show.legend = T) + 
  # geom_spatraster_contour(
  #   data = dat,
  #   color = "grey",
  #   linewidth = 0.1) +

  # Extend legend values to the break values above, froom 0 to 1
  scale_fill_viridis_d(drop = FALSE) +


  # Add stippling upper 70th percentile bound
  # geom_sf(data = stipple_points_sf, color = "white", alpha = 0.5, size = 0.5) +1  # Legend settings
  # guides(fill = guide_legend(nrow = 1, title = "Ensemble species turnover", title.position = "top")) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.title.align = 0.5,
    legend.key.size = unit(10, 'cm'), # change legend key size
    legend.key.height = unit(0.3, 'cm'), # change legend key height
    legend.key.width = unit(0.3, 'cm'),
    legend.title = element_text(size = 10), # change legend title font size
    legend.text = element_text(size = 8),

    # Plot background color
    panel.background = element_rect(fill = "white")) +

    # Axis labels
    xlab("") +
    ylab("") 


### =======================================
### Visualize annual biogeography HSI/PA
### ======================================
## HSI
# Directories
wd_in <- "/home/deriksson/Projects/Notion_DE/Code/Ensembles/Annual_individual_ensemble_biogeographies/hsi/"
# wd_in <- "/home/deriksson/Projects/Notion_DE/Code/Ensembles/Annual_individual_ensemble_biogeographies/pa/"
wd_out <- "/home/deriksson/Projects/Notion_DE/Code/Ensembles/Visualizations/hsi/"
# wd_out <- "/home/deriksson/Projects/Notion_DE/Code/Ensembles/Visualizations/pa/"

# # Load table for boxplot of TSS scores
# # Vector of filenames
# filenames <- list.files( paste0("/home/deriksson/Projects/Notion_DE/Code/6_SDM_fit/3_Output/Total_dataset/") )

# # Open file
# l <- list()
# for(i in seq_along(filenames)){

#     # Print progress
#     print( i )

#     # Open file
#     d <- read.csv( paste0("/home/deriksson/Projects/Notion_DE/Code/6_SDM_fit/3_Output/Total_dataset/", filenames[i]) )
#     d$id <- as.factor(rep(1:5, (nrow(d)/5) ))
#     d$id <- paste0( d$id, filenames[i] )

#     # Add to list
#     l[[i]] <- d

# } # end of loop across filenames

# Merge list into dataframe
# df_eval <- do.call( "rbind.fill", l )


# Load mean and sd
filenames <- list.files(wd_in)
filenames <- grep(".grd", filenames, value = TRUE)

l <- list()
l_linePlots <- list()
for(f in seq_along(filenames)){
    d <- brick( paste0(wd_in, filenames[f]) )

    # Load number of observations and observations
    # Directories
    obs <- get(load("/home/deriksson/Projects/Notion_DE/Code/4_Generate_absences/1_Output/Total_dataset/pres_abs,cr_bg_nonov(T_MLD1).RData"))
    obs <- obs[[ gsub("_annual_ensemble_mean_sd_individualEnsembleMembers.grd", "", filenames[f]) ]]
    obs <- rasterFromXYZ( obs[, c("x", "y", "obs")] )
    names(obs) <- "obs"

    # Merge with existing 
    obs <- resample(obs, d, method = "bilinear")
    d <- stack(d, obs)

    # Mean annual hsi
    dat1 <- as.data.frame(d[["Annual.ensemble.mean"]], xy = TRUE)
    dat1 <- rast(dat1)
    # Set reference system
    crs(dat1) <- 'epsg:4326'
    #
    df_hsi <- as.data.frame( dat1, xy = TRUE )
    sf_points_hsi <- st_as_sf(df_hsi, coords = c("x", "y"))
    # Define the CRS
    crs <- st_crs("+proj=longlat +datum=WGS84")
    # Set CRS for sf_points
    st_crs(sf_points_hsi) <- crs

    # Observations
    df_obs <- as.data.frame( d[["obs"]], xy = TRUE )
    df_obs <- df_obs[which(df_obs$obs == 1), ]
    sf_points_obs <- st_as_sf(df_obs, coords = c("x", "y"))
    # Define the CRS
    crs <- st_crs("+proj=longlat +datum=WGS84")
    # Set CRS for sf_points
    st_crs(sf_points_obs) <- crs

    # Plot world map mean hsi
    global_map <- ggplot(sf.ocean_pacific) + 
            geom_sf(fill = "grey") +
                    
            # Plot countour plot & and contour line settings
            geom_spatraster_contour_filled(data = dat1,breaks = seq(0, 1, 0.05), show.legend = FALSE) +
            ggtitle(paste0(gsub("_annual_ensemble_mean_sd_individualEnsembleMembers.grd", "", filenames[f]))) +  
            
            # Extend legend values to the break values above, froom 0 to 1
            scale_fill_viridis_d(drop = FALSE) +
            
            # Add observations
            geom_sf(data = sf_points_obs, aes(size = obs), color = "red", size = 1, alpha = 0.4) +
            coord_sf() +
            # Add annotations
            # geom_text(x = 0, y = 0, label = paste0("Number of observations: ", sum(sf_points_obs$obs)), hjust = 2.2, vjust = 18, color = "red", alpha = 0.4) +
            # geom_text(x = 0, y = 0, label = "Number of ensemble members: 90", hjust = 1.8, vjust = 21, color = "red", alpha = 0.4) +
            theme(
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              legend.position = "bottom",
              legend.title.align = 0.5,
              legend.key.size = unit(10, 'cm'), # change legend key size
              legend.key.height = unit(0.3, 'cm'), # change legend key height
              legend.key.width = unit(0.3, 'cm'),
              legend.title = element_text(size = 10), # change legend title font size
              legend.text = element_text(size = 8),
              plot.title = element_text(size = 40),

              # Plot background color
              panel.background = element_rect(fill = "white")) 

# Save maps in list
l[[f]] <- global_map


# Save lineplots -----------------------------
# Line plot ensemble members
test <- d[[ grep( "SDM", names(d) ) ]]
test <- as.data.frame( test, xy = TRUE )
# Aggregate
test <- test[, -1  ]
test <- aggregate( .~ y, data = test, mean )

# From wide to longformat
test <- gather(test, sdm, hsi, SDM_1:names(test)[length(names(test))], factor_key = TRUE)

# Aggregate to get mean line
mean_line <- aggregate( hsi ~ y, data = test, FUN = mean )

# Plot lineplot
lineplot <- ggplot() +
    geom_path(data = test, aes(x = hsi, y = y, color = sdm), alpha = 0.3, show.legend = F) +
    geom_path(data = mean_line, aes(x = hsi, y = y), color = "black", show.legend = F, size = 1.5) +
    guides(color = FALSE) +
    ggtitle( gsub("_annual_ensemble_mean_sd_individualEnsembleMembers.grd", "", filenames[f]) ) +
    xlab( "Mean habitat suitability" ) +
    ylab("Latitude [°]") +
    theme(
      
      axis.text.x = element_text(size = 32, angle = 90, vjust = 0.5),
      axis.text.y = element_text(size = 32),
      axis.title.x = element_text(size = 32),
      axis.title.y = element_text(size = 32),
      plot.title = element_text(size = 40),
      panel.background = element_rect(fill = "white")
    ) +
    scale_y_continuous(breaks = seq(min(test$y), max(test$y), by = 5)) 


l_linePlots[[f]] <- lineplot

} # end loop across filenames


## Multiple plot
# Arrange
g2 <- ggarrange(
    l[[1]], l[[2]], l[[3]],
    l[[4]], l[[5]], l[[6]],
    l[[7]], l[[8]], l[[9]],
    l[[10]], l[[11]], l[[12]],
    l[[13]], l[[14]], l[[15]], 
    ncol = 3, nrow = 5, common.legend = TRUE
    )

g3 <- ggarrange(
    l_linePlots[[1]], l_linePlots[[2]], l_linePlots[[3]],
    l_linePlots[[4]], l_linePlots[[5]], l_linePlots[[6]],
    l_linePlots[[7]], l_linePlots[[8]], l_linePlots[[9]],
    l_linePlots[[10]], l_linePlots[[11]], l_linePlots[[12]],
    l_linePlots[[13]], l_linePlots[[14]], l_linePlots[[15]],
    ncol = 5, nrow = 3, common.legend = TRUE, common.axis = TRUE
    )



# ggsave
ggsave(
    g2,
                filename = paste0(wd_out, "Annual_mean_biogeographies.png"),
                height = 20,
                width = 20,
                bg = "white"
                )

# ggsave
ggsave(
    g3,
                filename = paste0(wd_out, "Annual_mean_latitudinalGradients.png"),
                height = 20,
                width = 15,
                bg = "white"
                )

