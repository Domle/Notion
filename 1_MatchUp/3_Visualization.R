## Visualize outputs from this section. This is also a sanity to spot files that might have 
## been compiled wrong in script 1. For script 2 we will plot the presences in temperature-nitrate space.

## Dominic Eriksson
## Environmental Physics Group
## ETH Zurich
## Switzerland

## dominic.eriksson@usys.ethz.ch, 19th of June 2024


## Load libraries
lib_vec <- c("raster", "ggplot2", "dplyr", "ggpubr", "gridExtra",
  "betapart", "rnaturalearth", "terra", "tidyterra", "sf")

# Install and load packages
lapply(lib_vec, function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
})


### ===============================================
### Pacific Centered Global Map - Prep
### ===============================================
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


### ======================================================================
### 1_Output
### ======================================================================
# Directory
wd_in <- "/home/deriksson/Projects/Notion_DE/Code/GitHub_PublicRepo/1_MatchUp/1_Output/"

# Filenamnes
filenames <- list.files(wd_in)

# Load data
l <- lapply( paste0(wd_in, filenames), read.csv )

# Loop through and save plots in list
list_plots <- list()
for(f in seq_along(filenames)){

  # Print progress
  print( paste0("Plotting file ", f, ", out of ", length(filenames), ".") )

  # Subset list
  df <- l[[f]]

  #  Extract columns of interest
  df <- df[, c("x", "y", "occurrenceStatus")]

  # Convert to spatial object
  df_sf <- sf::st_as_sf( df, coords = c("x", "y"), crs = 4326 )

  ## Plot  
  # Plot pacific centered worldmap with presence and absences (across all taxa)
  gg_plot <-  ggplot(sf.ocean_pacific) + 
    geom_sf(fill = "grey") +
    # Add presences and absences
    geom_sf(data = df_sf, aes(color = occurrenceStatus), alpha = 0.5, size = 0.5) + # Legend settings
    ggtitle(filenames[f])

  # Save in list
  list_plots[[f]] <- gg_plot
} # close loop through filenames

# Arrange the plots in a 2-column, 3-row grid
grid.arrange(grobs = list_plots, ncol = 2, nrow = 3)

### --------------------------------------------------------------------------
### 2_Output
### --------------------------------------------------------------------------
# Directory
wd_in <- "/home/deriksson/Projects/Notion_DE/Code/GitHub_PublicRepo/1_MatchUp/2_Output/"

# Get filenames
filenames <- list.files(wd_in)

# Load data
l <- lapply( paste0(wd_in, filenames), read.csv )

# Loop through and save plots in list
list_plots <- list()
for(f in seq_along(filenames)){

  # Print progress
  print( paste0("Plotting file ", f, ", out of ", length(filenames), ".") )

  # Subset list
  df <- l[[f]]

  # Extract columns of interest
  df <- df[, c("T", "N", "occurrenceStatus")]
  # Remove absences
  df <- df[which(df$occurrenceStatus == "PRESENT"), ]

  # Plot
  gg_plot <- ggplot(data = df) +
    geom_point(aes( x = T, y = N )) + 
    ggtitle(filenames[f])

  # Save in list
  list_plots[[f]] <- gg_plot
} # close loop across filenames

# Multiplot
grid.arrange(grobs = list_plots, ncol = 2, nrow = 3)


### --------------------------------------------------------------------------
### END
### --------------------------------------------------------------------------