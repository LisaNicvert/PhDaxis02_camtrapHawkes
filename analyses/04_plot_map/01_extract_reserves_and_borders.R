# Header #############################################################
# 
# Author: Lisa Nicvert
# Email:  lisa.nicvert@univ-lyon1.fr
# 
# Date: 2022-10-14
#
# Script Description: Allows to import the shapefiles for each reserve


# Import libraries --------------------------------------------------------
# Load functions from main folder (for the require function)
library(camtrapHawkes)

packages <- c("here", "dplyr", 
              "sf", 
              "sp", "rgdal", 
              "ggsn", "osmdata")
base::lapply(packages, require)


# Define output folder ----------------------------------------------------
outputs_path <- here("outputs/04_plot_map")

# Get countries -----------------------------------------------------------

# Define bounding box
SA_dim <- sf::st_bbox(c(xmin = 14, xmax = 38.5, 
                        ymin = -35, ymax = -22),  
                      crs = CRS("+init=epsg:4326"))

# Get countries
q <- opq(bbox =  SA_dim, timeout = 60*3) %>%
  add_osm_feature(key = 'boundary', value = 'administrative') %>%
  add_osm_feature(key = 'admin_level', value = 2)
countries_res <- osmdata_sf(q) # Takes about 2 mins
countries <- countries_res$osm_multipolygons

countries <- countries %>% select(name)

# Get reserves -----------------------------------------------------------
bbox_dim <- sf::st_bbox(c(xmin = 26, xmax = 32, 
                          ymin = -28, ymax = -22),  
                        crs = CRS("+init=epsg:4326"))

q <- opq(bbox =  bbox_dim, timeout = 60*3) %>%
  add_osm_feature(key = 'boundary', value = 'protected_area')
parks_res <- osmdata_sf(q)

# Get relevant parks
parks <- parks_res$osm_multipolygons %>% 
  filter(name %in% c("Klaserie Private Nature Reserve",
                     "Kapama Game Reserve",
                     "Thornybush Game Reserve",
                     "Umbabat Nature Reserve",
                     "Kruger National Park",
                     "Madikwe Game Reserve",
                     "Pilanesberg National Park",
                     "Somkhanda Game Reserve",
                     "Balule Nature Reserve",
                     "Venetia Limpopo Nature Reserve")) %>% 
  select(name)


# Write to files ----------------------------------------------------------
st_write(parks, file.path(outputs_path, "reserves/reserves.shp"),
         append = FALSE,
         delete_layer = TRUE)
st_write(countries, file.path(outputs_path, "southern_africa/southern_africa.shp"),
         append = FALSE,
         delete_layer = TRUE)
