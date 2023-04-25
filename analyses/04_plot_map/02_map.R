# Header #############################################################
# 
# Author: Lisa Nicvert
# Email:  lisa.nicvert@univ-lyon1.fr
# 
# Date: 2022-10-14
#
# Script Description: Plots a map with the cameras that were used 
# to analyze the data.


# Import libraries -----------------------------------------------------
# Load functions from main folder (for the require function)
library(camtrapHawkes)

packages <- c("here", "dplyr", "sf", "sp", "ggsn")
base::lapply(packages, require)

# Read data ------------------------------------------------------------

# --- Countries borders
southern_africa <- st_read(here("outputs/04_plot_map/southern_africa/southern_africa.shp"))

# --- Parks
parks <- read_sf(here("outputs/04_plot_map/reserves/reserves.shp"))

# --- Figures
figures_path <- here("figures/04_plot_map")

# Select useful parks -------------------------------------------
# Remove Kapama
parks <- parks %>% filter(name != "Kapama Game Reserve")

# Select relevant columns
parks <- parks %>%
  select(name)

# Add reserve name
pnames <- unique(parks$name)
reserve <- c("KRU", "APN", "APN", "MAD",
             "PLN", "APN", "VEN", "SOM",
             "APN")
names(reserve) <- pnames

parks <- parks %>% 
  mutate(code = reserve[name])

# Prepare data --------------------------------------------------------

## Merge APN reserves ------
parks_APN <- parks %>% 
  filter(code == "APN")
parks_APN <- st_union(parks_APN)

# Convert to dataframe and add relevant columns
parks_APN <- st_sf(data.frame(name = "APN",
                              code = "APN",
                              geom = parks_APN))

## Replace multiple polygons with APN ------
parks <- parks %>%
  filter(code != "APN") %>%
  bind_rows(parks_APN)

# Arrange data
parks <- parks %>% arrange(name)

## Create map labels ------
sites_label <- st_centroid(parks)

# Create labels for sites
# Get first word
label <- regmatches(x =  sites_label$name, 
                    regexpr(pattern = "^[a-zA-Z]+", sites_label$name))
sites_label$label <- label

## Projection ------
CRS("+init=epsg:2046") # Inspiration of the string used for CRS
projstring <- "+proj=tmerc +axis=enu +lat_0=0 +lon_0=15 +k=1 +x_0=0 +y_0=0
+ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

init_crs <- CRS("+init=epsg:4326")

southern_africa <- southern_africa %>% 
  st_set_crs(init_crs) %>% 
  st_transform(projstring)

sites_label <- sites_label %>%
  st_set_crs(init_crs) %>% 
  st_transform(projstring)

parks <- parks %>%
  st_set_crs(init_crs) %>%
  st_transform(projstring)

# Plot -----------------------------------------------------------------

## Colors -----
blue <- "#365f95"
fill <- "antiquewhite2"
border <- "grey42"
park_fill <- "#bdde86"
park_border <- "darkolivegreen"

## Southern Africa inset -----
# --- Limits for South Africa inset
SA_dim <- sf::st_bbox(c(xmin = 14, xmax = 38.5, 
                        ymin = -35, ymax = -22),  
                      crs = init_crs)

SA_dim = SA_dim %>%
  st_as_sfc() %>%
  st_transform(crs = projstring) %>%
  st_bbox()

# --- Limits for map (rectangle on inset)
map_dim <- st_bbox(parks)
names(map_dim) <- c("xmin", "ymin", "xmax", "ymax")

map_dim["xmax"] <- map_dim["xmax"] + 300000
map_dim["xmin"] <- map_dim["xmin"] - 20000

SA_inset <- ggplotGrob(
  ggplot(southern_africa) +
  geom_sf(data = southern_africa, 
          lwd = .2, 
          fill = fill,
          col = border) +
  blank() +
  theme(axis.title = element_blank(),
        axis.text = element_text(size=15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black", 
                                    linewidth = 1.1, linetype = "solid"),
        panel.background = element_rect(fill = blue)) +
    annotate(geom = "rect",
             ymin = map_dim["ymin"], ymax = map_dim["ymax"],
             xmin = map_dim["xmin"], xmax = map_dim["xmax"],
             colour = "black", fill = NA,
             linewidth = 1) +
    coord_sf(xlim = SA_dim[c("xmin", "xmax")], 
             ylim = SA_dim[c("ymin", "ymax")],
             crs = projstring)
  )


## Main map -----

# --- Nudge formatting
nudge_x <- c(-1, 0.5, 0.3,
             0.3, -2, 0.2)*10^5
# "APN" "KRU" "MAD" "PLN" "SOM" "VEN" 


(p <- ggplot(parks) +
    # Plot countries
    geom_sf(data = southern_africa, lwd = .4, 
            fill = fill,
            col = border) +
    # Plot parks
    geom_sf(data = parks, 
            fill = park_fill,
            col = park_border,
            show.legend = FALSE,
            linewidth = 0.4) +
    # Plot labels
    geom_sf_text(data = sites_label,
               aes(label = label),
               show.legend = FALSE,
               size = 18/.pt,
               hjust = 0,
               nudge_x = nudge_x) +
    # Map limits
    coord_sf(xlim = map_dim[c("xmin", "xmax")],
             ylim = map_dim[c("ymin", "ymax")],
             crs = projstring) +
    # Scale
    ggsn::scalebar(location = "topright",
                   transform = FALSE,
                   dist_unit = "km",
                   st.size = 6, 
                   st.dist = .04,
                   border.size = 0.5,
                   x.min = map_dim["xmin"],
                   x.max = map_dim["xmax"] - 30000,
                   y.min = map_dim["ymin"],
                   y.max = map_dim["ymax"],
                   dist = 100) +
    # North arrow
    north(location = "topleft",
          x.min = map_dim["xmin"], x.max = map_dim["xmax"],
          y.min = map_dim["ymin"], y.max = map_dim["ymax"],
          symbol = 12) +
    # Cosmetics
    xlab("Longitude") +
    ylab("Latitude") +
    theme(axis.title = element_blank(),
          axis.text = element_text(size = 18, color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = blue),
          panel.border = element_rect(colour = "black", fill=NA,
                                      linewidth = 0.8))
    )


## Final plot ----------------------------
(p + annotation_custom(grob = SA_inset, 
                       xmin = map_dim["xmin"] - 100000, 
                       xmax =  map_dim["xmin"] + 400000,
                       ymin = map_dim["ymin"], 
                       ymax = map_dim["ymin"]+200000))

ggsave(file.path(figures_path, "map.jpeg"), 
       width =  (map_dim["xmax"] - map_dim["xmin"])/100000,
       height = (map_dim["ymax"] - map_dim["ymin"])/100000,
       dpi = 600)
