# Header #############################################################
# 
# Author: Lisa Nicvert
# Email:  lisa.nicvert@univ-lyon1.fr
# 
# Date: 2023-03-22
#
# Script Description: this scripts installs all the necessary dependencies to
# run this project.



# Packages needed for the whole project -----------------------------------


## Install UnitEvents ------------------------------------------------------
install.packages("igraph")
# Warning, to install UnitEvents you need to install dependencies outside of R
# see instrictions at https://sourcesup.renater.fr/frs/?group_id=3267&release_id=3729
install.packages('UnitEvents_0.0.8.tar.gz', 
                 type = 'source', 
                 repos = NULL, INSTALL_opts = c('--no-lock'))

## Install Quarto ----------------------------------------------------------
# Needed to render Quarto documents
install.packages("quarto")


# Install devtools --------------------------------------------------------
# Needed for install_local
install.packages("devtools")

## Install camtrapHawkes ---------------------------------------------------
# # This is the custom package in the R/ subdirectory
devtools::install_local(upgrade = "never")

# Packages needed per analyses --------------------------------------------
# These packages will be installed only if they are not already installed

library(camtrapHawkes) # To load custom 'require' function

## 01_hawkes_process_simulation_example ------------------------------------
packages <- c("here", "RColorBrewer", "ggplot2")
base::lapply(packages, require)


## 02_simulation_interevent_times ------------------------------------------
packages <- c("here", "foreach", "doParallel",
              "magrittr", "dplyr", "tidyr")
base::lapply(packages, require)

packages <- c("here", "dplyr", "tidyr", "ggplot2")
base::lapply(packages, require)


## 03_simulation_hawkes ----------------------------------------------------
packages <- c("here", "foreach", "doParallel",
              "gridExtra", 
              "tidyr", "tibble", "dplyr", 
              "ggplot2")
base::lapply(packages, require)


## 04_plot_map -------------------------------------------------------------
packages <- c("here", "dplyr", 
              "sf", "sp", "rgdal", "ggsn", "osmdata")
base::lapply(packages, require)

packages <- c("here", "dplyr", "sf", "sp", "ggsn")
base::lapply(packages, require)


## 05_example_real_data ----------------------------------------------------
packages <- c("here", "dplyr", "lubridate", "ggplot2")
base::lapply(packages, require)

packages <- c("here", "lubridate",
              "ggplot2")
base::lapply(packages, require)


## 06_circadian_rhythm -----------------------------------------------------
packages <- c("here", "NHPoisson",
              "ggplot2", "tidyr", "dplyr")
base::lapply(packages, require)
