# Header #############################################################
# 
# Author: Lisa Nicvert
# Email:  lisa.nicvert@univ-lyon1.fr
# 
# Date: 2023-04-17
#
# Script Description: run all the analyses (except the analyses from
# 02_simulation_interevent_times) to check that it works.

library(here)
library(quarto)

# 01_hawkes_process_simulation_example ------------------------------------
quarto::quarto_render(here("analyses/01_hawkes_process_simulation_example/hawkes_process_simulation_example.qmd"))

# 02_simulation_interevent_times ------------------------------------------
quarto::quarto_render(here("analyses/02_simulation_interevent_times/02_analyze_interevent_times.qmd"))

# 03_simulation_hawkes ----------------------------------------------------
quarto::quarto_render(here("analyses/03_simulation_hawkes/simulation_inference_hawkes.qmd"))

# 04_plot_map -------------------------------------------------------------
source(here("analyses/04_plot_map/01_extract_reserves_and_borders.R"))
source(here("analyses/04_plot_map/02_map.R"))

# 05_example_real_data ----------------------------------------------------
source(here("analyses/05_example_real_data/01_prepare_data.R"))
quarto::quarto_render(here("analyses/05_example_real_data/02_infer.qmd"))

# 06_circadian_rhythm -----------------------------------------------------
quarto::quarto_render(here("analyses/06_circadian_rhythm/circadian_rhythm.qmd"))
