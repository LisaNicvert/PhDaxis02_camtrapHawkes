# Header #############################################################
# 
# Author: Lisa Nicvert
# Email:  lisa.nicvert@univ-lyon1.fr
# 
# Date: 2024-05-19
#
# Script Description: plot species abundances

library(here)

library(ggplot2)
library(dplyr)
library(ggtext)

# Parameters --------------------------------------------------------------
input_path <- here("outputs/05_example_real_data")
figures_path <- here("figures/05_example_real_data")

silhouettes_path <- here("data/species_silhouettes/")


# Read data ---------------------------------------------------------------
# --- Silhouettes path
silhouettes_labs <- c(impala = paste0("<img src='", silhouettes_path, 
                                      "impala.png'", " height='20'/>", 
                                      "<br>impala"),
                      kudu = paste0("<img src='", silhouettes_path, 
                                    "kudu.png'", " height='20'/>", 
                                    "<br>kudu"),
                      lion = paste0("<img src='", silhouettes_path, 
                                    "lion.png'", " height='20'/>", 
                                    "<br>lion"),
                      wildebeestblue = paste0("<img src='", silhouettes_path, 
                                              "wildebeest.png'", " height='20'/>",
                                              "<br>wildebeest"),
                      zebraburchells = paste0("<img src='", silhouettes_path, 
                                              "zebra.png'", " height='20'/>",
                                              "<br>zebra"))

dat <- read.csv(file.path(input_path, "cleaned_data.csv"))

# Plot --------------------------------------------------------------------
dat <- dat |> 
  group_by(snapshotName) |> 
  summarize(n = n())

ggplot(dat, aes(y = reorder(snapshotName, n),
                x = n)) +
  geom_col() +
  geom_text(aes(label = n), 
            hjust = 1,
            nudge_x = -0.05,
            col = "white") +
  theme_linedraw(base_size = 13) +
  scale_x_log10() +
  scale_y_discrete(labels = silhouettes_labs) +
  xlab("Capture count") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_markdown(color = "black",
                                       size = 10,
                                       halign = 0.5,
                                       hjust = 0.5))
# Save figure
ggsave(file.path(figures_path, "spp_count.jpeg"),
       bg = "white",
       width = 8, height = 4,
       dpi = 300)
           