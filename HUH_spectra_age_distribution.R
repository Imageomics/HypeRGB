################################################################################
# Packages
################################################################################
library("tidyverse")
library("ggplot2")
library("dplyr")
library("tidyr")
library("purrr")
library("stringr")
library("readr")
library("viridis")

################################################################################
# Import data
################################################################################
metadata     = read.csv("FloraPalooza2026/Data/metadata.csv")
huh_metadata = read.csv("FloraPalooza2026/Data/HUH_spectra_metadata/DMWhiteHUHspec1_sp25leaf560_metadata.csv")
huh_spec     = read.csv("FloraPalooza2026/Data/DMWhiteHUHspec1_sp25leaf560_resampled_5nm.csv")

# Rename Species column
huh_spec$Species = huh_spec$scientificName

###################
# Specimen age
###################
age_plot = ggplot(huh_spec, aes(x = Age)) +
  geom_histogram(fill = "lightblue", color = "black") +
  labs(title = "Absolute age", 
       x = "Age", y = "Count") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid = element_blank()
  )

ggsave("FloraPalooza2026/Figs/Absolute_age_of_huh_specimens.pdf", plot = age_plot, width = 8.5, height = 5, dpi = 600)
ggsave("FloraPalooza2026/Figs/Absolute_age_of_huh_specimens.png", plot = age_plot, width = 8.5, height = 5, dpi = 600)


################################################################################

################################################################################
# Distribution map of herbarium specimens
################################################################################
# Packages
library("sf")
library("dplyr")
library("ggplot2")
library("rnaturalearth")
library("rnaturalearthdata")
library("ggspatial")
library("viridis")

################################################################################
# Convert occurrence records to spatial data
################################################################################

sol_spat = metadata %>%
  filter(
    !is.na(decimalLatitude),
    !is.na(decimalLongitude)
  ) %>%
  st_as_sf(
    coords = c("decimalLongitude", "decimalLatitude"),
    crs = 4326
  )

################################################################################
# Get USA base map and remove Alaska and Hawaii
################################################################################

usa = ne_countries(country = "United States of America",
  scale = "medium",
  returnclass = "sf"
)

# Remove Alaska and Hawaii
# Contiguous USA bounding box
contiguous_usa = st_crop(usa,
  xmin = -125, xmax = -66,
  ymin = 24, ymax = 50)

################################################################################
# Keep only points within USA
################################################################################

sol_sf_usa = sol_spat[usa, , op = st_within]

################################################################################
# Zoom to specimen locations
################################################################################

# Get bounding box of specimen locations
points_bbox = st_bbox(sol_sf_usa)

# Add buffer around points
buffer = 2

xlim = c(points_bbox["xmin"] - buffer, points_bbox["xmax"] + buffer)

ylim = c(points_bbox["ymin"] - buffer, points_bbox["ymax"] + buffer)

################################################################################
# Plot distribution by species
################################################################################

usa_map = ggplot() +
  # USA
  geom_sf(data = contiguous_usa, fill = "white", color = "black", linewidth = 0.4) +
  # Specimen locations
  geom_sf(data = sol_sf_usa, aes(color = scientificName),
    alpha = 0.7, size = 0.8) +
  
  # Species colors
  scale_color_viridis_d(name = "Species", option = "D") +
  
  # North arrow
  # annotation_north_arrow(location = "tr", which_north = "true", 
  #   style = north_arrow_fancy_orienteering,
  #   pad_x = unit(0.2, "cm"), pad_y = unit(0.5, "cm")) +
  # 
  coord_sf(clip = "on") +
  labs(title = "Distribution of Herbarium Specimens",
       x = "Longitude", y = "Latitude") +
    
  theme_minimal() +
  
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
    #panel.grid = element_rect(),
    legend.position = "none",
    #legend.text = element_text(size = 9, face = "italic"),
    plot.margin = margin(10, 30, 10, 10)
  )

ggsave("FloraPalooza2026/Figs/Distribution of HUH Specimens.pdf", usa_map, height = 6, width = 8.5, dpi = 600)
ggsave("FloraPalooza2026/Figs/Distribution of HUH Specimens.png", usa_map, height = 6, width = 8.5, dpi = 600)


################################################################################
# Distribution map - New England
################################################################################

new_eng_map = ggplot() +
  # USA / New England base map
  geom_sf(data = usa, fill = "white", color = "black", linewidth = 0.4) +
  
  # Specimen locations
  geom_sf(data = sol_sf_usa, aes(color = scientificName),
          alpha = 0.7, size = 1) +
  
  # Species colors
  scale_color_viridis_d(option = "D") +
  
  # Zoom to New England
  coord_sf(xlim = c(-74, -66.5),
           ylim = c(40.5, 47.5),
           expand = FALSE) +
  
  theme_minimal() +
  
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
    panel.grid = element_blank(),
    legend.position = "none"
  ) +
  
  labs(title = "Distribution of Herbarium Specimens",
    x = "Latitude",
    y = "Longitude"
  )

ggsave("FloraPalooza2026/Figs/Distribution in New England.pdf", new_eng_map, height = 6, width = 6.5, dpi = 600)
ggsave("FloraPalooza2026/Figs/Distribution in New England.png", new_eng_map, height = 6, width = 6.5, dpi = 600)
