################################################################################
#                             Script README!
# Spectral plots for Pressed-leaf (black backgroung), Pseudo-herbarium (paper background), 
# and Herbarium (black and paper backgrounds) and their coefficient of variation (CV).
################################################################################
# Package(s)
################################################################################
library("spectrolab")
library("tidyverse")
library("ggplot2")
library("dplyr")
library("tidyr")
library("purrr")
library("stringr")
library("readr")
library("viridis")

################################################################################
# Import Herbarium Spectra
################################################################################
huh_spec       = read_csv("FloraPalooza2026/Data/DMWhiteHUHspec1_sp25leaf560_resampled_5nm_updated.csv")
huh_norm_spec  = read.csv("FloraPalooza2026/Data/DMWhiteHUHspec1_sp25leaf560_resampled_5nm_norm_updated.csv")

huh_spec$Species = huh_spec$scientificName
huh_norm_spec$Species = huh_norm_spec$scientificName

# Season: Early = March - May; Peak = June - September; Late = October - November
# To assign the startDayOfYear/dayOfCollection to months, I used https://www.scp.byu.edu/docs/doychart.html

# Define season using month column
huh_spec = huh_spec %>%
  mutate(
    Month = as.integer(month),
    Leaf_Dev = case_when(
      Month %in% 3:5   ~ "Young",
      Month %in% 6:9   ~ "Mature",
      Month %in% 10:11 ~ "Old_Senescent",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Leaf_Dev))

# Pivot - Build spectra long data
huh_long = huh_spec %>%
  pivot_longer(
    cols = where(is.numeric) & matches("^X\\d+$"),
    names_to = "Wavelength",
    values_to = "Reflectance"
  ) %>%
  mutate(Wavelength = as.numeric(sub("^X", "", Wavelength)))

# Mean spectra per Species × Leaf_Dev
huh_mean = huh_long %>%
  group_by(Species, Leaf_Dev, Wavelength) %>%
  summarise(
    Reflectance = mean(Reflectance, na.rm = TRUE),
    .groups = "drop"
  )

# Spectra Plot - Mean Spectra with 95% CI
##########################################################
huh_season_ci = huh_long %>%
  #filter(Species %in% spp_keep) %>%
  group_by(Species, Leaf_Dev, Wavelength) %>%
  summarise(
    mean_R = mean(Reflectance, na.rm = TRUE),
    sd_R   = sd(Reflectance, na.rm = TRUE),
    n      = sum(!is.na(Reflectance)),
    se_R   = sd_R / sqrt(n),
    ci_lo  = mean_R - 1.96 * se_R,
    ci_hi  = mean_R + 1.96 * se_R,
    .groups = "drop"
  )

# Plot
huh_plot = ggplot(huh_season_ci, aes(Wavelength, mean_R, color = Leaf_Dev, fill = Leaf_Dev)) +
  geom_line(linewidth = 0.8) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi),
              alpha = 0.2, color = NA) +
  facet_wrap(~ Species, scales = "fixed") +
  scale_color_manual(
    values = c(Young = "lightgreen", Mature  = "forestgreen", Old_Senescent  = "darkorange"),
    breaks = c("Young", "Mature", "Old_Senescent")) +
  scale_fill_manual(
    values = c(Young = "lightgreen", Mature  = "forestgreen", Old_Senescent  = "darkorange"),
    breaks = c("Young", "Mature", "Old_Senescent")) +
  theme_bw() +
  labs(title = "HU Herbarium - White et al. 2025",
       x = "Wavelength (nm)",
       y = "Reflectance (%)") + 
  theme_bw(base_size = 10) +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "italic", size = 8))


ggsave("FloraPalooza2026/Figs/HUH_spectra_leaf_dev.pdf", 
       plot = huh_plot, width = 14, height = 11.5, dpi = 600)
ggsave("FloraPalooza2026/Figs/HUH_spectra_leaf_dev.png", 
       plot = huh_plot, width = 14, height = 11.5, dpi = 600)


########################################################
# (2) CV spectra
########################################################
huh_cv = huh_long %>%
  group_by(Species, Wavelength) %>%
  summarise(
    mean_R = mean(Reflectance, na.rm = TRUE),
    sd_R   = sd(Reflectance, na.rm = TRUE),
    CV     = sd_R / mean_R,
    .groups = "drop"
  )

# Define species colors
colors = viridis(length(unique(huh_cv$Species)))

huh_cv_plot = ggplot(huh_cv, aes(Wavelength, CV, color = Species)) +
  geom_line(linewidth = 0.9) +
  theme_bw() +
  scale_color_manual(values = colors) + 
  labs(
    title = "HU Herbarium - White et al. 2025",
    x = "Wavelength (nm)",
    y = "Coefficient of Variation (CV)") + 
  theme(
    #axis.title.x = element_blank(),
    #axis.text.x  = element_blank(),
    #axis.ticks.x = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 8, face = "italic"),
    legend.title = element_text(size = 10, face = "bold")
  )

ggsave("FloraPalooza2026/Figs/Herbarium_CV_plots.pdf", plot = huh_cv_plot, width = 14, height = 8.5, dpi = 600)
ggsave("FloraPalooza2026/Figs/Herbarium_CV_plots.png", plot = huh_cv_plot, width = 14, height = 8.5, dpi = 600)


################################################################################
# Species mean spectral profiles
################################################################################
# Import data
################################################################################
metadata     = read.csv("FloraPalooza2026/Data/metadata.csv")
huh_metadata = read.csv("FloraPalooza2026/Data/HUH_spectra_metadata/DMWhiteHUHspec1_sp25leaf560_metadata.csv")
huh_spec     = read.csv("FloraPalooza2026/Data/DMWhiteHUHspec1_sp25leaf560_resampled_5nm.csv")

# Species name
huh_spec$Species = huh_spec$scientificName

# Convert spectral data from wide to long format
huh_long = huh_spec %>%
  pivot_longer(
    cols = matches("^X\\d+$"),
    names_to = "Wavelength",
    values_to = "Reflectance"
  ) %>%
  mutate(
    Wavelength = as.numeric(sub("^X", "", Wavelength))
  )

# Calculate mean spectrum for each species
huh_spec_mean = huh_long %>%
  group_by(Species, Wavelength) %>%
  summarise(
    mean_R = mean(Reflectance, na.rm = TRUE),
    .groups = "drop"
  )

# Colors
colors = viridis(length(unique(huh_spec_mean$Species)))

# Plot mean spectra
huh_spec_plot = ggplot(huh_spec_mean, aes(x = Wavelength, y = mean_R, color = Species)) +
  geom_line(linewidth = 0.9) +
  theme_bw() +
  scale_color_manual(values = colors) +
  labs(
    title = "HU Herbarium - White et al. 2025",
    x = "Wavelength (nm)",
    y = "Mean Reflectance",
    color = "Species") +
  theme(
    legend.position = "none",
    #legend.text = element_text(size = 9, face = "italic"),
    #legend.title = element_text(size = 10, face = "bold")
  )

ggsave("FloraPalooza2026/Figs/Herbarium_spec_mean.pdf", plot = huh_spec_plot, width = 9, height = 6.5, dpi = 600)
ggsave("FloraPalooza2026/Figs/Herbarium_spec_mean.png", plot = huh_spec_plot, width = 9, height = 6.5, dpi = 600)

