################################################################################
# Packages
################################################################################
library("spectrolab")
library("shiny")
library("stringr")
library("dplyr")
library("readr")
library("tools")

################################################################################
# Import data
################################################################################

huh_metadata = read.csv("FloraPalooza2026/Data/HUH_spectra_metadata/DMWhiteHUHspec1_sp25leaf560_metadata.csv")
huh_spec     = spectrolab::read_spectra(path = "FloraPalooza2026/Data/HUH_spectra_metadata/DMWhiteHUHspec1_raw_spectra_files/", extract_metadata = TRUE)

# Get filenames
# Specify the path to the directory containing spectra files
data_dir  = "FloraPalooza2026/Data/HUH_spectra_metadata/DMWhiteHUHspec1_raw_spectra_files/"

# List all the files in the directory
file_list = list.files(path = data_dir, full.names = TRUE)

# Extract base filenames
#names(huh_spec) = basename(file_list)

################################################################################
# Process spec
################################################################################

##########################
# Match_sensor_overlap
##########################
## Check the spectra and confirm the splicing bands
#spectrolab::plot_interactive(huh_spec)

# Match the same region with MAINE spectra
# However, Dawson et al. (2025) matched 991.3, 1902.5
spec_match = spectrolab::match_sensors(x = huh_spec, splice_at = c(990, 1900), interpolate = c(5,1))

####################
# Fix names
####################
#spec = lapply(spec_match, function(x){
#  y = x
#  names(y) = str_extract(names(y), "[^_]+")
#  y
#})

####################
# Trim
####################

# Define 5 nm wavelength grid from 450 to 2400 nm
# To correspond with Dawson et al. 2025
lambda <- seq(450, 2400, by = 5)

# Create FWHM values for the new wavelength grid
fwhm <- spectrolab::make_fwhm(spec_match, lambda)

# Resample spectra to the 5 nm grid
spec_match <- resample(
  spec      = spec_match,
  new_bands = lambda,
  fwhm      = fwhm
)

###############################
# Normalize spectra
###############################
spec_match_norm = spectrolab::normalize(spec_match)

###############################
# Match metadata and spectra

metadata = huh_metadata %>%
  mutate(file_base = tools::file_path_sans_ext(filename))

# Assign base filenames as names
names(spec_match) = tools::file_path_sans_ext(basename(file_list))

names(spec_match_norm) = tools::file_path_sans_ext(basename(file_list))

####################################
# Convert spectra to data.frame
####################################
spec_df = as.data.frame(spec_match)
spec_df$file_base = tools::file_path_sans_ext(names(spec_match))

spec_norm_df = as.data.frame(spec_match_norm)
spec_norm_df$file_base = tools::file_path_sans_ext(names(spec_match_norm))

##########################################
# Aggregate spectra by filename (mean)
##########################################
spec_mean_df = spec_df %>%
  group_by(file_base) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop")

spec_mean_norm_df = spec_norm_df %>%
  group_by(file_base) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop")

###############################
# Match metadata
###############################
full_spec_df = metadata %>%
  inner_join(spec_mean_df, by = "file_base")

full_spec_norm_df = metadata %>%
  inner_join(spec_mean_norm_df, by = "file_base")


###############################
# Export to CSV
###############################

write.csv(full_spec_df, "FloraPalooza2026/Data/DMWhiteHUHspec1_sp25leaf560_resampled_5nm.csv", row.names = FALSE)

write.csv(full_spec_norm_df, "FloraPalooza2026/Data/DMWhiteHUHspec1_sp25leaf560_resampled_norm_5nm.csv", row.names = FALSE)

