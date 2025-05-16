#### packages #### 
library(vegan)
library(ggplot2)


#### subsections functional types data prepping ####


func_type_subsection_summary <- readRDS("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/01 Vegetation changes Kobbefjord/data/nmds_nero/nmds_nero/data/func_type_subsection_summary.rds")

names(func_type_subsection_summary)


func_subsection_wide <- func_type_subsection_summary |> 
  pivot_wider(names_from = ecoveg_gfc, values_from = func_plot_fraction, values_fill = 0)

head(func_subsection_wide)

#### meta data for nmds ####

# Extract species data (gramnoid, herb, etc.)
species_data <- func_subsection_wide[, 6:11]

# Metadata for grouping (year and veg_type)
metadata <- func_subsection_wide[, c("year", "veg_type")]

#### nmds ana ####

set.seed(123)  # For reproducibility
nmds_result <- metaMDS(species_data, distance = "bray", trymax = 1000)
