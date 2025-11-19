#### packages ####

library(tidyverse)

#### nmds with subsections (5 plots) ####

nero_subsection_summary <- readRDS("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/01 Vegetation changes Kobbefjord/data/nmds_nero/nmds_nero/data/nero_subsection_summary.rds")
nero_subsection_summary <- as.data.frame(readRDS("~/nmds_nero/data/nero_subsection_summary.rds"))
  

subsection_wide <- nero_subsection_summary |> 
  pivot_wider(names_from = species, values_from = species_plot_fraction, values_fill = 0 )

subsection_wide_2007 <- subsection_wide |> 
  filter(year == 2007) |> 
  select(-c(vt_section_id,year,total_occurrences,plots_containing_species,n_unique_plot_id))

# Metadata (includes veg_type)
meta_data <- subsection_wide_2007 |> select(veg_type)

# Community matrix (exclude veg_type)
community_matrix <- subsection_wide_2007 |> select(-veg_type)




