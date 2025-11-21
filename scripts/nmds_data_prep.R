library(tidyverse)
library(vegan)

# Summarize duplicate rows before pivoting ###################################################################################
merged_data <- readRDS("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/01 Vegetation changes Kobbefjord/data/nero_analysis/data/merged_data.rds")

names(merged_data)

merged_data <- merged_data |>
  filter(presence == 1) |> 
  filter(!is.na(raunkiaer_value)) |>   # Remove rows with missing values
  filter(taxon_code != "rock", veg_type != "saltmarsh")

str(merged_data$presence)

merged_data_summary <- merged_data %>%
  group_by(plot_id, year, veg_type, taxon_code) %>%
  summarize(presence = max(presence), .groups = 'drop') %>%
  unite("plot_year_vt", plot_id, year,veg_type, sep = "_", remove = FALSE)

community_matrix <- merged_data_summary %>%
  select(plot_year_vt, taxon_code, presence) %>%
  pivot_wider(names_from = taxon_code, values_from = presence, values_fill = 0) %>%
  column_to_rownames(var = "plot_year_vt")

########## nmds ############ 

# Set seed for reproducibility
set.seed(42)

# Run NMDS on presence/absence community matrix
nmds_result <- metaMDS(community_matrix, k = 3, trymax = 100)

# Check stress value and results
print(nmds_result)

# Extract NMDS scores for plots
nmds_scores <- as.data.frame(scores(nmds_result))

# Add plot_year_vt as row names for easier merging with metadata if needed
nmds_scores$plot_year_vt <- rownames(nmds_scores)

# •	Stress < 0.05: Excellent representation, very reliable ordination
# •	Stress < 0.1: Good representation, generally acceptable
# •	Stress < 0.2: Fair or usable representation, but some distortion expected
# •	Stress > 0.2: Poor representation, ordination may not be reliable or interpretable
