#### packages ####

library(tidyverse)
library(vegan)
library(dplyr)
library(ggplot2)
library(purrr)
library(lme4)
library(lmerTest)
library(ggeffects)

#### importing data ####

merged_data <- readRDS("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/01 Vegetation changes Kobbefjord/data/nero_analysis/data/merged_data.rds") |> 
  group_by(year, vt_section) |> 
  rename(section = vt_section) |> 
  mutate(no_plots = n_distinct(plot_id)) |> 
  ungroup() |> 
  filter(taxon_code != "rock", no_plots > 1)

pft_sec_long <- merged_data |>
  group_by(year, section, veg_type, taxon_code, no_plots, ecoveg_sgfc) |>
  reframe(
    count = n(),
    fraction_sec = as.double(count) / as.double(no_plots),
    .groups = "drop"
  ) |>
  filter(veg_type != "saltmarsh", no_plots != 1) |> 
  distinct(year, section, veg_type, ecoveg_sgfc, no_plots, fraction_sec, taxon_code)

func_sec_frac_sum <- pft_sec_long |> 
  group_by(year, section, veg_type, ecoveg_sgfc) |> 
  reframe(frac_sum = sum(fraction_sec))

func_sec_wide <- func_sec_frac_sum |> 
  pivot_wider(names_from = ecoveg_sgfc, values_from = frac_sum, values_fill = 0) |> 
  rename(graminoid = herb_graminoid, bryophyte = non_vascular_bryophyte, decidous = shrub_decidous, evergreen = shrub_evergreen, lichen = non_vascular_lichen)

community_matrix <- func_sec_wide |> 
  unite("sec_year_vt", section, year,veg_type, sep = "_", remove = FALSE) |> 
  column_to_rownames(var = "sec_year_vt") |> 
  select(-year, -section, -veg_type)

#### NMDS ######################################################################

# Set seed for reproducibility
set.seed(42)

# Run NMDS on presence/absence community matrix
nmds_result <- metaMDS(community_matrix, k = 3, trymax = 100)

print(nmds_result)

nmds_scores <- as.data.frame(scores(nmds_result, display = "sites"))

nrow(community_matrix)
nrow(nmds_scores) 

# Extract NMDS scores for plots
nmds_scores <- as.data.frame(scores(nmds_result, display = "sites"))
nmds_scores$sub_year_vt <- rownames(nmds_scores)
