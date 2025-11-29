#### packages ####

library(tidyverse)
library(vegan)
library(dplyr)
library(ggplot2)
library(purrr)

#### importing data ####

merged_data <- readRDS("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/01 Vegetation changes Kobbefjord/data/nero_analysis/data/merged_data.rds") |> 
  group_by(year, subsection) |> 
  mutate(no_plots = n_distinct(plot_id)) |> 
  ungroup() |> 
  filter(taxon_code != "rock")

species_sub_long <- merged_data |>
  group_by(year, subsection, veg_type, taxon_code, no_plots, ecoveg_sgfc) |>
  reframe(
    count = n(),
    fraction_sub = as.double(count) / as.double(no_plots),
    .groups = "drop"
  ) |>
  filter(veg_type != "saltmarsh", no_plots != 1) |> 
  distinct(year, subsection, veg_type, ecoveg_sgfc, taxon_code, no_plots, fraction_sub)

func_sub_frac_sum <- species_sub_long |> 
  group_by(year, subsection, veg_type, ecoveg_sgfc) |> 
  reframe(frac_sum = sum(fraction_sub))

func_sub_wide <- func_sub_frac_sum |> 
  pivot_wider(names_from = ecoveg_sgfc, values_from = frac_sum, values_fill = 0)


  
