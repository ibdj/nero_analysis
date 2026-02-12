
#### loading packages ####

library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggplot2)
library(ggeffects)
library(vegan)
library(codyn) # turnover calculations
library(multcomp)
library(multcompView)
library(patchwork)

#### importing data ####
merged_data <- readRDS("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/01 Vegetation changes Kobbefjord/data/nero_analysis/data/merged_data.rds") |> 
  group_by(year, subsection) |> 
  mutate(no_plots = n_distinct(plot_id)) |> 
  ungroup()

pivot <- merged_data |> 
dplyr::select(year, vt_section,subsection, plot, plot_id, no_plots, veg_type, func_type, ecoveg_gfc, ecoveg_sgfc,presence, taxon_code) |> 
distinct(year, vt_section,subsection, plot, plot_id, no_plots, veg_type, presence, taxon_code)

meta <- merged_data |> 
  dplyr::select(taxon_code, func_type, ecoveg_gfc, ecoveg_sgfc) |> 
  distinct()

pivot_wide <- pivot |> 
  pivot_wider(names_from = taxon_code, values_from = presence, values_fill = 0 )

species_long <- pivot_wide |> 
 pivot_longer(cols = 8:ncol(pivot_wide), names_to = "taxon_code", values_to = "occurrence")

species_sub <- species_long |> 
  group_by(year,subsection,no_plots,veg_type, taxon_code) |> 
  summarise(count = sum(occurrence)) |> 
  mutate(fraction = count/no_plots)
  

species_sub_salgla <- species_sub |> 
  filter(taxon_code == "salgla")
