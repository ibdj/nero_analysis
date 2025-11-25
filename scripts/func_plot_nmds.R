#### loading packages ####

library(tidyverse)
library(lme4)
library(lmerTest)
library(ggplot2)
library(ggeffects)
library(vegan)
library(codyn) # turnover calculations

#### loading data ################################

merged_data <- readRDS("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/01 Vegetation changes Kobbefjord/data/nero_analysis/data/merged_data.rds") |> 
  filter(veg_type != "saltmarsh", taxon_code != "rock") |> 
  select(year, plot_id, veg_type, taxon_code, presence, func_type, ecoveg_gfc, ecoveg_sgfc)

func_count_species <- merged_data |> 
  group_by(year, plot_id) |> 
  summarise(total_species = n_distinct(taxon_code))

func_count <- merged_data |> 
  group_by(year, plot_id, ecoveg_sgfc) |> 
  summarise(count_func = n_distinct(taxon_code))

func_plot <- func_count |> 
  left_join(func_count_species, by = c("year", "plot_id")) |> 
  mutate(frac = round(count_func/total_species,2))

unique(func_plot)

write_rds(func_plot, "data/func_plot.rds")


