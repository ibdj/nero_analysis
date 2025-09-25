#### nmds with 10 agreegate ####
library(tidyverse)

merged_data <- readRDS("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/01 Vegetation changes Kobbefjord/data/nero_analysis/data/merged_data.rds") |> 
  select(year,vt_section,plot_id,veg_type,taxon_code,presence,species,ecoveg_gfc)

### calculating percent pr plot ####
plots_per_section <- merged_data %>%
  distinct(year, vt_section, plot_id) %>%
  group_by(year, vt_section) %>%
  summarise(n_plots_total = n(), .groups = "drop")

presence_per_plot <- merged_data %>%
  group_by(year, vt_section, plot_id, species,veg_type) %>%
  summarise(presence = as.integer(any(presence == 1)), .groups = "drop")

species_percent <- presence_per_plot %>%
  group_by(year, vt_section, species,veg_type) %>%
  summarise(n_plots_present = sum(presence), .groups = "drop") %>%
  left_join(plots_per_section, by = c("year", "vt_section")) %>%
  mutate(
    perc_plots = 100 * n_plots_present / n_plots_total
  ) %>%
  arrange(year, vt_section, desc(perc_plots))


#### pivot for nmds ####

# make a sample_id
species_percent2 <- species_percent %>%
  mutate(sample_id = paste0(year, "-",veg_type,"_sec", vt_section))

# 1) metadata: one row per sample (year x vt_section)
metadata <- species_percent2 %>%
  distinct(sample_id, year, vt_section, veg_type, n_plots_total) %>%
  arrange(year, vt_section)

# 2) species matrix wide: rows = sample_id, cols = species
species_wide <- species_percent2 %>%
  select(sample_id, species, perc_plots) %>%
  pivot_wider(
    names_from = species,
    values_from = perc_plots,
    values_fill = 0
  ) %>%
  arrange(sample_id)

# convert to matrix and set rownames
species_mat <- species_wide %>%
  column_to_rownames(var = "sample_id") %>%   # needs tibble::column_to_rownames or use base R
  as.matrix()

# ensure numeric (NAs filled by values_fill above, but double-check)
mode(species_mat) <- "numeric"

# quick checks
dim(species_mat)         # samples x species
head(metadata)

