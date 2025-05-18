#### packages #### 
library(vegan)
library(ggplot2)


#### subsections functional types data prepping ####


func_type_subsection_summary <- readRDS("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/01 Vegetation changes Kobbefjord/data/nmds_nero/nmds_nero/data/func_type_subsection_summary.rds")

names(func_type_subsection_summary)


func_subsection_wide <- func_type_subsection_summary |> 
  pivot_wider(names_from = ecoveg_gfc, values_from = func_plot_fraction, values_fill = 0)

head(func_subsection_wide)

#### relative abundance of functional types copse #####

unique_years <- sort(unique(func_type_subsection_summary$year))

ggplot(func_type_subsection_summary |> filter(veg_type == "copse"), 
       aes(x = year, y = func_plot_fraction, color = ecoveg_gfc)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_continuous(breaks = unique_years)

#### relative abundance of functional types heath #####

ggplot(func_type_subsection_summary |> filter(veg_type == "heath"), 
       aes(x = year, y = func_plot_fraction, color = ecoveg_gfc)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_continuous(breaks = unique_years)

#### relative abundance of functional types fen #####

ggplot(func_type_subsection_summary |> filter(veg_type == "fen"), 
       aes(x = year, y = func_plot_fraction, color = ecoveg_gfc)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_continuous(breaks = unique_years)

#### relative abundance of functional types copse #####

ggplot(func_type_subsection_summary |> filter(veg_type == "copse"), 
       aes(x = year, y = func_plot_fraction, color = ecoveg_gfc)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_continuous(breaks = unique_years)

#### slope for all veg_types and growth forms #####
library(dplyr)
library(broom)

results <- func_type_subsection_summary %>%
  group_by(veg_type, ecoveg_gfc) %>%
  do(tidy(lm(func_plot_fraction ~ year, data = .))) %>%
  filter(term == "year") %>%
  mutate(significant = p.value < 0.05)

view(results)
#### nmds meta data for nmds ####

# Extract species data (gramnoid, herb, etc.)
species_data <- func_subsection_wide[, 6:11]

# Metadata for grouping (year and veg_type)
metadata <- func_subsection_wide[, c("year", "veg_type")]

#### nmds ana ####

set.seed(123)  # For reproducibility
nmds_result <- metaMDS(species_data, distance = "bray", trymax = 100)

#### nmds results #####

# For newer vegan versions (>2.6-2):
scores_df <- as.data.frame(scores(nmds_result)$sites)
scores_df$year <- metadata$year
scores_df$veg_type <- metadata$veg_type
