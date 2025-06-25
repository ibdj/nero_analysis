#### packages #### 
library(vegan)
library(ggplot2)
library(tidyverse)

#### subsections functional types data prepping ####


func_type_subsection_summary <- readRDS("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/01 Vegetation changes Kobbefjord/data/nmds_nero/nmds_nero/data/func_type_subsection_summary.rds")

names(func_type_subsection_summary)


func_subsection_wide <- func_type_subsection_summary |> 
  pivot_wider(names_from = ecoveg_gfc, values_from = func_plot_fraction, values_fill = 0)

head(func_subsection_wide)

#### relative abundance of functional types copse #####

unique_years <- sort(unique(func_type_subsection_summary$year))

veg_types <- unique(func_type_subsection_summary$veg_type)

for (veg in veg_types) {
  p <- ggplot(func_type_subsection_summary |> filter(veg_type == veg), 
              aes(x = year, y = func_plot_fraction, color = ecoveg_gfc)) +
    geom_point() +
    geom_smooth(method = "lm") +
    scale_x_continuous(breaks = unique_years) +
    ggtitle(paste("veg_type:", veg))
  print(p)
}

### facet wrap ###

unique_years <- sort(unique(func_type_subsection_summary$year))

ggplot(func_type_subsection_summary, 
       aes(x = year, y = func_plot_fraction, color = ecoveg_gfc)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_continuous(breaks = unique_years) +
  facet_wrap(~veg_type)

### facet wrap with sigificance lines ####
plot_data <- func_type_subsection_summary %>%
  left_join(
    results %>% select(veg_type, ecoveg_gfc, significant),
    by = c("veg_type", "ecoveg_gfc")
  )

ggplot(plot_data, 
       aes(x = year, y = func_plot_fraction, color = ecoveg_gfc)) +
  geom_point() +
  geom_smooth(
    method = "lm",
    aes(linetype = significant),  # Map linetype to significance
    se = TRUE,                    # Keep confidence intervals if desired
    size = 0.7                    # Adjust line thickness
  ) +
  scale_x_continuous(breaks = unique_years) +
  scale_linetype_manual(
    name = "Significance",
    values = c("TRUE" = "solid", "FALSE" = "dashed"),  # Define line styles
    labels = c("TRUE" = "p < 0.05", "FALSE" = "NS")
  ) +
  facet_wrap(~veg_type) +
  theme_minimal()


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

#### looking at frequency of functional types ####
plot_data <- func_type_subsection_summary %>%
  mutate(year_fct = factor(year))

plot_data <- func_type_subsection_summary %>%
  mutate(year_fct = factor(year))

ggplot(plot_data, aes(x = func_plot_fraction, fill = year_fct)) +
  geom_histogram(
    position = position_dodge2(preserve = "single"),  # Align bins across years
    bins = 15,                                         # Adjust based on your data
    color = "white",                                   # Outline color
    alpha = 0.8                                       
  ) +
  facet_wrap(~ecoveg_gfc, scales = "free_y") +
  labs(
    x = "func_plot_fraction",
    y = "Frequency",
    fill = "Year",
    title = "Yearly func_plot_fraction Distributions by Growth Form"
  ) +
  scale_fill_viridis_d() +  # Colorblind-friendly palette
  theme_minimal()

#### frequency of functional groups pr vegetation type #####
# Get unique veg_types
veg_types <- unique(func_type_subsection_summary$veg_type)

# Loop through each veg_type
for (v in veg_types) {
  
  # Filter data for current veg_type
  plot_data <- func_type_subsection_summary %>%
    filter(veg_type == v) %>%
    mutate(year = factor(year))  # Convert year to factor for discrete grouping
  
  # Skip empty data (if any veg_type has no observations)
  if (nrow(plot_data) == 0) next
  
  # Create plot
  p <- ggplot(plot_data, aes(x = func_plot_fraction, fill = year)) +
    geom_histogram(
      position = position_dodge2(preserve = "single"),  # Align bins
      bins = 15,                       
      color = "white",                # Bar outline color
      alpha = 0.8
    ) +
    facet_wrap(~ecoveg_gfc, scales = "free_y") +  # Facet by growth form
    labs(
      title = paste("veg_type:", v),
      x = "func_plot_fraction",
      y = "Frequency",
      fill = "Year"
    ) +
    scale_fill_viridis_d() +  # Color palette
    theme_minimal()
}
  # Print plot (or save)
  print(p)
  