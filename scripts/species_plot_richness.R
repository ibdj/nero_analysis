library(ggplot2)
library(ggpmisc)
library(tidyverse)
library(dplyr)
library(broom)
library(readxl)

#### diversity, reading data and mutating ####

merged_data <- readRDS("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/01 Vegetation changes Kobbefjord/data/nmds_nero/nmds_nero/data/merged_data.rds")

View(merged_data)

df <- merged_data |>  
  group_by(year, date, plot_id, veg_type) |> 
  summarize(unique_taxa = n_distinct(taxon_code))

mean_species_pr_plot <- df |> 
  group_by(year) |> 
  summarise(mean_no_species = mean(unique_taxa))

#### stats and labels needed for the plotting ####
all_years <- unique(df$year)

label_df <- df |> 
  group_by(veg_type) |> 
  summarise(
    n = n(),
    x = max(year, na.rm = TRUE),
    y = max(unique_taxa, na.rm = TRUE)
  )
label_df$n_label <- paste0("n = ", label_df$n)

model_stats <- df |>
  group_by(veg_type) |> 
  summarise(
    model = list(lm(unique_taxa ~ year)),
    .groups = 'drop'
  ) |> 
  rowwise() |> 
  mutate(
    slope_p = tidy(model)$p.value[2],  # Slope p-value
    linetype = ifelse(slope_p < 0.05, "solid", "dashed")
  ) |> 
  select(veg_type, linetype)

# 3. MERGE WITH PLOT DATA
plot_data <- df |>  
  left_join(model_stats, by = "veg_type")

# 4. CREATE N-LABELS
label_df <- df |> 
  group_by(veg_type) |> 
  summarise(n = n()) |> 
  mutate(n_label = paste0("n = ", n))

#### plot  ####

ggplot(plot_data, aes(x = year, y = unique_taxa, color = veg_type)) +
  geom_point(alpha = 0.5) +
  geom_smooth(
    method = "lm",
    aes(linetype = linetype),
    se = TRUE,
    size = 0.8
  ) +
  scale_linetype_identity() +
  stat_poly_eq(
    aes(label = paste0(..p.value.label..,"~", ..rr.label..,"~",..eq.label..)),
    formula = y ~ x,
    parse = TRUE,
    label.y = 0.95,
    label.x = 0.05,
    vjust = 1,
    hjust = 0,
    size = 3
  ) +
  facet_wrap(~veg_type) +
  geom_text(
    data = label_df,
    aes(x = Inf, y = Inf, label = n_label),
    inherit.aes = FALSE,
    hjust = 1.2,
    vjust = 1.2,
    size = 3
  ) +
  scale_x_continuous(breaks = all_years) +
  coord_cartesian(clip = "off") +
  theme(legend.position = "none")  # This removes all legends


#### all vegetation type together #####

ggplot(df, aes(x = year, y = unique_taxa))+
  geom_point(size = 3, color = "darkgreen", alpha = 0.5)+
  geom_smooth(method = "lm")+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~")),
    formula = y ~ x,
    parse = TRUE
  )+
  scale_x_continuous(breaks = all_years)+
  labs(x = "", y = "Number of species pr plot")+
  geom_point(data = mean_species_pr_plot, aes(x = year, y = mean_no_species), shape = 17, size = 5, alpha = 0.5)
  
#### end of script ####


#### make the plot only for non-shrub species ####

df_non_shrub <- merged_data |>  
  filter(func_type != "shrub") |>
  group_by(year, date, plot_id, veg_type) |>
  summarize(unique_taxa = n_distinct(taxon_code))

all_years_non_shrub <- unique(df_non_shrub$year)

label_df_non_shrub <- df_non_shrub  |> 
  group_by(veg_type)  |> 
  summarise(
    n = n(),
    x = max(year, na.rm = TRUE),
    y = max(unique_taxa, na.rm = TRUE)
  )
label_df_non_shrub$n_label <- paste0("n = ", label_df$n)

model_stats_non_shrub <- df_non_shrub |> 
  group_by(veg_type) |> 
  summarise(
    model = list(lm(unique_taxa ~ year)),
    .groups = 'drop'
  ) |> 
  rowwise() |> 
  mutate(
    slope_p = tidy(model)$p.value[2],  # Slope p-value
    linetype = ifelse(slope_p < 0.05, "solid", "dashed")
  ) |> 
  select(veg_type, linetype)

# 3. MERGE WITH PLOT DATA
plot_data_non_shrub <- df_non_shrub |>  
  left_join(model_stats_non_shrub, by = "veg_type")

# 4. CREATE N-LABELS
label_df_non_shrub <- df_non_shrub |> 
  group_by(veg_type) |> 
  summarise(n = n()) |> 
  mutate(n_label = paste0("n = ", n))

# plot  #

ggplot(plot_data_non_shrub, aes(x = year, y = unique_taxa, color = veg_type)) +
  geom_point(alpha = 0.5) +
  geom_smooth(
    method = "lm",
    aes(linetype = linetype),
    se = TRUE,
    size = 0.8
  ) +
  scale_linetype_identity() +
  stat_poly_eq(
    aes(label = paste0(..p.value.label..,"~", ..rr.label..,"~",..eq.label..)),
    formula = y ~ x,
    parse = TRUE,
    label.y = 0.95,
    label.x = 0.05,
    vjust = 1,
    hjust = 0,
    size = 3
  ) +
  facet_wrap(~veg_type) +
  geom_text(
    data = label_df_non_shrub,
    aes(x = Inf, y = Inf, label = n_label),
    inherit.aes = FALSE,
    hjust = 1.2,
    vjust = 1.2,
    size = 3
  ) +
  scale_x_continuous(breaks = all_years_non_shrub) +
  coord_cartesian(clip = "off") +
  theme(legend.position = "none")  # This removes all legends

#### plotting all veg types with only shrub species #####

df_shrub <- merged_data |>  
  filter(func_type == "shrub") |>
  group_by(year, date, plot_id, veg_type) |>
  summarize(unique_taxa = n_distinct(taxon_code))

all_years_shrub <- unique(df_shrub$year)

label_df_shrub <- df_shrub  |> 
  group_by(veg_type)  |> 
  summarise(
    n = n(),
    x = max(year, na.rm = TRUE),
    y = max(unique_taxa, na.rm = TRUE)
  )
label_df_shrub$n_label <- paste0("n = ", label_df$n)

model_stats_shrub <- df_shrub |> 
  group_by(veg_type) |> 
  summarise(
    model = list(lm(unique_taxa ~ year)),
    .groups = 'drop'
  ) |> 
  rowwise() |> 
  mutate(
    slope_p = tidy(model)$p.value[2],  # Slope p-value
    linetype = ifelse(slope_p < 0.05, "solid", "dashed")
  ) |> 
  select(veg_type, linetype)

# 3. MERGE WITH PLOT DATA
plot_data_shrub <- df_shrub |>  
  left_join(model_stats_shrub, by = "veg_type")

# 4. CREATE N-LABELS
label_df_shrub <- df_shrub |> 
  group_by(veg_type) |> 
  summarise(n = n()) |> 
  mutate(n_label = paste0("n = ", n))

# plot  #

ggplot(plot_data_shrub, aes(x = year, y = unique_taxa, color = veg_type)) +
  geom_point(alpha = 0.5) +
  geom_smooth(
    method = "lm",
    aes(linetype = linetype),
    se = TRUE,
    size = 0.8
  ) +
  scale_linetype_identity() +
  stat_poly_eq(
    aes(label = paste0(..p.value.label..,"~", ..rr.label..,"~",..eq.label..)),
    formula = y ~ x,
    parse = TRUE,
    label.y = 0.95,
    label.x = 0.05,
    vjust = 1,
    hjust = 0,
    size = 3
  ) +
  facet_wrap(~veg_type) +
  geom_text(
    data = label_df_shrub,
    aes(x = Inf, y = Inf, label = n_label),
    inherit.aes = FALSE,
    hjust = 1.2,
    vjust = 1.2,
    size = 3
  ) +
  scale_x_continuous(breaks = all_years_shrub) +
  coord_cartesian(clip = "off") +
  theme(legend.position = "none")  # This removes all legends

