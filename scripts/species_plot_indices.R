##################################################### loading packages

library(tidyverse)
library(lme4)
library(lmerTest)
library(ggplot2)
library(ggeffects)
library(vegan)

######################### loading data ################################
species_plot <- readRDS("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/01 Vegetation changes Kobbefjord/data/nero_analysis/data/merged_data.rds") |> 
  select(year, vt_section, plot, plot_id, veg_type, taxon_code) |> 
  filter(taxon_code != "rock")

write_rds(species_plot, "data/species_plot.rds")

species_plot <- read_rds("data/species_plot.rds")

richness.df <- species_plot |> 
  group_by(year, vt_section, plot, plot_id, veg_type) |> 
  summarise(richness = n_distinct(taxon_code), .groups = "drop")

names(richness.df)

######################### mixed linear modelling species plot richness ####

m_richness <- lmer(richness ~ year + (1|plot_id), data = richness.df)
summary(m_richness)

######################### visualisering plot species richness ####


# Get model predictions (with CI)
pred_df <- ggeffects::ggpredict(m_richness, terms = "year")

shannon_df |>
  ggplot(aes(x = factor(year), y = H)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(
    width = 0.15,       # slight horizontal jitter
    alpha = 0.35,       # light transparency to avoid overplotting
    size = 1.2,         # point size
    color = "grey30"
  ) +
  labs(
    x = "Year",
    y = "Shannon diversity (H)",
    title = "Shannon diversity across years"
  ) +
  theme_minimal(base_size = 13)


######################### EVENNESS #########################

# evenness not applicable for pres/abs data
evenness_species <- species_plot |>
  # count abundance per taxon in each plot-year
  count(year, vt_section, plot, plot_id, veg_type, taxon_code, name = "abundance") |>
  
  # compute diversity per plot-year
  group_by(year, vt_section, plot, plot_id, veg_type) |>
  summarise(
    H = diversity(abundance, index = "shannon"),
    S = n(),                          # number of taxa
    J = ifelse(S > 1, H / log(S), NA_real_),
    .groups = "drop"
  )
head(evenness_species)

species_plot %>% 
  count(year, plot_id, taxon_code) %>% 
  count(year, plot_id) %>% 
  summarise(
    S = n(),
    H = diversity(rep(1, n()), index = "shannon"),
    J = H / log(S)
  )

######################### SHANNON #########################

valid_levels <- c("2007", "2012", "2017", "2022")

shannon_df <- species_plot |>
  group_by(year, plot_id, veg_type) |>
  summarise(
    S = n(),                                 # richness = count of taxa
    H = diversity(rep(1, S), index = "shannon"),  # Shannon index
    .groups = "drop"
  ) |> 
  filter(!is.na(year)) |> 
  mutate(year = factor(as.character(year), levels = valid_levels))

head(shannon_df)

m_shannon <- lmer(H ~ year + (1 | plot_id), data = shannon_df)
summary(m_shannon)

######################### SHANNON visualisation #########################

# Only predict at actual survey years
pred_shannon_plot <- pred_shannon_plot %>%
  filter(!is.na(x)) %>%
  mutate(year = factor(as.character(x), levels = valid_levels))

pred_shannon_plot <- pred_shannon |>
  mutate(year = factor(x))

shannon_df |>
  ggplot(aes(x = year, y = H)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.35, size = 1.1, color = "grey30") +
  geom_ribbon(
    data = pred_shannon_plot,
    aes(x = year, ymin = conf.low, ymax = conf.high, group = 1),
    inherit.aes = FALSE,
    fill = "blue",
    alpha = 0.15
  ) +
  geom_line(
    data = pred_shannon_plot,
    aes(x = year, y = predicted, group = 1),
    inherit.aes = FALSE,
    color = "blue",
    linewidth = 1.2
  ) +
  labs(
    x = "Year",
    y = "Shannon diversity (H)",
    title = "Shannon diversity across years (observed + model predictions)"
  ) +
  theme_minimal(base_size = 13)

######################### TURNOVER #########################


