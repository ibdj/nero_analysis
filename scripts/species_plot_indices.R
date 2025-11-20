##################################################### loading packages

library(tidyverse)
library(lme4)
library(lmerTest)
library(ggplot2)
library(ggeffects)
library(vegan)
library(codyn) # turnover calculations

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

# Plot
ggplot(richness.df, aes(x = year, y = richness)) +
  geom_jitter(aes(group = plot_id), width = 0.2, alpha = 0.2, color = "#01ad7f") +
  geom_boxplot(aes(group = factor(year)), outlier.shape = NA, alpha = 0.5, color = "gray30", width = 0.6) +
  geom_line(
    data = as.data.frame(pred_df),
    aes(x = x, y = predicted),
    color = "#004d38",
    linewidth = 1.2,
    inherit.aes = FALSE
  ) +
  geom_ribbon(
    data = as.data.frame(pred_df),
    aes(x = x, ymin = conf.low, ymax = conf.high),
    fill = "#004d38",
    alpha = 0.2,
    inherit.aes = FALSE
  ) +
  scale_x_continuous(breaks = c(2007, 2012, 2017, 2022)) +
#  theme_bw(base_size = 13) +
  labs(
    x = "Year",
    y = "Species richness per plot",
    title = "Change in species richness over time",
    subtitle = paste0("Linear mixed model (p = ", 
                      formatC(summary(m_richness)$coefficients["year", "Pr(>|t|)"], format = "e", digits = 2),")")
    )


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



shannon_df <- species_plot |>
  group_by(year, plot_id, veg_type) |>
  summarise(
    S = n(),                                 # richness = count of taxa
    H = diversity(rep(1, S), index = "shannon"),  # Shannon index
    .groups = "drop"
  ) |> 
  filter(!is.na(year))

head(shannon_df)

m_shannon <- lmer(H ~ year + (1 | plot_id), data = shannon_df)
summary(m_shannon)

######################### SHANNON visualisation #########################

# Recreate predicted values only at the actual survey years
pred_shannon <- ggpredict(
  m_shannon,
  terms = c("year [2007,2012,2017,2022]")  # specify the years you want predictions for
)

# Convert x to factor for plotting
pred_shannon_plot <- pred_shannon %>% 
  filter(!is.na(x)) %>%
  mutate(x = factor(x))

shannon_df |>
  ggplot(aes(x = factor(year), y = H)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.7) +
  geom_jitter(
    width = 0.1,
    alpha = 0.35,
    size = 1.1,
    color = "grey30"
  ) +
  geom_ribbon(
    data = pred_shannon_plot,
    aes(x = factor(x), ymin = conf.low, ymax = conf.high, group = 1),
    inherit.aes = FALSE,
    fill = "blue",
    alpha = 0.15
  ) +
  geom_line(
    data = pred_shannon_plot,
    aes(x = factor(x), y = predicted, group = 1),
    inherit.aes = FALSE,
    color = "blue",
    linewidth = 1.2
  ) +
  scale_x_discrete(na.translate = FALSE) +   # <- hide NA tick
  labs(
    x = "Year",
    y = "Shannon diversity (H)",
    title = "Shannon diversity across years (observed + model predictions)"
  ) +
  theme_minimal(base_size = 13)


######################### TURNOVER #########################

# Prepare data with native pipe and remove duplicates
species_long <- species_plot |>
  select(year, plot_id, taxon_code) |>
  mutate(abundance = 1) |>
  distinct()

# Check number of unique years per plot
plot_years <- species_long |>
  group_by(plot_id) |>
  summarise(n_years = n_distinct(year), .groups = "drop")

# Filter plots with at least 2 years (turnover requires temporal comparison)
valid_plots <- plot_years |>
  filter(n_years > 1) |>
  pull(plot_id)

# Filter data to valid plots only
species_long_filtered <- species_long |>
  filter(plot_id %in% valid_plots)

# Calculate turnover per plot between years
turnover_results <- turnover(
  df = species_long_filtered,
  time.var = "year",
  species.var = "taxon_code",
  abundance.var = "abundance",
  replicate.var = "plot_id",
  metric = "total"
)

head(turnover_results)


######################### TURNOVER visualisation ####
# 2. Fit linear mixed model: turnover by year with random intercept for plot
m_turnover <- lmer(total ~ year + (1 | plot_id), data = turnover_results)

# 3. Generate predicted values at observed years
pred_turnover <- ggpredict(m_turnover, terms = c("year"))

# 4. Prepare predicted data for plotting
pred_turnover_plot <- pred_turnover |>
  mutate(year = as.numeric(x))  # convert x (character) to numeric for plotting

# 5. Boxplot + jitter for observed turnover, line + ribbon for predicted
ggplot(turnover_results, aes(x = factor(year), y = total)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.35, size = 1.1, color = "grey30") +
  geom_ribbon(
    data = pred_turnover_plot,
    aes(x = factor(year), ymin = conf.low, ymax = conf.high, group = 1),
    inherit.aes = FALSE,
    fill = "blue",
    alpha = 0.15
  ) +
  geom_line(
    data = pred_turnover_plot,
    aes(x = factor(year), y = predicted, group = 1),
    inherit.aes = FALSE,
    color = "blue",
    linewidth = 1.2
  ) +
  labs(
    x = "Year",
    y = "Turnover",
    title = "Species turnover across years (observed + predicted)"
  ) +
  theme_minimal(base_size = 13)

######################### combined plots ####

library(dplyr)
library(ggplot2)
library(ggeffects)

# 1. Prepare predicted data frames with numeric year and original scale values

# Richness predictions
pred_richness <- ggpredict(m_richness, terms = "year") |>
  as.data.frame() |>
  mutate(year = as.numeric(x),
         metric = "Richness",
         value = predicted,
         conf.low = conf.low,
         conf.high = conf.high) |>
  select(year, value, conf.low, conf.high, metric)

# Shannon predictions
pred_shannon <- ggpredict(m_shannon, terms = c("year [2007,2012,2017,2022]")) |>
  as.data.frame() |>
  mutate(year = as.numeric(x),
         metric = "Shannon",
         value = predicted,
         conf.low = conf.low,
         conf.high = conf.high) |>
  select(year, value, conf.low, conf.high, metric)

# Turnover predictions
pred_turnover <- ggpredict(m_turnover, terms = c("year")) |>
  as.data.frame() |>
  mutate(year = as.numeric(x),
         metric = "Turnover",
         value = predicted,
         conf.low = conf.low,
         conf.high = conf.high) |>
  select(year, value, conf.low, conf.high, metric)

# Combine predicted data
pred_all <- bind_rows(pred_richness, pred_shannon, pred_turnover)

# 2. Prepare observed data frames with scaling

# helper that returns both scaled values and the parameters used
scale_with_params <- function(x) {
  xmin <- min(x, na.rm = TRUE)
  xmax <- max(x, na.rm = TRUE)
  list(
    scaled = (x - xmin) / (xmax - xmin),
    min = xmin,
    max = xmax
  )
}

### 1. SCALE OBSERVED DATA BY METRIC, STORING MIN/MAX

## Richness
rich_scale <- scale_with_params(richness.df$richness)
richness_obs_scaled <- richness.df |>
  mutate(
    value  = rich_scale$scaled,
    metric = "Richness"
  ) |>
  select(year, plot_id, value, metric)

## Shannon
shan_scale <- scale_with_params(shannon_df$H)
shannon_obs_scaled <- shannon_df |>
  mutate(
    value  = shan_scale$scaled,
    metric = "Shannon"
  ) |>
  select(year, plot_id, value, metric)

## Turnover
turn_scale <- scale_with_params(turnover_results$total)
turnover_obs_scaled <- turnover_results |>
  mutate(
    value  = turn_scale$scaled,
    metric = "Turnover"
  ) |>
  select(year, plot_id, value, metric) |> 
  mutate(plot_id = as.double(plot_id))

observed_all <- bind_rows(richness_obs_scaled,
                          shannon_obs_scaled,
                          turnover_obs_scaled)

### 2. SCALE PREDICTIONS USING THE SAME MIN/MAX PER METRIC

## Richness preds
pred_richness <- ggpredict(m_richness, terms = "year") |>
  as.data.frame() |>
  mutate(
    year  = as.numeric(x),
    metric = "Richness",
    value = (predicted - rich_scale$min) /
      (rich_scale$max - rich_scale$min),
    conf.low  = (conf.low - rich_scale$min) /
      (rich_scale$max - rich_scale$min),
    conf.high = (conf.high - rich_scale$min) /
      (rich_scale$max - rich_scale$min)
  ) |>
  select(year, value, conf.low, conf.high, metric)

## Shannon preds
pred_shannon <- ggpredict(m_shannon, terms = "year") |>
  as.data.frame() |>
  mutate(
    year  = as.numeric(x),
    metric = "Shannon",
    value = (predicted - shan_scale$min) /
      (shan_scale$max - shan_scale$min),
    conf.low  = (conf.low - shan_scale$min) /
      (shan_scale$max - shan_scale$min),
    conf.high = (conf.high - shan_scale$min) /
      (shan_scale$max - shan_scale$min)
  ) |>
  select(year, value, conf.low, conf.high, metric)

## Turnover preds
pred_turnover <- ggpredict(m_turnover, terms = "year") |>
  as.data.frame() |>
  mutate(
    year  = as.numeric(x),
    metric = "Turnover",
    value = (predicted - turn_scale$min) /
      (turn_scale$max - turn_scale$min),
    conf.low  = (conf.low - turn_scale$min) /
      (turn_scale$max - turn_scale$min),
    conf.high = (conf.high - turn_scale$min) /
      (turn_scale$max - turn_scale$min)
  ) |>
  select(year, value, conf.low, conf.high, metric)

pred_all <- bind_rows(pred_richness, pred_shannon, pred_turnover)

### 3. PLOT
pos_box  <- position_dodge(width = 0.8)
pos_jitt <- position_jitterdodge(
  jitter.width = 0.15,
  jitter.height = 0,
  dodge.width  = 0.8
)

ggplot() +
  # one box per (year × metric)
  geom_boxplot(
    data = observed_all,
    aes(
      x     = factor(year),
      y     = value,
      fill  = metric,
      group = interaction(metric, year)
    ),
    alpha         = 0.5,
    outlier.shape = NA,
    position      = pos_box
  ) +
  # points aligned with their own box
  geom_point(
    data = observed_all,
    aes(
      x     = factor(year),
      y     = value,
      color = metric
    ),
    size     = 1,
    alpha    = 0.35,
    position = pos_jitt
  ) +
  # predictions (no dodging: one line per metric)
  geom_ribbon(
    data = pred_all,
    aes(
      x     = factor(year),
      ymin  = conf.low,
      ymax  = conf.high,
      fill  = metric,
      group = metric
    ),
    alpha = 0.15
  ) +
  geom_line(
    data = pred_all,
    aes(
      x     = factor(year),
      y     = value,
      color = metric,
      group = metric
    ),
    linewidth = 1.1
  ) +
  scale_x_discrete(breaks = c("2007","2012","2017","2022")) +
  labs(
    x     = "Year",
    y     = "Scaled value (0–1, within metric)",
    title = "Relative change in richness, Shannon diversity and turnover over time (pr plot)",
    fill  = "Metric",
    color = "Metric"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")
