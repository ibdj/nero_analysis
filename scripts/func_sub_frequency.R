#### packages ########################################################################################################################

library(tidyverse)
library(ggplot2)
library(purrr)
library(lme4)
library(broom.mixed)
library(emmeans)

#### import data ######################################################################################################################
# import data 

merged_data <- readRDS("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/01 Vegetation changes Kobbefjord/data/nero_analysis/data/merged_data.rds") |> 
  group_by(year, subsection) |> 
  mutate(no_plots = n_distinct(plot_id)) |> 
  ungroup() |> 
  filter(taxon_code != "rock")

names(merged_data)
#### adjust data ###################################################################################################################################

subsection_agreg_func_type <- merged_data |>
  group_by(year, subsection, ecoveg_sgfc, no_plots) |>
  summarise(frec = n(), .groups = "drop") |>
  mutate(frec_per_subsec = frec / no_plots)

#### visualising ####################################################################################################################################

# Recalculate summary statistics based on normalized frequency
summary_stats <- subsection_agreg_func_type |>
  group_by(year, ecoveg_sgfc) |>
  summarise(
    mean_frec = mean(frec_per_subsec),
    median_frec = median(frec_per_subsec),
    .groups = "drop"
  )

# Plot using normalized frequency per plot
ggplot(subsection_agreg_func_type, aes(x = frec_per_subsec, color = factor(year))) +
  geom_histogram(aes(fill = factor(year)), binwidth = 0.1, position = "dodge", alpha = 0.4) +
  geom_freqpoly(binwidth = 0.1, size = 1) +  # line outline of frequencies
  geom_vline(data = summary_stats, aes(xintercept = mean_frec, color = factor(year)),
             linetype = "dashed", size = 1) +
  facet_wrap(~ ecoveg_sgfc, scales = "free_y") +
  labs(title = "Normalized Frequency Distribution with Line Outlines and Mean Frequencies",
       x = "Frequency per Plot (frec_per_subsec)",
       y = "Number of Occurrences",
       fill = "Year",
       color = "Year") +
  theme_minimal()


#### lmm for fpts #########
models_by_group <- subsection_agreg_func_type |>
  group_by(ecoveg_sgfc) |>
  group_split() |>
  map(~ lmer(frec_per_subsec ~ year + (1 | subsection), data = .x))

# Extract summaries
map(models_by_group, summary)


#### visualisering lmms #####################################################

# Extract fixed effects estimates and t-values and determine significance per model
model_summaries <- map_df(models_by_group, ~ {
  s <- summary(.x)
  data.frame(
    ecoveg_sgfc = NA,  # will overwrite later
    intercept = s$coefficients[1, 1],
    year_estimate = s$coefficients[2, 1],
    year_tvalue = s$coefficients[2, 3],
    significant = abs(s$coefficients[2, 3]) > 2
  )
})

# Add functional group names to model summaries (in same order)
model_summaries$ecoveg_sgfc <- unique(subsection_agreg_func_type$ecoveg_sgfc)

# Prepare prediction  range of years per functional group
pred_data <- subsection_agreg_func_type |>
  group_by(ecoveg_sgfc) |>
  summarise(min_year = min(year), max_year = max(year)) |>
  rowwise() |>
  mutate(year = list(seq(min_year, max_year, by = 1))) |>
  unnest(year)

# Add model predictions with fixed effects and line type (significance)
predictions <- model_summaries |>
  left_join(pred_data, by = "ecoveg_sgfc") |>
  mutate(
    pred_frec = intercept + year_estimate * year,
    line_type = ifelse(significant, "solid", "dashed")
  )

# Plot observed points and predicted trend lines with line type mapped
ggplot(subsection_agreg_func_type, aes(x = year, y = frec_per_subsec)) +  # remove color = factor(year)
  geom_boxplot(alpha = 0.4, outlier.shape = NA, aes(group = factor(year))) +  # group boxplot by year
  geom_jitter(width = 0.1, alpha = 0.3, color = "black") +  # points in black without color grouping
  geom_line(data = predictions, aes(x = year, y = pred_frec, linetype = line_type), 
            size = 1.2, color = "black", inherit.aes = FALSE) +
  facet_wrap(~ ecoveg_sgfc, scales = "free_y") +
  scale_linetype_identity() +  # for solid/dashed lines
  labs(
    title = "LMM Predicted Trends with Significance (solid = significant)",
    x = "Year", y = "Frequency per Plot",
    linetype = "Significance"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")  # Adjust or remove legend as needed


#### lmm with year as factor ALL SHRUB ####

subsection_func_shrub <- merged_data |>
  group_by(year, subsection, ecoveg_gfc, no_plots) |>
  summarise(frec = n(), .groups = "drop") |>
  mutate(frec_per_subsec = frec / no_plots)|> 
  filter(ecoveg_gfc == "shrub")

data_model_shrub <- subsection_func_shrub

model_shrub <- lmer(data = data_model_shrub, frec_per_subsec ~ factor(year) + (1 | subsection))

# Extract summaries
model_shrub

model_shrub_f <- lmer(data = data_model_shrub, frec_per_subsec ~ factor(year) + (1 | subsection))

# Extract summaries
summary(model_shrub_f)

shrub_f_emmeans <- emmeans(model_shrub_f, ~ year) |>
  pairs()

shrub_f_emmeans

#### all shrub model visualisation ####

# Plot observed points and predicted trend lines with line type mapped
ggplot(subsection_func_shrub, aes(x = year, y = frec_per_subsec)) +  # remove color = factor(year)
  geom_boxplot(alpha = 0.4, outlier.shape = NA, aes(group = factor(year))) + # group boxplot by year
  geom_jitter(width = 0.1, alpha = 0.3, color = "darkgreen") + # points in black without color grouping
  theme_minimal()+
  scale_x_continuous(breaks = c(2007, 2012, 2017, 2022)) +
  labs(
    title = "LMM, year as factor",
    x = "Year", y = "Frequency per subsection"
  )

emm <- emmeans(model_shrub_f, ~ year)
emm_df <- as.data.frame(emm)

ggplot(emm_df, aes(x = year, y = emmean)) +
  geom_point(size = 2, color = "darkgreen") +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.3, color = "darkgreen") +
  labs(
    x = "Year",
    y = "Estimated fractional shrub abundance"
  ) +
  scale_x_continuous(breaks = c(2007, 2012, 2017,2022)) +
  theme_classic()

#### all shrub viasualising both the model and original data #####

# Ensure year is treated consistently
emm_df$year <- factor(emm_df$year)
subsection_agreg_func_type$year <- factor(subsection_agreg_func_type$year)

ggplot(emm_df, aes(x = year, y = emmean)) +
  
  # Raw data FIRST (background layer)
  geom_jitter(
    data = subsection_agreg_func_type,
    aes(x = year, y = frec_per_subsec),
    width = 0.15,
    alpha = 0.25,
    color = "black"
  ) +
  
  # Model estimates (foreground)
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.25,
    color = "darkgreen",
    linewidth = 0.8
  ) +
  geom_point(
    size = 2.5,
    color = "darkgreen"
  ) +
  
  labs(
    x = "Year",
    y = "Model-estimated fractional shrub abundance"
  ) +
  
  theme_classic()


#### ONLY decidious shrub ######################################################

subsection_func_d_shrub <- merged_data |>
  group_by(year, subsection, ecoveg_sgfc, no_plots) |>
  summarise(frec = n(), .groups = "drop") |>
  mutate(frec_per_subsec = frec / no_plots)|> 
  filter(ecoveg_sgfc == "shrub_decidous")

data_model_d_shrub <- subsection_func_d_shrub

model_d_shrub_f <- lmer(data = data_model_d_shrub, frec_per_subsec ~ factor(year) + (1 | subsection))

# Extract summaries
summary(model_d_shrub_f)

d_shrub_f_emmeans <- emmeans(model_d_shrub_f, ~ year) |>
  pairs()

d_shrub_f_emmeans


#### only dicidous shrub visualisation ####
emm_dshrub <- emmeans(model_d_shrub_f, ~ year)
emm_df_dshrub <- as.data.frame(emm_dshrub)

ggplot(emm_df_dshrub, aes(x = year, y = emmean)) +
  geom_point(size = 2, color = "darkgreen") +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.3, color = "darkgreen") +
  labs(
    x = "Year",
    y = "Estimated fractional shrub abundance"
  ) +
  scale_x_continuous(breaks = c(2007, 2012, 2017,2022)) +
  theme_classic()
