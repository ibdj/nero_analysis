#######################################################################################################################################

library(tidyverse)
library(ggplot2)
library(purrr)
library(lme4)
library(broom.mixed)

#######################################################################################################################################
# import data 

merged_data <- readRDS("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/01 Vegetation changes Kobbefjord/data/nero_analysis/data/merged_data.rds") |> 
  group_by(year, subsection) |> 
  mutate(no_plots = n_distinct(plot_id)) |> 
  ungroup()

names(merged_data)
#######################################################################################################################################
# adjust data

subsection_agreg_func_type <- merged_data |>
  group_by(year, subsection, ecoveg_sgfc, no_plots) |>
  summarise(frec = n(), .groups = "drop") |>
  mutate(frec_per_plot = frec / no_plots)

#######################################################################################################################################
# visualising

# Recalculate summary statistics based on normalized frequency
summary_stats <- subsection_agreg_func_type |>
  group_by(year, ecoveg_sgfc) |>
  summarise(
    mean_frec = mean(frec_per_plot),
    median_frec = median(frec_per_plot),
    .groups = "drop"
  )


# Plot using normalized frequency per plot
ggplot(subsection_agreg_func_type, aes(x = frec_per_plot, color = factor(year))) +
  geom_histogram(aes(fill = factor(year)), binwidth = 0.1, position = "dodge", alpha = 0.4) +
  geom_freqpoly(binwidth = 0.1, size = 1) +  # line outline of frequencies
  geom_vline(data = summary_stats, aes(xintercept = mean_frec, color = factor(year)),
             linetype = "dashed", size = 1) +
  facet_wrap(~ ecoveg_sgfc, scales = "free_y") +
  labs(title = "Normalized Frequency Distribution with Line Outlines and Mean Frequencies",
       x = "Frequency per Plot (frec_per_plot)",
       y = "Number of Occurrences",
       fill = "Year",
       color = "Year") +
  theme_minimal()


#############
models_by_group <- subsection_agreg_func_type |>
  group_by(ecoveg_sgfc) |>
  group_split() |>
  map(~ lmer(frec_per_plot ~ year + (1 | subsection), data = .x))

# Extract summaries
map(models_by_group, summary)


#########################################################

library(dplyr)
library(purrr)
library(ggplot2)
library(lme4)

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
ggplot(subsection_agreg_func_type, aes(x = year, y = frec_per_plot)) +  # remove color = factor(year)
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

