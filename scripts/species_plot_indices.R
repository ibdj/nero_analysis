##################################################### loading packages

library(tidyverse)
library(lme4)
library(lmerTest)
library(ggplot2)
library(ggeffects)

######################### loading data ################################
merged_data <- readRDS("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/01 Vegetation changes Kobbefjord/data/nero_analysis/data/merged_data.rds")

names(merged_data)

richness.df <- merged_data |> 
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
  geom_jitter(aes(group = plot_id), width = 0.3, alpha = 0.2, color = "#01ad7f") +
  geom_boxplot(aes(group = factor(year)), outlier.shape = NA, fill = "lightblue", alpha = 0.5, color = "gray30", width = 0.7) +
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
  theme_bw(base_size = 13) +
  labs(
    x = "Year",
    y = "Species richness per plot",
    title = "Change in species richness over time",
    subtitle = paste0("Linear mixed model (p = ", 
                      formatC(summary(m_richness)$coefficients["year", "Pr(>|t|)"], format = "e", digits = 2),")")
  )



######################### EVENNESS #########################

######################### SHANNON #########################

######################### TURNOVER #########################


