##################################################### loading packages

library(tidyverse)
library(lme4)
library(lmerTest)
library(ggplot2)
library(ggeffects)

##################################################### loading data
merged_data <- readRDS("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/01 Vegetation changes Kobbefjord/data/nero_analysis/data/merged_data.rds")

names(merged_data)

richness.df <- merged_data |> 
  group_by(year, vt_section, plot, plot_id, veg_type) |> 
  summarise(richness = n_distinct(taxon_code), .groups = "drop")

names(richness.df)

##################################################### mixed linear modelling


m1 <- lmer(richness ~ year + (1|plot_id), data = richness.df)
summary(m1)



######################### centering year
# richness.df <- richness.df |> mutate(year_c = year - mean(year))
# m2 <- lmer(richness ~ year_c + (1|plot_id), data = richness.df)
# summary(m2)
# 
# ######################### all plots increase/decrease equally, allow random slopes
# m3 <- lmer(richness ~ year_c + (year_c|plot_id), data = richness.df)
# summary(m3)

##################################################### visualisering


# Get model predictions (with CI)
pred_df <- ggeffects::ggpredict(m1, terms = "year")

# Plot
ggplot(richness.df, aes(x = year, y = richness)) +
  geom_jitter(aes(group = plot_id), width = 0.3, alpha = 0.2, color = "#01ad7f") +
  geom_line(
    data = as.data.frame(pred_df),
    aes(x = x, y = predicted),
    color = "#004d38",
    size = 1.2,
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
    subtitle = paste0("Linear mixed model: p = ", 
                      formatC(summary(m1)$coefficients["year", "Pr(>|t|)"], format = "e", digits = 2))
  )
