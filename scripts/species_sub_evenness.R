library(tidyverse)
library(vegan)
library(ggeffects)  # for model predictions


#### subsection evenness (Pielou) ####

# importing data #
merged_data <- readRDS("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/01 Vegetation changes Kobbefjord/data/nero_analysis/data/merged_data.rds") |> 
  group_by(year, subsection) |> 
  mutate(no_plots = n_distinct(plot_id)) |> 
  ungroup()

names(merged_data)

abundance.df <- merged_data |> 
  group_by(year, subsection, veg_type, taxon_code) |> 
  reframe(no_plots = mean(no_plots),
    count = n(),
            frac = count/no_plots) 
            

# Generate model predictions
pred_even <- ggpredict(m_even, terms = "year")

evenness_df |>
  ggplot(aes(x = year, y = J)) +
  geom_point(alpha = 0.4, color = "grey30") +
  geom_line(aes(group = subsection), alpha = 0.2, color = "grey60") +
  geom_line(
    data = pred_even,
    aes(x = x, y = predicted),
    color = "green",
    linewidth = 1.2,
    inherit.aes = FALSE
  ) +
  geom_ribbon(
    data = pred_even,
    aes(x = x, ymin = conf.low, ymax = conf.high),
    alpha = 0.2,
    fill = "green",
    inherit.aes = FALSE
  ) +
  labs(
    x = "Year",
    y = "Pielou's evenness (J)",
    title = "Temporal trend in evenness across plots"
  ) +
  scale_x_continuous(
    breaks = c(2007, 2012, 2017, 2022),
    labels = c("2007", "2012", "2017", "2022")
  ) +
  theme_minimal(base_size = 13)
library(lme4)

m_even <- lmer(J ~ year + (1 | subsection), data = evenness_df)
summary(m_even)

plot(m_even)       # residuals vs fitted
qqnorm(resid(m_even))
qqline(resid(m_even))

m_even2 <- lmer(J ~ year * veg_type + (1 | subsection), data = evenness_df)
summary(m_even2)

anova(m_even, m_even2)


########################### library(lme4)
library(ggeffects)

m_shannon <- lmer(H ~ year + (1 | subsection), data = evenness_df)
summary(m_shannon)

pred_shannon <- ggpredict(m_shannon, terms = "year")

evenness_df |>
  ggplot(aes(x = year, y = H)) +
  geom_point(alpha = 0.4, color = "grey30") +
  #geom_line(aes(group = subsection), alpha = 0.2, color = "grey60") +
  geom_line(data = pred_shannon, aes(x = x, y = predicted), color = "blue", linewidth = 0.8,linetype = "dashed"  ) +
  geom_ribbon(
    data = pred_shannon,
    aes(x = x, ymin = conf.low, ymax = conf.high, y = predicted),
    alpha = 0.2, fill = "blue", inherit.aes = FALSE
  )+
  scale_x_continuous(
    breaks = c(2007, 2012, 2017, 2022),
    labels = c("2007", "2012", "2017", "2022")
  ) +
  labs(
    x = "Year",
    y = "Shannon diversity (H')",
    title = "Temporal trend in Shannon diversity across plots"
  ) +
  theme_minimal(base_size = 13)
