#### loading packages ####

library(tidyverse)
library(lme4)
library(lmerTest)
library(ggplot2)
library(ggeffects)
library(vegan)
library(codyn) # turnover calculations

#### loading data ################################

merged_data <- readRDS("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/01 Vegetation changes Kobbefjord/data/nero_analysis/data/merged_data.rds") |> 
  group_by(year, vt_section) |> 
  mutate(no_plots = n_distinct(plot_id)) |> 
  rename(section = vt_section) |> 
  ungroup()

species_sec_long <- merged_data |>
  group_by(year, section, veg_type, taxon_code, no_plots) |>
  summarize(
    count = n(),
    fraction_sec = as.double(count) / as.double(no_plots),
    .groups = "drop"
  ) |>
  filter(veg_type != "saltmarsh", no_plots != 1) |> 
  distinct(year, section, veg_type, taxon_code, no_plots, fraction_sec)

write_rds(species_sec_long, "data/species_sec_long.rds")

richness_sec_df <- species_sec_long |> 
  group_by(year, section, veg_type) |> 
  summarise(richness = n_distinct(taxon_code), .groups = "drop")

names(richness_sec_df)

#### mixed linear modelling species plot richness ####

m_richness_sec <- lmer(richness ~ year + (1|section), data = richness_sec_df)
summary(m_richness_sec)

#### visualisering plot species richness ####


# Get model predictions for richhness
pred_richness <- ggeffects::ggpredict(m_richness_sec, terms = "year")

# Plot observed and predicted evenness
ggplot(richness_sec_df, aes(x = year, y = richness)) +
  geom_jitter(aes(group = section), width = 0.2, alpha = 0.2, color = "#076834") +
  geom_boxplot(aes(group = factor(year)), outlier.shape = NA, alpha = 0.5, color = "gray30", width = 0.6) +
  geom_line(
    data = as.data.frame(pred_richness),
    aes(x = x, y = predicted),
    linetype = "dashed",
    color = "#076834",
    linewidth = 1.2,
    inherit.aes = FALSE
  ) +
  geom_ribbon(
    data = as.data.frame(pred_richness),
    aes(x = x, ymin = conf.low, ymax = conf.high),
    fill = "#076834",
    alpha = 0.2,
    inherit.aes = FALSE
  ) +
  scale_x_continuous(breaks = c(2007, 2012, 2017, 2022)) +
  labs(
    x = "Year",
    y = "Richness (species pr subsection)",
    title = "Change in richness over time",
    subtitle = paste0("Linear mixed model (p = ",
                      formatC(summary(m_richness_sec)$coefficients["year", "Pr(>|t|)"], digits = 4, format = "f"), ")")
  )

#### EVENNESS #########################
str(species_sec_long)

species_sec_long <- species_sec_long |> 
  filter(taxon_code != "rock") |> 
  mutate(abundance = fraction_sec) |> 
  distinct()

evenness_species <- species_sec_long |> 
  # compute diversity per plot-year
  group_by(year, section, veg_type) |>
  summarise(
    H = diversity(abundance, index = "shannon"),
    S = n(),                          # number of taxa
    J = ifelse(S > 1, H / log(S), NA_real_),
    .groups = "drop"
  )

head(evenness_species)

#### evenness model ############################################################


m_evenness_sec <- lmer(J ~ year + (1|section), data = evenness_species)
summary(m_evenness_sec)

#### evenness visualisation ############################################################
# Get model predictions (with CI) for evenness
pred_evenness <- ggeffects::ggpredict(m_evenness_sec, terms = "year")

# Plot observed and predicted evenness
ggplot(evenness_species, aes(x = year, y = J)) +
  geom_boxplot(aes(group = factor(year)), outlier.shape = NA, alpha = 0.5, color = "gray30", width = 0.6) +
  geom_line(
    data = as.data.frame(pred_evenness),
    aes(x = x, y = predicted),
    color = "#076834",
    linewidth = 1.2,
    inherit.aes = FALSE
  ) +
  geom_ribbon(
    data = as.data.frame(pred_evenness),
    aes(x = x, ymin = conf.low, ymax = conf.high),
    fill = "#076834",
    alpha = 0.2,
    inherit.aes = FALSE
  ) +
  geom_jitter(aes(group = section), width = 0.2, alpha = 0.2, color = "#076834") +
  scale_x_continuous(breaks = c(2007, 2012, 2017, 2022)) +
  labs(
    x = "Year",
    y = "Evenness (Pielou's J)",
    title = "Change in evenness over time for section",
    subtitle = paste0("Linear mixed model (p = ",
                      formatC(summary(m_evenness_sec)$coefficients["year", "Pr(>|t|)"], digits = 8, format = "f"), ")")
  )

#### SHANNON #########################

m_shannon <- lmer(H ~ year + (1 | section), data = evenness_species)
summary(m_shannon)

#### SHANNON visualisation #########################

# Get model predictions
pred_shannon <- ggpredict(m_shannon, terms = "year") |> as.data.frame()

# Extract p-value for year effect from model summary
pval_shannon <- summary(m_shannon)$coefficients["year", "Pr(>|t|)"]

# Get model predictions (with CI) for Shannon diversity
pred_shannon <- ggeffects::ggpredict(m_shannon, terms = "year")

# Plot observed and predicted Shannon diversity
ggplot(evenness_species, aes(x = year, y = H)) +
  geom_jitter(aes(group = section), width = 0.2, alpha = 0.2, color = "#076834") +
  geom_boxplot(aes(group = factor(year)), outlier.shape = NA, alpha = 0.5, color = "gray30", width = 0.6) +
  geom_line(
    data = as.data.frame(pred_shannon),
    aes(x = x, y = predicted),
    color = "#076834",
    linewidth = 1.2,
    inherit.aes = FALSE
  ) +
  geom_ribbon(
    data = as.data.frame(pred_shannon),
    aes(x = x, ymin = conf.low, ymax = conf.high),
    fill = "#076834",
    alpha = 0.2,
    inherit.aes = FALSE
  ) +
  scale_x_continuous(breaks = c(2007, 2012, 2017, 2022)) +
  labs(
    x = "Year",
    y = "Shannon diversity (H)",
    title = "Change in Shannon diversity over time",
    subtitle = paste0("Linear mixed model (p = ",
                      formatC(summary(m_shannon)$coefficients["year", "Pr(>|t|)"], digits = 4, format = "f"), ")")
  )

