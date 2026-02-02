#### loading packages ####

library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggplot2)
library(ggeffects)
library(vegan)
library(codyn) # turnover calculations
library(multcomp)
library(multcompView)
library(patchwork)

#to do
# finished indicies plot
# finished report text for all nmds
# get overview of all the plots in the paper

#### loading data ################################

merged_data <- readRDS("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/01 Vegetation changes Kobbefjord/data/nero_analysis/data/merged_data.rds") |> 
  group_by(year, subsection) |> 
  mutate(no_plots = n_distinct(plot_id)) |> 
  ungroup()

species_sub_long <- merged_data |>
  group_by(year, subsection, veg_type, taxon_code, no_plots) |>
  summarize(
    count = n(),
    fraction_sub = as.double(count) / as.double(no_plots),
    .groups = "drop"
  ) |>
  filter(veg_type != "saltmarsh", no_plots != 1) |> 
  distinct(year, subsection, veg_type, taxon_code, no_plots, fraction_sub)

write_rds(species_sub_long, "data/species_sub_long.rds")

richness_sub_df <- species_sub_long |> 
  group_by(year, subsection, veg_type) |> 
  summarise(richness = n_distinct(taxon_code), .groups = "drop")

names(richness_sub_df)

species_sub_long <- species_sub_long |> 
  filter(taxon_code != "rock") |> 
  mutate(abundance = fraction_sub) |> 
  distinct()

evenness_species <- species_sub_long |> 

# compute diversity per plot-year
  group_by(year, subsection, veg_type) |>
  summarise(
    H = diversity(abundance, index = "shannon"),
    S = n(),                          # number of taxa
    J = ifelse(S > 1, H / log(S), NA_real_),
    .groups = "drop"
  )

head(evenness_species)

#### richness mlm ####

m_richness_sub <- lmer(richness ~ year + (1|subsection), data = richness_sub_df)
summary(m_richness_sub)

# year as factor and pairwise comparison
model_richness_sub_f <- lmer(
  richness ~ factor(year) + (1 | subsection),
  data = richness_sub_df
)

summary(model_richness_sub_f)

richness_emmeans_f <- emmeans(model_richness_sub_f, ~ year) 

richness_pairs_f <- emmeans(model_richness_sub_f, ~ year) |> 
  pairs()

richness_emmeans_f
richness_pairs_f

richness_pairs_df <- as.data.frame(richness_pairs_f)
richness_pairs_df
#### richness factor visualisation #####

richness_sub_df$year <- factor(richness_sub_df$year)

richness_emmeans_f_df <- as.data.frame(richness_emmeans_f)

# Get emmeans as data frame
richness_emm_df <- as.data.frame(richness_emmeans_f)

richness_emm_df <- richness_emm_df[order(richness_emm_df$year), ]

years <- c("2007", "2012", "2017", "2022")

richness_sub_df$year  <- factor(richness_sub_df$year,  levels = years)
richness_emm_df$year  <- factor(richness_emm_df$year,  levels = years)

richness_cld <- cld(
  richness_emmeans_f,
  adjust = "sidak",
  Letters = letters
)
richness_cld$year <- factor(richness_cld$year, levels = years)

plot_richness <- ggplot() +
  # Raw data
  geom_jitter(
    data = richness_sub_df,
    aes(x = year, y = richness),
    width = 0.15,
    alpha = 0.3,
    color = "black"
  ) +
  # Model-based confidence intervals
  geom_errorbar(
    data = richness_emm_df,
    aes(x = year, ymin = lower.CL, ymax = upper.CL),
    width = 0.25,
    linewidth = 0.9,
    color = "darkgreen"
  ) +
  # Model-based means
  geom_point(
    data = richness_emm_df,
    aes(x = year, y = emmean),
    size = 2.8,
    color = "darkgreen"
  ) +
  labs(
    x = "Year",
    y = "Species richness"
  ) +
  #scale_x_discontinuous(breaks = c(2007,2012,2017,2022))+
  theme_classic()  +
    geom_text(
    data = richness_cld,
    aes(x = year, label = .group),
    y = 25,  # Tune this
    size = 4, color = "darkblue"
  )

plot_richness

#### richness visualisation ####

# Get model predictions for richhness

years <- c("2007", "2012", "2017", "2022")

richness_sub_df$year  <- factor(richness_sub_df$year,  levels = years)
richness_emm_df$year  <- factor(richness_emm_df$year,  levels = years)
pred_richness <- ggeffects::ggpredict(m_richness_sub, terms = "year")

# Plot observed and predicted evenness

ggplot(richness_sub_df, aes(x = year, y = richness)) +
  geom_jitter(aes(group = subsection), width = 0.2, alpha = 0.2, color = "#076834") +
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
  #scale_x_continuous(breaks = c(2007, 2012, 2017, 2022)) +
  labs(
    x = "Year",
    y = "Richness (species pr subsection)",
    title = "Change in richness over time",
    subtitle = paste0("Linear mixed model (p = ",
                      formatC(summary(m_richness_sub)$coefficients["year", "Pr(>|t|)"], digits = 4, format = "f"), ")")
  )+
  theme_minimal()

#### EVENNESS #########################
str(species_sub_long)

species_sub_long <- species_sub_long |> 
  filter(taxon_code != "rock") |> 
  mutate(abundance = fraction_sub) |> 
  distinct()
  
evenness_species <- species_sub_long |> 
  # compute diversity per plot-year
  group_by(year, subsection, veg_type) |>
  summarise(
    H = diversity(abundance, index = "shannon"),
    S = n(),                          # number of taxa
    J = ifelse(S > 1, H / log(S), NA_real_),
    .groups = "drop"
  )

head(evenness_species)

#### evenness mlm ############################################################

m_evenness_sub <- lmer(J ~ year + (1|subsection), data = evenness_species)
summary(m_evenness_sub)

# year as factor and pairwise comparison
m_evenness_sub_f <- lmer(
  J ~ factor(year) + (1 | subsection),
  data = evenness_species
)

summary(m_evenness_sub_f)

evenness_emmeans <- emmeans(m_evenness_sub_f, ~ year) |>
  pairs()

evenness_emmeans
#### evenness factor visualisation #############################################

evenness_emm <- emmeans(m_evenness_sub_f, ~ year)
evenness_emm_df <- as.data.frame(evenness_emm)

years <- c("2007", "2012", "2017", "2022")

evenness_species$year <- factor(evenness_species$year, levels = years)
evenness_emm_df$year  <- factor(evenness_emm_df$year,  levels = years)

evenness_emm_df <- evenness_emm_df[order(evenness_emm_df$year), ]

evenness_cld <- cld(
  evenness_emm,
  adjust = "tukey",
  Letters = letters
)

evenness_cld$year <- factor(evenness_cld$year, 
                            levels = levels(evenness_emm_df$year))

evenness_cld$year
evenness_emm_df$year

plot_evenness <- ggplot() +
  geom_jitter(data = evenness_species, aes(x = year, y = J), 
              width = 0.15, alpha = 0.3, color = "black") +
  geom_errorbar(data = evenness_emm_df, 
                aes(x = year, ymin = lower.CL, ymax = upper.CL),
                width = 0.25, linewidth = 0.9, color = "darkblue") +
  geom_point(data = evenness_emm_df, 
             aes(x = year, y = emmean), 
             size = 2.8, color = "darkblue") +
  # FIXED: aes(x = year, ...) INSIDE geom_text
  geom_text(
    data = evenness_cld,
    aes(x = year, label = .group),
    y = 1,  # Tune this
    size = 4, color = "darkblue"
  ) +
  labs(x = "Year", y = "Pielou's evenness (J)") +
  theme_classic()

plot_evenness


#### evenness visualisation ############################################################

# Get model predictions (with CI) for evenness
pred_evenness <- ggeffects::ggpredict(m_evenness_sub, terms = "year")

# Plot observed and predicted evenness
ggplot(evenness_species, aes(x = year, y = J)) +
  geom_boxplot(aes(group = factor(year)), outlier.shape = NA, alpha = 0.5, color = "gray30", width = 0.6) +
  geom_line(
    data = as.data.frame(pred_evenness),
    aes(x = x, y = predicted),
    color = "#017fad",
    linewidth = 1.2,
    inherit.aes = FALSE
  ) +
  geom_ribbon(
    data = as.data.frame(pred_evenness),
    aes(x = x, ymin = conf.low, ymax = conf.high),
    fill = "#017fad",
    alpha = 0.2,
    inherit.aes = FALSE
  ) +
  geom_jitter(aes(group = subsection), width = 0.2, alpha = 0.2, color = "#017fad") +
  scale_x_continuous(breaks = c(2007, 2012, 2017, 2022)) +
  labs(
    x = "Year",
    y = "Evenness (Pielou's J)",
    title = "Change in evenness over time",
    subtitle = paste0("Linear mixed model (p = ",
                      formatC(summary(m_evenness_sub)$coefficients["year", "Pr(>|t|)"], format = "e", digits = 2), ")")
  )


#### SHANNON #########################

shannon_df <- evenness_species

m_shannon_sub <- lmer(H ~ year + (1 | subsection), data = evenness_species)
summary(m_shannon_sub)

# year as factor and pairwise comparison
model_shannon_sub_f <- lmer(
  H ~ factor(year) + (1 | subsection),
  data = evenness_species
)

summary(model_shannon_sub_f)

shannon_emmeans <- emmeans(model_shannon_sub_f, ~ year) |>
  pairs()

shannon_emmeans

#### SHANNON FACTOR visualisation #########################

# 1. Get emmeans (you have pairs, but need emmeans for plotting)
shannon_emmeans_only <- emmeans(model_shannon_sub_f, ~ year)
shannon_emm_df <- as.data.frame(shannon_emmeans_only)

# Make year a factor in your data if not already
evenness_species$year <- factor(evenness_species$year)

year_levels <- shannon_emm_df$year

# STEP 1: Marginal means (4 rows)
shannon_emmeans_marginals <- emmeans(model_shannon_sub_f, ~ year)

# STEP 2: CLD on marginal means (4 rows, with letters)
shannon_cld <- cld(shannon_emmeans_marginals, adjust = "sidak", Letters = letters)

# STEP 3: Convert marginals to df with x_pos (4 rows)
shannon_emm_df <- as.data.frame(shannon_emmeans_marginals) |>
  mutate(x_pos = 1:4) |>
  arrange(year)

shannon_cld$x_pos <- shannon_emm_df$x_pos


plot_shannon <- ggplot() +
  geom_jitter(data = evenness_species, aes(x = x_pos, y = H), 
              width = 0.15, alpha = 0.3, color = "black") +
  geom_errorbar(data = shannon_emm_df, 
                aes(x = x_pos, ymin = lower.CL, ymax = upper.CL),
                width = 0.25, linewidth = 0.9, color = "darkgreen") +
  geom_point(data = shannon_emm_df, 
             aes(x = x_pos, y = emmean), 
             size = 2.8, color = "darkgreen") +
  geom_text(data = shannon_cld,
            aes(x = x_pos, label = .group),
            y = 3,  # Above max upper.CL
            size = 4, color = "darkgreen") +
  scale_x_continuous(name = "Year", breaks = 1:4, labels = c("2007", "2012", "2017", "2022")) +
  labs(y = "Shannon diversity (H)") +
  theme_classic()

plot_shannon
#### SHANNON visualisation #########################

# Get model predictions
pred_shannon <- ggpredict(m_shannon, terms = "year") |> as.data.frame()

# Extract p-value for year effect from model summary
pval_shannon <- summary(m_shannon)$coefficients["year", "Pr(>|t|)"]

# Get model predictions (with CI) for Shannon diversity
pred_shannon <- ggeffects::ggpredict(m_shannon, terms = "year")

# Plot observed and predicted Shannon diversity
ggplot(evenness_species, aes(x = year, y = H)) +
  geom_jitter(aes(group = subsection), width = 0.2, alpha = 0.2, color = "#017fad") +
  geom_boxplot(aes(group = factor(year)), outlier.shape = NA, alpha = 0.5, color = "gray30", width = 0.6) +
  geom_line(
    data = as.data.frame(pred_shannon),
    aes(x = x, y = predicted),
    color = "#017fad",
    linewidth = 1.2,
    inherit.aes = FALSE
  ) +
  geom_ribbon(
    data = as.data.frame(pred_shannon),
    aes(x = x, ymin = conf.low, ymax = conf.high),
    fill = "#017fad",
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


#### TURNOVER #########################

# Check number of unique years per plot
plot_years <- species_sub_long |>
  group_by(subsection) |>
  summarise(n_years = n_distinct(year), .groups = "drop")

# Filter plots with at least 2 years (turnover requires temporal comparison)
valid_plots <- plot_years |>
  filter(n_years > 1) |>
  pull(subsection)

# Filter data to valid plots only
species_long_filtered <- species_sub_long |>
  filter(subsection %in% valid_plots)

# Calculate turnover per plot between years
turnover_results <- turnover(
  df = species_long_filtered,
  time.var = "year",
  species.var = "taxon_code",
  abundance.var = "abundance",
  replicate.var = "subsection",
  metric = "total"
)

head(turnover_results)

#### turnover mlm ##############################################################

m_turnover_sub <- lmer(total ~ year + (1 | subsection), data = turnover_results)

# year as factor and pairwise comparison
m_turnover_sub_f <- lmer(
  total ~ factor(year) + (1 | subsection),
  data = turnover_results
)

summary(m_turnover_sub_f)

turnover_emmeans <- emmeans(m_turnover_sub_f, ~ year) |>
  pairs()

turnover_emmeans

#### TURNOVER visualisation ####################################################
# 2. Fit linear mixed model: turnover by year with random intercept for plot


# 3. Generate predicted values at observed years
pred_turnover <- ggpredict(m_turnover, terms = c("year"))

# 4. Prepare predicted data for plotting
pred_turnover_plot <- pred_turnover |>
  mutate(year = as.numeric(x))  # convert x (character) to numeric for plotting

# 5. Boxplot + jitter for observed turnover, line + ribbon for predicted
ggplot(turnover_results, aes(x = factor(year), y = total)) + 
  geom_boxplot(aes(group = factor(year)), outlier.shape = NA, alpha = 0.5, color = "gray30", width = 0.4) +
  geom_jitter(aes(group = subsection), width = 0.15, alpha = 0.2, color = "#017fad") +
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
    title = "Species turnover across years (observed + predicted)",
    subtitle = paste0(
      "Linear mixed model (p = ",
      formatC(summary(m_turnover)$coefficients["year", "Pr(>|t|)"], format = "e", digits = 2),
      ")"
  )) 

plot(resid(m_turnover))   
qqnorm(resid(m_turnover)); qqline(resid(m_turnover))    # Normality
qqnorm(ranef(m_turnover)$group$`(Intercept)`); qqline(ranef(m_turnover)$group$`(Intercept)`)

#### combined plots ####

# 1. Prepare predicted data frames with numeric year and original scale values

# Richness predictions
pred_richness <- ggpredict(m_richness_sub, terms = "year") |>
  as.data.frame() |>
  mutate(year = as.numeric(x),
         metric = "Richness",
         value = predicted,
         conf.low = conf.low,
         conf.high = conf.high) |>
  select(year, value, conf.low, conf.high, metric)

# Shannon predictions
pred_shannon <- ggpredict(m_shannon_sub, terms = c("year [2007,2012,2017,2022]")) |>
  as.data.frame() |>
  mutate(year = as.numeric(x),
         metric = "Shannon",
         value = predicted,
         conf.low = conf.low,
         conf.high = conf.high) |>
  select(year, value, conf.low, conf.high, metric)

# Turnover predictions
pred_turnover <- ggpredict(m_turnover_sub, terms = c("year")) |>
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
rich_scale <- scale_with_params(richness_sub_df$richness)
richness_obs_scaled <- richness_sub_df |>
  mutate(
    value  = rich_scale$scaled,
    metric = "Richness"
  ) |>
  select(year, subsection, value, metric) |> 
  mutate(subsection = as.double(subsection))

## Shannon
shan_scale <- scale_with_params(shannon_df$H)
shannon_obs_scaled <- shannon_df |>
  mutate(
    value  = shan_scale$scaled,
    metric = "Shannon"
  ) |>
  select(year, subsection, value, metric) |> 
  mutate(subsection = as.double(subsection))

## Turnover
turn_scale <- scale_with_params(turnover_results$total)
turnover_obs_scaled <- turnover_results |>
  mutate(
    value  = turn_scale$scaled,
    metric = "Turnover"
  ) |>
  select(year, subsection, value, metric) |> 
  mutate(subsection = as.double(subsection))

observed_all <- bind_rows(richness_obs_scaled,
                          shannon_obs_scaled,
                          turnover_obs_scaled)

### 2. SCALE PREDICTIONS USING THE SAME MIN/MAX PER METRIC

## Richness preds
pred_richness <- ggpredict(m_richness_sub, terms = "year") |>
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
    linewidth = 1.1, alpha = 0.3
  ) +
  scale_x_discrete(breaks = c("2007","2012","2017","2022")) +
  labs(
    x     = "Year",
    y     = "Scaled value (0–1, within metric)",
    title = "",
    #titel: Relative change in richness, Shannon diversity and turnover over time (pr plot)
    fill  = "Metric",
    color = "Metric"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "inside", legend.position.inside = c(0.2, 1), legend.direction = "horizontal", legend.title = element_blank())


#### combined factor plots #######

# 2x2 grid (top row: Shannon + Richness, bottom row: Evenness + empty)
# WRAP everything in plot_layout()
( plot_shannon  | plot_richness ) /
  ( plot_evenness | plot_spacer() ) +
  
  plot_layout(
    guides  = "collect",
    heights = c(1, 1),
    widths  = c(1, 1)
  ) &
  
  theme(legend.position = "bottom")
  