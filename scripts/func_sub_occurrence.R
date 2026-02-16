
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
library(cld)

#### importing data ####
merged_data <- readRDS("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/01 Vegetation changes Kobbefjord/data/nero_analysis/data/merged_data.rds") |> 
  group_by(year, subsection) |> 
  mutate(no_plots = n_distinct(plot_id)) |> 
  ungroup()

pivot <- merged_data |> 
  dplyr::select(year, vt_section,subsection, plot, plot_id, no_plots, veg_type, func_type, ecoveg_gfc, ecoveg_sgfc,presence, taxon_code) |> 
  distinct(year, vt_section,subsection, plot, plot_id, no_plots, veg_type, presence, taxon_code)

meta <- merged_data |> 
  dplyr::select(taxon_code, func_type, ecoveg_gfc, ecoveg_sgfc) |> 
  distinct()

pivot_wide <- pivot |> 
  pivot_wider(names_from = taxon_code, values_from = presence, values_fill = 0 )

species_long <- pivot_wide |> 
  pivot_longer(cols = 8:ncol(pivot_wide), names_to = "taxon_code", values_to = "occurrence")

species_sub <- species_long |> 
  group_by(year,subsection,no_plots,veg_type, taxon_code) |> 
  summarise(count = sum(occurrence)) |> 
  mutate(fraction = count/no_plots) |> 
  left_join(meta, by = "taxon_code")

names(species_sub)

#### model all shrub ####

func_summed_allshrub <- species_sub |> 
  group_by(year,subsection, veg_type, fraction, func_type) |> 
  reframe(sum_frac = sum(fraction)) |> 
  filter(func_type == "shrub")

summary(func_summed_allshrub)

model_sum_frac_allshrub <- lmer(
  sum_frac ~ factor(year) + (1 | subsection),
  data = func_summed_allshrub
)

summary(model_sum_frac_allshrub)

# Assuming your model object is called `mod_abund`
pairwise_sum_frac_allshrub <- model_sum_frac_allshrub |>
  emmeans(~ factor(year)) |>
  pairs(adjust = "holm")   # Holm correction for multiple testing for 6 tests

pairwise_sum_frac_allshrub

#### plotting summed fraction of all shrub #####

year_emm_allshrub <- model_sum_frac_allshrub |>
  emmeans(~ factor(year)) |>
  as.data.frame()      # turn the result into a plain data frame

head(year_emm_allshrub)

# cld panel 

# 1️⃣  Get the marginal means (already stored in `year_emm_sal`)
emm_years_allshrub <- model_sum_frac_allshrub |>
  emmeans(~ factor(year))

# 2️⃣  Compute the CLD letters (Holm‑adjusted pairwise tests)
cld_years_allshrub <- emm_years_allshrub |>
  cld(method = "holm", Letters = letters)   # you can change `Letters` if you prefer

cld_years_allshrub

# 1️⃣  Keep only the columns we need from the CLD output
cld_labels_allshrub <- cld_years_allshrub |>
  dplyr::select(year, .group)

# 2️⃣  Join the letters to the emmeans data frame
year_emm_lab_allshrub <- year_emm_allshrub |>
  dplyr::left_join(cld_labels, by = "year")

p <- ggplot(data = func_summed_allshrub,
            aes(x = factor(year), y = sum_frac)) +
  
  # raw data points (jittered)
  geom_jitter(width = 0.15,
              height = 0,
              alpha = 0.25,
              colour = "gray60") +
  
  # model‑based means (black dots)
  geom_point(data = year_emm_allshrub,
             aes(x = factor(year), y = emmean),
             colour = "black",
             size = 3) +
  
  # *** corrected error‑bar layer ***
  geom_errorbar(data = year_emm_allshrub,
                aes(x = factor(year),
                    ymin = lower.CL ,
                    ymax = upper.CL),
                width = 0.2,
                colour = "black",
                inherit.aes = FALSE) +   # <‑ prevent inheritance of y = abundance
  
  labs(x = "Year",
       y = "Frequency of occurrence") +
  # **Add the CLD letters**
  geom_text(data = year_emm_lab_allshrub,
            aes(x = factor(year),
                y = 1,   # a little above the CI bar
                label = .group),
            vjust = 0,
            colour = "darkgreen",
            size = 5) +
  theme_minimal()

print(p)


#### model deciduous ####

func_summed_deciduous <- species_sub |> 
  group_by(year,subsection, veg_type, fraction, ecoveg_sgfc) |> 
  reframe(sum_frac = sum(fraction)) |> 
  filter(ecoveg_sgfc == "shrub_decidous")

model_sum_frac_deciduous <- lmer(
  sum_frac ~ factor(year) + (1 | subsection),
  data = func_summed_deciduous
)

summary(model_sum_frac_deciduous)

# Assuming your model object is called `mod_abund`
pairwise_sum_frac_deciduous <- model_sum_frac_deciduous |>
  emmeans(~ factor(year)) |>
  pairs(adjust = "holm")   # Holm correction for multiple testing for 6 tests

pairwise_sum_frac_deciduous

#### plotting summed fraction deciduous #####

year_emm_deciduous <- model_sum_frac_deciduous |>
  emmeans(~ factor(year)) |>
  as.data.frame()      # turn the result into a plain data frame

head(year_emm_deciduous)

# cld panel 

# 1️⃣  Get the marginal means (already stored in `year_emm_sal`)
emm_years_deciduous <- model_sum_frac_deciduous |>
  emmeans(~ factor(year))

# 2️⃣  Compute the CLD letters (Holm‑adjusted pairwise tests)
cld_years_deciduous <- emm_years_deciduous |>
  cld(method = "holm", Letters = letters)   # you can change `Letters` if you prefer

cld_years_deciduous

# 1️⃣  Keep only the columns we need from the CLD output
cld_labels_deciduous <- cld_years_deciduous |>
  dplyr::select(year, .group)

# 2️⃣  Join the letters to the emmeans data frame
year_emm_lab_deciduous <- year_emm_deciduous |>
  dplyr::left_join(cld_labels, by = "year")

p <- ggplot(data = func_summed_deciduous,
            aes(x = factor(year), y = sum_frac)) +
  
  # raw data points (jittered)
  geom_jitter(width = 0.15,
              height = 0,
              alpha = 0.25,
              colour = "gray60") +
  
  # model‑based means (black dots)
  geom_point(data = year_emm_deciduous,
             aes(x = factor(year), y = emmean),
             colour = "darkgreen",
             size = 3) +
  
  # *** corrected error‑bar layer ***
  geom_errorbar(data = year_emm_deciduous,
                aes(x = factor(year),
                    ymin = lower.CL ,
                    ymax = upper.CL),
                width = 0.2,
                colour = "darkgreen",
                inherit.aes = FALSE) +   # <‑ prevent inheritance of y = abundance
  
  labs(x = "Year",
       y = "Frequency of occurrence") +
  # **Add the CLD letters**
  geom_text(data = year_emm_lab_deciduous,
            aes(x = factor(year),
                y = 0.75,   # a little above the CI bar
                label = .group),
            vjust = 0,
            colour = "darkgreen",
            size = 5) +
  theme_minimal()

print(p)

#### model evergreen ####

func_summed_evergreen <- species_sub |> 
  group_by(year,subsection, veg_type, fraction, ecoveg_sgfc) |> 
  reframe(sum_frac = sum(fraction)) |> 
  filter(ecoveg_sgfc == "shrub_evergreen")

model_sum_frac_evergreen <- lmer(
  sum_frac ~ factor(year) + (1 | subsection),
  data = func_summed_evergreen
)

summary(model_sum_frac_evergreen)

# Assuming your model object is called `mod_abund`
pairwise_sum_frac_evergreen <- model_sum_frac_evergreen |>
  emmeans(~ factor(year)) |>
  pairs(adjust = "holm")   # Holm correction for multiple testing for 6 tests

pairwise_sum_frac_evergreen

#### plotting summed fraction evergreen #####

year_emm_evergreen <- model_sum_frac_evergreen |>
  emmeans(~ factor(year)) |>
  as.data.frame()      # turn the result into a plain data frame

head(year_emm_evergreen)

# cld panel 

# 1️⃣  Get the marginal means (already stored in `year_emm_sal`)
emm_years_evergreen <- model_sum_frac_evergreen |>
  emmeans(~ factor(year))

# 2️⃣  Compute the CLD letters (Holm‑adjusted pairwise tests)
cld_years_evergreen <- emm_years_evergreen |>
  cld(method = "holm", Letters = letters)   # you can change `Letters` if you prefer

cld_years_evergreen

# 1️⃣  Keep only the columns we need from the CLD output
cld_labels_evergreen <- cld_years_evergreen |>
  dplyr::select(year, .group)

# 2️⃣  Join the letters to the emmeans data frame
year_emm_lab_evergreen <- year_emm_evergreen |>
  dplyr::left_join(cld_labels, by = "year")

p <- ggplot(data = func_summed_evergreen,
            aes(x = factor(year), y = sum_frac)) +
  
  # raw data points (jittered)
  geom_jitter(width = 0.15,
              height = 0,
              alpha = 0.1,
              colour = "darkgray") +
  
  # model‑based means (black dots)
  geom_point(data = year_emm_evergreen,
             aes(x = factor(year), y = emmean),
             colour = "darkred",
             size = 3) +
  
  # *** corrected error‑bar layer ***
  geom_errorbar(data = year_emm_evergreen,
                aes(x = factor(year),
                    ymin = lower.CL ,
                    ymax = upper.CL),
                width = 0.2,
                colour = "darkred",
                inherit.aes = FALSE) +   # <‑ prevent inheritance of y = abundance
  
  labs(x = "Year",
       y = "Frequency of occurrence") +
  # **Add the CLD letters**
  geom_text(data = year_emm_lab_evergreen,
            aes(x = factor(year),
                y = 0.75,   # a little above the CI bar
                label = .group),
            vjust = 0,
            colour = "darkred",
            size = 5) +
  theme_minimal()

print(p)

#### combined plit#####

# Combine the raw‑point data
raw_long <- bind_rows(
  evergreen = func_summed_evergreen,
  deciduous = func_summed_deciduous,
  allshrub  = func_summed_allshrub,
  .id = "func_group"          # adds a column called func_group
) %>% 
  mutate(func_group = factor(func_group,
                             levels = c("evergreen", "deciduous", "allshrub")))

# Combine the estimated means and confidence limits
emm_long <- bind_rows(
  evergreen = year_emm_evergreen,
  deciduous = year_emm_deciduous,
  allshrub  = year_emm_allshrub,
  .id = "func_group"
) %>%
  mutate(
    func_group = factor(func_group,
                        levels = c("evergreen", "deciduous", "allshrub"))
  )
emm_long

# combine the three CLD tables into one long data frame
cld_long <- bind_rows(
  evergreen = year_emm_lab_evergreen,
  deciduous = year_emm_lab_deciduous,
  allshrub  = year_emm_lab_allshrub,
  .id = "func_group"
) |>
  # make sure the group identifier is a factor with the desired order
  mutate(
    func_group = factor(
      func_group,
      levels = c("evergreen", "deciduous", "allshrub")
    )
  )

p <-  ggplot(data = raw_long, aes(x = factor(year), y = sum_frac, colour = func_group)) +
  geom_jitter(
    width = 0.15,
    height = 0,
    alpha = 0.15,
    size = 1.2
  ) +
  scale_colour_manual(
    name = "Functional group",
    values = c(
      evergreen = "darkred",
      deciduous = "darkgreen",
      allshrub  = "blue"
    )
  )

str(raw_long)
p

p <- p +
  # ----- model‑based means (larger points) -----
geom_point(
  data = emm_long,
  aes(x = factor(year), y = emmean, colour = func_group),
  size = 3,
  shape = 19,
  inherit.aes = FALSE          # <-- do NOT pull in y = sum_frac
) +
  # ----- confidence‑interval error bars -----
geom_errorbar(
  data = emm_long,
  aes(
    x = factor(year),
    ymin = lower.CL,
    ymax = upper.CL,
    colour = func_group
  ),
  width = 0.2,
  size = 0.8,
  inherit.aes = FALSE          # <-- prevents the stray sum_frac mapping
)
p

p <- p +
  geom_text(
    data = cld_long,
    aes(
      x = factor(year),          # same x‑position as the other layers
      y = 0.85,                  # place the label just above the error bar;
      # adjust this value later if needed
      label = .group,            # the CLD letter (e.g., "a", "b")
      colour = func_group        # keep the same colour coding
    ),
    vjust = 0,                   # anchor the text at the bottom of the label
    size = 5,
    inherit.aes = FALSE          # crucial – don’t inherit y = sum_frac
  )
p
