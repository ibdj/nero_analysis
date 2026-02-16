
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

func_summed_allshrub <- species_sub |> 
  group_by(year,subsection, veg_type, fraction, func_type) |> 
  reframe(sum_frac = sum(fraction))

#### model ####

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
cld_labels <- cld_years_allshrub |>
  dplyr::select(year, .group)

# 2️⃣  Join the letters to the emmeans data frame
year_emm_lab <- year_emm_allshrub |>
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
                    ymin = asymp.LCL ,
                    ymax = asymp.UCL),
                width = 0.2,
                colour = "black",
                inherit.aes = FALSE) +   # <‑ prevent inheritance of y = abundance
  
  labs(x = "Year",
       y = "Frequency of occurrence") +
  # **Add the CLD letters**
  geom_text(data = year_emm_lab,
            aes(x = factor(year),
                y = 4,   # a little above the CI bar
                label = .group),
            vjust = 0,
            colour = "darkgreen",
            size = 5) +
  theme_minimal()

print(p)


