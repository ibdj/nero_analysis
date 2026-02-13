
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
  mutate(fraction = count/no_plots)
  

species_sub_salgla <- species_sub |> 
  filter(taxon_code == "salgla")


#### salix model frequency of occurrence ####


# Fit the model
mod_abund_sal <- lmer(
  fraction ~ factor(year) + (1 | subsection),
  data = species_sub_salgla,          
  REML = FALSE
)

# Show the fitted object (so you can verify it was created)
mod_abund_sal

# Assuming your model object is called `mod_abund`
pairwise_years <- mod_abund_sal |>
  emmeans(~ factor(year)) |>
  pairs(adjust = "holm")   # Holm correction for multiple testing

pairwise_years

#### salix plot frequency of occurrence ####
year_emm_sal <- mod_abund_sal |>
  emmeans(~ factor(year)) |>
  as.data.frame()      # turn the result into a plain data frame

head(year_emm_sal)

# cld panel 

# 1️⃣  Get the marginal means (already stored in `year_emm_sal`)
emm_years <- mod_abund_sal |>
  emmeans(~ factor(year))

# 2️⃣  Compute the CLD letters (Holm‑adjusted pairwise tests)
cld_years <- emm_years |>
  cld(method = "holm", Letters = letters)   # you can change `Letters` if you prefer

cld_years

# 1️⃣  Keep only the columns we need from the CLD output
cld_labels <- cld_years |>
  dplyr::select(year, .group)

# 2️⃣  Join the letters to the emmeans data frame
year_emm_lab <- year_emm_sal |>
  dplyr::left_join(cld_labels, by = "year")

p <- ggplot(data = species_sub_salgla,
            aes(x = factor(year), y = fraction)) +
  
  # raw data points (jittered)
  geom_jitter(width = 0.15,
              height = 0,
              alpha = 0.25,
              colour = "gray60") +
  
  # model‑based means (black dots)
  geom_point(data = year_emm_sal,
             aes(x = factor(year), y = emmean),
             colour = "black",
             size = 3) +
  
  # *** corrected error‑bar layer ***
  geom_errorbar(data = year_emm_sal,
                aes(x = factor(year),
                    ymin = lower.CL,
                    ymax = upper.CL),
                width = 0.2,
                colour = "black",
                inherit.aes = FALSE) +   # <‑ prevent inheritance of y = abundance
  
  labs(x = "Year",
       y = "Frequency of occurrence") +
  # **Add the CLD letters**
  geom_text(data = year_emm_lab,
            aes(x = factor(year),
                y = 0.75,   # a little above the CI bar
                label = .group),
            vjust = 0,
            colour = "red",
            size = 5) +
  theme_minimal()

print(p)



#### betula model ####

species_sub_betnan <- species_sub |> 
  filter(taxon_code == "betnan")
# Fit the model
mod_abund_betnan <- lmer(
  fraction ~ factor(year) + (1 | subsection),
  data = species_sub_betnan,          
  REML = FALSE
)

# Show the fitted object (so you can verify it was created)
mod_abund_betnan

# Assuming your model object is called `mod_abund`
pairwise_years_betnan <- mod_abund_betnan |>
  emmeans(~ factor(year)) |>
  pairs(adjust = "holm")   # Holm correction for multiple testing

pairwise_years_betnan

#### betula plot freq ####
year_emm_bet <- mod_abund_betnan |>
  emmeans(~ factor(year)) |>
  as.data.frame()      # turn the result into a plain data frame

head(year_emm_bet)

# cld panel 

# 1️⃣  Get the marginal means (already stored in `year_emm_sal`)
emm_years_bet <- mod_abund_betnan |>
  emmeans(~ factor(year))

# 2️⃣  Compute the CLD letters (Holm‑adjusted pairwise tests)
cld_years_bet <- emm_years_bet |>
  cld(method = "holm", Letters = letters)   # you can change `Letters` if you prefer

cld_years_bet

# 1️⃣  Keep only the columns we need from the CLD output
cld_labels_bet <- cld_years_bet |>
  dplyr::select(year, .group)

# 2️⃣  Join the letters to the emmeans data frame
year_emm_lab_bet <- year_emm_bet |>
  dplyr::left_join(cld_labels_bet, by = "year")

p <- ggplot(data = species_sub_betnan,
            aes(x = factor(year), y = fraction)) +
  
  # raw data points (jittered)
  geom_jitter(width = 0.15,
              height = 0,
              alpha = 0.25,
              colour = "gray60") +
  
  # model‑based means (black dots)
  geom_point(data = year_emm_bet,
             aes(x = factor(year), y = emmean),
             colour = "black",
             size = 3) +
  
  # *** corrected error‑bar layer ***
  geom_errorbar(data = year_emm_bet,
                aes(x = factor(year),
                    ymin = lower.CL,
                    ymax = upper.CL),
                width = 0.2,
                colour = "black",
                inherit.aes = FALSE) +   # <‑ prevent inheritance of y = abundance
  
  labs(x = "Year",
       y = "Frequency of occurrence") +
  # **Add the CLD letters**
  geom_text(data = year_emm_lab_bet,
            aes(x = factor(year),
                y = 0.75,   # a little above the CI bar
                label = .group),
            vjust = 0,
            colour = "red",
            size = 5) +
  theme_minimal()

print(p)




#### empetrum model ####

species_sub_empnig <- species_sub |> 
  filter(taxon_code == "empnig")
# Fit the model
mod_abund_empnig <- lmer(
  fraction ~ factor(year) + (1 | subsection),
  data = species_sub_empnig,          
  REML = FALSE
)

# Show the fitted object (so you can verify it was created)
mod_abund_empnig

# Assuming your model object is called `mod_abund`
pairwise_years_empnig <- mod_abund_empnig |>
  emmeans(~ factor(year)) |>
  pairs(adjust = "holm")   # Holm correction for multiple testing

pairwise_years_empnig

#### betula plot freq ####
year_emm_empnig <- mod_abund_empnig |>
  emmeans(~ factor(year)) |>
  as.data.frame()      # turn the result into a plain data frame

head(year_emm_empnig)

# cld panel 

# 1️⃣  Get the marginal means (already stored in `year_emm_sal`)
emm_years_empnig <- mod_abund_empnig |>
  emmeans(~ factor(year))

# 2️⃣  Compute the CLD letters (Holm‑adjusted pairwise tests)
cld_years_empnig <- emm_years_empnig |>
  cld(method = "holm", Letters = letters)   # you can change `Letters` if you prefer

cld_years_empnig

# 1️⃣  Keep only the columns we need from the CLD output
cld_labels_empnig <- cld_years_empnig |>
  dplyr::select(year, .group)

# 2️⃣  Join the letters to the emmeans data frame
year_emm_lab_empnig <- year_emm_empnig |>
  dplyr::left_join(cld_labels_empnig, by = "year")

p <- ggplot(data = species_sub_empnig,
            aes(x = factor(year), y = fraction)) +
  
  # raw data points (jittered)
  geom_jitter(width = 0.15,
              height = 0,
              alpha = 0.25,
              colour = "gray60") +
  
  # model‑based means (black dots)
  geom_point(data = year_emm_empnig,
             aes(x = factor(year), y = emmean),
             colour = "black",
             size = 3) +
  
  # *** corrected error‑bar layer ***
  geom_errorbar(data = year_emm_empnig,
                aes(x = factor(year),
                    ymin = lower.CL,
                    ymax = upper.CL),
                width = 0.2,
                colour = "black",
                inherit.aes = FALSE) +   # <‑ prevent inheritance of y = abundance
  
  labs(x = "Year",
       y = "Frequency of occurrence") +
  # **Add the CLD letters**
  geom_text(data = year_emm_lab_empnig,
            aes(x = factor(year),
                y = 0.75,   # a little above the CI bar
                label = .group),
            vjust = 0,
            colour = "red",
            size = 5) +
  theme_minimal()

print(p)

