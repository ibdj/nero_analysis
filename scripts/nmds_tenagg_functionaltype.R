#### nmds with 10 agreegate and functional type ####
####################################### packages ############################################
# Install if not already installed
if (!requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("pmartinezarbizu/pairwiseAdonis")

library(tidyverse)
library(viridis)
library(vegan)
library(remotes)
library(pairwiseAdonis)


####################################### loading data ############################################

merged_data <- readRDS("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/01 Vegetation changes Kobbefjord/data/nero_analysis/data/merged_data.rds") |> 
  select(year,vt_section,plot_id,veg_type,taxon_code,presence,species,ecoveg_gfc) |> 
  filter(veg_type != "saltmarsh")

summary(merged_data)

names(merged_data)
#"year"       "vt_section" "plot_id"    "veg_type"   "taxon_code" "presence"   "species"    "ecoveg_gfc"

tenaggr_functype <- merged_data |> 
  group_by(year, vt_section, veg_type, ecoveg_gfc) |> 
  summarise(count = n())

tenaggr_functype_wide <- tenaggr_functype |> 
  pivot_wider(names_from = "ecoveg_gfc", values_from = "count", values_fill = 0)

head(tenaggr_functype_wide)

####################################### nmds ############################################

nmds_all <- tenaggr_functype_wide |>
  ungroup() |>
  select(-year, -vt_section, -veg_type, -`NA`) |>
  mutate(across(everything(), as.numeric)) |>
  metaMDS(distance = "bray", k = 2, trymax = 100)

nmds_all

# Call:
#   metaMDS(comm = mutate(select(ungroup(tenaggr_functype_wide),      -year, -vt_section, -veg_type, -`NA`), across(everything(),      as.numeric)), distance = "bray", k = 2, trymax = 100) 
# 
# global Multidimensional Scaling using monoMDS
# 
# Data:     wisconsin(mutate(select(ungroup(tenaggr_functype_wide), -year, -vt_section, -veg_type, -`NA`), across(everything(), as.numeric))) 
# Distance: bray 
# 
# Dimensions: 2 
# Stress:     0.1427319 
# Stress type 1, weak ties
# Best solution was not repeated after 100 tries
# The best solution was from try 0 (metric scaling or null solution)
# Scaling: centring, PC rotation, halfchange scaling 
# Species: expanded scores based on ‘wisconsin(mutate(select(ungroup(tenaggr_functype_wide), -year, -vt_section, -veg_type, -`NA`), across(everything(), as.numeric)))’ 

####################################### nmds meta data ############################################
nmds_points <- as.data.frame(nmds_all$points) |>
  bind_cols(tenaggr_functype_wide |> select(year, vt_section, veg_type))


####################################### convex hulls ############################################
hulls <- nmds_points |>
  group_by(year) |>
  slice(chull(MDS1, MDS2))

####################################### centroids ############################################
centroids <- nmds_points |>
  group_by(year) |>
  summarise(
    MDS1 = mean(MDS1),
    MDS2 = mean(MDS2)
  )
####################################### calculating loadings ############################################
species_scores <- as.data.frame(scores(nmds_all, "species"))
species_scores$variable <- rownames(species_scores)

species_scores <- species_scores |>
  mutate(
    MDS1 = NMDS1 * 1.2,
    MDS2 = NMDS2 * 1.2
  )

####################################### plotting ############################################
ggplot(nmds_points, aes(x = MDS1, y = MDS2, colour = factor(year))) +
  geom_point(size = 2, alpha = 0.7) +
  geom_polygon(
    data = hulls,
    aes(fill = factor(year), colour = NULL),
    alpha = 0.2
  ) +
  geom_point(
    data = centroids,
    aes(x = MDS1, y = MDS2, fill = factor(year)),
    shape = 21,
    colour = "black",
    size = 4
  ) +
  theme_minimal() +
  labs(
    title = "NMDS of Functional Types by Year with Centroids",
    colour = "Year",
    fill = "Year"
  )

####################################### inspecting functional type grouping pr veg_type ############################################

ggplot(nmds_points, aes(x = MDS1, y = MDS2)) +
  geom_point(aes(colour = veg_type), size = 2, alpha = 0.8) +
  geom_polygon(
    data = hulls,
    aes(group = year, fill = factor(year)),
    alpha = 0.15,
    colour = "grey40"
  ) +
  geom_point(
    data = centroids,
    aes(x = MDS1, y = MDS2, fill = factor(year)),
    shape = 21,
    colour = "black",
    size = 4
  ) +
  theme_minimal() +
  labs(
    title = "NMDS of Functional Types",
    subtitle = "Points coloured by vegetation type; hulls and centroids by year",
    colour = "Vegetation type",
    fill = "Year"
  )

####################################### inspecting functional type grouping pr veg_type WITH LOADINGS ############################################

ggplot(nmds_points, aes(x = MDS1, y = MDS2)) +
  # Points by veg_type
  geom_point(aes(colour = veg_type), size = 2, alpha = 0.8) +
  
  # Convex hulls and centroids by year
  geom_polygon(
    data = hulls,
    aes(group = year, fill = factor(year)),
    alpha = 0.15,
    colour = "grey40"
  ) +
  geom_point(
    data = centroids,
    aes(x = MDS1, y = MDS2, fill = factor(year)),
    shape = 21,
    colour = "black",
    size = 4
  ) +
  
  # Loadings (arrows)
  geom_segment(
    data = species_scores,
    aes(x = 0, y = 0, xend = MDS1, yend = MDS2),
    arrow = arrow(length = unit(0.2, "cm")),
    colour = "black"
  ) +
  
  # Loadings labels
  geom_text(
    data = species_scores,
    aes(x = MDS1, y = MDS2, label = variable),
    size = 3,
    hjust = 0.5,
    vjust = -0.6
  ) +
  
  theme_minimal() +
  labs(
    title = "NMDS of Functional Types with Loadings",
    subtitle = "Points coloured by vegetation type; hulls and centroids by year",
    colour = "Vegetation type",
    fill = "Year"
  )
