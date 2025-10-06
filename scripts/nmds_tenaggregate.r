#### nmds with 10 agreegate ####
library(tidyverse)
library(viridis)

merged_data <- readRDS("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/01 Vegetation changes Kobbefjord/data/nero_analysis/data/merged_data.rds") |> 
  select(year,vt_section,plot_id,veg_type,taxon_code,presence,species,ecoveg_gfc)

####################################### calculating percent pr plot ############################################
plots_per_section <- merged_data |>
  distinct(year, vt_section, plot_id) |> 
  group_by(year, vt_section) |>
  summarise(n_plots_total = n(), .groups = "drop")

presence_per_plot <- merged_data |>
  group_by(year, vt_section, plot_id, species,veg_type) |>
  summarise(presence = as.integer(any(presence == 1)), .groups = "drop")

species_percent <- presence_per_plot |>
  group_by(year, vt_section, species,veg_type) |>
  summarise(n_plots_present = sum(presence), .groups = "drop") |>
  left_join(plots_per_section, by = c("year", "vt_section")) |>
  mutate(
    perc_plots = 100 * n_plots_present / n_plots_total
  ) |>
  arrange(year, vt_section, desc(perc_plots))


################################################ pivot for nmds ####################################################

# make a sample_id
species_percent2 <- species_percent |>
  mutate(sample_id = paste0(year, "-",veg_type,"_sec", vt_section))

# 1) metadata: one row per sample (year x vt_section)
metadata <- species_percent2 |>
  distinct(sample_id, year, vt_section, veg_type, n_plots_total) |>
  arrange(year, vt_section)

# 2) species matrix wide: rows = sample_id, cols = species
species_wide <- species_percent2 |>
  select(sample_id, species, perc_plots) |>
  pivot_wider(
    names_from = species,
    values_from = perc_plots,
    values_fill = 0
  ) |>
  arrange(sample_id)

# convert to matrix and set rownames
species_mat <- species_wide |>
  column_to_rownames(var = "sample_id") |>   # needs tibble::column_to_rownames or use base R
  as.matrix()

# ensure numeric (NAs filled by values_fill above, but double-check)
mode(species_mat) <- "numeric"

# quick checks
dim(species_mat)         # samples x species
head(metadata)

######################################## actual nmds ########################################

library(vegan)

set.seed(42)
nmds_bray <- metaMDS(species_mat, distance = "bray", k = 3, trymax = 200)

# extract sample scores
scores_samples <- scores(nmds_bray, display = "sites") |>
  as.data.frame() |>
  tibble::rownames_to_column("sample_id")

# join with metadata
nmds_data <- scores_samples |>
  left_join(metadata, by = "sample_id")

head(nmds_data)

######################################## visualising for all veg types ########################################
# 2D NMDS plot colored by year
nmds_data |>
  ggplot(aes(x = NMDS1, y = NMDS2, color = factor(year))) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(color = "Year") +
  ggtitle("NMDS of vegetation % occupancy (Bray-Curtis)")

# Ensure data is sorted for geom_path
nmds_data <- nmds_data |>
  arrange(vt_section, year)

# Define color palette and shapes manually for years
year_colors <- c("2007" = "#E69F00",  # orange
                 "2012" = "#56B4E9",  # blue
                 "2017" = "#009E73",  # green
                 "2022" = "#D55E00")  # red

year_shapes <- c("2007" = 16,  # circle
                 "2012" = 15,  # square
                 "2017" = 17,  # triangle
                 "2022" = 18)  # diamond

nmds_data |>
  ggplot(aes(x = NMDS1, y = NMDS2,
             shape = factor(year),
             color = factor(year))) +
  geom_point(size = 3) +
  geom_path(aes(group = vt_section), color = "gray50", linetype = "dashed") +
  scale_color_manual(values = year_colors) +
  scale_shape_manual(values = year_shapes) +
  theme_minimal() +
  labs(shape = "Year", color = "Year") +
  ggtitle("NMDS trajectories of vegetation sections over time (all vegetation types)")
######################################## visualising for all veg types WITH ELIPSES ########################################

# Ensure data is sorted for geom_path
nmds_data <- nmds_data |>
  arrange(vt_section, year)

nmds_data |>
  ggplot(aes(x = NMDS1, y = NMDS2,
             color = factor(year),
             shape = factor(year))) +
  geom_point(size = 3) +
  stat_ellipse(type = "t", linetype = "solid", linewidth = 1) +
  scale_color_manual(values = year_colors) +
  scale_shape_manual(values = year_shapes) +
  theme_minimal() +
  labs(shape = "Year", color = "Year") +
  ggtitle("NMDS ordination with 95% confidence ellipses per year")

######################################## visualising for all veg types WITH LOADINGS ########################################

# Fit species vectors
species_fit <- envfit(nmds_bray, species_mat, permutations = 999)

# Extract coordinates of the fitted vectors
species_scores <- as.data.frame(scores(species_fit, display = "vectors"))

# Clean up to only keep numeric axes
species_scores <- species_scores[, c("NMDS1", "NMDS2")]

# Add species names
species_scores$species <- rownames(species_scores)

# Base NMDS plot with ellipses
p <- nmds_data |>
  ggplot(aes(x = NMDS1, y = NMDS2, color = factor(year), shape = factor(year))) +
  geom_point(size = 3) +
  stat_ellipse(type = "t", linewidth = 1) +
  scale_color_manual(values = year_colors) +
  scale_shape_manual(values = year_shapes) +
  theme_minimal() +
  labs(shape = "Year", color = "Year") +
  ggtitle("NMDS with species vectors (envfit loadings)")

# Add species arrows and labels
p +
  geom_segment(
    data = species_scores,
    aes(x = 0, y = 0, xend = NMDS1 * 0.7, yend = NMDS2 * 0.7),
    arrow = arrow(length = unit(0.2, "cm")),
    color = "black",
    alpha = 0.7,
    inherit.aes = FALSE  # very important so it doesn't mix with NMDS point data
  ) +
  geom_text(
    data = species_scores,
    aes(x = NMDS1 * 0.8, y = NMDS2 * 0.8, label = species),
    size = 3,
    color = "black",
    alpha = 0.7,
    hjust = 0.5,
    inherit.aes = FALSE
  )

############################################ for each vegetation type ############################################

# Shapes for years
year_shapes <- c("2007" = 16, "2012" = 15, "2017" = 17, "2022" = 18)

# Ensure chronological order for paths
nmds_data <- nmds_data |>
  arrange(vt_section, year)

# Get all unique vegetation types
veg_types <- unique(nmds_data$veg_type)

# Loop over veg_types and plot
for (veg in veg_types) {
  plot_data <- nmds_data |>
    subset(veg_type == veg)
  
  p <- plot_data |>
    ggplot(aes(x = NMDS1, y = NMDS2, shape = factor(year))) +
    geom_point(size = 3, color = "black") +             # black points for clarity
    geom_path(aes(group = vt_section), color = "gray50", linetype = "dashed") +
    scale_shape_manual(values = year_shapes) +
    theme_minimal() +
    labs(shape = "Year") +
    ggtitle(paste("NMDS trajectories - Vegetation type:", veg))
  
  print(p)  # prints plot to R graphics window
}

# Define shapes for years
year_shapes <- c("2007" = 16,  # circle
                 "2012" = 15,  # square
                 "2017" = 17,  # triangle
                 "2022" = 18)  # diamond

nmds_data |>
  ggplot(aes(x = NMDS1, y = NMDS2,
             shape = factor(year),
             color = veg_type)) +
  geom_point(size = 3) +
  geom_path(aes(group = vt_section), color = "gray50", linetype = "dashed") +
  scale_shape_manual(values = year_shapes) +
  scale_color_viridis(discrete = TRUE, option = "D") +  # color-blind-friendly
  theme_minimal() +
  labs(shape = "Year", color = "Vegetation type") +
  ggtitle("NMDS trajectories of vegetation sections over time")

#### quantifying direction change ####

# Ensure data is sorted by section and year
nmds_data <- nmds_data |>
  arrange(vt_section, year)

# Compute first-to-last year vector for each section
vectors <- nmds_data |>
  group_by(vt_section) |>
  summarise(
    start_NMDS1 = first(NMDS1),
    start_NMDS2 = first(NMDS2),
    end_NMDS1   = last(NMDS1),
    end_NMDS2   = last(NMDS2),
    .groups = "drop"
  ) |>
  mutate(
    delta_NMDS1 = end_NMDS1 - start_NMDS1,
    delta_NMDS2 = end_NMDS2 - start_NMDS2,
    magnitude   = sqrt(delta_NMDS1^2 + delta_NMDS2^2),
    angle       = atan2(delta_NMDS2, delta_NMDS1)  # in radians
  )
vectors

library(circular)

angles <- circular(vectors$angle)
mean_direction <- mean(angles)             # average vector direction
result <- rayleigh.test(angles)

library(circular)
library(ggplot2)

# angles is already a circular object
# convert radians to degrees for easier interpretation
angles_deg <- as.numeric(angles) * 180 / pi

# simple histogram
tibble::tibble(angle_deg = angles_deg) |>
  ggplot(aes(x = angle_deg)) +
  geom_histogram(binwidth = 20, fill = "skyblue", color = "black") +
  coord_polar(start = 0) +  # polar coordinates for circular plot
  theme_minimal() +
  labs(title = "Distribution of NMDS change directions", x = "Direction (degrees)", y = "Count")

summary(vectors$magnitude)
hist(vectors$magnitude)

### plots for each vegetation type ###
arrows_mid <- nmds_data |>
  arrange(vt_section, year) |>
  group_by(vt_section) |>
  mutate(
    NMDS1_next = lead(NMDS1),
    NMDS2_next = lead(NMDS2)
  ) |>
  filter(!is.na(NMDS1_next)) |>
  mutate(
    # midpoint coordinates
    NMDS1_mid = (NMDS1 + NMDS1_next)/2,
    NMDS2_mid = (NMDS2 + NMDS2_next)/2,
    # tiny offset for arrow direction
    dx = (NMDS1_next - NMDS1) * 0.05,  
    dy = (NMDS2_next - NMDS2) * 0.05
  ) |>
  ungroup()

library(ggplot2)
library(grid)

year_shapes <- c("2007" = 16, "2012" = 15, "2017" = 17, "2022" = 18)
veg_types <- unique(nmds_data$veg_type)

for (veg in veg_types) {
  plot_points <- nmds_data |>
    subset(veg_type == veg)
  
  plot_arrows <- arrows_mid |>
    subset(veg_type == veg)
  
  p <- plot_points |>
    ggplot(aes(x = NMDS1, y = NMDS2, shape = factor(year))) +
    geom_point(size = 3, color = "black") +
    # full dashed line
    geom_path(aes(group = vt_section), color = "gray50", linetype = "dashed") +
    # small arrow at midpoint
    geom_segment(data = plot_arrows,
                 aes(x = NMDS1_mid, y = NMDS2_mid,
                     xend = NMDS1_mid + dx, yend = NMDS2_mid + dy),
                 arrow = arrow(type = "open", length = unit(0.15, "inches")),
                 color = "gray50") +
    scale_shape_manual(values = year_shapes) +
    theme_minimal() +
    labs(shape = "Year") +
    ggtitle(paste("NMDS trajectories - Vegetation type:", veg))
  
  print(p)
}
