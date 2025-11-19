#### nmds with 10 agreegate ####
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
  select(year,vt_section,plot_id,veg_type,taxon_code,presence,species,ecoveg_gfc,ecoveg_sgfc) |> 
  filter(veg_type != "saltmarsh")

summary(merged_data)

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
  ggtitle("NMDS trajectories of vegetation sections over time (all vegetation types !saltmarsh)")
######################################## visualising for all veg types WITH ELIPSES AND CENTERS ########################################

# Compute centroid (mean NMDS position) per year
centroids <- nmds_data |>
  group_by(year) |>
  summarise(
    NMDS1 = mean(NMDS1, na.rm = TRUE),
    NMDS2 = mean(NMDS2, na.rm = TRUE)
  )


# Plot with ellipses + centroids
nmds_data |>
  ggplot(aes(x = NMDS1, y = NMDS2,
             color = factor(year),
             shape = factor(year))) +
  geom_point(size = 3) +
  stat_ellipse(type = "t", linetype = "solid", linewidth = 1) +
  # add centroids
  geom_point(
    data = centroids,
    aes(x = NMDS1, y = NMDS2, color = factor(year)),
    shape = 4,        # cross
    size = 5,
    stroke = 2
  ) +
  scale_color_manual(values = year_colors) +
  scale_shape_manual(values = year_shapes) +
  theme_minimal() +
  labs(shape = "Year", color = "Year") +
  ggtitle("NMDS ordination with 95% confidence ellipses and centroids per year")

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

############################################ permanova + pairvise adonis for all plots ##############################################################################

# Calculate the dissimilarity matrix (same as used for NMDS)
bray_dist <- vegdist(species_mat, method = "bray")

# Run PERMANOVA
permanova_result <- adonis2(bray_dist ~ year, data = nmds_data, permutations = 999)

# Show results
print(permanova_result)

# > print(permanova_result)
# Permutation test for adonis under reduced model
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = bray_dist ~ year, data = nmds_data, permutations = 999)
# Df SumOfSqs      R2     F Pr(>F)
# Model      1    0.071 0.00172 0.498  0.855
# Residual 289   41.001 0.99828             
# Total    290   41.071 1.00000     


# Then you can run it like:
pairwise_ad <- pairwiseAdonis::pairwise.adonis(
  x = species_mat,              # or a dist object like vegdist(species_mat)
  factors = nmds_data$year,
  sim.method = "bray",
  p.adjust.m = "bonferroni",
  perm = 999
)

print(pairwise_ad)

# pairs Df  SumsOfSqs   F.Model          R2 p.value p.adjusted sig
# 1 2007 vs 2012  1 0.13086565 0.9235351 0.006285821   0.498          1    
# 2 2007 vs 2017  1 0.05800062 0.4285292 0.002967068   0.898          1    
# 3 2007 vs 2022  1 0.11326349 0.8260349 0.005743292   0.570          1    
# 4 2012 vs 2017  1 0.08307035 0.5641961 0.003902737   0.810          1    
# 5 2012 vs 2022  1 0.03846340 0.2579899 0.001800876   0.979          1    
# 6 2017 vs 2022  1 0.09878380 0.6922204 0.004885380   0.711          1  
############################################ convex hulls visualising ##############################################################################

library(ggplot2)
library(dplyr)

# Function to get convex hull for each group (here: year)
hulls <- nmds_data |>
  group_by(year) |>
  slice(chull(NMDS1, NMDS2))  # convex hull points per year

ggplot(nmds_data, aes(x = NMDS1, y = NMDS2, color = factor(year))) +
  geom_point(size = 3, alpha = 0.7) +
  geom_polygon(data = hulls, aes(fill = factor(year), group = year), 
               alpha = 0.2, color = NA) +
  geom_path(data = hulls, aes(group = year), linewidth = 0, color = "black", alpha = 0.5) +
  scale_color_manual(values = year_colors) +
  scale_fill_manual(values = year_colors) +
  theme_minimal() +
  labs(color = "Year", fill = "Year") +
  ggtitle("NMDS convex hulls by year — visualising community spread")

############################################ convex hulls calculating ##############################################################################
library(geometry)  # for convhulln and convhulln function

# Function to compute 2D convex hull area
hull_area <- function(df) {
  pts <- df[, c("NMDS1", "NMDS2")]
  if (nrow(pts) < 3) return(NA)  # too few points
  hull <- chull(pts)
  hull_coords <- pts[hull, ]
  # polygon area using the shoelace formula
  area <- abs(sum(hull_coords$NMDS1 * dplyr::lead(hull_coords$NMDS2) - 
                    dplyr::lead(hull_coords$NMDS1) * hull_coords$NMDS2, 
                  na.rm = TRUE)) / 2
  return(area)
}

# Apply per year
hull_areas <- nmds_data |>
  group_by(year) |>
  summarise(
    n_points = n(),
    hull_area = hull_area(cur_data())
  )

print(hull_areas)

plot(hull_areas)

############################################ convex hulls sizes pr veg type pr year #############################################################
hull_areas_veg <- nmds_data |>
  group_by(veg_type, year) |>
  summarise(
    n_points = n(),
    hull_area = hull_area(cur_data())
  )

print(hull_areas_veg)

ggplot(hull_areas_veg, aes(x = as.numeric(as.character(year)), y = hull_area, group = veg_type, color = veg_type)) +
  geom_line() +
  geom_point() +
  facet_wrap(~veg_type, scales = "free_y") +
  theme_minimal() +
  labs(x = "Year", y = "Convex hull area", 
       title = "Convex hull areas per vegetation type over time")

############################################ change in dispersion #############################################################

bray_dist <- vegdist(species_mat, method = "bray")
bd <- betadisper(bray_dist, nmds_data$year)
anova(bd)
TukeyHSD(bd)

# Analysis of Variance Table
# 
# Response: Distances
# Df Sum Sq  Mean Sq F value Pr(>F)
# Groups      3 0.0056 0.001870  0.1194 0.9486
# Residuals 283 4.4305 0.015655  


####################################################################################################################################################
############################################ for each vegetation type ##############################################################################

# Shapes for years
year_shapes <- c("2007" = 16, "2012" = 15, "2017" = 17, "2022" = 18)

# Ensure chronological order for paths
nmds_data <- nmds_data |>
  arrange(vt_section, year)

# Get all unique vegetation types
veg_types <- unique(nmds_data$veg_type)

# Function to find convex hull points
find_hull <- function(df) df[chull(df$NMDS1, df$NMDS2), ]

hulls <- plot_data |>
  group_by(year) |>
  group_modify(~find_hull(.x))

plots <- list()

# Loop over veg types
for (veg in veg_types) {
  
  plot_data <- nmds_data |>
    filter(veg_type == veg)
  
  # Compute convex hulls per year
  hulls <- plot_data |>
    group_by(year) |>
    group_modify(~find_hull(.x))
  
  # Compute centroids per year
  centroids <- plot_data |>
    group_by(year) |>
    summarise(
      NMDS1 = mean(NMDS1, na.rm = TRUE),
      NMDS2 = mean(NMDS2, na.rm = TRUE)
    )
  
  # Plot
  p <- ggplot(plot_data, aes(x = NMDS1, y = NMDS2,
                             shape = factor(year),
                             color = factor(year))) +
    geom_polygon(data = hulls,
                 aes(fill = factor(year), group = year),
                 color = NA, alpha = 0.2) +
    geom_point(size = 3) +
    #stat_ellipse(type = "t", linewidth = 1) +
    geom_point(data = centroids,
               aes(x = NMDS1, y = NMDS2, color = factor(year)),
               shape = 4, size = 4, stroke = 1.5) +
    scale_shape_manual(values = year_shapes) +
    scale_color_manual(values = year_colors) +
    scale_fill_manual(values = year_colors) +
    theme_minimal() +
    labs(shape = "Year", color = "Year", fill = "Year") +
    ggtitle(paste("NMDS with convex hulls -", veg))
  
  plots[[veg]] <- p
}
plots[[veg]]

library(patchwork)

# Arrange all six plots (2 rows × 3 columns, for example)
combined_plot <- wrap_plots(plots, ncol = 2)

combined_plot

####################################################################################################################################################
######################################### pairwise adonis vegetation types #########################################################################

# Example assumptions:
# df = your main data frame
# df$vegtype = vegetation type factor
# df$year = factor with survey years
# comm = community matrix (species abundances) with same row order as df

# Run PERMANOVA by vegetation type
results_list <- list()

for (vt in unique(df$vegtype)) {
  sub_df <- df %>% filter(vegtype == vt)
  sub_comm <- comm[rownames(comm) %in% rownames(sub_df), ]
  
  # Ensure there are at least 2 years with >1 sample each
  year_counts <- table(sub_df$year)
  valid_years <- names(year_counts[year_counts > 1])
  
  if (length(valid_years) < 2) {
    message(paste("Skipping", vt, "- not enough samples per year"))
    next
  }
  
  # Subset to valid years
  sub_df <- sub_df[sub_df$year %in% valid_years, ]
  sub_comm <- sub_comm[rownames(sub_comm) %in% rownames(sub_df), ]
  
  # Pairwise comparisons
  year_pairs <- combn(valid_years, 2, simplify = FALSE)
  
  for (pair in year_pairs) {
    sub_pair_df <- sub_df[sub_df$year %in% pair, ]
    sub_pair_comm <- sub_comm[rownames(sub_comm) %in% rownames(sub_pair_df), ]
    
    if (nrow(sub_pair_df) < 4) next  # skip if too few total samples
    
    perm <- adonis2(sub_pair_comm ~ year, data = sub_pair_df, method = "bray", permutations = 999)
    
    results_list[[paste(vt, pair[1], pair[2], sep = "_")]] <- data.frame(
      vegtype = vt,
      year1 = pair[1],
      year2 = pair[2],
      F_value = perm$F[1],
      R2 = perm$R2[1],
      p_value = perm$`Pr(>F)`[1]
    )
  }
}

pairwise_results <- bind_rows(results_list)
pairwise_results
####################################################################################################################################################
######################################### quantifying direction change ##############################################################################

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
