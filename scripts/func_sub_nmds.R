#### packages ####

library(tidyverse)
library(vegan)
library(dplyr)
library(ggplot2)
library(purrr)
library(lme4)
library(lmerTest)
library(ggeffects)
library(patchwork)

#### importing data ####

merged_data <- readRDS("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/01 Vegetation changes Kobbefjord/data/nero_analysis/data/merged_data.rds") |> 
  group_by(year, subsection) |> 
  mutate(no_plots = n_distinct(plot_id)) |> 
  ungroup() |> 
  filter(taxon_code != "rock", no_plots > 1)

species_sub_long <- merged_data |>
  group_by(year, subsection, veg_type, taxon_code, no_plots, ecoveg_sgfc) |>
  reframe(
    count = n(),
    fraction_sub = as.double(count) / as.double(no_plots),
    .groups = "drop"
  ) |>
  filter(veg_type != "saltmarsh", no_plots != 1) |> 
  distinct(year, subsection, veg_type, ecoveg_sgfc, taxon_code, no_plots, fraction_sub)

unique(func_sub_frac_sum$ecoveg_sgfc)

func_sub_frac_sum <- species_sub_long |> 
  group_by(year, subsection, veg_type, ecoveg_sgfc) |> 
  reframe(frac_sum = sum(fraction_sub))

func_sub_wide <- func_sub_frac_sum |> 
  pivot_wider(names_from = ecoveg_sgfc, values_from = frac_sum, values_fill = 0) |> 
  rename(graminoid = herb_graminoid, bryophyte = non_vascular_bryophyte, decidous = shrub_decidous, evergreen = shrub_evergreen, lichen = non_vascular_lichen)

community_matrix <- func_sub_wide |> 
  unite("sub_year_vt", subsection, year,veg_type, sep = "_", remove = FALSE) |> 
  column_to_rownames(var = "sub_year_vt") |> 
  select(-year, -subsection, -veg_type)

#### shrub fraction ####

names(func_sub_frac_sum)

lmm_func_sub_shrub <- lmer(frac_sum ~ year + (1 | subsection), data = func_sub_frac_sum |> filter(ecoveg_sgfc == "shrub_decidous"))
summary(lmm_func_sub_shrub)

lmm_func_sub_shrub_factor <- lmer(frac_sum ~ factor(year) + (1 | subsection), data = func_sub_frac_sum |> filter(ecoveg_sgfc == "shrub_decidous"))
summary(lmm_func_sub_shrub_factor)


pred_frac_sum <- ggpredict(lmm_func_sub_shrub, terms = c("year"))

pred_frac_sum_plot <- pred_frac_sum  |> 
  mutate(year_num = as.numeric(x))

ggplot(func_sub_frac_sum |> filter(ecoveg_sgfc == "shrub_decidous"), aes(x = factor(year), y = frac_sum)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.35, size = 1.1, color = "#076834") +
  geom_ribbon(
    data = pred_frac_sum_plot,
    aes(x = factor(x), ymin = conf.low, ymax = conf.high, group = 1),
    inherit.aes = FALSE,
    fill = "#076834",
    alpha = 0.15
  ) +
  geom_line(
    data = pred_frac_sum,
    aes(x = factor(x), y = predicted, group = 1),
    inherit.aes = FALSE,
    color = "#076834",
    size = 1.2
  ) +
  labs(
    x = "Year",
    y = "Summed fraction of deciduous shrub",
    title = "Summed of deciduous shrub in subsection"
  ) +
  theme_minimal(base_size = 13)

#### NMDS ######################################################################

# Set seed for reproducibility
set.seed(42)

# Run NMDS on presence/absence community matrix
nmds_result <- metaMDS(community_matrix, k = 3, trymax = 100)

print(nmds_result)

nmds_scores <- as.data.frame(scores(nmds_result, display = "sites"))

nrow(community_matrix)
nrow(nmds_scores) 

# Extract NMDS scores for plots
nmds_scores <- as.data.frame(scores(nmds_result, display = "sites"))
nmds_scores$sub_year_vt <- rownames(nmds_scores)

#### NMDS visualisation ######################################################################## 

# Add plot_year_vt as a column to nmds_scores for merging
nmds_scores$sub_year_vt <- rownames(nmds_scores)

plot_metadata <- merged_data |>
  mutate(sub_year_vt = paste(subsection, year, veg_type, sep = "_")) |>
  select(sub_year_vt, year, veg_type) |> 
  distinct()

# Merge NMDS site scores with metadata using plot_year_vt as key
# Convert year to factor (if not already)
nmds_plot_data <- nmds_scores |>
  left_join(plot_metadata, by = "sub_year_vt") |>
  mutate(year = factor(year))  # convert year to factor for plotting

# Define function to safely compute convex hull
find_hull <- function(df) {
  df <- na.omit(df)  # remove NA rows
  if(nrow(df) < 3) {  # chull requires at least 3 points
    return(df)
  } else {
    df[chull(df$NMDS1, df$NMDS2), ]
  }
}

# Calculate hulls grouping strictly by year
hulls <- nmds_plot_data |>
  group_by(year) |>
  do(find_hull(.))

hulls_closed <- hulls |>
  group_by(year) |>
  do({
    df <- .
    if(nrow(df) > 1) {
      df <- rbind(df, df[1, ])  # Add first row at the end to close polygon
    }
    df
  })

centroids <- nmds_plot_data |>
  group_by(year) |>
  summarize(centroid_NMDS1 = mean(NMDS1),
            centroid_NMDS2 = mean(NMDS2))

# Generate NMDS plot
ggplot(nmds_plot_data, aes(x = NMDS1, y = NMDS2, shape = veg_type, color = factor(year))) +
  geom_polygon(data = hulls_closed, aes(fill = factor(year), group = factor(year)), alpha = 0.2, color = NA) +
  geom_path(data = hulls_closed, aes(group = factor(year)), color = "white", size = 0.2) +
  geom_point(size = 3, alpha = 0.7) +
  geom_point(data = centroids, aes(x = centroid_NMDS1, y = centroid_NMDS2, fill = factor(year)),
             shape = 21, size = 5, color = "black") +
  theme_minimal() +
  labs(shape = "Vegetation Type", color = "Year", fill = "Year", title = "NMDS Ordination with Year Centroids") +
  theme(legend.position = "right")

ggplot(nmds_plot_data, aes(x = NMDS1, y = NMDS2, color = factor(year))) +
  # centroids
  geom_point(
    data = centroids,
    aes(x = centroid_NMDS1, y = centroid_NMDS2, fill = factor(year)),
    shape = 21, size = 5, color = "black", show.legend = FALSE
  ) +
  geom_text(
    data = centroids,
    aes(x = centroid_NMDS1, y = centroid_NMDS2, label = year),
    vjust = -1.2, size = 4, color = "black", show.legend = FALSE
  ) +
  theme_minimal() +
  labs(
    title = "NMDS Ordination with Year Centroids"
  ) +
  guides(
    color = "none",
    fill  = "none"
  ) +
  theme(legend.position = "right")

#### NMDS year-vegtype #####

# Function to compute convex hull safely
find_hull <- function(df) {
  df <- na.omit(df)
  if (nrow(df) < 3) return(df)
  df[chull(df$NMDS1, df$NMDS2), ]
}

# Compute hulls per year and veg_type
hulls <- nmds_plot_data %>%
  group_by(year, veg_type) %>%
  do(find_hull(.))

# Close hull polygons
hulls_closed <- hulls %>%
  group_by(year, veg_type) %>%
  do({
    df <- .
    if (nrow(df) > 1) df <- rbind(df, df[1, ])
    df
  })

# Compute centroids per year and veg_type
centroids <- nmds_plot_data %>%
  group_by(year, veg_type) %>%
  summarize(centroid_NMDS1 = mean(NMDS1),
            centroid_NMDS2 = mean(NMDS2),
            .groups = "drop")

# Create separate plots for each year
years <- unique(nmds_plot_data$year)

plots <- lapply(years, function(y) {
  df <- nmds_plot_data %>% filter(year == y)
  hull_df <- hulls_closed %>% filter(year == y)
  cent_df <- centroids %>% filter(year == y)
  
  ggplot(df, aes(x = NMDS1, y = NMDS2, shape = veg_type, color = veg_type)) +
    geom_polygon(data = hull_df, aes(fill = veg_type, group = veg_type), alpha = 0.2, color = NA) +
    geom_path(data = hull_df, aes(group = veg_type), color = "white", size = 0.3) +
    geom_point(size = 3, alpha = 0.7) +
    geom_point(data = cent_df, aes(x = centroid_NMDS1, y = centroid_NMDS2, fill = veg_type),
               shape = 21, size = 5, color = "black") +
    geom_text(data = cent_df, aes(x = centroid_NMDS1, y = centroid_NMDS2, label = veg_type),
              vjust = -1.2, size = 3, color = "black") +
    theme_minimal() +
    labs(title = paste("NMDS for Year", y),
         shape = "Vegetation Type", color = "Vegetation Type", fill = "Vegetation Type") +
    theme(legend.position = "right")
})

# Combine plots in a grid
wrap_plots(plots, ncol = 2)

#### NMDS stats BETADISPER ######################################################################## 
# within group variation testing. Done before adonis, because significant within group variation can make the interpretation of an Adonis result more complex.

# Generate Bray-Curtis distance matrix (dist object)
bray_dist <- vegdist(community_matrix, method = "bray")

# Prepare your metadata ensuring it matches community_matrix rows exactly
metadata <- merged_data %>%
  mutate(sub_year_vt = paste(subsection, year, veg_type, sep = "_")) %>%
  filter(sub_year_vt %in% rownames(community_matrix)) %>%
  arrange(match(sub_year_vt, rownames(community_matrix))) |> 
  distinct(sub_year_vt, year, veg_type)

# Check row counts
nrow(metadata)
nrow(community_matrix)
attr(bray_dist, "Labels")   # sample names used in bray_dist

# Verify matching identifiers
all(rownames(community_matrix) == metadata$sub_year_vt)  # should be TRUE
all(attr(bray_dist, "Labels") == rownames(community_matrix))  # should be TRUE

# Calculate Bray-Curtis distance on presence/absence community matrix
bray_dist <- vegdist(community_matrix, method = "bray")

# Perform betadisper to calculate multivariate dispersions by year
dispersion <- betadisper(bray_dist, group = as.factor(metadata$year))

# Test for significant differences in dispersion among years (permutation test)
permutest_result <- permutest(dispersion)

print(permutest_result)

# Optional: plot the dispersions
plot(dispersion, hull = TRUE, ellipse = TRUE, main = "Beta Dispersion by Year")

#### disp distance test ###################################################
dispersion      # your betadisper result
distances <- dispersion$distances  # vector of distances
groups <- metadata_agg$year

disp_df <- data.frame(
  year = groups,
  distance = as.numeric(dispersion$distances)
)

disp_df$sub_year_vt <- metadata_agg$sub_year_vt

disp_df <- disp_df %>%
  left_join(metadata_agg %>% select(sub_year_vt, subsection, veg_type), by = "sub_year_vt")

# Example model: distance explained by year, with random intercept for section
names(disp_df)

str(disp_df$year)
lmm_func_sub <- lmer(distance ~ year + (1 | subsection), data = disp_df)

summary(lmm_func_sub)

# Generate predicted values across observed years
pred_distance <- ggpredict(lmm_func_sub, terms = c("year"))

# Prepare predicted data for plotting: convert x (year) to numeric if needed
pred_distance_plot <- pred_distance %>%
  mutate(year_num = as.numeric(x))  # x is factor or character representing years

# Plot observed distances (boxplot + jitter) + predicted line + confidence ribbon
ggplot(disp_df, aes(x = factor(year), y = distance)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.35, size = 1.1, color = "#076834") +
  geom_ribbon(
    data = pred_distance_plot,
    aes(x = factor(x), ymin = conf.low, ymax = conf.high, group = 1),
    inherit.aes = FALSE,
    fill = "#076834",
    alpha = 0.15
  ) +
  geom_line(
    data = pred_distance_plot,
    aes(x = factor(x), y = predicted, group = 1),
    inherit.aes = FALSE,
    color = "#076834",
    size = 1.2
  ) +
  labs(
    x = "Year",
    y = "Distance to centroid",
    title = "Distance to year centroid (observed + predicted)"
  ) +
  theme_minimal(base_size = 13)
#### disp distance test FACTOR ####

lmm_func_sub_factor <- lmer(distance ~ factor(year) + (1 | subsection), data = disp_df)
summary(lmm_func_sub)

#### pairwise PERMANOVA ###################################################

# Prepare dataset, same as before
metadata_agg <- species_sub_long |> 
  mutate(sub_year_vt = paste(subsection, year, veg_type, sep = "_")) |> 
  distinct(sub_year_vt, year, veg_type, subsection) |> 
  arrange(match(sub_year_vt, rownames(community_matrix)))

nrow(metadata_agg)

stopifnot(nrow(metadata_agg) == nrow(community_matrix))
stopifnot(all(metadata_agg$sub_year_vt == rownames(community_matrix)))

# Function for pairwise PERMANOVA between levels of a factor
pairwise_adonis <- function(dist_matrix, factors, permutations = 999) {
  levels <- unique(factors)
  results <- list()
  for(i in 1:(length(levels)-1)) {
    for(j in (i+1):length(levels)) {
      group1 <- levels[i]
      group2 <- levels[j]
      idx <- which(factors %in% c(group1, group2))
      sub_dist <- as.dist(as.matrix(dist_matrix)[idx, idx])
      sub_factors <- droplevels(factors[idx])
      ad <- adonis2(sub_dist ~ sub_factors, permutations = permutations)
      res <- data.frame(
        group1 = group1,
        group2 = group2,
        F.Model = ad$F[1],
        R2 = ad$R2[1],
        p.value = ad$`Pr(>F)`[1]
      )
      results[[paste(group1, group2, sep = "_vs_")]] <- res
    }}
  do.call(rbind, results)
}

# Run pairwise PERMANOVA on years
pairwise_results <- pairwise_adonis(bray_dist, as.factor(metadata_agg$year))
print(pairwise_results)

#### NMDS for each vegetation type ######################################################################
#for each vegetation type

# Get unique vegetation types
veg_types <- unique(func_sub_wide$veg_type)

# Function to run NMDS and generate plot for one veg_type
run_nmds_for_veg <- function(vt) {
  cat("Processing vegetation type:", vt, "\n")
  
  # Filter metadata and community matrix rows for this vegetation type
  plot_metadata_vt <- func_sub_wide %>%
    filter(veg_type == vt) %>%
    mutate(sub_year_vt = paste(subsection, year, veg_type, sep = "_")) %>%
    distinct(sub_year_vt, year, veg_type)
  
  # Filter community matrix by matching rows (plot_year_vt)
  rows_to_keep <- rownames(community_matrix) %in% plot_metadata_vt$sub_year_vt
  community_matrix_vt <- community_matrix[rows_to_keep, ]
  
  # Run NMDS on subset
  set.seed(42)
  nmds_result <- metaMDS(community_matrix_vt, k = 3, trymax = 100)
  cat("Stress:", nmds_result$stress, "\n")
  
  # Extract site scores
  nmds_scores <- as.data.frame(scores(nmds_result, display = "sites"))
  nmds_scores$sub_year_vt <- rownames(nmds_scores)
  
  # Merge with metadata
  nmds_plot_data <- nmds_scores %>%
    left_join(plot_metadata_vt, by = "sub_year_vt") %>%
    mutate(year = factor(year))
  
  # Convex hull function
  find_hull <- function(df) {
    df <- na.omit(df)
    if(nrow(df) < 3) return(df)
    df[chull(df$NMDS1, df$NMDS2), ]
  }
  
  # Hulls by year
  hulls <- nmds_plot_data %>%
    group_by(year) %>%
    do(find_hull(.))
  
  hulls_closed <- hulls %>%
    group_by(year) %>%
    do({
      df <- .
      if(nrow(df) > 1) df <- rbind(df, df[1, ])
      df
    })
  
  # Calculate centroids
  centroids <- nmds_plot_data %>%
    group_by(year) %>%
    summarize(centroid_NMDS1 = mean(NMDS1),
              centroid_NMDS2 = mean(NMDS2))
  
  # Generate plot
  p <- ggplot(nmds_plot_data, aes(x = NMDS1, y = NMDS2, color = factor(year), shape = factor(year))) +
    geom_polygon(data = hulls_closed, aes(fill = factor(year), group = factor(year)), alpha = 0.2, color = NA) +
    geom_path(data = hulls_closed, aes(group = factor(year)), color = "white", size = 0.2) +
    geom_point(size = 3, alpha = 0.7) +                                # points colored and shaped by year
    geom_point(data = centroids,
               aes(x = centroid_NMDS1, y = centroid_NMDS2),
               shape = 21, size = 7, color = "black", fill = NA, stroke = 1) +
    geom_point(data = centroids, 
               aes(x = centroid_NMDS1, y = centroid_NMDS2, 
                   color = factor(year), shape = factor(year), fill = factor(year)),
               size = 5, stroke = 1) +  # stroke controls border thickness
    theme_minimal() +
    labs(shape = "Year", color = "Year", fill = "Year",
         title = paste("NMDS Year Centroids for", vt)) +
    theme(legend.position = "inside", legend.justification.inside = c(0, 0))
  print(p)
  
  return(list(nmds_result = nmds_result, plot = p))
}

# Run for all veg_types
results_list <- map(veg_types, run_nmds_for_veg)
names(results_list) <- veg_types
# Extract just the ggplot objects
plots <- lapply(results_list, `[[`, "plot")

# Wrap them into a grid with 3 columns
wrap_plots(plots, ncol = 3)
