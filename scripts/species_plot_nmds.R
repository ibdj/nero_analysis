#### packages ####

library(tidyverse)
library(viridis)
library(vegan)
library(remotes)
library(pairwiseAdonis)
library(patchwork)
library(dplyr)
library(ggplot2)
library(purrr)

#### importing data ##########################################################
merged_data <- readRDS("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/01 Vegetation changes Kobbefjord/data/nero_analysis/data/merged_data.rds")

names(merged_data)

merged_data <- merged_data |>
  filter(presence == 1) |> 
  filter(!is.na(raunkiaer_value)) |>   # Remove rows with missing values
  filter(taxon_code != "rock", veg_type != "saltmarsh")

str(merged_data$presence)

merged_data_summary <- merged_data |>
  group_by(plot_id, year, veg_type, taxon_code) |>
  summarize(presence = max(presence), .groups = 'drop') |>
  unite("plot_year_vt", plot_id, year,veg_type, sep = "_", remove = FALSE)

community_matrix <- merged_data_summary |>
  select(plot_year_vt, taxon_code, presence) |>
  pivot_wider(names_from = taxon_code, values_from = presence, values_fill = 0) |>
  column_to_rownames(var = "plot_year_vt")

community_matrix <- community_matrix |> 
  filter(if_any(where(is.numeric), ~ . != 0))

#### nmds ############ 

# Set seed for reproducibility
set.seed(42)

# Run NMDS on presence/absence community matrix
nmds_result <- metaMDS(community_matrix, k = 3, trymax = 100)

# Check stress value and results
print(nmds_result)

nmds_scores <- as.data.frame(scores(nmds_result, display = "sites"))

nrow(community_matrix)
nrow(nmds_scores) 

site_scores <- as.data.frame(scores(nmds_result, display = "sites"))    # 2190 rows
species_scores <- as.data.frame(scores(nmds_result, display = "species"))  # 91 rows
site_scores$plot_year_vt <- rownames(site_scores)

# Extract NMDS scores for plots
nmds_scores <- as.data.frame(scores(nmds_result, display = "sites"))
nmds_scores$plot_year_vt <- rownames(nmds_scores)

# Add plot_year_vt as row names for easier merging with metadata if needed
nmds_scores$plot_year_vt <- rownames(nmds_scores)

# •	Stress < 0.05: Excellent representation, very reliable ordination
# •	Stress < 0.1: Good representation, generally acceptable
# •	Stress < 0.2: Fair or usable representation, but some distortion expected
# •	Stress > 0.2: Poor representation, ordination may not be reliable or interpretable

#### NMDS visualisation ######################################################################## 

# Add plot_year_vt as a column to nmds_scores for merging
nmds_scores$plot_year_vt <- rownames(nmds_scores)

plot_metadata <- merged_data_summary %>%
  mutate(plot_year_vt = paste(plot_id, year, veg_type, sep = "_")) %>%
  select(plot_year_vt, year, veg_type)

# Merge NMDS site scores with metadata using plot_year_vt as key
# Convert year to factor (if not already)
nmds_plot_data <- nmds_scores %>%
  left_join(plot_metadata, by = "plot_year_vt") %>%
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

ggplot(nmds_plot_data, aes(x = NMDS1, y = NMDS2, shape = veg_type, color = factor(year))) +
  geom_point(data = centroids, aes(x = centroid_NMDS1, y = centroid_NMDS2, fill = factor(year)),
             shape = 21, size = 5, color = "black") +
  theme_minimal() +
  labs(shape = "Vegetation Type", color = "Year", fill = "Year", title = "NMDS Ordination with Year Centroids") +
  theme(legend.position = "right")

#### NMDS stats BETADISPER ######################################################################## 
# within group variation testing. Done before adonis, because significant within group variation can make the interpretation of an Adonis result more complex.

# Generate Bray-Curtis distance matrix (dist object)
bray_dist <- vegdist(community_matrix, method = "bray")

# Prepare your metadata ensuring it matches community_matrix rows exactly
metadata <- merged_data_summary %>%
  mutate(plot_year_vt = paste(plot_id, year, veg_type, sep = "_")) %>%
  filter(plot_year_vt %in% rownames(community_matrix)) %>%
  arrange(match(plot_year_vt, rownames(community_matrix))) |> 
  distinct(plot_year_vt, year, veg_type)

# Check row counts
nrow(metadata)
nrow(community_matrix)
attr(bray_dist, "Labels")   # sample names used in bray_dist

# Verify matching identifiers
all(rownames(community_matrix) == metadata$plot_year_vt)  # should be TRUE
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

#### disp distance test !!!! NOT ADJUSTED!! ###################################################
dispersion      # your betadisper result
distances <- dispersion$distances  # vector of distances
groups <- metadata_agg$year

disp_df <- data.frame(
  year = groups,
  distance = as.numeric(dispersion$distances)
)

disp_df$plot_year_vt <- metadata_agg$plot_year_vt

disp_df <- disp_df %>%
  left_join(metadata_agg %>% select(plot_year_vt, plot_id, veg_type), by = "plot_year_vt")

# Example model: distance explained by year, with random intercept for section

str(disp_df$year)
lmm_func_plot <- lmer(distance ~ year.x + (1 | plot_id), data = disp_df)

summary(lmm_func_plot)

# Generate predicted values across observed years
pred_distance <- ggpredict(lmm_func_plot, terms = c("year.x"))

# Prepare predicted data for plotting: convert x (year) to numeric if needed
pred_distance_plot <- pred_distance %>%
  mutate(year_num = as.numeric(x))  # x is factor or character representing years

# Plot observed distances (boxplot + jitter) + predicted line + confidence ribbon
ggplot(disp_df, aes(x = factor(year.x), y = distance)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.35, size = 1.1, color = "#076834") +
  geom_ribbon(
    data = pred_distance_plot,!!!
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
#### ADONIS ###################################################
library(vegan)
library(dplyr)

# Ensure metadata matches community_matrix rows exactly
metadata_agg <- merged_data_summary %>%
  mutate(plot_year_vt = paste(plot_id, year, veg_type, sep = "_")) %>%
  distinct(plot_year_vt, year, veg_type) %>%
  arrange(match(plot_year_vt, rownames(community_matrix)))

# Confirm alignment before running PERMANOVA
stopifnot(all(rownames(community_matrix) == metadata_agg$plot_year_vt))

# Run PERMANOVA testing for differences in community composition between years
adonis_result <- adonis2(community_matrix ~ year, data = metadata_agg, method = "bray")

print(adonis_result)

#### pairwise PERMANOVA ###################################################

library(vegan)
library(dplyr)

# Prepare dataset, same as before
metadata_agg <- merged_data_summary %>%
  mutate(plot_year_vt = paste(plot_id, year, veg_type, sep = "_")) %>%
  distinct(plot_year_vt, year, veg_type) %>%
  arrange(match(plot_year_vt, rownames(community_matrix)))

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
    }
  }
  do.call(rbind, results)
}

# Run pairwise PERMANOVA on years
pairwise_results <- pairwise_adonis(bray_dist, as.factor(metadata_agg$year))
print(pairwise_results)





#### NMDS for each vegetation type ######################################################################

# Get unique vegetation types
veg_types <- unique(merged_data_summary$veg_type)

# Function to run NMDS and generate plot for one veg_type
run_nmds_for_veg <- function(vt) {
  cat("Processing vegetation type:", vt, "\n")
  
  # Filter metadata and community matrix rows for this vegetation type
  plot_metadata_vt <- merged_data_summary %>%
    filter(veg_type == vt) %>%
    mutate(plot_year_vt = paste(plot_id, year, veg_type, sep = "_")) %>%
    distinct(plot_year_vt, year, veg_type)
  
  # Filter community matrix by matching rows (plot_year_vt)
  rows_to_keep <- rownames(community_matrix) %in% plot_metadata_vt$plot_year_vt
  community_matrix_vt <- community_matrix[rows_to_keep, ]
  
  # Run NMDS on subset
  set.seed(42)
  nmds_result <- metaMDS(community_matrix_vt, k = 3, trymax = 100)
  cat("Stress:", nmds_result$stress, "\n")
  
  # Extract site scores
  nmds_scores <- as.data.frame(scores(nmds_result, display = "sites"))
  nmds_scores$plot_year_vt <- rownames(nmds_scores)
  
  # Merge with metadata
  nmds_plot_data <- nmds_scores %>%
    left_join(plot_metadata_vt, by = "plot_year_vt") %>%
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
