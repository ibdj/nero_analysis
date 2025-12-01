#### packages ####

library(tidyverse)
library(vegan)
library(dplyr)
library(ggplot2)
library(purrr)

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

#### pairwise PERMANOVA ###################################################

# Prepare dataset, same as before
metadata_agg <- species_sub_long |> 
  mutate(sub_year_vt = paste(subsection, year, veg_type, sep = "_")) |> 
  distinct(sub_year_vt, year, veg_type) |> 
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
veg_types <- unique(species_sub_long$veg_type)

# Function to run NMDS and generate plot for one veg_type
run_nmds_for_veg <- function(vt) {
  cat("Processing vegetation type:", vt, "\n")
  
  # Filter metadata and community matrix rows for this vegetation type
  plot_metadata_vt <- species_sub_long %>%
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
  p <- ggplot(nmds_plot_data, aes(x = NMDS1, y = NMDS2, shape = veg_type, color = factor(year))) +
    geom_polygon(data = hulls_closed, aes(fill = factor(year), group = factor(year)), alpha = 0.2, color = NA) +
    geom_path(data = hulls_closed, aes(group = factor(year)), color = "white", size = 0.2) +
    geom_point(size = 3, alpha = 0.7) +
    geom_point(data = centroids, aes(x = centroid_NMDS1, y = centroid_NMDS2, fill = factor(year)),
               shape = 21, size = 5, color = "black") +
    theme_minimal() +
    labs(shape = "Vegetation Type", color = "Year", fill = "Year",
         title = paste("NMDS Ordination with Year Centroids for", vt)) +
    theme(legend.position = "right")
  
  print(p)
  
  return(list(nmds_result = nmds_result, plot = p))
}

# Run for all veg_types
results_list <- map(veg_types, run_nmds_for_veg)
names(results_list) <- veg_types