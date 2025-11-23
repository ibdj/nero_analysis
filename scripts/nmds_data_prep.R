library(tidyverse)
library(vegan)

#### Summarize duplicate rows before pivoting ##########################################################
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

nrow(community_matrix)  # Should be 91
nrow(nmds_scores) 

site_scores <- as.data.frame(scores(nmds_result, display = "sites"))    # 2190 rows
species_scores <- as.data.frame(scores(nmds_result, display = "species"))  # 91 rows
site_scores$plot_year_vt <- rownames(site_scores)

# Extract NMDS scores for plots
nmds_scores <- as.data.frame(scores(nmds_result))

# Add plot_year_vt as row names for easier merging with metadata if needed
nmds_scores$plot_year_vt <- rownames(nmds_scores)

# •	Stress < 0.05: Excellent representation, very reliable ordination
# •	Stress < 0.1: Good representation, generally acceptable
# •	Stress < 0.2: Fair or usable representation, but some distortion expected
# •	Stress > 0.2: Poor representation, ordination may not be reliable or interpretable

#### NMDS visualisation ######################################################################## 

# Add plot_year_vt as a column to nmds_scores for merging
nmds_scores$plot_year_vt <- rownames(nmds_scores)

# Merge NMDS site scores with metadata using plot_year_vt as key
# Convert year to factor (if not already)
nmds_plot_data <- nmds_plot_data |>
  mutate(year = factor(year))

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


#### NMDS stats BETADISPER ######################################################################## 
# within group variation testing. Done before adonis, because significant within group variation can make the interpretation of an Adonis result more complex.

# Check row counts
nrow(metadata)
nrow(community_matrix)
attr(bray_dist, "Labels")   # sample names used in bray_dist

# Verify matching identifiers
all(rownames(community_matrix) == metadata$plot_year_vt)  # should be TRUE
all(attr(bray_dist, "Labels") == rownames(community_matrix))  # should be TRUE

# Ensure your metadata matches the rows of community_matrix
metadata <- merged_data_summary |>
  mutate(plot_year_vt = paste(plot_id, year, veg_type, sep = "_")) |>
  filter(plot_year_vt %in% rownames(community_matrix)) |>
  arrange(match(plot_year_vt, rownames(community_matrix)))

# Calculate Bray-Curtis distance on presence/absence community matrix
bray_dist <- vegdist(community_matrix, method = "bray")

# Perform betadisper to calculate multivariate dispersions by year
dispersion <- betadisper(bray_dist, group = as.factor(metadata$year))

# Test for significant differences in dispersion among years (permutation test)
permutest_result <- permutest(dispersion)

print(permutest_result)

# Optional: plot the dispersions
plot(dispersion, hull = TRUE, ellipse = TRUE, main = "Beta Dispersion by Year")

