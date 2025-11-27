#### loading packages ####

library(tidyverse)
library(lme4)
library(lmerTest)
library(ggplot2)
library(ggeffects)
library(vegan)
library(codyn) # turnover calculations

#### loading data ################################

merged_data <- readRDS("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/01 Vegetation changes Kobbefjord/data/nero_analysis/data/merged_data.rds") |> 
  filter(veg_type != "saltmarsh", taxon_code != "rock") |> 
  select(year, plot_id, veg_type, taxon_code, presence, func_type, ecoveg_gfc, ecoveg_sgfc)

func_count_species <- merged_data |> 
  group_by(year, plot_id, veg_type) |> 
  summarise(total_species = n_distinct(taxon_code))

func_count <- merged_data |> 
  group_by(year, plot_id, ecoveg_sgfc) |> 
  summarise(count_func = n_distinct(taxon_code))

func_plot <- func_count |> 
  left_join(func_count_species, by = c("year", "plot_id")) |> 
  mutate(frac = round(count_func/total_species,2))

unique(func_plot)

write_rds(func_plot, "data/func_plot.rds")

func_plot_wide <- func_plot |> 
  select(year,plot_id,veg_type, ecoveg_sgfc,frac) |> 
  pivot_wider(names_from = "ecoveg_sgfc", values_from = "frac", values_fill = 0)

community_matrix <- func_plot_wide |> 
  unite("plot_year_vt", plot_id, year,veg_type, sep = "_", remove = FALSE) |> 
  column_to_rownames(var = "plot_year_vt") |> 
  select(-plot_id,-year,-veg_type)

### NMDS ###########################

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
nmds_scores$plot_year_vt <- rownames(nmds_scores)

#### NMDS visualisation ######################################################################## 

# Add plot_year_vt as a column to nmds_scores for merging
nmds_scores$plot_year_vt <- rownames(nmds_scores)

plot_metadata <- func_plot |>
  mutate(plot_year_vt = paste(plot_id, year, veg_type, sep = "_")) |>
  select(plot_year_vt, year, veg_type) |> 
  distinct()

# Merge NMDS site scores with metadata using plot_year_vt as key
# Convert year to factor (if not already)
nmds_plot_data <- nmds_scores |>
  left_join(plot_metadata, by = "plot_year_vt") |>
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
  labs(shape = "Vegetation Type", color = "Year", fill = "Year", title = "NMDS Ordination with Year Centroids (PFT)") +
  theme(legend.position = "right")

ggplot(nmds_plot_data, aes(x = NMDS1, y = NMDS2, shape = veg_type, color = factor(year))) +
  geom_point(data = centroids, aes(x = centroid_NMDS1, y = centroid_NMDS2, fill = factor(year)),
             shape = 21, size = 5, color = "black") +
  theme_minimal() +
  labs(shape = "Vegetation Type", color = "Year", fill = "Year", title = "NMDS Ordination with Year Centroids (PFT)") +
  theme(legend.position = "right")

#### NMDS stats BETADISPER ######################################################################## 
# within group variation testing. Done before adonis, because significant within group variation can make the interpretation of an Adonis result more complex.

# Generate Bray-Curtis distance matrix (dist object)
bray_dist <- vegdist(community_matrix, method = "bray")

# Prepare your metadata ensuring it matches community_matrix rows exactly
metadata <- func_plot %>%
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

#### disp distance test ###################################################
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

#### pairwise PERMANOVA distance#### pairwise PERMANOVA ###################################################

# Prepare dataset, same as before
metadata_agg <-  func_plot |> 
  mutate(plot_year_vt = paste(plot_id, year, veg_type, sep = "_")) |> 
  distinct(plot_year_vt, year, veg_type) |> 
  arrange(match(plot_year_vt, rownames(community_matrix)))

nrow(metadata_agg)

stopifnot(nrow(metadata_agg) == nrow(community_matrix))
stopifnot(all(metadata_agg$plot_year_vt == rownames(community_matrix)))

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
  
