#### packages ####

library(tidyverse)
library(viridis)
library(vegan)
library(remotes)
library(pairwiseAdonis)
library(patchwork)

#### importing data ####
  
merged_data <- readRDS("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/01 Vegetation changes Kobbefjord/data/nero_analysis/data/merged_data.rds") |> 
  filter(veg_type != "saltmarsh", taxon_code != "rock") |> 
  group_by(year, subsection) |> 
  mutate(no_plots = n_distinct(plot_id)) |> 
  ungroup()

species_sub_long <- merged_data |>
  group_by(year, subsection, veg_type, taxon_code, no_plots) |>
  summarize(
    count = n(),
    fraction_sub = as.double(count) / as.double(no_plots),
    .groups = "drop"
  ) |>
  filter(veg_type != "saltmarsh", no_plots != 1) |> 
  distinct(year, subsection, veg_type, taxon_code, no_plots, fraction_sub)

species_sub_wide <- species_sub_long |> 
  pivot_wider(
    names_from = taxon_code, 
    values_from = fraction_sub, 
    values_fill = 0
  )

write_rds(species_sub_wide, "data/species_sub_wide.rds")

names(species_sub_wide)

community_matrix <- read_rds("data/species_sub_wide.rds") |> 
  unite("sub_year_vt", subsection, year,veg_type, sep = "_", remove = FALSE) |> 
  column_to_rownames(var = "sub_year_vt") |> 
  select(-year, -subsection, -veg_type, -no_plots)

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

ggplot(nmds_plot_data, aes(x = NMDS1, y = NMDS2, shape = veg_type, color = factor(year))) +
  geom_point(data = centroids, aes(x = centroid_NMDS1, y = centroid_NMDS2, fill = factor(year)),
             shape = 21, size = 5, color = "black") +
  theme_minimal() +
  labs(shape = "Vegetation Type", color = "Year", fill = "Year", title = "NMDS Ordination with Year Centroids") +
  theme(legend.position = "right")

#### NMDS stats BETADISPER ################################################################### 
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

#### testing mean distance within vegtype ####

veg_types <- unique(merged_data_summary$veg_type)

test_dist_to_centroid_by_veg <-
  veg_types |>
  purrr::map_dfr(\(vt) {
    cat("Testing vegetation type:", vt, "\n")
    
    # 1) Metadata for this veg type
    plot_metadata_vt <-
      merged_data_summary |>
      filter(veg_type == vt) |>
      mutate(plot_year_vt = paste(plot_id, year, veg_type, sep = "_")) |>
      distinct(plot_year_vt, year, veg_type)
    
    # 2) Subset community matrix
    rows_to_keep <- rownames(community_matrix) %in% plot_metadata_vt$plot_year_vt
    community_matrix_vt <- community_matrix[rows_to_keep, , drop = FALSE]
    
    # Skip if too few rows
    if (nrow(community_matrix_vt) < 3) {
      return(tibble(
        veg_type = vt,
        p_value  = NA_real_,
        note     = "too few samples for NMDS"
      ))
    }
    
    # 3) NMDS for this veg type
    set.seed(42)
    nmds_result <- metaMDS(community_matrix_vt, k = 3, trymax = 100, trace = FALSE)
    
    # 4) Site scores and merge with metadata
    nmds_scores <-
      scores(nmds_result, display = "sites") |>
      as.data.frame() |>
      tibble::rownames_to_column("plot_year_vt")
    
    nmds_plot_data <-
      nmds_scores |>
      left_join(plot_metadata_vt, by = "plot_year_vt") |>
      mutate(year = factor(year))
    
    # 5) Centroids per year within this veg type
    centroids <-
      nmds_plot_data |>
      group_by(year, veg_type) |>
      summarise(
        centroid_NMDS1 = mean(NMDS1),
        centroid_NMDS2 = mean(NMDS2),
        .groups = "drop"
      )
    
    # 6) Distances to centroid
    nmds_with_centroid <-
      nmds_plot_data |>
      left_join(centroids, by = c("year", "veg_type")) |>
      mutate(
        dist_to_centroid = sqrt(
          (NMDS1 - centroid_NMDS1)^2 +
            (NMDS2 - centroid_NMDS2)^2
        )
      )
    
    # 7) One value per plot_year_vt for ANOVA
    dist_summary_vt <-
      nmds_with_centroid |>
      distinct(plot_year_vt, year, dist_to_centroid)
    
    # Guard against degenerate cases
    if (n_distinct(dist_summary_vt$year) < 2) {
      return(tibble(
        veg_type = vt,
        p_value  = NA_real_,
        note     = "only one year present"
      ))
    }
    
    dist_aov_vt <- aov(dist_to_centroid ~ factor(year), data = dist_summary_vt)
    p_val <- summary(dist_aov_vt)[[1]]["factor(year)", "Pr(>F)"]
    
    tibble(
      veg_type = vt,
      p_value  = p_val,
      note     = NA_character_
    )
  })

test_dist_to_centroid_by_veg

#plotting the distances
all_distances <-
  purrr::map_dfr(
    names(results_list),
    \(vt) {
      results_list[[vt]]$dist_data |>
        mutate(veg_type = vt)
    }
  )

ggplot(nmds_with_centroid,
       aes(x = factor(year), y = dist_to_centroid)) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1.5, color = "#076834") +
  geom_boxplot(width = 0.3, outlier.shape = NA, alpha = 0.1) +
  facet_wrap(~ veg_type, scales = "free_y") +
  labs(
    x = "Year",
    y = "Distance to centroid",
    title = "Distances to year centroid within each vegetation type"
  ) +
  theme_minimal()
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
    labs(title = paste("NMDS for Year", y)) +
    #labs(title = paste("NMDS for Year", y),
    #     shape = "Vegetation Type", color = "Vegetation Type", fill = "Vegetation Type") +
    theme(legend.position = "none")
})

# Combine plots in a grid
wrap_plots(plots, ncol = 2)  

#### NMDS year-vegtype mean distance ####

# centroids: year, veg_type, centroid_NMDS1, centroid_NMDS2

centroids |>
  summarise(
    n = n(),
    n_NA_x = sum(is.na(centroid_NMDS1)),
    n_NA_y = sum(is.na(centroid_NMDS2))
  )

centroid_pairs <-
  centroids |>
  # self-join on year: all veg_type pairs within a year
  inner_join(
    centroids,
    by = "year",
    suffix = c("1", "2")
  ) |>
  # drop symmetric duplicates and self pairs
  filter(as.character(veg_type1) < as.character(veg_type2)) |>
  mutate(
    dist = sqrt(
      (centroid_NMDS11 - centroid_NMDS12)^2 +
        (centroid_NMDS21 - centroid_NMDS22)^2
    )
  ) |>
  select(
    year,
    veg_type1,
    veg_type2,
    dist
  )

centroid_dist_summary <-
  centroid_pairs |>
  group_by(year) |>
  summarise(
    mean_dist = mean(dist),
    median_dist = median(dist),
    sd_dist = sd(dist),
    n_pairs = n(),
    .groups = "drop"
  )

centroid_dist_summary

ggplot(centroid_pairs, aes(x = factor(year), y = dist)) +
  geom_boxplot(width = 0.3, outlier.shape = NA, alpha = 0.7, size = 0.3) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 3, color = "#076834") +
  labs(
    x = "Year",
    y = "Pairwise distance between veg-type centroids",
    title = "Distances among vegetation-type centroids by year"
  )

centroid_dist_summary |>
  mutate(year_num = as.numeric(as.character(year))) |>
  with(cor.test(year_num, mean_dist, method = "kendall"))

centroid_pairs |>
  mutate(year = factor(year)) |>
  aov(dist ~ year, data = _) |>
  summary()
