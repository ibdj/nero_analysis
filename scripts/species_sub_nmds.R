#### packages ####
library(remotes)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

library(tidyverse)
library(viridis)
library(vegan)
library(remotes)
library(pairwiseAdonis)
library(patchwork)
library(lme4)
library(lmerTest)
devtools::install_github("rvlenth/emmeans")
library(emmeans)


#### importing data ############################################################
  
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

community_matrix <- species_sub_wide |> 
  unite("sub_year_vt", subsection, year, veg_type, sep = "_", remove = FALSE) |> 
  as.data.frame() |> 
  column_to_rownames(var = "sub_year_vt") |> 
  dplyr::select(-year, -subsection, -veg_type, -no_plots)

meta <- merged_data |>
  mutate(sub_year_vt = paste(subsection, year, veg_type, sep = "_")) |>
  dplyr::select(sub_year_vt, year, veg_type, subsection) |> 
  distinct()

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


#### NMDS visualisation ########################################################

# Add plot_year_vt as a column to nmds_scores for merging
nmds_scores$sub_year_vt <- rownames(nmds_scores)

# Merge NMDS site scores with metadata using plot_year_vt as key
# Convert year to factor (if not already)
nmds_plot_data <- nmds_scores |>
  left_join(meta, by = "sub_year_vt") |>
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
  geom_path(data = hulls_closed, aes(group = factor(year)), color = "white", linewidth = 0.2) +
  geom_point(size = 3, alpha = 0.7) +
  geom_point(data = centroids, aes(x = centroid_NMDS1, y = centroid_NMDS2, fill = factor(year)),
             shape = 21, size = 5, color = "black") +
  theme_minimal() +
  labs(shape = "Vegetation Type", color = "Year", fill = "Year", title = "NMDS Ordination with Year Centroids") +
  theme(legend.position = "right")

ggplot(nmds_plot_data, aes(x = NMDS1, y = NMDS2, color = factor(year))) +
  geom_point(
    data = centroids,
    aes(x = centroid_NMDS1, y = centroid_NMDS2, fill = factor(year)),
    shape = 21, size = 5, color = "black", show.legend = FALSE
  ) +
  geom_text(
    data = centroids,
    aes(x = centroid_NMDS1, y = centroid_NMDS2, label = year),
    vjust = -1.2, hjust = -0.2, size = 3, color = "black", show.legend = FALSE
  ) +
  theme_minimal() +
  labs(shape = "Vegetation Type", title = "NMDS Ordination with Year Centroids") +
  guides(color = "none", fill = "none") +
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


#### mean distance to centroid #################################################

# this is for all veg types, but calculations have been done with respect to each subsecstions vegetation type

centroids <- nmds_plot_data %>%
  # 1. Remove unwanted vegetation type
  filter(veg_type != "saltmarsh") %>%               # <-- change if spelling differs
  # 2. Group by the two factors that define each centroid
  group_by(year, veg_type) %>%
  # 3. Compute the mean of each NMDS axis for the group
  summarise(
    cent_NMDS1 = mean(NMDS1, na.rm = TRUE),
    cent_NMDS2 = mean(NMDS2, na.rm = TRUE),
    cent_NMDS3 = mean(NMDS3, na.rm = TRUE),
    .groups = "drop"
  )

# Quick peek at the result
print(centroids)

# 1️⃣  Left‑join so every row keeps its original data
nmds_with_dist <- nmds_plot_data %>%
  left_join(centroids, by = c("year", "veg_type"))

# 2️⃣  Verify the join succeeded (optional sanity check)
if(any(is.na(nmds_with_dist$cent_NMDS1))) {
  warning("Some rows did not find a matching centroid – check year/veg_type combos.")
}

# 3️⃣  Euclidean distance in 3‑dimensional NMDS space
nmds_with_dist <- nmds_with_dist %>%
  mutate(
    dist_to_centroid = sqrt(
      (NMDS1 - cent_NMDS1)^2 +
        (NMDS2 - cent_NMDS2)^2 +
        (NMDS3 - cent_NMDS3)^2
    )
  )

# Quick glimpse of the new column
head(nmds_with_dist[, c("sub_year_vt", "year", "veg_type",
                        "dist_to_centroid")])

nmds_with_dist <- nmds_plot_data %>%
  # Keep every column from nmds_plot_data (the left table)
  left_join(centroids, by = c("year", "veg_type"))

# checking that all wors have a value
summary(nmds_with_dist$dist_to_centroid)

n_unique <- n_distinct(paste0(nmds_with_dist$year, "_", nmds_with_dist$veg_type))
cat("Unique (year × veg_type) combos:", n_unique, "\n")
n_unique

cat("Unique (year × veg_type) combos:", n_unique, "\n")

nmds_with_dist <- nmds_plot_data %>%
  left_join(centroids, by = c("year", "veg_type")) %>%   # keep veg_type
  mutate(
    dist_to_centroid = sqrt(
      (NMDS1 - cent_NMDS1)^2 +
        (NMDS2 - cent_NMDS2)^2 +
        (NMDS3 - cent_NMDS3)^2
    )
  )

# 4️⃣  Summary table: mean ± SD of distance for each (year, veg_type)
dispersion_summary <- nmds_with_dist %>%
  group_by(year, veg_type) %>%                # keep the two‑factor structure
  summarise(
    n          = n(),
    mean_dist  = mean(dist_to_centroid, na.rm = TRUE),
    sd_dist    = sd(dist_to_centroid,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(year, veg_type)

# View the table
print(dispersion_summary)

# Are there any NAs in the distance column?
sum(is.na(nmds_with_dist$dist_to_centroid))

# Does every (year, veg_type) combo have at least one observation?
nmds_with_dist %>%
  count(year, veg_type) %>%
  filter(n == 0)

# baseline model (M1)
model_dist_to_centroid <- lmer(dist_to_centroid ~ year + (1 | subsection),
           data = nmds_with_dist,
           REML = FALSE)   # ML for model comparison
summary(model_dist_to_centroid) 

model_dist_to_centroid2 <- lmer(dist_to_centroid ~ year + (1 | subsection) + (1 | veg_type),
           data = nmds_with_dist,
           REML = FALSE)
summary(model_dist_to_centroid2)

anova(model_dist_to_centroid, model_dist_to_centroid2)

emm <- emmeans(model_dist_to_centroid, ~ year)
pairs(emm, adjust = "tukey")

#### plotting distance to respective centroid #####

emm_year <- emmeans(model_dist_to_centroid, ~ year)

# Generate the CLD (default uses Tukey letters)
cld_tbl <- cld(emm_year, Letters = letters, adjust = "sidak") %>%
  as.data.frame() %>%
  rename(
    year = year,
    est  = emmean,
    se   = SE,
    lower = lower.CL,
    upper = upper.CL,
    cld  = .group   # column that holds the letters
  ) %>%
  mutate(year = factor(year))   # keep factor ordering consistent


## 2.  Get marginal means (EMMs) and their SEs
emm_year <- emmeans(model_dist_to_centroid, ~ year) %>%   # model object name
  as.data.frame() %>%                                   # turn into plain df
  rename(
    est   = emmean,      # estimated mean distance
    se    = SE,          # standard error
    lower = lower.CL,    # lower 95 % CI (optional)
    upper = upper.CL     # upper 95 % CI (optional)
  ) %>%
  mutate(year = factor(year))   # ensure year is treated as categorical

## 3.  Prepare the raw data for plotting
raw_plot_df <- nmds_with_dist |> 
  dplyr::select(sub_year_vt, year, dist_to_centroid) |> 
  mutate(year = factor(year))   # match factor levels used in emm_year


## 4.  Build the ggplot

p <- ggplot() +
  
  ## 4a. Jittered raw observations (semi‑transparent)
  geom_jitter(data = raw_plot_df,
              aes(x = year, y = dist_to_centroid),
              width = 0.15,               # horizontal jitter amount
              alpha = 0.35,               # point transparency
              size  = 1.5,
              colour = "#555555") +      # neutral grey
  
  ## 4b. Model‑estimated means (solid larger points)
  geom_point(data = emm_year,
             aes(x = year, y = est),
             colour = "#D55E00",        # a strong orange for visibility
             size   = 3) +
  
  ## 4c. Error bars (SE)
  geom_errorbar(data = emm_year,
                aes(x = year,
                    ymin = est - se,
                    ymax = est + se),
                width = 0.2,
                colour = "#D55E00",
                size   = 0.9) +
  
  ## 4d. Optional: 95 % CI bars (uncomment if you prefer CIs)
  # geom_errorbar(data = emm_year,
  #               aes(x = year,
  #                   ymin = lower,
  #                   ymax = upper),
  #               width = 0.2,
  #               colour = "#D55E00",
  #               size   = 0.7,
  #               linetype = "dashed") +
  
  ## 4e. Axis labels, theme, etc.
  labs(
    x = "",
    y = "Distance to vegetation-type centroid",
    title = ""
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 0, vjust = 0.7)
  )+
  geom_text(data = cld_tbl,
            aes(x = year,
                y = 1.25,   # small offset above the CI bar
                label = cld),
            colour = "black",
            vjust = 0,           # align bottom of text with the offset
            size = 5) 


## 5.  Print / save the figure
print(p)                     # display in RStudio / notebook
ggsave("fig_distance_vs_year.png", p,
       width = 7, height = 5, dpi = 300)

#### NMDS for each vegetation type #############################################

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

#### testing mean distance EACH vegtype ####################################

## 1️⃣  Global NMDS site scores + metadat

# 1.1  Pull the three NMDS axes for every plot
site_scores <- scores(nmds_result, display = "sites") |>
  as.data.frame() |>
  tibble::rownames_to_column("sub_year_vt")

# 1.2  Join the metadata (veg_type, year, etc.)
site_scores <- site_scores |>
  left_join(meta, by = "sub_year_vt")

# 1.3  Keep only the columns we need
site_scores <- scores(nmds_result, display = "sites") |>
  as.data.frame() |>
  tibble::rownames_to_column("sub_year_vt") |>
  left_join(meta, by = "sub_year_vt") |>
  dplyr::select(
    sub_year_vt,
    veg_type,
    year,
    subsection,          # <<< keep this column
    NMDS1, NMDS2, NMDS3 # keep the three NMDS axes
  )

# 1. Compute global centroids (veg_type × year)
centroids_all <- site_scores |>
  dplyr::group_by(veg_type, year) |>
  dplyr::summarise(
    cen1 = mean(NMDS1),
    cen2 = mean(NMDS2),
    cen3 = mean(NMDS3),
    .groups = "drop"
  )

# 2. Join and create the distance column (single operation)
site_scores <- site_scores |>
  dplyr::left_join(centroids_all,
                   by = c("veg_type", "year")) |>
  dplyr::mutate(
    dist_to_centroid = sqrt(
      (NMDS1 - cen1)^2 +
        (NMDS2 - cen2)^2 +
        (NMDS3 - cen3)^2
    )
  ) |>
  # <‑‑ make sure it is numeric
  dplyr::mutate(dist_to_centroid = as.numeric(dist_to_centroid)) |>
  dplyr::select(-cen1, -cen2, -cen3)
##  Global y‑axis limits (computed once)
# Raw distances from the full data set (used for the jitter points)
global_y_max <- site_scores |>
  dplyr::summarise(max_dist = max(dist_to_centroid, na.rm = TRUE)) |>
  dplyr::pull(max_dist)

# Upper bound for the error‑bar + mean point (adds a small cushion)
global_y_upper <- global_y_max * 1.10   # 10 % headroom for the CLD letters
global_y_limits <- c(0, global_y_upper)   # lower bound = 0 (distances cannot be negative)

## 2️⃣  Loop over vegetation types

veg_vec   <- unique(site_scores$veg_type)
plot_list <- vector("list", length(veg_vec))
mlm_table <- tibble()   # collect fixed‑effect summaries
pairwise_table  <- tibble()   # pairwise year contrasts for all veg types

for (i in seq_along(veg_vec)) {
  vt <- veg_vec[i]
  
  ## --------------------------------------------------------------
  ## 2.1  Subset data (year stays a factor for the model)
  ## --------------------------------------------------------------
  vt_dat <- site_scores |>
    dplyr::filter(veg_type == vt) |>
    dplyr::mutate(year = factor(year))
  
  ## --------------------------------------------------------------
  ## 2.2  Fit mixed‑effects model (year fixed, subsection random)
  ## --------------------------------------------------------------
  mod <- lmerTest::lmer(
    dist_to_centroid ~ year + (1 | subsection),
    data = vt_dat,
    REML = FALSE
  )
  
  ## --------------------------------------------------------------
  ## 2.3  EMMeans (needed for both CLD and pairwise contrasts)
  ## --------------------------------------------------------------
  emm <- emmeans::emmeans(mod, ~ year)
  
  ## --------------------------------------------------------------
  ## 2.4  CLD letters (Sidak‑adjusted) – for the plot
  ## --------------------------------------------------------------
  cld_df <- multcomp::cld(
    emm,
    adjust = "sidak",
    Letters = letters,
    sort = FALSE
  ) |> as_tibble() |> dplyr::select(year, .group)
  
  ## --------------------------------------------------------------
  ## 2.5  Pairwise year contrasts (the numbers behind the CLD)
  ## --------------------------------------------------------------
  pairwise_cmp <- emmeans::contrast(
    emm,
    method = "pairwise",
    adjust = "sidak"
  ) |> as_tibble() |> dplyr::rename(
    contrast = contrast,
    estimate = estimate,
    SE       = SE,
    df       = df,
    t_ratio  = t.ratio,
    p_adj    = p.value
  ) |> dplyr::mutate(veg_type = vt)
  
  ## --------------------------------------------------------------
  ## 2.6  Append pairwise results to the master table
  ## --------------------------------------------------------------
  pairwise_table <- dplyr::bind_rows(pairwise_table, pairwise_cmp)
  
  ## --------------------------------------------------------------
  ## 2.7  Plot‑ready summary (mean ± SE) + attach CLD letters
  ## --------------------------------------------------------------
  plot_stats <- vt_dat |>
    dplyr::group_by(year) |>
    dplyr::summarise(
      mean_dist = mean(dist_to_centroid),
      se_dist   = sd(dist_to_centroid) / sqrt(n()),
      .groups = "drop"
    ) |>
    dplyr::left_join(cld_df, by = "year")
  
  centroid_cld <- plot_stats |>
    dplyr::mutate(
      y_pos = global_y_upper * 0.92   # fixed vertical position inside panel
    ) |>
    dplyr::select(year, .group, y_pos)
  
  ## --------------------------------------------------------------
  ## 2.8  Build the facet plot (all visual settings unchanged)
  ## --------------------------------------------------------------
  p <- ggplot(vt_dat, aes(x = year, y = dist_to_centroid)) +
    geom_jitter(width = 0.15, alpha = 0.4, size = 1.2,
                colour = "black") +
    
    geom_errorbar(
      data = plot_stats,
      aes(y = mean_dist,
          ymin = mean_dist - se_dist,
          ymax = mean_dist + se_dist),
      width = 0.2, colour = "#076834", size = 0.8
    ) +
    
    geom_point(
      data = plot_stats,
      aes(y = mean_dist),
      colour = "#076834", size = 3
    ) +
    
    geom_text(
      data = centroid_cld,
      aes(x = year, y = y_pos, label = .group),
      vjust = 0, size = 4, colour = "#076834"
    ) +
    
    labs(
      title = vt,      # plain external title
      x = "", y = ""
    ) +
    
    coord_cartesian(ylim = global_y_limits, clip = "off") +
    
    theme_minimal() +
    
    theme(
      plot.title = element_text(face = "plain", hjust = 0.5),
      text = element_text(size = 8)
    )
  
  ## --------------------------------------------------------------
  ## 2.9  Store the plot
  ## --------------------------------------------------------------
  plot_list[[i]] <- p
  
  ## --------------------------------------------------------------
  ## 2.10  Extract fixed‑effect summary (includes p‑values)
  ## --------------------------------------------------------------
  fixed_tab <- summary(mod)$coefficients |>
    as_tibble(rownames = "term") |>
    dplyr::rename(
      estimate   = Estimate,
      std.error  = `Std. Error`,
      statistic  = `t value`,
      p.value    = `Pr(>|t|)`
    ) |>
    dplyr::mutate(
      veg_type = vt,
      effect   = "fixed"
    )
  
  ## --------------------------------------------------------------
  ## 2.11  Append to the master MLM table
  ## --------------------------------------------------------------
  mlm_table <- dplyr::bind_rows(mlm_table, fixed_tab)
}

# Show the figure
print(grid_plot)

## 3.1  Combined figure
grid_plot <- wrap_plots(plot_list, ncol = 2) &
  theme(plot.margin = margin(5, 5, 5, 5))


## 3.2  MLM results table (six columns you asked for)

mlm_results <- mlm_table |>
  dplyr::select(
    veg_type,
    term,
    estimate,
    std.error,
    statistic,
    p.value
  ) |>
  dplyr::rename(
    SE      = std.error,
    t_value = statistic,
    p_val   = p.value
  )


## 3.3  Pairwise year‑wise contrast table (supplemental)

pairwise_results <- pairwise_table |>
  dplyr::select(
    veg_type,
    contrast,
    estimate,
    SE,
    df,
    t_ratio,
    p_adj
  ) |>
  dplyr::arrange(veg_type, contrast) |>
  dplyr::mutate(
    estimate = round(estimate, 4),
    SE       = round(SE, 4),
    t_ratio  = round(t_ratio, 2),
    p_adj    = signif(p_adj, 3)   # scientific notation for tiny p‑values
  )

view(pairwise_table)

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
lims_x <- range(nmds_plot_data$NMDS1)
lims_y <- range(nmds_plot_data$NMDS2)

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
    annotate("text", x = -1.5, y = -1, 
             label = paste(y), 
             hjust = 0, vjust = 0,
             size = 4) +
    theme_minimal() +
    coord_cartesian(xlim = lims_x, ylim = lims_y)+
    theme(legend.position = "none")
})

# Combine plots in a grid
wrap_plots(plots, ncol = 2)  

#### NMDS year-vegtype centroid distance ####

# centroid distance construction
centroid_pairs <-
  centroids |>
  inner_join(
    centroids,
    by = "year",
    suffix = c("1", "2")
  ) |>
  filter(as.character(veg_type1) < as.character(veg_type2)) |>
  mutate(
    dist = sqrt(
      (centroid_NMDS11 - centroid_NMDS12)^2 +
        (centroid_NMDS21 - centroid_NMDS22)^2
    ),
    veg_pair = interaction(veg_type1, veg_type2, drop = TRUE)
  ) |>
  dplyr::select(year, veg_type1, veg_type2, veg_pair, dist)

# Keep the exploratory summaries (but reframe them)
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

# mixed model
m_convergence <- lmer(
  dist ~ factor(year) + (1 | veg_pair),
  data = centroid_pairs
)

summary(m_convergence)
anova(m_convergence)

#convergence

emm_year <- emmeans(m_convergence, ~ year)
emm_year
pairs(emm_year)

emm_df <- as.data.frame(emm_year)

emm_df <- as.data.frame(emm_year) %>%
  mutate(year = factor(year, levels = levels(emm_year@grid$year)),
         x_pos = as.integer(year))


#cld
centroid_cld <- cld(emm_year,
                    adjust = "sidak",      # same adjustment you used for pairs()
                    Letters = letters) %>%
  as.data.frame() %>%
  # keep the same ordering as the plot
  mutate(year = factor(year, levels = levels(emm_year@grid$year)),
         x_pos = as.integer(year))

#ggplot
y_max <- max(emm_df$upper.CL, na.rm = TRUE)
y_max <- y_max * 1.05

ggplot() +
  
  ## 4.1  Raw pairwise distances (jittered points)
  geom_jitter(data = centroid_pairs,
              aes(x = as.integer(factor(year)), y = dist),
              width = 0.15, height = 0,
              alpha = 0.4, colour = "gray30") +
  
  ## 4.2  Model‑based means with 95 % CIs
  geom_errorbar(data = emm_df,
                aes(x = x_pos, ymin = lower.CL, ymax = upper.CL),
                width = 0.2, colour = "steelblue") +
  
  geom_point(data = emm_df,
             aes(x = x_pos, y = emmean),
             size = 2.8, colour = "steelblue") +
  
  ## 4.3  Significance letters (CLD)
  geom_text(data = centroid_cld,
            aes(x = x_pos, label = .group),
            # place just above the highest CI bar
            y = 1.45,
            vjust = 0, size = 4, colour = "steelblue") +
  
  ## 4.4  Axis formatting (four years)
  scale_x_continuous(name = "Year",
                     breaks = 1:4,
                     labels = c("2007","2012","2017","2022")) +
  
  labs(y = "Estimated mean centroid distance") +
  
  theme_minimal()

#Diagnostics
plot(m_convergence)
qqnorm(residuals(m_convergence))
qqline(residuals(m_convergence))

