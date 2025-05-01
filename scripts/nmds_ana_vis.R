# nmds analysis

library(vegan)

#List of number of plots in each veg type
# heath 8076
# copse 1875
# snowbed 1765
# fen 1738
# abrasion 284
# salt marsh 95

#### heath ######################################################################################################################################################################
# Perform NMDS using jaccard dissimilarity (because now presence-absence data)
head(heath_nmds_clean)
heath_nmds_result <- metaMDS(heath_nmds_clean, distance = "jaccard", k = 3, trymax = 100)

# Check the stress value
print(paste("Stress value:", heath_nmds_result$stress))


# Extract NMDS site scores (coordinates for plots)
nmds_scores <- as.data.frame(scores(heath_nmds_result)$sites)


# Add plot IDs as a column
nmds_scores$plot_id <- rownames(nmds_scores)

# View the extracted scores
print(head(nmds_scores))

# Remove prefix "X" from plot_id in nmds_scores
nmds_scores$plot_id <- gsub("^X", "", nmds_scores$plot_id)
print(head(nmds_scores))

# Convert plot_id in both data frames to character (if not already done)
nmds_scores$plot_id <- as.numeric(nmds_scores$plot_id)
heath_metadata$plot_id <- as.numeric(heath_metadata$plot_id)

# Perform the join again
nmds_plot_data <- nmds_scores %>%
  left_join(heath_metadata, by = "plot_id") |> 
  filter(is.na(year))

nmds_plot_data <- merge(nmds_scores, heath_nmds[, c("plot_id", "year")], by = "plot_id")


head(nmds_plot_data)

# View the merged data
print(head(nmds_plot_data))


#### copse #################
# Perform NMDS using Jaccard dissimilarity (because now presence-absence data)
head(copse_nmds_clean)
copse_nmds_result <- metaMDS(copse_nmds_clean, distance = "jaccard", k = 3, trymax = 100)

# Check the stress value
print(paste("Stress value:", copse_nmds_result$stress))


# Extract NMDS site scores (coordinates for plots)
nmds_scores <- as.data.frame(scores(copse_nmds_result)$sites)


# Add plot IDs as a column
nmds_scores$plot_id <- rownames(nmds_scores)

# View the extracted scores
print(head(nmds_scores))

# Remove prefix "X" from plot_id in nmds_scores
nmds_scores$plot_id <- gsub("^X", "", nmds_scores$plot_id)
print(head(nmds_scores))

# Convert plot_id in both data frames to character (if not already done)
nmds_scores$plot_id <- as.numeric(nmds_scores$plot_id)
metadata$plot_id <- as.numeric(metadata$plot_id)

head(metadata)

# Perform the join again

metadata <- merged_data  |> 
  select(plot_id, year) |> 
  distinct()  # Ensure unique rows for each plot_id and year
head(metadata)

nmds_scores <- nmds_scores |>
  mutate(plot_id = as.numeric(plot_id))

head(nmds_scores)
head(metadata)

nmds_plot_data <- nmds_scores %>%
  left_join(metadata, by = "plot_id") |> 
  drop_na()


# View the merged data
print(head(nmds_plot_data))

#### plotting nmds same for all ######################################################

library(ggplot2)

# Plot NMDS results
ggplot(nmds_plot_data, aes(x = NMDS1, y = NMDS2, color = as.factor(year))) +
  geom_point(size = 4) +  # Add points
  geom_text(aes(label = plot_id), vjust = -1.5, size = 3) +  # Add plot labels
  labs(
    title = "NMDS of copse",
    x = "NMDS1",
    y = "NMDS2",
    color = "Year"
  ) +
  geom_line(aes(group = plot_id), color = "grey50", linetype = "dashed") +  # Add connecting lines
  stat_ellipse(aes(group = as.factor(year), fill = as.factor(year)), 
               alpha = 0.3, type = "norm") +  # Add ellipses
  theme_minimal()

hulls <- nmds_plot_data |>
  group_by(year) |>
  slice(chull(NMDS1, NMDS2)) |>
  ungroup()



# Plot NMDS with convex hulls
ggplot(nmds_plot_data, aes(x = NMDS1, y = NMDS2, color = as.factor(year))) +
  geom_point(size = 3) +  # Plot points
  geom_text(aes(label = plot_id), vjust = -1.5, size = 3) +  # Add plot labels
  geom_polygon(data = hulls, aes(group = year, fill = as.factor(year)), alpha = 0.3) +  # Add convex hulls
  labs(
    title = "NMDS with Convex Hulls by Year",
    x = "NMDS Dimension 1",
    y = "NMDS Dimension 2",
    color = "Year",
    fill = "Year"
  ) +
  theme_minimal()

