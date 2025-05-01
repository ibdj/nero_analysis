# Load necessary libraries
# Load necessary libraries
library(tidyverse) # For reshaping data


# Filter data for a specific year (e.g., 2007)  
data_2007 <- merged_data  |>  filter(year == 2007)
data_2012 <- merged_data  |>  filter(year == 2012)

# Transform long format to wide format
wide_data <- data_2007 |> 
  select(plot_id, taxon_code, presence) |>  # Select relevant columns
  pivot_wider(names_from = taxon_code, values_from = presence, values_fill = 0) # Fill missing values with 0

# Inspect the wide-format data
head(wide_data)

# Perform DIANA clustering
diana_result <- diana(clustering_data)

# Plot the dendrogram to visualize the clustering structure
plot(diana_result, main = "DIANA Clustering Dendrogram")

# Cut the dendrogram into k clusters (e.g., k = 2 or k = 3 based on dendrogram inspection)
clusters <- cutree(diana_result, k = 2) # Example with k = 2

# Add cluster labels to the original dataset
data_2007$cluster <- clusters

# Inspect the clustered data
head(data_2007)
