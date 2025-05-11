library(tidyverse)

# Summarize duplicate rows before pivoting ###################################################################################
names(merged_data)

nmds_data <- merged_data |>
  filter(presence == 1) |> 
  select(year, plot_id, veg_type, taxon_code, raunkiaer_value, presence) |> 
  filter(!is.na(raunkiaer_value)) |>   # Remove rows with missing values
  group_by(year, plot_id, veg_type, taxon_code, presence) |> 
  summarize(raunkiaer_value = sum(raunkiaer_value), .groups = "drop")  # Summarize duplicates

# see if any plots dont have data
check_empty_plots <- nmds_data |> 
  group_by(year, plot_id) |> 
  summarize(sum = sum(presence)) |> 
  filter(sum == 0)


names(nmds_data)

# Pivot the data to create a wide matrix
nmds_matrix <- nmds_data |> 
  pivot_wider(
    names_from = taxon_code,
    values_from = presence,
    values_fill = 0  # Fill missing values with 0
  ) |> 
  ungroup() |>   # Remove grouping
  select(-raunkiaer_value)
  
empty <- nmds_matrix |> 
  filter(rowSums(across(4:ncol(.))) == 0)

# Convert to a numeric matrix for NMDS
nmds_matrix <- as.data.frame(nmds_matrix)
rownames(nmds_matrix) <- nmds_matrix$plot_id
nmds_matrix <- nmds_matrix |>  select(-plot_id, -veg_type, -year)  # Exclude non-numeric columns


#### checking for empty rows ####
empty <- nmds_matrix |> 
  filter(rowSums(across(4:ncol(.))) == 0)

#### substting for heath ############################################################################################################

# Filter the data for a single vegetation type
heath_nmds <- merged_data |> 
  filter(veg_type == "heath")

print(head(heath_nmds))

# Summarize duplicate rows before pivoting
heath_nmds_summary <- heath_nmds |> 
  group_by(year, plot_id, veg_type, taxon_code) |> 
  summarize(presence = sum(presence), .groups = "drop")

# Pivot the data to create a wide matrix
heath_nmds_matrix <- heath_nmds_summary |> 
  pivot_wider(
    names_from = taxon_code,
    values_from = presence,
    values_fill = 0  # Fill missing values with 0
  )

# Regenerate metadata (plot_id and year) from heath_nmds
heath_metadata <- heath_nmds |> 
  select(plot_id, year) |> 
  distinct()  # Ensure unique rows for each plot_id and year

# View the regenerated metadata
print(head(heath_metadata))

# View the resulting matrix
print(heath_nmds_matrix)

# Ensure only numeric columns are included in the NMDS matrix
heath_nmds_numeric <- heath_nmds_matrix |> 
  select(-plot_id, -veg_type, -year)  # Remove non-numeric columns

# Convert to a numeric matrix
heath_nmds_numeric <- as.matrix(heath_nmds_numeric)

# Set row names to plot IDs for identification in NMDS
rownames(heath_nmds_numeric) <- heath_nmds_matrix$plot_id

# Check for rows with all zero values
empty_rows <- rowSums(heath_nmds_numeric) == 0
print(paste("Number of empty rows:", sum(empty_rows)))

# Check for negative values
negative_values <- heath_nmds_numeric < 0
print(paste("Number of negative entries:", sum(negative_values)))

# Remove rows with all zeros
heath_nmds_clean <- heath_nmds_numeric[!empty_rows, ]

# Ensure there are no negative values (replace negatives with 0, if appropriate)
heath_nmds_clean[heath_nmds_clean < 0] <- 0

# Check for rows with all zeros
empty_rows <- rowSums(heath_nmds_clean) == 0
print(paste("Number of empty rows:", sum(empty_rows)))

# Remove rows with all zeros
heath_nmds_clean <- heath_nmds_clean[!empty_rows, ]

# Verify the cleaned matrix
print(dim(heath_nmds_clean))

# Check for missing values
missing_values <- apply(heath_nmds_clean, 1, function(row) any(is.na(row)))
print(paste("Number of rows with missing values:", sum(missing_values)))

# Remove rows with missing values
heath_nmds_clean <- heath_nmds_clean[!missing_values, ]

# Verify the cleaned matrix again
print(dim(heath_nmds_clean))


#### substting for copse ##################################################################################################################################

# Filter the data for a single vegetation type
copse_nmds <- merged_data |> 
  filter(veg_type == "copse")

print(head(copse_nmds))

# Summarize duplicate rows before pivoting
copse_nmds_summary <- copse_nmds |> 
  group_by(year, plot_id, veg_type, taxon_code) |> 
  summarize(presence = sum(presence), .groups = "drop")

# Pivot the data to create a wide matrix
copse_nmds_matrix <- copse_nmds_summary |> 
  pivot_wider(
    names_from = taxon_code,
    values_from = presence,
    values_fill = 0  # Fill missing values with 0
  )

# Regenerate metadata (plot_id and year) from heath_nmds
copse_metadata <- copse_nmds |> 
  select(plot_id, year) |> 
  distinct()  # Ensure unique rows for each plot_id and year

# View the regenerated metadata
print(head(copse_metadata))

# View the resulting matrix
print(copse_nmds_matrix)

# Ensure only numeric columns are included in the NMDS matrix
copse_nmds_numeric <- copse_nmds_matrix |> 
  select(-plot_id, -veg_type, -year)  # Remove non-numeric columns

# Convert to a numeric matrix
copse_nmds_numeric <- as.matrix(copse_nmds_numeric)

# Set row names to plot IDs for identification in NMDS
rownames(copse_nmds_numeric) <- copse_nmds_matrix$plot_id

# Check for rows with all zero values
empty_rows <- rowSums(copse_nmds_numeric) == 0
print(paste("Number of empty rows:", sum(empty_rows)))

# Check for negative values
negative_values <- copse_nmds_numeric < 0
print(paste("Number of negative entries:", sum(negative_values)))

# Remove rows with all zeros
copse_nmds_clean <- heath_nmds_numeric[!empty_rows, ]

# Ensure there are no negative values (replace negatives with 0, if appropriate)
copse_nmds_clean[heath_nmds_clean < 0] <- 0

# Check for rows with all zeros
empty_rows <- rowSums(copse_nmds_clean) == 0
print(paste("Number of empty rows:", sum(empty_rows)))

# Remove rows with all zeros
copse_nmds_clean <- heath_nmds_clean[!empty_rows, ]

# Verify the cleaned matrix
print(dim(copse_nmds_clean))

# Check for missing values
missing_values <- apply(copse_nmds_clean, 1, function(row) any(is.na(row)))
print(paste("Number of rows with missing values:", sum(missing_values)))

# Remove rows with missing values
copse_nmds_clean <- copse_nmds_clean[!missing_values, ]

# Verify the cleaned matrix again
print(dim(copse_nmds_clean))
