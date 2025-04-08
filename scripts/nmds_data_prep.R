library(dplyr)
library(tidyr)

# Summarize duplicate rows before pivoting
nmds_data <- merged_data %>%
  select(year, plot_id, veg_type, taxon_code, raunkiaer_value) %>%
  filter(!is.na(raunkiaer_value)) %>%  # Remove rows with missing values
  group_by(year, plot_id, veg_type, taxon_code) %>%
  summarize(raunkiaer_value = sum(raunkiaer_value), .groups = "drop")  # Summarize duplicates

# Pivot the data to create a wide matrix
nmds_matrix <- nmds_data %>%
  pivot_wider(
    names_from = taxon_code,
    values_from = raunkiaer_value,
    values_fill = 0  # Fill missing values with 0
  ) %>%
  ungroup()  # Remove grouping

# Convert to a numeric matrix for NMDS
nmds_matrix <- as.data.frame(nmds_matrix)
rownames(nmds_matrix) <- nmds_matrix$plot_id
nmds_matrix <- nmds_matrix %>% select(-plot_id, -veg_type, -year)  # Exclude non-numeric columns
nmds_matrix <- as.matrix(nmds_matrix)
