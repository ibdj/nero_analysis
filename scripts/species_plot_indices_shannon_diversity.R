### packages      ################################################################################################################

library(tidyverse)
library(vegan)
library(ggplot)

#### loading data ################################################################################################################

# Load necessary libraries
library(dplyr)
library(readr)

# Define the path to the subfolder containing your CSV files
data_folder <- "~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/01 Vegetation changes Kobbefjord/data/nmds_nero/nmds_nero/data"

# List all files in the subfolder matching the pattern 'data_[YEAR]_data.csv'
file_list <- list.files(path = data_folder, pattern = "data_\\d{4}_data\\.csv$", full.names = TRUE)

# Function to read and standardize column types
read_and_standardize <- function(file) {
  # Read the file with explicit delimiter (adjust if needed)
  df <- read_csv(file)
  
  # Standardize the "cfr" column (convert to numeric)
  if ("cfr" %in% colnames(df)) {
    df <- df %>%
      mutate(cfr = as.character(cfr),
             raunkiaer_value = as.numeric(raunkiaer_value),
             fertile = as.numeric(fertile),
             veg_type = as.factor(veg_type)) 
  }
  
  return(df)
}

# Read and merge all CSV files into one data frame
merged_data <- file_list |> 
  lapply(read_and_standardize) |>   # Apply standardization function to each file
  bind_rows() |>                    # Combine them into one data frame
  filter(raunkiaer_value != -9999) |>  # Remove invalid values
  mutate(presence = as.numeric(ifelse(raunkiaer_value > 0, 1, 0)))  # Convert to presence-absence

clean <- merged_data |>
  filter(presence == 1) |> 
  select(year, plot_id, veg_type, taxon_code, raunkiaer_value, presence) |> 
  filter(!is.na(raunkiaer_value)) |>   # Remove rows with missing values
  group_by(year, plot_id, veg_type, taxon_code, presence) |> 
  summarize(raunkiaer_value = sum(raunkiaer_value), .groups = "drop") |> 
  select(-raunkiaer_value)# Summarize duplicates
  
matrix <- clean |> 
  pivot_wider(
    names_from = taxon_code,
    values_from = presence,
    values_fill = 0  # Fill missing values with 0
  ) |> 
  ungroup()   # Remove grouping
  
head(matrix)
### shannon index   ################################################################################################################

library(vegan)

# Identify species columns (assuming first three are metadata)
species_cols <- setdiff(names(matrix), c("year", "plot_id", "veg_type"))

# Calculate Shannon index
matrix$shannon_index <- diversity(matrix[, species_cols], index = "shannon")

ggplot(matrix, aes(x = year, y = shannon_index))+
  geom_point()+
  geom_smooth(method = lm)

model <- lm(shannon_index ~ year, data = matrix)
p_value <- summary(model)$coefficients[2, 4]
if (p_value < 0.05) {
  print(paste("Significant trend (p =", round(p_value, 4), "): Shannon diversity changes over time."))
} else {
  print(paste("No significant trend (p =", round(p_value, 4), "): No evidence of change."))
}

### shannon index for each vegetation type (all data)##############################################################################################

library(vegan)
library(ggplot2)
library(dplyr)
library(patchwork)


# 1. Split data by vegetation type and create models (with completeness check)
veg_results <- matrix %>% 
  group_by(veg_type) %>% 
  group_modify(~ {
    if(nrow(.x) >= 2) {  # Need at least 2 data points for regression
      model <- lm(shannon_index ~ year, data = .x)
      p_value <- summary(model)$coefficients[2, 4]
      data.frame(p_value = p_value)
    } else {
      data.frame(p_value = NA_real_)
    }
  }) %>% 
  complete(veg_type)  # Ensure all veg_types are represented

# 2. Create plots with conditional p-value display
plots <- lapply(unique(matrix$veg_type), function(veg) {
  subset_data <- filter(matrix, veg_type == veg)
  veg_p <- veg_results$p_value[veg_results$veg_type == veg]
  # Set linetype based on significance
  ltype <- ifelse(!is.na(veg_p) && veg_p < 0.05, "solid", "dotted")
  
  ggplot(subset_data, aes(x = year, y = shannon_index)) +
    geom_point(color = "darkgreen") +
    geom_smooth(color = "black",
      linewidth = 0.5,
      method = "lm",
      data = subset_data %>% filter(shannon_index >= 0.5),
      linetype = ltype
    ) +
    scale_x_continuous(breaks = c(2007, 2012, 2017, 2022)) +
    labs(
      title = paste(veg),
      subtitle = ifelse(is.na(veg_p), "Insufficient data",
                        paste("p-value:", round(veg_p, 4)))
    ) +
    theme(plot.title = element_text(size = 12), plot.subtitle = element_text(size = 8))
})


# 3. Print all plots in a 2x3 grid (REPLACE THIS SECTION)
wrap_plots(plots, ncol = 3)  # Shows all plots in one window

# 4. Display significance summary (with completeness)
veg_results <- veg_results %>%
  mutate(significant = case_when(
    is.na(p_value) ~ "Insufficient data",
    p_value < 0.05 ~ "Yes",
    TRUE ~ "No"
  ))

print(veg_results)

### shannon index cheching with out low values #########
# 1. Create filtered dataset FIRST
filtered_matrix <- matrix %>% filter(shannon_index >= 0.5)

# 2. Split data by vegetation type and create models
veg_results <- filtered_matrix %>% 
  group_by(veg_type) %>%
  group_modify(~ {
    if(nrow(.x) >= 2) {  # Minimum 2 points for regression
      model <- lm(shannon_index ~ year, data = .x)
      p_value <- summary(model)$coefficients[2, 4]
      data.frame(p_value = p_value)
    } else {
      data.frame(p_value = NA_real_)
    }
  }) %>% 
  complete(veg_type)  # Keep all veg_types even if filtered out

# 3. Create plots using FILTERED data
plots <- lapply(unique(matrix$veg_type), function(veg) {
  # Use filtered data for plotting
  subset_filtered <- filter(filtered_matrix, veg_type == veg)
  
  ggplot() +
    # Plot filtered points
    geom_point(data = subset_filtered, aes(x = year, y = shannon_index)) +
    # Plot trend line only if filtered data exists
    {if(nrow(subset_filtered) >= 2) geom_smooth(
      data = subset_filtered,
      aes(x = year, y = shannon_index),
      method = lm
    )} +
    # Set y-axis limits to exclude values <0.5
    coord_cartesian(ylim = c(0.5, max(matrix$shannon_index))) +
    labs(
      title = paste("Vegetation type:", veg),
      subtitle = ifelse(
        is.na(veg_results$p_value[veg_results$veg_type == veg]),
        "Insufficient data after filtering",
        paste("p-value:", round(veg_results$p_value[veg_results$veg_type == veg], 4))
      )
    ) +
    theme_minimal()
})

# 4. Print all plots
invisible(lapply(plots, print))

# 5. Display significance summary
veg_results <- veg_results %>%
  mutate(significant = case_when(
    is.na(p_value) ~ "Insufficient data",
    p_value < 0.05 ~ "Yes",
    TRUE ~ "No"
  ))

print(veg_results)
