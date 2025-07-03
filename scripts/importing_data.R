# data prepping

# Load necessary libraries
library(dplyr)
library(readr)
library(readxl)

# Define the path to the subfolder containing your CSV files
data_folder <- "data"

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

# read the functional types files #
eco_veg_growth_forms <- read_excel("data/eco_veg_growth_forms.xlsx")

# Read and merge all CSV files into one data frame
merged_data <- file_list |> 
  lapply(read_and_standardize) |>   # Apply standardization function to each file
  bind_rows() |>                    # Combine them into one data frame
  filter(raunkiaer_value != -9999) |>  # Remove invalid values
  mutate(presence = ifelse(raunkiaer_value > 0, 1, 0)) |> # Convert to presence-absence
  left_join(eco_veg_growth_forms, by = "taxon_code") |> 
  select(!validation,)

inspect <- merged_data |> 
  filter(year == 2022, 
         plot_id == 3.02 )

# View the merged data
print(merged_data)
merged_data$veg_type <- as.factor(merged_data$veg_type)
merged_data$taxon_code <- as.factor(merged_data$taxon_code)
merged_data$species <- as.factor(merged_data$species)
merged_data$func_type <- as.factor(merged_data$func_type)
merged_data$ecoveg_gfc <- as.factor(merged_data$ecoveg_gfc)
merged_data$ecoveg_sgfc <- as.factor(merged_data$ecoveg_sgfc)
summary(merged_data)

write_rds(merged_data, "data/merged_data.rds")

this.path::this.path()
