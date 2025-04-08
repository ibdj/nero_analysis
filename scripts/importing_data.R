# data prepping

# Load necessary libraries
library(dplyr)
library(readr)


# Load necessary libraries
library(dplyr)
library(readr)

# Define the path to the subfolder containing your CSV files
data_folder <- "data"

# List all files in the subfolder matching the pattern 'data_[YEAR]_data.csv'
#file_list <- list.files(path = data_folder, pattern = "data_2022_data.csv", full.names = TRUE)
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
             fertile = as.numeric(fertile) ) 
  }
  
  return(df)
}

# Read and merge all CSV files into one data frame
merged_data <- file_list %>%
  lapply(read_and_standardize) %>%  # Apply standardization function to each file
  bind_rows()                       # Combine them into one data frame

# View the merged data
print(merged_data)

test_file <- "data/data_2022_data.csv"
test_data <- read_delim(test_file, delim = ",")  # Adjust `delim` if needed
print(head(test_data))  # View the first few rows
print(colnames(test_data))  # Check column names
