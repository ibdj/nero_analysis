# data prepping

#######################################################################################################################################
# Load necessary libraries
library(dplyr)
library(readr)
library(readxl)
#######################################################################################################################################
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
#######################################################################################################################################
# read the functional types files #
eco_veg_growth_forms <- read_excel("data/eco_veg_growth_forms.xlsx")
View(eco_veg_growth_forms)

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

merged_data<- merged_data |> 
  filter(!is.na(func_type)) |>  
  mutate(subsection = case_when(
    plot >= 1  & plot <= 5  ~ paste(vt_section, "1", sep = "."),
    plot >= 6  & plot <= 10 ~ paste(vt_section, "2", sep = "."),
    plot >= 11 & plot <= 15 ~ paste(vt_section, "3", sep = "."),
    plot >= 16 & plot <= 20 ~ paste(vt_section, "4", sep = "."),
    TRUE ~ NA_character_
  ))
merged_data$subsection <- as.factor(merged_data$subsection)

write_rds(merged_data, "data/merged_data.rds")
View(merged_data)

this.path::this.path()

merged_data[is.na(merged_data$ecoveg_sgfc),]

#### nero subsection summmary ####

nero_subsection_summary1 <- merged_data |>
  mutate(subsection = case_when(
    plot >= 1  & plot <= 5  ~ paste(vt_section, "1", sep = "."),
    plot >= 6  & plot <= 10 ~ paste(vt_section, "2", sep = "."),
    plot >= 11 & plot <= 15 ~ paste(vt_section, "3", sep = "."),
    plot >= 16 & plot <= 20 ~ paste(vt_section, "4", sep = "."),
    TRUE ~ NA_character_
  ))

summary(nero_subsection_summary1)  

#### subsection functional type ####

nero_subsection_pft <- nero_subsection_summary1 |>
  group_by(year, subsection, veg_type, func_type) |>
  summarise(
    n_plots = n_distinct(plot_id),     # unique plots in group
    count = n(),                       # row count (species/records) per group
    pft_pr_plot = count/n_plots,       # per-plot value, as before
    .groups = 'drop'                   # ensures a clean, ungrouped result
  )
  
write_rds(nero_subsection_pft, "data/nero_subsection_pft.rds")

#### subsection species #####

# Step 1: Calculate n_plots by year, subsection and veg_type (no species)
plots_summary <- nero_subsection_summary1 |> 
  group_by(year, subsection, veg_type) |> 
  summarise(n_plots = n_distinct(plot_id), .groups = 'drop')

# Step 2: Summarise by species (with other vars)
nero_subsection_species <- nero_subsection_summary1 |> 
  group_by(year, subsection, veg_type, species) |> 
  summarise(
    count = n(),     # records per group (e.g. per species)
    .groups = 'drop'
  ) |> 
  # Step 3: Join n_plots onto this summary
  left_join(plots_summary, by = c("year", "subsection", "veg_type")) |> 
  mutate(fraction = count / n_plots)

summary(nero_subsection_species)
write_rds(nero_subsection_species, "data/nero_subsection_species.rds")
