# data prepping
# Load necessary libraries ###########################################################################################
library(dplyr)
library(readr)
library(readxl)

#### Define the path to the subfolder containing your CSV files ###########################################################################
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

#### read the functional types files ##############################################
eco_veg_growth_forms <- read_excel("data/eco_veg_growth_forms.xlsx")
#View(eco_veg_growth_forms)

# Read and merge all CSV files into one data frame
data_raw <- file_list |> 
  lapply(read_and_standardize) |>   # Apply standardization function to each file
  bind_rows() |>                    # Combine them into one data frame
  filter(raunkiaer_value != -9999) |>  # Remove invalid values
  mutate(presence = ifelse(raunkiaer_value > 0, 1, 0)) |> # Convert to presence-absence
  dplyr::select(-validation,-remarks,-cfr,-fertile,-raunkiaer_value) |> 
  dplyr::distinct()

length(unique(data_raw$plot_id))
length(unique(data_raw$year))
length(unique(data_raw$year))*length(unique(data_raw$plot_id))
length(unique(data_raw$taxon_code))
length(unique(data_raw$vt_section))

n_years  <- dplyr::n_distinct(data_raw$year)
n_plots  <- dplyr::n_distinct(data_raw$plot_id)
n_taxa   <- dplyr::n_distinct(data_raw$taxon_code)

n_years * n_plots * n_taxa

data_wide <- data_raw |> 
  pivot_wider(names_from = "taxon_code", values_from = "presence", values_fill = 0, values_fn = max) |> 
  dplyr::distinct()

data_long <- data_wide |> 
  pivot_longer(cols = carbig:viospe, names_to = "taxon_code", values_to = "presence") |> 
  dplyr::distinct() |> 
  left_join(eco_veg_growth_forms, by = "taxon_code")

inspect <- merged_data |> 
  filter(year == 2022, 
         plot_id == 3.02 )

# View the merged data
print(data_long)
data_long$veg_type <- as.factor(data_long$veg_type)
data_long$taxon_code <- as.factor(data_long$taxon_code)
data_long$species <- as.factor(data_long$species)
data_long$func_type <- as.factor(data_long$func_type)
data_long$ecoveg_gfc <- as.factor(data_long$ecoveg_gfc)
data_long$ecoveg_sgfc <- as.factor(data_long$ecoveg_sgfc)
summary(data_long)

data_long<- data_long |> 
#  filter(!is.na(func_type)) |>  
  mutate(subsection = case_when(
    plot >= 1  & plot <= 5  ~ paste(vt_section, "1", sep = "."),
    plot >= 6  & plot <= 10 ~ paste(vt_section, "2", sep = "."),
    plot >= 11 & plot <= 15 ~ paste(vt_section, "3", sep = "."),
    plot >= 16 & plot <= 20 ~ paste(vt_section, "4", sep = "."),
    TRUE ~ NA_character_
  ))

summary(data_long)
subsection <- data_long |> 
  dplyr::select(year,vt_section,veg_type,taxon_code,species,presence,func_type,ecoveg_gfc,ecoveg_sgfc)

data_long$subsection <- as.factor(data_long$subsection)

write_rds(data_long, "data/data_long.rds")
View(data_long)

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
