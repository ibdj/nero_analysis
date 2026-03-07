#### arthropod data ####

library(readxl)
library(tidyverse)

#### import data ####

Biobasis_Nuuk_Arthropod <- read_excel("~/Library/CloudStorage/OneDrive-GrønlandsNaturinstitut/General - BioBasis/03_GEM_Database/Datafiler excel/Biobasis_Nuuk_Arthropod.xlsx")
names(Biobasis_Nuuk_Arthropod)

# > names(art)
# [1] "Date"      "Hour"      "Trap"      "Plot"      "Sorter"    "Snow_A"    "Snow_B"    "Snow_C"    "Snow_D"   
# [10] "Snow_W"    "Days_A"    "Days_B"    "Days_C"    "Days_D"    "Days_W"    "Phylum"    "Order"     "Family"   
# [19] "Genus"     "Species"   "Species_A" "Species_B" "Species_C" "Species_D" "Species_W" "Latitude"  "Longitude"
# [28] "Remarks" 

base_taxonomy <- Biobasis_Nuuk_Arthropod |> 
  group_by(Phylum, Order, Family, Genus, Species) |> 
  summarise(rows = n())

base_taxonomy_year <- Biobasis_Nuuk_Arthropod |> 
  #mutate(Order = str_replace_all(Order, "^(Aranea|Aranae)$", "Araneae")) |> 
  mutate(Year = year(Date)) |> 
  group_by(Phylum, Order, Family, Genus, Species, Year) |> 
  summarise(rows = n()) |> 
  pivot_wider(values_from = rows, names_from = Year)


base_taxonomy_year_edit <- base_taxonomy_year |> 
mutate(Order = str_replace_all(Order, "^(Aranea|Aranae)$", "Araneae")) |> 

  
#Aranae #Aranea > #Araneae
#Diapridae > #Diapriidae Diapridae
#Mycetophiliidae
#Nematoda > Phyllum or order?

#### basic split #####

base_split1_art <- Biobasis_Nuuk_Arthropod |> 
  filter(Trap == "Art") |> 
  select(Date, Hour, Plot, Sorter, Snow_A, Snow_B, Snow_C, Snow_D, Days_A, Days_B, Days_C, Days_D, Phylum, Order, Family, Genus, Species, Count_A, Count_B, Count_C, Count_D, Latitude, Longitude)

base_split1_wart <- Biobasis_Nuuk_Arthropod |> 
  filter(Trap == "WArt") |> 
  select(Date, Hour, Plot, Sorter, Snow_W, Days_W, Phylum, Order, Family, Genus, Species, Count_W, Latitude, Longitude)



#### pivot all ####
rearrange1 <- Biobasis_Nuuk_Arthropod |>
  rename_with(~ gsub("^Species_", "Count_", .x), starts_with("Species_")) |>
  pivot_longer(
    cols = matches("^(Snow|Days|Count)_[A-Z]$"),
    names_to = c(".value", "trap_letter"),
    names_sep = "_"
  )

rearrange2 <- rearrange1 |> 
  mutate(
    Trap_id = if_else(
      Trap == "WArt",
      paste0(Trap, Plot),
      paste0(Trap, Plot, trap_letter)
    )
  ) |> 
  filter(Count > -1) |> 
  mutate(Year = year(Date)) |> 
  select(Year, Date, Hour, Trap, Plot, trap_letter, Trap_id, Sorter, Snow, Days, Phylum, Order, Family, Genus, Species, Count, Latitude, Longitude, Remarks)
  
rearrange2$Trap_id <- as.factor(rearrange2$Trap_id)

summary(rearrange2)

unique(rearrange2$Trap_id)

rearrange3_filtered <- rearrange2 |> 
  filter(Count > -1) |> 
  mutate(Order = str_replace_all(Order, "^(Aranea|Aranae)$", "Araneae"),
         Week = week(Date), 
         sample_week = floor_date(Date, unit = "week"))

#Aranae #Aranea > #Araneae
#Diapridae > #Diapriidae
#Mycetophiliidae
#Nematoda > Phyllum or order?

stat <- rearrange3_filtered |>
  distinct(Year, Phylum, Order, Family, Genus, Species) |>
  count(Year, name = "n_unique_taxa")

taxonomy <- rearrange3_filtered |> 
  group_by(Phylum, Order, Family, Genus, Species, Year) |> 
  reframe(count = n())

taxonomy_pivot <- taxonomy |>
  pivot_wider(
    names_from = Year,
    values_from = count,
    values_fill = NA
  )

#### pivoting wider ####

rearrange4 <- rearrange3_filtered |>
  complete(
    nesting(sample_week, Trap_id),
    nesting(Phylum, Order, Family, Genus, Species),
    fill = list(Count = 0)
  ) |>
  group_by(sample_week, Trap_id) |>
  fill(Snow, Date, Days, Year, Trap, Plot, Hour, Sorter,
       trap_letter, Latitude, Longitude, Remarks,
       .direction = "downup") |>
  ungroup()

unique(rearrange4$Trap)


taxonomy4 <- rearrange4 |> 
  group_by(Phylum, Order, Family, Genus, Species) |> 
  summarise(count = n())

#### splitting ####

split_wart <- rearrange4 |> 
  filter(Trap == "WArt")

split_0_wart <- split_wart |> 
  group_by(Phylum, Order, Family, Genus, Species) |> 
  summarise(sum = sum(Count))

split_art <- rearrange4 |> 
  filter(Trap == "Art")

split_0_art <- split_art |> 
  group_by(Phylum, Order, Family, Genus, Species) |> 
  summarise(sum = sum(Count))

stat_split <- rearrange4 |> 
  group_by(Trap) |> 
  summarise(count = n())

dplyr::anti_join(split_0_art |> select(!sum), split_0_wart |> select(!sum)) #: Returns all rows from df1 that do not have a match in df2.
dplyr::intersect(df1, df2): Returns rows that are present in both data frames.

rearrange4 |>
  group_by(sample_week, Trap_id) |>
  summarise(all_na = all(is.na(Trap)), .groups = "drop") |>
  filter(all_na)

rearrange4 |>
  count(sample_week, Trap_id, ) |>
  summarise(min_rows = min(n), max_rows = max(n))
