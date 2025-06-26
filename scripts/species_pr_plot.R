library(ggplot2)
library(ggpmisc)

#### diversity ####

merged_data <- readRDS("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/01 Vegetation changes Kobbefjord/data/nmds_nero/nmds_nero/data/merged_data.rds")

df <- merged_data |>  
  group_by(year, date, plot_id, veg_type) |> 
  summarize(unique_taxa = n_distinct(taxon_code))

ggplot(df, aes(x = year, y = unique_taxa, color = veg_type))+
  geom_point()+
  geom_smooth(method = "lm")+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~")),
    formula = y ~ x,
    parse = TRUE
  )+
  facet_grid(~veg_type)
