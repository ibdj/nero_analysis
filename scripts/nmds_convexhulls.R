library(ggplot2)
library(dplyr)

# Create a data frame for convex hulls
hulls <- nmds_plot_data |>
  group_by(year) |>
  slice(chull(NMDS1, NMDS2)) |>
  ungroup()

hulls

# Plot NMDS with convex hulls
ggplot(nmds_plot_data, aes(x = NMDS1, y = NMDS2, color = as.factor(year))) +
  geom_point(size = 3) +  # Plot points
  geom_polygon(data = hulls, aes(group = year, fill = as.factor(year)), alpha = 0.3) +  # Add convex hulls
  labs(
    title = "NMDS with Convex Hulls by Year",
    x = "NMDS Dimension 1",
    y = "NMDS Dimension 2",
    color = "Year",
    fill = "Year"
  ) +
  theme_minimal()

install.packages("sp")  # Install if not already installed
library(sp)

# Calculate convex hull areas for each year
hull_areas <- nmds_plot_data |>
  group_by(year) |>
  summarise(
    hull_area = {
      hull_points <- slice(nmds_plot_data, chull(NMDS1, NMDS2)) |> select(NMDS1, NMDS2)
      polygon <- Polygon(as.matrix(hull_points), hole = FALSE)
      polygon@area
    }
  )

# View the calculated hull areas
print(hull_areas)
