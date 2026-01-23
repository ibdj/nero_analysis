library(vegan)


# Looking at dissimilariry between years

# Calculate pairwise Jaccard dissimilarities for each year ###
dissimilarity_matrices <- lapply(yearly_data, function(data) {
  # Remove non-species columns (e.g., plot_id, year)
  data_numeric <- data |>
    select(-plot_id, -year) |>
    mutate(across(everything(), ~ ifelse(is.na(as.numeric(.)), 0, as.numeric(.)))) |>  # Replace non-numeric values with 0
    as.matrix()
  
  # Calculate Jaccard dissimilarity
  vegdist(data_numeric, method = "jaccard")
})

# Calculate mean pairwise dissimilarity for each year
mean_dissimilarities <- sapply(dissimilarity_matrices, function(matrix) {
  mean(as.numeric(matrix), na.rm = TRUE)
})

# Convert the result to a data frame for easier handling
dissimilarity_summary <- tibble(
  year = as.numeric(names(mean_dissimilarities)),
  mean_dissimilarity = mean_dissimilarities
)

# View the summary
print(dissimilarity_summary)
 
ggplot(dissimilarity_summary, aes(x = year, y = mean_dissimilarity))+
  geom_point()+
  geom_smooth(method = lm)
