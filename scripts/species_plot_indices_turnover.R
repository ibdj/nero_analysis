library(dplyr)
library(vegan)

library(dplyr)
library(vegan)

library(dplyr)
library(vegan)

# Identify species columns (exclude metadata)
species_cols <- setdiff(names(filtered_matrix), c("year", "plot_id", "veg_type"))

turnover_results <- filtered_matrix |> filter(veg_type == "heath") %>%
  group_by(plot_id) %>%
  arrange(year, .by_group = TRUE) %>%
  group_modify(~ {
    n <- nrow(.x)
    if (n < 2) {
      # Return a tibble with plot_id and NAs for other columns
      return(tibble(plot_id = .x$plot_id[1], year1 = NA_real_, year2 = NA_real_, turnover = NA_real_))
    }
    sp_data <- as.matrix(.x[, species_cols])
    dist_mat <- as.matrix(vegdist(sp_data, method = "jaccard", binary = TRUE))
    turnovers <- sapply(1:(n - 1), function(i) dist_mat[i, i + 1])
    tibble(
      plot_id = .x$plot_id[1],
      year1 = .x$year[1:(n - 1)],
      year2 = .x$year[2:n],
      turnover = turnovers
    )
  }) %>%
  ungroup() %>%
  filter(!is.na(turnover))  # Remove rows where turnover couldn't be calculated

turnover_results <- turnover_results %>%
  mutate(
    period = case_when(
      year1 == 2007 & year2 == 2012 ~ 1,
      year1 == 2012 & year2 == 2017 ~ 2,
      year1 == 2017 & year2 == 2022 ~ 3,
      TRUE ~ NA_integer_  # Assign NA to non-matching pairs
    )
  )

ggplot(turnover_results, aes(x = period, y = turnover))+
  geom_point(color = "darkgreen")+
  geom_smooth(method = lm)


ggplot(turnover_results |> filter(!is.na(period)), aes(x = factor(period), y = turnover)) +
  geom_boxplot(aes(group = factor(period))) +  # Force grouping by period
  geom_jitter(width = 0.05, alpha = 0.5) +     # Better than geom_point for overlapping values
  scale_x_discrete(
    labels = c("1" = "2007-2012", "2" = "2012-2017", "3" = "2017-2022")
  ) +
  labs(x = "", y = "Turnover (Jaccard dissimilarity)") +
  theme_minimal()

turnover_results$period <- factor(turnover_results$period)

# Fit the model
model <- lm(turnover ~ period, data = turnover_results)

# View the summary
summary(model)

#### basic stats ######

stats <- turnover_results %>%
  group_by(period) %>%
  summarise(
    mean_turnover = mean(turnover, na.rm = TRUE),
    n_zero = sum(turnover == 0, na.rm = TRUE),
    n_total = n_distinct(plot_id),
    percent_zero = 100 * n_zero / n_total
  )



summary(turnover_results)
