# Load packages
library(tidyverse)
library(epiworldR)

# Load test data
test_data <- read_csv("test_data_for_r_validation.csv")

# Configuration
n_days <- 60       # Sequence length

# Display summary of test data
cat("=== Test Data Summary ===\n")
cat("Number of samples:", nrow(test_data), "\n\n")
cat("Simulation parameter ranges:\n")
cat("  Population size (N):", min(test_data$N), "to", max(test_data$N), "\n")
cat("  Incubation days:", round(min(test_data$incub), 2), "to",
    round(max(test_data$incub), 2), "\n")
cat("  Recovery rate:", round(min(test_data$recov), 4), "to",
    round(max(test_data$recov), 4), "\n")
cat("  Prevalence:", round(min(test_data$prevalence), 6), "to",
    round(max(test_data$prevalence), 6), "\n\n")

# Function to run SEIR simulation and extract E->I transitions
simulate_seir <- function(ptran, crate, n_days, N, incub, recov, prev, nsims = 100) {

  model <- ModelSEIRCONN(
    name = "Test validation",
    n = N,
    prevalence = prev,
    contact_rate = crate,
    transmission_rate = ptran,
    incubation_days = incub,
    recovery_rate = recov
  )

  saver <- make_saver("transition")

  run_multiple(
    model,
    ndays = n_days,
    nsims = nsims,
    saver = saver,
    nthreads = 8
  )

  sim_results <- run_multiple_get_results(model, nthreads = 8)

  # Extract E->I transitions
  transitions <- sim_results$transition %>%
    filter(from == "Exposed", to == "Infected") %>%
    group_by(sim_num, date) %>%
    summarise(counts = sum(counts), .groups = "drop")

  # Get median trajectory
  median_trajectory <- transitions %>%
    group_by(date) %>%
    summarise(
      median_count = median(counts),
      lower_ci = quantile(counts, 0.025),
      upper_ci = quantile(counts, 0.975),
      .groups = "drop"
    ) %>%
    pull(median_count)

  return(median_trajectory)
}

# Validate a single sample
validate_sample <- function(sample_id, test_data) {
  row <- test_data[sample_id, ]

  # Extract ALL sample-specific parameters
  N <- row$N
  incub <- row$incub
  recov <- row$recov
  prev <- row$prevalence  # Use the actual prevalence from the test set!

  # Parse incidence sequence (this is the OBSERVED data from simulation)
  observed_incidence <- as.numeric(strsplit(row$incidence_sequence, ",")[[1]])[-1]

  cat("\n=== Sample", sample_id, "===\n")
  cat("Simulation params: N =", N, "| incub =", round(incub, 2),
      "days | recov =", round(recov, 4), "| prev =", round(prev, 6), "\n")
  cat("Actual params: ptran =", round(row$actual_ptran, 4),
      "| crate =", round(row$actual_crate, 3),
      "| R0 =", round(row$actual_R0, 2), "\n")
  cat("Predicted params: ptran =", round(row$pred_ptran, 4),
      "| crate =", round(row$pred_crate, 3),
      "| R0 =", round(row$pred_R0, 2), "\n")

  # Simulate ONLY with PREDICTED parameters
  pred_trajectory <- simulate_seir(
    ptran = row$pred_ptran,
    crate = row$pred_crate,
    n_days = n_days,
    N = N,
    incub = incub,
    recov = recov,
    prev = prev,  # Use the actual prevalence
    nsims = 100
  )

  # Calculate error per day
  mae_per_day <- abs(pred_trajectory - observed_incidence)
  mae_pred <- mean(mae_per_day)

  cat("MAE (predicted params vs observed):", round(mae_pred, 2), "cases/day\n")
  cat("Error range per day:", round(min(mae_per_day), 2), "to",
      round(max(mae_per_day), 2), "cases\n")

  # Plot
  plot_df <- tibble(
    day = 1:n_days,
    observed = observed_incidence,
    predicted_params = pred_trajectory
  )

  p <- ggplot(plot_df, aes(x = day)) +
    geom_line(aes(y = observed, color = "Observed (from actual params)"),
              linewidth = 1.2) +
    geom_point(aes(y = observed, color = "Observed (from actual params)"),
               size = 2, alpha = 0.6) +
    geom_line(aes(y = predicted_params, color = "Simulated (predicted params)"),
              linewidth = 1, linetype = "dashed") +
    scale_color_manual(values = c(
      "Observed (from actual params)" = "blue",
      "Simulated (predicted params)" = "red"
    )) +
    labs(
      title = paste("Sample", sample_id, "- Model Validation"),
      subtitle = sprintf(
        "N=%d, incub=%.1f days, recov=%.3f, prev=%.4f\nActual: ptran=%.4f, crate=%.3f, R0=%.2f | Predicted: ptran=%.4f, crate=%.3f, R0=%.2f",
        N, incub, recov, prev,
        row$actual_ptran, row$actual_crate, row$actual_R0,
        row$pred_ptran, row$pred_crate, row$pred_R0
      ),
      x = "Day",
      y = "Daily E→I Transitions",
      color = ""
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      plot.subtitle = element_text(size = 9)
    )

  print(p)

  return(list(mae = mae_pred, mae_per_day = mae_per_day))
}

# Validate first 5 samples
results <- list()
for (i in 1:min(5, nrow(test_data))) {
  results[[i]] <- validate_sample(i, test_data)
}

# Compute overall statistics across ALL test samples
cat("\n=== Computing Overall Statistics Across All Samples ===\n")

all_maes <- numeric(nrow(test_data))
mae_per_day_matrix <- matrix(NA, nrow = nrow(test_data), ncol = n_days)

for (i in 1:nrow(test_data)) {
  row <- test_data[i, ]

  # Parse observed incidence (remove first day)
  observed_incidence <- as.numeric(strsplit(row$incidence_sequence, ",")[[1]])[-1]

  # Use the actual prevalence from the test set
  prev <- row$prevalence

  # Simulate with predicted parameters
  pred_trajectory <- simulate_seir(
    ptran = row$pred_ptran,
    crate = row$pred_crate,
    n_days = n_days,
    N = row$N,
    incub = row$incub,
    recov = row$recov,
    prev = prev,  # Use actual prevalence
    nsims = 100
  )

  # Calculate MAE per day for this sample
  mae_per_day <- abs(pred_trajectory - observed_incidence)

  # Store in matrix (each row is a sample, each column is a day)
  mae_per_day_matrix[i, ] <- mae_per_day

  # Calculate overall MAE for this sample
  all_maes[i] <- mean(mae_per_day)

  if (i %% 100 == 0) {
    cat("Processed", i, "samples...\n")
  }
}

# Create dataframe with mean MAE per day across all samples
mae_per_day_df <- tibble(
  day = 1:n_days,
  mean_mae = colMeans(mae_per_day_matrix, na.rm = TRUE),
  median_mae = apply(mae_per_day_matrix, 2, median, na.rm = TRUE),
  sd_mae = apply(mae_per_day_matrix, 2, sd, na.rm = TRUE),
  q25_mae = apply(mae_per_day_matrix, 2, quantile, 0.25, na.rm = TRUE),
  q75_mae = apply(mae_per_day_matrix, 2, quantile, 0.75, na.rm = TRUE),
  min_mae = apply(mae_per_day_matrix, 2, min, na.rm = TRUE),
  max_mae = apply(mae_per_day_matrix, 2, max, na.rm = TRUE)
)

# Display the dataframe
cat("\n=== MAE Per Day Summary ===\n")
print(mae_per_day_df)

# Save the dataframe
write_csv(mae_per_day_df, "mae_per_day_summary.csv")
cat("\nSaved MAE per day summary to 'mae_per_day_summary.csv'\n")

# Plot mean MAE per day
p_mae_per_day <- ggplot(mae_per_day_df, aes(x = day)) +
  geom_ribbon(aes(ymin = q25_mae, ymax = q75_mae), fill = "lightblue", alpha = 0.5) +
  geom_line(aes(y = mean_mae, color = "Mean MAE"), linewidth = 1.2) +
  geom_line(aes(y = median_mae, color = "Median MAE"), linewidth = 1, linetype = "dashed") +
  scale_color_manual(values = c("Mean MAE" = "blue", "Median MAE" = "red")) +
  labs(
    title = "Mean Absolute Error per Day Across All Test Samples",
    subtitle = paste0(nrow(test_data), " samples, shaded region = IQR (25th-75th percentile)"),
    x = "Day",
    y = "Mean Absolute Error (cases)",
    color = ""
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

print(p_mae_per_day)

# Save the plot
ggsave("mae_per_day_plot.png", p_mae_per_day, width = 10, height = 6, dpi = 300)

cat("\n=== Overall Test Set Performance ===\n")
cat("Number of test samples:", nrow(test_data), "\n\n")

cat("Parameter prediction errors:\n")
cat("  ptran MAE:", round(mean(abs(test_data$pred_ptran - test_data$actual_ptran)), 4), "\n")
cat("  crate MAE:", round(mean(abs(test_data$pred_crate - test_data$actual_crate)), 3), "\n")
cat("  R0 MAE:", round(mean(abs(test_data$pred_R0 - test_data$actual_R0)), 2), "\n\n")

cat("Incidence prediction performance:\n")
cat("  Mean MAE across all samples:", round(mean(all_maes), 2), "cases/day\n")
cat("  Median MAE:", round(median(all_maes), 2), "cases/day\n")
cat("  MAE range:", round(min(all_maes), 2), "to", round(max(all_maes), 2), "\n\n")

cat("MAE per day statistics:\n")
cat("  Overall mean MAE per day:", round(mean(mae_per_day_df$mean_mae), 2), "cases\n")
cat("  Days with highest mean error:\n")
top_error_days <- mae_per_day_df %>% arrange(desc(mean_mae)) %>% head(5)
print(top_error_days)

cat("\n  Days with lowest mean error:\n")
low_error_days <- mae_per_day_df %>% arrange(mean_mae) %>% head(5)
print(low_error_days)
