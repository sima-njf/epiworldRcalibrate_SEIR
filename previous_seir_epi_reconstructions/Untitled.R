# =============================================================================
# SEIR Model Validation with parallel::mclapply
# Cluster: notchpeak | Account/Partition: vegayon-np | 10 cores | 6 hrs
# =============================================================================

library(epiworldR)
library(data.table)
library(parallel)
library(tidyverse)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
n_days <- 60
nsims  <- 100
ncores <- 10

# ---------------------------------------------------------------------------
# Load test data
# ---------------------------------------------------------------------------
test_data <- data.table::fread("test_data_for_r_validation.csv")

cat("=== Test Data Summary ===\n")
cat("Number of samples:", nrow(test_data), "\n")
cat("Simulation parameter ranges:\n")
cat("  Population size (N):", min(test_data$N), "to", max(test_data$N), "\n")
cat("  Incubation days:", round(min(test_data$incub), 2), "to",
    round(max(test_data$incub), 2), "\n")
cat("  Recovery rate:", round(min(test_data$recov), 4), "to",
    round(max(test_data$recov), 4), "\n")
cat("  Prevalence:", round(min(test_data$prevalence), 6), "to",
    round(max(test_data$prevalence), 6), "\n\n")

# ---------------------------------------------------------------------------
# Function: run nsims SEIR simulations for ONE test sample
# Returns median trajectory and MAE vs observed
# ---------------------------------------------------------------------------
run_seir_validation_single <- function(sample_idx, test_data, n_days, nsims, ncores) {

  row <- test_data[sample_idx, ]

  # Parse observed incidence (drop day 0)
  observed_incidence <- as.numeric(strsplit(row$incidence_sequence, ",")[[1]])[-1]

  # Generate seeds for reproducibility
  set.seed(sample_idx * 1000L)
  seeds <- sample(1:1000000, nsims)

  # Run simulations in parallel
  incidence_list <- parallel::mclapply(1:nsims, FUN = function(i) {

    set.seed(seeds[i])

    # Create SEIR model with PREDICTED parameters
    m <- epiworldR::ModelSEIRCONN(
      name              = "Validation",
      n                 = row$N,
      prevalence        = row$prevalence,
      contact_rate      = row$pred_crate,
      incubation_days   = row$incub,
      transmission_rate = row$pred_ptran,
      recovery_rate     = row$recov
    )

    # Turn off verbose output
    epiworldR::verbose_off(m)

    # Run the simulation
    epiworldR::run(m, ndays = n_days)

    # Extract incidence data
    incidence <- epiworldR::plot_incidence(m, plot = FALSE)
    incidence_dt <- data.table::as.data.table(incidence)

    # Extract infected counts (skip day 0)
    infected_vector <- incidence_dt$Infected[2:(n_days + 1)]
    return(infected_vector)

  }, mc.cores = ncores)

  # Convert list to matrix (each row = one simulation, each column = one day)
  infected_matrix <- do.call(rbind, incidence_list)

  # Get median trajectory across all simulations
  median_trajectory <- apply(infected_matrix, 2, median)

  # Calculate per-day MAE
  mae_per_day <- abs(median_trajectory - observed_incidence)

  return(list(
    sample_id         = sample_idx,
    mae               = mean(mae_per_day),
    mae_per_day       = mae_per_day,
    median_trajectory = median_trajectory
  ))
}

# ---------------------------------------------------------------------------
# Run validation across ALL test samples
# ---------------------------------------------------------------------------
N_SAMPLES <- nrow(test_data)

cat("==============================================\n")
cat("Running validation for", N_SAMPLES, "samples\n")
cat("Each sample:", nsims, "simulations on", ncores, "cores\n")
cat("==============================================\n")

all_results <- vector("list", N_SAMPLES)

for (i in 1:N_SAMPLES) {

  all_results[[i]] <- run_seir_validation_single(
    sample_idx = i,
    test_data  = test_data,
    n_days     = n_days,
    nsims      = nsims,
    ncores     = ncores
  )

  if (i %% 100 == 0) {
    cat("Processed", i, "/", N_SAMPLES, "samples...\n")
  }
}

cat("\n==============================================\n")
cat("All samples processed!\n")
cat("==============================================\n\n")

# ---------------------------------------------------------------------------
# Collect & aggregate results
# ---------------------------------------------------------------------------
all_maes           <- vapply(all_results, function(r) r$mae, numeric(1))
mae_per_day_matrix <- do.call(rbind, lapply(all_results, function(r) r$mae_per_day))
saveRDS(all_results, file = "all_results.rds")
cat("Saved all_results to 'all_results.rds'\n")
# Per-day summary table
mae_per_day_df <- tibble(
  day        = 1:n_days,
  mean_mae   = colMeans(mae_per_day_matrix, na.rm = TRUE),
  median_mae = apply(mae_per_day_matrix, 2, median,   na.rm = TRUE),
  sd_mae     = apply(mae_per_day_matrix, 2, sd,       na.rm = TRUE),
  q25_mae    = apply(mae_per_day_matrix, 2, quantile, 0.25, na.rm = TRUE),
  q75_mae    = apply(mae_per_day_matrix, 2, quantile, 0.75, na.rm = TRUE),
  min_mae    = apply(mae_per_day_matrix, 2, min,      na.rm = TRUE),
  max_mae    = apply(mae_per_day_matrix, 2, max,      na.rm = TRUE)
)

cat("\n=== MAE Per Day Summary ===\n")
print(mae_per_day_df)
max(mae_per_day_df)
data.table::fwrite(mae_per_day_df, file = "mae_per_day_summary.csv")
cat("\nSaved MAE per day summary to 'mae_per_day_summary.csv'\n")

# ---------------------------------------------------------------------------
# Plot: Mean MAE per day
# ---------------------------------------------------------------------------
p_mae_per_day <- ggplot(mae_per_day_df, aes(x = day)) +
  geom_ribbon(aes(ymin = q25_mae, ymax = q75_mae),
              fill = "lightblue", alpha = 0.5) +
  geom_line(aes(y = mean_mae,   color = "Mean MAE"),   linewidth = 1.2) +
  geom_line(aes(y = median_mae, color = "Median MAE"),
            linewidth = 1, linetype = "dashed") +
  scale_color_manual(values = c("Mean MAE" = "blue", "Median MAE" = "red")) +
  labs(
    title    = "Mean Absolute Error per Day Across All Test Samples",
    subtitle = paste0(N_SAMPLES,
                      " samples, shaded region = IQR (25th-75th percentile)"),
    x     = "Day",
    y     = "Mean Absolute Error (cases)",
    color = ""
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

print(p_mae_per_day)
ggsave("mae_per_day_plot.png", p_mae_per_day, width = 10, height = 6, dpi = 300)

# ---------------------------------------------------------------------------
# Validation plots for first 5 samples
# ---------------------------------------------------------------------------
for (i in 1:min(5, N_SAMPLES)) {

  row <- test_data[i, ]
  observed_incidence <- as.numeric(strsplit(row$incidence_sequence, ",")[[1]])[-1]
  pred_trajectory    <- all_results[[i]]$median_trajectory

  plot_df <- tibble(
    day              = 1:n_days,
    observed         = observed_incidence,
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
      "Simulated (predicted params)"  = "red"
    )) +
    labs(
      title    = paste("Sample", i, "- Model Validation"),
      subtitle = sprintf(
        "N=%d, incub=%.1f, recov=%.3f, prev=%.4f\nActual: ptran=%.4f, crate=%.3f, R0=%.2f | Predicted: ptran=%.4f, crate=%.3f, R0=%.2f",
        row$N, row$incub, row$recov, row$prevalence,
        row$actual_ptran, row$actual_crate, row$actual_R0,
        row$pred_ptran, row$pred_crate, row$pred_R0
      ),
      x     = "Day",
      y     = "Daily E→I Transitions",
      color = ""
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", plot.subtitle = element_text(size = 9))

  print(p)
}

# ---------------------------------------------------------------------------
# Final summary
# ---------------------------------------------------------------------------
cat("\n==============================================\n")
cat("Overall Test Set Performance\n")
cat("==============================================\n")
cat("Number of test samples:", N_SAMPLES, "\n\n")

cat("Parameter prediction errors:\n")
cat("  ptran MAE:", round(mean(abs(test_data$pred_ptran - test_data$actual_ptran)), 4), "\n")
cat("  crate MAE:", round(mean(abs(test_data$pred_crate - test_data$actual_crate)), 3), "\n")
cat("  R0 MAE:",    round(mean(abs(test_data$pred_R0    - test_data$actual_R0)),    2), "\n\n")

cat("Incidence prediction performance:\n")
cat("  Mean MAE across all samples:", round(mean(all_maes), 2), "cases/day\n")
cat("  Median MAE:", round(median(all_maes), 2), "cases/day\n")
cat("  MAE range:", round(min(all_maes), 2), "to", round(max(all_maes), 2), "\n\n")

cat("MAE per day statistics:\n")
cat("  Overall mean MAE per day:", round(mean(mae_per_day_df$mean_mae), 2), "cases\n")

cat("\n  Days with highest mean error:\n")
print(mae_per_day_df %>% arrange(desc(mean_mae)) %>% head(5))

cat("\n  Days with lowest mean error:\n")
print(mae_per_day_df %>% arrange(mean_mae) %>% head(5))

cat("\n==============================================\n")
cat("Validation complete!\n")
cat("Matrix dimensions:", dim(mae_per_day_matrix), "\n")
cat("==============================================\n")

