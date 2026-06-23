# =====================================================
# Load packages
# =====================================================
library(tidyverse)
library(epiworldR)

# =====================================================
# Load test data
# =====================================================
test_data    <- read_csv("test_data_for_r_validation.csv")
df_incidence <- read_csv("test_data_incidence.csv")

# Parse incidence into a list of integer vectors
incidence_list <- lapply(1:nrow(df_incidence), function(i) {
  as.integer(df_incidence[i, -1])  # exclude sample_id column
})

# Configuration
n_days <- 60

# =====================================================
# Display summary of test data
# =====================================================
cat("=== Test Data Summary ===\n")
cat("Number of samples:", nrow(test_data), "\n\n")
cat("Column names:", paste(names(test_data), collapse = ", "), "\n\n")
cat("Simulation parameter ranges:\n")
cat("  Population size (N):",  min(test_data$N), "to", max(test_data$N), "\n")
cat("  Incubation days:",      round(min(test_data$incub), 2), "to",
    round(max(test_data$incub), 2), "\n")
cat("  Recovery rate:",        round(min(test_data$recov), 4), "to",
    round(max(test_data$recov), 4), "\n")
cat("  Prevalence:",           round(min(test_data$prevalence), 6), "to",
    round(max(test_data$prevalence), 6), "\n\n")

# =====================================================
# Add derived columns
# =====================================================
# beta = ptran * crate  (both predicted and actual)
test_data <- test_data %>%
  mutate(
    pred_beta   = pred_ptran   * pred_crate,
    actual_beta = actual_ptran * actual_crate
  )

cat("=== Derived beta (first 5 rows) ===\n")
test_data %>%
  select(actual_ptran, actual_crate, actual_beta, actual_R0,
         pred_ptran,   pred_crate,   pred_beta,   pred_R0) %>%
  head(5) %>%
  print()

# =====================================================
# Function: Run SEIR simulation
# =====================================================
simulate_seir <- function(ptran, crate, incub, n_days, N, recov, prev,
                          nsims = 100) {
  model <- ModelSEIRCONN(
    name              = "Validation model",
    n                 = N,
    prevalence        = prev,
    contact_rate      = crate,
    transmission_rate = ptran,
    incubation_days   = incub,
    recovery_rate     = recov
  )

  saver <- make_saver("transition")

  run_multiple(
    model,
    ndays    = n_days,
    nsims    = nsims,
    saver    = saver,
    nthreads = 8
  )

  sim_results <- run_multiple_get_results(model, nthreads = 8)

  # Extract E→I transitions
  transitions <- sim_results$transition %>%
    filter(from == "Exposed", to == "Infected") %>%
    group_by(sim_num, date) %>%
    summarise(counts = sum(counts), .groups = "drop")

  # Median trajectory across simulations
  median_trajectory <- transitions %>%
    group_by(date) %>%
    summarise(median_count = median(counts), .groups = "drop") %>%
    arrange(date) %>%
    pull(median_count)

  return(median_trajectory)
}

# =====================================================
# Function: Validate one sample
# =====================================================
validate_sample <- function(sample_id, test_data, incidence_list) {

  row   <- test_data[sample_id, ]
  N     <- row$N
  incub <- row$incub
  recov <- row$recov
  prev  <- row$prevalence

  observed_incidence <- incidence_list[[sample_id]]

  # R0 consistency check: R0 = beta / recov = ptran * crate / recov
  pred_R0_check   <- row$pred_beta   / recov
  actual_R0_check <- row$actual_beta / recov

  cat("\n=== Sample", sample_id, "===\n")
  cat(sprintf("  N = %d | incub = %.2f days | recov = %.4f | prev = %.6f\n",
              N, incub, recov, prev))
  cat(sprintf("  Actual  : ptran = %.4f | crate = %.3f | beta = %.4f | R0 = %.2f\n",
              row$actual_ptran, row$actual_crate, row$actual_beta, row$actual_R0))
  cat(sprintf("  Pred    : ptran = %.4f | crate = %.3f | beta = %.4f | R0 = %.2f\n",
              row$pred_ptran,   row$pred_crate,   row$pred_beta,   row$pred_R0))
  cat(sprintf("  R0 check (beta/recov): pred = %.4f | actual = %.4f\n",
              pred_R0_check, actual_R0_check))

  # Run simulation with predicted parameters
  pred_trajectory <- simulate_seir(
    ptran  = row$pred_ptran,
    crate  = row$pred_crate,
    incub  = incub,
    n_days = n_days,
    N      = N,
    recov  = recov,
    prev   = prev,
    nsims  = 100
  )

  # Remove day 0 (epiworldR includes day 0)
  pred_trajectory <- pred_trajectory[-1]

  # Enforce equal length
  min_len            <- min(length(pred_trajectory), length(observed_incidence))
  pred_trajectory    <- pred_trajectory[1:min_len]
  observed_incidence <- observed_incidence[1:min_len]

  # Compute MAE
  mae_per_day <- abs(pred_trajectory - observed_incidence)
  mae_pred    <- mean(mae_per_day)

  cat(sprintf("  MAE: %.2f cases/day | range: %.2f – %.2f\n",
              mae_pred, min(mae_per_day), max(mae_per_day)))

  # Plot
  plot_df <- tibble(
    day       = 1:min_len,
    observed  = observed_incidence,
    predicted = pred_trajectory
  )

  p <- ggplot(plot_df, aes(x = day)) +
    geom_line(aes(y = observed,  color = "Observed"),
              linewidth = 1.2) +
    geom_point(aes(y = observed, color = "Observed"),
               size = 2, alpha = 0.6) +
    geom_line(aes(y = predicted, color = "Simulated (predicted params)"),
              linewidth = 1, linetype = "dashed") +
    scale_color_manual(values = c(
      "Observed"                     = "blue",
      "Simulated (predicted params)" = "red"
    )) +
    labs(
      title    = paste("Sample", sample_id, "— Model Validation"),
      subtitle = sprintf(
        "N=%d | incub=%.1f days | recov=%.3f | prev=%.4f\nActual: ptran=%.4f  crate=%.3f  beta=%.4f  R0=%.2f\nPred:   ptran=%.4f  crate=%.3f  beta=%.4f  R0=%.2f",
        N, incub, recov, prev,
        row$actual_ptran, row$actual_crate, row$actual_beta, row$actual_R0,
        row$pred_ptran,   row$pred_crate,   row$pred_beta,   row$pred_R0
      ),
      x     = "Day",
      y     = "Daily E→I Transitions",
      color = ""
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      plot.subtitle   = element_text(size = 9)
    )

  print(p)

  return(list(mae = mae_pred, mae_per_day = mae_per_day))
}

# =====================================================
# Validate first 10 samples
# =====================================================
results <- list()

for (i in 1:min(10, nrow(test_data))) {
  results[[i]] <- validate_sample(i, test_data, incidence_list)
}

# =====================================================
# Summary of results
# =====================================================
cat("\n=== Overall Validation Summary ===\n")
mae_values <- sapply(results, function(r) r$mae)
for (i in seq_along(mae_values)) {
  cat(sprintf("  Sample %d: %.2f cases/day\n", i, mae_values[i]))
}
cat(sprintf("\nMean MAE : %.2f cases/day\n", mean(mae_values)))
cat(sprintf("Min  MAE : %.2f cases/day\n", min(mae_values)))
cat(sprintf("Max  MAE : %.2f cases/day\n", max(mae_values)))
