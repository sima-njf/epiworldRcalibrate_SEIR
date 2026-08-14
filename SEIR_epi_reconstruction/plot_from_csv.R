# =============================================================================
#  plot_from_csv.R  — regenerate all comparison plots from saved CSVs
#
#  Use this when compare_all_methods.R was interrupted OR you just want
#  to regenerate plots without re-running SEIR.
#
#  Priority:
#    1. Loads comparison_5method_metrics.csv if it exists (fastest)
#    2. Otherwise reconstructs metrics from individual method summary CSVs
#       (BiLSTM curve metrics will be NA but param-error plots still work)
#
#  Run:
#    Rscript plot_from_csv.R
# =============================================================================

PROJECT_DIR <- path.expand(
  "~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction"
)

source(file.path(PROJECT_DIR, "seir_common.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(cowplot)
})

BILSTM_WINDOW <- "late_365d"
METHOD_ORDER  <- SEIR_METHOD_ORDER
METHOD_COLORS <- SEIR_METHOD_COLORS[names(SEIR_METHOD_COLORS) != "Oracle"]
PLOTS_DIR     <- file.path(PROJECT_DIR, "plots", "comparison")
dir.create(PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Load or reconstruct metrics_all
# =============================================================================

metrics_path <- file.path(PROJECT_DIR, "comparison_5method_metrics.csv")

if (file.exists(metrics_path)) {

  cat("Loading comparison_5method_metrics.csv...\n")
  metrics_all <- read.csv(metrics_path)

} else {

  cat("comparison_5method_metrics.csv not found — reconstructing from method summaries.\n")

  abc_summary    <- read.csv(file.path(PROJECT_DIR, "abc_seir_summary.csv"))
  abcsmc_summary <- read.csv(file.path(PROJECT_DIR, "abcsmc_seir_summary.csv"))
  nm_summary     <- read.csv(file.path(PROJECT_DIR, "nm_seir_summary.csv"))
  de_summary     <- read.csv(file.path(PROJECT_DIR, "de_seir_summary.csv"))
  actual         <- read.csv(file.path(PROJECT_DIR, SEIR_TEST_PARAMS_FILE))
  preds_all      <- read.csv(file.path(PROJECT_DIR, "test_bilstm_predictions_tuned.csv"))

  want_cols <- c("sim_idx", "method", "pred_beta", "pred_R0",
                 "true_beta", "true_R0", "mae", "smape", "pearson_r",
                 "peak_day_err", "rel_peak_err", "err_beta_pct", "err_R0_pct",
                 "time_sec", "cpu_sec", "n_evals", "n_model_runs",
                 "converged", "dieout_pred")

  clean <- function(df, mth) {
    df <- df[df$converged %in% TRUE & !is.na(df$pred_beta), ]
    df$method <- mth
    df[, intersect(want_cols, names(df))]
  }

  # BiLSTM: pred params only; derive param errors from actual
  bl_raw <- preds_all[preds_all$window == BILSTM_WINDOW, ] |>
    left_join(actual |> select(sim_idx, true_beta = beta, true_R0 = R0),
              by = "sim_idx") |>
    transmute(
      sim_idx      = sim_idx,
      method       = "BiLSTM",
      pred_beta    = beta_pred,
      pred_R0      = R0_pred,
      true_beta    = true_beta,
      true_R0      = true_R0,
      err_beta_pct = abs(beta_pred - true_beta) / true_beta * 100,
      err_R0_pct   = abs(R0_pred   - true_R0)   / true_R0   * 100
    )

  metrics_all <- bind_rows(
    clean(abc_summary,    "ABC"),
    clean(abcsmc_summary, "ABC-SMC"),
    clean(nm_summary,     "NelderMead"),
    clean(de_summary,     "DE"),
    bl_raw
  )

  cat(sprintf("Reconstructed %d rows. Note: BiLSTM curve metrics (mae/smape/pearson_r) are NA.\n",
              nrow(metrics_all)))
}

metrics_all <- metrics_all |>
  mutate(method = factor(method, levels = METHOD_ORDER))

n_sims <- n_distinct(metrics_all$sim_idx)
cat(sprintf("Plotting %d unique sims, %d methods.\n",
            n_sims, n_distinct(metrics_all$method)))

# =============================================================================
# Cost table
# =============================================================================

cost_path <- file.path(PROJECT_DIR, "comparison_5method_cost.csv")

if (file.exists(cost_path)) {
  cost_tbl <- read.csv(cost_path) |>
    mutate(method = factor(method, levels = METHOD_ORDER))
} else {
  # Load summaries if not already in memory
  if (!exists("abc_summary"))    abc_summary    <- read.csv(file.path(PROJECT_DIR, "abc_seir_summary.csv"))
  if (!exists("abcsmc_summary")) abcsmc_summary <- read.csv(file.path(PROJECT_DIR, "abcsmc_seir_summary.csv"))
  if (!exists("nm_summary"))     nm_summary     <- read.csv(file.path(PROJECT_DIR, "nm_seir_summary.csv"))
  if (!exists("de_summary"))     de_summary     <- read.csv(file.path(PROJECT_DIR, "de_seir_summary.csv"))

  cost_tbl <- bind_rows(
    abc_summary    |> transmute(sim_idx, method = "ABC",        time_sec, cpu_sec, n_evals, n_model_runs),
    abcsmc_summary |> transmute(sim_idx, method = "ABC-SMC",    time_sec, cpu_sec, n_evals, n_model_runs),
    nm_summary     |> transmute(sim_idx, method = "NelderMead", time_sec, cpu_sec, n_evals, n_model_runs),
    de_summary     |> transmute(sim_idx, method = "DE",         time_sec, cpu_sec, n_evals, n_model_runs)
  ) |>
    group_by(method) |>
    summarise(mean_wall_sec   = mean(time_sec,     na.rm = TRUE),
              median_wall_sec = median(time_sec,   na.rm = TRUE),
              mean_cpu_sec    = mean(cpu_sec,      na.rm = TRUE),
              mean_evals      = mean(n_evals,      na.rm = TRUE),
              mean_model_runs = mean(n_model_runs, na.rm = TRUE),
              .groups = "drop")

  timing_path  <- file.path(PROJECT_DIR, "bilstm_timing.csv")
  bl_mean_wall <- if (file.exists(timing_path))
    mean(read.csv(timing_path)$mean_time_sec, na.rm = TRUE) else NA_real_

  cost_tbl <- cost_tbl |>
    bind_rows(data.frame(method = "BiLSTM", mean_wall_sec = bl_mean_wall,
                         median_wall_sec = NA_real_, mean_cpu_sec = NA_real_,
                         mean_evals = 0, mean_model_runs = 0,
                         stringsAsFactors = FALSE)) |>
    mutate(method = factor(method, levels = METHOD_ORDER)) |>
    arrange(method)
}

bl_mean_wall <- cost_tbl$mean_wall_sec[cost_tbl$method == "BiLSTM"]

# =============================================================================
# Summary
# =============================================================================

comparison_summary <- metrics_all |>
  group_by(method) |>
  summarise(
    n_sims            = n(),
    mean_mae          = mean(mae,          na.rm = TRUE),
    median_mae        = median(mae,        na.rm = TRUE),
    mean_smape        = mean(smape,        na.rm = TRUE),
    median_smape      = median(smape,      na.rm = TRUE),
    mean_pearson_r    = mean(pearson_r,    na.rm = TRUE),
    mean_peak_day_err = mean(peak_day_err, na.rm = TRUE),
    mean_rel_peak_err = mean(rel_peak_err, na.rm = TRUE),
    mean_err_beta_pct = mean(err_beta_pct, na.rm = TRUE),
    mean_err_R0_pct   = mean(err_R0_pct,   na.rm = TRUE),
    n_dieout          = sum(dieout_pred,   na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(cost_tbl |> select(method, mean_wall_sec, mean_model_runs),
            by = "method")

cat("\n========== 5-Method Comparison Summary ==========\n")
print(as.data.frame(comparison_summary), digits = 4, row.names = FALSE)

write.csv(comparison_summary,
          file.path(PROJECT_DIR, "comparison_5method_summary.csv"),
          row.names = FALSE)

# =============================================================================
# Plot 1: all metrics boxplots
# =============================================================================

metrics_long <- metrics_all |>
  select(sim_idx, method, any_of(c("mae", "smape", "pearson_r", "peak_day_err",
                                    "rel_peak_err", "err_beta_pct", "err_R0_pct"))) |>
  pivot_longer(-c(sim_idx, method), names_to = "metric", values_to = "value") |>
  filter(!is.na(value))

metric_labels <- c(
  mae          = "MAE (incidence)",
  smape        = "sMAPE (%)",
  pearson_r    = "Pearson r (shape)",
  peak_day_err = "Peak Timing Error (days)",
  rel_peak_err = "Relative Peak Error (%)",
  err_beta_pct = "beta Error (%)",
  err_R0_pct   = "R0 Error (%)"
)
metrics_long$metric <- factor(metric_labels[metrics_long$metric],
                              levels = unname(metric_labels))

p1 <- ggplot(metrics_long, aes(x = method, y = value, fill = method)) +
  geom_boxplot(alpha = 0.65, outlier.size = 0.7) +
  geom_jitter(width = 0.15, size = 1.6, alpha = 0.55, color = "black") +
  facet_wrap(~ metric, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = METHOD_COLORS) +
  labs(title    = "5-Method Calibration Comparison (SEIR)",
       subtitle = sprintf("N=%d sims | E->I incidence | days 2-365", n_sims),
       x = NULL, y = "Value", fill = "Method") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold"), strip.text = element_text(size = 9),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 7))

ggsave(file.path(PLOTS_DIR, "comparison_5method_all_metrics.png"),
       p1, width = 16, height = 9, dpi = 150)
cat("Saved: comparison_5method_all_metrics.png\n")

# =============================================================================
# Plot 2: sMAPE violin + boxplot
# =============================================================================

smape_data <- metrics_all |> filter(!is.na(smape))
if (nrow(smape_data) > 0) {
  p2 <- ggplot(smape_data, aes(x = method, y = smape, fill = method)) +
    geom_violin(alpha = 0.55, scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.10, outlier.size = 0.7, alpha = 0.85) +
    geom_jitter(width = 0.07, size = 1.5, alpha = 0.5, color = "black") +
    scale_fill_manual(values = METHOD_COLORS) +
    labs(title    = "sMAPE Distribution by Calibration Method (SEIR)",
         subtitle = sprintf("N=%d sims | E->I incidence | days 2-365", n_sims),
         x = NULL, y = "sMAPE (%)", fill = "Method") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "none",
          plot.title = element_text(face = "bold"))

  ggsave(file.path(PLOTS_DIR, "comparison_5method_smape.png"),
         p2, width = 8, height = 5, dpi = 150)
  cat("Saved: comparison_5method_smape.png\n")
}

# =============================================================================
# Plot 3: predicted vs true beta / R0
# =============================================================================

param_scatter <- metrics_all |>
  select(sim_idx, method, true_beta, true_R0, pred_beta, pred_R0) |>
  filter(!is.na(pred_beta))

p3a <- ggplot(param_scatter,
              aes(x = true_beta, y = pred_beta, color = method, shape = method)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_color_manual(values = METHOD_COLORS) +
  labs(title = "Predicted vs True beta",
       x = "True beta", y = "Predicted beta",
       color = "Method", shape = "Method") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

p3b <- ggplot(param_scatter,
              aes(x = true_R0, y = pred_R0, color = method, shape = method)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_color_manual(values = METHOD_COLORS) +
  labs(title = "Predicted vs True R0",
       x = "True R0", y = "Predicted R0",
       color = "Method", shape = "Method") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

ggsave(file.path(PLOTS_DIR, "comparison_5method_param_scatter.png"),
       plot_grid(p3a, p3b, ncol = 2), width = 13, height = 5, dpi = 150)
cat("Saved: comparison_5method_param_scatter.png\n")

# =============================================================================
# Plot 4: accuracy vs cost
# =============================================================================

acc_cost <- comparison_summary |>
  select(method, mean_smape, mean_wall_sec) |>
  filter(!is.na(mean_wall_sec) & !is.na(mean_smape))

if (nrow(acc_cost) > 0) {
  p5 <- ggplot(acc_cost, aes(x = mean_wall_sec, y = mean_smape,
                             color = method, label = method)) +
    geom_point(size = 4) +
    geom_text(vjust = -1.1, size = 3.5, show.legend = FALSE) +
    scale_x_log10() +
    scale_color_manual(values = METHOD_COLORS) +
    labs(title    = "Accuracy vs Computation Cost per Calibrated Simulation",
         subtitle = if (is.na(bl_mean_wall))
           "BiLSTM: amortised inference (run generate_tuned_predictions.py for timing)"
         else
           sprintf("BiLSTM inference: %.1f ms/sim (amortised)",
                   bl_mean_wall * 1000),
         x = "Mean wall time per sim (s, log scale)",
         y = "Mean sMAPE (%)", color = "Method") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"))

  ggsave(file.path(PLOTS_DIR, "comparison_5method_accuracy_vs_cost.png"),
         p5, width = 8, height = 5.5, dpi = 150)
  cat("Saved: comparison_5method_accuracy_vs_cost.png\n")
}

# =============================================================================
# Plot 5: mean metric bars
# =============================================================================

mean_metrics <- comparison_summary |>
  select(method,
         "MAE"                = mean_mae,
         "sMAPE (%)"          = mean_smape,
         "1 - Pearson r"      = mean_pearson_r,
         "Peak Timing (days)" = mean_peak_day_err,
         "Rel. Peak Err (%)"  = mean_rel_peak_err,
         "beta Err (%)"       = mean_err_beta_pct,
         "R0 Err (%)"         = mean_err_R0_pct) |>
  mutate(`1 - Pearson r` = 1 - `1 - Pearson r`) |>
  pivot_longer(-method, names_to = "metric", values_to = "value") |>
  filter(!is.na(value))

p6 <- ggplot(mean_metrics, aes(x = method, y = value, fill = method)) +
  geom_col(alpha = 0.75, width = 0.6) +
  facet_wrap(~ metric, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = METHOD_COLORS) +
  labs(title    = "Mean Metric Comparison: 5 Methods (SEIR)",
       subtitle = "Lower is better for all metrics except Pearson r",
       x = NULL, y = "Mean value", fill = "Method") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 7),
        strip.text = element_text(size = 9))

ggsave(file.path(PLOTS_DIR, "comparison_5method_mean_bar.png"),
       p6, width = 16, height = 8, dpi = 150)
cat("Saved: comparison_5method_mean_bar.png\n")

# =============================================================================
# Plot 6: timing bar
# =============================================================================

timing_long <- cost_tbl |>
  filter(!is.na(mean_wall_sec)) |>
  mutate(label = ifelse(method == "BiLSTM",
                        sprintf("%.3f s\n(amortised)", mean_wall_sec),
                        sprintf("%.1f s", mean_wall_sec)))

p7 <- ggplot(timing_long, aes(x = method, y = mean_wall_sec, fill = method)) +
  geom_col(alpha = 0.80, width = 0.6) +
  geom_text(aes(label = label), vjust = -0.4, size = 3.2) +
  scale_y_log10() +
  scale_fill_manual(values = METHOD_COLORS) +
  labs(title    = "Computation Time per Calibrated Simulation",
       subtitle = "Log scale | BiLSTM = forward pass only (training is amortised)",
       x = NULL, y = "Mean wall time (s, log scale)", fill = "Method") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "none",
        plot.title = element_text(face = "bold"))

ggsave(file.path(PLOTS_DIR, "comparison_5method_timing.png"),
       p7, width = 7, height = 5, dpi = 150)
cat("Saved: comparison_5method_timing.png\n")

cat(sprintf("\nAll plots saved to: %s\n", PLOTS_DIR))
