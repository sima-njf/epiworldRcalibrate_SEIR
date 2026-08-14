# =============================================================================
#  fulltest_submit.R  — hypertuned BiLSTM (Optuna)
#
#  STEP 1 — Submit:   Rscript fulltest_submit.R
#  STEP 2 — Collect:  Rscript fulltest_submit.R --collect
#
#  Evaluates the hypertuned BiLSTM on ALL test simulations x ALL windows.
#  One SLURM task per sim_idx. Single ModelSEIRCONN run per evaluation.
#
#  PREREQUISITE: run generate_tuned_predictions.py first to create
#    test_bilstm_predictions_tuned.csv
#
#  NOTE: Day 1 of simulation is excluded from all plots.
# =============================================================================

library(slurmR)
library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)

# =============================================================================
# Config
# =============================================================================
PROJECT_DIR <- path.expand(
  "~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction"
)

SCRATCH <- file.path(
  "/scratch/general/vast",
  Sys.getenv("USER"),
  "slurmR"
)

OUT_DIR <- file.path(PROJECT_DIR, "plots", "fulltest")
dir.create(path.expand(SCRATCH), recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR,              recursive = TRUE, showWarnings = FALSE)

SLURM_OPTS <- list(
  account         = "vegayon-np",
  partition       = "vegayon-np",
  `cpus-per-task` = 1,
  `mem-per-cpu`   = "4G",
  time            = "02:00:00"
)

# =============================================================================
# Detect run mode
# =============================================================================

args         <- commandArgs(trailingOnly = TRUE)
collect_only <- "--collect" %in% args

# =============================================================================
# STEP 1 — Submit
# =============================================================================

if (!collect_only) {

  source(path.expand(file.path(PROJECT_DIR, "fulltest_worker.R")))

  actual    <- read.csv(path.expand(file.path(PROJECT_DIR, "test_actual_parameters (2).csv")))
  preds_all <- read.csv(path.expand(file.path(PROJECT_DIR, "test_bilstm_predictions_tuned.csv")))

  sim_ids <- intersect(unique(preds_all$sim_idx), actual$sim_idx)
  cat(sprintf("Total sims to evaluate: %d\n", length(sim_ids)))
  cat(sprintf("Windows in predictions: %s\n",
              paste(sort(unique(preds_all$window)), collapse = ", ")))

  job <- Slurm_lapply(
    X          = as.list(sim_ids),
    FUN        = eval_one_sim,
    njobs      = min(length(sim_ids), 100),
    mc.cores   = 1,
    job_name   = "seir_fulltest_tuned",
    plan       = "submit",
    sbatch_opt = SLURM_OPTS,
    export     = c("eval_one_sim"),
    tmp_path   = SCRATCH
  )

  cat("\nJobs submitted. When finished, run:\n")
  cat("  Rscript fulltest_submit.R --collect\n\n")
  quit(save = "no")
}

# =============================================================================
# STEP 2 — Collect
# =============================================================================

cat("Collecting results...\n")

job_path <- file.path(path.expand(SCRATCH), "seir_fulltest_tuned")

results_raw <- Slurm_collect(
  read_slurm_job(job_path),
  any. = TRUE
)
results_raw <- Filter(Negate(is.null), results_raw)
cat(sprintf("Non-null results: %d\n", length(results_raw)))

perf_df    <- bind_rows(lapply(results_raw, `[[`, "metrics")) |>
  mutate(regime = factor(regime, levels = c("early", "mid", "late")))

mae_day_df <- bind_rows(lapply(results_raw, `[[`, "mae_day")) |>
  mutate(regime = factor(regime, levels = c("early", "mid", "late")))

cat(sprintf("Metrics rows : %d  (sims: %d, windows: %d)\n",
            nrow(perf_df),
            length(unique(perf_df$sim_idx)),
            length(unique(perf_df$window))))
cat(sprintf("MAE-day rows : %d\n", nrow(mae_day_df)))

write.csv(perf_df,
          file.path(OUT_DIR, "fulltest_all_metrics.csv"),
          row.names = FALSE)
write.csv(mae_day_df,
          file.path(OUT_DIR, "fulltest_mae_per_day.csv"),
          row.names = FALSE)
cat("Saved: fulltest_all_metrics.csv\n")
cat("Saved: fulltest_mae_per_day.csv\n")

# =============================================================================
# STEP 3 — Summary tables
# =============================================================================

LENGTHS <- sort(unique(perf_df$win_len))

summary_win <- perf_df |>
  group_by(window, regime, win_len) |>
  summarise(
    n_sims        = n(),
    mae_beta      = mean(err_beta,      na.rm = TRUE),
    mae_R0        = mean(err_R0,        na.rm = TRUE),
    curve_mae     = mean(curve_mae,     na.rm = TRUE),
    obs_curve_mae = mean(obs_curve_mae, na.rm = TRUE),
    peak_day_err  = mean(peak_day_err,  na.rm = TRUE),
    peak_size_err = mean(peak_size_err, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(regime, win_len)

cat("\n=== Performance by window (all test sims, tuned model) ===\n")
summary_win |>
  mutate(across(where(is.numeric) & !c(n_sims, win_len), round, 3)) |>
  print(n = 40)

write.csv(summary_win,
          file.path(OUT_DIR, "fulltest_summary_by_window.csv"),
          row.names = FALSE)

# Per-day MAE averaged across all sims (exclude day 1)
mae_day_avg <- mae_day_df |>
  filter(day > 1) |>
  group_by(window, regime, win_len, day) |>
  summarise(
    mean_mae      = mean(mae,      na.rm = TRUE),
    sd_mae        = sd(mae,        na.rm = TRUE),
    mean_act_inc  = mean(act_inc,  na.rm = TRUE),
    mean_pred_inc = mean(pred_inc, na.rm = TRUE),
    mean_obs_inc  = mean(obs_inc,  na.rm = TRUE),
    .groups = "drop"
  )

write.csv(mae_day_avg,
          file.path(OUT_DIR, "fulltest_mae_per_day_avg.csv"),
          row.names = FALSE)
cat("Saved: fulltest_mae_per_day_avg.csv\n")

# =============================================================================
# Plots (all starting from day 2)
# =============================================================================

regime_colors <- c(early = "#E65100", mid = "#1565C0", late = "#2E7D32")

# -- Plot 1: Per-day MAE curves — one panel per window -----------------------
p1 <- ggplot(mae_day_avg, aes(x = day)) +
  geom_ribbon(
    aes(ymin = pmax(mean_mae - sd_mae, 0), ymax = mean_mae + sd_mae, fill = regime),
    alpha = 0.18
  ) +
  geom_line(aes(y = mean_mae, color = regime), linewidth = 0.9) +
  facet_wrap(~ window, scales = "free_y", ncol = 6) +
  scale_color_manual(values = regime_colors) +
  scale_fill_manual(values  = regime_colors) +
  scale_y_continuous(labels = comma) +
  labs(
    title    = sprintf("Per-day Incidence MAE — All %d Test Sims [tuned BiLSTM]",
                       length(unique(perf_df$sim_idx))),
    subtitle = "Mean ± SD | Days 2-365 | Each panel = one prediction window",
    x        = "Day",
    y        = "Mean |act incidence − pred incidence|",
    color    = "Regime", fill = "Regime"
  ) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "bottom",
        plot.title       = element_text(face = "bold"),
        strip.text       = element_text(size = 8))

ggsave(file.path(OUT_DIR, "fulltest_perday_mae_by_window.png"),
       p1, width = 20, height = 10, dpi = 150)
cat("Saved: fulltest_perday_mae_by_window.png\n")

# -- Plot 2: Mean incidence curves (actual vs predicted) per window ----------
curves_long <- mae_day_avg |>
  select(window, regime, win_len, day, mean_act_inc, mean_pred_inc, mean_obs_inc) |>
  pivot_longer(c(mean_act_inc, mean_pred_inc, mean_obs_inc),
               names_to = "source", values_to = "incidence") |>
  mutate(source = recode(source,
                         mean_act_inc  = "SEIR - Actual Params",
                         mean_pred_inc = "SEIR - BiLSTM Predicted",
                         mean_obs_inc  = "Observed ABM"))

p2 <- ggplot(curves_long,
             aes(x = day, y = incidence, color = source, linetype = source)) +
  geom_line(linewidth = 0.85, alpha = 0.95) +
  facet_wrap(~ window, scales = "free_y", ncol = 6) +
  scale_color_manual(values = c(
    "Observed ABM"           = "black",
    "SEIR - Actual Params"   = "#1565C0",
    "SEIR - BiLSTM Predicted"= "#C62828"
  )) +
  scale_linetype_manual(values = c(
    "Observed ABM"           = "solid",
    "SEIR - Actual Params"   = "dashed",
    "SEIR - BiLSTM Predicted"= "solid"
  )) +
  scale_y_continuous(labels = comma) +
  labs(
    title    = sprintf("Mean Incidence Curves — All %d Test Sims [tuned BiLSTM]",
                       length(unique(perf_df$sim_idx))),
    subtitle = "Lines = mean across all test sims | Days 2-365 | Each panel = one prediction window",
    x        = "Day",
    y        = "Mean daily incidence",
    color    = NULL, linetype = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "bottom",
        legend.key.width = unit(2, "cm"),
        plot.title       = element_text(face = "bold"),
        strip.text       = element_text(size = 8))

ggsave(file.path(OUT_DIR, "fulltest_mean_incidence_by_window.png"),
       p2, width = 20, height = 10, dpi = 150)
cat("Saved: fulltest_mean_incidence_by_window.png\n")

# -- Plot 3: Overall curve MAE vs window length ------------------------------
p3 <- ggplot(summary_win,
             aes(x = win_len, y = curve_mae, color = regime, group = regime)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  scale_x_log10(breaks = LENGTHS) +
  scale_color_manual(values = regime_colors) +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Overall Curve MAE vs Window Length [tuned BiLSTM]",
    x     = "Window length (days, log scale)",
    y     = "Mean curve MAE (act vs pred incidence)",
    color = "Regime"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title       = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "fulltest_curve_mae_vs_window.png"),
       p3, width = 10, height = 5, dpi = 150)
cat("Saved: fulltest_curve_mae_vs_window.png\n")

# -- Plot 4: Parameter errors vs window length --------------------------------
param_long <- summary_win |>
  select(win_len, regime, mae_beta, mae_R0) |>
  pivot_longer(c(mae_beta, mae_R0), names_to = "param", values_to = "mae") |>
  mutate(param = recode(param, mae_beta = "MAE beta", mae_R0 = "MAE R0"))

p4 <- ggplot(param_long,
             aes(x = win_len, y = mae, color = regime, group = regime)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  facet_wrap(~ param, scales = "free_y", ncol = 2) +
  scale_x_log10(breaks = LENGTHS) +
  scale_color_manual(values = regime_colors) +
  labs(
    title = "Parameter Error vs Window Length [tuned BiLSTM]",
    x     = "Window length (days, log scale)",
    y     = "Mean MAE",
    color = "Regime"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title       = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "fulltest_param_error_vs_window.png"),
       p4, width = 10, height = 5, dpi = 150)
cat("Saved: fulltest_param_error_vs_window.png\n")

# -- Plot 5: R0 error violin --------------------------------------------------
p5 <- ggplot(perf_df,
             aes(x = factor(win_len), y = err_R0, fill = regime)) +
  geom_violin(alpha = 0.6, scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.08, outlier.size = 0.3, alpha = 0.8,
               position = position_dodge(0.9)) +
  facet_wrap(~ regime, ncol = 3) +
  scale_fill_manual(values = regime_colors) +
  labs(
    title = "R0 Error Distribution by Window Length and Regime [tuned BiLSTM]",
    x     = "Window length (days)",
    y     = "|Predicted R0 - Actual R0|",
    fill  = "Regime"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position  = "none",
        panel.grid.minor = element_blank(),
        plot.title       = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "fulltest_R0_error_violin.png"),
       p5, width = 13, height = 5, dpi = 150)
cat("Saved: fulltest_R0_error_violin.png\n")

# =============================================================================
# Done
# =============================================================================

cat(sprintf("\nAll outputs saved under:\n  %s\n", OUT_DIR))
