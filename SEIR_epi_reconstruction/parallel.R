# =============================================================================
#  PART 4 — Overall performance on ALL test data x ALL 18 windows (Slurm)
#
#  Usage (on login node):
#    Rscript part4_slurm_submit.R
# =============================================================================

library(slurmR)
library(dplyr)

# -- Config -------------------------------------------------------------------
DATA_DIR     <- "~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction"
SCRATCH_DIR  <- file.path(DATA_DIR, "part4_scratch")
NSIMS_WORKER <- 100L
NDAYS        <- 365L
CONTACT_RATE <- 1

CURVE_WINDOWS <- c("early_015d", "early_060d", "early_180d",
                   "mid_090d",   "mid_180d",   "late_180d",  "late_365d")

LENGTHS <- c(15, 30, 60, 90, 180, 365)

make_windows <- function(t_max = 365, lengths = LENGTHS) {
  out <- list()
  for (L in lengths) {
    if (L > t_max) next
    out[[sprintf("early_%03dd", L)]] <- c(start = 0L,                        len = L)
    out[[sprintf("mid_%03dd",   L)]] <- c(start = max(0L, as.integer((t_max - L) / 2)), len = L)
    out[[sprintf("late_%03dd",  L)]] <- c(start = max(0L, t_max - L),        len = L)
  }
  out
}
TEST_WINDOWS <- make_windows()

dir.create(path.expand(SCRATCH_DIR), recursive = TRUE, showWarnings = FALSE)

# -- Load IDs (needed here to size the array) ---------------------------------
library(slurmR)
library(dplyr)

DATA_DIR    <- "~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction"
SCRATCH_DIR <- file.path(DATA_DIR, "part4_scratch")

actual     <- read.csv(path.expand(file.path(DATA_DIR, "test_actual_parameters.csv")))
preds_all  <- read.csv(path.expand(file.path(DATA_DIR, "test_bilstm_predictions.csv")))
common_ids <- intersect(unique(preds_all$sim_idx), actual$sim_idx)

slurm_opts <- list(
  account         = "vegayon-np",
  partition       = "vegayon-np",
  time            = "10:00:00",
  `mem-per-cpu`   = "4G",
  `cpus-per-task` = 1
)

# Redefine the function fresh in this session
eval_one_sim <- function(sid) {

  sid <- as.integer(sid)   # fix: coerce list element to scalar

  suppressPackageStartupMessages({
    library(epiworldR)
    library(dplyr)
    library(data.table)
  })

  DATA_DIR      <- "~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction"
  NDAYS         <- 365L
  CONTACT_RATE  <- 1
  NSIMS_WORKER  <- 100L
  CURVE_WINDOWS <- c("early_015d", "early_060d", "early_180d",
                     "mid_090d",   "mid_180d",   "late_180d",  "late_365d")
  LENGTHS       <- c(15, 30, 60, 90, 180, 365)

  make_windows <- function(t_max = 365, lengths = LENGTHS) {
    out <- list()
    for (L in lengths) {
      if (L > t_max) next
      out[[sprintf("early_%03dd", L)]] <- c(start = 0L,                                   len = L)
      out[[sprintf("mid_%03dd",   L)]] <- c(start = max(0L, as.integer((t_max - L) / 2)), len = L)
      out[[sprintf("late_%03dd",  L)]] <- c(start = max(0L, t_max - L),                   len = L)
    }
    out
  }
  TEST_WINDOWS <- make_windows()

  actual    <- read.csv(path.expand(file.path(DATA_DIR, "test_actual_parameters.csv")))
  preds_all <- read.csv(path.expand(file.path(DATA_DIR, "test_bilstm_predictions.csv")))

  np      <- reticulate::import("numpy")
  inc_raw <- as.matrix(
    np$load(path.expand(file.path(DATA_DIR, "test_incidence_raw.npy"))))

  .run_multi_worker <- function(n, prevalence, beta, incub, recov, seed,
                                ndays = NDAYS, nsims = NSIMS_WORKER,
                                cr = CONTACT_RATE) {
    n_int <- max(as.integer(round(n)), 10L)
    prev  <- max(min(prevalence, 1.0), 1.0 / n_int)
    model <- ModelSEIRCONN(
      name              = "sim",
      n                 = n_int,
      prevalence        = prev,
      contact_rate      = cr,
      transmission_rate = max(beta / cr, 1e-6),
      incubation_days   = max(incub, 1.0),
      recovery_rate     = max(recov, 1e-6)
    )
    verbose_off(model)
    saver <- make_saver("transition")
    run_multiple(model, ndays = ndays + 1L, nsims = nsims,
                 saver = saver, nthreads = 1L)
    res <- run_multiple_get_results(model, nthreads = 1L,
                                    freader = data.table::fread)
    res$transition |>
      filter(from == "Exposed", to == "Infected", date > 0) |>
      group_by(date) |>
      summarise(lower = quantile(counts, 0.025),
                med   = quantile(counts, 0.500),
                upper = quantile(counts, 0.975),
                .groups = "drop") |>
      right_join(data.frame(date = 1:ndays), by = "date") |>
      arrange(date) |>
      mutate(across(c(lower, med, upper), \(x) replace_na(x, 0)))
  }

  act     <- actual[actual$sim_idx == sid, ]
  mat_row <- which(actual$sim_idx == sid)
  if (nrow(act) == 0 || length(mat_row) == 0) return(NULL)

  n          <- act$n[1];         recov      <- act$recov[1]
  incub      <- act$incub[1];     prevalence <- act$prevalence[1]
  obs_full   <- as.numeric(inc_raw[mat_row, ])

  act_q <- tryCatch(
    .run_multi_worker(n, prevalence, act$beta[1], incub, recov, seed = sid),
    error = function(e) NULL)
  if (is.null(act_q)) return(NULL)

  metric_rows <- list()
  curve_rows  <- list()

  for (win_tag in names(TEST_WINDOWS)) {
    w       <- TEST_WINDOWS[[win_tag]]
    regime  <- sub("_.*", "", win_tag)
    win_len <- w["len"]
    t_start <- w["start"] + 1L
    t_end   <- min(t_start + win_len - 1L, 365L)

    prd <- preds_all[preds_all$sim_idx == sid & preds_all$window == win_tag, ]
    if (nrow(prd) == 0) next

    pred_q <- tryCatch(
      .run_multi_worker(n, prevalence, prd$beta_pred[1], incub, recov, seed = sid),
      error = function(e) NULL)
    if (is.null(pred_q)) next

    metric_rows[[win_tag]] <- data.frame(
      sim_idx        = sid, window = win_tag, regime = regime, win_len = win_len,
      true_beta      = act$beta[1],       pred_beta      = prd$beta_pred[1],
      true_R0        = act$R0[1],         pred_R0        = prd$R0_pred[1],
      err_beta       = abs(prd$beta_pred[1] - act$beta[1]),
      err_R0         = abs(prd$R0_pred[1]   - act$R0[1]),
      curve_mae      = mean(abs(act_q$med - pred_q$med)),
      peak_day_act   = which.max(act_q$med), peak_day_pred = which.max(pred_q$med),
      peak_day_err   = abs(which.max(pred_q$med) - which.max(act_q$med)),
      peak_size_act  = max(act_q$med),       peak_size_pred = max(pred_q$med),
      peak_size_err  = abs(max(pred_q$med) - max(act_q$med)),
      stringsAsFactors = FALSE)

    if (win_tag %in% CURVE_WINDOWS) {
      curve_rows[[win_tag]] <- data.frame(
        day = 1:NDAYS, sim_idx = sid, window = win_tag,
        win_len = win_len, regime = regime, obs = obs_full,
        act_lower = act_q$lower, act_med = act_q$med, act_upper = act_q$upper,
        pred_lower = pred_q$lower, pred_med = pred_q$med, pred_upper = pred_q$upper,
        in_window = (1:NDAYS) >= t_start & (1:NDAYS) <= t_end,
        stringsAsFactors = FALSE)
    }
  }

  list(metrics = bind_rows(metric_rows),
       curves  = bind_rows(curve_rows))
}

# Now submit — slurmR will serialize THIS version of eval_one_sim
step1 <- Slurm_lapply(
  X          = common_ids,
  FUN        = eval_one_sim,
  job_name   = "part4_array",
  njobs      = 100,
  overwrite  = TRUE,
  plan       = "submit",
  tmp_path   = path.expand(SCRATCH_DIR),
  sbatch_opt = slurm_opts,
  export     = c("eval_one_sim")
)
# =============================================================================
# STEP 1: Submit array job — one Slurm task per sim_idx
# =============================================================================

cat("Step 1 submitted. Waiting for jobs to finish before collecting...\n")
# slurmR::wait_slurm_job(step1)   # uncomment if you want to block here

# =============================================================================
# STEP 2: Collect results (run after jobs finish)
# =============================================================================

# Re-read the job object by path so this section can be run independently:
#   step1 <- read_slurm_job(path.expand(file.path(SCRATCH_DIR, "part4_array")))

results_raw <- Slurm_collect(
  read_slurm_job(path.expand(file.path(SCRATCH_DIR, "part4_array"))),
  any. = TRUE
)
results_raw <- Filter(Negate(is.null), results_raw)
cat(sprintf("Collected %d non-null results.\n", length(results_raw)))

# =============================================================================
# STEP 3: Merge + plots (runs on the login node after collect)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

perf_df <- bind_rows(lapply(results_raw, `[[`, "metrics")) |>
  mutate(regime = factor(regime, levels = c("early", "mid", "late")))

curves_df <- bind_rows(lapply(results_raw, `[[`, "curves")) |>
  mutate(window = factor(window, levels = CURVE_WINDOWS),
         regime = factor(regime, levels = c("early", "mid", "late")))

cat(sprintf("Metrics: %d rows | Curves: %d rows\n",
            nrow(perf_df), nrow(curves_df)))

# -- Save CSV -----------------------------------------------------------------
out_csv <- path.expand(file.path(DATA_DIR, "part4_overall_performance.csv"))
write.csv(perf_df, out_csv, row.names = FALSE)
cat(sprintf("Saved: %s\n", out_csv))

# -- Console summary ----------------------------------------------------------
cat("\n=== Performance by window (averaged across all sims) ===\n")
perf_df |>
  group_by(window, regime, win_len) |>
  summarise(mae_beta     = mean(err_beta),
            mae_R0       = mean(err_R0),
            curve_mae    = mean(curve_mae),
            peak_day_err = mean(peak_day_err),
            n            = n(),
            .groups      = "drop") |>
  arrange(regime, win_len) |>
  mutate(across(where(is.numeric), round, 3)) |>
  print(n = 40)

# =========================================================================
# Plots
# =========================================================================
regime_colors <- c(early = "#E65100", mid = "#1565C0", late = "#2E7D32")

summary_win <- perf_df |>
  group_by(window, regime, win_len) |>
  summarise(across(c(err_beta, err_R0, curve_mae, peak_day_err, peak_size_err),
                   mean, .names = "mean_{.col}"),
            .groups = "drop")

# Plot 1: error vs window length
metric_long <- summary_win |>
  select(win_len, regime, mean_err_beta, mean_err_R0, mean_curve_mae) |>
  pivot_longer(starts_with("mean_"), names_to = "metric", values_to = "value") |>
  mutate(metric = recode(metric,
                         mean_err_beta  = "MAE beta",
                         mean_err_R0    = "MAE R0",
                         mean_curve_mae = "Curve MAE"))

pW1 <- ggplot(metric_long,
              aes(x = win_len, y = value, color = regime, group = regime)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.5) +
  facet_wrap(~ metric, scales = "free_y", ncol = 3) +
  scale_x_log10(breaks = LENGTHS) +
  scale_color_manual(values = regime_colors) +
  labs(title = sprintf(
    "Part 4: Error vs Window Length (N=%d sims, %d runs each)",
    length(unique(perf_df$sim_idx)), NSIMS_WORKER),
    x = "Window length (days, log scale)", y = "Mean error", color = "Regime") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

ggsave(path.expand(file.path(DATA_DIR, "part4_error_vs_window.png")),
       pW1, width = 14, height = 5, dpi = 150)

# Plot 2: peak errors
peak_long <- summary_win |>
  select(win_len, regime, mean_peak_day_err, mean_peak_size_err) |>
  pivot_longer(starts_with("mean_peak"), names_to = "metric", values_to = "value") |>
  mutate(metric = recode(metric,
                         mean_peak_day_err  = "Peak day error (days)",
                         mean_peak_size_err = "Peak size error"))

pW2 <- ggplot(peak_long,
              aes(x = win_len, y = value, color = regime, group = regime)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.5) +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  scale_x_log10(breaks = LENGTHS) +
  scale_color_manual(values = regime_colors) +
  labs(title = "Part 4: Peak Prediction Error vs Window Length",
       x = "Window length (days, log scale)", y = "Mean error", color = "Regime") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

ggsave(path.expand(file.path(DATA_DIR, "part4_peak_vs_window.png")),
       pW2, width = 10, height = 5, dpi = 150)

# Plot 3: param scatter
windows_show <- CURVE_WINDOWS
sub_df <- perf_df |>
  filter(window %in% windows_show) |>
  mutate(window = factor(window, levels = windows_show))

pW3 <- (ggplot(sub_df, aes(x = true_beta, y = pred_beta, color = regime)) +
          geom_abline(slope = 1, intercept = 0,
                      linetype = "dashed", color = "gray50") +
          geom_point(alpha = 0.25, size = 0.8) +
          facet_wrap(~ window, ncol = 3, scales = "free") +
          scale_color_manual(values = regime_colors) +
          labs(title = "beta: predicted vs actual",
               x = "Actual", y = "Predicted", color = "Regime") +
          theme_bw(base_size = 9) +
          theme(panel.grid.minor = element_blank(),
                strip.text = element_text(size = 8))) /
  (ggplot(sub_df, aes(x = true_R0, y = pred_R0, color = regime)) +
     geom_abline(slope = 1, intercept = 0,
                 linetype = "dashed", color = "gray50") +
     geom_point(alpha = 0.25, size = 0.8) +
     facet_wrap(~ window, ncol = 3, scales = "free") +
     scale_color_manual(values = regime_colors) +
     labs(title = "R0: predicted vs actual",
          x = "Actual", y = "Predicted", color = "Regime") +
     theme_bw(base_size = 9) +
     theme(panel.grid.minor = element_blank(),
           strip.text = element_text(size = 8)))

ggsave(path.expand(file.path(DATA_DIR, "part4_param_scatter_by_window.png")),
       pW3, width = 13, height = 10, dpi = 150)

# Plot 4: R0 error violin
pW4 <- ggplot(perf_df,
              aes(x = factor(win_len), y = err_R0, fill = regime)) +
  geom_violin(alpha = 0.6, scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.08, outlier.size = 0.3, alpha = 0.8,
               position = position_dodge(0.9)) +
  facet_wrap(~ regime, ncol = 3) +
  scale_fill_manual(values = regime_colors) +
  labs(title = "Part 4: R0 Error Distribution by Window Length and Regime",
       x = "Window length (days)",
       y = "|Predicted R0 - Actual R0|",
       fill = "Regime") +
  theme_bw(base_size = 11) +
  theme(legend.position  = "none",
        panel.grid.minor = element_blank(),
        plot.title       = element_text(face = "bold"))

ggsave(path.expand(file.path(DATA_DIR, "part4_R0_error_violin.png")),
       pW4, width = 13, height = 5, dpi = 150)

# Plot 5: Aggregated incidence curves
cat("\nBuilding Plot 5...\n")

agg_curves <- curves_df |>
  group_by(window, win_len, regime, day) |>
  summarise(
    obs_med  = median(obs),
    obs_lo   = quantile(obs,      0.05),
    obs_hi   = quantile(obs,      0.95),
    act_med  = median(act_med),
    act_lo   = quantile(act_med,  0.05),
    act_hi   = quantile(act_med,  0.95),
    pred_med = median(pred_med),
    pred_lo  = quantile(pred_med, 0.05),
    pred_hi  = quantile(pred_med, 0.95),
    .groups  = "drop"
  )

rect_agg <- curves_df |>
  group_by(window) |>
  summarise(xmin = min(day[in_window]),
            xmax = max(day[in_window]),
            .groups = "drop")

pW5 <- ggplot(agg_curves, aes(x = day)) +
  geom_rect(data = rect_agg,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "gold", alpha = 0.20, inherit.aes = FALSE) +
  geom_ribbon(aes(ymin = obs_lo,  ymax = obs_hi),
              fill = "gray40", alpha = 0.12) +
  geom_line(aes(y = obs_med, color = "Observed (ABM)"),
            linewidth = 0.6, alpha = 0.7) +
  geom_ribbon(aes(ymin = act_lo,  ymax = act_hi),
              fill = "#1976D2", alpha = 0.18) +
  geom_line(aes(y = act_med,  color = "SEIR - Actual Params"),
            linewidth = 1.0) +
  geom_ribbon(aes(ymin = pred_lo, ymax = pred_hi),
              fill = "#D32F2F", alpha = 0.18) +
  geom_line(aes(y = pred_med, color = "SEIR - BiLSTM Predicted"),
            linewidth = 1.0) +
  facet_wrap(~ window, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c(
    "Observed (ABM)"          = "gray30",
    "SEIR - Actual Params"    = "#1976D2",
    "SEIR - BiLSTM Predicted" = "#D32F2F")) +
  labs(
    title = sprintf(
      "Part 4 Plot 5: Aggregated Incidence across ALL %d Val Sims (%d runs each)",
      length(unique(curves_df$sim_idx)), NSIMS_WORKER),
    subtitle = paste0(
      "Lines = median across sims | Bands = 5th-95th pct across sims",
      " | Gold = observation window"),
    x = "Day", y = "Daily Incidence", color = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position  = "bottom",
        legend.key.width = unit(1.5, "cm"),
        strip.text       = element_text(size = 9),
        panel.grid.minor = element_blank(),
        plot.title       = element_text(face = "bold", size = 11))

ggsave(path.expand(file.path(DATA_DIR, "part4_aggregated_curves.png")),
       pW5, width = 14, height = 9, dpi = 150)

cat("\nPart 4 complete. Saved:\n")
cat("  part4_overall_performance.csv\n")
cat("  part4_error_vs_window.png\n")
cat("  part4_peak_vs_window.png\n")
cat("  part4_param_scatter_by_window.png\n")
cat("  part4_R0_error_violin.png\n")
cat("  part4_aggregated_curves.png\n")
