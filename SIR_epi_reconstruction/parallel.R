# =============================================================================
#  PART 4 (SIR) — Overall performance on ALL test data x ALL 18 windows (Slurm)
#
#  Usage (on login node):
#    Rscript part4_sir_slurm_submit.R
# =============================================================================

library(slurmR)
library(dplyr)

# -- Config -------------------------------------------------------------------
DATA_DIR     <- "~/sima/epiworldRcalibrate_SEIR/SIR_epi_reconstruction"
SCRATCH_DIR  <- file.path(DATA_DIR, "part4_scratch")
NSIMS_WORKER <- 100L
NDAYS        <- 365L

CURVE_WINDOWS <- c("early_015d", "early_060d", "early_180d",
                   "mid_090d",   "mid_180d",   "late_180d",  "late_365d")

LENGTHS <- c(15, 30, 60, 90, 180, 365)

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

dir.create(path.expand(SCRATCH_DIR), recursive = TRUE, showWarnings = FALSE)

# -- Load IDs -----------------------------------------------------------------
actual    <- read.csv(path.expand(file.path(DATA_DIR, "test_actual_parameters.csv")))
preds_all <- read.csv(path.expand(file.path(DATA_DIR, "test_bilstm_predictions.csv")))
common_ids <- intersect(unique(preds_all$sim_idx), actual$sim_idx)
cat(sprintf("Submitting array over %d sim IDs\n", length(common_ids)))

# -- Slurm options ------------------------------------------------------------
slurm_opts <- list(
  account         = "vegayon-np",
  partition       = "vegayon-np",
  time            = "10:00:00",
  `mem-per-cpu`   = "4G",
  `cpus-per-task` = 1
)

# =============================================================================
# Core worker function — one sim_idx at a time
# =============================================================================

eval_one_sim_sir <- function(sid) {

  # FIX: slurmR wraps each X element in a list; coerce to scalar
  sid <- as.integer(sid)

  suppressPackageStartupMessages({
    library(epiworldR)
    library(dplyr)
    library(data.table)
  })

  # Re-declare all constants (never rely on closure serialisation)
  DATA_DIR      <- "~/sima/epiworldRcalibrate_SEIR/SIR_epi_reconstruction"
  NDAYS         <- 365L
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

  # -- Load data per worker ---------------------------------------------------
  actual    <- read.csv(path.expand(file.path(DATA_DIR, "test_actual_parameters.csv")))
  preds_all <- read.csv(path.expand(file.path(DATA_DIR, "test_bilstm_predictions.csv")))

  np      <- reticulate::import("numpy")
  inc_raw <- as.matrix(
    np$load(path.expand(file.path(DATA_DIR, "test_incidence_raw.npy"))))

  # -- run_multiple helper (SIR: S -> I transitions) -------------------------
  .run_multi_worker <- function(n, prevalence, ptran, crate, recov, seed,
                                ndays = NDAYS, nsims = NSIMS_WORKER) {
    n_int <- max(as.integer(round(n)), 10L)
    prev  <- max(min(prevalence, 1.0), 1.0 / n_int)

    model <- ModelSIRCONN(
      name              = "sim",
      n                 = n_int,
      prevalence        = prev,
      contact_rate      = max(crate, 1e-6),
      transmission_rate = max(ptran, 1e-6),
      recovery_rate     = max(recov, 1e-6)
    )
    verbose_off(model)

    saver <- make_saver("transition")
    run_multiple(model, ndays = ndays, nsims = nsims,
                 saver = saver, nthreads = 1L)

    res <- run_multiple_get_results(model, nthreads = 1L,
                                    freader = data.table::fread)

    res$transition |>
      filter(from == "Susceptible", to == "Infected", date > 0) |>
      group_by(date) |>
      summarise(
        lower = quantile(counts, 0.025),
        med   = quantile(counts, 0.500),
        upper = quantile(counts, 0.975),
        .groups = "drop"
      ) |>
      right_join(data.frame(date = 1:ndays), by = "date") |>
      arrange(date) |>
      mutate(across(c(lower, med, upper), \(x) replace_na(x, 0)))
  }

  # -- Lookup this sim --------------------------------------------------------
  act     <- actual[actual$sim_idx == sid, ]
  mat_row <- which(actual$sim_idx == sid)
  if (nrow(act) == 0 || length(mat_row) == 0) return(NULL)

  n          <- act$n[1]
  recov      <- act$recov[1]
  prevalence <- act$prevalence[1]
  obs_full   <- as.numeric(inc_raw[mat_row, ])

  act_q <- tryCatch(
    .run_multi_worker(n, prevalence,
                      ptran = act$ptran[1], crate = act$crate[1],
                      recov = recov, seed = sid),
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

    prd <- preds_all[preds_all$sim_idx == sid &
                       preds_all$window  == win_tag, ]
    if (nrow(prd) == 0) next

    pred_q <- tryCatch(
      .run_multi_worker(n, prevalence,
                        ptran = prd$ptran_pred[1], crate = prd$crate_pred[1],
                        recov = recov, seed = sid),
      error = function(e) NULL)
    if (is.null(pred_q)) next

    metric_rows[[win_tag]] <- data.frame(
      sim_idx        = sid,
      window         = win_tag,
      regime         = regime,
      win_len        = win_len,
      true_ptran     = act$ptran[1],
      pred_ptran     = prd$ptran_pred[1],
      true_crate     = act$crate[1],
      pred_crate     = prd$crate_pred[1],
      true_R0        = act$R0[1],
      pred_R0        = prd$R0_pred[1],
      err_ptran      = abs(prd$ptran_pred[1] - act$ptran[1]),
      err_crate      = abs(prd$crate_pred[1] - act$crate[1]),
      err_R0         = abs(prd$R0_pred[1]    - act$R0[1]),
      curve_mae      = mean(abs(act_q$med    - pred_q$med)),
      peak_day_act   = which.max(act_q$med),
      peak_day_pred  = which.max(pred_q$med),
      peak_day_err   = abs(which.max(pred_q$med) - which.max(act_q$med)),
      peak_size_act  = max(act_q$med),
      peak_size_pred = max(pred_q$med),
      peak_size_err  = abs(max(pred_q$med) - max(act_q$med)),
      stringsAsFactors = FALSE
    )

    if (win_tag %in% CURVE_WINDOWS) {
      curve_rows[[win_tag]] <- data.frame(
        day        = 1:NDAYS,
        sim_idx    = sid,
        window     = win_tag,
        win_len    = win_len,
        regime     = regime,
        obs        = obs_full,
        act_lower  = act_q$lower,
        act_med    = act_q$med,
        act_upper  = act_q$upper,
        pred_lower = pred_q$lower,
        pred_med   = pred_q$med,
        pred_upper = pred_q$upper,
        in_window  = (1:NDAYS) >= t_start & (1:NDAYS) <= t_end,
        stringsAsFactors = FALSE
      )
    }
  }

  list(metrics = bind_rows(metric_rows),
       curves  = bind_rows(curve_rows))
}

# =============================================================================
# STEP 1: Submit array job
# =============================================================================

cat("Submitting SLURM array job for Part 4 (SIR)...\n")

step1 <- Slurm_lapply(
  X          = common_ids,
  FUN        = eval_one_sim_sir,
  job_name   = "part4_sir_array",
  njobs      = 100,
  overwrite  = TRUE,
  plan       = "submit",
  tmp_path   = path.expand(SCRATCH_DIR),
  sbatch_opt = slurm_opts,
  export     = c("eval_one_sim_sir")
)

cat("Step 1 submitted. Check status with: squeue -u $USER\n")
cat("Once all jobs finish, run the collect section below.\n")

# =============================================================================
# STEP 2: Collect
#   To re-run independently:
#   step1 <- read_slurm_job(path.expand(file.path(SCRATCH_DIR, "part4_sir_array")))
# =============================================================================

results_raw <- Slurm_collect(
  read_slurm_job(path.expand(file.path(SCRATCH_DIR, "part4_sir_array"))),
  any. = TRUE
)
results_raw <- Filter(Negate(is.null), results_raw)

# Each chunk task returns a list of per-sid results; flatten one level
results_flat <- Filter(Negate(is.null), unlist(results_raw, recursive = FALSE))
cat(sprintf("Collected %d non-null sim results.\n", length(results_flat)))

# =============================================================================
# STEP 3: Merge + plots
# =============================================================================

suppressPackageStartupMessages({
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

perf_df <- bind_rows(lapply(results_flat, `[[`, "metrics")) |>
  mutate(regime = factor(regime, levels = c("early", "mid", "late")))

curves_df <- bind_rows(lapply(results_flat, `[[`, "curves")) |>
  mutate(window = factor(window, levels = CURVE_WINDOWS),
         regime = factor(regime, levels = c("early", "mid", "late")))

cat(sprintf("Metrics: %d rows | Curves: %d rows\n",
            nrow(perf_df), nrow(curves_df)))

# -- Save CSV -----------------------------------------------------------------
out_csv <- path.expand(file.path(DATA_DIR, "part4_sir_overall_performance.csv"))
write.csv(perf_df, out_csv, row.names = FALSE)
cat(sprintf("Saved: %s\n", out_csv))

# -- Console summary ----------------------------------------------------------
cat("\n=== Performance by window (averaged across all sims) ===\n")
perf_df |>
  group_by(window, regime, win_len) |>
  summarise(mae_ptran    = mean(err_ptran),
            mae_crate    = mean(err_crate),
            mae_R0       = mean(err_R0),
            curve_mae    = mean(curve_mae),
            peak_day_err = mean(peak_day_err),
            n            = n(),
            .groups      = "drop") |>
  arrange(regime, win_len) |>
  mutate(across(where(is.numeric), round, 3)) |>
  print(n = 40)

# =========================================================================
# Plots 1-4 (summary plots)
# =========================================================================
regime_colors <- c(early = "#E65100", mid = "#1565C0", late = "#2E7D32")

summary_win <- perf_df |>
  group_by(window, regime, win_len) |>
  summarise(across(c(err_ptran, err_crate, err_R0,
                     curve_mae, peak_day_err, peak_size_err),
                   mean, .names = "mean_{.col}"),
            .groups = "drop")

# Plot 1: error vs window length
metric_long <- summary_win |>
  select(win_len, regime, mean_err_ptran, mean_err_crate,
         mean_err_R0, mean_curve_mae) |>
  pivot_longer(starts_with("mean_"), names_to = "metric", values_to = "value") |>
  mutate(metric = recode(metric,
                         mean_err_ptran = "MAE ptran",
                         mean_err_crate = "MAE crate",
                         mean_err_R0    = "MAE R0",
                         mean_curve_mae = "Curve MAE"))

pW1 <- ggplot(metric_long,
              aes(x = win_len, y = value, color = regime, group = regime)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.5) +
  facet_wrap(~ metric, scales = "free_y", ncol = 4) +
  scale_x_log10(breaks = LENGTHS) +
  scale_color_manual(values = regime_colors) +
  labs(title = sprintf(
    "Part 4 (SIR): Error vs Window Length (N=%d sims, %d runs each)",
    length(unique(perf_df$sim_idx)), NSIMS_WORKER),
    x = "Window length (days, log scale)", y = "Mean error", color = "Regime") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

ggsave(path.expand(file.path(DATA_DIR, "part4_sir_error_vs_window.png")),
       pW1, width = 16, height = 5, dpi = 150)

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
  labs(title = "Part 4 (SIR): Peak Prediction Error vs Window Length",
       x = "Window length (days, log scale)", y = "Mean error", color = "Regime") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

ggsave(path.expand(file.path(DATA_DIR, "part4_sir_peak_vs_window.png")),
       pW2, width = 10, height = 5, dpi = 150)

# Plot 3: param scatter (ptran + crate + R0)
sub_df <- perf_df |>
  filter(window %in% CURVE_WINDOWS) |>
  mutate(window = factor(window, levels = CURVE_WINDOWS))

pW3 <- (ggplot(sub_df, aes(x = true_ptran, y = pred_ptran, color = regime)) +
          geom_abline(slope = 1, intercept = 0,
                      linetype = "dashed", color = "gray50") +
          geom_point(alpha = 0.25, size = 0.8) +
          facet_wrap(~ window, ncol = 3, scales = "free") +
          scale_color_manual(values = regime_colors) +
          labs(title = "ptran: predicted vs actual",
               x = "Actual", y = "Predicted", color = "Regime") +
          theme_bw(base_size = 9) +
          theme(panel.grid.minor = element_blank(),
                strip.text = element_text(size = 8))) /
  (ggplot(sub_df, aes(x = true_crate, y = pred_crate, color = regime)) +
     geom_abline(slope = 1, intercept = 0,
                 linetype = "dashed", color = "gray50") +
     geom_point(alpha = 0.25, size = 0.8) +
     facet_wrap(~ window, ncol = 3, scales = "free") +
     scale_color_manual(values = regime_colors) +
     labs(title = "crate: predicted vs actual",
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

ggsave(path.expand(file.path(DATA_DIR, "part4_sir_param_scatter_by_window.png")),
       pW3, width = 13, height = 14, dpi = 150)

# Plot 4: R0 error violin
pW4 <- ggplot(perf_df,
              aes(x = factor(win_len), y = err_R0, fill = regime)) +
  geom_violin(alpha = 0.6, scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.08, outlier.size = 0.3, alpha = 0.8,
               position = position_dodge(0.9)) +
  facet_wrap(~ regime, ncol = 3) +
  scale_fill_manual(values = regime_colors) +
  labs(title = "Part 4 (SIR): R0 Error Distribution by Window Length and Regime",
       x = "Window length (days)",
       y = "|Predicted R0 - Actual R0|",
       fill = "Regime") +
  theme_bw(base_size = 11) +
  theme(legend.position  = "none",
        panel.grid.minor = element_blank(),
        plot.title       = element_text(face = "bold"))

ggsave(path.expand(file.path(DATA_DIR, "part4_sir_R0_error_violin.png")),
       pW4, width = 13, height = 5, dpi = 150)

# =========================================================================
# Plot 5: Aggregated incidence curves (median + bands across ALL sims)
# =========================================================================
cat("\nBuilding Plot 5 (aggregated curves)...\n")

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
  geom_line(aes(y = act_med,  color = "SIR - Actual Params"),
            linewidth = 1.0) +
  geom_ribbon(aes(ymin = pred_lo, ymax = pred_hi),
              fill = "#D32F2F", alpha = 0.18) +
  geom_line(aes(y = pred_med, color = "SIR - BiLSTM Predicted"),
            linewidth = 1.0) +
  facet_wrap(~ window, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c(
    "Observed (ABM)"         = "gray30",
    "SIR - Actual Params"    = "#1976D2",
    "SIR - BiLSTM Predicted" = "#D32F2F")) +
  labs(
    title = sprintf(
      "Part 4 (SIR) Plot 5: Aggregated Incidence across ALL %d Val Sims (%d runs each)",
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

ggsave(path.expand(file.path(DATA_DIR, "part4_sir_aggregated_curves.png")),
       pW5, width = 14, height = 9, dpi = 150)

# =========================================================================
# Plot 6: Worst-performing individual sim panels per window
#         Top 20 highest curve_mae, one PNG per CURVE_WINDOW
# =========================================================================
cat("\nBuilding Plot 6 (worst-performing individual curves per window)...\n")

N_WORST <- 20

for (wt in CURVE_WINDOWS) {

  # top N_WORST sims by curve_mae for this window
  worst_ids <- perf_df |>
    filter(window == wt) |>
    arrange(desc(curve_mae)) |>
    slice_head(n = N_WORST) |>
    pull(sim_idx)

  plot_data <- curves_df |>
    filter(window == wt, sim_idx %in% worst_ids) |>
    left_join(
      perf_df |>
        filter(window == wt) |>
        select(sim_idx, true_ptran, pred_ptran, true_crate, pred_crate,
               true_R0, pred_R0, curve_mae),
      by = "sim_idx"
    ) |>
    mutate(
      panel = sprintf(
        "sim %d | MAE=%.1f\nptran: %.3f->%.3f | R0: %.2f->%.2f",
        sim_idx, curve_mae, true_ptran, pred_ptran, true_R0, pred_R0),
      panel = factor(panel, levels = unique(panel[order(-curve_mae)]))
    )

  # per-panel window shading
  rect_df <- plot_data |>
    group_by(panel) |>
    summarise(xmin = min(day[in_window]),
              xmax = max(day[in_window]),
              .groups = "drop")

  p6 <- ggplot(plot_data, aes(x = day)) +
    geom_rect(data = rect_df,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = "gold", alpha = 0.15, inherit.aes = FALSE) +
    geom_ribbon(aes(ymin = act_lower,  ymax = act_upper),
                fill = "#1976D2", alpha = 0.15) +
    geom_line(aes(y = act_med,  color = "SIR - Actual Params"),
              linewidth = 0.7) +
    geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper),
                fill = "#D32F2F", alpha = 0.15) +
    geom_line(aes(y = pred_med, color = "SIR - BiLSTM Predicted"),
              linewidth = 0.7, linetype = "dashed") +
    geom_line(aes(y = obs, color = "Observed (ABM)"),
              linewidth = 0.5, linetype = "dotdash") +
    facet_wrap(~ panel, scales = "free_y", ncol = 4) +
    scale_color_manual(values = c(
      "Observed (ABM)"         = "black",
      "SIR - Actual Params"    = "#1976D2",
      "SIR - BiLSTM Predicted" = "#D32F2F")) +
    labs(
      title    = sprintf("Part 4 (SIR) Plot 6: Worst %d Sims — Window: %s", N_WORST, wt),
      subtitle = "Ranked by highest curve MAE | Gold = observation window | Dashed = BiLSTM predicted",
      x = "Day", y = "Daily Incidence", color = NULL) +
    theme_bw(base_size = 9) +
    theme(legend.position  = "bottom",
          legend.key.width = unit(1.5, "cm"),
          strip.text       = element_text(size = 6.5),
          panel.grid.minor = element_blank(),
          plot.title       = element_text(face = "bold"))

  out_file <- path.expand(file.path(DATA_DIR,
                                    sprintf("part4_sir_worst_curves_%s.png", wt)))
  ggsave(out_file, p6, width = 18, height = 14, dpi = 150)
  cat(sprintf("  Saved: part4_sir_worst_curves_%s.png\n", wt))
}

cat("\nPart 4 (SIR) complete. Saved:\n")
cat("  part4_sir_overall_performance.csv\n")
cat("  part4_sir_error_vs_window.png\n")
cat("  part4_sir_peak_vs_window.png\n")
cat("  part4_sir_param_scatter_by_window.png\n")
cat("  part4_sir_R0_error_violin.png\n")
cat("  part4_sir_aggregated_curves.png\n")
cat("  part4_sir_worst_curves_<window>.png  (7 files, one per CURVE_WINDOW)\n")
