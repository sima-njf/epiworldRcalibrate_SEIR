# =============================================================================
#  SIR Calibration Validation — with run_multiple (uncertainty bands)
#
#  Part 1 — Test set: run ModelSIRCONN with actual vs BiLSTM-predicted
#            parameters, compare incidence curves + 95% CI, per-day MAE
#
#  Part 2 — Random parameters pipeline:
#            sample random params → simulate incidence → calibrate_sir()
#            → simulate with predicted params → plot comparison with CI
# =============================================================================

library(epiworldR)
library(ggplot2)
library(dplyr)
library(tidyr)

# ── Paths ─────────────────────────────────────────────────────────────────────
MODEL_DIR <- "~/sima/epiworldRcalibrate_SEIR/SIR_epi_reconstruction/model"
DATA_DIR  <- "~/sima/epiworldRcalibrate_SEIR/SIR_epi_reconstruction"

NDAYS    <- 365
N_SAMPLE <- 20
WINDOW   <- "late_365d"
NSIMS    <- 500    # simulations per ModelSIRCONN run (increase for smoother CIs)
NTHREADS <- 4      # parallel threads

# ── Load BiLSTM model ─────────────────────────────────────────────────────────
init_bilstm_model(MODEL_DIR)

# =============================================================================
# Helper: run ModelSIRCONN NSIMS times, return median + 95% CI per day
# =============================================================================
run_sir_multi <- function(n, prevalence, ptran, crate, recov,
                          ndays = NDAYS, nsims = NSIMS, nthreads = NTHREADS,
                          seed = 1L) {
  n    <- max(as.integer(round(n)), 10L)
  prev <- max(min(prevalence, 1.0), 1.0 / n)

  model <- ModelSIRCONN(
    name              = "sim",
    n                 = n,
    prevalence        = prev,
    contact_rate      = max(crate, 1e-6),
    transmission_rate = max(ptran, 1e-6),
    recovery_rate     = max(recov, 1e-6)
  )
  verbose_off(model)

  saver <- make_saver("transition")
  run_multiple(model, ndays = ndays, nsims = nsims,
               saver = saver, nthreads = nthreads)

  res <- run_multiple_get_results(model, nthreads = nthreads,
                                  freader = data.table::fread)

  res$transition |>
    filter(from == "Susceptible", to == "Infected") |>
    group_by(date) |>
    summarise(
      lower = quantile(counts, 0.025),
      med   = quantile(counts, 0.500),
      upper = quantile(counts, 0.975),
      .groups = "drop"
    ) |>
    right_join(data.frame(date = 0:(ndays - 1)), by = "date") |>
    arrange(date) |>
    mutate(across(c(lower, med, upper), \(x) replace_na(x, 0)))
}

# =============================================================================
# PART 1 — Test set: actual vs predicted incidence
# =============================================================================

cat("\n-- Part 1: loading test data --\n")
actual    <- read.csv(path.expand(file.path(DATA_DIR, "test_actual_parameters.csv")))
preds_all <- read.csv(path.expand(file.path(DATA_DIR, "test_bilstm_predictions.csv")))

np      <- reticulate::import("numpy")
inc_raw <- np$load(path.expand(file.path(DATA_DIR, "test_incidence_raw.npy")))

preds <- preds_all |> filter(window == WINDOW)
cat(sprintf("  Val sims: %d  |  Window: %s  |  NSIMS per run: %d\n",
            nrow(actual), WINDOW, NSIMS))

# ── Sample test simulations ───────────────────────────────────────────────────
set.seed(42)
sample_ids <- sample(intersect(preds$sim_idx, actual$sim_idx), N_SAMPLE)
cat(sprintf("Running %d sims x 2 models x %d runs each...\n", N_SAMPLE, NSIMS))

curves_list  <- vector("list", N_SAMPLE)
mae_day_list <- vector("list", N_SAMPLE)

for (k in seq_along(sample_ids)) {
  sid <- sample_ids[k]
  cat(sprintf("  [%d/%d] sim_idx = %d\n", k, N_SAMPLE, sid))

  act <- actual |> filter(sim_idx == sid)
  prd <- preds  |> filter(sim_idx == sid)

  n          <- act$n[1]
  recov      <- act$recov[1]
  prevalence <- act$prevalence[1]

  mat_row <- which(actual$sim_idx == sid)
  obs_inc <- as.numeric(inc_raw[mat_row, ])

  # ModelSIRCONN x NSIMS with ACTUAL parameters
  act_q <- run_sir_multi(
    n = n, prevalence = prevalence,
    ptran = act$ptran[1], crate = act$crate[1], recov = recov, seed = sid)

  # ModelSIRCONN x NSIMS with BiLSTM PREDICTED parameters
  pred_q <- run_sir_multi(
    n = n, prevalence = prevalence,
    ptran = prd$ptran_pred[1], crate = prd$crate_pred[1], recov = recov,
    seed = sid)

  curves_list[[k]] <- data.frame(
    day        = 1:NDAYS,
    sim_idx    = sid,
    panel      = sprintf(
      "sim %d\nptran: %.3f->%.3f | crate: %.2f->%.2f | R0: %.2f->%.2f",
      sid,
      act$ptran[1], prd$ptran_pred[1],
      act$crate[1], prd$crate_pred[1],
      act$R0[1],    prd$R0_pred[1]
    ),
    obs        = obs_inc,
    act_lower  = act_q$lower,
    act_med    = act_q$med,
    act_upper  = act_q$upper,
    pred_lower = pred_q$lower,
    pred_med   = pred_q$med,
    pred_upper = pred_q$upper
  )

  mae_day_list[[k]] <- data.frame(
    day     = 1:NDAYS,
    sim_idx = sid,
    mae     = abs(act_q$med - pred_q$med)
  )
}

curves_df  <- bind_rows(curves_list)
mae_day_df <- bind_rows(mae_day_list)

# ── Plot 1: incidence curves with CI ribbons ──────────────────────────────────
p1 <- ggplot(curves_df, aes(x = day)) +
  geom_ribbon(aes(ymin = act_lower,  ymax = act_upper),
              fill = "#1976D2", alpha = 0.15) +
  geom_line(aes(y = act_med,  color = "SIR - Actual Params"), linewidth = 0.7) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper),
              fill = "#D32F2F", alpha = 0.15) +
  geom_line(aes(y = pred_med, color = "SIR - BiLSTM Predicted"),
            linewidth = 0.7, linetype = "dashed") +
  geom_line(aes(y = obs,      color = "Observed (ABM)"),
            linewidth = 0.6, linetype = "dotdash") +
  facet_wrap(~ panel, scales = "free_y", ncol = 4) +
  scale_color_manual(values = c(
    "Observed (ABM)"        = "black",
    "SIR - Actual Params"   = "#1976D2",
    "SIR - BiLSTM Predicted"= "#D32F2F")) +
  labs(title    = paste0("SIR Incidence - Actual vs BiLSTM Predicted (", WINDOW, ")"),
       subtitle = paste0("Ribbons = 95% CI across ", NSIMS, " stochastic runs"),
       x = "Day", y = "Daily Incidence", color = NULL) +
  theme_bw(base_size = 9) +
  theme(legend.position  = "bottom",
        legend.key.width = unit(1.8, "cm"),
        strip.text       = element_text(size = 6.5),
        panel.grid.minor = element_blank(),
        plot.title       = element_text(face = "bold"))

ggsave(path.expand(file.path(DATA_DIR, "part1_incidence_curves.png")),
       p1, width = 18, height = 14, dpi = 150)
print(p1)

# ── Plot 2: per-day MAE (on medians) ─────────────────────────────────────────
mae_avg <- mae_day_df |>
  group_by(day) |>
  summarise(mean_mae = mean(mae), sd_mae = sd(mae), .groups = "drop")

p2 <- ggplot(mae_avg, aes(x = day)) +
  geom_ribbon(aes(ymin = pmax(mean_mae - sd_mae, 0),
                  ymax = mean_mae + sd_mae),
              fill = "#D32F2F", alpha = 0.15) +
  geom_line(aes(y = mean_mae), color = "#D32F2F", linewidth = 0.9) +
  labs(title    = "Per-day MAE between Actual-Param and Predicted-Param SIR (medians)",
       subtitle = sprintf("Mean +/- SD across %d sampled simulations", N_SAMPLE),
       x = "Day", y = "MAE (incidence)") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

ggsave(path.expand(file.path(DATA_DIR, "part1_per_day_mae.png")),
       p2, width = 10, height = 5, dpi = 150)
print(p2)

cat("\n-- Per-day MAE summary --\n")
mae_avg |>
  summarise(overall_mean = mean(mean_mae),
            peak_day     = day[which.max(mean_mae)],
            peak_mae     = max(mean_mae)) |>
  with({
    cat(sprintf("  Overall mean MAE : %.2f\n", overall_mean))
    cat(sprintf("  Worst day        : day %d  (MAE = %.2f)\n", peak_day, peak_mae))
  })

# =============================================================================
# PART 2 — Random parameters pipeline
# =============================================================================

cat("\n-- Part 2: random parameters pipeline --\n")

N_RANDOM  <- 15
NDAYS_RAN <- 365
set.seed(123)

random_params <- data.frame(
  n          = sample(5000:10000, N_RANDOM, replace = TRUE),
  prevalence = runif(N_RANDOM, min = 0.03, max = 0.35),
  ptran      = runif(N_RANDOM, min = 0.01, max = 0.15),
  crate      = runif(N_RANDOM, min = 5,    max = 15),
  recov      = runif(N_RANDOM, min = 0.07, max = 0.25)
)
random_params$R0 <- random_params$ptran * random_params$crate / random_params$recov

cat("Random parameters sampled:\n")
print(round(random_params, 4))

random_curves <- vector("list", N_RANDOM)

for (k in seq_len(N_RANDOM)) {
  p <- random_params[k, ]
  cat(sprintf("\n  [%d/%d]  n=%d  ptran=%.3f  crate=%.2f  recov=%.3f  R0=%.2f\n",
              k, N_RANDOM, p$n, p$ptran, p$crate, p$recov, p$R0))

  # Step 1: simulate with TRUE params x NSIMS
  true_q <- run_sir_multi(
    n = p$n, prevalence = p$prevalence,
    ptran = p$ptran, crate = p$crate, recov = p$recov,
    ndays = NDAYS_RAN, seed = k)

  # Step 2: feed the MEDIAN into BiLSTM
  predicted <- calibrate_sir(
    daily_cases     = true_q$med,
    population_size = p$n,
    recovery_rate   = p$recov
  )
  cat(sprintf("    True      : ptran=%.4f  crate=%.4f  R0=%.4f\n",
              p$ptran, p$crate, p$R0))
  cat(sprintf("    Predicted : ptran=%.4f  crate=%.4f  R0=%.4f\n",
              predicted["ptran"], predicted["crate"], predicted["R0"]))

  # Step 3: simulate with PREDICTED params x NSIMS
  pred_q <- run_sir_multi(
    n = p$n, prevalence = p$prevalence,
    ptran = predicted["ptran"], crate = predicted["crate"], recov = p$recov,
    ndays = NDAYS_RAN, seed = k)

  random_curves[[k]] <- data.frame(
    day        = 1:NDAYS_RAN,
    panel      = sprintf(
      "sim %d | R0: %.2f->%.2f\nptran: %.3f->%.3f | crate: %.2f->%.2f",
      k, p$R0, predicted["R0"],
      p$ptran, predicted["ptran"],
      p$crate, predicted["crate"]
    ),
    true_lower = true_q$lower,
    true_med   = true_q$med,
    true_upper = true_q$upper,
    pred_lower = pred_q$lower,
    pred_med   = pred_q$med,
    pred_upper = pred_q$upper
  )
}

random_df <- bind_rows(random_curves)

# ── Plot 3: random pipeline comparison with CI ribbons ────────────────────────
p3 <- ggplot(random_df, aes(x = day)) +
  geom_ribbon(aes(ymin = true_lower, ymax = true_upper),
              fill = "#1976D2", alpha = 0.15) +
  geom_line(aes(y = true_med, color = "SIR - True (Random) Params"),
            linewidth = 0.8) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper),
              fill = "#D32F2F", alpha = 0.15) +
  geom_line(aes(y = pred_med, color = "SIR - BiLSTM Predicted"),
            linewidth = 0.8, linetype = "dashed") +
  facet_wrap(~ panel, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c(
    "SIR - True (Random) Params" = "#1976D2",
    "SIR - BiLSTM Predicted"     = "#D32F2F")) +
  labs(title    = "Part 2 - Random Parameters: True vs BiLSTM-Predicted SIR Curves",
       subtitle = paste0("Ribbons = 95% CI across ", NSIMS,
                         " stochastic runs | BiLSTM input = median of true runs"),
       x = "Day", y = "Daily Incidence", color = NULL) +
  theme_bw(base_size = 9) +
  theme(legend.position  = "bottom",
        legend.key.width = unit(1.8, "cm"),
        strip.text       = element_text(size = 7),
        panel.grid.minor = element_blank(),
        plot.title       = element_text(face = "bold"))

ggsave(path.expand(file.path(DATA_DIR, "part2_random_pipeline.png")),
       p3, width = 15, height = 16, dpi = 150)
print(p3)

# ── Part 2 summary ────────────────────────────────────────────────────────────
cat("\n-- Part 2 summary --\n")
random_df |>
  group_by(panel) |>
  summarise(
    mae_medians    = mean(abs(true_med - pred_med)),
    peak_true      = max(true_med),
    peak_predicted = max(pred_med),
    peak_day_true  = which.max(true_med),
    peak_day_pred  = which.max(pred_med),
    .groups = "drop"
  ) |>
  mutate(peak_day_err = abs(peak_day_true - peak_day_pred)) |>
  print()

cat("\nSaved:\n")
cat("  part1_incidence_curves.png\n")
cat("  part1_per_day_mae.png\n")
cat("  part2_random_pipeline.png\n")

# =============================================================================
# PART 3 — Window analysis: inference across all 18 windows
#           (6 lengths x 3 regimes: early / mid / late)
#
# For each window we:
#   1. Slice the observed incidence to that window
#   2. Run calibrate_sir() on the slice
#   3. Run ModelSIRCONN x NSIMS with predicted params over the FULL 365 days
#   4. Compare recovered full curve to actual-param SIR (median)
#
# Plots produced:
#   A. Parameter error (R0, ptran, crate) vs window length — by regime
#   B. Recovered-curve MAE vs window length — by regime
#   C. Example curves for one sim across all window lengths — early / mid / late
#   D. "Can we predict the peak?" — predicted peak day & size vs actual, by window
# =============================================================================

cat("\n-- Part 3: window analysis --\n")

# ── Window definitions (match training exactly) ───────────────────────────────
LENGTHS <- c(15, 30, 60, 90, 180, 365)
make_windows <- function(t_max = 365, lengths = LENGTHS) {
  out <- list()
  for (L in lengths) {
    if (L > t_max) next
    out[[sprintf("early_%03dd", L)]] <- c(start = 0,                   len = L)
    mid_s <- max(0L, as.integer((t_max - L) / 2))
    out[[sprintf("mid_%03dd",   L)]] <- c(start = mid_s,               len = L)
    late_s <- max(0L, t_max - L)
    out[[sprintf("late_%03dd",  L)]] <- c(start = late_s,              len = L)
  }
  out
}
TEST_WINDOWS <- make_windows()

# ── Sample a smaller set for Part 3 (expensive: 18 windows x N_P3 sims x NSIMS runs)
N_P3 <- 10    # number of val sims to use in Part 3
set.seed(7)
p3_ids <- sample(intersect(preds$sim_idx, actual$sim_idx), N_P3)
cat(sprintf("Part 3: %d sims x %d windows x %d runs each\n",
            N_P3, length(TEST_WINDOWS), NSIMS))

# ── Pre-compute actual-param SIR quantiles for these sims (full 365 days) ────
cat("  Pre-computing actual-param baselines...\n")
act_baseline <- lapply(p3_ids, function(sid) {
  act <- actual |> filter(sim_idx == sid)
  q   <- run_sir_multi(n = act$n[1], prevalence = act$prevalence[1],
                       ptran = act$ptran[1], crate = act$crate[1],
                       recov = act$recov[1], ndays = 365, seed = sid)
  list(sid = sid, q = q,
       true_ptran = act$ptran[1], true_crate = act$crate[1],
       true_R0    = act$R0[1],    recov = act$recov[1],
       n = act$n[1], prevalence = act$prevalence[1])
})
names(act_baseline) <- as.character(p3_ids)

# ── Main loop: for each window, calibrate and recover ─────────────────────────
window_rows   <- list()   # parameter errors + curve MAE
example_curves <- list()  # full curves for plotting (subset of windows)

EXAMPLE_SID <- p3_ids[1]  # use first sim for example curve plots

for (win_tag in names(TEST_WINDOWS)) {
  w      <- TEST_WINDOWS[[win_tag]]
  t_start <- w["start"] + 1L   # R is 1-indexed
  t_end   <- t_start + w["len"] - 1L
  regime  <- sub("_.*", "", win_tag)
  win_len <- w["len"]

  cat(sprintf("  Window: %-14s  [day %3d - %3d]\n",
              win_tag, t_start, t_end))

  for (bl in act_baseline) {
    sid <- bl$sid
    mat_row <- which(actual$sim_idx == sid)
    obs_full <- as.numeric(inc_raw[mat_row, ])          # full 365-day observed
    obs_win  <- obs_full[t_start:min(t_end, 365)]       # sliced window

    # calibrate from this window slice
    pred <- tryCatch(
      calibrate_sir(daily_cases     = obs_win,
                    population_size = bl$n,
                    recovery_rate   = bl$recov),
      error = function(e) NULL)
    if (is.null(pred)) next

    # recover full 365-day curve with predicted params
    pred_q <- run_sir_multi(n = bl$n, prevalence = bl$prevalence,
                            ptran = pred["ptran"], crate = pred["crate"],
                            recov = bl$recov, ndays = 365, seed = sid)

    # parameter errors
    window_rows[[length(window_rows) + 1]] <- data.frame(
      window  = win_tag,
      regime  = regime,
      win_len = win_len,
      sim_idx = sid,
      err_ptran = abs(pred["ptran"] - bl$true_ptran),
      err_crate = abs(pred["crate"] - bl$true_crate),
      err_R0    = abs(pred["R0"]    - bl$true_R0),
      pred_ptran = pred["ptran"],
      pred_crate = pred["crate"],
      pred_R0    = pred["R0"],
      true_ptran = bl$true_ptran,
      true_crate = bl$true_crate,
      true_R0    = bl$true_R0,
      # curve MAE (predicted median vs actual median, full 365 days)
      curve_mae  = mean(abs(pred_q$med - bl$q$med)),
      # peak comparison
      pred_peak_day  = which.max(pred_q$med),
      pred_peak_size = max(pred_q$med),
      true_peak_day  = which.max(bl$q$med),
      true_peak_size = max(bl$q$med),
      stringsAsFactors = FALSE
    )

    # store example curves for one sim
    if (sid == EXAMPLE_SID &&
        win_tag %in% c("early_015d", "early_030d", "early_060d",
                       "mid_090d",  "mid_180d",
                       "late_365d")) {
      example_curves[[win_tag]] <- data.frame(
        day        = 1:365,
        win_tag    = win_tag,
        regime     = regime,
        win_len    = win_len,
        obs        = obs_full,
        act_med    = bl$q$med,
        act_lower  = bl$q$lower,
        act_upper  = bl$q$upper,
        pred_med   = pred_q$med,
        pred_lower = pred_q$lower,
        pred_upper = pred_q$upper,
        # shade the window used
        in_window  = (1:365) >= t_start & (1:365) <= min(t_end, 365)
      )
    }
  }
}

win_df  <- bind_rows(window_rows)
ex_df   <- bind_rows(example_curves)

# ── Summarise across sims ─────────────────────────────────────────────────────
win_summary <- win_df |>
  group_by(window, regime, win_len) |>
  summarise(across(c(err_ptran, err_crate, err_R0, curve_mae,
                     pred_peak_day, pred_peak_size,
                     true_peak_day, true_peak_size),
                   mean, .names = "mean_{.col}"),
            .groups = "drop") |>
  mutate(peak_day_err  = abs(mean_pred_peak_day  - mean_true_peak_day),
         peak_size_err = abs(mean_pred_peak_size - mean_true_peak_size),
         regime = factor(regime, levels = c("early", "mid", "late")))

# ── Plot A: Parameter error vs window length ──────────────────────────────────
param_long <- win_summary |>
  select(win_len, regime, mean_err_ptran, mean_err_crate, mean_err_R0) |>
  pivot_longer(cols = starts_with("mean_err"),
               names_to = "param", values_to = "mae") |>
  mutate(param = recode(param,
                        mean_err_ptran = "ptran",
                        mean_err_crate = "crate",
                        mean_err_R0    = "R0"))

pA <- ggplot(param_long,
             aes(x = win_len, y = mae, color = regime, group = regime)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  facet_wrap(~ param, scales = "free_y", ncol = 3) +
  scale_x_log10(breaks = LENGTHS) +
  scale_color_manual(values = c(early = "#E65100", mid = "#1565C0",
                                late  = "#2E7D32")) +
  labs(title    = "Part 3A: Parameter Error vs Window Length",
       subtitle = sprintf("Mean MAE across %d val sims", N_P3),
       x = "Window length (days, log scale)", y = "Mean MAE", color = "Regime") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

ggsave(path.expand(file.path(DATA_DIR, "part3A_param_error_vs_window.png")),
       pA, width = 13, height = 5, dpi = 150)
print(pA)

# ── Plot B: Recovered-curve MAE vs window length ──────────────────────────────
pB <- ggplot(win_summary,
             aes(x = win_len, y = mean_curve_mae,
                 color = regime, group = regime)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  scale_x_log10(breaks = LENGTHS) +
  scale_color_manual(values = c(early = "#E65100", mid = "#1565C0",
                                late  = "#2E7D32")) +
  labs(title    = "Part 3B: Recovered-Curve MAE vs Window Length",
       subtitle = "MAE between predicted-param and actual-param SIR medians (full 365 days)",
       x = "Window length (days, log scale)", y = "Mean curve MAE", color = "Regime") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

ggsave(path.expand(file.path(DATA_DIR, "part3B_curve_mae_vs_window.png")),
       pB, width = 9, height = 5, dpi = 150)
print(pB)

# ── Plot C: Example recovered curves across window lengths for one sim ─────────
if (nrow(ex_df) > 0) {
  ex_df <- ex_df |>
    mutate(label = sprintf("%s (%d days)", win_tag, win_len),
           label = factor(label, levels = unique(label[order(win_len)])))

  # Pre-compute shading rectangles: one row per window
  rect_df <- ex_df |>
    group_by(win_tag, label) |>
    summarise(xmin = min(day[in_window]),
              xmax = max(day[in_window]),
              .groups = "drop")

  pC <- ggplot(ex_df, aes(x = day)) +
    # shade the observation window
    geom_rect(data = rect_df,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = "gold", alpha = 0.15, inherit.aes = FALSE) +
    # actual param CI + median
    geom_ribbon(aes(ymin = act_lower, ymax = act_upper),
                fill = "#1976D2", alpha = 0.15) +
    geom_line(aes(y = act_med, color = "Actual params"),
              linewidth = 0.7) +
    # predicted param CI + median
    geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper),
                fill = "#D32F2F", alpha = 0.15) +
    geom_line(aes(y = pred_med, color = "Predicted params"),
              linewidth = 0.7, linetype = "dashed") +
    # observed ABM
    geom_line(aes(y = obs, color = "Observed (ABM)"),
              linewidth = 0.5, linetype = "dotdash") +
    facet_wrap(~ label, ncol = 3, scales = "free_y") +
    scale_color_manual(values = c("Observed (ABM)"  = "black",
                                  "Actual params"   = "#1976D2",
                                  "Predicted params" = "#D32F2F")) +
    labs(title    = sprintf("Part 3C: Recovered Curves for sim %d across Window Lengths", EXAMPLE_SID),
         subtitle = "Gold shading = observation window used for calibration",
         x = "Day", y = "Daily Incidence", color = NULL) +
    theme_bw(base_size = 9) +
    theme(legend.position  = "bottom",
          strip.text       = element_text(size = 8),
          panel.grid.minor = element_blank(),
          plot.title       = element_text(face = "bold"))

  ggsave(path.expand(file.path(DATA_DIR, "part3C_example_windows.png")),
         pC, width = 14, height = 10, dpi = 150)
  print(pC)
}

# ── Plot D: Can we predict the peak? ─────────────────────────────────────────
pD1 <- ggplot(win_summary,
              aes(x = win_len, y = peak_day_err,
                  color = regime, group = regime)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.5) +
  scale_x_log10(breaks = LENGTHS) +
  scale_color_manual(values = c(early = "#E65100", mid = "#1565C0",
                                late  = "#2E7D32")) +
  labs(title = "Part 3D: Peak Day Error vs Window Length",
       x = "Window length (days, log scale)",
       y = "Mean |predicted peak day - actual peak day|",
       color = "Regime") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

pD2 <- ggplot(win_summary,
              aes(x = win_len, y = peak_size_err,
                  color = regime, group = regime)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.5) +
  scale_x_log10(breaks = LENGTHS) +
  scale_color_manual(values = c(early = "#E65100", mid = "#1565C0",
                                late  = "#2E7D32")) +
  labs(title = "Part 3D: Peak Size Error vs Window Length",
       x = "Window length (days, log scale)",
       y = "Mean |predicted peak size - actual peak size|",
       color = "Regime") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

library(patchwork)
pD <- pD1 / pD2
ggsave(path.expand(file.path(DATA_DIR, "part3D_peak_prediction.png")),
       pD, width = 10, height = 9, dpi = 150)
print(pD)

# ── Summary table ─────────────────────────────────────────────────────────────
cat("\n-- Part 3 summary table --\n")
win_summary |>
  select(window, regime, win_len,
         mean_err_R0, mean_curve_mae, peak_day_err, peak_size_err) |>
  arrange(regime, win_len) |>
  mutate(across(where(is.numeric), round, 3)) |>
  print(n = 40)

cat("\nSaved:\n")
cat("  part3A_param_error_vs_window.png\n")
cat("  part3B_curve_mae_vs_window.png\n")
cat("  part3C_example_windows.png\n")
cat("  part3D_peak_prediction.png\n")
