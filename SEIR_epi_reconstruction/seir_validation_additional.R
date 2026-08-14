# =============================================================================
#  SEIR Calibration Validation — hypertuned BiLSTM (Optuna)
#
#  Part 1 — Test set: run ModelSEIRCONN with actual vs BiLSTM-predicted
#            parameters, compare incidence curves + 95% CI, per-day MAE
#
#  Part 2 — Random parameters pipeline:
#            sample random params -> simulate incidence -> calibrate_seir()
#            -> simulate with predicted params -> plot comparison with CI
#
#  Part 3 — Window analysis: inference across all 18 windows
#            (6 lengths x 3 regimes: early / mid / late)
#
#  NOTE: Day 1 of simulation is excluded from all plots (initialization effect).
#
#  SEIR known inputs : n, recov, incub, prevalence
#  SEIR predicted    : beta, R0
#  ModelSEIRCONN     : contact_rate=CONTACT_RATE, transmission_rate=beta/CONTACT_RATE
# =============================================================================

library(epiworldR)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# -- Source calibrate_seir if not already loaded ------------------------------
if (!exists("init_bilstm_model")) {
  source(path.expand("~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction/calibrate_seir.R"))
}

# -- Paths --------------------------------------------------------------------
MODEL_DIR <- "~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction/model"
DATA_DIR  <- "~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction"

NDAYS    <- 365
N_SAMPLE <- 20
WINDOW       <- "late_365d"
CONTACT_RATE <- 1
NSIMS    <- 500
NTHREADS <- 4

# -- Load SEIR BiLSTM model ---------------------------------------------------
init_bilstm_model(MODEL_DIR)

# =============================================================================
# Helper: run ModelSEIRCONN NSIMS times, return median + 95% CI per day
# =============================================================================
run_seir_multi <- function(n, prevalence, beta, incub, recov,
                           ndays = NDAYS, nsims = NSIMS, nthreads = NTHREADS,
                           seed = 1L) {
  n    <- max(as.integer(round(n)), 10L)
  prev <- max(min(prevalence, 1.0), 1.0 / n)

  model <- ModelSEIRCONN(
    name              = "sim",
    n                 = n,
    prevalence        = prev,
    contact_rate      = CONTACT_RATE,
    transmission_rate = max(beta / CONTACT_RATE, 1e-6),
    incubation_days   = max(incub, 1.0),
    recovery_rate     = max(recov, 1e-6)
  )
  verbose_off(model)

  saver <- make_saver("transition")
  run_multiple(model, ndays = ndays + 1L, nsims = nsims,
               saver = saver, nthreads = nthreads)

  res <- run_multiple_get_results(model, nthreads = nthreads,
                                  freader = data.table::fread)

  res$transition |>
    filter(from == "Exposed", to == "Infected", date > 0) |>
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

# =============================================================================
# PART 1 -- Test set: actual vs predicted incidence
# =============================================================================

cat("\n-- Part 1: loading test data --\n")
actual    <- read.csv(path.expand(file.path(DATA_DIR, "test_actual_parameters.csv")))
preds_all <- read.csv(path.expand(file.path(DATA_DIR, "test_bilstm_predictions_tuned.csv")))

np      <- reticulate::import("numpy")
inc_raw <- np$load(path.expand(file.path(DATA_DIR, "test_incidence_raw.npy")))

preds <- preds_all |> filter(window == WINDOW)
cat(sprintf("  Val sims: %d  |  Window: %s  |  NSIMS per run: %d\n",
            nrow(actual), WINDOW, NSIMS))

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
  incub      <- act$incub[1]
  prevalence <- act$prevalence[1]

  mat_row <- which(actual$sim_idx == sid)
  obs_inc <- as.numeric(inc_raw[mat_row, ])

  act_q <- run_seir_multi(
    n = n, prevalence = prevalence,
    beta = act$beta[1], incub = incub, recov = recov, seed = sid)

  pred_q <- run_seir_multi(
    n = n, prevalence = prevalence,
    beta = prd$beta_pred[1], incub = incub, recov = recov, seed = sid)

  curves_list[[k]] <- data.frame(
    day        = 1:NDAYS,
    sim_idx    = sid,
    panel      = sprintf(
      "sim %d\nbeta: %.3f->%.3f | R0: %.2f->%.2f | incub=%.1fd",
      sid,
      act$beta[1], prd$beta_pred[1],
      act$R0[1],   prd$R0_pred[1],
      incub
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

# -- Plot 1: incidence curves (day > 1) ---------------------------------------
p1 <- ggplot(curves_df |> filter(day > 1), aes(x = day)) +
  geom_ribbon(aes(ymin = act_lower, ymax = act_upper),
              fill = "#1976D2", alpha = 0.20) +
  geom_line(aes(y = act_med, color = "SEIR - Actual Params"),
            linewidth = 1.1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper),
              fill = "#D32F2F", alpha = 0.20) +
  geom_line(aes(y = pred_med, color = "SEIR - BiLSTM Predicted"),
            linewidth = 1.1) +
  geom_line(aes(y = obs, color = "Observed (ABM)"),
            linewidth = 0.5, alpha = 0.7) +
  facet_wrap(~ panel, scales = "free_y", ncol = 4) +
  scale_color_manual(
    values = c("Observed (ABM)"          = "black",
               "SEIR - Actual Params"    = "#1976D2",
               "SEIR - BiLSTM Predicted" = "#D32F2F"),
    guide  = guide_legend(override.aes = list(linewidth = c(0.5, 1.1, 1.1),
                                              alpha     = c(0.7, 1.0, 1.0)))) +
  labs(title    = paste0("SEIR Incidence - Actual vs BiLSTM Predicted [tuned] (", WINDOW, ")"),
       subtitle = paste0("Shaded bands = 95% CI across ", NSIMS,
                         " runs | Days 2-365 shown"),
       x = "Day", y = "Daily Incidence (new exposures)", color = NULL) +
  theme_bw(base_size = 9) +
  theme(legend.position  = "bottom",
        legend.key.width = unit(1.8, "cm"),
        strip.text       = element_text(size = 6.5),
        panel.grid.minor = element_blank(),
        plot.title       = element_text(face = "bold"))

ggsave(path.expand(file.path(DATA_DIR, "part1_seir_incidence_curves.png")),
       p1, width = 18, height = 14, dpi = 150)
print(p1)

# -- Plot 2: per-day MAE (day > 1) --------------------------------------------
mae_avg <- mae_day_df |>
  filter(day > 1) |>
  group_by(day) |>
  summarise(mean_mae = mean(mae), sd_mae = sd(mae), .groups = "drop")

p2 <- ggplot(mae_avg, aes(x = day)) +
  geom_ribbon(aes(ymin = pmax(mean_mae - sd_mae, 0),
                  ymax = mean_mae + sd_mae),
              fill = "#D32F2F", alpha = 0.15) +
  geom_line(aes(y = mean_mae), color = "#D32F2F", linewidth = 0.9) +
  labs(title    = "Per-day MAE: Actual-Param vs Predicted-Param SEIR (medians)",
       subtitle = sprintf("Mean +/- SD across %d sampled simulations | Days 2-365 shown",
                          N_SAMPLE),
       x = "Day", y = "MAE (incidence)") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

ggsave(path.expand(file.path(DATA_DIR, "part1_seir_per_day_mae.png")),
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
# PART 2 -- Random parameters pipeline
# =============================================================================

cat("\n-- Part 2: random parameters pipeline --\n")

N_RANDOM  <- 15
NDAYS_RAN <- 365
set.seed(123)

random_params <- data.frame(
  n          = sample(5000:10000, N_RANDOM, replace = TRUE),
  prevalence = runif(N_RANDOM, min = 0.03, max = 0.35),
  beta       = runif(N_RANDOM, min = 0.05, max = 0.50),
  recov      = runif(N_RANDOM, min = 0.07, max = 0.25),
  incub      = runif(N_RANDOM, min = 2.0,  max = 10.0)
)
random_params$R0 <- random_params$beta / random_params$recov

cat("Random parameters sampled:\n")
print(round(random_params, 4))

random_curves <- vector("list", N_RANDOM)

for (k in seq_len(N_RANDOM)) {
  rp <- random_params[k, ]
  cat(sprintf("\n  [%d/%d]  n=%d  beta=%.3f  recov=%.3f  incub=%.1f  R0=%.2f\n",
              k, N_RANDOM, rp$n, rp$beta, rp$recov, rp$incub, rp$R0))

  true_q <- run_seir_multi(
    n = rp$n, prevalence = rp$prevalence,
    beta = rp$beta, incub = rp$incub, recov = rp$recov,
    ndays = NDAYS_RAN, seed = k)

  predicted <- calibrate_seir(
    daily_cases     = true_q$med,
    population_size = rp$n,
    recovery_rate   = rp$recov,
    incubation_days = rp$incub
  )
  cat(sprintf("    True      : beta=%.4f  R0=%.4f\n", rp$beta, rp$R0))
  cat(sprintf("    Predicted : beta=%.4f  R0=%.4f\n",
              predicted["beta"], predicted["R0"]))

  pred_q <- run_seir_multi(
    n = rp$n, prevalence = rp$prevalence,
    beta  = predicted["beta"],
    incub = rp$incub,
    recov = rp$recov,
    ndays = NDAYS_RAN, seed = k)

  random_curves[[k]] <- data.frame(
    day        = 1:NDAYS_RAN,
    panel      = sprintf(
      "sim %d | R0: %.2f->%.2f\nbeta: %.3f->%.3f | incub=%.1fd",
      k, rp$R0, predicted["R0"],
      rp$beta, predicted["beta"], rp$incub
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

# -- Summary CSV --------------------------------------------------------------
random_summary <- random_df |>
  filter(day > 1) |>
  group_by(panel) |>
  summarise(
    mae_medians    = mean(abs(true_med - pred_med)),
    peak_true      = max(true_med),
    peak_predicted = max(pred_med),
    peak_day_true  = which.max(true_med),
    peak_day_pred  = which.max(pred_med),
    .groups = "drop"
  ) |>
  mutate(peak_day_err = abs(peak_day_true - peak_day_pred))

write.csv(random_summary,
          path.expand(file.path(DATA_DIR, "part2_random_pipeline_summary.csv")),
          row.names = FALSE)

# -- Plot 3: random pipeline (day > 1) ----------------------------------------
p3 <- ggplot(random_df |> filter(day > 1), aes(x = day)) +
  geom_ribbon(aes(ymin = true_lower, ymax = true_upper),
              fill = "#1976D2", alpha = 0.20) +
  geom_line(aes(y = true_med, color = "SEIR - True (Random) Params"),
            linewidth = 1.1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper),
              fill = "#D32F2F", alpha = 0.20) +
  geom_line(aes(y = pred_med, color = "SEIR - BiLSTM Predicted"),
            linewidth = 1.1) +
  facet_wrap(~ panel, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c(
    "SEIR - True (Random) Params" = "#1976D2",
    "SEIR - BiLSTM Predicted"     = "#D32F2F")) +
  labs(title    = "Part 2 - Random Parameters: True vs BiLSTM-Predicted SEIR Curves [tuned]",
       subtitle = paste0("Shaded bands = 95% CI across ", NSIMS,
                         " runs | Days 2-365 shown"),
       x = "Day", y = "Daily Incidence (new exposures)", color = NULL) +
  theme_bw(base_size = 9) +
  theme(legend.position  = "bottom",
        legend.key.width = unit(1.8, "cm"),
        strip.text       = element_text(size = 7),
        panel.grid.minor = element_blank(),
        plot.title       = element_text(face = "bold"))

ggsave(path.expand(file.path(DATA_DIR, "part2_seir_random_pipeline.png")),
       p3, width = 15, height = 16, dpi = 150)
print(p3)

cat("\n-- Part 2 summary --\n")
print(random_summary)

# =============================================================================
# PART 3 -- Window analysis: inference across all 18 windows
# =============================================================================

cat("\n-- Part 3: window analysis --\n")

LENGTHS <- c(15, 30, 60, 90, 180, 365)
make_windows <- function(t_max = 365, lengths = LENGTHS) {
  out <- list()
  for (L in lengths) {
    if (L > t_max) next
    out[[sprintf("early_%03dd", L)]] <- c(start = 0L, len = L)
    mid_s <- max(0L, as.integer((t_max - L) / 2))
    out[[sprintf("mid_%03dd",   L)]] <- c(start = mid_s, len = L)
    late_s <- max(0L, t_max - L)
    out[[sprintf("late_%03dd",  L)]] <- c(start = late_s, len = L)
  }
  out
}
TEST_WINDOWS <- make_windows()

N_P3 <- 10
set.seed(7)
p3_ids <- sample(intersect(preds$sim_idx, actual$sim_idx), N_P3)
cat(sprintf("Part 3: %d sims x %d windows x %d runs each\n",
            N_P3, length(TEST_WINDOWS), NSIMS))

cat("  Pre-computing actual-param baselines...\n")
act_baseline <- lapply(p3_ids, function(sid) {
  act <- actual |> filter(sim_idx == sid)
  q   <- run_seir_multi(n = act$n[1], prevalence = act$prevalence[1],
                        beta  = act$beta[1], incub = act$incub[1],
                        recov = act$recov[1], ndays = 365, seed = sid)
  list(sid        = sid,       q          = q,
       true_beta  = act$beta[1], true_R0  = act$R0[1],
       recov      = act$recov[1], incub   = act$incub[1],
       n          = act$n[1],   prevalence = act$prevalence[1])
})
names(act_baseline) <- as.character(p3_ids)

window_rows    <- list()
example_curves <- list()
EXAMPLE_SID    <- p3_ids[1]

for (win_tag in names(TEST_WINDOWS)) {
  w       <- TEST_WINDOWS[[win_tag]]
  t_start <- w["start"] + 1L
  t_end   <- t_start + w["len"] - 1L
  regime  <- sub("_.*", "", win_tag)
  win_len <- w["len"]

  cat(sprintf("  Window: %-14s  [day %3d - %3d]\n", win_tag, t_start, t_end))

  for (bl in act_baseline) {
    sid      <- bl$sid
    mat_row  <- which(actual$sim_idx == sid)
    obs_full <- as.numeric(inc_raw[mat_row, ])
    obs_win  <- obs_full[t_start:min(t_end, 365)]

    pred <- tryCatch(
      calibrate_seir(daily_cases     = obs_win,
                     population_size = bl$n,
                     recovery_rate   = bl$recov,
                     incubation_days = bl$incub),
      error = function(e) NULL)
    if (is.null(pred)) next

    pred_q <- run_seir_multi(n = bl$n, prevalence = bl$prevalence,
                             beta  = pred["beta"],
                             incub = bl$incub,
                             recov = bl$recov,
                             ndays = 365, seed = sid)

    window_rows[[length(window_rows) + 1]] <- data.frame(
      window   = win_tag, regime  = regime, win_len = win_len, sim_idx = sid,
      err_beta = abs(pred["beta"] - bl$true_beta),
      err_R0   = abs(pred["R0"]   - bl$true_R0),
      pred_beta = pred["beta"],   pred_R0 = pred["R0"],
      true_beta = bl$true_beta,   true_R0 = bl$true_R0,
      curve_mae      = mean(abs(pred_q$med - bl$q$med)),
      pred_peak_day  = which.max(pred_q$med),
      pred_peak_size = max(pred_q$med),
      true_peak_day  = which.max(bl$q$med),
      true_peak_size = max(bl$q$med),
      stringsAsFactors = FALSE
    )

    if (sid == EXAMPLE_SID &&
        win_tag %in% c("early_015d", "early_030d", "early_060d",
                       "mid_090d",   "mid_180d",   "late_365d")) {
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
        in_window  = (1:365) >= t_start & (1:365) <= min(t_end, 365)
      )
    }
  }
}

win_df <- bind_rows(window_rows)
ex_df  <- bind_rows(example_curves)

win_summary <- win_df |>
  group_by(window, regime, win_len) |>
  summarise(across(c(err_beta, err_R0, curve_mae,
                     pred_peak_day, pred_peak_size,
                     true_peak_day, true_peak_size),
                   mean, .names = "mean_{.col}"),
            .groups = "drop") |>
  mutate(peak_day_err  = abs(mean_pred_peak_day  - mean_true_peak_day),
         peak_size_err = abs(mean_pred_peak_size - mean_true_peak_size),
         regime = factor(regime, levels = c("early", "mid", "late")))

# -- Plot A: Parameter error vs window length ---------------------------------
param_long <- win_summary |>
  select(win_len, regime, mean_err_beta, mean_err_R0) |>
  pivot_longer(cols = starts_with("mean_err"),
               names_to = "param", values_to = "mae") |>
  mutate(param = recode(param, mean_err_beta = "beta", mean_err_R0 = "R0"))

pA <- ggplot(param_long,
             aes(x = win_len, y = mae, color = regime, group = regime)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.5) +
  facet_wrap(~ param, scales = "free_y", ncol = 2) +
  scale_x_log10(breaks = LENGTHS) +
  scale_color_manual(values = c(early = "#E65100", mid = "#1565C0",
                                late  = "#2E7D32")) +
  labs(title    = "Part 3A: Parameter Error vs Window Length (SEIR, tuned)",
       subtitle = sprintf("Mean MAE across %d val sims", N_P3),
       x = "Window length (days, log scale)", y = "Mean MAE", color = "Regime") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

ggsave(path.expand(file.path(DATA_DIR, "part3A_seir_param_error_vs_window.png")),
       pA, width = 10, height = 5, dpi = 150)
print(pA)

# -- Plot B: Recovered-curve MAE vs window length ----------------------------
pB <- ggplot(win_summary,
             aes(x = win_len, y = mean_curve_mae, color = regime, group = regime)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.5) +
  scale_x_log10(breaks = LENGTHS) +
  scale_color_manual(values = c(early = "#E65100", mid = "#1565C0",
                                late  = "#2E7D32")) +
  labs(title    = "Part 3B: Recovered-Curve MAE vs Window Length (SEIR, tuned)",
       subtitle = "MAE between predicted-param and actual-param SEIR medians (full 365 days)",
       x = "Window length (days, log scale)", y = "Mean curve MAE", color = "Regime") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

ggsave(path.expand(file.path(DATA_DIR, "part3B_seir_curve_mae_vs_window.png")),
       pB, width = 9, height = 5, dpi = 150)
print(pB)

# -- Plot C: Example recovered curves (day > 1) -------------------------------
if (nrow(ex_df) > 0) {
  ex_df <- ex_df |>
    filter(day > 1) |>
    mutate(label = sprintf("%s (%d days)", win_tag, win_len),
           label = factor(label, levels = unique(label[order(win_len)])))

  rect_df <- ex_df |>
    group_by(win_tag, label) |>
    summarise(xmin = min(day[in_window]),
              xmax = max(day[in_window]),
              .groups = "drop")

  pC <- ggplot(ex_df, aes(x = day)) +
    geom_rect(data = rect_df,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = "gold", alpha = 0.18, inherit.aes = FALSE) +
    geom_ribbon(aes(ymin = act_lower,  ymax = act_upper),
                fill = "#1976D2", alpha = 0.20) +
    geom_line(aes(y = act_med,  color = "Actual params"),  linewidth = 1.0) +
    geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper),
                fill = "#D32F2F", alpha = 0.20) +
    geom_line(aes(y = pred_med, color = "Predicted params"), linewidth = 1.0) +
    geom_line(aes(y = obs, color = "Observed (ABM)"),
              linewidth = 0.5, alpha = 0.7) +
    facet_wrap(~ label, ncol = 3, scales = "free_y") +
    scale_color_manual(
      values = c("Observed (ABM)"   = "black",
                 "Actual params"    = "#1976D2",
                 "Predicted params" = "#D32F2F"),
      guide  = guide_legend(override.aes = list(linewidth = c(0.5, 1.0, 1.0),
                                                alpha     = c(0.7, 1.0, 1.0)))) +
    labs(title    = sprintf("Part 3C: Recovered SEIR Curves for sim %d across Window Lengths [tuned]",
                            EXAMPLE_SID),
         subtitle = "Gold shading = observation window | Shaded bands = 95% CI | Days 2-365 shown",
         x = "Day", y = "Daily Incidence (new exposures)", color = NULL) +
    theme_bw(base_size = 9) +
    theme(legend.position  = "bottom",
          strip.text       = element_text(size = 8),
          panel.grid.minor = element_blank(),
          plot.title       = element_text(face = "bold"))

  ggsave(path.expand(file.path(DATA_DIR, "part3C_seir_example_windows.png")),
         pC, width = 14, height = 10, dpi = 150)
  print(pC)
}

# -- Plot D: Peak prediction --------------------------------------------------
pD1 <- ggplot(win_summary,
              aes(x = win_len, y = peak_day_err, color = regime, group = regime)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.5) +
  scale_x_log10(breaks = LENGTHS) +
  scale_color_manual(values = c(early = "#E65100", mid = "#1565C0",
                                late  = "#2E7D32")) +
  labs(title = "Part 3D: Peak Day Error vs Window Length (SEIR, tuned)",
       x = "Window length (days, log scale)",
       y = "Mean |predicted peak day - actual peak day|",
       color = "Regime") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

pD2 <- ggplot(win_summary,
              aes(x = win_len, y = peak_size_err, color = regime, group = regime)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.5) +
  scale_x_log10(breaks = LENGTHS) +
  scale_color_manual(values = c(early = "#E65100", mid = "#1565C0",
                                late  = "#2E7D32")) +
  labs(title = "Part 3D: Peak Size Error vs Window Length (SEIR, tuned)",
       x = "Window length (days, log scale)",
       y = "Mean |predicted peak size - actual peak size|",
       color = "Regime") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

pD <- pD1 / pD2
ggsave(path.expand(file.path(DATA_DIR, "part3D_seir_peak_prediction.png")),
       pD, width = 10, height = 9, dpi = 150)
print(pD)

cat("\n-- Part 3 summary table --\n")
win_summary |>
  select(window, regime, win_len,
         mean_err_beta, mean_err_R0, mean_curve_mae,
         peak_day_err, peak_size_err) |>
  arrange(regime, win_len) |>
  mutate(across(where(is.numeric), round, 3)) |>
  print(n = 40)

cat("\nSaved:\n")
cat("  part1_seir_incidence_curves.png\n")
cat("  part1_seir_per_day_mae.png\n")
cat("  part2_seir_random_pipeline.png\n")
cat("  part2_random_pipeline_summary.csv\n")
cat("  part3A_seir_param_error_vs_window.png\n")
cat("  part3B_seir_curve_mae_vs_window.png\n")
cat("  part3C_seir_example_windows.png\n")
cat("  part3D_seir_peak_prediction.png\n")
