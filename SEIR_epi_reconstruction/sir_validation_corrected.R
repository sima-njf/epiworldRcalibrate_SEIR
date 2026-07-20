# =============================================================================
#  SIR Calibration Validation — with run_multiple (uncertainty bands)
#
#  FIX: run_sir_multi now uses make_saver("transition") + filters
#       from=="Susceptible", to=="Infected", date>0
#       so the single-run blue line and the ribbon both use the same
#       metric and both drop the day-0 initialization spike.
# =============================================================================

library(epiworldR)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(data.table)
library(reticulate)

# -- Paths --------------------------------------------------------------------
MODEL_DIR <- "~/sima/epiworldRcalibrate_SEIR/SIR_epi_reconstruction/model"
DATA_DIR  <- "~/sima/epiworldRcalibrate_SEIR/SIR_epi_reconstruction"

NDAYS    <- 365
N_SAMPLE <- 20
WINDOW   <- "late_365d"
NSIMS    <- 1000
NTHREADS <- 10

# -- Load SIR BiLSTM model ----------------------------------------------------
source("~/sima/epiworldRcalibrate_SEIR/SIR_epi_reconstruction/calibrate_sir.R")
init_bilstm_model(MODEL_DIR)

# =============================================================================
# Helper: run ModelSIRCONN NSIMS times, return median + 95% CI per day
#
# FIX applied here:
#   1. run ndays+1 so we can drop the day-0 initialization spike
#   2. filter(from=="Susceptible", to=="Infected", date>0)
#      — same metric as what ModelSIRCONN reports for incidence
#   3. right_join to dates 1:ndays (not 0:(ndays-1))
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
  run_multiple(model, ndays = ndays + 1L, nsims = nsims,
               saver = saver, nthreads = nthreads)

  res <- run_multiple_get_results(model, nthreads = nthreads,
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

# =============================================================================
# PART 1 — Test set: actual vs predicted incidence
# =============================================================================

cat("\n-- Part 1: loading test data --\n")
actual    <- read.csv(path.expand(file.path(DATA_DIR, "test_actual_parameters.csv")))
preds_all <- read.csv(path.expand(file.path(DATA_DIR, "test_bilstm_predictions.csv")))

np      <- reticulate::import("numpy")
inc_raw <- as.matrix(np$load(path.expand(
  file.path(DATA_DIR, "test_incidence_raw.npy"))))

preds <- preds_all |> filter(window == WINDOW)
cat(sprintf("  Val sims: %d  |  Window: %s  |  NSIMS: %d\n",
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

# -- Plot 1: incidence curves with CI ribbons ---------------------------------
p1 <- ggplot(curves_df, aes(x = day)) +
  geom_ribbon(aes(ymin = act_lower,  ymax = act_upper),
              fill = "#1976D2", alpha = 0.20) +
  geom_line(aes(y = act_med,  color = "SIR - Actual Params"),
            linewidth = 1.1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper),
              fill = "#D32F2F", alpha = 0.20) +
  geom_line(aes(y = pred_med, color = "SIR - BiLSTM Predicted"),
            linewidth = 1.1) +
  geom_line(aes(y = obs,      color = "Observed (ABM)"),
            linewidth = 0.5, alpha = 0.7) +
  facet_wrap(~ panel, scales = "free_y", ncol = 4) +
  scale_color_manual(
    values = c("Observed (ABM)"         = "black",
               "SIR - Actual Params"    = "#1976D2",
               "SIR - BiLSTM Predicted" = "#D32F2F"),
    guide  = guide_legend(override.aes = list(
      linewidth = c(0.5, 1.1, 1.1), alpha = c(0.7, 1.0, 1.0)))) +
  labs(title    = paste0("SIR Incidence - Actual vs BiLSTM Predicted (", WINDOW, ")"),
       subtitle = paste0("Shaded bands = 95% CI across ", NSIMS, " stochastic runs"),
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

# -- Plot 2: per-day MAE (on medians) ----------------------------------------
mae_avg <- mae_day_df |>
  group_by(day) |>
  summarise(mean_mae = mean(mae), sd_mae = sd(mae), .groups = "drop")

p2 <- ggplot(mae_avg, aes(x = day)) +
  geom_ribbon(aes(ymin = pmax(mean_mae - sd_mae, 0),
                  ymax = mean_mae + sd_mae),
              fill = "#D32F2F", alpha = 0.15) +
  geom_line(aes(y = mean_mae), color = "#D32F2F", linewidth = 0.9) +
  labs(title    = "Per-day MAE: Actual-Param vs Predicted-Param SIR (medians)",
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

cat("Random parameters:\n"); print(round(random_params, 4))

random_curves <- vector("list", N_RANDOM)

for (k in seq_len(N_RANDOM)) {
  rp <- random_params[k, ]
  cat(sprintf("\n  [%d/%d] n=%d ptran=%.3f crate=%.2f recov=%.3f R0=%.2f\n",
              k, N_RANDOM, rp$n, rp$ptran, rp$crate, rp$recov, rp$R0))

  true_q <- run_sir_multi(
    n = rp$n, prevalence = rp$prevalence,
    ptran = rp$ptran, crate = rp$crate, recov = rp$recov,
    ndays = NDAYS_RAN, seed = k)

  predicted <- calibrate_sir(
    daily_cases     = true_q$med,
    population_size = rp$n,
    recovery_rate   = rp$recov)

  cat(sprintf("    True      : ptran=%.4f  crate=%.4f  R0=%.4f\n",
              rp$ptran, rp$crate, rp$R0))
  cat(sprintf("    Predicted : ptran=%.4f  crate=%.4f  R0=%.4f\n",
              predicted["ptran"], predicted["crate"], predicted["R0"]))

  pred_q <- run_sir_multi(
    n = rp$n, prevalence = rp$prevalence,
    ptran = predicted["ptran"], crate = predicted["crate"], recov = rp$recov,
    ndays = NDAYS_RAN, seed = k)

  random_curves[[k]] <- data.frame(
    day        = 1:NDAYS_RAN,
    panel      = sprintf(
      "sim %d | R0: %.2f->%.2f\nptran: %.3f->%.3f | crate: %.2f->%.2f",
      k, rp$R0, predicted["R0"],
      rp$ptran, predicted["ptran"],
      rp$crate, predicted["crate"]),
    true_lower = true_q$lower,  true_med  = true_q$med,  true_upper = true_q$upper,
    pred_lower = pred_q$lower,  pred_med  = pred_q$med,  pred_upper = pred_q$upper
  )
}

random_df <- bind_rows(random_curves)

p3 <- ggplot(random_df, aes(x = day)) +
  geom_ribbon(aes(ymin = true_lower, ymax = true_upper),
              fill = "#1976D2", alpha = 0.20) +
  geom_line(aes(y = true_med, color = "SIR - True (Random) Params"),
            linewidth = 0.8) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper),
              fill = "#D32F2F", alpha = 0.20) +
  geom_line(aes(y = pred_med, color = "SIR - BiLSTM Predicted"),
            linewidth = 0.8) +
  facet_wrap(~ panel, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c(
    "SIR - True (Random) Params" = "#1976D2",
    "SIR - BiLSTM Predicted"     = "#D32F2F")) +
  labs(title    = "Part 2 - Random Parameters: True vs BiLSTM-Predicted SIR Curves",
       subtitle = paste0("Shaded bands = 95% CI across ", NSIMS,
                         " runs | BiLSTM input = median of true runs"),
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

