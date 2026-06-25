# =============================================================================
#  SEIR Calibration Validation
#
#  Part 1 — Test set: run ModelSEIRCONN with actual vs BiLSTM-predicted
#            parameters, compare incidence curves, compute per-day MAE
#
#  Part 2 — Random parameters pipeline:
#            sample random params -> simulate incidence -> calibrate_seir()
#            -> simulate with predicted params -> plot comparison
#
#  SEIR known inputs : n, recov, incub, prevalence
#  SEIR predicted    : beta, R0
#  ModelSEIRCONN     : contact_rate=1, transmission_rate=beta
# =============================================================================

library(epiworldR)
library(ggplot2)
library(dplyr)
library(tidyr)

# ── Paths ─────────────────────────────────────────────────────────────────────
MODEL_DIR <- "~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction/model"
DATA_DIR  <- "~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction"

NDAYS    <- 365
N_SAMPLE <- 20
WINDOW   <- "late_365d"
CONTACT_RATE=1
init_bilstm_model(MODEL_DIR)

# =============================================================================
# PART 1 — Test set: actual vs predicted incidence
# =============================================================================

cat("\n-- Part 1: loading test data --\n")

actual    <- read.csv(path.expand(file.path(DATA_DIR, "test_actual_parameters.csv")))
preds_all <- read.csv(path.expand(file.path(DATA_DIR, "test_bilstm_predictions.csv")))

# FIX 1: variable must be np (not p), used as np$load below
np      <- reticulate::import("numpy")
inc_raw <- np$load(path.expand(file.path(DATA_DIR, "test_incidence_raw.npy")))

preds <- preds_all %>% filter(window == WINDOW)
cat(sprintf("  Val sims: %d  |  Window: %s\n", nrow(actual), WINDOW))

# ── Helper: run ModelSEIRCONN ─────────────────────────────────────────────────
run_seir_incidence <- function(n, prevalence, beta, incub, recov,
                               ndays = NDAYS, seed = 1L) {
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
  run(model, ndays = ndays + 1L, seed = seed)   # +1 to compensate for dropping day 0

  inc <- plot_incidence(model, plot = FALSE)[, 1]
  inc <- inc[-1]                                 # drop day 0 initialization spike
  if (length(inc) < ndays) inc <- c(inc, rep(0L, ndays - length(inc)))
  as.numeric(inc[1:ndays])
}

# ── Sample test simulations ───────────────────────────────────────────────────
set.seed(42)
sample_ids <- sample(intersect(preds$sim_idx, actual$sim_idx), N_SAMPLE)
cat(sprintf("Running SEIR simulations for %d sample sims...\n", N_SAMPLE))

curves_list  <- vector("list", N_SAMPLE)
mae_day_list <- vector("list", N_SAMPLE)

for (k in seq_along(sample_ids)) {
  sid <- sample_ids[k]
  cat(sprintf("  [%d/%d] sim_idx = %d\n", k, N_SAMPLE, sid))

  act <- actual %>% filter(sim_idx == sid)
  prd <- preds  %>% filter(sim_idx == sid)

  n          <- act$n[1]
  recov      <- act$recov[1]
  incub      <- act$incub[1]
  prevalence <- act$prevalence[1]

  mat_row <- which(actual$sim_idx == sid)
  obs_inc <- as.numeric(inc_raw[mat_row, ])

  act_inc <- run_seir_incidence(
    n = n, prevalence = prevalence,
    beta = act$beta[1], incub = incub, recov = recov, seed = sid)

  pred_inc <- run_seir_incidence(
    n = n, prevalence = prevalence,
    beta = prd$beta_pred[1], incub = incub, recov = recov, seed = sid)

  mae_per_day <- abs(act_inc - pred_inc)

  curves_list[[k]] <- data.frame(
    day     = 1:NDAYS,
    sim_idx = sid,
    panel   = sprintf(
      "sim %d\nbeta: %.3f->%.3f | R0: %.2f->%.2f | incub=%.1fd",
      sid,
      act$beta[1], prd$beta_pred[1],
      act$R0[1],   prd$R0_pred[1],
      incub
    ),
    Observed_ABM          = obs_inc,
    SEIR_Actual_Params    = act_inc,
    SEIR_BiLSTM_Predicted = pred_inc,
    MAE_per_day           = mae_per_day
  )

  mae_day_list[[k]] <- data.frame(
    day     = 1:NDAYS,
    sim_idx = sid,
    mae     = mae_per_day
  )
}

curves_df  <- bind_rows(curves_list)
mae_day_df <- bind_rows(mae_day_list)

# ── Plot 1: incidence curves ──────────────────────────────────────────────────
curves_long <- curves_df %>%
  select(day, sim_idx, panel,
         Observed_ABM, SEIR_Actual_Params, SEIR_BiLSTM_Predicted) %>%
  pivot_longer(
    cols      = c(Observed_ABM, SEIR_Actual_Params, SEIR_BiLSTM_Predicted),
    names_to  = "source",
    values_to = "incidence"
  ) %>%
  mutate(source = factor(source,
                         levels = c("Observed_ABM", "SEIR_Actual_Params", "SEIR_BiLSTM_Predicted"),
                         labels = c("Observed (ABM)", "SEIR - Actual Params", "SEIR - BiLSTM Predicted")))

p1 <- ggplot(curves_long,
             aes(x = day, y = incidence, color = source, linetype = source)) +
  geom_line(linewidth = 0.7, alpha = 0.9) +
  facet_wrap(~ panel, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c(
    "Observed (ABM)"        = "black",
    "SEIR - Actual Params"  = "#1976D2",
    "SEIR - BiLSTM Predicted" = "#D32F2F")) +
  scale_linetype_manual(values = c(
    "Observed (ABM)"          = "solid",
    "SEIR - Actual Params"    = "dashed",
    "SEIR - BiLSTM Predicted" = "dotdash")) +
  labs(title    = paste0("SEIR Incidence - Actual vs BiLSTM Predicted (", WINDOW, ")"),
       subtitle = "incub and recov are known; only beta is predicted",
       x = "Day", y = "Daily Incidence", color = NULL, linetype = NULL) +
  theme_bw(base_size = 9) +
  theme(legend.position  = "bottom",
        legend.key.width = unit(1.8, "cm"),
        strip.text       = element_text(size = 7),
        panel.grid.minor = element_blank(),
        plot.title       = element_text(face = "bold"))

ggsave(path.expand(file.path(DATA_DIR, "part1_seir_incidence_curves.png")),
       p1, width = 15, height = 13, dpi = 150)
print(p1)

# ── Plot 2: per-day MAE ───────────────────────────────────────────────────────
mae_avg <- mae_day_df %>%
  group_by(day) %>%
  summarise(mean_mae = mean(mae), sd_mae = sd(mae), .groups = "drop")

p2 <- ggplot(mae_avg, aes(x = day)) +
  geom_ribbon(aes(ymin = pmax(mean_mae - sd_mae, 0),
                  ymax = mean_mae + sd_mae),
              fill = "#D32F2F", alpha = 0.15) +
  geom_line(aes(y = mean_mae), color = "#D32F2F", linewidth = 0.9) +
  labs(title    = "Per-day MAE: Actual-Param vs Predicted-Param SEIR",
       subtitle = sprintf("Mean +/- SD across %d sampled simulations", N_SAMPLE),
       x = "Day", y = "MAE (incidence)") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

ggsave(path.expand(file.path(DATA_DIR, "part1_seir_per_day_mae.png")),
       p2, width = 10, height = 5, dpi = 150)
print(p2)

cat("\n-- Per-day MAE summary --\n")
mae_avg %>%
  summarise(overall_mean = mean(mean_mae),
            peak_day     = day[which.max(mean_mae)],
            peak_mae     = max(mean_mae)) %>%
  with({
    cat(sprintf("  Overall mean MAE : %.2f\n", overall_mean))
    cat(sprintf("  Worst day        : day %d  (MAE = %.2f)\n", peak_day, peak_mae))
  })

# =============================================================================
# PART 2 — Random parameters pipeline
# =============================================================================

cat("\n-- Part 2: random parameters pipeline --\n")

N_RANDOM <- 9
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

# FIX: loop variable renamed to rp to avoid clash with np (numpy) above
for (k in seq_len(N_RANDOM)) {
  rp <- random_params[k, ]
  cat(sprintf("\n  [%d/%d]  n=%d  beta=%.3f  recov=%.3f  incub=%.1f  R0=%.2f\n",
              k, N_RANDOM, rp$n, rp$beta, rp$recov, rp$incub, rp$R0))

  true_inc <- run_seir_incidence(
    n = rp$n, prevalence = rp$prevalence,
    beta = rp$beta, incub = rp$incub, recov = rp$recov,
    ndays = NDAYS, seed = k)

  predicted <- calibrate_seir(
    daily_cases     = true_inc,
    population_size = rp$n,
    recovery_rate   = rp$recov,
    incubation_days = rp$incub
  )
  cat(sprintf("    True      : beta=%.4f  R0=%.4f\n", rp$beta, rp$R0))
  cat(sprintf("    Predicted : beta=%.4f  R0=%.4f\n",
              predicted["beta"], predicted["R0"]))

  pred_inc <- run_seir_incidence(
    n = rp$n, prevalence = rp$prevalence,
    beta  = predicted["beta"],
    incub = rp$incub,
    recov = rp$recov,
    ndays = NDAYS, seed = k)

  random_curves[[k]] <- data.frame(
    day   = 1:NDAYS,
    panel = sprintf(
      "sim %d | R0: %.2f->%.2f\nbeta: %.3f->%.3f | incub=%.1fd",
      k, rp$R0, predicted["R0"],
      rp$beta, predicted["beta"], rp$incub
    ),
    True_Params      = true_inc,
    Predicted_Params = pred_inc
  )
}

random_df <- bind_rows(random_curves)

# ── Plot 3: random pipeline comparison ───────────────────────────────────────
random_long <- random_df %>%
  pivot_longer(cols = c(True_Params, Predicted_Params),
               names_to = "source", values_to = "incidence") %>%
  mutate(source = factor(source,
                         levels = c("True_Params", "Predicted_Params"),
                         labels = c("SEIR - True (Random) Params", "SEIR - BiLSTM Predicted")))

p3 <- ggplot(random_long,
             aes(x = day, y = incidence, color = source, linetype = source)) +
  geom_line(linewidth = 0.8, alpha = 0.9) +
  facet_wrap(~ panel, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c(
    "SEIR - True (Random) Params" = "#1976D2",
    "SEIR - BiLSTM Predicted"     = "#D32F2F")) +
  scale_linetype_manual(values = c(
    "SEIR - True (Random) Params" = "solid",
    "SEIR - BiLSTM Predicted"     = "dashed")) +
  labs(title    = "Part 2 - Random Parameters: True vs BiLSTM-Predicted SEIR Curves",
       subtitle = "Blue = SEIR with random true params | Red = SEIR with BiLSTM output",
       x = "Day", y = "Daily Incidence", color = NULL, linetype = NULL) +
  theme_bw(base_size = 9) +
  theme(legend.position  = "bottom",
        legend.key.width = unit(1.8, "cm"),
        strip.text       = element_text(size = 7),
        panel.grid.minor = element_blank(),
        plot.title       = element_text(face = "bold"))

ggsave(path.expand(file.path(DATA_DIR, "part2_seir_random_pipeline.png")),
       p3, width = 15, height = 13, dpi = 150)
print(p3)

cat("\n-- Part 2 summary --\n")
random_df %>%
  group_by(panel) %>%
  summarise(
    mae_incidence  = mean(abs(True_Params - Predicted_Params)),
    peak_true      = max(True_Params),
    peak_predicted = max(Predicted_Params),
    peak_day_true  = which.max(True_Params),
    peak_day_pred  = which.max(Predicted_Params),
    .groups = "drop"
  ) %>%
  mutate(peak_day_err = abs(peak_day_true - peak_day_pred)) %>%
  print()

cat("\nSaved:\n")
cat("  part1_seir_incidence_curves.png\n")
cat("  part1_seir_per_day_mae.png\n")
cat("  part2_seir_random_pipeline.png\n")
