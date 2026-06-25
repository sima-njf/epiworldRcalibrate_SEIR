# =============================================================================
#  SIR Calibration Validation
#
#  Part 1 — Test set: run ModelSIRCONN with actual vs BiLSTM-predicted
#            parameters, compare incidence curves, compute per-day MAE
#
#  Part 2 — Random parameters pipeline:
#            sample random params → simulate incidence → calibrate_sir()
#            → simulate with predicted params → plot comparison
# =============================================================================

# install.packages(c("epiworldR", "RcppCNPy", "ggplot2", "dplyr", "tidyr"))
library(epiworldR)
library(RcppCNPy)
library(ggplot2)
library(dplyr)
library(tidyr)

# ── Paths ─────────────────────────────────────────────────────────────────────
MODEL_DIR <- "~/sima/epiworldRcalibrate_SEIR/SIR_epi_reconstruction/model"   # <-- UPDATE
DATA_DIR  <- "~/sima/epiworldRcalibrate_SEIR/SIR_epi_reconstruction"                # <-- UPDATE (test CSVs & npy)

NDAYS     <- 365
N_SAMPLE  <- 20 # number of test sims to plot (perfect for 3×3 grid)
WINDOW    <- "late_365d"   # which prediction window to use from predictions CSV

# ── Load BiLSTM model ─────────────────────────────────────────────────────────

init_bilstm_model(MODEL_DIR)

# =============================================================================
# PART 1 — Test set: actual vs predicted incidence
# =============================================================================

cat("\n── Part 1: loading test data ──\n")
actual    <- read.csv(file.path(DATA_DIR, "test_actual_parameters.csv"))
preds_all <- read.csv(file.path(DATA_DIR, "test_bilstm_predictions.csv"))
np      <- reticulate::import("numpy")
inc_raw <- np$load(path.expand(file.path(DATA_DIR, "test_incidence_raw.npy")))
# inc_raw is now a proper R matrix with correct values
preds <- preds_all %>% filter(window == WINDOW)
cat(sprintf("  Val sims: %d  |  Using window: %s\n", nrow(actual), WINDOW))

# ── Helper: run ModelSIRCONN and return incidence vector ──────────────────────
run_sir_incidence <- function(n, prevalence, ptran, crate, recov,
                              ndays = NDAYS, seed = 1L) {
  n    <- max(as.integer(round(n)), 10L)
  prev <- max(min(prevalence, 1.0), 1.0 / n)

  model <- ModelSIRCONN(
    name              = "sim",
    n                 = n,
    prevalence        = prev,
    contact_rate      = max(crate,  1e-6),
    transmission_rate = max(ptran,  1e-6),
    recovery_rate     = max(recov,  1e-6)
  )
  verbose_off(model)
  run(model, ndays = ndays, seed = seed)

  inc <- plot_incidence(model, plot = FALSE)[, 1]
  # Pad or trim to exactly ndays values
  if (length(inc) < ndays) inc <- c(inc, rep(0L, ndays - length(inc)))
  as.numeric(inc[1:ndays])
}

# ── Sample test simulations ───────────────────────────────────────────────────
set.seed(42)
sample_ids <- sample(intersect(preds$sim_idx, actual$sim_idx), N_SAMPLE)

cat(sprintf("Running simulations for %d sample sims...\n", N_SAMPLE))

# Collect per-sim curves and per-day MAE
curves_list  <- vector("list", N_SAMPLE)
mae_day_list <- vector("list", N_SAMPLE)

for (k in seq_along(sample_ids)) {
  sid <- sample_ids[k]
  cat(sprintf("  [%d/%d] sim_idx = %d\n", k, N_SAMPLE, sid))

  act <- actual %>% filter(sim_idx == sid)
  prd <- preds  %>% filter(sim_idx == sid)

  n          <- act$n[1]
  recov      <- act$recov[1]
  prevalence <- act$prevalence[1]

  # Observed incidence from the saved ABM data
  mat_row  <- which(actual$sim_idx == sid)
  obs_inc  <- as.numeric(inc_raw[mat_row, ])

  # ModelSIRCONN with ACTUAL parameters
  act_inc <- run_sir_incidence(
    n = n, prevalence = prevalence,
    ptran = act$ptran[1], crate = act$crate[1], recov = recov, seed = sid)

  # ModelSIRCONN with BiLSTM PREDICTED parameters
  pred_inc <- run_sir_incidence(
    n = n, prevalence = prevalence,
    ptran = prd$ptran_pred[1], crate = prd$crate_pred[1], recov = recov, seed = sid)

  # Per-day MAE between actual-param SIR and predicted-param SIR
  mae_per_day <- abs(act_inc - pred_inc)

  curves_list[[k]] <- data.frame(
    day     = 1:NDAYS,
    sim_idx = sid,
    panel   = sprintf(
      "sim %d\nptran: %.3f→%.3f | crate: %.2f→%.2f | R0: %.2f→%.2f",
      sid,
      act$ptran[1], prd$ptran_pred[1],
      act$crate[1], prd$crate_pred[1],
      act$R0[1],    prd$R0_pred[1]
    ),
    Observed_ABM          = obs_inc,
    SIR_Actual_Params     = act_inc,
    SIR_BiLSTM_Predicted  = pred_inc,
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

# ── Plot 1: incidence curves ─────────────────────────────────────────────────
curves_long <- curves_df %>%
  select(day, sim_idx, panel, Observed_ABM,
         SIR_Actual_Params, SIR_BiLSTM_Predicted) %>%
  pivot_longer(
    cols      = c(Observed_ABM, SIR_Actual_Params, SIR_BiLSTM_Predicted),
    names_to  = "source",
    values_to = "incidence"
  ) %>%
  mutate(source = factor(source,
                         levels = c("Observed_ABM", "SIR_Actual_Params", "SIR_BiLSTM_Predicted"),
                         labels = c("Observed (ABM)", "SIR — Actual Params", "SIR — BiLSTM Predicted")))

p1 <- ggplot(curves_long,
             aes(x = day, y = incidence, color = source, linetype = source)) +
  geom_line(linewidth = 0.7, alpha = 0.9) +
  facet_wrap(~ panel, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c(
    "Observed (ABM)"         = "black",
    "SIR — Actual Params"    = "#1976D2",
    "SIR — BiLSTM Predicted" = "#D32F2F")) +
  scale_linetype_manual(values = c(
    "Observed (ABM)"         = "solid",
    "SIR — Actual Params"    = "dashed",
    "SIR — BiLSTM Predicted" = "dotdash")) +
  labs(title    = paste0("SIR Incidence — Actual vs BiLSTM Predicted (", WINDOW, ")"),
       x = "Day", y = "Daily Incidence", color = NULL, linetype = NULL) +
  theme_bw(base_size = 9) +
  theme(legend.position  = "bottom",
        legend.key.width = unit(1.8, "cm"),
        strip.text       = element_text(size = 7),
        panel.grid.minor = element_blank(),
        plot.title       = element_text(face = "bold"))

ggsave(file.path(DATA_DIR, "part1_incidence_curves.png"),
       p1, width = 15, height = 13, dpi = 150)
print(p1)

# ── Plot 2: per-day MAE ───────────────────────────────────────────────────────
# Average MAE across all sampled sims
mae_avg <- mae_day_df %>%
  group_by(day) %>%
  summarise(mean_mae = mean(mae), sd_mae = sd(mae), .groups = "drop")

p2 <- ggplot(mae_avg, aes(x = day)) +
  geom_ribbon(aes(ymin = pmax(mean_mae - sd_mae, 0),
                  ymax = mean_mae + sd_mae),
              fill = "#D32F2F", alpha = 0.15) +
  geom_line(aes(y = mean_mae), color = "#D32F2F", linewidth = 0.9) +
  labs(title    = "Per-day MAE between Actual-Param and Predicted-Param SIR",
       subtitle = sprintf("Mean ± SD across %d sampled simulations", N_SAMPLE),
       x = "Day", y = "MAE (incidence)") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

ggsave(file.path(DATA_DIR, "part1_per_day_mae.png"),
       p2, width = 10, height = 5, dpi = 150)
print(p2)

# ── Per-sim per-day MAE table ─────────────────────────────────────────────────
cat("\n── Per-day MAE summary (mean across all sampled sims) ──\n")
mae_summary <- mae_avg %>%
  summarise(
    overall_mean_mae = mean(mean_mae),
    peak_day         = day[which.max(mean_mae)],
    peak_mae         = max(mean_mae)
  )
cat(sprintf("  Overall mean MAE : %.2f\n",  mae_summary$overall_mean_mae))
cat(sprintf("  Worst day        : day %d  (MAE = %.2f)\n",
            mae_summary$peak_day, mae_summary$peak_mae))

# =============================================================================
# PART 2 — Random parameters pipeline
# =============================================================================

cat("\n── Part 2: random parameters pipeline ──\n")

# ── Define plausible SIR parameter ranges ────────────────────────────────────
# (Based on ranges seen in theta_use_sir.csv)
N_RANDOM    <- 15     # number of random sims (3×3 grid)
NDAYS_RAN   <- 365
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

# ── Run pipeline for each random sim ─────────────────────────────────────────
random_curves <- vector("list", N_RANDOM)

for (k in seq_len(N_RANDOM)) {
  p <- random_params[k, ]
  cat(sprintf("\n  [%d/%d]  n=%d  ptran=%.3f  crate=%.2f  recov=%.3f  R0=%.2f\n",
              k, N_RANDOM, p$n, p$ptran, p$crate, p$recov, p$R0))

  # Step 1: simulate incidence with TRUE (random) parameters
  true_inc <- run_sir_incidence(
    n = p$n, prevalence = p$prevalence,
    ptran = p$ptran, crate = p$crate, recov = p$recov,
    ndays = NDAYS_RAN, seed = k)

  # Step 2: feed simulated incidence into BiLSTM
  predicted <- calibrate_sir(
    daily_cases     = true_inc,
    population_size = p$n,
    recovery_rate   = p$recov
  )
  cat(sprintf("    True   : ptran=%.4f  crate=%.4f  R0=%.4f\n",
              p$ptran, p$crate, p$R0))
  cat(sprintf("    Predicted: ptran=%.4f  crate=%.4f  R0=%.4f\n",
              predicted["ptran"], predicted["crate"], predicted["R0"]))

  # Step 3: simulate incidence with PREDICTED parameters
  pred_inc <- run_sir_incidence(
    n = p$n, prevalence = p$prevalence,
    ptran = predicted["ptran"], crate = predicted["crate"], recov = p$recov,
    ndays = NDAYS_RAN, seed = k)

  random_curves[[k]] <- data.frame(
    day   = 1:NDAYS_RAN,
    panel = sprintf(
      "sim %d | R0: %.2f→%.2f\nptran: %.3f→%.3f | crate: %.2f→%.2f",
      k, p$R0, predicted["R0"],
      p$ptran, predicted["ptran"],
      p$crate, predicted["crate"]
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
                         labels = c("SIR — True (Random) Params", "SIR — BiLSTM Predicted")))

p3 <- ggplot(random_long,
             aes(x = day, y = incidence, color = source, linetype = source)) +
  geom_line(linewidth = 0.8, alpha = 0.9) +
  facet_wrap(~ panel, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c(
    "SIR — True (Random) Params" = "#1976D2",
    "SIR — BiLSTM Predicted"     = "#D32F2F")) +
  scale_linetype_manual(values = c(
    "SIR — True (Random) Params" = "solid",
    "SIR — BiLSTM Predicted"     = "dashed")) +
  labs(title    = "Part 2 — Random Parameters: True vs BiLSTM-Predicted SIR Curves",
       subtitle = "Blue = SIR simulated with random true params | Red = SIR with BiLSTM output",
       x = "Day", y = "Daily Incidence", color = NULL, linetype = NULL) +
  theme_bw(base_size = 9) +
  theme(legend.position  = "bottom",
        legend.key.width = unit(1.8, "cm"),
        strip.text       = element_text(size = 7),
        panel.grid.minor = element_blank(),
        plot.title       = element_text(face = "bold"))

ggsave(file.path(DATA_DIR, "part2_random_pipeline.png"),
       p3, width = 15, height = 13, dpi = 150)
print(p3)

# ── Part 2 summary ────────────────────────────────────────────────────────────
cat("\n── Part 2 summary ──\n")
random_df %>%
  group_by(panel) %>%
  summarise(
    mae_incidence   = mean(abs(True_Params - Predicted_Params)),
    peak_true       = max(True_Params),
    peak_predicted  = max(Predicted_Params),
    peak_day_true   = which.max(True_Params),
    peak_day_pred   = which.max(Predicted_Params),
    .groups = "drop"
  ) %>%
  mutate(peak_day_err = abs(peak_day_true - peak_day_pred)) %>%
  print()

cat("\nSaved:\n")
cat("  part1_incidence_curves.png\n")
cat("  part1_per_day_mae.png\n")
cat("  part2_random_pipeline.png\n")
