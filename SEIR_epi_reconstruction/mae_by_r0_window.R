# =============================================================================
#  mae_by_r0_window.R
#
#  Per-day incidence MAE grouped by true R0 bin and observation window length.
#
#  Layout: facet_grid(R0 bin ~ regime), colored lines per window length.
#          Y-scales are FIXED across all panels so trends are comparable.
#          Day 1 excluded from plots.
#
#  R0 bins:
#    "R0 ~ 1.5"  : true R0 in [1.2, 1.8)
#    "R0 ~ 2.0"  : true R0 in [1.8, 2.5)
#    "R0 ~ 3.0"  : true R0 in [2.5, 3.8)
#
#  Windows shown: 15d, 30d, 60d, 180d  x  early / mid / late
#
#  Usage: Rscript mae_by_r0_window.R
#  Prerequisite: generate_tuned_predictions.py must have been run first.
# =============================================================================

suppressPackageStartupMessages({
  library(epiworldR)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
DATA_DIR     <- path.expand("~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction")
CONTACT_RATE <- 1
NDAYS        <- 365L
N_PER_BIN    <- 30     # sims sampled from each R0 bin

# All three regimes x four lengths = 12 windows
WINDOWS_SHOW <- c(
  "early_015d", "early_030d", "early_060d", "early_180d",
  "mid_015d",   "mid_030d",   "mid_060d",   "mid_180d",
  "late_015d",  "late_030d",  "late_060d",  "late_180d"
)
WIN_LABELS <- c(
  "015d" = "15 days", "030d" = "30 days",
  "060d" = "60 days", "180d" = "180 days"
)

R0_BINS <- list(
  "R0 ~ 1.5" = c(1.2, 1.8),
  "R0 ~ 2.0" = c(1.8, 2.5),
  "R0 ~ 3.0" = c(2.5, 3.8)
)

set.seed(42)

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
preds_file <- file.path(DATA_DIR, "test_bilstm_predictions_tuned.csv")
if (!file.exists(preds_file)) {
  # fall back to old predictions if tuned not yet generated
  preds_file <- file.path(DATA_DIR, "test_bilstm_predictions.csv")
  cat("NOTE: using old predictions (run generate_tuned_predictions.py for tuned model)\n")
}

actual    <- read.csv(file.path(DATA_DIR, "test_actual_parameters.csv"))
preds_all <- read.csv(preds_file)

cat(sprintf("Loaded %d test sims, %d prediction rows\n",
            nrow(actual), nrow(preds_all)))

# Keep only windows of interest
preds_all <- preds_all |> filter(window %in% WINDOWS_SHOW)

# ---------------------------------------------------------------------------
# Assign R0 bins
# ---------------------------------------------------------------------------
assign_bin <- function(R0_vals) {
  out <- rep(NA_character_, length(R0_vals))
  for (nm in names(R0_BINS)) {
    lo <- R0_BINS[[nm]][1]; hi <- R0_BINS[[nm]][2]
    out[R0_vals >= lo & R0_vals < hi] <- nm
  }
  out
}

actual <- actual |> mutate(r0_bin = assign_bin(R0))
cat("Sims per R0 bin:\n")
print(table(actual$r0_bin, useNA = "ifany"))

# ---------------------------------------------------------------------------
# Sample sims from each bin
# ---------------------------------------------------------------------------
sampled_ids <- actual |>
  filter(!is.na(r0_bin)) |>
  group_by(r0_bin) |>
  slice_sample(n = N_PER_BIN) |>
  ungroup() |>
  select(sim_idx, r0_bin)

cat(sprintf("\nUsing %d sims total (%d per bin)\n", nrow(sampled_ids), N_PER_BIN))

# ---------------------------------------------------------------------------
# Helper: single deterministic SEIR run -> incidence vector (days 1-365)
# ---------------------------------------------------------------------------
run_seir_single <- function(n, prevalence, beta, incub, recov, seed = 1L) {
  n_int <- max(as.integer(round(n)), 10L)
  prev  <- max(min(prevalence, 1.0), 1.0 / n_int)
  model <- ModelSEIRCONN(
    name              = "sim",
    n                 = n_int,
    prevalence        = prev,
    contact_rate      = CONTACT_RATE,
    transmission_rate = max(beta / CONTACT_RATE, 1e-6),
    incubation_days   = max(incub, 1.0),
    recovery_rate     = max(recov, 1e-6)
  )
  verbose_off(model)
  run(model, ndays = NDAYS + 1L, seed = as.integer(seed))
  inc <- plot_incidence(model, plot = FALSE)[, 2]
  inc <- inc[-1]   # drop day-0 init spike
  if (length(inc) < NDAYS) inc <- c(inc, rep(0L, NDAYS - length(inc)))
  as.numeric(inc[1:NDAYS])
}

# ---------------------------------------------------------------------------
# Main loop: for each sampled sim x window, compute per-day |act - pred|
# ---------------------------------------------------------------------------
n_total <- nrow(sampled_ids)
rows    <- vector("list", n_total * length(WINDOWS_SHOW))
idx     <- 1L

cat("\nRunning SEIR simulations...\n")

for (i in seq_len(n_total)) {
  sid     <- sampled_ids$sim_idx[i]
  r0_bin  <- sampled_ids$r0_bin[i]
  act     <- actual |> filter(sim_idx == sid)

  if ((i %% 10) == 0)
    cat(sprintf("  [%d/%d]\n", i, n_total))

  act_inc <- tryCatch(
    run_seir_single(act$n, act$prevalence, act$beta, act$incub, act$recov, seed = sid),
    error = function(e) NULL
  )
  if (is.null(act_inc)) next

  for (win_tag in WINDOWS_SHOW) {
    prd <- preds_all |> filter(sim_idx == sid, window == win_tag)
    if (nrow(prd) == 0) next

    pred_inc <- tryCatch(
      run_seir_single(act$n, act$prevalence, prd$beta_pred[1],
                      act$incub, act$recov, seed = sid),
      error = function(e) NULL
    )
    if (is.null(pred_inc)) next

    rows[[idx]] <- data.frame(
      sim_idx = sid,
      r0_bin  = r0_bin,
      window  = win_tag,
      day     = 1:NDAYS,
      mae     = abs(act_inc - pred_inc),
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L
  }
}

mae_df <- bind_rows(rows[seq_len(idx - 1L)])
cat(sprintf("Done. %d sim x window x day rows.\n", nrow(mae_df)))

# ---------------------------------------------------------------------------
# Average per-day MAE across sims in each (R0 bin, window)
# ---------------------------------------------------------------------------
mae_avg <- mae_df |>
  filter(day > 1) |>
  mutate(
    regime  = sub("_.*", "", window),                        # early/mid/late
    len_tag = sub("^[a-z]+_", "", window)                    # 015d/030d/...
  ) |>
  group_by(r0_bin, regime, len_tag, day) |>
  summarise(
    mean_mae = mean(mae, na.rm = TRUE),
    sd_mae   = sd(mae,   na.rm = TRUE),
    n_sims   = n(),
    .groups  = "drop"
  ) |>
  mutate(
    r0_bin  = factor(r0_bin, levels = names(R0_BINS)),
    regime  = factor(regime, levels = c("early", "mid", "late")),
    win_lab = factor(unname(WIN_LABELS[len_tag]),
                     levels = c("15 days", "30 days", "60 days", "180 days"))
  )

# ---------------------------------------------------------------------------
# Plot: fixed Y-scale across all R0 panels
# ---------------------------------------------------------------------------
win_colors <- c(
  "15 days"  = "#E65100",
  "30 days"  = "#1565C0",
  "60 days"  = "#2E7D32",
  "180 days" = "#6A1B9A"
)

p <- ggplot(mae_avg, aes(x = day, color = win_lab, fill = win_lab)) +
  geom_ribbon(aes(ymin = pmax(mean_mae - sd_mae, 0),
                  ymax = mean_mae + sd_mae),
              alpha = 0.12, color = NA) +
  geom_line(aes(y = mean_mae), linewidth = 0.9) +
  facet_grid(r0_bin ~ regime, scales = "fixed") +
  scale_color_manual(values = win_colors, name = "Window length") +
  scale_fill_manual(values  = win_colors, name = "Window length") +
  labs(
    title    = "Per-day Incidence MAE by R0 Group, Regime, and Window Length",
    subtitle = sprintf(
      "Mean ± SD across ~%d sims per R0 bin | Days 2-365 | Fixed Y-scale",
      N_PER_BIN
    ),
    x = "Day",
    y = "Mean |Actual incidence − Predicted incidence|"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "bottom",
    legend.key.width = unit(1.5, "cm"),
    strip.text       = element_text(size = 11, face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold")
  )

out_path <- file.path(DATA_DIR, "mae_by_r0_window.png")
ggsave(out_path, p, width = 15, height = 10, dpi = 150)
cat(sprintf("\nSaved: %s\n", out_path))
print(p)

# ---------------------------------------------------------------------------
# Also save a summary table
# ---------------------------------------------------------------------------
summary_tbl <- mae_avg |>
  group_by(r0_bin, regime, win_lab) |>
  summarise(
    overall_mae  = mean(mean_mae),
    peak_mae_day = day[which.max(mean_mae)],
    peak_mae     = max(mean_mae),
    .groups = "drop"
  ) |>
  mutate(across(where(is.numeric), \(x) round(x, 3)))

cat("\n=== Summary: Mean MAE per R0 bin x window ===\n")
print(summary_tbl, n = 40)

write.csv(summary_tbl,
          file.path(DATA_DIR, "mae_by_r0_window_summary.csv"),
          row.names = FALSE)
cat("Saved: mae_by_r0_window_summary.csv\n")

# ---------------------------------------------------------------------------
# Parameter MAE: beta and R0 errors directly from predictions
# ---------------------------------------------------------------------------
param_mae <- preds_all |>
  filter(sim_idx %in% sampled_ids$sim_idx) |>
  left_join(sampled_ids, by = "sim_idx") |>
  mutate(
    err_beta = abs(beta_pred - beta_true),
    err_R0   = abs(R0_pred   - R0_true),
    regime   = factor(sub("_.*", "", window), levels = c("early", "mid", "late")),
    len_tag  = sub("^[a-z]+_", "", window),
    win_lab  = factor(unname(WIN_LABELS[len_tag]),
                      levels = c("15 days", "30 days", "60 days", "180 days")),
    r0_bin   = factor(r0_bin, levels = names(R0_BINS))
  ) |>
  group_by(r0_bin, regime, win_lab) |>
  summarise(
    mae_beta = mean(err_beta, na.rm = TRUE),
    mae_R0   = mean(err_R0,   na.rm = TRUE),
    sd_beta  = sd(err_beta,   na.rm = TRUE),
    sd_R0    = sd(err_R0,     na.rm = TRUE),
    .groups  = "drop"
  )

param_long <- bind_rows(
  param_mae |> mutate(param = "MAE \u03b2 (beta)", mae = mae_beta, sd = sd_beta),
  param_mae |> mutate(param = "MAE R\u2080",        mae = mae_R0,   sd = sd_R0)
) |>
  select(r0_bin, regime, win_lab, param, mae, sd) |>
  mutate(param = factor(param, levels = c("MAE \u03b2 (beta)", "MAE R\u2080")))

regime_colors <- c(early = "#E65100", mid = "#1565C0", late = "#2E7D32")

p_param <- ggplot(param_long,
                  aes(x = win_lab, y = mae, color = regime, group = regime)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = pmax(mae - sd, 0), ymax = mae + sd),
                width = 0.2, alpha = 0.6) +
  facet_grid(r0_bin ~ param, scales = "free_y") +
  scale_color_manual(values = regime_colors, name = "Regime") +
  labs(
    title    = "Parameter MAE (\u03b2 and R\u2080) by R0 Group, Regime, and Window Length",
    subtitle = sprintf("Mean \u00b1 SD across ~%d sims per R0 bin", N_PER_BIN),
    x = "Window length",
    y = "Mean |Predicted \u2212 True|"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "bottom",
    strip.text       = element_text(size = 11, face = "bold"),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(angle = 20, hjust = 1),
    plot.title       = element_text(face = "bold")
  )

out_param <- file.path(DATA_DIR, "mae_by_r0_window_params.png")
ggsave(out_param, p_param, width = 10, height = 10, dpi = 150)
cat(sprintf("Saved: %s\n", out_param))

# ---------------------------------------------------------------------------
# Combined plot: incidence MAE (top) + parameter MAE (bottom)
# ---------------------------------------------------------------------------
p_combined <- p / p_param +
  plot_layout(heights = c(2, 1.5)) +
  plot_annotation(
    title   = "BiLSTM Calibration Performance by R0 Group",
    subtitle = sprintf(
      "Top: per-day incidence MAE (Days 2-365) | Bottom: parameter errors | ~%d sims per R0 bin",
      N_PER_BIN
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10)
    )
  )

out_combined <- file.path(DATA_DIR, "mae_by_r0_window_combined.png")
ggsave(out_combined, p_combined, width = 15, height = 18, dpi = 150)
cat(sprintf("Saved: %s\n", out_combined))
