# =============================================================================
#  SEIR Calibration Validation — All Windows + Plot Folder Version
#
#  Part 1:
#    For every available prediction window:
#      - run ModelSEIRCONN with actual parameters
#      - run ModelSEIRCONN with BiLSTM-predicted parameters
#      - compare incidence curves
#      - compute per-day MAE
#      - save plots inside DATA_DIR/plots/<window>/
#
#  Part 2:
#    Random parameter pipeline
#    - simulate random true SEIR curve
#    - calibrate_seir()
#    - simulate with predicted parameters
#    - save plots inside DATA_DIR/plots/random_pipeline/
#
#  SEIR known inputs : n, recov, incub, prevalence
#  SEIR predicted    : beta, R0
#  ModelSEIRCONN     : contact_rate = 1, transmission_rate = beta
# =============================================================================

library(epiworldR)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(reticulate)
library(stringr)

# =============================================================================
# Paths and settings
# =============================================================================

MODEL_DIR <- "~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction/model"
DATA_DIR  <- "~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction"

PLOT_DIR <- path.expand(file.path(DATA_DIR, "plots"))

dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

NDAYS        <- 365
N_SAMPLE     <- 9     # readable plots; change to 20 if you want more panels
N_RANDOM     <- 9
CONTACT_RATE <- 1

init_bilstm_model(MODEL_DIR)

# =============================================================================
# Global plot theme
# =============================================================================

big_theme <- theme_bw(base_size = 15) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 13),
    legend.title = element_blank(),
    legend.key.width = unit(2.5, "cm"),

    strip.text = element_text(size = 11, face = "bold"),
    strip.background = element_rect(fill = "grey95", color = "grey70"),

    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 13),

    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),

    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25),

    plot.margin = margin(10, 10, 10, 10)
  )

# =============================================================================
# Helper: clean names for folders/files
# =============================================================================

safe_name <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("[^A-Za-z0-9_\\-]+", "_")
}

# =============================================================================
# Helper: extract numeric length from window name
# examples: early_015d -> 15, mid_180d -> 180, late_365d -> 365
# =============================================================================

get_window_length <- function(window_name) {
  as.integer(str_extract(window_name, "\\d+"))
}

get_window_regime <- function(window_name) {
  str_extract(window_name, "^[A-Za-z]+")
}

# =============================================================================
# Helper: order windows nicely
# =============================================================================

order_windows <- function(windows) {
  data.frame(
    window = windows,
    regime = get_window_regime(windows),
    length = get_window_length(windows),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      regime_order = case_when(
        regime == "early" ~ 1L,
        regime == "mid"   ~ 2L,
        regime == "late"  ~ 3L,
        TRUE              ~ 99L
      )
    ) %>%
    arrange(regime_order, length) %>%
    pull(window)
}

# =============================================================================
# Helper: run ModelSEIRCONN and return daily incidence
# =============================================================================

run_seir_incidence <- function(n, prevalence, beta, incub, recov,
                               ndays = NDAYS, seed = 1L) {

  n <- max(as.integer(round(n)), 10L)
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

  run(
    model,
    ndays = ndays + 1L,
    seed  = as.integer(seed)
  )

  inc <- plot_incidence(model, plot = FALSE)[, 2]

  # Drop day 0 initialization spike
  inc <- inc[-1]

  if (length(inc) < ndays) {
    inc <- c(inc, rep(0L, ndays - length(inc)))
  }

  as.numeric(inc[1:ndays])
}

# =============================================================================
# Load test data
# =============================================================================

cat("\n-- Loading test data --\n")

actual <- read.csv(
  path.expand(file.path(DATA_DIR, "test_actual_parameters.csv"))
)

preds_all <- read.csv(
  path.expand(file.path(DATA_DIR, "test_bilstm_predictions.csv"))
)

if (!"window" %in% names(preds_all)) {
  stop("The file test_bilstm_predictions.csv must contain a column named 'window'.")
}

np <- reticulate::import("numpy")

inc_raw <- np$load(
  path.expand(file.path(DATA_DIR, "test_incidence_raw.npy"))
)

inc_raw <- as.matrix(inc_raw)

WINDOWS <- order_windows(unique(preds_all$window))

cat("\nDetected windows:\n")
print(WINDOWS)

# =============================================================================
# Select the same sample simulations across all windows when possible
# =============================================================================

common_ids_all <- Reduce(
  intersect,
  split(preds_all$sim_idx, preds_all$window)
)

common_ids_all <- intersect(common_ids_all, actual$sim_idx)

if (length(common_ids_all) == 0) {
  stop("No common sim_idx values across all windows and actual parameter file.")
}

set.seed(42)

N_SAMPLE_USE <- min(N_SAMPLE, length(common_ids_all))

sample_ids_global <- sample(common_ids_all, N_SAMPLE_USE)

cat(sprintf("\nUsing %d shared sample simulations across windows:\n", N_SAMPLE_USE))
print(sample_ids_global)

# =============================================================================
# Cache actual-parameter simulations and observed ABM incidence
# These do not depend on prediction window, so we compute them once.
# =============================================================================

cat("\n-- Precomputing observed and actual-parameter SEIR curves --\n")

actual_cache <- list()

for (sid in sample_ids_global) {

  cat(sprintf("  Actual cache for sim_idx = %d\n", sid))

  act <- actual %>%
    filter(sim_idx == sid)

  if (nrow(act) == 0) {
    warning(sprintf("Skipping sim_idx %d because actual row is missing.", sid))
    next
  }

  n          <- act$n[1]
  recov      <- act$recov[1]
  incub      <- act$incub[1]
  prevalence <- act$prevalence[1]

  mat_row <- which(actual$sim_idx == sid)

  obs_inc <- as.numeric(inc_raw[mat_row, ])

  act_inc <- run_seir_incidence(
    n          = n,
    prevalence = prevalence,
    beta       = act$beta[1],
    incub      = incub,
    recov      = recov,
    seed       = sid
  )

  actual_cache[[as.character(sid)]] <- list(
    act        = act,
    obs_inc    = obs_inc,
    act_inc    = act_inc,
    n          = n,
    recov      = recov,
    incub      = incub,
    prevalence = prevalence
  )
}

# =============================================================================
# PART 1 — Loop over all windows
# =============================================================================

all_window_summaries <- list()
all_window_mae_days  <- list()

for (WINDOW in WINDOWS) {

  cat("\n============================================================\n")
  cat(sprintf("Running window: %s\n", WINDOW))
  cat("============================================================\n")

  WINDOW_SAFE <- safe_name(WINDOW)

  WINDOW_DIR <- file.path(PLOT_DIR, WINDOW_SAFE)

  INDIV_DIR <- file.path(WINDOW_DIR, "individual_incidence_plots")

  dir.create(WINDOW_DIR, showWarnings = FALSE, recursive = TRUE)
  dir.create(INDIV_DIR, showWarnings = FALSE, recursive = TRUE)

  preds <- preds_all %>%
    filter(window == WINDOW)

  common_ids <- intersect(preds$sim_idx, sample_ids_global)

  if (length(common_ids) == 0) {
    warning(sprintf("Skipping window %s because no sampled sim_idx values are available.", WINDOW))
    next
  }

  curves_list  <- vector("list", length(common_ids))
  mae_day_list <- vector("list", length(common_ids))

  for (k in seq_along(common_ids)) {

    sid <- common_ids[k]

    cat(sprintf("  [%d/%d] sim_idx = %d\n", k, length(common_ids), sid))

    prd <- preds %>%
      filter(sim_idx == sid)

    if (nrow(prd) == 0) {
      warning(sprintf("Skipping sim_idx %d in window %s because prediction row is missing.", sid, WINDOW))
      next
    }

    cache <- actual_cache[[as.character(sid)]]

    if (is.null(cache)) {
      warning(sprintf("Skipping sim_idx %d because actual cache is missing.", sid))
      next
    }

    act <- cache$act

    pred_inc <- run_seir_incidence(
      n          = cache$n,
      prevalence = cache$prevalence,
      beta       = prd$beta_pred[1],
      incub      = cache$incub,
      recov      = cache$recov,
      seed       = sid
    )

    mae_per_day <- abs(cache$act_inc - pred_inc)

    curves_list[[k]] <- data.frame(
      day     = 1:NDAYS,
      sim_idx = sid,
      window  = WINDOW,
      regime  = get_window_regime(WINDOW),
      win_len = get_window_length(WINDOW),
      panel   = sprintf(
        "sim %d\nbeta %.3f -> %.3f | R0 %.2f -> %.2f | incub %.1fd",
        sid,
        act$beta[1],
        prd$beta_pred[1],
        act$R0[1],
        prd$R0_pred[1],
        cache$incub
      ),
      Observed_ABM          = cache$obs_inc,
      SEIR_Actual_Params    = cache$act_inc,
      SEIR_BiLSTM_Predicted = pred_inc,
      MAE_per_day           = mae_per_day
    )

    mae_day_list[[k]] <- data.frame(
      day     = 1:NDAYS,
      sim_idx = sid,
      window  = WINDOW,
      regime  = get_window_regime(WINDOW),
      win_len = get_window_length(WINDOW),
      mae     = mae_per_day
    )
  }

  curves_df <- bind_rows(curves_list)
  mae_day_df <- bind_rows(mae_day_list)

  if (nrow(curves_df) == 0) {
    warning(sprintf("No curves produced for window %s.", WINDOW))
    next
  }

  # ===========================================================================
  # Prepare long format
  # ===========================================================================

  curves_long <- curves_df %>%
    select(
      day,
      sim_idx,
      window,
      regime,
      win_len,
      panel,
      Observed_ABM,
      SEIR_Actual_Params,
      SEIR_BiLSTM_Predicted
    ) %>%
    pivot_longer(
      cols = c(
        Observed_ABM,
        SEIR_Actual_Params,
        SEIR_BiLSTM_Predicted
      ),
      names_to = "source",
      values_to = "incidence"
    ) %>%
    mutate(
      source = factor(
        source,
        levels = c(
          "Observed_ABM",
          "SEIR_Actual_Params",
          "SEIR_BiLSTM_Predicted"
        ),
        labels = c(
          "Observed ABM",
          "SEIR - Actual Params",
          "SEIR - BiLSTM Predicted"
        )
      )
    )

  # ===========================================================================
  # Plot 1 — Combined incidence curves for this window
  # ===========================================================================

  p1_big <- ggplot(
    curves_long,
    aes(x = day, y = incidence, color = source, linetype = source)
  ) +
    geom_line(linewidth = 1.15, alpha = 0.95) +
    facet_wrap(~ panel, scales = "free_y", ncol = 2) +
    scale_color_manual(values = c(
      "Observed ABM" = "black",
      "SEIR - Actual Params" = "#1565C0",
      "SEIR - BiLSTM Predicted" = "#C62828"
    )) +
    scale_linetype_manual(values = c(
      "Observed ABM" = "solid",
      "SEIR - Actual Params" = "dashed",
      "SEIR - BiLSTM Predicted" = "solid"
    )) +
    scale_y_continuous(labels = comma) +
    labs(
      title = paste0("SEIR Incidence Curves — ", WINDOW),
      subtitle = "Black = observed ABM | Blue dashed = SEIR with actual beta | Red = SEIR with BiLSTM-predicted beta",
      x = "Day",
      y = "Daily incidence"
    ) +
    big_theme

  ggsave(
    filename = file.path(WINDOW_DIR, paste0("part1_incidence_curves_", WINDOW_SAFE, ".png")),
    plot     = p1_big,
    width    = 18,
    height   = 26,
    dpi      = 300
  )

  print(p1_big)

  # ===========================================================================
  # Plot 1B — One plot per simulation for this window
  # ===========================================================================

  for (sid in unique(curves_long$sim_idx)) {

    df_one <- curves_long %>%
      filter(sim_idx == sid)

    p_one <- ggplot(
      df_one,
      aes(x = day, y = incidence, color = source, linetype = source)
    ) +
      geom_line(linewidth = 1.35, alpha = 0.95) +
      scale_color_manual(values = c(
        "Observed ABM" = "black",
        "SEIR - Actual Params" = "#1565C0",
        "SEIR - BiLSTM Predicted" = "#C62828"
      )) +
      scale_linetype_manual(values = c(
        "Observed ABM" = "solid",
        "SEIR - Actual Params" = "dashed",
        "SEIR - BiLSTM Predicted" = "solid"
      )) +
      scale_y_continuous(labels = comma) +
      labs(
        title = paste("SEIR Incidence Comparison —", WINDOW, "— sim", sid),
        subtitle = unique(df_one$panel),
        x = "Day",
        y = "Daily incidence"
      ) +
      big_theme

    ggsave(
      filename = file.path(
        INDIV_DIR,
        paste0("sim_", sid, "_", WINDOW_SAFE, "_incidence.png")
      ),
      plot   = p_one,
      width  = 12,
      height = 7,
      dpi    = 300
    )
  }

  # ===========================================================================
  # Plot 2 — Per-day MAE for this window
  # ===========================================================================

  mae_avg <- mae_day_df %>%
    group_by(day, window, regime, win_len) %>%
    summarise(
      mean_mae = mean(mae, na.rm = TRUE),
      sd_mae   = sd(mae, na.rm = TRUE),
      .groups  = "drop"
    )

  p2_big <- ggplot(mae_avg, aes(x = day)) +
    geom_ribbon(
      aes(
        ymin = pmax(mean_mae - sd_mae, 0),
        ymax = mean_mae + sd_mae
      ),
      fill = "#C62828",
      alpha = 0.18
    ) +
    geom_line(
      aes(y = mean_mae),
      color = "#C62828",
      linewidth = 1.4
    ) +
    scale_y_continuous(labels = comma) +
    labs(
      title = paste0("Per-day MAE — ", WINDOW),
      subtitle = sprintf("Mean ± SD across %d sampled simulations", length(unique(mae_day_df$sim_idx))),
      x = "Day",
      y = "Mean absolute error in daily incidence"
    ) +
    big_theme +
    theme(legend.position = "none")

  ggsave(
    filename = file.path(WINDOW_DIR, paste0("part1_per_day_mae_", WINDOW_SAFE, ".png")),
    plot     = p2_big,
    width    = 13,
    height   = 7,
    dpi      = 300
  )

  print(p2_big)

  # ===========================================================================
  # Window summary
  # ===========================================================================

  mae_summary <- mae_avg %>%
    summarise(
      window       = first(window),
      regime       = first(regime),
      win_len      = first(win_len),
      n_sims       = length(unique(mae_day_df$sim_idx)),
      overall_mean = mean(mean_mae, na.rm = TRUE),
      peak_day     = day[which.max(mean_mae)],
      peak_mae     = max(mean_mae, na.rm = TRUE),
      .groups      = "drop"
    )

  print(mae_summary)

  write.csv(
    mae_summary,
    file.path(WINDOW_DIR, paste0("mae_summary_", WINDOW_SAFE, ".csv")),
    row.names = FALSE
  )

  write.csv(
    mae_avg,
    file.path(WINDOW_DIR, paste0("mae_per_day_", WINDOW_SAFE, ".csv")),
    row.names = FALSE
  )

  all_window_summaries[[WINDOW]] <- mae_summary
  all_window_mae_days[[WINDOW]]  <- mae_avg
}

# =============================================================================
# Save combined summaries across all windows
# =============================================================================

all_window_summary_df <- bind_rows(all_window_summaries)
all_window_mae_day_df <- bind_rows(all_window_mae_days)

write.csv(
  all_window_summary_df,
  file.path(PLOT_DIR, "all_windows_mae_summary.csv"),
  row.names = FALSE
)

write.csv(
  all_window_mae_day_df,
  file.path(PLOT_DIR, "all_windows_mae_per_day.csv"),
  row.names = FALSE
)

# =============================================================================
# Plot 4 — Summary: overall MAE by window length and regime
# =============================================================================

if (nrow(all_window_summary_df) > 0) {

  all_window_summary_df <- all_window_summary_df %>%
    mutate(
      regime = factor(regime, levels = c("early", "mid", "late")),
      window_label = paste0(regime, "_", win_len, "d")
    )

  p_mae_window <- ggplot(
    all_window_summary_df,
    aes(x = win_len, y = overall_mean, color = regime, group = regime)
  ) +
    geom_line(linewidth = 1.2, alpha = 0.95) +
    geom_point(size = 3) +
    scale_x_continuous(
      breaks = sort(unique(all_window_summary_df$win_len))
    ) +
    scale_y_continuous(labels = comma) +
    labs(
      title = "Overall Incidence MAE Across Prediction Windows",
      subtitle = "Lower values mean the BiLSTM-predicted parameter curve is closer to the actual-parameter SEIR curve",
      x = "Window length used by BiLSTM prediction",
      y = "Overall mean MAE",
      color = "Window regime"
    ) +
    big_theme +
    theme(
      legend.title = element_text(size = 13, face = "bold")
    )

  ggsave(
    filename = file.path(PLOT_DIR, "all_windows_overall_mae_by_length.png"),
    plot     = p_mae_window,
    width    = 12,
    height   = 7,
    dpi      = 300
  )

  print(p_mae_window)

  # ---------------------------------------------------------------------------
  # Bar plot version: easier to compare all window labels
  # ---------------------------------------------------------------------------

  p_mae_bar <- all_window_summary_df %>%
    mutate(
      window_label = factor(
        window,
        levels = WINDOWS
      )
    ) %>%
    ggplot(
      aes(x = window_label, y = overall_mean, fill = regime)
    ) +
    geom_col(alpha = 0.85) +
    scale_y_continuous(labels = comma) +
    labs(
      title = "Overall MAE for Each Window",
      subtitle = "Each bar summarizes the mean per-day incidence MAE for one prediction window",
      x = "Window",
      y = "Overall mean MAE",
      fill = "Regime"
    ) +
    big_theme +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_text(size = 13, face = "bold")
    )

  ggsave(
    filename = file.path(PLOT_DIR, "all_windows_overall_mae_barplot.png"),
    plot     = p_mae_bar,
    width    = 15,
    height   = 8,
    dpi      = 300
  )

  print(p_mae_bar)
}

# =============================================================================
# PART 2 — Random parameters pipeline
# =============================================================================

cat("\n-- Part 2: random parameters pipeline --\n")

RANDOM_DIR <- file.path(PLOT_DIR, "random_pipeline")
RANDOM_INDIV_DIR <- file.path(RANDOM_DIR, "individual_random_pipeline_plots")

dir.create(RANDOM_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RANDOM_INDIV_DIR, showWarnings = FALSE, recursive = TRUE)

set.seed(123)

random_params <- data.frame(
  n          = sample(5000:10000, N_RANDOM, replace = TRUE),
  prevalence = runif(N_RANDOM, min = 0.03, max = 0.35),
  beta       = runif(N_RANDOM, min = 0.05, max = 0.50),
  recov      = runif(N_RANDOM, min = 0.07, max = 0.25),
  incub      = runif(N_RANDOM, min = 2.0,  max = 10.0)
)

random_params$R0 <- random_params$beta / random_params$recov

cat("\nRandom parameters sampled:\n")
print(round(random_params, 4))

random_curves <- vector("list", N_RANDOM)

for (k in seq_len(N_RANDOM)) {

  rp <- random_params[k, ]

  cat(sprintf(
    "\n  [%d/%d] n=%d beta=%.3f recov=%.3f incub=%.1f R0=%.2f\n",
    k,
    N_RANDOM,
    rp$n,
    rp$beta,
    rp$recov,
    rp$incub,
    rp$R0
  ))

  true_inc <- run_seir_incidence(
    n          = rp$n,
    prevalence = rp$prevalence,
    beta       = rp$beta,
    incub      = rp$incub,
    recov      = rp$recov,
    ndays      = NDAYS,
    seed       = k
  )

  predicted <- calibrate_seir(
    daily_cases     = true_inc,
    population_size = rp$n,
    recovery_rate   = rp$recov,
    incubation_days = rp$incub
  )

  beta_pred <- as.numeric(predicted["beta"])
  R0_pred   <- as.numeric(predicted["R0"])

  cat(sprintf("    True      : beta=%.4f R0=%.4f\n", rp$beta, rp$R0))
  cat(sprintf("    Predicted : beta=%.4f R0=%.4f\n", beta_pred, R0_pred))

  pred_inc <- run_seir_incidence(
    n          = rp$n,
    prevalence = rp$prevalence,
    beta       = beta_pred,
    incub      = rp$incub,
    recov      = rp$recov,
    ndays      = NDAYS,
    seed       = k
  )

  random_curves[[k]] <- data.frame(
    day = 1:NDAYS,
    sim_idx = k,
    panel = sprintf(
      "sim %d | R0 %.2f -> %.2f\nbeta %.3f -> %.3f | incub %.1fd",
      k,
      rp$R0,
      R0_pred,
      rp$beta,
      beta_pred,
      rp$incub
    ),
    True_Params      = true_inc,
    Predicted_Params = pred_inc
  )
}

random_df <- bind_rows(random_curves)

random_long <- random_df %>%
  pivot_longer(
    cols = c(True_Params, Predicted_Params),
    names_to = "source",
    values_to = "incidence"
  ) %>%
  mutate(
    source = factor(
      source,
      levels = c("True_Params", "Predicted_Params"),
      labels = c(
        "SEIR - True Random Params",
        "SEIR - BiLSTM Predicted"
      )
    )
  )

# =============================================================================
# Plot 5 — Random parameter pipeline combined plot
# =============================================================================

p_random_big <- ggplot(
  random_long,
  aes(x = day, y = incidence, color = source, linetype = source)
) +
  geom_line(linewidth = 1.2, alpha = 0.95) +
  facet_wrap(~ panel, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c(
    "SEIR - True Random Params" = "#1565C0",
    "SEIR - BiLSTM Predicted" = "#C62828"
  )) +
  scale_linetype_manual(values = c(
    "SEIR - True Random Params" = "solid",
    "SEIR - BiLSTM Predicted" = "solid"
  )) +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Random SEIR Parameters: True vs BiLSTM-Predicted Curves",
    subtitle = "Blue = SEIR with true random parameters | Red = SEIR with BiLSTM-predicted parameters",
    x = "Day",
    y = "Daily incidence"
  ) +
  big_theme

ggsave(
  filename = file.path(RANDOM_DIR, "part2_random_pipeline_combined.png"),
  plot     = p_random_big,
  width    = 16,
  height   = 18,
  dpi      = 300
)

print(p_random_big)

# =============================================================================
# Plot 5B — One random-pipeline plot per simulation
# =============================================================================

for (sid in unique(random_long$sim_idx)) {

  df_one <- random_long %>%
    filter(sim_idx == sid)

  p_one_random <- ggplot(
    df_one,
    aes(x = day, y = incidence, color = source, linetype = source)
  ) +
    geom_line(linewidth = 1.35, alpha = 0.95) +
    scale_color_manual(values = c(
      "SEIR - True Random Params" = "#1565C0",
      "SEIR - BiLSTM Predicted" = "#C62828"
    )) +
    scale_linetype_manual(values = c(
      "SEIR - True Random Params" = "solid",
      "SEIR - BiLSTM Predicted" = "solid"
    )) +
    scale_y_continuous(labels = comma) +
    labs(
      title = paste("Random Parameter Pipeline — sim", sid),
      subtitle = unique(df_one$panel),
      x = "Day",
      y = "Daily incidence"
    ) +
    big_theme

  ggsave(
    filename = file.path(
      RANDOM_INDIV_DIR,
      paste0("random_sim_", sid, "_pipeline.png")
    ),
    plot   = p_one_random,
    width  = 12,
    height = 7,
    dpi    = 300
  )
}

# =============================================================================
# Part 2 summary
# =============================================================================

part2_summary <- random_df %>%
  group_by(panel) %>%
  summarise(
    mae_incidence  = mean(abs(True_Params - Predicted_Params), na.rm = TRUE),
    peak_true      = max(True_Params, na.rm = TRUE),
    peak_predicted = max(Predicted_Params, na.rm = TRUE),
    peak_day_true  = which.max(True_Params),
    peak_day_pred  = which.max(Predicted_Params),
    .groups = "drop"
  ) %>%
  mutate(
    peak_day_err = abs(peak_day_true - peak_day_pred)
  )

print(part2_summary)

write.csv(
  part2_summary,
  file.path(RANDOM_DIR, "part2_random_pipeline_summary.csv"),
  row.names = FALSE
)

# =============================================================================
# Final saved output message
# =============================================================================

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")

cat("\nAll plots saved under:\n")
cat(PLOT_DIR, "\n")

cat("\nMain outputs:\n")
cat("  plots/<window>/part1_incidence_curves_<window>.png\n")
cat("  plots/<window>/part1_per_day_mae_<window>.png\n")
cat("  plots/<window>/individual_incidence_plots/\n")
cat("  plots/all_windows_mae_summary.csv\n")
cat("  plots/all_windows_mae_per_day.csv\n")
cat("  plots/all_windows_overall_mae_by_length.png\n")
cat("  plots/all_windows_overall_mae_barplot.png\n")
cat("  plots/random_pipeline/part2_random_pipeline_combined.png\n")
cat("  plots/random_pipeline/individual_random_pipeline_plots/\n")
cat("  plots/random_pipeline/part2_random_pipeline_summary.csv\n")

cat("\nFinished successfully.\n")
