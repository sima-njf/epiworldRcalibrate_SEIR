# calibrate_real_seir.R
# Calibrates SEIR parameters on two real datasets and plots reconstructions.
#
# Methods:
#   1. Nelder-Mead optimization  (implemented here via optim())
#   2. Bernardo BiLSTM           (pre-computed by bernardo_predict_real.py)
#
# Prerequisites — run in order:
#   Rscript real_data/prepare_utah_covid.R
#   Rscript real_data/prepare_measles.R
#   python  real_data/bernardo_predict_real.py
#
# Output:
#   real_data/plots/covid_reconstruction.png
#   real_data/plots/covid_bilstm_windows.png
#   real_data/plots/measles_reconstruction.png
#   real_data/plots/measles_bilstm_windows.png
#   real_data/real_data_results.csv
#
# Usage:
#   Rscript real_data/calibrate_real_seir.R

suppressPackageStartupMessages({
  library(epiworldR)
  library(data.table)
  library(ggplot2)
  library(dplyr)
})

PROJECT_DIR <- path.expand(
  "~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction"
)
REAL_DIR  <- file.path(PROJECT_DIR, "real_data")
PLOTS_DIR <- file.path(REAL_DIR, "plots")
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)

source(file.path(PROJECT_DIR, "seir_common.R"))
seir_set_libpath()

# =============================================================================
# Nelder-Mead calibration on a real incidence vector
# =============================================================================

nm_calibrate_real <- function(obs_inc, known,
                               n_reps   = 200L,
                               nthreads = 4L,
                               maxit    = 2000L) {
  n_days  <- length(obs_inc)
  obj_fn  <- function(z) {
    theta <- seir_from_unc(z)
    if (!seir_feasible(theta, known)) return(1e10)
    p   <- seir_resolve(theta, known)
    sim <- tryCatch(
      seir_run_single(p, ndays = n_days),
      error = function(e) rep(NA_real_, n_days)
    )
    if (any(is.na(sim)) || all(sim == 0)) return(1e8)
    seir_smape(obs_inc, sim)
  }

  # Try multiple starting points spanning the beta search space
  init_betas <- c(0.2, 0.5, 1.0, 2.0, 5.0, 10.0)
  best <- NULL
  for (b0 in init_betas) {
    b0  <- pmin(pmax(b0, SEIR_BOUNDS$beta[1]), SEIR_BOUNDS$beta[2])
    z0  <- seir_to_unc(b0)
    res <- tryCatch(
      optim(z0, obj_fn, method = "Nelder-Mead",
            control = list(maxit = maxit, reltol = 1e-7)),
      error = function(e) NULL
    )
    if (!is.null(res) && (is.null(best) || res$value < best$value))
      best <- res
  }
  if (is.null(best)) {
    cat("  NM calibration failed (all starts returned NULL).\n")
    return(NULL)
  }

  theta_best <- seir_from_unc(best$par)
  p_best     <- seir_resolve(theta_best, known)
  beta_best  <- as.numeric(theta_best["beta"])
  R0_best    <- beta_best / known$recov

  # Ensemble CI with best parameters
  ci <- tryCatch(
    seir_run_multi_ci(p_best, ndays = n_days,
                      nreps = n_reps, nthreads = nthreads),
    error = function(e) {
      cat("  NM ensemble failed:", conditionMessage(e), "\n")
      NULL
    }
  )

  list(beta      = beta_best,
       R0        = R0_best,
       p         = p_best,
       ci        = ci,
       smape     = best$value,
       converged = (best$convergence == 0))
}

# =============================================================================
# Plot: observed + NM ribbon + BiLSTM ribbon
# =============================================================================

plot_reconstruction <- function(obs_inc, nm_ci, lstm_ci,
                                dates = NULL, title = "", out_png = NULL) {
  n  <- length(obs_inc)
  xv <- if (!is.null(dates)) dates else seq_len(n)

  df_obs <- data.frame(x = xv, y = obs_inc)

  g <- ggplot()

  if (!is.null(nm_ci)) {
    df_nm <- data.frame(x = xv, lower = nm_ci$lower,
                        med = nm_ci$med, upper = nm_ci$upper)
    g <- g +
      geom_ribbon(data = df_nm,
                  aes(x = x, ymin = lower, ymax = upper),
                  fill = "#6A1B9A", alpha = 0.18) +
      geom_line(data = df_nm,
                aes(x = x, y = med, colour = "NelderMead"), linewidth = 1)
  }

  if (!is.null(lstm_ci)) {
    df_lstm <- data.frame(x = xv, lower = lstm_ci$lower,
                          med = lstm_ci$med, upper = lstm_ci$upper)
    g <- g +
      geom_ribbon(data = df_lstm,
                  aes(x = x, ymin = lower, ymax = upper),
                  fill = "#1565C0", alpha = 0.18) +
      geom_line(data = df_lstm,
                aes(x = x, y = med, colour = "BiLSTM"), linewidth = 1)
  }

  g <- g +
    geom_point(data = df_obs, aes(x = x, y = y),
               colour = "black", size = 0.8, alpha = 0.6) +
    geom_line(data = df_obs, aes(x = x, y = y),
              colour = "black", linewidth = 0.35, alpha = 0.5) +
    scale_colour_manual(
      values = c(NelderMead = "#6A1B9A", BiLSTM = "#1565C0"),
      name   = "Method"
    ) +
    labs(title    = title,
         subtitle = "Black = observed; ribbon = 95% SEIR ensemble CI",
         x        = if (!is.null(dates)) "Date" else "Day",
         y        = "Daily cases") +
    theme_bw(base_size = 12) +
    theme(legend.position = "top")

  if (!is.null(out_png)) {
    ggsave(out_png, g, width = 9, height = 5, dpi = 150)
    cat(sprintf("  Saved: %s\n", basename(out_png)))
  }
  invisible(g)
}

# =============================================================================
# Plot: BiLSTM R0 predictions across windows
# =============================================================================

plot_bilstm_windows <- function(preds, nm_R0 = NA_real_,
                                title = "", out_png = NULL) {
  preds <- preds %>%
    mutate(window = factor(window, levels = window[order(t_start, win_len)]))

  g <- ggplot(preds, aes(x = window, y = R0_pred, fill = regime)) +
    geom_col(colour = "white", width = 0.7) +
    scale_fill_brewer(palette = "Set2", name = "Regime") +
    labs(title = title,
         x     = "Window",
         y     = "Predicted R0 (Bernardo BiLSTM)") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

  if (!is.na(nm_R0)) {
    g <- g +
      geom_hline(yintercept = nm_R0, linetype = "dashed",
                 colour = "#6A1B9A", linewidth = 0.9) +
      annotate("text", x = 1, y = nm_R0 * 1.05,
               label = sprintf("NM R0=%.2f", nm_R0),
               colour = "#6A1B9A", hjust = 0, size = 3.5)
  }

  if (!is.null(out_png)) {
    ggsave(out_png, g, width = 10, height = 5, dpi = 150)
    cat(sprintf("  Saved: %s\n", basename(out_png)))
  }
  invisible(g)
}

# =============================================================================
# Process one dataset end-to-end
# =============================================================================

process_dataset <- function(label, inc_csv, meta_csv, lstm_csv,
                             prefer_window = "late_365d",
                             nm_nreps = 200L, nm_threads = 4L) {

  cat(sprintf("\n=== %s ===\n", label))

  missing_files <- c(inc_csv, meta_csv)[!file.exists(c(inc_csv, meta_csv))]
  if (length(missing_files)) {
    cat("  Missing files — skipping:\n")
    for (f in missing_files) cat("    ", f, "\n")
    return(NULL)
  }

  inc_df  <- read.csv(inc_csv)
  meta    <- read.csv(meta_csv)[1L, ]
  obs     <- as.numeric(inc_df$daily_cases)
  n_days  <- length(obs)
  dates   <- if ("date" %in% names(inc_df)) as.Date(inc_df$date) else NULL

  n_pop      <- as.numeric(meta$n_pop)
  incub_days <- as.numeric(meta$incub_days)
  recov_rate <- as.numeric(meta$recov_rate)
  init <- seir_initial_conditions(obs, n_pop, incub_days, recov_rate)

  known <- list(n          = n_pop,
                prevalence = init$prevalence,
                initial_infected_fraction = init$initial_infected_fraction,
                incub      = incub_days,
                recov      = recov_rate)

  cat(sprintf(paste0(
      "  %d days | n=%g, initial E=%.1f, initial I=%.1f, ",
      "incub=%.0f d, recov=%.3f/d\n"),
      n_days, n_pop, init$initial_exposed, init$initial_infected,
      incub_days, recov_rate))

  # ── Nelder-Mead ─────────────────────────────────────────────────────────────
  cat("  Running NM calibration (multi-start)...\n")
  t0_nm  <- proc.time()[["elapsed"]]
  nm_fit <- nm_calibrate_real(obs, known, n_reps = nm_nreps,
                               nthreads = nm_threads)
  nm_sec <- proc.time()[["elapsed"]] - t0_nm

  if (!is.null(nm_fit)) {
    cat(sprintf("  NM : beta=%.4f  R0=%.3f  SMAPE=%.1f%%  (%.1f s)\n",
        nm_fit$beta, nm_fit$R0, nm_fit$smape, nm_sec))
  }

  # ── BiLSTM (Bernardo) — pick best window ─────────────────────────────────────
  lstm_beta <- NA_real_; lstm_R0 <- NA_real_
  lstm_ci   <- NULL
  preds_df  <- NULL

  if (file.exists(lstm_csv)) {
    preds_df <- read.csv(lstm_csv)

    # Prefer an explicitly requested prediction. If it is unavailable, fall
    # back to the longest window that starts with the reconstructed series;
    # this avoids silently substituting a late, shorter slice merely because
    # it happens to be listed first among equally long windows.
    anchored <- preds_df[preds_df$t_start == 0, ]

    if (prefer_window %in% preds_df$window) {
      row_lstm <- preds_df[preds_df$window == prefer_window, ][1L, ]
    } else if (nrow(anchored) > 0) {
      row_lstm <- anchored[which.max(anchored$win_len), ][1L, ]
    } else {
      row_lstm <- preds_df[which.max(preds_df$win_len), ][1L, ]
    }
    lstm_beta <- row_lstm$beta_pred
    lstm_R0   <- row_lstm$R0_pred
    cat(sprintf("  BiLSTM (%s): beta=%.4f  R0=%.3f\n",
        row_lstm$window, lstm_beta, lstm_R0))

    p_lstm <- seir_resolve(c(beta = lstm_beta), known)
    lstm_ci <- tryCatch(
      seir_run_multi_ci(p_lstm, ndays = n_days,
                        nreps = nm_nreps, nthreads = nm_threads),
      error = function(e) {
        cat("  BiLSTM SEIR ensemble failed:", conditionMessage(e), "\n")
        NULL
      }
    )
  } else {
    cat(sprintf("  BiLSTM predictions not found (%s)\n", basename(lstm_csv)))
    cat("  -> Run: python real_data/bernardo_predict_real.py\n")
  }

  # ── Reconstruction plot ──────────────────────────────────────────────────────
  slug <- gsub("[^A-Za-z0-9]", "_", tolower(label))
  plot_reconstruction(
    obs_inc = obs,
    nm_ci   = if (!is.null(nm_fit)) nm_fit$ci else NULL,
    lstm_ci = lstm_ci,
    dates   = dates,
    title   = sprintf("%s — SEIR Reconstruction", label),
    out_png = file.path(PLOTS_DIR, sprintf("%s_reconstruction.png", slug))
  )

  # ── BiLSTM window R0 plot ────────────────────────────────────────────────────
  if (!is.null(preds_df)) {
    plot_bilstm_windows(
      preds   = preds_df,
      nm_R0   = if (!is.null(nm_fit)) nm_fit$R0 else NA_real_,
      title   = sprintf("%s — Bernardo BiLSTM R0 by window", label),
      out_png = file.path(PLOTS_DIR, sprintf("%s_bilstm_windows.png", slug))
    )
  }

  data.frame(
    dataset    = label,
    n_pop      = n_pop,
    incub_days = incub_days,
    recov_rate = recov_rate,
    n_days     = n_days,
    nm_beta    = if (!is.null(nm_fit)) nm_fit$beta  else NA_real_,
    nm_R0      = if (!is.null(nm_fit)) nm_fit$R0    else NA_real_,
    nm_smape   = if (!is.null(nm_fit)) nm_fit$smape else NA_real_,
    nm_time_s  = nm_sec,
    lstm_beta  = lstm_beta,
    lstm_R0    = lstm_R0,
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# Run both datasets
# =============================================================================

results <- list(

  process_dataset(
    label          = "Utah COVID-19 (Wave 1)",
    inc_csv        = file.path(REAL_DIR, "utah_covid_wave1.csv"),
    meta_csv       = file.path(REAL_DIR, "utah_covid_meta.csv"),
    lstm_csv       = file.path(REAL_DIR, "bernardo_real_covid_predictions.csv"),
    prefer_window  = "late_365d",  # full first-wave window
    nm_nreps       = 200L,
    nm_threads     = 4L
  ),

  process_dataset(
    label          = "Hagelloch Measles 1861",
    inc_csv        = file.path(REAL_DIR, "measles_hagelloch_incidence.csv"),
    meta_csv       = file.path(REAL_DIR, "measles_hagelloch_meta.csv"),
    lstm_csv       = file.path(REAL_DIR, "bernardo_real_measles_predictions.csv"),
    prefer_window  = "full_087d",  # full 87-day outbreak, anchored at day 1
    nm_nreps       = 300L,         # more reps for small stochastic population
    nm_threads     = 4L
  )
)

# =============================================================================
# Save summary
# =============================================================================

res_df <- do.call(rbind, Filter(Negate(is.null), results))

if (nrow(res_df) > 0) {
  write.csv(res_df, file.path(REAL_DIR, "real_data_results.csv"),
            row.names = FALSE)
  cat("\n=== Final parameter estimates ===\n")
  print(res_df[, c("dataset","nm_beta","nm_R0","nm_smape","lstm_beta","lstm_R0")],
        row.names = FALSE)
  cat(sprintf("\nResults : real_data/real_data_results.csv\n"))
  cat(sprintf("Plots   : real_data/plots/\n"))
} else {
  cat("\nNo results to save — check that input CSVs exist.\n")
}
