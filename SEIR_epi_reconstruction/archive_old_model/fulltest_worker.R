# =============================================================================
#  fulltest_worker.R
#  Worker function: evaluate one sim_idx across ALL available windows.
#  Returns:
#    $metrics   — one row per (sim x window): parameter + curve summary errors
#    $mae_day   — one row per (sim x window x day): per-day |act - pred|
#
#  Sourced by fulltest_submit.R before Slurm_lapply().
# =============================================================================

eval_one_sim <- function(sim_idx) {

  suppressPackageStartupMessages({
    library(epiworldR)
    library(dplyr)
    library(reticulate)
  })

  PROJECT_DIR <- path.expand("~/DeepIMC")
  DATA_DIR    <- PROJECT_DIR

  NDAYS        <- 365L
  CONTACT_RATE <- 1

  # ---------------------------------------------------------------------------
  # Load test data
  # ---------------------------------------------------------------------------
  actual    <- read.csv(path.expand(file.path(DATA_DIR, "test_actual_parameters (2).csv")))
  preds_all <- read.csv(path.expand(file.path(DATA_DIR, "test_bilstm_predictions (2).csv")))

  np      <- reticulate::import("numpy")
  inc_raw <- as.matrix(
    np$load(path.expand(file.path(DATA_DIR, "test_incidence_raw (2).npy")))
  )

  # ---------------------------------------------------------------------------
  # Helper: single ModelSEIRCONN run → daily incidence vector
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
  # Look up this sim
  # ---------------------------------------------------------------------------
  sid     <- as.integer(sim_idx)
  act     <- actual[actual$sim_idx == sid, ]
  mat_row <- which(actual$sim_idx == sid)

  if (nrow(act) == 0 || length(mat_row) == 0) return(NULL)

  n          <- act$n[1]
  recov      <- act$recov[1]
  incub      <- act$incub[1]
  prevalence <- act$prevalence[1]

  # Actual-parameter SEIR curve (computed once, reused for all windows)
  act_inc <- tryCatch(
    run_seir_single(n, prevalence, act$beta[1], incub, recov, seed = sid),
    error = function(e) NULL
  )
  if (is.null(act_inc)) return(NULL)

  # Observed ABM incidence
  obs_inc <- as.numeric(inc_raw[mat_row, ])

  # ---------------------------------------------------------------------------
  # Loop over all windows that have a prediction for this sim
  # ---------------------------------------------------------------------------
  sim_preds <- preds_all[preds_all$sim_idx == sid, ]
  windows   <- unique(sim_preds$window)

  metric_rows  <- vector("list", length(windows))
  mae_day_rows <- vector("list", length(windows))

  for (i in seq_along(windows)) {

    win_tag <- windows[i]
    prd     <- sim_preds[sim_preds$window == win_tag, ]
    if (nrow(prd) == 0) next

    pred_inc <- tryCatch(
      run_seir_single(n, prevalence, prd$beta_pred[1], incub, recov, seed = sid),
      error = function(e) NULL
    )
    if (is.null(pred_inc)) next

    regime  <- sub("_.*", "", win_tag)
    win_len <- as.integer(regmatches(win_tag, regexpr("\\d+", win_tag)))

    mae_vec <- abs(act_inc - pred_inc)

    # Summary metrics (one row per sim x window)
    metric_rows[[i]] <- data.frame(
      sim_idx        = sid,
      window         = win_tag,
      regime         = regime,
      win_len        = win_len,
      true_beta      = act$beta[1],
      pred_beta      = prd$beta_pred[1],
      true_R0        = act$R0[1],
      pred_R0        = prd$R0_pred[1],
      err_beta       = abs(prd$beta_pred[1] - act$beta[1]),
      err_R0         = abs(prd$R0_pred[1]   - act$R0[1]),
      curve_mae      = mean(mae_vec),
      obs_curve_mae  = mean(abs(obs_inc - pred_inc)),
      peak_day_act   = which.max(act_inc),
      peak_day_pred  = which.max(pred_inc),
      peak_day_err   = abs(which.max(pred_inc) - which.max(act_inc)),
      peak_size_act  = max(act_inc),
      peak_size_pred = max(pred_inc),
      peak_size_err  = abs(max(pred_inc) - max(act_inc)),
      stringsAsFactors = FALSE
    )

    # Per-day incidence + MAE (one row per sim x window x day)
    mae_day_rows[[i]] <- data.frame(
      sim_idx  = sid,
      window   = win_tag,
      regime   = regime,
      win_len  = win_len,
      day      = 1:NDAYS,
      act_inc  = act_inc,
      pred_inc = pred_inc,
      obs_inc  = obs_inc,
      mae      = mae_vec,
      stringsAsFactors = FALSE
    )
  }

  list(
    metrics = dplyr::bind_rows(metric_rows),
    mae_day = dplyr::bind_rows(mae_day_rows)
  )
}
