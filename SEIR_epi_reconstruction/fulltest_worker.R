# =============================================================================
#  fulltest_worker.R  — hypertuned BiLSTM
#  Worker function: evaluate one sim_idx across ALL available windows.
#  Returns:
#    $metrics   — one row per (sim x window): parameter + curve summary errors
#    $mae_day   — one row per (sim x window x day): per-day |act - pred|
#
#  Sourced by fulltest_submit.R before Slurm_lapply().
#  NOTE: day 1 is included in raw data but excluded from plots in submit script.
#
#  INCIDENCE DEFINITION: Exposed -> Infected (E->I) transitions, consistent
#  with seir_common.R and the training data. plot_incidence() is NOT used.
# =============================================================================

eval_one_sim <- function(sim_idx) {

  suppressPackageStartupMessages({
    library(epiworldR)
    library(dplyr)
    library(reticulate)
  })

  PROJECT_DIR <- path.expand(
    "~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction"
  )

  source(file.path(PROJECT_DIR, "seir_common.R"))
  seir_set_libpath()

  NDAYS <- SEIR_NDAYS

  # ---------------------------------------------------------------------------
  # Helper: single SEIR run -> E->I daily incidence (days 1-365)
  # Uses seir_build and seir_incidence_e2i from seir_common.R for consistency
  # with training data and all other calibration scripts.
  # beta is decomposed as: crate = max(beta, 1.0), ptran = beta / crate
  # This ensures ptran <= 1 while preserving the effective beta.
  # ---------------------------------------------------------------------------
  run_seir_e2i <- function(n, prevalence, beta, incub, recov, seed = 1L) {
    crate <- max(beta, 1.0)
    ptran <- beta / crate
    p <- list(
      n          = n,
      prevalence = prevalence,
      crate      = crate,
      ptran      = ptran,
      incub      = incub,
      recov      = recov
    )
    seir_run_single(p, ndays = NDAYS, seed = as.integer(seed))
  }

  # ---------------------------------------------------------------------------
  # Load test data (tuned model predictions)
  # ---------------------------------------------------------------------------
  actual    <- read.csv(path.expand(file.path(PROJECT_DIR, SEIR_TEST_PARAMS_FILE)))
  preds_all <- read.csv(path.expand(file.path(PROJECT_DIR, "test_bilstm_predictions_tuned.csv")))

  np      <- reticulate::import("numpy")
  inc_raw <- as.matrix(
    np$load(path.expand(file.path(PROJECT_DIR, SEIR_TEST_INCIDENCE_NPY)))
  )

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

  act_inc <- tryCatch(
    run_seir_e2i(n, prevalence, act$beta[1], incub, recov, seed = sid),
    error = function(e) NULL
  )
  if (is.null(act_inc)) return(NULL)

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
      run_seir_e2i(n, prevalence, prd$beta_pred[1], incub, recov, seed = sid),
      error = function(e) NULL
    )
    if (is.null(pred_inc)) next

    regime  <- sub("_.*", "", win_tag)
    win_len <- as.integer(regmatches(win_tag, regexpr("\\d+", win_tag)))

    mae_vec <- abs(act_inc - pred_inc)

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
