# =============================================================================
#  seir_common.R — shared helpers for SEIR calibration
#
#  Sourced by every script and every SLURM worker.
#
#  CALIBRATED PARAMETER: beta (= contact_rate * transmission_rate).
#    The simulator (ModelSEIRCONN) needs separate crate and ptran, so beta is
#    decomposed canonically:
#        crate = max(beta, 1.0)
#        ptran = beta / crate           => ptran is always in (0, 1]
#    This gives crate * ptran == beta exactly, and ptran never exceeds 1.
#
#  KNOWN (not calibrated): n, prevalence, incub, recov.
#
#  DERIVED QUANTITIES
#    R0 = beta / recov
#
#  TEST DATA (4000 sets, held-out split):
#    test_actual_parameters (2).csv   — columns: sim_idx, n, recov, incub,
#                                                 prevalence, beta, R0
#    test_incidence_raw (2).npy       — shape (n_test, 365), E->I incidence
#
#  INCIDENCE: Exposed -> Infected (E->I) transitions, days 1..ndays.
# =============================================================================

SEIR_LIBPATH <- "/uufs/chpc.utah.edu/common/home/u1418987/R/x86_64-pc-linux-gnu-library/4.4"

seir_set_libpath <- function() {
  .libPaths(c(SEIR_LIBPATH, .libPaths()))
  invisible(NULL)
}

# =============================================================================
# Canonical data file paths
# =============================================================================

SEIR_TEST_PARAMS_FILE   <- "test_actual_parameters (2).csv"
SEIR_TEST_INCIDENCE_NPY <- "test_incidence_raw (2).npy"

# =============================================================================
# Calibrated parameter
# =============================================================================

SEIR_PARS <- c("beta")   # single effective transmission rate

# Prior / search bounds
SEIR_BOUNDS <- list(
  beta  = c(0.001, 20.0),   # = crate * ptran
  recov = c(0.02,   0.60),
  R0    = c(0.20,   15.0)
)

SEIR_NDAYS <- 365L

seir_bounds_mat <- function(pars = SEIR_PARS) {
  lo <- vapply(pars, function(p) SEIR_BOUNDS[[p]][1], numeric(1))
  hi <- vapply(pars, function(p) SEIR_BOUNDS[[p]][2], numeric(1))
  list(lo = unname(lo), hi = unname(hi))
}

seir_in_box <- function(theta, pars = SEIR_PARS) {
  b <- seir_bounds_mat(pars)
  all(theta >= b$lo & theta <= b$hi)
}

# Unconstrained reparameterisation (logit on the box) for Nelder-Mead
seir_to_unc <- function(theta, pars = SEIR_PARS) {
  b <- seir_bounds_mat(pars)
  u <- (as.numeric(theta) - b$lo) / (b$hi - b$lo)
  u <- pmin(pmax(u, 1e-6), 1 - 1e-6)
  log(u / (1 - u))
}

seir_from_unc <- function(z, pars = SEIR_PARS) {
  b <- seir_bounds_mat(pars)
  u <- 1 / (1 + exp(-as.numeric(z)))
  setNames(as.numeric(b$lo + u * (b$hi - b$lo)), pars)
}

# =============================================================================
# Canonical beta -> (crate, ptran) decomposition
# =============================================================================

seir_beta_to_cp <- function(beta) {
  beta  <- as.numeric(beta)
  crate <- max(beta, 1.0)
  ptran <- beta / crate   # always in (0, 1]
  c(crate = crate, ptran = ptran)
}

# =============================================================================
# Assemble a full parameter list from calibrated + known values
# =============================================================================

seir_resolve <- function(theta, known, pars = SEIR_PARS) {
  p  <- as.list(known)
  th <- as.numeric(theta)
  for (i in seq_along(pars)) p[[pars[i]]] <- th[i]

  # Parameterisation: calibrating beta directly -> decompose to crate / ptran
  if ("beta" %in% pars) {
    cp      <- seir_beta_to_cp(p$beta)
    p$crate <- cp["crate"]
    p$ptran <- cp["ptran"]
  }

  # Old parameterisation B (kept for backwards compat): calibrating R0 directly
  if ("R0" %in% pars && "crate" %in% names(p)) {
    p$ptran <- p$R0 * p$recov / max(p$crate, 1e-12)
  }

  p
}

seir_feasible <- function(theta, known, pars = SEIR_PARS) {
  if (!seir_in_box(theta, pars)) return(FALSE)
  p <- seir_resolve(theta, known)
  is.finite(p$crate) && p$crate > 0 &&
    is.finite(p$ptran) && p$ptran > 0 && p$ptran <= 1 &&
    is.finite(p$recov) && p$recov > 0
}

# Derive beta and R0 from any parameter list that has either
#   (crate, ptran, recov)  OR  (beta, recov)
seir_derive <- function(p) {
  beta <- if (!is.null(p$crate) && !is.null(p$ptran)) {
    p$crate * p$ptran
  } else {
    as.numeric(p$beta)
  }
  list(beta = beta, R0 = beta / p$recov)
}

# =============================================================================
# Model construction and simulation
# =============================================================================

seir_build <- function(p) {
  n_int <- max(as.integer(round(p$n)), 10L)

  # Floating-point trap (found 2026-08-17 on n=187, prevalence=1/187 exactly):
  # 1/187 * 187 == 0.99999999999999988898 in double precision, not 1.0. If
  # ModelSEIRCONN converts prevalence*n to an integer seed count by
  # truncation rather than rounding, that silently seeds ZERO initial
  # infections instead of 1 -- a completely inert simulation (no epidemic
  # possible; get_hist_transition_matrix() then only ever shows the
  # "Susceptible" state, since nothing ever transitions, and
  # seir_incidence_e2i()/seir_run_multi_ci() error out with "No Exposed ->
  # Infected transitions"). seed_n/n_int reproduces the exact same trap
  # (it's the same value), so nudge prev just above the seed_n/n_int
  # boundary by a relative epsilon -- comfortably larger than double
  # precision error, negligible next to the seed count itself -- so
  # prev*n_int safely truncates (or rounds) back to seed_n instead of
  # seed_n - 1.
  seed_n <- max(1L, as.integer(round(min(p$prevalence, 1.0) * n_int)))
  prev   <- min((seed_n / n_int) * (1 + 1e-6), 1.0)

  m <- epiworldR::ModelSEIRCONN(
    name              = "sim",
    n                 = n_int,
    prevalence        = prev,
    contact_rate      = max(p$crate, 1e-6),
    transmission_rate = max(min(p$ptran, 1.0), 1e-8),
    incubation_days   = max(p$incub, 1.0),
    recovery_rate     = max(min(p$recov, 1.0), 1e-6)
  )
  epiworldR::verbose_off(m)
  m
}

seir_incidence_e2i <- function(model, ndays) {
  tm <- as.data.frame(epiworldR::get_hist_transition_matrix(model))
  cc <- seir_transition_cols(tm)

  from <- as.character(tm[[cc$from]])
  to   <- as.character(tm[[cc$to]])
  date <- as.integer(tm[[cc$date]])
  cnt  <- as.numeric(tm[[cc$counts]])

  is_e <- grepl("^expos",  from, ignore.case = TRUE)
  is_i <- grepl("^infect", to,   ignore.case = TRUE)
  if (!any(is_e) || !any(is_i)) {
    stop("No Exposed -> Infected transitions. States: ",
         paste(unique(c(from, to)), collapse = ", "))
  }

  keep <- is_e & is_i & date >= 1L & date <= ndays
  out  <- numeric(ndays)
  if (any(keep)) {
    agg <- tapply(cnt[keep], date[keep], sum)
    out[as.integer(names(agg))] <- as.numeric(agg)
  }
  out
}

seir_run_single <- function(p, ndays = SEIR_NDAYS, seed = NULL) {
  m <- seir_build(p)
  if (is.null(seed)) {
    epiworldR::run(m, ndays = ndays + 1L)
  } else {
    epiworldR::run(m, ndays = ndays + 1L, seed = as.integer(seed))
  }
  seir_incidence_e2i(m, ndays)
}

seir_run_multi_ci <- function(p, ndays = SEIR_NDAYS, nreps = 300L,
                              nthreads = 1L, seed = NULL) {

  m     <- seir_build(p)
  saver <- epiworldR::make_saver("transition")

  if (is.null(seed)) {
    epiworldR::run_multiple(m, ndays = ndays + 1L, nsims = nreps,
                            saver = saver, nthreads = nthreads)
  } else {
    epiworldR::run_multiple(m, ndays = ndays + 1L, nsims = nreps,
                            saver = saver, nthreads = nthreads,
                            seed = as.integer(seed))
  }

  # NOTE: nthreads=1 here (not the simulation's nthreads) -- run_multiple_get_results()
  # reads per-rep CSVs via parallel::makeCluster()/parLapply() when nthreads > 1, and
  # those cluster workers are fresh R processes that don't inherit this session's
  # custom .libPaths()/LD_LIBRARY_PATH, so data.table::fread() fails inside them with
  # "could not find function %chin%" (bug found 2026-08-17). Serial reading is fine
  # here since it only runs once per ensemble, not per simulation.
  res <- epiworldR::run_multiple_get_results(m, nthreads = 1L,
                                             freader = data.table::fread)
  tr  <- data.table::as.data.table(res$transition)

  tc <- seir_transition_cols(tr)
  data.table::setnames(tr, c(tc$from, tc$to, tc$date, tc$counts),
                       c("from", "to", "date", "counts"))

  tr[, keep_e2i := grepl("^expos", from, ignore.case = TRUE) &
                   grepl("^infect", to,  ignore.case = TRUE)]
  if (!any(tr$keep_e2i)) {
    stop("No E->I transitions in run_multiple output. States: ",
         paste(unique(c(tr$from, tr$to)), collapse = ", "))
  }

  id_col <- intersect(c("sim_num", "nsim", "sim", "run", "id"), names(tr))
  id_col <- if (length(id_col) > 0) id_col[1] else NA_character_

  if (is.na(id_col)) {
    warning("No per-run id column; quantiles may be biased by missing zeros.")
    q <- tr[keep_e2i == TRUE & date >= 1L & date <= ndays,
            .(lower = stats::quantile(counts, 0.025, names = FALSE),
              med   = stats::quantile(counts, 0.500, names = FALSE),
              upper = stats::quantile(counts, 0.975, names = FALSE)), by = date]
  } else {
    all_ids <- unique(tr[[id_col]])
    data.table::setnames(tr, id_col, "rep_id")
    tr <- tr[keep_e2i == TRUE & date >= 1L & date <= ndays,
             .(rep_id, date, counts)]
    grid <- data.table::CJ(rep_id = all_ids, date = seq_len(ndays))
    full <- merge(grid, tr, by = c("rep_id", "date"), all.x = TRUE)
    full[is.na(counts), counts := 0]
    q <- full[, .(lower = stats::quantile(counts, 0.025, names = FALSE),
                  med   = stats::quantile(counts, 0.500, names = FALSE),
                  upper = stats::quantile(counts, 0.975, names = FALSE)),
              by = date]
  }

  q <- merge(data.table::data.table(date = seq_len(ndays)), q,
             by = "date", all.x = TRUE)
  for (cl in c("lower", "med", "upper")) {
    data.table::set(q, which(is.na(q[[cl]])), cl, 0)
  }
  as.data.frame(q[order(date)])
}

# =============================================================================
# Transition-table column mapping
# =============================================================================

seir_transition_cols <- function(df) {
  nm <- tolower(names(df))
  pick <- function(cands) {
    i <- which(nm %in% cands)[1]
    if (is.na(i)) NA_character_ else names(df)[i]
  }
  cols <- list(
    from   = pick(c("from", "state_from", "source", "origin")),
    to     = pick(c("to",   "state_to",   "target", "dest", "destination")),
    date   = pick(c("date", "day", "time", "t")),
    counts = pick(c("counts", "count", "n", "value", "freq"))
  )
  if (any(vapply(cols, is.na, logical(1)))) {
    stop("Cannot map transition columns. Available: ",
         paste(names(df), collapse = ", "))
  }
  cols
}

# =============================================================================
# Metrics
# =============================================================================

seir_smape <- function(obs, pred) {
  mean(2 * abs(obs - pred) / (abs(obs) + abs(pred) + 1e-6)) * 100
}

seir_cor <- function(obs, pred) {
  if (stats::sd(obs) == 0 || stats::sd(pred) == 0) return(NA_real_)
  suppressWarnings(stats::cor(obs, pred))
}

# One row per (sim, method).
# true_p : list with at minimum (beta, recov).  May also have (crate, ptran).
# pred_p : same structure for predicted values.
# known  : list with (n, prevalence, incub).  recov if not in true_p.
seir_metrics <- function(sim_idx, method, obs_cmp, pred_cmp,
                         true_p, pred_p, known,
                         converged = TRUE, time_sec = NA_real_,
                         cpu_sec = NA_real_, n_evals = NA_integer_,
                         n_model_runs = NA_integer_,
                         pars = SEIR_PARS) {

  td <- seir_derive(true_p)
  pd <- if (converged) seir_derive(pred_p) else list(beta = NA_real_, R0 = NA_real_)

  out <- data.frame(
    sim_idx    = sim_idx,
    method     = method,
    converged  = converged,
    calibrated = paste(pars, collapse = "+"),
    n          = known$n,
    prevalence = known$prevalence,
    incub      = known$incub,
    stringsAsFactors = FALSE
  )

  # Report true / pred / error for each calibrated parameter
  for (nm in pars) {
    tv <- if (!is.null(true_p[[nm]])) true_p[[nm]] else NA_real_
    pv <- if (converged && !is.null(pred_p[[nm]])) pred_p[[nm]] else NA_real_
    out[[paste0("true_", nm)]]        <- tv
    out[[paste0("pred_", nm)]]        <- pv
    out[[paste0("err_", nm, "_pct")]] <- abs(pv - tv) / (abs(tv) + 1e-9) * 100
  }

  # Always report derived beta and R0 regardless of parameterisation
  out$true_beta    <- td$beta
  out$pred_beta    <- pd$beta
  out$err_beta_pct <- abs(pd$beta - td$beta) / (td$beta + 1e-9) * 100
  out$true_R0      <- td$R0
  out$pred_R0      <- pd$R0
  out$err_R0_pct   <- abs(pd$R0 - td$R0) / (td$R0 + 1e-9) * 100

  if (converged && !any(is.na(pred_cmp))) {
    out$mae            <- mean(abs(obs_cmp - pred_cmp))
    out$rmse           <- sqrt(mean((obs_cmp - pred_cmp)^2))
    out$smape          <- seir_smape(obs_cmp, pred_cmp)
    out$pearson_r      <- seir_cor(obs_cmp, pred_cmp)
    out$peak_day_true  <- which.max(obs_cmp)  + 1L
    out$peak_day_pred  <- which.max(pred_cmp) + 1L
    out$peak_day_err   <- abs(which.max(pred_cmp) - which.max(obs_cmp))
    out$peak_size_true <- max(obs_cmp)
    out$peak_size_pred <- max(pred_cmp)
    out$rel_peak_err   <- abs(max(pred_cmp) - max(obs_cmp)) / (max(obs_cmp) + 1e-6) * 100
    out$dieout_pred    <- as.integer(max(pred_cmp) <= 1)
  } else {
    for (cl in c("mae","rmse","smape","pearson_r","peak_day_true","peak_day_pred",
                 "peak_day_err","peak_size_true","peak_size_pred","rel_peak_err",
                 "dieout_pred")) out[[cl]] <- NA
  }

  out$time_sec     <- time_sec
  out$cpu_sec      <- cpu_sec
  out$n_evals      <- n_evals
  out$n_model_runs <- n_model_runs
  out
}

# =============================================================================
# Plot constants
# =============================================================================

SEIR_METHOD_ORDER  <- c("ABC", "ABC-SMC", "NelderMead", "DE", "BiLSTM")
SEIR_METHOD_COLORS <- c(
  Oracle     = "#388E3C",
  ABC        = "#E65100",
  `ABC-SMC`  = "#D4537E",
  NelderMead = "#6A1B9A",
  DE         = "#00695C",
  BiLSTM     = "#1565C0"
)

seir_mkey <- function(m) gsub("[^A-Za-z0-9]", "", m)
