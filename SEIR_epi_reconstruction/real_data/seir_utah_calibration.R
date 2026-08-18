################################################################################
#
#  seir_utah_calibration.R
#
#  Calibrates a stochastic SEIR agent-based model (epiworldR::ModelSEIRCONN)
#  to Utah COVID-19 daily case data, by two methods:
#
#     Method A  Fast direct search on R0        (grid scan + golden section)
#     Method B  ABC via likelihood-free MCMC    (epiworldR::LFMCMC)
#
#  plus a BiLSTM prediction on the SAME window (bernardo_predict_current61.py).
#
#  Everything is in this one file. The block marked "USER SETTINGS" is the only
#  part you should normally need to edit.
#
#  Corrected 2026-08-17 from a draft with 4 issues found before first run:
#    1. SEIR_PARAM labels didn't match installed epiworldR (verified directly
#       against summary(ModelSEIRCONN(...)): "Prob. Transmission", not
#       "Transmission rate"; "Prob. Recovery", not "Recovery rate";
#       "Avg. Incubation days", not "Incubation days"). Would have failed
#       loudly at the first seir_set_param() call in Method A -- the script's
#       own fail-loud design caught this before it could produce silent
#       garbage, but it still needed a human/LLM to fix the label.
#    2. N_MODEL=5000 was a 4th independently hard-coded population value,
#       duplicating (and sitting right at the boundary of) what's already in
#       seir_scale_config.json's n_bilstm=8000. Now reads from that file.
#    3. INCUBATION_DAYS=3 didn't match the rest of the project (Utah SEIR
#       meta uses 5 days, consistent with published COVID incubation
#       estimates ~5 days). Changed to 5 for consistency, not because 3 was
#       provably wrong -- there was just no reason to introduce a 3rd value.
#    4. BILSTM_CSV pointed at bernardo_real_covid_predictions.csv, which was
#       built from wave1.csv (2020-03-18 to 2021-03-17). This script fits
#       utah_covid_data's last 61 days, which as packaged today is
#       2025-03-07 to 2025-05-06 -- a totally different era. Comparing the
#       two would silently mix periods. Now points at a fresh prediction on
#       the actual fitting window (bernardo_predict_current61.py).
#
#  Usage:  Rscript seir_utah_calibration.R
#          (or source() it in an interactive session)
#
################################################################################


## =============================================================================
## 0. LIBRARIES
## =============================================================================

suppressPackageStartupMessages({
  library(epiworldR)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(data.table)
  library(jsonlite)
})

set.seed(122)
MODEL_SEED <- 122L


## =============================================================================
## 1. USER SETTINGS -- every assumption lives here, and only here
## =============================================================================
##
## READ THIS FIRST.
##
## utah_covid_data has exactly two columns: Date and Daily.Cases. It contains
## NO population size, NO recovery rate, NO incubation period and NO reporting
## fraction. An SEIR model needs all four. They are modelling CHOICES, not
## measurements, and every number this script prints is conditional on them.
##
## Previously these were scattered across three scripts: N = 20000 in the ABC
## script, N = 5000 in the vignette, recovery 1/7 in both but 1/10 in
## utah_covid_meta.csv. That last inconsistency is why nm_R0 and lstm_R0 in
## real_data_results.csv were never comparable -- the gap between them was
## arithmetic, not methodology.
## -----------------------------------------------------------------------------

# ---- Population --------------------------------------------------------------

# True population at risk (Utah statewide).
UTAH_POPULATION <- 3.34e6

# Simulated population, read from the single shared config also used by
# calibrate_real_5method.R / bernardo_predict_real.py (n_bilstm) -- NOT
# hard-coded here. A pretrained network's predictions are only valid at the
# population it was trained on; see seir_scale_config.json for the trained
# range this value must stay inside.
.scale_cfg <- jsonlite::fromJSON(file.path("real_data", "seir_scale_config.json"))
N_MODEL <- as.integer(.scale_cfg$n_bilstm)

# Fraction of true infections that appear in the data as confirmed cases.
# Scales epidemic size but NOT R0.
#
# The window this script fits (see section 6) is utah_covid_data's last 61
# days, which as packaged today is 2025-03-07 to 2025-05-06 -- the
# home-testing era, where official case counts capture a much smaller share
# of infections than early-pandemic PCR-based surveillance did. 0.25 (a
# reasonable estimate for 2020) is very likely too high for this period; 0.10
# is used as a more defensible default, but this is the single most
# uncertain assumption in this file and should be revisited if the fitting
# window changes.
REPORTING_FRACTION <- 0.10


# ---- Disease durations -------------------------------------------------------

# Mean infectious period, days.  gamma = 1/INFECTIOUS_PERIOD_DAYS
#
#   *** R0 IS LINEAR IN THIS NUMBER ***
#   R0 = contact_rate * transmission_rate / gamma
#   Using 10 days instead of 7 inflates every R0 by 43% with no change to fit.
#
# 7 days matches the SIR-side calibration (vignette + ABC script) so SIR and
# SEIR R0 estimates on this project stay comparable.
INFECTIOUS_PERIOD_DAYS <- 7

# Mean latent period, days.  sigma = 1/INCUBATION_DAYS
# This is the SEIR-specific one. Note that sigma does NOT appear in R0 -- a
# latent stage does not change how many people one case infects, only how fast.
# It still moves your calibrated R0; see section 5.
#
# 5 days matches utah_covid_meta.csv (incub_days=5) and published COVID
# incubation estimates (~5 days median, Lauer et al. 2020).
INCUBATION_DAYS <- 5


# ---- Analysis window and run control -----------------------------------------

N_DAYS_WINDOW  <- 61L    # days of data to use
RUN_ABC        <- TRUE   # set FALSE to skip Method B (it is the slow one)
N_SIMS_UQ      <- 2000L  # simulations for the uncertainty ribbon
N_THREADS      <- 4L
ABC_SAMPLES    <- 3000L
ABC_BURNIN     <- 1500L
# RELATIVE tolerance; see kernel_fun in section 7. 0.25 gave only 2.5%
# acceptance on the current window (script's own check flags <10% as too
# low) -- likely because "day of peak" (one of 6 summary stats, unweighted)
# is a small integer where even a 1-2 day miss is a large relative error,
# dominating the combined distance. Loosened epsilon and tightened the
# proposal step (below) rather than reweighting the summary stats, since
# that's a smaller, more easily-verified change.
ABC_EPSILON    <- 0.45

# Transmission probability per contact. Held fixed because crate and ptran are
# only jointly identifiable through their product; the contact rate is then
# derived from crate = R0 * gamma / ptran (the same soft constraint the BiLSTM
# uses). Change it and the contact rate rescales, but R0 does not.
PTRAN_FIXED <- 0.05

# BiLSTM prediction on the SAME window this script fits (see header note #4
# above -- the old bernardo_real_covid_predictions.csv is from 2020-21 and is
# not comparable to the 2025 window fit here).
BILSTM_CSV     <- file.path("real_data", "bernardo_real_covid_current61_predictions.csv")
BILSTM_WINDOW  <- "current_061d"

OUT_DIR <- "seir_output"


# ---- Derived (do not edit) ---------------------------------------------------

RECOVERY_RATE   <- 1 / INFECTIOUS_PERIOD_DAYS   # gamma
INCUBATION_RATE <- 1 / INCUBATION_DAYS          # sigma

# Burn-in days discarded from the front of each simulation. ModelSEIRCONN seeds
# all initial agents as EXPOSED, so the first few days are an artefact of the
# initial condition rather than epidemic dynamics. See section 4.
SEED_BURN_IN <- as.integer(ceiling(INCUBATION_DAYS))

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)


print_config <- function() {
  cat("\n--- Active assumptions (none of these come from the data) ---\n")
  cat(sprintf("  True population      : %s\n", format(UTAH_POPULATION, big.mark = ",")))
  cat(sprintf("  Simulated population : %s (from seir_scale_config.json)\n", format(N_MODEL, big.mark = ",")))
  cat(sprintf("  Scale factor         : %.3e\n", N_MODEL / UTAH_POPULATION))
  cat(sprintf("  Reporting fraction   : %.2f\n", REPORTING_FRACTION))
  cat(sprintf("  Latent period        : %g days (sigma = %.4f)\n",
              INCUBATION_DAYS, INCUBATION_RATE))
  cat(sprintf("  Infectious period    : %g days (gamma = %.4f)\n",
              INFECTIOUS_PERIOD_DAYS, RECOVERY_RATE))
  cat(sprintf("  Seed burn-in         : %d days\n", SEED_BURN_IN))
  cat("-------------------------------------------------------------\n\n")
}


## =============================================================================
## 2. SCALING HELPERS
## =============================================================================
##
## Why this section exists.
##
## The simulation has N_MODEL agents, not 3.34 million. Over a 61-day window the
## raw Utah case count can exceed the entire simulated population. An SEIR model
## cannot infect more people than it has: the susceptible pool empties, the
## objective surface goes flat, and the optimiser returns whatever value sits at
## the edge of its search range. The symptom is a fitted curve that saturates
## early plus an implausibly large R0.
##
##   >>> If you shrink the population, you must shrink the incidence too. <<<
##
## R0 is per-capita, so it is invariant under this transformation and remains
## interpretable at true population scale.
## -----------------------------------------------------------------------------

scale_factor <- function() N_MODEL / UTAH_POPULATION

#' Confirmed cases (real world) -> new infections (simulated population)
to_model_scale <- function(daily_cases) {
  stopifnot(is.numeric(daily_cases), all(daily_cases >= 0, na.rm = TRUE))
  (daily_cases / REPORTING_FRACTION) * scale_factor()
}

#' New infections (simulated population) -> confirmed cases (real world)
to_case_scale <- function(model_counts) {
  model_counts / scale_factor() * REPORTING_FRACTION
}

#' Refuse to proceed if the scaled epidemic cannot fit in the population.
check_scaling <- function(model_incidence, max_attack = 0.60) {
  total <- sum(model_incidence)
  ar    <- total / N_MODEL
  cat(sprintf("  Scaling check: %.0f scaled infections over %d days in n=%d",
              total, length(model_incidence), N_MODEL))
  cat(sprintf(" (attack rate %.1f%%)\n", 100 * ar))

  if (ar > max_attack) {
    stop(sprintf(paste0(
      "Scaled epidemic needs an attack rate of %.0f%%, above the %.0f%% ceiling.\n",
      "  Calibration would saturate at the parameter bounds.\n",
      "  Fix by raising N_MODEL, raising REPORTING_FRACTION, or shortening the window."),
      100 * ar, 100 * max_attack))
  }
  if (ar < 0.001) {
    warning(sprintf(paste0(
      "Scaled attack rate is only %.3f%%. Stochastic extinction is likely; ",
      "consider raising N_MODEL."), 100 * ar))
  }
  invisible(ar)
}


## =============================================================================
## 3. R0, THE GROWTH RATE, AND WHY SEIR DIFFERS FROM SIR
## =============================================================================
##
## R0 has the SAME formula in both models: crate * ptran / gamma. What differs
## is how an OBSERVED growth rate r maps onto R0:
##
##     SIR :  R0 = 1 + r/gamma
##     SEIR:  R0 = (1 + r/sigma) * (1 + r/gamma)
##
## The latent period delays transmission, so SEIR needs a LARGER R0 to
## reproduce the same observed growth. If your SEIR and SIR fits disagree on
## R0, that is structural, not a bug. Report which model each came from.
##
## The growth-rate estimate below needs no simulation at all, which makes it a
## useful independent target: if a calibrated R0 lands far from it, suspect the
## fit rather than the data.
## -----------------------------------------------------------------------------

r_to_R0_seir <- function(r, sigma = INCUBATION_RATE, gamma = RECOVERY_RATE) {
  (1 + r / sigma) * (1 + r / gamma)
}

r_to_R0_sir <- function(r, gamma = RECOVERY_RATE) {
  1 + r / gamma
}

#' Exponential growth rate from the early phase, by log-linear regression.
growth_rate <- function(incidence, days = 14L) {
  k  <- min(days, length(incidence))
  y  <- incidence[seq_len(k)]
  t  <- seq_len(k)
  ok <- y > 0
  if (sum(ok) < 3L) return(NA_real_)
  unname(coef(lm(log(y[ok]) ~ t[ok]))[2])
}

#' Symmetric mean absolute percentage error.
smape <- function(obs, sim) {
  d <- (abs(obs) + abs(sim)) / 2
  mean(ifelse(d == 0, 0, abs(obs - sim) / d)) * 100
}


## =============================================================================
## 4. INITIAL CONDITIONS AND PARAMETER LABELS
## =============================================================================

#' Initial prevalence for SEIR, derived from an INCIDENCE series.
#'
#' Two corrections over the old `incidence_vec[1] / N`:
#'
#'   (a) Incidence is not prevalence. The number currently infectious is roughly
#'       incidence x infectious period.
#'   (b) SEIR has an invisible compartment. At the start of an outbreak roughly
#'       incidence x latent period people are already exposed but not yet
#'       infectious. They never appear in case data, but the model needs them.
#'
#' So the seed is incidence x (latent + infectious), about 10x day-one incidence
#' with the default settings -- not 1x. Seeding too small forces the calibration
#' to compensate with an inflated transmission estimate.
initial_prevalence <- function(model_incidence, n_smooth = 3L) {
  k    <- min(n_smooth, length(model_incidence))
  seed <- mean(model_incidence[seq_len(k)]) *
            (INFECTIOUS_PERIOD_DAYS + INCUBATION_DAYS)
  min(max(seed / N_MODEL, 1 / N_MODEL), 0.5)
}

# set_param() matches on the model's parameter LABEL. Verified 2026-08-17
# against the installed epiworldR by calling summary() on a ModelSEIRCONN
# instance directly -- the draft's labels ("Transmission rate", "Recovery
# rate", "Incubation days") were all wrong except contact_rate, and would
# have failed loudly (by design) at the first seir_set_param() call.
SEIR_PARAM <- list(
  contact_rate      = "Contact rate",
  transmission_rate = "Prob. Transmission",
  incubation_days   = "Avg. Incubation days",
  recovery_rate     = "Prob. Recovery"
)

#' set_param() that fails loudly rather than silently on an unknown label.
seir_set_param <- function(model, key, value) {
  label <- SEIR_PARAM[[key]]
  if (is.null(label)) stop("Unknown parameter key: ", key)
  ok <- tryCatch({ set_param(model, label, value); TRUE },
                 error = function(e) FALSE)
  if (!ok) {
    stop(sprintf(paste0(
      "epiworldR rejected parameter label '%s'.\n",
      "  Run summary(model) to see the labels your version uses, then update\n",
      "  SEIR_PARAM in section 4."), label))
  }
  invisible(TRUE)
}


## =============================================================================
## 5. WHICH TRANSITION IS A "CASE"?
## =============================================================================
##
## SIR offers one candidate, S -> I. SEIR offers two, and they are NOT
## interchangeable:
##
##     S -> E   new infection            unobservable -- nobody is tested on
##                                       the day they are infected
##     E -> I   onset of infectiousness  observable  -- this is what
##                                       surveillance approximates
##
## The two series are separated by roughly one latent period. Fitting reported
## cases against S -> E shifts the whole curve about INCUBATION_DAYS days early,
## and the optimiser absorbs the mismatch by inflating transmission.
##
## We use E -> I everywhere.
## -----------------------------------------------------------------------------

CASE_FROM <- "Exposed"
CASE_TO   <- "Infected"

#' Daily E->I counts from a single run, aligned to the observed window.
case_series <- function(model, n_obs) {
  inc <- plot_incidence(model, plot = FALSE)
  if (!"Infected" %in% names(inc)) {
    stop("plot_incidence() has no 'Infected' column; available: ",
         paste(names(inc), collapse = ", "))
  }
  y <- as.numeric(inc$Infected)          # new entries into I, i.e. E -> I

  # Drop the burn-in: all seed agents start Exposed, so the opening days
  # reflect the initial condition rather than transmission.
  if (length(y) > SEED_BURN_IN) y <- y[-seq_len(SEED_BURN_IN)]

  # Force the length to match. Relying on run(ndays) happening to line up with
  # the observed window silently corrupts phase statistics rather than erroring.
  if (length(y) > n_obs)      y <- y[seq_len(n_obs)]
  else if (length(y) < n_obs) y <- c(y, rep(0, n_obs - length(y)))
  y
}


## =============================================================================
## 6. LOAD DATA
## =============================================================================

print_config()
cat("=== Loading data ===\n")

data("utah_covid_data", package = "epiworldRcalibrate")

last_date  <- max(utah_covid_data$Date, na.rm = TRUE)
covid_data <- utah_covid_data |>
  filter(Date > (last_date - N_DAYS_WINDOW)) |>
  arrange(Date)

stopifnot("Daily.Cases" %in% names(covid_data))
if (nrow(covid_data) != N_DAYS_WINDOW) {
  stop(sprintf(paste0(
    "Expected %d rows, got %d -- the series has gaps. Reindex onto a complete\n",
    "  daily date sequence first; a missing day shifts every phase statistic."),
    N_DAYS_WINDOW, nrow(covid_data)))
}

cases_observed <- as.numeric(covid_data$Daily.Cases)
stopifnot(!any(is.na(cases_observed)))

cat(sprintf("  %s to %s\n", min(covid_data$Date), max(covid_data$Date)))
cat(sprintf("  Total confirmed cases : %s\n",
            format(sum(cases_observed), big.mark = ",")))
cat(sprintf("  Mean daily cases      : %.1f\n", mean(cases_observed)))

incidence_vec <- to_model_scale(cases_observed)
check_scaling(incidence_vec)
n_obs       <- length(incidence_vec)
model_ndays <- n_obs + SEED_BURN_IN - 1L

# Model-free reference, computed before any fitting.
r_hat <- growth_rate(incidence_vec, days = 14L)
if (is.finite(r_hat)) {
  cat(sprintf("\n  Early growth rate r   : %.4f /day (doubling %.1f d)\n",
              r_hat, log(2) / r_hat))
  cat(sprintf("  Implied R0 under SEIR : %.2f   <- structural target\n",
              r_to_R0_seir(r_hat)))
  cat(sprintf("  Implied R0 under SIR  : %.2f   (reference only)\n",
              r_to_R0_sir(r_hat)))
}


## =============================================================================
## 7. BUILD THE SEIR MODEL
## =============================================================================

cat("\n=== Building SEIR model ===\n")

init_prev <- initial_prevalence(incidence_vec)
cat(sprintf("  Initial prevalence: %.5f (%.0f agents seeded as Exposed)\n",
            init_prev, init_prev * N_MODEL))

model_seir <- ModelSEIRCONN(
  name              = "Utah COVID-19 (SEIR)",
  n                 = N_MODEL,
  prevalence        = init_prev,
  contact_rate      = 5.0,               # placeholder, calibrated below
  transmission_rate = PTRAN_FIXED,
  incubation_days   = INCUBATION_DAYS,
  recovery_rate     = RECOVERY_RATE
)
verbose_off(model_seir)

cat("\n  Verify these labels match SEIR_PARAM in section 4:\n")
print(summary(model_seir))


## =============================================================================
## 8. METHOD A -- FAST DIRECT CALIBRATION
## =============================================================================
##
## With durations fixed and ptran fixed, one scalar remains: R0. The contact
## rate follows from crate = R0 * gamma / ptran.
##
## A note on the optimiser. ?optim states that one-dimensional Nelder-Mead is
## unreliable -- the simplex degenerates to two points and collapses into
## whichever basin the start lands in. That is why the old code needed six
## restarts and why its `convergence` flag could not be trusted. A coarse grid
## scan followed by golden-section refinement on the bracketing interval is both
## correct and faster for a scalar problem.
## -----------------------------------------------------------------------------

cat("\n=== Method A: direct R0 search ===\n")

objective_R0 <- function(R0, n_rep = 3L) {
  crate <- R0 * RECOVERY_RATE / PTRAN_FIXED
  seir_set_param(model_seir, "contact_rate",      crate)
  seir_set_param(model_seir, "transmission_rate", PTRAN_FIXED)
  # Average over replicates: a single stochastic run makes the objective noisy
  # enough to trap any optimiser in spurious local minima.
  mean(replicate(n_rep, {
    run(model_seir, ndays = model_ndays)
    smape(incidence_vec, case_series(model_seir, n_obs))
  }))
}

t0_direct <- proc.time()[["elapsed"]]

R0_grid <- seq(0.5, 8.0, by = 0.25)
cat(sprintf("  Coarse scan over %d values of R0...\n", length(R0_grid)))
grid_vals <- vapply(R0_grid, objective_R0, numeric(1))

j       <- which.min(grid_vals)
bracket <- c(R0_grid[max(j - 1L, 1L)],
             R0_grid[min(j + 1L, length(R0_grid))])

# Warn if the minimum sits on an edge -- that usually means the scaling is off
# or the search range is too narrow, not that the answer is extreme.
if (j == 1L || j == length(R0_grid)) {
  warning(sprintf(paste0(
    "Grid minimum is at the boundary (R0 = %.2f). Widen R0_grid, or check the\n",
    "  scaling assumptions -- a saturated epidemic pins R0 at the edge."),
    R0_grid[j]))
}

opt_direct <- optimize(objective_R0, interval = bracket, tol = 1e-3)

R0_direct    <- opt_direct$minimum
crate_direct <- R0_direct * RECOVERY_RATE / PTRAN_FIXED
smape_direct <- opt_direct$objective
time_direct  <- proc.time()[["elapsed"]] - t0_direct

cat(sprintf("  R0           : %.3f\n",   R0_direct))
cat(sprintf("  Contact rate : %.3f\n",   crate_direct))
cat(sprintf("  SMAPE        : %.1f%%\n", smape_direct))
cat(sprintf("  Time         : %.1f s\n", time_direct))


## =============================================================================
## 9. METHOD B -- ABC VIA LIKELIHOOD-FREE MCMC
## =============================================================================
##
## Durations are FIXED here, not sampled. With a single incidence curve,
## (crate, ptran, gamma, sigma) are only weakly identifiable: R0 is well
## determined but the individual parameters trade off against each other, which
## is exactly what a wandering trace looks like. SEIR makes this worse, since
## sigma is close to unidentifiable from incidence alone. Fixing the durations
## leaves two free parameters and a well-behaved posterior -- and makes the
## result comparable with Method A and with the BiLSTM, both of which also
## assume them.
##
## Set ESTIMATE_DURATIONS <- TRUE if you deliberately want the 4-parameter
## version, but then do not compare its R0 against the other methods.
## -----------------------------------------------------------------------------

ESTIMATE_DURATIONS <- FALSE

abc_result <- NULL

if (RUN_ABC) {

  cat("\n=== Method B: ABC (LFMCMC) ===\n")

  simulation_fun <- function(params, lfmcmc_obj) {
    seir_set_param(model_seir, "contact_rate",      params[1])
    seir_set_param(model_seir, "transmission_rate", params[2])
    if (ESTIMATE_DURATIONS) {
      seir_set_param(model_seir, "recovery_rate",   params[3])
      seir_set_param(model_seir, "incubation_days", params[4])
    }
    run(model_seir, ndays = model_ndays)
    case_series(model_seir, n_obs)
  }

  # Phase windows derived from the series length. The old code hard-coded
  # 1:20 / 21:40 / 41:60, which silently dropped day 61 and broke outright on
  # any other window length.
  summary_fun <- function(data, lfmcmc_obj) {
    n     <- length(data)
    phase <- cut(seq_len(n), breaks = 3, labels = FALSE)
    c(sum(data),               # total cases
      max(data),               # peak height
      which.max(data),         # day of peak
      mean(data[phase == 1]),  # early
      mean(data[phase == 2]),  # middle
      mean(data[phase == 3]))  # late
  }

  proposal_fun <- function(old_params, lfmcmc_obj) {
    # Step sizes halved from the original (0.10, 0.15) alongside the epsilon
    # increase above -- see ABC_EPSILON comment. Smaller steps mean each
    # proposal is more likely to land within the (now looser) tolerance.
    new_crate <- exp(log(old_params[1]) + rnorm(1, sd = 0.05))
    new_ptran <- plogis(qlogis(old_params[2]) + rnorm(1, sd = 0.08))
    if (ESTIMATE_DURATIONS) {
      return(c(new_crate, new_ptran,
               plogis(qlogis(old_params[3]) + rnorm(1, sd = 0.025)),
               exp(log(old_params[4]) + rnorm(1, sd = 0.10))))
    }
    c(new_crate, new_ptran)
  }

  # Scale-free discrepancy. The old kernel weighted by 1/(|obs|+1), leaving the
  # distance in squared-case units, so epsilon had no interpretable scale and
  # had to be tuned by trial and error. Relative errors make epsilon a
  # fractional tolerance: 0.25 means "within roughly 25% on the statistics".
  kernel_fun <- function(simulated_stat, observed_stat, epsilon, lfmcmc_obj) {
    rel <- (simulated_stat - observed_stat) / pmax(abs(observed_stat), 1e-8)
    exp(-sum(rel^2) / (2 * epsilon^2))
  }

  lfmcmc_obj <- LFMCMC(model_seir)
  lfmcmc_obj <- set_simulation_fun(lfmcmc_obj, simulation_fun)
  lfmcmc_obj <- set_summary_fun(lfmcmc_obj, summary_fun)
  lfmcmc_obj <- set_proposal_fun(lfmcmc_obj, proposal_fun)
  lfmcmc_obj <- set_kernel_fun(lfmcmc_obj, kernel_fun)
  lfmcmc_obj <- set_observed_data(lfmcmc_obj, incidence_vec)

  init_params <- if (ESTIMATE_DURATIONS) {
    c(crate_direct, PTRAN_FIXED, RECOVERY_RATE, INCUBATION_DAYS)
  } else {
    c(crate_direct, PTRAN_FIXED)   # start from the Method A estimate
  }

  cat(sprintf("  Observed summary stats: %s\n",
              paste(round(summary_fun(incidence_vec, NULL), 3), collapse = ", ")))
  cat(sprintf("  %d samples, %d burn-in, epsilon %.2f (relative), %d free params\n",
              ABC_SAMPLES, ABC_BURNIN, ABC_EPSILON, length(init_params)))

  t0_abc <- proc.time()[["elapsed"]]
  run_lfmcmc(
    lfmcmc      = lfmcmc_obj,
    params_init = init_params,
    n_samples   = ABC_SAMPLES,
    epsilon     = ABC_EPSILON,
    seed        = MODEL_SEED
  )
  time_abc <- proc.time()[["elapsed"]] - t0_abc

  accepted <- as.matrix(get_all_accepted_params(lfmcmc_obj))

  # Real acceptance rate. The old code used nrow(accepted)/n_samples, but the
  # chain records one row per iteration whether or not the proposal was
  # accepted -- so it always reported ~100%. Count steps that actually moved.
  acceptance_rate <- mean(rowSums(abs(diff(accepted))) > 0) * 100

  cat(sprintf("  Acceptance rate: %.1f%%", acceptance_rate))
  if (acceptance_rate < 10) {
    cat("  <- too low: raise ABC_EPSILON or shrink the proposal sd\n")
  } else if (acceptance_rate > 60) {
    cat("  <- too high: the chain is barely constrained, lower ABC_EPSILON\n")
  } else {
    cat("\n")
  }

  post <- accepted[(ABC_BURNIN + 1L):nrow(accepted), , drop = FALSE]

  abc_crate <- median(post[, 1])
  abc_ptran <- median(post[, 2])
  abc_recov <- if (ESTIMATE_DURATIONS) median(post[, 3]) else RECOVERY_RATE
  abc_incub <- if (ESTIMATE_DURATIONS) median(post[, 4]) else INCUBATION_DAYS

  # Propagate the joint posterior rather than combining marginal quantiles,
  # which would overstate the interval.
  recov_draws <- if (ESTIMATE_DURATIONS) post[, 3] else RECOVERY_RATE
  R0_draws    <- (post[, 1] * post[, 2]) / recov_draws
  abc_R0      <- median(R0_draws)
  abc_R0_ci   <- quantile(R0_draws, c(0.025, 0.975))

  cat(sprintf("  Contact rate : %.4f [%.4f, %.4f]\n", abc_crate,
              quantile(post[, 1], 0.025), quantile(post[, 1], 0.975)))
  cat(sprintf("  Trans. prob  : %.4f [%.4f, %.4f]\n", abc_ptran,
              quantile(post[, 2], 0.025), quantile(post[, 2], 0.975)))
  cat(sprintf("  R0           : %.3f [%.3f, %.3f]\n",
              abc_R0, abc_R0_ci[1], abc_R0_ci[2]))
  cat(sprintf("  Time         : %.1f s\n", time_abc))

  abc_result <- list(
    crate = abc_crate, ptran = abc_ptran, recov = abc_recov,
    incub = abc_incub, R0 = abc_R0, R0_ci = abc_R0_ci,
    posterior = post, acceptance_rate = acceptance_rate, time = time_abc
  )

  # ---- Trace plots: one panel per parameter --------------------------------
  # The old version overlaid a contact rate of ~5 and rates of ~0.05 on a shared
  # y-axis, so the rate traces were flat lines at the bottom of the plot and
  # showed nothing at all about mixing.
  par_names <- if (ESTIMATE_DURATIONS) {
    c("Contact rate", "Transmission rate", "Recovery rate", "Incubation days")
  } else {
    c("Contact rate", "Transmission rate")
  }
  par_cols <- c("tomato", "darkgreen", "steelblue", "darkorange")

  png(file.path(OUT_DIR, "abc_traces.png"),
      width = 900, height = 220 * length(par_names), res = 110)
  op <- par(mfrow = c(length(par_names), 1), mar = c(4, 4, 2, 1))
  for (jj in seq_along(par_names)) {
    plot(accepted[, jj], type = "l", lwd = 1, col = par_cols[jj],
         main = par_names[jj], xlab = "Step", ylab = "Value")
    abline(v = ABC_BURNIN, lty = 2, col = "grey40")
  }
  par(op)
  dev.off()
  cat("  Saved: abc_traces.png\n")

} else {
  cat("\n=== Method B skipped (RUN_ABC = FALSE) ===\n")
}


## =============================================================================
## 10. BERNARDO BiLSTM PREDICTION (same window as sections 6-9)
## =============================================================================
##
## Note: the BiLSTM shipped with epiworldRcalibrate (calibrate_sir()) was
## trained on SIRCONN trajectories -- do NOT use it here, it answers a
## different question. The prediction below is from the project's own
## SEIR-trained Bernardo network (bernardo_predict_current61.py), run on the
## exact same 61-day window fit in sections 6-9.
## -----------------------------------------------------------------------------

bilstm_R0 <- NA_real_

if (file.exists(BILSTM_CSV)) {
  cat("\n=== BiLSTM prediction (same window) ===\n")
  preds <- read.csv(BILSTM_CSV)

  row_p <- if (BILSTM_WINDOW %in% preds$window) {
    preds[preds$window == BILSTM_WINDOW, ][1L, ]
  } else {
    fb <- preds[which.max(preds$win_len), ][1L, ]
    # Log the fallback. Silently substituting a different window is how a
    # reported R0 ends up describing a period nobody intended.
    cat(sprintf("  NOTE: window '%s' not found; using '%s' (%d d) instead.\n",
                BILSTM_WINDOW, fb$window, fb$win_len))
    fb
  }
  bilstm_R0 <- row_p$R0_pred
  cat(sprintf("  R0 (%s): %.3f\n", row_p$window, bilstm_R0))
} else {
  cat("\n  BiLSTM predictions not found at:", BILSTM_CSV, "\n")
  cat("  -> Run: python real_data/bernardo_predict_current61.py\n")
}


## =============================================================================
## 11. UNCERTAINTY QUANTIFICATION
## =============================================================================
##
## The ribbon below reflects SIMULATION variance -- stochastic variation in who
## contacts whom -- not parameter uncertainty. Method A gives a point estimate,
## so that is all it can express. The ABC posterior in section 9 is what carries
## parameter uncertainty; the R0 credible interval printed there is the honest
## statement of it.
## -----------------------------------------------------------------------------

cat("\n=== Uncertainty quantification ===\n")

R0_final    <- if (!is.null(abc_result)) abc_result$R0    else R0_direct
crate_final <- if (!is.null(abc_result)) abc_result$crate else crate_direct
ptran_final <- if (!is.null(abc_result)) abc_result$ptran else PTRAN_FIXED

seir_set_param(model_seir, "contact_rate",      crate_final)
seir_set_param(model_seir, "transmission_rate", ptran_final)

cat(sprintf("  Running %d simulations at R0 = %.3f...\n", N_SIMS_UQ, R0_final))

saver <- make_saver("transition")
run_multiple(model_seir, ndays = model_ndays, nsims = N_SIMS_UQ,
             saver = saver, nthreads = N_THREADS)
sim_results <- run_multiple_get_results(model_seir, nthreads = 1L,
                                        freader = data.table::fread)
# nthreads=1 here (not N_THREADS): run_multiple_get_results() reads per-rep
# CSVs via a parallel::makeCluster() when nthreads > 1, and those cluster
# workers are fresh R processes that don't inherit this session's custom
# .libPaths()/LD_LIBRARY_PATH, so data.table::fread() fails inside them with
# "could not find function %chin%" (same bug found and fixed 2026-08-17 in
# seir_common.R's seir_run_multi_ci()). Serial reading is fine here since it
# runs once, not per simulation.

quantiles_df <- sim_results$transition |>
  filter(from == CASE_FROM, to == CASE_TO) |>   # E -> I, not S -> E
  group_by(date) |>
  summarize(lower_ci = quantile(counts, 0.025),
            median   = quantile(counts, 0.500),
            upper_ci = quantile(counts, 0.975),
            .groups  = "drop") |>
  arrange(date) |>
  slice(-seq_len(SEED_BURN_IN))                 # drop burn-in, as in the fit

# Align lengths explicitly. Assuming one row per observed day and letting R
# recycle would shift the entire series without any warning.
k <- min(nrow(quantiles_df), n_obs)

plot_df <- tibble(
  Date     = covid_data$Date[seq_len(k)],
  observed = cases_observed[seq_len(k)],
  # Back-transform to the observed case scale. Plotting simulated counts
  # against raw case counts compares two different units.
  lower    = to_case_scale(quantiles_df$lower_ci[seq_len(k)]),
  median   = to_case_scale(quantiles_df$median[seq_len(k)]),
  upper    = to_case_scale(quantiles_df$upper_ci[seq_len(k)])
)

coverage <- mean(plot_df$observed >= plot_df$lower &
                 plot_df$observed <= plot_df$upper) * 100
rmse     <- sqrt(mean((plot_df$median - plot_df$observed)^2))
mae      <- mean(abs(plot_df$median - plot_df$observed))

cat(sprintf("  95%% coverage : %.1f%% (%d of %d days)\n", coverage,
            sum(plot_df$observed >= plot_df$lower &
                plot_df$observed <= plot_df$upper), nrow(plot_df)))
cat(sprintf("  RMSE         : %.1f cases/day\n", rmse))
cat(sprintf("  MAE          : %.1f cases/day\n", mae))


## =============================================================================
## 12. PLOTS
## =============================================================================

g_fit <- ggplot(plot_df, aes(x = Date)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#C62828", alpha = 0.25) +
  geom_line(aes(y = median,   colour = "SEIR median"), linewidth = 1.1) +
  geom_line(aes(y = observed, colour = "Observed"),    linewidth = 1.1) +
  geom_point(aes(y = observed, colour = "Observed"), size = 1.6) +
  scale_colour_manual(values = c("SEIR median" = "#C62828",
                                 "Observed"    = "#1565C0"), name = NULL) +
  labs(
    title    = "Utah COVID-19: observed vs calibrated SEIR",
    subtitle = sprintf(
      "R₀ = %.2f | latent %g d, infectious %g d (assumed) | coverage %.0f%% | MAE %.0f cases/day",
      R0_final, INCUBATION_DAYS, INFECTIOUS_PERIOD_DAYS, coverage, mae),
    x = "Date", y = "Daily confirmed cases"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom",
        plot.title      = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "seir_fit.png"), g_fit,
       width = 10, height = 5.5, dpi = 150)
cat("  Saved: seir_fit.png\n")

# Objective surface -- shows at a glance whether the minimum is a real interior
# optimum or the optimiser hitting a wall.
g_obj <- ggplot(data.frame(R0 = R0_grid, smape = grid_vals), aes(R0, smape)) +
  geom_line(colour = "grey40") +
  geom_point(size = 1.4, colour = "grey30") +
  geom_vline(xintercept = R0_direct, colour = "#C62828", linetype = "dashed") +
  annotate("text", x = R0_direct, y = max(grid_vals),
           label = sprintf(" direct R0 = %.2f", R0_direct),
           hjust = 0, colour = "#C62828", size = 3.5) +
  labs(title = "Objective surface over R0",
       subtitle = "A minimum at the edge of the range means the fit saturated",
       x = "R0", y = "SMAPE (%)") +
  theme_minimal(base_size = 12)

ggsave(file.path(OUT_DIR, "objective_surface.png"), g_obj,
       width = 8, height = 4.5, dpi = 150)
cat("  Saved: objective_surface.png\n")


## =============================================================================
## 13. SUMMARY AND SAVE
## =============================================================================

results <- data.frame(
  method = c("Growth rate (model-free)",
             "Direct search (Method A)",
             "ABC / LFMCMC (Method B)",
             "Bernardo BiLSTM"),
  R0 = c(if (is.finite(r_hat)) r_to_R0_seir(r_hat) else NA_real_,
         R0_direct,
         if (!is.null(abc_result)) abc_result$R0 else NA_real_,
         bilstm_R0),
  R0_lower = c(NA_real_, NA_real_,
               if (!is.null(abc_result)) abc_result$R0_ci[1] else NA_real_,
               NA_real_),
  R0_upper = c(NA_real_, NA_real_,
               if (!is.null(abc_result)) abc_result$R0_ci[2] else NA_real_,
               NA_real_),
  time_s = c(NA_real_, time_direct,
             if (!is.null(abc_result)) abc_result$time else NA_real_,
             NA_real_),
  stringsAsFactors = FALSE
)

cat("\n")
cat("################################################################\n")
cat("  RESULTS\n")
cat("################################################################\n\n")
print(results, row.names = FALSE, digits = 4)

cat("\n  Fit (final parameters):\n")
cat(sprintf("    95%% coverage : %.1f%%\n", coverage))
cat(sprintf("    RMSE         : %.1f cases/day\n", rmse))
cat(sprintf("    MAE          : %.1f cases/day\n", mae))

cat("\n  All R0 values above are conditional on an assumed infectious period\n")
cat(sprintf("  of %g days and latent period of %g days. R0 is linear in the\n",
            INFECTIOUS_PERIOD_DAYS, INCUBATION_DAYS))
cat("  former; change it and every number in this table moves proportionally.\n")

# Sensitivity: same data, same growth rate, different assumed durations.
sens <- expand.grid(latent     = c(2, 3, 5),
                    infectious = c(5, 7, 10)) |>
  mutate(R0 = round(r_to_R0_seir(r_hat, sigma = 1 / latent,
                                 gamma = 1 / infectious), 2),
         pct_vs_baseline = sprintf("%+.0f%%",
                                   100 * (R0 / r_to_R0_seir(r_hat) - 1)))

cat("\n  Sensitivity of R0 to the assumed durations (growth-rate based):\n\n")
print(sens, row.names = FALSE)

write.csv(results, file.path(OUT_DIR, "seir_results.csv"), row.names = FALSE)
write.csv(plot_df, file.path(OUT_DIR, "seir_fitted_curve.csv"), row.names = FALSE)
write.csv(sens,    file.path(OUT_DIR, "seir_sensitivity.csv"), row.names = FALSE)

saveRDS(
  list(
    results         = results,
    fit             = list(coverage_95 = coverage, rmse = rmse, mae = mae),
    direct          = list(R0 = R0_direct, crate = crate_direct,
                           smape = smape_direct, time = time_direct),
    abc             = abc_result,
    bilstm_R0       = bilstm_R0,
    growth_rate     = r_hat,
    case_transition = c(from = CASE_FROM, to = CASE_TO),
    # Store the assumptions alongside the estimates so a later comparison can
    # verify every method used the same ones.
    assumptions = list(
      n_model                = N_MODEL,
      true_population        = UTAH_POPULATION,
      reporting_fraction     = REPORTING_FRACTION,
      infectious_period_days = INFECTIOUS_PERIOD_DAYS,
      incubation_days        = INCUBATION_DAYS,
      recovery_rate          = RECOVERY_RATE,
      ptran_fixed            = PTRAN_FIXED,
      seed_burn_in           = SEED_BURN_IN,
      window_days            = N_DAYS_WINDOW
    )
  ),
  file.path(OUT_DIR, "seir_calibration.rds")
)

cat(sprintf("\n  Output written to: %s/\n", OUT_DIR))
cat("    seir_results.csv, seir_fitted_curve.csv, seir_sensitivity.csv\n")
cat("    seir_calibration.rds, seir_fit.png, objective_surface.png")
if (RUN_ABC) cat(", abc_traces.png")
cat("\n\n")
