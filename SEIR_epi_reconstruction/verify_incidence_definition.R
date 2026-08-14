# =============================================================================
#  verify_incidence_definition.R  — RUN THIS FIRST
#
#  Answers two questions that silently invalidate every downstream metric:
#
#    Q1. Does test_incidence_raw.npy record S->E (new exposures) or
#        E->I (new infectious)? Everything downstream now assumes E->I.
#
#    Q2. Does row i of test_incidence_raw.npy correspond to row i of
#        test_actual_parameters.csv?
#
#  Method: for a handful of sims, simulate with the TRUE parameters, extract
#  both conventions, and see which one matches the stored curve better in
#  peak timing / peak size / correlation. Because the ABM is stochastic the
#  match is never exact — we average over several replicates and compare the
#  two conventions relative to each other, not in absolute terms.
#
#  Usage: Rscript verify_incidence_definition.R
# =============================================================================

source("seir_common.R")
seir_set_libpath()

suppressPackageStartupMessages({
  library(epiworldR)
  library(data.table)
  library(reticulate)
})

PROJECT_DIR <- path.expand(
  "~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction"
)

N_CHECK <- 8L    # how many sims to check
N_REPS  <- 25L   # replicates per sim (median curve is the comparison target)

# =============================================================================
# Load
# =============================================================================

actual <- read.csv(file.path(PROJECT_DIR, "test_actual_parameters.csv"))

np_module <- reticulate::import("numpy")
inc_raw   <- as.matrix(
  np_module$load(path.expand(file.path(PROJECT_DIR, "test_incidence_raw.npy")))
)

cat(sprintf("actual        : %d rows\n", nrow(actual)))
cat(sprintf("incidence .npy: %d rows x %d cols\n", nrow(inc_raw), ncol(inc_raw)))

stopifnot(nrow(inc_raw) == nrow(actual))
stopifnot(ncol(inc_raw) == SEIR_NDAYS)

cat("\nParameter ranges in the test set:\n")
print(summary(actual[, intersect(c("beta", "recov", "incub", "n",
                                   "prevalence", "R0"), names(actual))]))

cat(sprintf("\nOptimiser bounds are [%.3f, %.3f].\n", SEIR_BETA_LO, SEIR_BETA_HI))
if (max(actual$beta) > SEIR_BETA_HI || min(actual$beta) < SEIR_BETA_LO) {
  cat("  *** WARNING: some true beta values fall OUTSIDE the bounds. ***\n")
  cat("  *** DE cannot possibly recover them. Widen SEIR_BETA_* in   ***\n")
  cat("  *** seir_common.R before running anything else.             ***\n")
} else {
  cat("  All true beta values are inside the bounds.\n")
}

# =============================================================================
# Both conventions from one model object
# =============================================================================

incidence_both <- function(model, ndays) {

  tm <- as.data.frame(epiworldR::get_hist_transition_matrix(model))
  cc <- seir_transition_cols(tm)   # handles state_from/state_to vs from/to

  from <- as.character(tm[[cc$from]])
  to   <- as.character(tm[[cc$to]])
  date <- as.integer(tm[[cc$date]])
  cnt  <- as.numeric(tm[[cc$counts]])

  grab <- function(f_pat, t_pat) {
    keep <- grepl(f_pat, from, ignore.case = TRUE) &
            grepl(t_pat, to,   ignore.case = TRUE) &
            date >= 1L & date <= ndays
    out <- numeric(ndays)
    if (any(keep)) {
      agg <- tapply(cnt[keep], date[keep], sum)
      out[as.integer(names(agg))] <- as.numeric(agg)
    }
    out
  }

  list(
    s2e = grab("^suscep", "^expos"),
    e2i = grab("^expos",  "^infect")
  )
}

# =============================================================================
# Check
# =============================================================================

set.seed(1)
rows <- sort(sample(seq_len(nrow(actual)), min(N_CHECK, nrow(actual))))

report <- vector("list", length(rows))

for (k in seq_along(rows)) {
  i   <- rows[k]
  s   <- actual[i, ]
  obs <- as.numeric(inc_raw[i, ])

  s2e_mat <- matrix(0, nrow = N_REPS, ncol = SEIR_NDAYS)
  e2i_mat <- matrix(0, nrow = N_REPS, ncol = SEIR_NDAYS)

  for (r in seq_len(N_REPS)) {
    m <- seir_build_model(s$n, s$prevalence, s$beta, s$incub, s$recov)
    run(m, ndays = SEIR_NDAYS + 1L, seed = as.integer(1000 * i + r))
    both <- incidence_both(m, SEIR_NDAYS)
    s2e_mat[r, ] <- both$s2e
    e2i_mat[r, ] <- both$e2i
  }

  s2e_med <- apply(s2e_mat, 2, median)
  e2i_med <- apply(e2i_mat, 2, median)

  report[[k]] <- data.frame(
    row            = i,
    sim_idx        = s$sim_idx,
    true_beta      = s$beta,
    obs_peak_day   = which.max(obs),
    obs_peak_size  = max(obs),
    obs_total      = sum(obs),
    s2e_peak_day   = which.max(s2e_med),
    e2i_peak_day   = which.max(e2i_med),
    s2e_peak_size  = max(s2e_med),
    e2i_peak_size  = max(e2i_med),
    s2e_total      = sum(s2e_med),
    e2i_total      = sum(e2i_med),
    cor_s2e        = seir_cor(obs, s2e_med),
    cor_e2i        = seir_cor(obs, e2i_med),
    peakday_err_s2e = abs(which.max(s2e_med) - which.max(obs)),
    peakday_err_e2i = abs(which.max(e2i_med) - which.max(obs)),
    stringsAsFactors = FALSE
  )

  cat(sprintf("  [%d/%d] row %d  obs peak day %d | S->E %d | E->I %d\n",
              k, length(rows), i, which.max(obs),
              which.max(s2e_med), which.max(e2i_med)))
}

rep_df <- do.call(rbind, report)

cat("\n========== Convention check ==========\n")
print(rep_df[, c("row", "sim_idx", "obs_peak_day",
                 "s2e_peak_day", "e2i_peak_day",
                 "cor_s2e", "cor_e2i")], digits = 3, row.names = FALSE)

cat(sprintf("\nMean |peak-day error| vs S->E : %.2f days\n",
            mean(rep_df$peakday_err_s2e)))
cat(sprintf("Mean |peak-day error| vs E->I : %.2f days\n",
            mean(rep_df$peakday_err_e2i)))
cat(sprintf("Mean correlation with S->E    : %.4f\n", mean(rep_df$cor_s2e, na.rm = TRUE)))
cat(sprintf("Mean correlation with E->I    : %.4f\n", mean(rep_df$cor_e2i, na.rm = TRUE)))

better <- if (mean(rep_df$peakday_err_e2i) <= mean(rep_df$peakday_err_s2e)) "E->I" else "S->E"
cat(sprintf("\n=> test_incidence_raw.npy looks like %s.\n", better))
if (better != "E->I") {
  cat("   *** The pipeline assumes E->I. Regenerate the .npy with the   ***\n")
  cat("   *** E->I convention, or change seir_incidence_e2i() to match. ***\n")
}

# =============================================================================
# Row alignment: correct pairing should beat shuffled pairing by a lot
# =============================================================================

cat("\n========== Row alignment check ==========\n")

aligned  <- mean(abs(rep_df$obs_peak_day - rep_df$e2i_peak_day))
shuffled <- mean(abs(rep_df$obs_peak_day - sample(rep_df$e2i_peak_day)))

cat(sprintf("Mean |peak-day diff|, aligned pairing : %.2f days\n", aligned))
cat(sprintf("Mean |peak-day diff|, shuffled pairing: %.2f days\n", shuffled))
if (aligned < shuffled) {
  cat("=> Rows appear correctly aligned.\n")
} else {
  cat("=> *** Rows do NOT appear aligned. Fix this before anything else. ***\n")
}

write.csv(rep_df, file.path(PROJECT_DIR, "verify_incidence_report.csv"),
          row.names = FALSE)
cat("\nSaved: verify_incidence_report.csv\n")
