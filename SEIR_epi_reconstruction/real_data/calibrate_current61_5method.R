# =============================================================================
#  real_data/calibrate_current61_5method.R
#
#  Same 5-method comparison as calibrate_real_5method.R (ABC, ABC-SMC,
#  Nelder-Mead, DE, BiLSTM), but on the CURRENT Utah COVID window (as
#  packaged today: 2025-03-07 to 2025-05-06), not the 2020-21 wave1 data --
#  and with EVAL_REPS overridden to 1000 (per request) instead of 2000.
#
#  Sources calibrate_real_5method.R with SKIP_WAVE1=1 SKIP_MEASLES=1 so only
#  the function definitions run, then calls process_dataset() itself for the
#  current-window files. This re-uses the exact same (already-fixed)
#  calibration code -- N_SCALED/N_BILSTM from seir_scale_config.json, Brent
#  for the 1-D NM case, nthreads=1 for run_multiple_get_results() -- rather
#  than duplicating it.
#
#  Prerequisites:
#    (utah_covid_current61.csv, utah_covid_current61_meta.csv,
#     bernardo_real_covid_current61_predictions.csv already generated)
#
#  Usage:
#    SKIP_WAVE1=1 SKIP_MEASLES=1 Rscript real_data/calibrate_current61_5method.R
# =============================================================================

Sys.setenv(SKIP_WAVE1 = "1", SKIP_MEASLES = "1")
source("real_data/calibrate_real_5method.R")

# Override the run settings for this window (functions look these up from
# the global env at call time, so reassigning here after sourcing is safe).
WINDOWS   <- c(61L)   # the current window is exactly 61 days, no early/mid/late split
EVAL_REPS <- 1000L    # per request: 1000 SEIRconn reps per method, not 2000

process_dataset(
  label      = "Utah COVID-19 Current (2025)",
  inc_csv    = file.path(REAL_DIR, "utah_covid_current61.csv"),
  meta_csv   = file.path(REAL_DIR, "utah_covid_current61_meta.csv"),
  bilstm_csv = file.path(REAL_DIR, "bernardo_real_covid_current61_predictions.csv")
)

cat("\nDone. Plots saved to real_data/plots/, params to real_data/utah_covid_19_current_2025_5method_params.csv\n")
