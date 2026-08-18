# =============================================================================
#  real_data/calibrate_measles1861_5method.R
#
#  Same 5-method comparison (ABC, ABC-SMC, Nelder-Mead, DE, BiLSTM) as
#  calibrate_real_5method.R, run on the EpiEstim Measles1861 smoke-test
#  dataset. Per instructions: use the FULL 48-day curve, not the usual
#  60/90/180/365-day windows (this outbreak is only 48 days long).
#
#  Sources calibrate_real_5method.R with SKIP_WAVE1=1 SKIP_MEASLES=1 so only
#  the function definitions run, then calls process_dataset() itself for
#  this dataset's files -- reuses the exact same (already-fixed) calibration
#  code: N_SCALED from seir_scale_config.json, no obs-scaling round-trip,
#  prevalence seeded as obs[1]/n_use, fixed ABC-SMC generation schedule.
#
#  Prerequisites:
#    Rscript real_data/prepare_measles1861_epiestim.R
#    python  real_data/bernardo_predict_measles1861.py
#
#  Usage:
#    SKIP_WAVE1=1 SKIP_MEASLES=1 Rscript real_data/calibrate_measles1861_5method.R
# =============================================================================

Sys.setenv(SKIP_WAVE1 = "1", SKIP_MEASLES = "1")
source("real_data/calibrate_real_5method.R")

# Full 48-day curve, not a sub-window -- this outbreak is only 48 days long.
WINDOWS   <- c(48L)
EVAL_REPS <- 1000L

# ABC-SMC budget: n=187 evaluates ~1ms/call (vs ~0.1-0.3s at N=7000), so the
# wave1/current61-tuned SMC_ATT_MULT=50/SMC_BUD=20000 stalls out too early
# here (diagnosed 2026-08-17: stalled gen 5/7 at the default budget). Cost
# still explodes each generation as eps nears the noise floor (553 -> 1258
# -> 3814 -> 11275 -> 23517 evals/gen observed even at 300x/300k), so rather
# than chasing full convergence, give it enough headroom to complete several
# real improvement rounds and rely on calibrate_smc()'s graceful-degradation
# fallback (returns the last completed generation) for whichever generation
# it stalls on.
SMC_NG       <- 10L
SMC_BUD      <- 100000L
SMC_ATT_MULT <- 100L

process_dataset(
  label      = "Measles 1861 (EpiEstim)",
  inc_csv    = file.path(REAL_DIR, "measles1861_epiestim.csv"),
  meta_csv   = file.path(REAL_DIR, "measles1861_epiestim_meta.csv"),
  bilstm_csv = file.path(REAL_DIR, "bernardo_real_measles1861_predictions.csv")
)

cat("\nDone. Plots saved to real_data/plots/, params to real_data/measles_1861__epiestim__5method_params.csv\n")
