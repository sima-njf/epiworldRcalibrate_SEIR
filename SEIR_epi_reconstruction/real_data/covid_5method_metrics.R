# =============================================================================
#  real_data/covid_5method_metrics.R
#
#  SIR-vignette-style writeup for the Utah COVID-19 5-method SEIR comparison:
#  builds the CI-ribbon plot + a coverage/RMSE/MAE table per method, reusing
#  the beta values already found by real_data/calibrate_real_5method.R (no
#  re-calibration -- just re-runs each method's 2000-rep CI ensemble).
#
#  Usage (submit via SLURM, not on the login node):
#    sbatch real_data/covid_5method_metrics_slurm.sh
# =============================================================================

suppressPackageStartupMessages({
  library(epiworldR)
  library(ggplot2)
  library(dplyr)
})

PROJECT_DIR <- path.expand("~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction")
REAL_DIR    <- file.path(PROJECT_DIR, "real_data")
PLOTS_DIR   <- file.path(REAL_DIR, "plots")

source(file.path(PROJECT_DIR, "seir_common.R"))
seir_set_libpath()

# N_SCALED and prevalence-seeding must match calibrate_real_5method.R exactly
# -- this script re-simulates the CI ensemble for the SAME already-found beta
# values, so using a different population/prevalence here would silently
# make the ensemble inconsistent with the params it's plotting against (this
# is exactly what happened 2026-08-17: this script kept its own hardcoded
# N_SCALED=100000 after calibrate_real_5method.R was fixed to N=5000, so it
# kept reproducing the sawtooth/quantization bug even after the "fix").
scale_cfg <- jsonlite::fromJSON(file.path(REAL_DIR, "seir_scale_config.json"))
N_SCALED  <- as.integer(scale_cfg$n_classical)
EVAL_REPS <- 2000L
NTHREADS  <- 12L
WIN_LEN   <- 60L

METHOD_ORDER  <- c("BiLSTM", "ABC", "ABC-SMC", "NelderMead", "DE")
METHOD_COLORS <- SEIR_METHOD_COLORS[METHOD_ORDER]

# -- Load data and the already-found beta values -----------------------------
inc_df <- read.csv(file.path(REAL_DIR, "utah_covid_wave1.csv"))
meta   <- read.csv(file.path(REAL_DIR, "utah_covid_meta.csv"))[1L, ]
obs    <- as.numeric(inc_df$daily_cases)
dates  <- as.Date(inc_df$date)
ndays  <- length(obs)

n_true <- as.numeric(meta$n_pop)
n_use  <- if (n_true > N_SCALED) N_SCALED else n_true

# obs is NOT scaled down by n_use/n_true -- see calibrate_real_5method.R's
# process_dataset() for the full rationale. prevalence is seeded directly
# from the observed window on the n_use scale (same convention as the SIR
# vignette), not from meta$prevalence (a fraction of the true population).
known <- list(n = n_use, prevalence = obs[1] / n_use,
              incub = as.numeric(meta$incub_days), recov = as.numeric(meta$recov_rate))

params_df <- read.csv(file.path(REAL_DIR, "utah_covid_19_wave_1_5method_params.csv"))
params_df <- params_df[params_df$win_len == WIN_LEN, ]

beta_to_p <- function(beta, known) {
  crate <- max(as.numeric(beta), 1.0)
  ptran <- as.numeric(beta) / crate
  list(n = known$n, prevalence = known$prevalence, incub = known$incub,
       recov = known$recov, crate = crate, ptran = ptran)
}

# -- Re-run each method's CI ensemble at its already-found beta --------------
results <- list()
for (i in seq_len(nrow(params_df))) {
  meth <- params_df$method[i]
  beta <- params_df$beta[i]
  cat(sprintf("%-12s beta=%.4f  running %d-rep ensemble...\n", meth, beta, EVAL_REPS))
  p  <- beta_to_p(beta, known)
  ci <- tryCatch(seir_run_multi_ci(p, ndays = ndays, nreps = EVAL_REPS, nthreads = NTHREADS),
                 error = function(e) { cat("  FAILED:", conditionMessage(e), "\n"); NULL })
  if (!is.null(ci)) {
    results[[meth]] <- list(beta = beta, R0 = beta / known$recov, ci = ci)
  }
}

# -- Per-method coverage / RMSE / MAE, computed over two periods -------------
#    "in_sample"  = the WIN_LEN calibration window itself (this is what the
#                    SIR vignette reports: fit measured on the days the model
#                    actually saw).
#    "full_wave"  = the full 365-day wave, i.e. how far a *constant-parameter*
#                    model fit to only the first 60 days extrapolates. Low
#                    coverage here is expected/correct, not a method failure:
#                    Utah's real transmission rate changed over the year
#                    (interventions, behavior, variants) while beta here is
#                    held fixed, so this number answers a different question
#                    than in_sample and the two should not be conflated.
compute_metrics <- function(idx, label) {
  rows <- list()
  for (meth in names(results)) {
    ci <- results[[meth]]$ci
    o  <- obs[idx]; lo <- ci$lower[idx]; me <- ci$med[idx]; up <- ci$upper[idx]
    in_ci <- o >= lo & o <= up
    rows[[meth]] <- data.frame(
      period = label, method = meth, beta = results[[meth]]$beta, R0 = results[[meth]]$R0,
      coverage_pct = mean(in_ci) * 100, rmse = sqrt(mean((me - o)^2)), mae = mean(abs(me - o)),
      n_in_ci = sum(in_ci), n_days = length(idx)
    )
  }
  bind_rows(rows)
}
metrics_df <- bind_rows(
  compute_metrics(seq_len(WIN_LEN), "in_sample_60d"),
  compute_metrics(seq_len(ndays),   "full_wave_365d")
) |> arrange(period, match(method, METHOD_ORDER))
write.csv(metrics_df, file.path(REAL_DIR, "utah_covid_5method_metrics.csv"), row.names = FALSE)

# Also save the raw per-day CI curves so metrics can be recomputed without
# re-running simulations.
ci_long <- bind_rows(lapply(names(results), function(meth) {
  data.frame(method = meth, date = dates, day = seq_len(ndays), observed = obs,
             lower = results[[meth]]$ci$lower, med = results[[meth]]$ci$med,
             upper = results[[meth]]$ci$upper)
}))
write.csv(ci_long, file.path(REAL_DIR, "utah_covid_5method_ci_curves.csv"), row.names = FALSE)

cat("\n========== Coverage / RMSE / MAE by method ==========\n")
print(metrics_df, row.names = FALSE, digits = 4)

# -- Combined plot: all 5 ribbons + observed data -----------------------------
ribbon_rows <- list()
for (meth in names(results)) {
  ci <- results[[meth]]$ci
  ribbon_rows[[meth]] <- data.frame(date = dates, lower = ci$lower, med = ci$med,
                                     upper = ci$upper, method = meth)
}
df_ribbons <- bind_rows(ribbon_rows) |> mutate(method = factor(method, levels = METHOD_ORDER))
df_obs     <- data.frame(date = dates, y = obs)

g <- ggplot() +
  geom_ribbon(data = df_ribbons, aes(x = date, ymin = lower, ymax = upper, fill = method), alpha = 0.18) +
  geom_line(data = df_ribbons, aes(x = date, y = med, colour = method), linewidth = 0.9) +
  geom_vline(xintercept = dates[WIN_LEN], linetype = "dashed", colour = "grey40", linewidth = 0.6) +
  annotate("text", x = dates[WIN_LEN], y = Inf, vjust = 1.5, hjust = -0.1, size = 3,
           colour = "grey40", label = sprintf("day %d (calibration cutoff)", WIN_LEN)) +
  geom_line(data = df_obs, aes(x = date, y = y), colour = "black", linewidth = 0.7) +
  geom_point(data = df_obs, aes(x = date, y = y), colour = "black", size = 0.6, alpha = 0.5) +
  scale_fill_manual(values = METHOD_COLORS, name = "Method") +
  scale_colour_manual(values = METHOD_COLORS, name = "Method") +
  labs(title = "Utah COVID-19 (365-day wave) — 5-method SEIR calibration, 60-day window",
       subtitle = sprintf("Black = observed daily cases | Ribbons = 95%% CI (%d SEIR runs, N=%s)",
                           EVAL_REPS, format(N_SCALED, big.mark = ",")),
       x = "Date", y = "Daily incidence") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"))

out_png <- file.path(PLOTS_DIR, "utah_covid_5method_metrics.png")
ggsave(out_png, g, width = 11, height = 6, dpi = 150)
cat(sprintf("\nSaved plot: %s\n", out_png))
cat(sprintf("Saved metrics: %s\n", file.path(REAL_DIR, "utah_covid_5method_metrics.csv")))
