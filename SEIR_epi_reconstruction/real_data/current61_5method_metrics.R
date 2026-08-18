# =============================================================================
#  real_data/current61_5method_metrics.R
#
#  Coverage/RMSE/MAE comparison for the 5-method current-window (2025) run,
#  same pattern as covid_5method_metrics.R but pointed at the current-61
#  dataset. Reuses the already-found beta values (no re-calibration), just
#  re-runs each method's 1000-rep CI ensemble and builds the comparison plot.
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
# -- see the comment in covid_5method_metrics.R / n_classical_comment in
# seir_scale_config.json for why (bug found 2026-08-17: this script kept its
# own hardcoded N_SCALED=100000 after the main pipeline was fixed to N=5000,
# so it kept reproducing the sawtooth/quantization bug even after the "fix").
scale_cfg <- jsonlite::fromJSON(file.path(REAL_DIR, "seir_scale_config.json"))
N_SCALED  <- as.integer(scale_cfg$n_classical)
EVAL_REPS <- 1000L
NTHREADS  <- 12L
WIN_LEN   <- 61L

METHOD_ORDER  <- c("BiLSTM", "ABC", "ABC-SMC", "NelderMead", "DE")
METHOD_COLORS <- SEIR_METHOD_COLORS[METHOD_ORDER]

inc_df <- read.csv(file.path(REAL_DIR, "utah_covid_current61.csv"))
meta   <- read.csv(file.path(REAL_DIR, "utah_covid_current61_meta.csv"))[1L, ]
obs    <- as.numeric(inc_df$daily_cases)
dates  <- as.Date(inc_df$date)
ndays  <- length(obs)

n_true <- as.numeric(meta$n_pop)
n_use  <- if (n_true > N_SCALED) N_SCALED else n_true

# obs is NOT scaled down by n_use/n_true; prevalence is seeded directly from
# the observed window on the n_use scale. See covid_5method_metrics.R.
known <- list(n = n_use, prevalence = obs[1] / n_use,
              incub = as.numeric(meta$incub_days), recov = as.numeric(meta$recov_rate))

params_df <- read.csv(file.path(REAL_DIR, "utah_covid_19_current__2025__5method_params.csv"))

beta_to_p <- function(beta, known) {
  crate <- max(as.numeric(beta), 1.0)
  ptran <- as.numeric(beta) / crate
  list(n = known$n, prevalence = known$prevalence, incub = known$incub,
       recov = known$recov, crate = crate, ptran = ptran)
}

results <- list()
for (i in seq_len(nrow(params_df))) {
  meth <- params_df$method[i]
  beta <- params_df$beta[i]
  cat(sprintf("%-12s beta=%.6f  running %d-rep ensemble...\n", meth, beta, EVAL_REPS))
  p  <- beta_to_p(beta, known)
  ci <- tryCatch(seir_run_multi_ci(p, ndays = ndays, nreps = EVAL_REPS, nthreads = NTHREADS),
                 error = function(e) { cat("  FAILED:", conditionMessage(e), "\n"); NULL })
  if (!is.null(ci)) {
    results[[meth]] <- list(beta = beta, R0 = beta / known$recov, ci = ci)
  }
}

metrics_rows <- list()
for (meth in names(results)) {
  ci <- results[[meth]]$ci
  in_ci <- obs >= ci$lower & obs <= ci$upper
  metrics_rows[[meth]] <- data.frame(
    method = meth, beta = results[[meth]]$beta, R0 = results[[meth]]$R0,
    coverage_pct = mean(in_ci) * 100, rmse = sqrt(mean((ci$med - obs)^2)),
    mae = mean(abs(ci$med - obs)), n_in_ci = sum(in_ci), n_days = ndays
  )
}
metrics_df <- bind_rows(metrics_rows) |> arrange(match(method, METHOD_ORDER))

timing_df <- read.csv(file.path(REAL_DIR, "utah_covid_current61_method_timing.csv"))
metrics_df <- merge(metrics_df, timing_df, by = "method", all.x = TRUE) |>
  arrange(match(method, METHOD_ORDER))

write.csv(metrics_df, file.path(REAL_DIR, "utah_covid_current61_5method_metrics.csv"), row.names = FALSE)
cat("\n========== Coverage / RMSE / MAE / timing, current 61-day window ==========\n")
print(metrics_df, row.names = FALSE, digits = 4)

ci_long <- bind_rows(lapply(names(results), function(meth) {
  data.frame(method = meth, date = dates, day = seq_len(ndays), observed = obs,
             lower = results[[meth]]$ci$lower, med = results[[meth]]$ci$med,
             upper = results[[meth]]$ci$upper)
}))
write.csv(ci_long, file.path(REAL_DIR, "utah_covid_current61_5method_ci_curves.csv"), row.names = FALSE)

df_ribbons <- bind_rows(lapply(names(results), function(meth) {
  ci <- results[[meth]]$ci
  data.frame(date = dates, lower = ci$lower, med = ci$med, upper = ci$upper, method = meth)
})) |> mutate(method = factor(method, levels = METHOD_ORDER))
df_obs <- data.frame(date = dates, y = obs)

g <- ggplot() +
  geom_ribbon(data = df_ribbons, aes(x = date, ymin = lower, ymax = upper, fill = method), alpha = 0.18) +
  geom_line(data = df_ribbons, aes(x = date, y = med, colour = method), linewidth = 0.9) +
  geom_line(data = df_obs, aes(x = date, y = y), colour = "black", linewidth = 0.7) +
  geom_point(data = df_obs, aes(x = date, y = y), colour = "black", size = 0.8, alpha = 0.6) +
  scale_fill_manual(values = METHOD_COLORS, name = "Method") +
  scale_colour_manual(values = METHOD_COLORS, name = "Method") +
  labs(title = "Utah COVID-19 (current, 2025-03-07 to 2025-05-06) — 5-method SEIR comparison",
       subtitle = sprintf("Black = observed daily cases | Ribbons = 95%% CI (%d SEIR runs, N=%s) | full window calibrated",
                           EVAL_REPS, format(N_SCALED, big.mark = ",")),
       x = "Date", y = "Daily incidence") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"))

out_png <- file.path(PLOTS_DIR, "utah_covid_current61_5method_metrics.png")
ggsave(out_png, g, width = 11, height = 6, dpi = 150)
cat(sprintf("\nSaved plot: %s\n", out_png))
