suppressPackageStartupMessages({library(epiworldR); library(ggplot2); library(dplyr); library(jsonlite)})
PROJECT_DIR <- path.expand("~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction")
REAL_DIR    <- file.path(PROJECT_DIR, "real_data")
PLOTS_DIR   <- file.path(REAL_DIR, "plots")
source(file.path(PROJECT_DIR, "seir_common.R")); seir_set_libpath()

scale_cfg <- jsonlite::fromJSON(file.path(REAL_DIR, "seir_scale_config.json"))
N_SCALED  <- as.integer(scale_cfg$n_classical)   # 7000, same scale as the rest of the pipeline
EVAL_REPS <- 1000L; NTHREADS <- 12L

# Old-model predictions from bernardo_predict_oldmodel.py (run 2026-08-17)
old_preds <- list(
  current61     = list(beta = 0.1462, R0 = 1.1208, csv = "utah_covid_current61.csv", win = 61L,
                        label = "Utah COVID-19 Current (2025)", slug = "utah_covid_19_current__2025_"),
  wave1_first60 = list(beta = 0.1951, R0 = 1.4830, csv = "utah_covid_wave1.csv", win = 60L,
                        label = "Utah COVID-19 Wave 1", slug = "utah_covid_19_wave_1")
)
new_beta <- list(current61 = 0.1352, wave1_first60 = 0.1208)

run_ci <- function(beta, known, ndays) {
  crate <- max(beta, 1.0); ptran <- beta / crate
  p <- list(n = known$n, prevalence = known$prevalence, incub = known$incub,
            recov = known$recov, crate = crate, ptran = ptran)
  seir_run_multi_ci(p, ndays = ndays, nreps = EVAL_REPS, nthreads = NTHREADS)
}

for (key in names(old_preds)) {
  cfg <- old_preds[[key]]
  inc_df <- read.csv(file.path(REAL_DIR, cfg$csv))
  obs_full <- as.numeric(inc_df$daily_cases)
  dates_full <- if ("date" %in% names(inc_df)) as.Date(inc_df$date) else NULL
  obs <- if (key == "wave1_first60") obs_full[1:60] else obs_full
  dates <- if (!is.null(dates_full) && key == "wave1_first60") dates_full[1:60] else dates_full
  ndays <- length(obs)
  known <- list(n = N_SCALED, prevalence = obs[1] / N_SCALED, incub = 5, recov = 1/7)

  cat(sprintf("\n=== %s ===\n", key))
  ci_new <- run_ci(new_beta[[key]], known, ndays)
  ci_old <- run_ci(cfg$beta, known, ndays)

  df <- bind_rows(
    data.frame(date = if (!is.null(dates)) dates else seq_len(ndays),
               lower = ci_new$lower, med = ci_new$med, upper = ci_new$upper,
               model = sprintf("Current tuned model (beta=%.4f, R0=%.3f)", new_beta[[key]], new_beta[[key]]/known$recov)),
    data.frame(date = if (!is.null(dates)) dates else seq_len(ndays),
               lower = ci_old$lower, med = ci_old$med, upper = ci_old$upper,
               model = sprintf("Archived old model (beta=%.4f, R0=%.3f)", cfg$beta, cfg$R0))
  )
  df_obs <- data.frame(date = if (!is.null(dates)) dates else seq_len(ndays), y = obs)

  rmse <- function(ci) sqrt(mean((ci$med - obs)^2))
  mae  <- function(ci) mean(abs(ci$med - obs))
  cov  <- function(ci) mean(obs >= ci$lower & obs <= ci$upper) * 100
  cat(sprintf("  Current tuned: RMSE=%.2f MAE=%.2f coverage=%.1f%%\n", rmse(ci_new), mae(ci_new), cov(ci_new)))
  cat(sprintf("  Archived old : RMSE=%.2f MAE=%.2f coverage=%.1f%%\n", rmse(ci_old), mae(ci_old), cov(ci_old)))

  g <- ggplot() +
    geom_ribbon(data = df, aes(x = date, ymin = lower, ymax = upper, fill = model), alpha = 0.25) +
    geom_line(data = df, aes(x = date, y = med, colour = model), linewidth = 1) +
    geom_line(data = df_obs, aes(x = date, y = y), colour = "black", linewidth = 0.7) +
    geom_point(data = df_obs, aes(x = date, y = y), colour = "black", size = 0.8, alpha = 0.6) +
    labs(title = sprintf("%s -- Archived old BiLSTM vs current tuned BiLSTM", cfg$label),
         subtitle = "Black = observed daily cases | Ribbons = 95% CI (1000 SEIR runs, N matches rest of pipeline)",
         x = "Date", y = "Daily incidence", colour = "", fill = "") +
    theme_bw(base_size = 12) + theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

  out_png <- file.path(PLOTS_DIR, sprintf("%s_bilstm_old_vs_new.png", key))
  ggsave(out_png, g, width = 11, height = 6, dpi = 150)
  cat(sprintf("  Saved: %s\n", out_png))
}
