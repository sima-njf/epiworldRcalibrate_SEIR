suppressPackageStartupMessages({library(epiworldR); library(ggplot2); library(dplyr)})
PROJECT_DIR <- path.expand("~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction")
REAL_DIR    <- file.path(PROJECT_DIR, "real_data")
PLOTS_DIR   <- file.path(REAL_DIR, "plots")
source(file.path(PROJECT_DIR, "seir_common.R")); seir_set_libpath()

EVAL_REPS <- 1000L; NTHREADS <- 12L
N <- 5000L
beta_old <- 0.1462   # archived model/model_bilstm (1).pt prediction on current61
recov <- 1/7

inc_df <- read.csv(file.path(REAL_DIR, "utah_covid_current61.csv"))
obs    <- as.numeric(inc_df$daily_cases)
dates  <- as.Date(inc_df$date)
ndays  <- length(obs)

crate <- max(beta_old, 1.0); ptran <- beta_old / crate
p <- list(n = N, prevalence = obs[1] / N, incub = 5, recov = recov, crate = crate, ptran = ptran)
cat(sprintf("N=%d  prevalence=%.5f (seed=%.1f)  beta=%.4f  R0=%.4f\n",
    N, p$prevalence, p$prevalence*N, beta_old, beta_old/recov))

ci <- seir_run_multi_ci(p, ndays = ndays, nreps = EVAL_REPS, nthreads = NTHREADS)

rmse <- sqrt(mean((ci$med - obs)^2)); mae <- mean(abs(ci$med - obs))
cov  <- mean(obs >= ci$lower & obs <= ci$upper) * 100
cat(sprintf("RMSE=%.2f  MAE=%.2f  coverage=%.1f%%\n", rmse, mae, cov))

df <- data.frame(date = dates, lower = ci$lower, med = ci$med, upper = ci$upper)
df_obs <- data.frame(date = dates, y = obs)

g <- ggplot() +
  geom_ribbon(data = df, aes(x = date, ymin = lower, ymax = upper), fill = "#D4537E", alpha = 0.25) +
  geom_line(data = df, aes(x = date, y = med), colour = "#D4537E", linewidth = 1.1) +
  geom_line(data = df_obs, aes(x = date, y = y), colour = "black", linewidth = 0.7) +
  geom_point(data = df_obs, aes(x = date, y = y), colour = "black", size = 0.8, alpha = 0.6) +
  labs(title = "Utah COVID-19 Current (2025) -- Archived old BiLSTM, N=5000",
       subtitle = sprintf("Black = observed | Pink = 95%% CI (%d SEIR runs) | beta=%.4f, R0=%.3f | RMSE=%.1f MAE=%.1f coverage=%.1f%%",
                           EVAL_REPS, beta_old, beta_old/recov, rmse, mae, cov),
       x = "Date", y = "Daily incidence") +
  theme_bw(base_size = 12) + theme(plot.title = element_text(face = "bold"))

out_png <- file.path(PLOTS_DIR, "current61_oldmodel_N5000.png")
ggsave(out_png, g, width = 11, height = 6, dpi = 150)
cat(sprintf("Saved: %s\n", out_png))
