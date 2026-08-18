suppressPackageStartupMessages({library(epiworldR); library(ggplot2); library(dplyr)})
PROJECT_DIR <- path.expand("~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction")
REAL_DIR    <- file.path(PROJECT_DIR, "real_data")
PLOTS_DIR   <- file.path(REAL_DIR, "plots")
source(file.path(PROJECT_DIR, "seir_common.R")); seir_set_libpath()

EVAL_REPS <- 500L; NTHREADS <- 12L
N_VALUES  <- c(5000L, 6000L, 7000L, 8000L, 9000L, 10000L)

# Predictions from the two models on current61, both already computed:
betas <- list(
  current_tuned = 0.1352,   # bernardo_model/model_bilstm_tuned.pt
  archived_old  = 0.1462    # model/model_bilstm (1).pt
)
recov <- 1/7

inc_df <- read.csv(file.path(REAL_DIR, "utah_covid_current61.csv"))
obs    <- as.numeric(inc_df$daily_cases)
dates  <- as.Date(inc_df$date)
ndays  <- length(obs)

run_ci <- function(beta, N) {
  crate <- max(beta, 1.0); ptran <- beta / crate
  p <- list(n = N, prevalence = obs[1] / N, incub = 5, recov = recov,
            crate = crate, ptran = ptran)
  seir_run_multi_ci(p, ndays = ndays, nreps = EVAL_REPS, nthreads = NTHREADS)
}

rows <- list()
cat(sprintf("%-16s %6s  %8s %8s %8s\n", "model", "N", "seed", "RMSE", "MAE"))
for (mname in names(betas)) {
  beta <- betas[[mname]]
  for (N in N_VALUES) {
    ci <- run_ci(beta, N)
    rmse <- sqrt(mean((ci$med - obs)^2)); mae <- mean(abs(ci$med - obs))
    cat(sprintf("%-16s %6d  %8.1f %8.2f %8.2f\n", mname, N, obs[1], rmse, mae))
    rows[[length(rows)+1]] <- data.frame(date = dates, lower = ci$lower, med = ci$med,
                                          upper = ci$upper, model = mname, N = factor(N))
  }
}
df <- bind_rows(rows)
df_obs <- data.frame(date = dates, y = obs)

g <- ggplot() +
  geom_line(data = df, aes(x = date, y = med, colour = N), linewidth = 0.8) +
  geom_line(data = df_obs, aes(x = date, y = y), colour = "black", linewidth = 0.9) +
  geom_point(data = df_obs, aes(x = date, y = y), colour = "black", size = 1) +
  facet_wrap(~model, ncol = 1) +
  labs(title = "Utah COVID-19 Current (2025) -- median simulated curve across N=5000..10000",
       subtitle = "Black = observed | Coloured lines = median of 500 SEIR runs at each N (prevalence re-seeded as obs[1]/N each time)",
       x = "Date", y = "Daily incidence", colour = "N") +
  theme_bw(base_size = 12) + theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

out_png <- file.path(PLOTS_DIR, "current61_bilstm_N_sweep.png")
ggsave(out_png, g, width = 11, height = 8, dpi = 150)
cat(sprintf("\nSaved: %s\n", out_png))
