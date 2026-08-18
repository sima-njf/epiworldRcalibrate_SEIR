suppressPackageStartupMessages(library(epiworldR))
setwd("~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction")
source("seir_common.R"); seir_set_libpath()

inc_df <- read.csv("real_data/measles1861_epiestim.csv")
obs <- as.numeric(inc_df$daily_cases)
known <- list(n = 187, prevalence = obs[1]/187, incub = 10, recov = 0.125)
win_len <- length(obs); ndays <- win_len

seir_summary_stats <- function(x, win_len) {
  n_phase <- max(1L, win_len %/% 3L)
  idx1 <- seq_len(min(n_phase, win_len))
  idx2 <- seq.int(length(idx1) + 1L, min(2L * n_phase, win_len))
  idx3 <- seq.int(min(2L * n_phase, win_len) + 1L, win_len)
  c(total = sum(x), peak = max(x), peak_day = which.max(x) / win_len,
    phase1 = mean(x[idx1]),
    phase2 = if (length(idx2)) mean(x[idx2]) else mean(x[idx1]),
    phase3 = if (length(idx3)) mean(x[idx3]) else mean(x[idx1]))
}
distance_fn <- function(theta, known, obs_win, win_len, ndays) {
  p     <- seir_resolve(theta, known)
  pred  <- seir_run_single(p, ndays = ndays)[seq_len(win_len)]
  s_obs <- seir_summary_stats(obs_win, win_len)
  s_sim <- seir_summary_stats(pred,    win_len)
  wts   <- 1 / (abs(s_obs) + 1)
  sqrt(mean(wts * (s_sim - s_obs)^2))
}

cat("--- New (summary-stat) objective: noise check ---\n")
for (beta_test in c(0.625, 3.8)) {  # R0 = 5, 30.4
  ds <- replicate(15, distance_fn(beta_test, known, obs, win_len, ndays))
  cat(sprintf("beta=%.3f (R0=%.1f): dist mean=%.3f sd=%.3f range=[%.3f,%.3f]\n",
      beta_test, beta_test/known$recov, mean(ds), sd(ds), min(ds), max(ds)))
}

# Also scan a broader range of R0 to see where the objective's minimum
# actually sits now
cat("\n--- Objective vs beta (10 reps averaged each) ---\n")
for (beta_test in c(0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0)) {
  ds <- replicate(10, distance_fn(beta_test, known, obs, win_len, ndays))
  cat(sprintf("beta=%.2f R0=%5.1f  mean_dist=%.3f\n", beta_test, beta_test/known$recov, mean(ds)))
}
