# =============================================================================
#  bernardo/compare_slurm.R  — 5-method SEIR comparison using Bernardo BiLSTM
#
#  Identical pipeline as the main compare_slurm.R but uses predictions from
#  the Bernardo model (bernardo_bilstm_predictions.csv).
#
#  STEP 1 — Submit:   Rscript bernardo/compare_slurm.R
#  STEP 2 — Collect:  Rscript bernardo/compare_slurm.R --collect
#
#  Prerequisites:
#    python bernardo/generate_bernardo_predictions.py   (run first)
#    abc_seir_summary.csv, abcsmc_seir_summary.csv,
#    nm_seir_summary.csv, de_seir_summary.csv          (from main pipeline)
# =============================================================================

PROJECT_DIR  <- path.expand(
  "~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction"
)
BERNARDO_DIR <- file.path(PROJECT_DIR, "bernardo")

SEIR_COMMON <- file.path(PROJECT_DIR, "seir_common.R")
source(SEIR_COMMON)
seir_set_libpath()

suppressPackageStartupMessages({
  library(slurmR)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(cowplot)
  library(reticulate)
})

SCRATCH <- file.path("/scratch/general/vast", Sys.getenv("USER"), "slurmR")
dir.create(path.expand(SCRATCH), recursive = TRUE, showWarnings = FALSE)

SLURM_OPTS <- list(
  account         = "vegayon-np",
  partition       = "vegayon-np",
  `cpus-per-task` = 10L,
  `mem-per-cpu`   = "4G",
  time            = "04:00:00"
)

NDAYS         <- SEIR_NDAYS
N_SAMPLE      <- 2000L
SAMPLE_SEED   <- 42L
USE_ALL_SIMS  <- FALSE
METHOD_REPS   <- 500L
NTHREADS      <- 10L
N_CURVE_PLOTS <- 6L
BILSTM_WINDOW <- "late_365d"

METHOD_ORDER  <- SEIR_METHOD_ORDER
METHOD_COLORS <- SEIR_METHOD_COLORS[names(SEIR_METHOD_COLORS) != "Oracle"]

PLOTS_DIR <- file.path(BERNARDO_DIR, "plots")
dir.create(PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Load data
# =============================================================================

actual    <- read.csv(file.path(PROJECT_DIR, SEIR_TEST_PARAMS_FILE))
preds_all <- read.csv(file.path(BERNARDO_DIR, "bernardo_bilstm_predictions.csv"))

np_module <- reticulate::import("numpy")
inc_raw   <- as.matrix(
  np_module$load(path.expand(file.path(PROJECT_DIR, SEIR_TEST_INCIDENCE_NPY)))
)
stopifnot(nrow(inc_raw) == nrow(actual), ncol(inc_raw) == NDAYS)

abc_summary    <- read.csv(file.path(PROJECT_DIR, "abc_seir_summary.csv"))
abcsmc_summary <- read.csv(file.path(PROJECT_DIR, "abcsmc_seir_summary.csv"))
nm_summary     <- read.csv(file.path(PROJECT_DIR, "nm_seir_summary.csv"))
de_summary     <- read.csv(file.path(PROJECT_DIR, "de_seir_summary.csv"))

bl_all <- preds_all |> filter(window == BILSTM_WINDOW)

timing_path <- file.path(BERNARDO_DIR, "bernardo_bilstm_timing.csv")
bl_timing   <- if (file.exists(timing_path)) read.csv(timing_path) else NULL

# =============================================================================
# Common sim set
# =============================================================================

conv_ids <- function(df) df$sim_idx[df$converged %in% TRUE & !is.na(df$pred_beta)]

common_ids <- Reduce(intersect, list(
  actual$sim_idx,
  conv_ids(abc_summary),
  conv_ids(abcsmc_summary),
  conv_ids(nm_summary),
  conv_ids(de_summary),
  bl_all$sim_idx
))

cat(sprintf("Sims with converged predictions from ALL 5 methods: %d\n",
            length(common_ids)))

if (USE_ALL_SIMS) {
  sampled_ids <- common_ids
} else {
  set.seed(SAMPLE_SEED)
  sampled_ids <- sort(sample(common_ids, min(N_SAMPLE, length(common_ids))))
}

test_params <- actual[actual$sim_idx %in% sampled_ids, ]
test_params <- test_params[order(test_params$sim_idx), ]
mat_rows    <- match(test_params$sim_idx, actual$sim_idx)
inc_matrix  <- inc_raw[mat_rows, , drop = FALSE]

cat(sprintf("Evaluating %d sims x 5 methods x %d reps\n",
            nrow(test_params), METHOD_REPS))

beta_tbl <- bind_rows(
  abc_summary    |> select(sim_idx, pred_beta, pred_R0) |> mutate(method = "ABC"),
  abcsmc_summary |> select(sim_idx, pred_beta, pred_R0) |> mutate(method = "ABC-SMC"),
  nm_summary     |> select(sim_idx, pred_beta, pred_R0) |> mutate(method = "NelderMead"),
  de_summary     |> select(sim_idx, pred_beta, pred_R0) |> mutate(method = "DE"),
  bl_all |> filter(sim_idx %in% sampled_ids) |>
    transmute(sim_idx, pred_beta = beta_pred, pred_R0 = R0_pred, method = "BiLSTM")
)

# =============================================================================
# Worker
# =============================================================================

compare_one_sim <- function(row_idx,
                            test_params, inc_matrix, beta_tbl,
                            seir_common, method_order,
                            ndays, method_reps, nthreads) {

  source(seir_common); seir_set_libpath()
  suppressPackageStartupMessages({
    library(epiworldR); library(dplyr)
  })

  s       <- test_params[row_idx, ]
  obs_inc <- as.numeric(inc_matrix[row_idx, ])
  obs_cmp <- obs_inc[-1]

  known  <- list(n = s$n, prevalence = s$prevalence, incub = s$incub)
  true_p <- list(beta = s$beta, recov = s$recov)

  beta_to_p_local <- function(beta) {
    crate <- max(as.numeric(beta), 1.0)
    ptran <- as.numeric(beta) / crate
    list(n = s$n, prevalence = s$prevalence, incub = s$incub, recov = s$recov,
         crate = crate, ptran = ptran)
  }

  sim_curves  <- data.frame(day = seq_len(ndays), sim_idx = s$sim_idx, obs = obs_inc)
  metric_rows <- list()

  for (mth in method_order) {
    row <- beta_tbl[beta_tbl$sim_idx == s$sim_idx & beta_tbl$method == mth, ]
    if (nrow(row) == 0 || is.na(row$pred_beta[1])) next

    p_run  <- beta_to_p_local(row$pred_beta[1])
    pred_p <- list(beta = row$pred_beta[1], recov = s$recov)

    q <- tryCatch(
      seir_run_multi_ci(p_run, ndays = ndays, nreps = method_reps,
                        nthreads = nthreads, seed = s$sim_idx),
      error = function(e) NULL
    )
    if (is.null(q)) next

    mr <- seir_metrics(
      sim_idx = s$sim_idx, method = mth,
      obs_cmp = obs_cmp, pred_cmp = q$med[-1],
      true_p = true_p, pred_p = pred_p, known = known,
      converged = TRUE, n_model_runs = method_reps
    )
    metric_rows[[length(metric_rows) + 1L]] <- mr

    mkey <- seir_mkey(mth)
    sim_curves[[paste0(mkey, "_lower")]] <- q$lower
    sim_curves[[paste0(mkey, "_med")]]   <- q$med
    sim_curves[[paste0(mkey, "_upper")]] <- q$upper
  }

  list(metrics = dplyr::bind_rows(metric_rows), curves = sim_curves)
}

# =============================================================================
# STEP 1 — Submit
# =============================================================================

args         <- commandArgs(trailingOnly = TRUE)
collect_only <- "--collect" %in% args

if (!collect_only) {
  n_sims <- nrow(test_params)
  cat(sprintf("Submitting %d SLURM jobs...\n", n_sims))

  job <- Slurm_lapply(
    X   = as.list(seq_len(n_sims)),
    FUN = function(i) {
      compare_one_sim(i, test_params, inc_matrix, beta_tbl,
                      SEIR_COMMON, METHOD_ORDER,
                      NDAYS, METHOD_REPS, NTHREADS)
    },
    njobs      = min(n_sims, 500),
    mc.cores   = 1,
    job_name   = "bernardo_cmp",
    plan       = "submit",
    sbatch_opt = SLURM_OPTS,
    export     = c("compare_one_sim",
                   "test_params", "inc_matrix", "beta_tbl",
                   "SEIR_COMMON", "METHOD_ORDER",
                   "NDAYS", "METHOD_REPS", "NTHREADS"),
    tmp_path   = SCRATCH
  )

  cat("\nJobs submitted. When done:\n")
  cat("  Rscript bernardo/compare_slurm.R --collect\n\n")
  quit(save = "no")
}

# =============================================================================
# STEP 2 — Collect
# =============================================================================

cat("Collecting results...\n")

job_path    <- file.path(path.expand(SCRATCH), "bernardo_cmp")
results_raw <- Slurm_collect(read_slurm_job(job_path), any. = TRUE)
results_raw <- Filter(Negate(is.null), results_raw)
cat(sprintf("Non-null results: %d / %d\n", length(results_raw), nrow(test_params)))

metrics_all <- bind_rows(lapply(results_raw, `[[`, "metrics")) |>
  mutate(method = factor(method, levels = METHOD_ORDER))

curves_df <- bind_rows(lapply(results_raw, `[[`, "curves"))

write.csv(metrics_all,
          file.path(BERNARDO_DIR, "bernardo_5method_metrics.csv"),
          row.names = FALSE)
cat("Saved: bernardo_5method_metrics.csv\n")

# =============================================================================
# Cost table
# =============================================================================

cost_tbl <- bind_rows(
  abc_summary    |> transmute(sim_idx, method = "ABC",        time_sec, cpu_sec, n_evals, n_model_runs),
  abcsmc_summary |> transmute(sim_idx, method = "ABC-SMC",    time_sec, cpu_sec, n_evals, n_model_runs),
  nm_summary     |> transmute(sim_idx, method = "NelderMead", time_sec, cpu_sec, n_evals, n_model_runs),
  de_summary     |> transmute(sim_idx, method = "DE",         time_sec, cpu_sec, n_evals, n_model_runs)
) |>
  group_by(method) |>
  summarise(mean_wall_sec   = mean(time_sec,     na.rm = TRUE),
            median_wall_sec = median(time_sec,   na.rm = TRUE),
            mean_cpu_sec    = mean(cpu_sec,      na.rm = TRUE),
            mean_evals      = mean(n_evals,      na.rm = TRUE),
            mean_model_runs = mean(n_model_runs, na.rm = TRUE),
            .groups = "drop")

bl_mean_wall <- if (!is.null(bl_timing))
  mean(bl_timing$mean_time_sec, na.rm = TRUE) else NA_real_

cost_tbl <- cost_tbl |>
  bind_rows(data.frame(method = "BiLSTM", mean_wall_sec = bl_mean_wall,
                       median_wall_sec = NA_real_, mean_cpu_sec = NA_real_,
                       mean_evals = 0, mean_model_runs = 0,
                       stringsAsFactors = FALSE)) |>
  mutate(method = factor(method, levels = METHOD_ORDER)) |>
  arrange(method)

write.csv(cost_tbl, file.path(BERNARDO_DIR, "bernardo_5method_cost.csv"),
          row.names = FALSE)

# =============================================================================
# Summary
# =============================================================================

comparison_summary <- metrics_all |>
  group_by(method) |>
  summarise(
    n_sims            = n(),
    mean_mae          = mean(mae,          na.rm = TRUE),
    median_mae        = median(mae,        na.rm = TRUE),
    mean_smape        = mean(smape,        na.rm = TRUE),
    median_smape      = median(smape,      na.rm = TRUE),
    mean_pearson_r    = mean(pearson_r,    na.rm = TRUE),
    mean_peak_day_err = mean(peak_day_err, na.rm = TRUE),
    mean_rel_peak_err = mean(rel_peak_err, na.rm = TRUE),
    mean_err_beta_pct = mean(err_beta_pct, na.rm = TRUE),
    mean_err_R0_pct   = mean(err_R0_pct,   na.rm = TRUE),
    n_dieout          = sum(dieout_pred,   na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(cost_tbl |> select(method, mean_wall_sec, mean_model_runs),
            by = "method")

cat("\n========== Bernardo 5-Method Comparison Summary ==========\n")
print(as.data.frame(comparison_summary), digits = 4, row.names = FALSE)

write.csv(comparison_summary,
          file.path(BERNARDO_DIR, "bernardo_5method_summary.csv"),
          row.names = FALSE)

# =============================================================================
# Plots
# =============================================================================

n_sims <- nrow(test_params)

# -- Plot 1: all metrics boxplots -------------------------------------------

metrics_long <- metrics_all |>
  select(sim_idx, method, mae, smape, pearson_r, peak_day_err,
         rel_peak_err, err_beta_pct, err_R0_pct) |>
  pivot_longer(-c(sim_idx, method), names_to = "metric", values_to = "value")

metric_labels <- c(
  mae          = "MAE (incidence)",
  smape        = "sMAPE (%)",
  pearson_r    = "Pearson r (shape)",
  peak_day_err = "Peak Timing Error (days)",
  rel_peak_err = "Relative Peak Error (%)",
  err_beta_pct = "beta Error (%)",
  err_R0_pct   = "R0 Error (%)"
)
metrics_long$metric <- factor(metric_labels[metrics_long$metric],
                              levels = unname(metric_labels))

p1 <- ggplot(metrics_long, aes(x = method, y = value, fill = method)) +
  geom_boxplot(alpha = 0.65, outlier.size = 0.7) +
  geom_jitter(width = 0.15, size = 1.6, alpha = 0.55, color = "black") +
  facet_wrap(~ metric, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = METHOD_COLORS) +
  labs(title    = "5-Method Calibration Comparison — Bernardo BiLSTM (SEIR)",
       subtitle = sprintf("N=%d sims | E->I incidence | days 2-365 | %d ensemble runs",
                          n_sims, METHOD_REPS),
       x = NULL, y = "Value", fill = "Method") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold"), strip.text = element_text(size = 9),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 7))

ggsave(file.path(PLOTS_DIR, "bernardo_5method_all_metrics.png"),
       p1, width = 16, height = 9, dpi = 150)
cat("\nSaved: bernardo_5method_all_metrics.png\n")

# -- Plot 2: sMAPE violin ---------------------------------------------------

p2 <- ggplot(metrics_all, aes(x = method, y = smape, fill = method)) +
  geom_violin(alpha = 0.55, scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.10, outlier.size = 0.7, alpha = 0.85) +
  geom_jitter(width = 0.07, size = 1.5, alpha = 0.5, color = "black") +
  scale_fill_manual(values = METHOD_COLORS) +
  labs(title    = "sMAPE Distribution — Bernardo BiLSTM (SEIR)",
       subtitle = sprintf("N=%d sims | E->I incidence | days 2-365 | %d ensemble runs",
                          n_sims, METHOD_REPS),
       x = NULL, y = "sMAPE (%)", fill = "Method") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "none",
        plot.title = element_text(face = "bold"))

ggsave(file.path(PLOTS_DIR, "bernardo_5method_smape.png"),
       p2, width = 8, height = 5, dpi = 150)
cat("Saved: bernardo_5method_smape.png\n")

# -- Plot 3: predicted vs true beta / R0 ------------------------------------

param_scatter <- metrics_all |>
  select(sim_idx, method, true_beta, true_R0, pred_beta, pred_R0)

p3a <- ggplot(param_scatter,
              aes(x = true_beta, y = pred_beta, color = method, shape = method)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_color_manual(values = METHOD_COLORS) +
  labs(title = "Predicted vs True beta (Bernardo BiLSTM)",
       x = "True beta", y = "Predicted beta", color = "Method", shape = "Method") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

p3b <- ggplot(param_scatter,
              aes(x = true_R0, y = pred_R0, color = method, shape = method)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_color_manual(values = METHOD_COLORS) +
  labs(title = "Predicted vs True R0 (Bernardo BiLSTM)",
       x = "True R0", y = "Predicted R0", color = "Method", shape = "Method") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

ggsave(file.path(PLOTS_DIR, "bernardo_5method_param_scatter.png"),
       plot_grid(p3a, p3b, ncol = 2), width = 13, height = 5, dpi = 150)
cat("Saved: bernardo_5method_param_scatter.png\n")

# -- Plot 4: example incidence curves ---------------------------------------

n_plot   <- min(N_CURVE_PLOTS, nrow(test_params))
plot_ids <- test_params$sim_idx[seq_len(n_plot)]

panel_lab <- metrics_all |>
  filter(sim_idx %in% plot_ids) |>
  left_join(select(test_params, sim_idx, R0), by = "sim_idx") |>
  group_by(sim_idx, R0) |>
  summarise(txt = paste(sprintf("%s=%.2f", method, pred_R0), collapse = "  "),
            .groups = "drop") |>
  mutate(panel = sprintf("sim %d | true R0=%.2f\n%s", sim_idx, R0, txt)) |>
  select(sim_idx, panel)

curves_plot <- curves_df |>
  filter(sim_idx %in% plot_ids, day > 1) |>
  left_join(panel_lab, by = "sim_idx")

key_map    <- setNames(METHOD_ORDER, seir_mkey(METHOD_ORDER))
curve_long <- curves_plot |>
  pivot_longer(cols = matches("_(lower|med|upper)$"),
               names_to = c("mkey", ".value"),
               names_pattern = "(.*)_(lower|med|upper)") |>
  mutate(method = factor(unname(key_map[mkey]), levels = METHOD_ORDER)) |>
  filter(!is.na(method))

p4 <- ggplot(curve_long, aes(x = day)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = method), alpha = 0.10, color = NA) +
  geom_line(aes(y = med, color = method), linewidth = 0.85) +
  geom_line(data = curves_plot, aes(y = obs),
            color = "black", linewidth = 0.6, alpha = 0.7) +
  facet_wrap(~ panel, scales = "free_y", ncol = 3) +
  scale_color_manual(values = METHOD_COLORS) +
  scale_fill_manual(values = METHOD_COLORS, guide = "none") +
  labs(title    = "SEIR E->I Incidence: Observed vs 5 Methods (Bernardo BiLSTM)",
       subtitle = sprintf("Black = observed | %d ensemble runs | days 2-365",
                          METHOD_REPS),
       x = "Day", y = "Daily Incidence (E -> I)", color = NULL) +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", legend.key.width = unit(1.4, "cm"),
        strip.text = element_text(size = 7.5), panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

ggsave(file.path(PLOTS_DIR, "bernardo_5method_incidence_curves.png"),
       p4, width = 16, height = 12, dpi = 150)
cat("Saved: bernardo_5method_incidence_curves.png\n")

# -- Plot 5: accuracy vs cost -----------------------------------------------

acc_cost <- comparison_summary |>
  select(method, mean_smape, mean_wall_sec) |>
  filter(!is.na(mean_wall_sec))

if (nrow(acc_cost) > 0) {
  p5 <- ggplot(acc_cost, aes(x = mean_wall_sec, y = mean_smape,
                             color = method, label = method)) +
    geom_point(size = 4) +
    geom_text(vjust = -1.1, size = 3.5, show.legend = FALSE) +
    scale_x_log10() +
    scale_color_manual(values = METHOD_COLORS) +
    labs(title    = "Accuracy vs Cost — Bernardo BiLSTM",
         subtitle = if (is.na(bl_mean_wall))
           "BiLSTM: amortised inference (run generate_bernardo_predictions.py for timing)"
         else
           sprintf("BiLSTM inference: %.1f ms/sim (amortised)", bl_mean_wall * 1000),
         x = "Mean wall time per sim (s, log scale)",
         y = "Mean sMAPE (%)", color = "Method") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

  ggsave(file.path(PLOTS_DIR, "bernardo_5method_accuracy_vs_cost.png"),
         p5, width = 8, height = 5.5, dpi = 150)
  cat("Saved: bernardo_5method_accuracy_vs_cost.png\n")
}

# -- Plot 6: mean metric bars -----------------------------------------------

mean_metrics <- comparison_summary |>
  select(method,
         "MAE"                = mean_mae,
         "sMAPE (%)"          = mean_smape,
         "1 - Pearson r"      = mean_pearson_r,
         "Peak Timing (days)" = mean_peak_day_err,
         "Rel. Peak Err (%)"  = mean_rel_peak_err,
         "beta Err (%)"       = mean_err_beta_pct,
         "R0 Err (%)"         = mean_err_R0_pct) |>
  mutate(`1 - Pearson r` = 1 - `1 - Pearson r`) |>
  pivot_longer(-method, names_to = "metric", values_to = "value")

p6 <- ggplot(mean_metrics, aes(x = method, y = value, fill = method)) +
  geom_col(alpha = 0.75, width = 0.6) +
  facet_wrap(~ metric, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = METHOD_COLORS) +
  labs(title    = "Mean Metric Comparison — Bernardo BiLSTM (SEIR)",
       subtitle = "Lower is better for all metrics except Pearson r",
       x = NULL, y = "Mean value", fill = "Method") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 7),
        strip.text = element_text(size = 9))

ggsave(file.path(PLOTS_DIR, "bernardo_5method_mean_bar.png"),
       p6, width = 16, height = 8, dpi = 150)
cat("Saved: bernardo_5method_mean_bar.png\n")

# -- Plot 7: timing bar -----------------------------------------------------

timing_long <- cost_tbl |>
  filter(!is.na(mean_wall_sec)) |>
  mutate(label = ifelse(method == "BiLSTM",
                        sprintf("%.3f s\n(amortised)", mean_wall_sec),
                        sprintf("%.1f s", mean_wall_sec)))

p7 <- ggplot(timing_long, aes(x = method, y = mean_wall_sec, fill = method)) +
  geom_col(alpha = 0.80, width = 0.6) +
  geom_text(aes(label = label), vjust = -0.4, size = 3.2) +
  scale_y_log10() +
  scale_fill_manual(values = METHOD_COLORS) +
  labs(title    = "Computation Time per Calibrated Simulation — Bernardo BiLSTM",
       subtitle = "Log scale | BiLSTM = forward pass only (training is amortised)",
       x = NULL, y = "Mean wall time (s, log scale)", fill = "Method") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "none",
        plot.title = element_text(face = "bold"))

ggsave(file.path(PLOTS_DIR, "bernardo_5method_timing.png"),
       p7, width = 7, height = 5, dpi = 150)
cat("Saved: bernardo_5method_timing.png\n")

cat(sprintf("\nAll Bernardo outputs saved under:\n  %s\n", BERNARDO_DIR))
