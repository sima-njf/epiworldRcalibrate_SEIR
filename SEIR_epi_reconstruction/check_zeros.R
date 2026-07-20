# ================================================================================
# ZERO-TAIL DIAGNOSTICS: why does calibrate_seir() still "recover" R0 from
# an all-zero input window? Two competing hypotheses:
#
#   H1 (mean-collapse):    with no curve signal, the model's R0 prediction
#                          regresses toward the marginal mean of R0 in the
#                          training distribution, regardless of the true R0.
#   H2 (metadata-leakage): the prediction is substantially driven by the
#                          explicit inputs (population_size, recovery_rate,
#                          incubation_days) rather than the (uninformative)
#                          incidence curve itself.
#
# This script:
#   1. Runs the zero-tail prediction for ALL candidate sims (not a demo subset),
#      deliberately oversampling LOW R0 cases per the advisor's suggestion.
#   2. Plots actual R0 vs zero-tail-predicted R0 (tests H1: near-zero slope /
#      correlation would support "prediction ignores the true value").
#   3. Regresses the zero-tail prediction on the metadata covariates
#      (recovery_rate, incubation_days, population_size) to see how much of
#      the zero-tail R0 they explain (tests H2: high R^2 here means metadata,
#      not "the mean", is driving the output).
#   4. Reports the training-set-like mean/median of actual R0 as the reference
#      point for H1.
# ================================================================================

library(epiworldR)
library(dplyr)
library(ggplot2)
library(broom)

MODEL_DIR <- "~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction/model"
DATA_DIR  <- "~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction"
TAIL_DAYS <- 25
NDAYS     <- 365

source("~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction/calibrate_seir.R")
init_bilstm_model(MODEL_DIR)

actual    <- read.csv(path.expand(file.path(DATA_DIR, "test_actual_parameters.csv")))
np        <- reticulate::import("numpy")
inc_raw   <- as.matrix(np$load(path.expand(file.path(DATA_DIR, "test_incidence_raw.npy"))))

candidate_ids <- unique(actual$sim_idx)

# -- Identify zero-tail sims (same rule as your demo) ----------------------------
zero_tail_ids <- Filter(function(sid) {
  mat_row <- which(actual$sim_idx == sid)
  if (length(mat_row) == 0) return(FALSE)
  obs <- as.numeric(inc_raw[mat_row, ])
  mean(obs[(NDAYS - TAIL_DAYS + 1):NDAYS]) < 1
}, candidate_ids)

cat(sprintf("Zero-tail sims available: %d / %d\n", length(zero_tail_ids), length(candidate_ids)))

# -- Reference point for H1: distribution of actual R0 among these sims ---------
r0_all <- actual %>% filter(sim_idx %in% zero_tail_ids) %>% pull(R0)
cat(sprintf("Actual R0 among zero-tail sims: mean=%.3f  median=%.3f  sd=%.3f  range=[%.2f, %.2f]\n",
            mean(r0_all), median(r0_all), sd(r0_all), min(r0_all), max(r0_all)))

# -- Stratified sample: oversample LOW R0 sims, per advisor's suggestion --------
# Low R0 defined as bottom tertile of the zero-tail-eligible sims' actual R0.
r0_low_cut <- quantile(r0_all, 1/3)
low_ids  <- zero_tail_ids[actual$R0[match(zero_tail_ids, actual$sim_idx)] <= r0_low_cut]
rest_ids <- setdiff(zero_tail_ids, low_ids)

set.seed(1)
N_LOW  <- min(60, length(low_ids))
N_REST <- min(60, length(rest_ids))
sample_ids <- c(sample(low_ids, N_LOW), sample(rest_ids, N_REST))
cat(sprintf("Sampling %d low-R0 sims (R0 <= %.2f) and %d others for the systematic test\n",
            N_LOW, r0_low_cut, N_REST))

# -- Run zero-tail prediction for every sampled sim ------------------------------
results <- vector("list", length(sample_ids))

for (k in seq_along(sample_ids)) {
  sid <- sample_ids[k]
  act <- actual[actual$sim_idx == sid, ]
  mat_row <- which(actual$sim_idx == sid)
  obs_full <- as.numeric(inc_raw[mat_row, ])
  tail_window <- obs_full[(NDAYS - TAIL_DAYS + 1):NDAYS]

  pred <- tryCatch(
    calibrate_seir(daily_cases     = tail_window,
                   population_size = act$n[1],
                   recovery_rate   = act$recov[1],
                   incubation_days = act$incub[1]),
    error = function(e) NULL)

  if (is.null(pred)) next

  results[[k]] <- data.frame(
    sim_idx        = sid,
    actual_R0      = act$R0[1],
    actual_beta    = act$beta[1],
    n              = act$n[1],
    recov          = act$recov[1],
    incub          = act$incub[1],
    pred_R0_zero   = pred["R0"],
    pred_beta_zero = pred["beta"],
    low_r0_group   = sid %in% low_ids
  )
}

res <- bind_rows(results)
write.csv(res, path.expand(file.path(DATA_DIR, "zero_tail_diagnostic_results.csv")), row.names = FALSE)
cat(sprintf("\nSaved per-sim results: %s  (n=%d)\n",
            file.path(DATA_DIR, "zero_tail_diagnostic_results.csv"), nrow(res)))

# ================================================================================
# TEST 1 (H1 -- mean collapse): actual R0 vs zero-tail predicted R0.
# If the model ignores the true value, points will be flat/near-horizontal
# around the marginal mean rather than following the y = x line.
# ================================================================================
cor_r0 <- cor(res$actual_R0, res$pred_R0_zero, use = "complete.obs")
cat(sprintf("\nCorrelation(actual R0, zero-tail predicted R0) = %.3f\n", cor_r0))
cat("  (Near 0  -> prediction carries almost no info about true R0: supports H1)\n")
cat("  (Near 1  -> curve isn't truly zero-information, or another leak: rules out H1 alone)\n")

p1 <- ggplot(res, aes(x = actual_R0, y = pred_R0_zero, color = low_r0_group)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = mean(r0_all), linetype = "dotted", color = "#1565C0") +
  geom_point(alpha = 0.7, size = 1.8) +
  scale_color_manual(values = c("FALSE" = "#455A64", "TRUE" = "#D32F2F"),
                     labels = c("FALSE" = "Other R0", "TRUE" = "Low R0 (bottom third)"),
                     name = NULL) +
  labs(title = "Zero-tail predicted R0 vs actual R0",
       subtitle = sprintf("Dashed = y=x (perfect recovery). Dotted blue = mean actual R0 (%.2f).\nCorrelation = %.3f",
                          mean(r0_all), cor_r0),
       x = "Actual R0", y = "Zero-tail predicted R0") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(path.expand(file.path(DATA_DIR, "zero_tail_diag_actual_vs_pred.png")),
       p1, width = 7, height = 6, dpi = 200)
cat("Saved plot: zero_tail_diag_actual_vs_pred.png\n")

# ================================================================================
# TEST 2 (H2 -- metadata leakage): regress zero-tail predicted R0 on the
# metadata covariates that ARE passed explicitly (population size, recovery
# rate, incubation days). High R^2 means these -- not "the mean" -- are
# driving the output when the curve itself is uninformative.
# ================================================================================
fit <- lm(pred_R0_zero ~ n + recov + incub, data = res)
cat("\n=== H2 test: does metadata alone predict the zero-tail output? ===\n")
print(summary(fit))
cat(sprintf("R^2 = %.3f  (high R^2 => metadata, not curve signal, drives the zero-tail prediction)\n",
            summary(fit)$r.squared))

# ================================================================================
# TEST 3: does the SAME metadata relationship also explain the ACTUAL R0?
# If metadata predicts actual_R0 almost as well as it predicts pred_R0_zero,
# that would mean the training data itself ties R0 to these covariates,
# and the model may have learned that shortcut instead of reading the curve.
# ================================================================================
fit_actual <- lm(actual_R0 ~ n + recov + incub, data = res)
cat("\n=== Does metadata ALSO explain the TRUE R0 in this sample? ===\n")
cat(sprintf("R^2 (metadata -> actual R0)    = %.3f\n", summary(fit_actual)$r.squared))
cat(sprintf("R^2 (metadata -> zero-tail R0) = %.3f\n", summary(fit)$r.squared))
cat("If these two R^2 values are similar and both high, the model may simply be\n")
cat("recovering the metadata->R0 relationship that exists in the DATA GENERATING\n")
cat("PROCESS itself, not doing anything with the incidence curve at all.\n")

# -- Low-R0-specific bias summary (what the advisor asked for) ------------------
low_res <- res %>% filter(low_r0_group)
cat(sprintf("\n=== Low-R0 subset (actual R0 <= %.2f), n=%d ===\n", r0_low_cut, nrow(low_res)))
cat(sprintf("Mean actual R0       : %.3f\n", mean(low_res$actual_R0)))
cat(sprintf("Mean zero-tail pred  : %.3f\n", mean(low_res$pred_R0_zero)))
cat(sprintf("Mean bias (pred-true): %+.3f  (positive = systematic OVERestimate for low R0)\n",
            mean(low_res$pred_R0_zero - low_res$actual_R0)))




# ================================================================================
# EXTRA VISUALIZATIONS for the zero-tail mean-collapse diagnosis.
# Run this AFTER zero_tail_diagnostics.R -- it reuses the `res` data frame
# already in your session (or reload it from the saved CSV).
# ================================================================================

library(ggplot2)
library(dplyr)

DATA_DIR <- "~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction"

# If starting a fresh session, uncomment:
# res <- read.csv(path.expand(file.path(DATA_DIR, "zero_tail_diagnostic_results.csv")))

res <- res %>% mutate(bias = pred_R0_zero - actual_R0)

# ================================================================================
# PLOT 1: Bias (predicted - actual) vs actual R0.
# This is the clearest single plot for the advisor: it isolates the error itself.
# A flat model with correlation ~0 produces a DOWNWARD-SLOPING line here even
# though the y-axis is "error" -- because pred is ~constant while actual varies,
# so bias = const - actual is mechanically anti-correlated with actual.
# ================================================================================
p_bias <- ggplot(res, aes(x = actual_R0, y = bias, color = low_r0_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(alpha = 0.7, size = 1.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.6) +
  scale_color_manual(values = c("FALSE" = "#455A64", "TRUE" = "#D32F2F"),
                     labels = c("FALSE" = "Other R0", "TRUE" = "Low R0 (bottom third)"),
                     name = NULL) +
  labs(title = "Prediction bias vs true R0",
       subtitle = "Zero line = no bias. Downward trend = systematic overestimate at low R0,\nunderestimate at high R0 -- the signature of collapse toward the mean.",
       x = "Actual R0", y = "Bias (predicted - actual)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(path.expand(file.path(DATA_DIR, "zt_plot1_bias_vs_actual.png")),
       p_bias, width = 7, height = 6, dpi = 200)

# ================================================================================
# PLOT 2: Distribution overlay -- actual R0 vs zero-tail predicted R0.
# Shows visually that predictions are squeezed into a narrow band around the
# mean, while actual R0 spans the full 1-5 range.
# ================================================================================
dist_df <- bind_rows(
  data.frame(R0 = res$actual_R0,    type = "Actual R0"),
  data.frame(R0 = res$pred_R0_zero, type = "Zero-tail predicted R0")
)

p_dist <- ggplot(dist_df, aes(x = R0, fill = type)) +
  geom_density(alpha = 0.45, color = NA) +
  geom_vline(xintercept = mean(res$actual_R0), linetype = "dashed", color = "#455A64") +
  scale_fill_manual(values = c("Actual R0" = "#1565C0",
                               "Zero-tail predicted R0" = "#D32F2F"), name = NULL) +
  labs(title = "Predicted R0 collapses into a narrow band",
       subtitle = "Actual R0 spans 1-5; zero-tail predictions cluster tightly near the sample mean",
       x = expression(R[0]), y = "Density") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(path.expand(file.path(DATA_DIR, "zt_plot2_distribution_overlay.png")),
       p_dist, width = 7, height = 5, dpi = 200)

# ================================================================================
# PLOT 3: Binned boxplot -- predicted R0 by TRUE-R0 quintile.
# If the model tracked the true value, boxes would step upward left to right.
# If it collapsed to the mean, boxes will look nearly identical across bins.
# ================================================================================
res_binned <- res %>%
  mutate(R0_bin = cut(actual_R0, breaks = quantile(actual_R0, probs = seq(0, 1, 0.2)),
                      include.lowest = TRUE,
                      labels = c("Q1 (lowest)", "Q2", "Q3", "Q4", "Q5 (highest)")))

p_box <- ggplot(res_binned, aes(x = R0_bin, y = pred_R0_zero, fill = R0_bin)) +
  geom_boxplot(alpha = 0.8, outlier.size = 1) +
  geom_hline(yintercept = mean(res$actual_R0), linetype = "dashed", color = "grey30") +
  labs(title = "Zero-tail predicted R0, grouped by TRUE R0 quintile",
       subtitle = "Flat boxes across quintiles = prediction does not track the true R0 at all",
       x = "True R0 quintile", y = "Zero-tail predicted R0") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

ggsave(path.expand(file.path(DATA_DIR, "zt_plot3_boxplot_by_quintile.png")),
       p_box, width = 7, height = 5, dpi = 200)

# ================================================================================
# PLOT 4: Beta version of Plot 1 (same story for the other calibrated parameter)
# ================================================================================
res <- res %>% mutate(bias_beta = pred_beta_zero - actual_beta)

p_beta <- ggplot(res, aes(x = actual_beta, y = pred_beta_zero, color = low_r0_group)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.7, size = 1.8) +
  scale_color_manual(values = c("FALSE" = "#455A64", "TRUE" = "#D32F2F"),
                     labels = c("FALSE" = "Other R0", "TRUE" = "Low R0 (bottom third)"),
                     name = NULL) +
  labs(title = "Zero-tail predicted beta vs actual beta",
       subtitle = "Same collapse pattern visible in the underlying beta parameter",
       x = "Actual beta", y = "Zero-tail predicted beta") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(path.expand(file.path(DATA_DIR, "zt_plot4_beta_actual_vs_pred.png")),
       p_beta, width = 7, height = 6, dpi = 200)

cat("Saved 4 plots to", DATA_DIR, ":\n")
cat("  zt_plot1_bias_vs_actual.png\n")
cat("  zt_plot2_distribution_overlay.png\n")
cat("  zt_plot3_boxplot_by_quintile.png\n")
cat("  zt_plot4_beta_actual_vs_pred.png\n")








# ================================================================================
# LOW-R0 ZERO-INPUT RECONSTRUCTION TEST
#
# Idea: pick sims with LOW actual R0. Force the input passed to calibrate_seir()
# to be ALL ZEROS (length TAIL_DAYS), regardless of what the real tail looked
# like. Get the predicted parameters. Then simulate the incidence curve TWICE:
#   (A) using the ACTUAL parameters (ground truth)
#   (B) using the PREDICTED parameters (from the all-zero input)
# and plot both curves together to see whether the zero-input prediction can
# still reconstruct anything resembling the true epidemic shape.
#
# Uses the same run_one_seir() / calibrate_seir() functions as your original
# zero_tail_demo.R script.
# ================================================================================

library(epiworldR)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# -- Config ---------------------------------------------------------------------
MODEL_DIR    <- "~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction/model"
DATA_DIR     <- "~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction"
CONTACT_RATE <- 1
NDAYS        <- 365
TAIL_DAYS    <- 25      # length of the all-zero input we force
N_LOW_DEMO   <- 6        # how many low-R0 sims to show
LOW_R0_CUT   <- 2.0      # "low R0" threshold -- adjust as needed

# -- Load model and data ---------------------------------------------------------
source("~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction/calibrate_seir.R")
init_bilstm_model(MODEL_DIR)

actual <- read.csv(path.expand(file.path(DATA_DIR, "test_actual_parameters.csv")))

# -- Single SEIR run helper (same as your original script) ----------------------
run_one_seir <- function(n, prevalence, beta, incub, recov,
                         ndays = NDAYS, seed = 1L) {
  n_int <- max(as.integer(round(n)), 10L)
  prev  <- max(min(prevalence, 1.0), 1.0 / n_int)
  model <- ModelSEIRCONN(
    name              = "sim",
    n                 = n_int,
    prevalence        = prev,
    contact_rate      = CONTACT_RATE,
    transmission_rate = max(beta / CONTACT_RATE, 1e-6),
    incubation_days   = max(incub, 1.0),
    recovery_rate     = max(recov, 1e-6)
  )
  verbose_off(model)
  run(model, ndays = ndays + 1L, seed = seed)
  inc <- plot_incidence(model, plot = FALSE)[["Infected"]]
  inc <- inc[-1]
  if (length(inc) < ndays) inc <- c(inc, rep(0L, ndays - length(inc)))
  as.numeric(inc[1:ndays])
}

# -- Pick LOW-R0 sims from the actual parameters ---------------------------------
low_r0_ids <- actual %>% filter(R0 <= LOW_R0_CUT) %>% pull(sim_idx)
cat(sprintf("Sims with actual R0 <= %.2f : %d available\n", LOW_R0_CUT, length(low_r0_ids)))

set.seed(42)
demo_ids <- sample(low_r0_ids, min(N_LOW_DEMO, length(low_r0_ids)))

# -- All-zero input, forced regardless of the sim's real tail --------------------
zero_input <- rep(0, TAIL_DAYS)

# -- Run the reconstruction test for each low-R0 sim -----------------------------
plot_list <- list()

for (sid in demo_ids) {

  act <- actual[actual$sim_idx == sid, ]
  n          <- act$n[1]
  recov      <- act$recov[1]
  incub      <- act$incub[1]
  prevalence <- act$prevalence[1]
  beta_act   <- act$beta[1]
  R0_act     <- act$R0[1]

  cat(sprintf("\nsim %d  |  actual R0 = %.3f  |  forcing ALL-ZERO input (%d days)\n",
              sid, R0_act, TAIL_DAYS))

  # Predict from the forced all-zero input, using the sim's real metadata
  pred <- tryCatch(
    calibrate_seir(daily_cases     = zero_input,
                   population_size = n,
                   recovery_rate   = recov,
                   incubation_days = incub),
    error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL })

  if (is.null(pred)) next

  cat(sprintf("  Actual   : beta=%.4f  R0=%.3f\n", beta_act, R0_act))
  cat(sprintf("  Zero-pred: beta=%.4f  R0=%.3f\n", pred["beta"], pred["R0"]))

  # -- Simulate the incidence curve TWICE -----------------------------------
  # (A) with the ACTUAL parameters
  curve_actual <- run_one_seir(n, prevalence, beta_act, incub, recov, seed = sid)
  # (B) with the PREDICTED parameters (from zero input)
  curve_pred   <- run_one_seir(n, prevalence, pred["beta"], incub, recov, seed = sid)

  lbl_act  <- sprintf("Actual: beta=%.3f, R0=%.2f", beta_act, R0_act)
  lbl_pred <- sprintf("Zero-input prediction: beta=%.3f, R0=%.2f",
                      pred["beta"], pred["R0"])

  curve_df <- data.frame(
    day    = 1:NDAYS,
    Actual = curve_actual,
    Zero_prediction = curve_pred
  ) %>%
    pivot_longer(-day, names_to = "source", values_to = "incidence") %>%
    mutate(source = recode(source,
                           Actual          = lbl_act,
                           Zero_prediction = lbl_pred))

  color_map <- setNames(c("#2E7D32", "#D32F2F"), c(lbl_act, lbl_pred))

  p <- ggplot(curve_df, aes(x = day, y = incidence, color = source)) +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = color_map) +
    labs(title = sprintf("sim %d  |  actual R0=%.2f  vs  zero-input predicted R0=%.2f",
                         sid, R0_act, pred["R0"]),
         x = "Day", y = "Incidence", color = NULL) +
    theme_bw(base_size = 9) +
    theme(legend.position  = "bottom",
          legend.direction = "vertical",
          legend.text      = element_text(size = 7),
          panel.grid.minor = element_blank(),
          plot.title       = element_text(size = 8, face = "bold"))

  plot_list[[as.character(sid)]] <- p
}

# -- Combine into one figure ------------------------------------------------------
final_plot <- patchwork::wrap_plots(plot_list, ncol = min(3, N_LOW_DEMO))

out_path <- path.expand(file.path(DATA_DIR, "low_r0_zero_input_reconstruction.png"))
ggsave(out_path, final_plot,
       width = 5 * min(3, N_LOW_DEMO), height = 4 * ceiling(N_LOW_DEMO / 3), dpi = 150)
cat(sprintf("\nSaved: %s\n", out_path))
print(final_plot)

cat("\n=== What to look for ===\n")
cat("If the RED (zero-input predicted) curve reaches a much higher peak, grows\n")
cat("faster, and dies out sooner than the GREEN (actual) curve, that visually\n")
cat("confirms the mean-collapse finding: the model reconstructs an epidemic\n")
cat("shaped like a TYPICAL (higher-R0) outbreak, not the true, slower, low-R0\n")
cat("epidemic that actually happened.\n")
