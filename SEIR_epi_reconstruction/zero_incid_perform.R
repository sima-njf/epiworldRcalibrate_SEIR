# =============================================================================
#  ZERO-TAIL DEMO
#
#  Illustrates the problem described by the advisor:
#  When the last N days of an incidence curve are all zeros,
#  passing that zero-tail window to calibrate_seir() gives
#  misleading parameter estimates.
#
#  Three panels per sim:
#    (A) Full 365-day incidence — the shaded box marks the zero-tail window
#    (B) Predicted parameters from zero-tail vs actual parameters
#    (C) SEIR curve recovered from zero-tail prediction vs actual SEIR curve
# =============================================================================

library(epiworldR)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# -- Config -------------------------------------------------------------------
MODEL_DIR    <- "~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction/model"
DATA_DIR     <- "~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction"
CONTACT_RATE <- 1
NDAYS        <- 365
TAIL_DAYS    <- 25      # how many days from the end to pass to calibrate_seir
N_DEMO       <- 6       # number of example sims to show

# -- Load model and data ------------------------------------------------------
source("~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction/calibrate_seir.R")   # or library(epiworldRcalibrate)
init_bilstm_model(MODEL_DIR)

actual    <- read.csv(path.expand(file.path(DATA_DIR, "test_actual_parameters.csv")))
preds_all <- read.csv(path.expand(file.path(DATA_DIR, "test_bilstm_predictions.csv")))

np      <- reticulate::import("numpy")
inc_raw <- as.matrix(np$load(path.expand(
  file.path(DATA_DIR, "test_incidence_raw.npy"))))

# -- Single SEIR run helper ---------------------------------------------------
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

# -- Pick sims where the epidemic has clearly stabilised by day (365-TAIL_DAYS)
# i.e. the zero-tail window is genuinely flat
set.seed(42)
candidate_ids <- intersect(unique(preds_all$sim_idx), actual$sim_idx)

# Find sims where last TAIL_DAYS are all near-zero
zero_tail_ids <- Filter(function(sid) {
  mat_row <- which(actual$sim_idx == sid)
  if (length(mat_row) == 0) return(FALSE)
  obs <- as.numeric(inc_raw[mat_row, ])
  tail_vals <- obs[(NDAYS - TAIL_DAYS + 1):NDAYS]
  mean(tail_vals) < 1   # average < 1 case per day in the tail
}, candidate_ids)

cat(sprintf("Sims with zero-tail (last %d days avg < 1): %d / %d\n",
            TAIL_DAYS, length(zero_tail_ids), length(candidate_ids)))

demo_ids <- sample(zero_tail_ids, min(N_DEMO, length(zero_tail_ids)))

# -- Run the demo for each sim ------------------------------------------------
plot_list <- list()

for (sid in demo_ids) {
  act     <- actual[actual$sim_idx == sid, ]
  mat_row <- which(actual$sim_idx == sid)
  obs_full <- as.numeric(inc_raw[mat_row, ])

  n          <- act$n[1]
  recov      <- act$recov[1]
  incub      <- act$incub[1]
  prevalence <- act$prevalence[1]

  # Window passed to model: LAST TAIL_DAYS (all zeros)
  tail_start  <- NDAYS - TAIL_DAYS + 1
  tail_window <- obs_full[tail_start:NDAYS]

  cat(sprintf("\nsim %d  |  tail window [day %d - %d]  |  ",
              sid, tail_start, NDAYS))
  cat(sprintf("sum of zeros: %d / %d  |  mean: %.3f\n",
              sum(tail_window == 0), TAIL_DAYS, mean(tail_window)))

  # Predict from zero-tail
  pred_zero <- tryCatch(
    calibrate_seir(daily_cases     = tail_window,
                   population_size = n,
                   recovery_rate   = recov,
                   incubation_days = incub),
    error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL })

  if (is.null(pred_zero)) next

  # Predict from full curve (for comparison — best-case)
  pred_full <- tryCatch(
    calibrate_seir(daily_cases     = obs_full,
                   population_size = n,
                   recovery_rate   = recov,
                   incubation_days = incub),
    error = function(e) NULL)

  cat(sprintf("  Actual     : beta=%.4f  R0=%.3f\n", act$beta[1], act$R0[1]))
  cat(sprintf("  Zero-tail  : beta=%.4f  R0=%.3f\n",
              pred_zero["beta"], pred_zero["R0"]))
  if (!is.null(pred_full))
    cat(sprintf("  Full curve : beta=%.4f  R0=%.3f\n",
                pred_full["beta"], pred_full["R0"]))

  # SEIR curves
  act_inc      <- run_one_seir(n, prevalence, act$beta[1], incub, recov, seed = sid)
  zero_inc     <- run_one_seir(n, prevalence, pred_zero["beta"], incub, recov, seed = sid)
  full_inc     <- if (!is.null(pred_full))
    run_one_seir(n, prevalence, pred_full["beta"], incub, recov, seed = sid)
  else rep(NA_real_, NDAYS)

  # ── Panel A: full observed curve with zero-tail shaded ─────────────────────
  df_obs <- data.frame(day = 1:NDAYS, incidence = obs_full)

  pA <- ggplot(df_obs, aes(x = day, y = incidence)) +
    annotate("rect",
             xmin = tail_start, xmax = NDAYS,
             ymin = -Inf, ymax = Inf,
             fill = "#FF6F00", alpha = 0.20) +
    annotate("text",
             x = tail_start + TAIL_DAYS / 2, y = max(obs_full) * 0.9,
             label = sprintf("passed to\nmodel\n(%d days)", TAIL_DAYS),
             size = 2.5, color = "#FF6F00", fontface = "bold") +
    geom_line(color = "black", linewidth = 0.7) +
    labs(title = sprintf("sim %d  |  actual: beta=%.3f, R0=%.2f",
                         sid, act$beta[1], act$R0[1]),
         x = "Day", y = "Observed incidence") +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(size = 8, face = "bold"))

  # ── Panel B: parameter comparison bar chart ────────────────────────────────
  param_df <- data.frame(
    param  = rep(c("beta", "R0"), 3),
    source = rep(c("Actual", "Zero-tail (25d)", "Full curve (365d)"), each = 2),
    value  = c(act$beta[1],          act$R0[1],
               pred_zero["beta"],    pred_zero["R0"],
               if (!is.null(pred_full)) c(pred_full["beta"], pred_full["R0"])
               else c(NA, NA))
  ) |>
    mutate(source = factor(source,
                           levels = c("Actual", "Zero-tail (25d)", "Full curve (365d)")))

  pB <- ggplot(param_df, aes(x = source, y = value, fill = source)) +
    geom_col(alpha = 0.85, width = 0.6) +
    geom_text(aes(label = round(value, 3)), vjust = -0.3, size = 2.5) +
    facet_wrap(~ param, scales = "free_y") +
    scale_fill_manual(values = c("Actual"           = "#2E7D32",
                                 "Zero-tail (25d)"  = "#D32F2F",
                                 "Full curve (365d)"= "#1565C0")) +
    labs(title = "Predicted parameters", x = NULL, y = "Value") +
    theme_bw(base_size = 9) +
    theme(legend.position = "none",
          axis.text.x     = element_text(size = 7, angle = 15, hjust = 1),
          panel.grid.minor = element_blank(),
          plot.title       = element_text(size = 8, face = "bold"))

  # ── Panel C: recovered SEIR curves ────────────────────────────────────────
  curve_df <- data.frame(
    day        = 1:NDAYS,
    Actual     = act_inc,
    Zero_tail  = zero_inc,
    Full_curve = full_inc
  ) |>
    pivot_longer(-day, names_to = "source", values_to = "incidence") |>
    mutate(source = recode(source,
                           Actual     = "Actual params",
                           Zero_tail  = "Zero-tail (25d) prediction",
                           Full_curve = "Full-curve prediction"
    )) |>
    filter(!is.na(incidence))

  # Build label string first (can't use sprintf() as a name directly)
  lbl_zero <- sprintf("Zero-tail (%dd): beta=%.3f, R0=%.2f",
                      TAIL_DAYS, pred_zero["beta"], pred_zero["R0"])
  lbl_full <- if (!is.null(pred_full))
    sprintf("Full curve: beta=%.3f, R0=%.2f",
            pred_full["beta"], pred_full["R0"])
  else "Full curve"
  lbl_act  <- sprintf("Actual: beta=%.3f, R0=%.2f", act$beta[1], act$R0[1])

  curve_df <- curve_df |>
    mutate(source = recode(source,
                           "Actual params"         = lbl_act,
                           "Zero-tail (25d) prediction" = lbl_zero,
                           "Full-curve prediction" = lbl_full
    ))

  color_map <- setNames(
    c("#2E7D32", "#D32F2F", "#1565C0"),
    c(lbl_act,   lbl_zero,  lbl_full))

  pC <- ggplot(curve_df, aes(x = day, y = incidence, color = source)) +
    annotate("rect",
             xmin = tail_start, xmax = NDAYS,
             ymin = -Inf, ymax = Inf,
             fill = "#FF6F00", alpha = 0.10) +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = color_map) +
    labs(title = "Recovered SEIR curves — parameters shown in legend",
         x = "Day", y = "Incidence", color = NULL) +
    theme_bw(base_size = 9) +
    theme(legend.position  = "bottom",
          legend.direction = "vertical",
          legend.text      = element_text(size = 7),
          legend.key.width = unit(1.2, "cm"),
          panel.grid.minor = element_blank(),
          plot.title       = element_text(size = 8, face = "bold"))

  plot_list[[as.character(sid)]] <- pA / pB / pC
}

# -- Combine all sims into one figure ----------------------------------------
# wrap_plots from patchwork: one column per sim
final_plot <- patchwork::wrap_plots(plot_list, ncol = min(3, N_DEMO))

out_path <- path.expand(file.path(DATA_DIR, "zero_tail_demo.png"))
ggsave(out_path,
       final_plot,
       width  = 5 * min(3, N_DEMO),
       height = 14,
       dpi    = 150)
cat(sprintf("\nSaved: %s\n", out_path))
print(final_plot)

# -- Console summary ----------------------------------------------------------
cat(sprintf("\n=== Zero-tail (%d days) prediction summary ===\n", TAIL_DAYS))
cat("The model receives only zeros. Typical behaviour:\n")
cat("  - beta is predicted near the lower end of its training range\n")
cat("  - R0 is similarly underestimated\n")
cat("  - The recovered SEIR curve is flat (no epidemic)\n")
cat("This illustrates WHY passing zero-tail windows is misleading:\n")
cat("  The MAE looks good (all zeros vs all zeros) but the\n")
cat("  epidemiologically meaningful part of the curve is missed.\n")
cat(sprintf("\nFix: trim incidence curves to the active period before\n"))
cat(sprintf("  passing to calibrate_seir() (e.g. stop when 7-day\n"))
cat(sprintf("  rolling mean < 1 case).\n"))
input_df <- data.frame(
  day       = tail_start:NDAYS,
  incidence = tail_window
)

p_input <- ggplot(input_df, aes(x = day, y = incidence)) +
  geom_col(fill = "#FF6F00", alpha = 0.85) +
  geom_text(aes(label = round(incidence, 1)),
            vjust = -0.3, size = 2.5) +
  labs(title    = sprintf("sim %d — input passed to calibrate_seir()", sid),
       subtitle = sprintf("Days %d – %d  |  mean = %.3f  |  sum = %.0f",
                          tail_start, NDAYS,
                          mean(tail_window), sum(tail_window)),
       x = "Day", y = "Daily incidence") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 9))

print(p_input)




# 25 zeros — simulating the zero-tail of an epidemic
zero_input <- rep(0, 25)

# Call calibrate_seir with example known parameters
result <- calibrate_seir(
  daily_cases     = zero_input,
  population_size = 7000,
  recovery_rate   = 0.15,
  incubation_days = 5
)

print(result)
