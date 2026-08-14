# =============================================================================
#  validate_methods.R  — multi-parameter recovery validation
#
#  Run AFTER all calibration scripts have collected.
#
#  Compares every method on FOUR quantities:
#     crate  — calibrated
#     ptran  — calibrated
#     beta   — DERIVED as crate * ptran
#     R0     — DERIVED as beta / recov  (or calibrated if recov is in SEIR_PARS)
#
#  THE CENTRAL QUESTION this script is built to answer:
#  mean SEIR dynamics depend on crate and ptran only through their product, so
#  the data constrain beta far better than either factor. Every method should
#  therefore recover beta well and crate/ptran poorly. The script quantifies
#  that gap (the "ridge ratio") and checks whether each method's uncertainty
#  reporting is honest about it.
#
#  Point estimators (NelderMead, DE, BiLSTM) return a confident number for
#  quantities the data cannot determine. ABC and ABC-SMC return a posterior
#  that says so. That contrast is the result.
# =============================================================================

PROJECT_DIR <- path.expand(
  "~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction"
)

source(file.path(PROJECT_DIR, "seir_common.R"))
seir_set_libpath()

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(cowplot)
})

# Method that the paired tests compare everything against. Set to "BiLSTM"
# once BiLSTM predictions exist for the current test set; until then ABC-SMC
# is the natural reference among the mechanistic methods.
REFERENCE_METHOD <- "ABC-SMC"
BILSTM_WINDOW    <- "late_365d"
N_BOOT           <- 2000L
QUANTITIES       <- c("crate", "ptran", "beta", "R0")
OUT_DIR          <- PROJECT_DIR

set.seed(42)

# =============================================================================
# Load
# =============================================================================

actual <- read.csv(file.path(PROJECT_DIR, "test_actual_parameters.csv"))
abc_s  <- read.csv(file.path(PROJECT_DIR, "abc_seir_summary.csv"))
smc_s  <- read.csv(file.path(PROJECT_DIR, "abcsmc_seir_summary.csv"))
nm_s   <- read.csv(file.path(PROJECT_DIR, "nm_seir_summary.csv"))
de_s   <- read.csv(file.path(PROJECT_DIR, "de_seir_summary.csv"))

bl_path <- file.path(PROJECT_DIR, "test_bilstm_predictions_tuned.csv")
bl_s    <- if (file.exists(bl_path)) {
  b <- read.csv(bl_path); b[b$window == BILSTM_WINDOW, ]
} else NULL

ok <- function(df) df[df$converged %in% TRUE, ]

cat(sprintf("Calibrated parameters: %s\n\n", paste(SEIR_PARS, collapse = ", ")))

# =============================================================================
# Tidy predictions
# =============================================================================

truth <- actual |>
  transmute(sim_idx,
            true_crate = crate, true_ptran = ptran, true_recov = recov,
            true_beta  = crate * ptran,
            true_R0    = crate * ptran / recov)

grab <- function(df, label) {
  ok(df) |>
    transmute(sim_idx, method = label,
              pred_crate, pred_ptran, pred_beta, pred_R0)
}

pred_list <- list(grab(abc_s, "ABC"), grab(smc_s, "ABC-SMC"),
                  grab(nm_s, "NelderMead"), grab(de_s, "DE"))

# BiLSTM: its heads predict beta and R0. crate/ptran are only available if the
# network was retrained with those as targets; absent, they stay NA and BiLSTM
# simply does not appear in the crate/ptran panels.
if (!is.null(bl_s)) {
  bl <- bl_s |>
    transmute(sim_idx, method = "BiLSTM",
              pred_crate = if ("crate_pred" %in% names(bl_s)) crate_pred else NA_real_,
              pred_ptran = if ("ptran_pred" %in% names(bl_s)) ptran_pred else NA_real_,
              pred_beta  = beta_pred,
              pred_R0    = R0_pred)
  pred_list <- c(pred_list, list(bl))
  if (all(is.na(bl$pred_crate))) {
    cat("NOTE: the BiLSTM predictions file has no crate/ptran columns.\n")
    cat("      BiLSTM is compared on beta and R0 only. To include it in the\n")
    cat("      crate/ptran panels, retrain with target_cols =\n")
    cat("      ['crate','ptran','R0'] and re-export predictions.\n\n")
  }
}

preds <- bind_rows(pred_list)

n_methods  <- length(unique(preds$method))
common_ids <- preds |> count(sim_idx) |> filter(n == n_methods) |> pull(sim_idx)

preds <- preds |>
  filter(sim_idx %in% common_ids) |>
  left_join(truth, by = "sim_idx") |>
  mutate(method = factor(method, levels = SEIR_METHOD_ORDER))

cat(sprintf("Paired validation set: %d sims x %d methods\n\n",
            length(common_ids), n_methods))

# =============================================================================
# Recovery statistics
# =============================================================================

ccc <- function(x, y) 2 * cov(x, y) / (var(x) + var(y) + (mean(x) - mean(y))^2)

agreement <- function(true, pred) {
  keep <- is.finite(true) & is.finite(pred)
  true <- true[keep]; pred <- pred[keep]
  if (length(true) < 3) {
    return(data.frame(n = length(true), bias = NA, mae = NA, rmse = NA,
                      mape = NA, median_ape = NA, pearson_r = NA, r2 = NA,
                      ccc = NA, slope = NA))
  }
  fit <- lm(pred ~ true)
  data.frame(
    n          = length(true),
    bias       = mean(pred - true),
    mae        = mean(abs(pred - true)),
    rmse       = sqrt(mean((pred - true)^2)),
    mape       = mean(abs(pred - true) / (abs(true) + 1e-9)) * 100,
    median_ape = median(abs(pred - true) / (abs(true) + 1e-9)) * 100,
    pearson_r  = cor(true, pred),
    r2         = 1 - sum((pred - true)^2) / sum((true - mean(true))^2),
    ccc        = ccc(true, pred),
    slope      = unname(coef(fit)[2])
  )
}

recovery <- bind_rows(lapply(QUANTITIES, function(q) {
  preds |>
    group_by(method) |>
    reframe(quantity = q,
            agreement(.data[[paste0("true_", q)]], .data[[paste0("pred_", q)]]))
})) |>
  mutate(quantity = factor(quantity, levels = QUANTITIES)) |>
  arrange(quantity, method)

cat("========== Parameter recovery ==========\n")
print(as.data.frame(recovery), digits = 4, row.names = FALSE)
write.csv(recovery, file.path(OUT_DIR, "validation_parameter_recovery.csv"),
          row.names = FALSE)

# =============================================================================
# The ridge: is beta recovered better than its own factors?
# =============================================================================

ridge <- recovery |>
  select(method, quantity, median_ape) |>
  pivot_wider(names_from = quantity, values_from = median_ape) |>
  mutate(
    ridge_ratio_crate = crate / beta,
    ridge_ratio_ptran = ptran / beta
  )

cat("\n========== Identifiability ridge ==========\n")
print(as.data.frame(ridge), digits = 4, row.names = FALSE)
cat("\nridge_ratio = median APE of the factor divided by median APE of beta.\n")
cat("Values >> 1 mean the method recovers the PRODUCT far better than the\n")
cat("FACTORS -- the expected signature of a non-identifiable direction.\n")
cat("A ratio near 1 for a point estimator is more suspicious than reassuring:\n")
cat("it usually means the objective happened to break the tie arbitrarily.\n")

write.csv(ridge, file.path(OUT_DIR, "validation_ridge_ratio.csv"), row.names = FALSE)

# =============================================================================
# Paired comparison against the reference method
# =============================================================================

paired_test <- function(q) {
  tcol <- paste0("true_", q); pcol <- paste0("pred_", q)

  err <- preds |>
    mutate(abs_err = abs(.data[[pcol]] - .data[[tcol]])) |>
    select(sim_idx, method, abs_err) |>
    pivot_wider(names_from = method, values_from = abs_err)

  if (!REFERENCE_METHOD %in% names(err)) return(NULL)
  ref <- err[[REFERENCE_METHOD]]
  if (all(is.na(ref))) return(NULL)

  others <- setdiff(names(err), c("sim_idx", REFERENCE_METHOD))

  do.call(rbind, lapply(others, function(m) {
    dd <- err[[m]] - ref
    dd <- dd[is.finite(dd)]
    if (length(dd) < 5) return(NULL)
    bs <- replicate(N_BOOT, mean(sample(dd, length(dd), replace = TRUE)))
    data.frame(
      quantity = q, baseline = m, reference = REFERENCE_METHOD,
      n_pairs = length(dd), mean_diff = mean(dd),
      boot_lo95 = unname(quantile(bs, 0.025)),
      boot_hi95 = unname(quantile(bs, 0.975)),
      pct_ref_better = 100 * mean(dd > 0),
      wilcoxon_p = suppressWarnings(wilcox.test(dd, mu = 0)$p.value),
      stringsAsFactors = FALSE)
  }))
}

paired <- bind_rows(lapply(QUANTITIES, paired_test))

if (!is.null(paired) && nrow(paired) > 0) {
  cat(sprintf("\n========== Paired comparison vs %s ==========\n", REFERENCE_METHOD))
  print(as.data.frame(paired), digits = 4, row.names = FALSE)
  cat(sprintf("\nmean_diff > 0 means the baseline errs MORE than %s.\n",
              REFERENCE_METHOD))
  cat("A bootstrap CI spanning 0 means the two are indistinguishable, which\n")
  cat("for a validation study is usually the target claim.\n")
  write.csv(paired, file.path(OUT_DIR, "validation_paired_comparison.csv"),
            row.names = FALSE)
}

# =============================================================================
# Credible interval coverage (Bayesian methods only)
# =============================================================================

cov_one <- function(df, label) {
  d <- ok(df)
  d <- d[d$sim_idx %in% common_ids, ]
  rows <- lapply(QUANTITIES, function(q) {
    lo <- paste0(q, "_lo_95"); hi <- paste0(q, "_hi_95")
    cv <- paste0(q, "_covered")
    if (!all(c(lo, hi) %in% names(d))) return(NULL)
    data.frame(method = label, quantity = q,
               n_sims = nrow(d),
               coverage_pct = 100 * mean(d[[cv]], na.rm = TRUE),
               mean_ci_width = mean(d[[hi]] - d[[lo]], na.rm = TRUE),
               mean_rel_width = mean((d[[hi]] - d[[lo]]) /
                                       (abs(d[[paste0("true_", q)]]) + 1e-9),
                                     na.rm = TRUE),
               stringsAsFactors = FALSE)
  })
  bind_rows(rows)
}

cov_tbl <- bind_rows(cov_one(abc_s, "ABC"), cov_one(smc_s, "ABC-SMC")) |>
  mutate(quantity = factor(quantity, levels = QUANTITIES)) |>
  arrange(quantity, method)

if (nrow(cov_tbl) > 0) {
  cat("\n========== 95% credible interval coverage ==========\n")
  print(as.data.frame(cov_tbl), digits = 4, row.names = FALSE)
  cat("\nNominal is 95%. The informative pattern here is WIDE intervals for\n")
  cat("crate and ptran alongside a NARROW interval for beta, all with good\n")
  cat("coverage. That is a correctly calibrated posterior reporting an\n")
  cat("identifiability limit. NelderMead, DE and BiLSTM cannot appear in\n")
  cat("this table at all -- that capability gap is a finding, not a caveat.\n")
  write.csv(cov_tbl, file.path(OUT_DIR, "validation_interval_coverage.csv"),
            row.names = FALSE)
}

# Posterior correlation between crate and ptran
if ("post_cor_crate_ptran" %in% names(smc_s)) {
  cor_tbl <- bind_rows(
    ok(abc_s) |> transmute(method = "ABC",     post_cor_crate_ptran,
                           cv_crate, cv_ptran, cv_beta),
    ok(smc_s) |> transmute(method = "ABC-SMC", post_cor_crate_ptran,
                           cv_crate, cv_ptran, cv_beta)
  ) |>
    group_by(method) |>
    summarise(mean_post_cor = mean(post_cor_crate_ptran, na.rm = TRUE),
              mean_cv_crate = mean(cv_crate, na.rm = TRUE),
              mean_cv_ptran = mean(cv_ptran, na.rm = TRUE),
              mean_cv_beta  = mean(cv_beta,  na.rm = TRUE), .groups = "drop")

  cat("\n========== Posterior geometry ==========\n")
  print(as.data.frame(cor_tbl), digits = 4, row.names = FALSE)
  cat("A strongly NEGATIVE cor(crate, ptran) with cv_beta << cv_crate is the\n")
  cat("ridge made visible: the posterior lies along crate*ptran = constant.\n")
  write.csv(cor_tbl, file.path(OUT_DIR, "validation_posterior_geometry.csv"),
            row.names = FALSE)
}

# =============================================================================
# Plots
# =============================================================================

long <- bind_rows(lapply(QUANTITIES, function(q) {
  preds |> transmute(sim_idx, method, quantity = q,
                     true = .data[[paste0("true_", q)]],
                     pred = .data[[paste0("pred_", q)]])
})) |>
  filter(is.finite(true), is.finite(pred)) |>
  mutate(quantity = factor(quantity, levels = QUANTITIES))

# Report exactly which method x quantity combinations have data, so a method
# that is silently absent from a panel is obvious rather than mysterious.
panel_n <- long |> count(method, quantity, name = "n") |>
  tidyr::complete(method, quantity, fill = list(n = 0L))
cat("\n========== Points available per panel ==========\n")
print(as.data.frame(tidyr::pivot_wider(panel_n, names_from = quantity,
                                       values_from = n)), row.names = FALSE)
cat("A zero means that method produced no prediction for that quantity and\n")
cat("its panel will be empty (BiLSTM has no crate/ptran output heads).\n")

p1 <- ggplot(long, aes(x = true, y = pred, color = method)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray45") +
  geom_point(size = 2.2, alpha = 0.75) +
  facet_grid(quantity ~ method, scales = "free") +
  scale_color_manual(values = SEIR_METHOD_COLORS, drop = FALSE) +
  labs(title = "Predicted vs true, all calibrated and derived quantities",
       subtitle = paste0("Tight along the diagonal for beta and loose for ",
                         "crate/ptran = the identifiability ridge. ",
                         "Empty panels = method has no output for that quantity."),
       x = "True", y = "Predicted") +
  theme_bw(base_size = 9) +
  theme(legend.position = "none", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        strip.text = element_text(size = 8))

ggsave(file.path(OUT_DIR, "validation_scatter_all.png"), p1,
       width = 14, height = 10, dpi = 150)
cat("\nSaved: validation_scatter_all.png\n")

ape <- long |> mutate(ape = pmax(100 * abs(pred - true) / (abs(true) + 1e-9), 1e-3))

# Violins are invisible at small n, so points are always drawn on top and each
# panel is annotated with its sample size.
n_lab <- ape |> count(method, quantity, name = "n") |>
  left_join(ape |> group_by(method, quantity) |>
              summarise(y = max(ape), .groups = "drop"),
            by = c("method", "quantity"))

p2 <- ggplot(ape, aes(x = method, y = ape, fill = method, color = method)) +
  geom_violin(alpha = 0.35, scale = "width", trim = TRUE, color = NA) +
  geom_boxplot(width = 0.14, outlier.shape = NA, alpha = 0.6, color = "gray25") +
  geom_jitter(width = 0.12, size = 2, alpha = 0.9) +
  geom_text(data = n_lab, aes(x = method, y = y, label = paste0("n=", n)),
            vjust = -0.6, size = 2.6, color = "gray25", inherit.aes = FALSE) +
  facet_wrap(~ quantity, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = SEIR_METHOD_COLORS, drop = FALSE) +
  scale_color_manual(values = SEIR_METHOD_COLORS, drop = FALSE) +
  scale_y_log10() +
  labs(title = "Absolute percentage error by method and quantity",
       subtitle = "Log scale, individual sims shown as points. Paired design: identical sims per method.",
       x = NULL, y = "APE (%)") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 35, hjust = 1, size = 7))

ggsave(file.path(OUT_DIR, "validation_ape_by_quantity.png"), p2,
       width = 14, height = 5, dpi = 150)
cat("Saved: validation_ape_by_quantity.png\n")

# Beta and R0 only: the quantities every method predicts, so all five appear
p2b <- ape |>
  filter(quantity %in% c("beta", "R0")) |>
  ggplot(aes(x = method, y = ape, fill = method, color = method)) +
  geom_violin(alpha = 0.35, scale = "width", trim = TRUE, color = NA) +
  geom_boxplot(width = 0.14, outlier.shape = NA, alpha = 0.6, color = "gray25") +
  geom_jitter(width = 0.12, size = 2.4, alpha = 0.9) +
  facet_wrap(~ quantity, scales = "free_y") +
  scale_fill_manual(values = SEIR_METHOD_COLORS, drop = FALSE) +
  scale_color_manual(values = SEIR_METHOD_COLORS, drop = FALSE) +
  scale_y_log10() +
  labs(title = "Beta and R0 recovery: all methods",
       subtitle = "The quantities every method predicts, including BiLSTM",
       x = NULL, y = "APE (%)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(OUT_DIR, "validation_ape_beta_R0.png"), p2b,
       width = 10, height = 5, dpi = 150)
cat("Saved: validation_ape_beta_R0.png\n")

# Estimates in the (crate, ptran) plane with the true ridge overlaid
ridge_plot_ids <- head(sort(common_ids), 9)
rp <- preds |> filter(sim_idx %in% ridge_plot_ids,
                      is.finite(pred_crate), is.finite(pred_ptran))

if (nrow(rp) > 0) {
  hyp <- rp |>
    distinct(sim_idx, true_beta) |>
    rowwise() |>
    reframe(sim_idx = sim_idx,
            crate = seq(SEIR_BOUNDS$crate[1], SEIR_BOUNDS$crate[2], length.out = 200),
            ptran = true_beta / seq(SEIR_BOUNDS$crate[1], SEIR_BOUNDS$crate[2],
                                    length.out = 200)) |>
    filter(ptran >= SEIR_BOUNDS$ptran[1], ptran <= SEIR_BOUNDS$ptran[2])

  p3 <- ggplot() +
    geom_line(data = hyp, aes(x = crate, y = ptran),
              color = "gray55", linetype = "dashed", linewidth = 0.6) +
    geom_point(data = rp, aes(x = pred_crate, y = pred_ptran, color = method),
               size = 2.2, alpha = 0.85) +
    geom_point(data = distinct(rp, sim_idx, true_crate, true_ptran),
               aes(x = true_crate, y = true_ptran), shape = 4, size = 3,
               stroke = 1.1, color = "black") +
    facet_wrap(~ sim_idx, scales = "free", ncol = 3) +
    scale_color_manual(values = SEIR_METHOD_COLORS) +
    labs(title = "Estimates in the (crate, ptran) plane",
         subtitle = paste0("Dashed = the true crate*ptran ridge, X = truth. ",
                           "Points scattered ALONG the ridge means the data ",
                           "cannot separate the factors."),
         x = "contact rate", y = "transmission rate", color = "Method") +
    theme_bw(base_size = 9) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"))

  ggsave(file.path(OUT_DIR, "validation_ridge_plane.png"), p3,
         width = 12, height = 10, dpi = 150)
  cat("Saved: validation_ridge_plane.png\n")
}

cat(sprintf("\nAll validation outputs saved under:\n  %s\n", OUT_DIR))
