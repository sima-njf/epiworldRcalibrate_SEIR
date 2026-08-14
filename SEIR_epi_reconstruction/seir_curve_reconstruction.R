# =============================================================================
#  seir_curve_reconstruction.R  — incidence reconstruction and curve plots
#
#  STEP 1 — Submit:   Rscript seir_curve_reconstruction.R
#  STEP 2 — Collect:  Rscript seir_curve_reconstruction.R --collect
#
#  For each test sim, re-simulates with each method's estimated (crate, ptran)
#  and compares the resulting ensemble against the observed ABM incidence.
#
#  ORACLE = the TRUE parameters, run ORACLE_REPS times. It is not a competing
#  method: it is the stochastic floor, the error a perfect calibrator would
#  still make because the observation is one realisation. It is drawn as a
#  reference and excluded from rankings. Every method also reports "excess"
#  metrics = its value minus the Oracle's on the same sim, which isolates the
#  part of the error actually attributable to calibration.
#
#  TWO sMAPE VARIANTS
#    smape_full  days 2..365
#    smape_info  days 2..last informative day, where the last informative day
#                is the final day whose trailing 7-day mean observed incidence
#                is at least 1 case (the same rule the BiLSTM notebook uses)
#
#  Why both: sMAPE charges the maximum 200% for any day where one curve is
#  zero and the other is not. With epidemics that burn out early on a 365-day
#  window, smape_full is dominated by WHEN THE EPIDEMIC ENDED rather than by
#  the fit during the epidemic. smape_info measures the epidemic phase. Report
#  smape_info as the headline and smape_full for completeness; if the two
#  disagree sharply, the dead tail was driving the number.
#
#  PREREQUISITES: the four *_seir_summary.csv files, each carrying pred_crate
#  and pred_ptran.
# =============================================================================

PROJECT_DIR <- path.expand(
  "~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction"
)

SEIR_COMMON <- file.path(PROJECT_DIR, "seir_common.R")
source(SEIR_COMMON)
seir_set_libpath()

suppressPackageStartupMessages({
  library(slurmR); library(dplyr); library(tidyr)
  library(ggplot2); library(scales); library(reticulate)
})

SCRATCH <- file.path("/scratch/general/vast", Sys.getenv("USER"), "slurmR")
OUT_DIR <- file.path(PROJECT_DIR, "plots", "reconstruction")
dir.create(path.expand(SCRATCH), recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR,              recursive = TRUE, showWarnings = FALSE)

SLURM_OPTS <- list(
  account         = "vegayon-np",
  partition       = "vegayon-np",
  `cpus-per-task` = 4,
  `mem-per-cpu`   = "8G",
  time            = "08:00:00"
)

NDAYS         <- SEIR_NDAYS
METHOD_REPS   <- 300L    # ensemble size per calibration method
ORACLE_REPS   <- 1000L   # ensemble size for the true-parameter reference
NTHREADS      <- 4L
N_RECON_SIMS  <- 200L    # subsample: coverage/bias curves converge well before 1000
RECON_SEED    <- 11L
N_CURVE_PLOTS <- 9L
BILSTM_WINDOW <- "late_365d"

# How to obtain (crate, ptran) for BiLSTM when it only predicts beta:
#   "median" — split beta at the median contact rate (see the note where this
#              is applied). Lets BiLSTM appear in every reconstruction plot.
#   "none"   — omit BiLSTM from the reconstruction entirely.
BILSTM_CRATE_MODE <- "median"

METHOD_LEVELS  <- c("Oracle", "ABC", "ABC-SMC", "NelderMead", "DE", "BiLSTM")
RANKED_METHODS <- setdiff(METHOD_LEVELS, "Oracle")

# =============================================================================
# Load
# =============================================================================

actual <- read.csv(file.path(PROJECT_DIR, "test_actual_parameters.csv"))

np      <- reticulate::import("numpy")
inc_raw <- as.matrix(
  np$load(path.expand(file.path(PROJECT_DIR, "test_incidence_raw.npy")))
)
stopifnot(nrow(inc_raw) == nrow(actual), ncol(inc_raw) == NDAYS)

read_ok <- function(f) {
  p <- file.path(PROJECT_DIR, f)
  if (!file.exists(p)) return(NULL)
  d <- read.csv(p)
  d[d$converged %in% TRUE & !is.na(d$pred_crate) & !is.na(d$pred_ptran), ]
}

abc_s <- read_ok("abc_seir_summary.csv")
smc_s <- read_ok("abcsmc_seir_summary.csv")
nm_s  <- read_ok("nm_seir_summary.csv")
de_s  <- read_ok("de_seir_summary.csv")

bl_p <- file.path(PROJECT_DIR, "test_bilstm_predictions_tuned.csv")
bl_s <- NULL
if (file.exists(bl_p)) {
  b <- read.csv(bl_p)
  b <- b[b$window == BILSTM_WINDOW, ]

  if (all(c("crate_pred", "ptran_pred") %in% names(b))) {
    bl_s <- b |> transmute(sim_idx, pred_crate = crate_pred, pred_ptran = ptran_pred)
    cat("BiLSTM: using its own crate/ptran output heads.\n")

  } else if ("beta_pred" %in% names(b) && BILSTM_CRATE_MODE != "none") {
    # -------------------------------------------------------------------------
    # The network predicts beta but not its factors, so it cannot be re-simulated
    # directly. Because MEAN SEIR dynamics depend on crate and ptran only through
    # their product, any pair with the correct product reproduces essentially the
    # same expected curve. We therefore place BiLSTM on the ridge at a canonical
    # contact rate:
    #
    #     crate = median crate across the test set   (a prior, not fitted)
    #     ptran = beta_pred / crate
    #
    # crate is raised where necessary to keep ptran a valid probability.
    #
    # LIMITATION, state it when reporting: the mean curve is unaffected by which
    # point on the ridge is chosen, but the NOISE STRUCTURE is not. crate is the
    # number of contact draws per step, so a wrong crate changes the variance and
    # therefore the width of BiLSTM's 95% band and its coverage. Its MAE, sMAPE,
    # bias and peak-timing numbers are comparable to the other methods; its
    # coverage and CI width carry this extra assumption.
    # -------------------------------------------------------------------------
    crate_ref <- median(actual$crate)
    ptran_hi  <- SEIR_BOUNDS$ptran[2]

    bl_s <- b |>
      transmute(sim_idx,
                pred_crate = pmax(crate_ref, beta_pred / ptran_hi),
                pred_ptran = beta_pred / pmax(crate_ref, beta_pred / ptran_hi))

    n_raised <- sum(bl_s$pred_crate > crate_ref + 1e-9)
    cat(sprintf(paste0("BiLSTM: no crate/ptran heads. Splitting beta_pred at a ",
                       "canonical crate = %.3f\n         (median of the test set); ",
                       "crate raised for %d sims to keep ptran <= %.4f.\n",
                       "         Mean curve is unaffected by the split; CI width ",
                       "and coverage carry this assumption.\n"),
                crate_ref, n_raised, ptran_hi))

  } else {
    cat("NOTE: BiLSTM predictions cannot be re-simulated; omitted.\n")
  }
}

for (nm in c("abc_s", "smc_s", "nm_s", "de_s")) {
  d <- get(nm)
  cat(sprintf("  %-6s : %s\n", nm,
              if (is.null(d)) "MISSING" else sprintf("%d converged sims", nrow(d))))
}

# Sims where every available method converged
id_sets <- Filter(Negate(is.null),
                  list(actual$sim_idx, abc_s$sim_idx, smc_s$sim_idx,
                       nm_s$sim_idx, de_s$sim_idx,
                       if (!is.null(bl_s)) bl_s$sim_idx else NULL))
common_ids <- Reduce(intersect, id_sets)

set.seed(RECON_SEED)
use_ids <- sort(sample(common_ids, min(N_RECON_SIMS, length(common_ids))))

test_params <- actual[actual$sim_idx %in% use_ids, ]
test_params <- test_params[order(test_params$sim_idx), ]
inc_matrix  <- inc_raw[match(test_params$sim_idx, actual$sim_idx), , drop = FALSE]

cat(sprintf("\nSims with all methods converged: %d | reconstructing %d\n",
            length(common_ids), nrow(test_params)))
cat(sprintf("Reps: Oracle %d, each method %d\n", ORACLE_REPS, METHOD_REPS))

args         <- commandArgs(trailingOnly = TRUE)
collect_only <- "--collect" %in% args

# =============================================================================
# Helpers shared by worker and collector
# =============================================================================

# Last day whose trailing k-day mean observed incidence is >= min_count.
# Matches the informative-window rule used in the BiLSTM notebook.
last_informative_day <- function(x, k = 7L, min_count = 1.0) {
  n <- length(x)
  if (n < k) return(n)
  cs <- c(0, cumsum(x))
  rm <- (cs[(k + 1):(n + 1)] - cs[1:(n - k + 1)]) / k
  idx <- which(rm >= min_count)
  if (length(idx) == 0) return(k)
  min(max(idx) + k - 1L, n)
}

# =============================================================================
# Worker
# =============================================================================

recon_one_sim <- function(row_idx, test_params, inc_matrix, seir_common,
                          abc_s, smc_s, nm_s, de_s, bl_s,
                          ndays, method_reps, oracle_reps, nthreads) {

  source(seir_common)
  seir_set_libpath()
  suppressPackageStartupMessages({
    library(epiworldR); library(data.table)
  })

  s       <- test_params[row_idx, ]
  sid     <- as.integer(s$sim_idx)
  obs_inc <- as.numeric(inc_matrix[row_idx, ])

  known <- list(n = s$n, prevalence = s$prevalence, incub = s$incub,
                recov = s$recov)

  # (crate, ptran) per method
  pick <- function(df) {
    if (is.null(df)) return(NULL)
    r <- df[df$sim_idx == sid, ]
    if (nrow(r) == 0) return(NULL)
    c(crate = r$pred_crate[1], ptran = r$pred_ptran[1])
  }

  methods <- list(Oracle = c(crate = s$crate, ptran = s$ptran))
  for (nm in c("ABC", "ABC-SMC", "NelderMead", "DE", "BiLSTM")) {
    src <- switch(nm, ABC = abc_s, `ABC-SMC` = smc_s,
                  NelderMead = nm_s, DE = de_s, BiLSTM = bl_s)
    v <- pick(src)
    if (!is.null(v)) methods[[nm]] <- v
  }

  out <- vector("list", length(methods))

  for (j in seq_along(methods)) {
    mth  <- names(methods)[j]
    th   <- methods[[j]]
    reps <- if (mth == "Oracle") oracle_reps else method_reps

    p <- known
    p$crate <- unname(th["crate"])
    p$ptran <- unname(th["ptran"])

    cat(sprintf("  [sim %d] %-11s crate=%.3f ptran=%.4f reps=%d\n",
                sid, mth, p$crate, p$ptran, reps))

    q <- tryCatch(
      seir_run_multi_ci(p, ndays = ndays, nreps = reps, nthreads = nthreads,
                        seed = sid),
      error = function(e) {
        cat(sprintf("  [sim %d] %s FAILED: %s\n", sid, mth, e$message))
        NULL
      })

    if (!is.null(q)) {
      out[[j]] <- data.frame(
        sim_idx = sid, method = mth,
        pred_crate = p$crate, pred_ptran = p$ptran,
        pred_beta = p$crate * p$ptran,
        n_reps = reps, day = q$date, obs_inc = obs_inc,
        lower = q$lower, med = q$med, upper = q$upper,
        stringsAsFactors = FALSE)
    }
  }

  do.call(rbind, Filter(Negate(is.null), out))
}

# =============================================================================
# STEP 1 — Submit
# =============================================================================

if (!collect_only) {
  n_sims <- nrow(test_params)
  cat(sprintf("\nSubmitting reconstruction for %d sims...\n", n_sims))

  job <- Slurm_lapply(
    X   = as.list(seq_len(n_sims)),
    FUN = function(i) {
      recon_one_sim(i, test_params, inc_matrix, SEIR_COMMON,
                    abc_s, smc_s, nm_s, de_s, bl_s,
                    NDAYS, METHOD_REPS, ORACLE_REPS, NTHREADS)
    },
    njobs      = min(n_sims, 100),
    mc.cores   = 1,
    job_name   = "seir_recon",
    plan       = "submit",
    sbatch_opt = SLURM_OPTS,
    export     = c("recon_one_sim", "test_params", "inc_matrix", "SEIR_COMMON",
                   "abc_s", "smc_s", "nm_s", "de_s", "bl_s",
                   "NDAYS", "METHOD_REPS", "ORACLE_REPS", "NTHREADS"),
    tmp_path   = SCRATCH
  )

  cat("\nJobs submitted. When finished, run:\n")
  cat("  Rscript seir_curve_reconstruction.R --collect\n\n")
  quit(save = "no")
}

# =============================================================================
# STEP 2 — Collect
# =============================================================================

cat("Collecting reconstruction results...\n")

job_path    <- file.path(path.expand(SCRATCH), "seir_recon")
results_raw <- Slurm_collect(read_slurm_job(job_path), any. = TRUE)
results_raw <- Filter(Negate(is.null), results_raw)
cat(sprintf("Non-null results: %d / %d\n", length(results_raw), nrow(test_params)))

curves_df <- do.call(rbind, results_raw) |>
  mutate(method = factor(method, levels = METHOD_LEVELS)) |>
  filter(!is.na(method))

write.csv(curves_df, file.path(OUT_DIR, "seir_reconstruction_curves.csv"),
          row.names = FALSE)
cat("Saved: seir_reconstruction_curves.csv\n")

# ---- Per-sim informative window --------------------------------------------

lid_tbl <- curves_df |>
  filter(method == "Oracle") |>
  group_by(sim_idx) |>
  summarise(lid = last_informative_day(obs_inc[order(day)]), .groups = "drop")

cat(sprintf("\nInformative window: median last day %d (range %d-%d of %d)\n",
            median(lid_tbl$lid), min(lid_tbl$lid), max(lid_tbl$lid), NDAYS))

analysis_df <- curves_df |>
  filter(day > 1) |>
  left_join(lid_tbl, by = "sim_idx") |>
  mutate(in_window = day <= lid)

# =============================================================================
# Per-sim, per-method statistics
# =============================================================================

smape_v <- function(o, p) mean(2 * abs(o - p) / (abs(o) + abs(p) + 1e-6)) * 100

sim_stats <- analysis_df |>
  group_by(sim_idx, method, pred_crate, pred_ptran, pred_beta, n_reps, lid) |>
  summarise(
    coverage       = mean(obs_inc >= lower & obs_inc <= upper),
    mean_bias      = mean(med - obs_inc),
    mae            = mean(abs(med - obs_inc)),
    rmse           = sqrt(mean((med - obs_inc)^2)),
    smape_full     = smape_v(obs_inc, med),
    smape_info     = smape_v(obs_inc[in_window], med[in_window]),
    pearson_r      = seir_cor(obs_inc, med),
    ci_width       = mean(upper - lower),
    peak_day_obs   = day[which.max(obs_inc)],
    peak_day_pred  = day[which.max(med)],
    peak_day_err   = abs(day[which.max(med)] - day[which.max(obs_inc)]),
    rel_peak_err   = abs(max(med) - max(obs_inc)) / (max(obs_inc) + 1e-6) * 100,
    .groups = "drop"
  )

# Excess over the Oracle floor
oracle_ref <- sim_stats |>
  filter(method == "Oracle") |>
  select(sim_idx, o_mae = mae, o_smape_full = smape_full,
         o_smape_info = smape_info, o_peak = peak_day_err, o_cov = coverage)

sim_stats <- sim_stats |>
  left_join(oracle_ref, by = "sim_idx") |>
  mutate(excess_mae        = mae        - o_mae,
         excess_smape_full = smape_full - o_smape_full,
         excess_smape_info = smape_info - o_smape_info,
         excess_peak_err   = peak_day_err - o_peak,
         coverage_gap      = o_cov      - coverage)

write.csv(sim_stats, file.path(OUT_DIR, "seir_reconstruction_stats.csv"),
          row.names = FALSE)
cat("Saved: seir_reconstruction_stats.csv\n")

method_summary <- sim_stats |>
  group_by(method) |>
  summarise(
    n_sims             = n(),
    n_reps             = first(n_reps),
    mean_coverage      = mean(coverage),
    mean_ci_width      = mean(ci_width),
    mean_bias          = mean(mean_bias),
    mean_mae           = mean(mae),
    median_mae         = median(mae),
    mean_excess_mae    = mean(excess_mae),
    mean_smape_full    = mean(smape_full),
    median_smape_full  = median(smape_full),
    mean_smape_info    = mean(smape_info),
    median_smape_info  = median(smape_info),
    mean_excess_smape_info = mean(excess_smape_info),
    mean_pearson_r     = mean(pearson_r, na.rm = TRUE),
    mean_peak_err      = mean(peak_day_err),
    mean_rel_peak      = mean(rel_peak_err),
    .groups = "drop"
  )

cat("\n========== Reconstruction summary ==========\n")
print(as.data.frame(method_summary), digits = 3)
cat("\nOracle = true parameters: the stochastic floor, not a competitor.\n")
cat("Rank methods by mean_excess_mae and mean_excess_smape_info (0 = as good\n")
cat("as knowing the truth). If Oracle coverage is far below 0.95, the observed\n")
cat("curves and the ensembles are not from the same process -- investigate\n")
cat("before trusting any coverage number.\n")
cat("\nsmape_full vs smape_info: a large gap means the dead tail (days after\n")
cat("the epidemic ended) is driving smape_full, not the epidemic-phase fit.\n")

if (!is.null(bl_s) && BILSTM_CRATE_MODE == "median") {
  cat("\nBiLSTM CAVEAT: it predicts beta only, so its (crate, ptran) split was\n")
  cat("assigned at the median contact rate rather than estimated. Mean-curve\n")
  cat("metrics (MAE, sMAPE, bias, peak timing) are directly comparable; CI\n")
  cat("width and coverage additionally depend on that assumed contact rate,\n")
  cat("since crate sets the number of contact draws and hence the variance.\n")
}

write.csv(method_summary,
          file.path(OUT_DIR, "seir_reconstruction_method_summary.csv"),
          row.names = FALSE)

# =============================================================================
# Plots
# =============================================================================

mcol <- SEIR_METHOD_COLORS
n_used <- length(unique(curves_df$sim_idx))

# -- 1. Example incidence curves ---------------------------------------------

plot_ids <- head(sort(unique(curves_df$sim_idx)), N_CURVE_PLOTS)

lab <- sim_stats |>
  filter(sim_idx %in% plot_ids, method == "Oracle") |>
  left_join(select(actual, sim_idx, crate, ptran, R0), by = "sim_idx") |>
  mutate(panel = sprintf("sim %d | R0=%.2f  crate=%.1f  ptran=%.3f",
                         sim_idx, R0, crate, ptran)) |>
  select(sim_idx, panel)

cp <- analysis_df |> filter(sim_idx %in% plot_ids) |> left_join(lab, by = "sim_idx")

p1 <- ggplot(cp, aes(x = day)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = method),
              alpha = 0.12, color = NA) +
  geom_line(aes(y = med, color = method), linewidth = 0.75) +
  geom_line(data = filter(cp, method == "Oracle"), aes(y = obs_inc),
            color = "black", linewidth = 0.45, alpha = 0.8) +
  facet_wrap(~ panel, scales = "free", ncol = 3) +
  scale_color_manual(values = mcol) + scale_fill_manual(values = mcol) +
  scale_y_continuous(labels = comma) +
  labs(title = "SEIR incidence reconstruction",
       subtitle = sprintf(paste0("Black = observed ABM | lines = ensemble median",
                                 " | bands = 95%% CI (Oracle %d reps, methods %d)"),
                          ORACLE_REPS, METHOD_REPS),
       x = "Day", y = "Daily incidence (E -> I)", color = NULL, fill = NULL) +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", legend.key.width = unit(1.3, "cm"),
        strip.text = element_text(size = 7.5), panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "recon_incidence_curves.png"), p1,
       width = 15, height = 11, dpi = 150)
cat("\nSaved: recon_incidence_curves.png\n")

# -- 2. sMAPE: full window vs informative window -----------------------------

sm <- sim_stats |>
  select(sim_idx, method, smape_full, smape_info) |>
  pivot_longer(c(smape_full, smape_info), names_to = "window", values_to = "smape") |>
  mutate(window = recode(window,
                         smape_full = "Days 2-365 (full)",
                         smape_info = "Days 2-last informative"))

or_sm <- sm |> filter(method == "Oracle") |>
  group_by(window) |> summarise(o = mean(smape), .groups = "drop")

p2 <- ggplot(filter(sm, method %in% RANKED_METHODS),
             aes(x = method, y = smape, fill = method)) +
  geom_hline(data = or_sm, aes(yintercept = o), linetype = "dashed",
             color = mcol[["Oracle"]], linewidth = 0.7) +
  geom_violin(alpha = 0.55, scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.size = 0.5, alpha = 0.85) +
  facet_wrap(~ window, scales = "free_y") +
  scale_fill_manual(values = mcol) +
  labs(title = "sMAPE distribution by method",
       subtitle = paste0("Dashed green = Oracle (true parameters). sMAPE charges ",
                         "200% for any day where one curve is zero and the other ",
                         "is not,\nso the full window is partly measuring when the ",
                         "epidemic ended rather than epidemic-phase fit."),
       x = NULL, y = "sMAPE (%)") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 8))

ggsave(file.path(OUT_DIR, "recon_smape.png"), p2, width = 11, height = 5.5, dpi = 150)
cat("Saved: recon_smape.png\n")

# -- 3. Per-day coverage ------------------------------------------------------

cov_day <- analysis_df |>
  group_by(method, day) |>
  summarise(coverage = 100 * mean(obs_inc >= lower & obs_inc <= upper),
            .groups = "drop")

p3 <- ggplot(cov_day, aes(x = day, y = coverage, color = method)) +
  geom_hline(yintercept = 95, linetype = "dashed", color = "gray45") +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = mcol) +
  scale_y_continuous(limits = c(0, 100), labels = function(x) paste0(x, "%")) +
  labs(title = "Per-day 95% CI coverage of the observed incidence",
       subtitle = sprintf("Across %d sims | Oracle %d reps, methods %d reps",
                          n_used, ORACLE_REPS, METHOD_REPS),
       x = "Day", y = "Coverage", color = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "recon_coverage.png"), p3, width = 12, height = 5, dpi = 150)
cat("Saved: recon_coverage.png\n")

# -- 4. Per-day bias and MAE --------------------------------------------------

day_stats <- analysis_df |>
  group_by(method, day) |>
  summarise(bias = mean(med - obs_inc), mae = mean(abs(med - obs_inc)),
            .groups = "drop")

p4 <- ggplot(day_stats, aes(x = day, color = method)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray45", linewidth = 0.4) +
  geom_line(aes(y = bias), linewidth = 0.8) +
  scale_color_manual(values = mcol) + scale_y_continuous(labels = comma) +
  labs(title = "Per-day bias: ensemble median minus observed",
       subtitle = sprintf("Mean across %d sims", n_used),
       x = "Day", y = "Bias", color = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "recon_bias_per_day.png"), p4, width = 12, height = 5, dpi = 150)

p5 <- ggplot(day_stats, aes(x = day, y = mae, color = method)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = mcol) + scale_y_continuous(labels = comma) +
  labs(title = "Per-day MAE: |ensemble median - observed|",
       subtitle = sprintf("Mean across %d sims", n_used),
       x = "Day", y = "MAE", color = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "recon_mae_per_day.png"), p5, width = 12, height = 5, dpi = 150)
cat("Saved: recon_bias_per_day.png / recon_mae_per_day.png\n")

# -- 4b. Per-day sMAPE --------------------------------------------------------
#
#  Two panels, because per-day sMAPE means different things depending on which
#  days are included:
#
#   "All sims"  -- every sim contributes on every day. Late in the series most
#      sims have BOTH curves at zero, which scores 0% and therefore DILUTES the
#      mean downward. A falling tail here does not mean the reconstruction is
#      improving; it means the epidemic is over for most sims.
#
#   "Within informative window" -- a sim contributes to day d only while its
#      observed incidence is still meaningful (trailing 7-day mean >= 1 case).
#      This measures epidemic-phase fit. Sample size shrinks with day, so days
#      backed by fewer than MIN_SIMS_PER_DAY sims are dropped as unreliable.

MIN_SIMS_PER_DAY <- 20L

smape_pt <- analysis_df |>
  mutate(sm = 2 * abs(obs_inc - med) / (abs(obs_inc) + abs(med) + 1e-6) * 100)

sd_all <- smape_pt |>
  group_by(method, day) |>
  summarise(mean_smape = mean(sm), sd_smape = sd(sm), n = n(), .groups = "drop") |>
  mutate(window = "All sims (days 2-365)")

sd_info <- smape_pt |>
  filter(in_window) |>
  group_by(method, day) |>
  summarise(mean_smape = mean(sm), sd_smape = sd(sm), n = n(), .groups = "drop") |>
  filter(n >= MIN_SIMS_PER_DAY) |>
  mutate(window = "Within informative window")

smape_day <- bind_rows(sd_all, sd_info) |>
  mutate(window = factor(window, levels = c("All sims (days 2-365)",
                                            "Within informative window")),
         lo = pmax(mean_smape - sd_smape, 0),
         hi = pmin(mean_smape + sd_smape, 200))

p5b <- ggplot(smape_day, aes(x = day, color = method, fill = method)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.10, color = NA) +
  geom_line(aes(y = mean_smape), linewidth = 0.85) +
  facet_wrap(~ window, scales = "free_x") +
  scale_color_manual(values = mcol) + scale_fill_manual(values = mcol) +
  coord_cartesian(ylim = c(0, 200)) +
  labs(title = "Per-day sMAPE: ensemble median vs observed ABM incidence",
       subtitle = sprintf(paste0("Mean +/- SD across %d test sims | Oracle %d reps, ",
                                 "methods %d reps\nRight panel drops days backed by ",
                                 "fewer than %d sims"),
                          n_used, ORACLE_REPS, METHOD_REPS, MIN_SIMS_PER_DAY),
       x = "Day", y = "sMAPE (%)", color = NULL, fill = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "recon_smape_per_day.png"), p5b,
       width = 13, height = 5.5, dpi = 150)
cat("Saved: recon_smape_per_day.png\n")

# Single-panel version matching the per-day MAE figure exactly
p5c <- ggplot(filter(smape_day, window == "All sims (days 2-365)"),
              aes(x = day, color = method, fill = method)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.10, color = NA) +
  geom_line(aes(y = mean_smape), linewidth = 0.85) +
  scale_color_manual(values = mcol) + scale_fill_manual(values = mcol) +
  coord_cartesian(ylim = c(0, 200)) +
  labs(title = "Per-day sMAPE: ensemble median vs observed ABM incidence",
       subtitle = sprintf("Mean +/- SD across %d test sims | days 2-365 | methods %d reps",
                          n_used, METHOD_REPS),
       x = "Day", y = "sMAPE (%)", color = NULL, fill = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "recon_smape_per_day_single.png"), p5c,
       width = 12, height = 5, dpi = 150)
cat("Saved: recon_smape_per_day_single.png\n")

write.csv(smape_day, file.path(OUT_DIR, "seir_reconstruction_smape_per_day.csv"),
          row.names = FALSE)

# -- 5. Summary bars ----------------------------------------------------------

bars <- method_summary |>
  select(method,
         `Mean MAE` = mean_mae,
         `Mean sMAPE (informative)` = mean_smape_info,
         `Mean 95% CI coverage` = mean_coverage,
         `Mean peak timing error (days)` = mean_peak_err) |>
  pivot_longer(-method, names_to = "metric", values_to = "value") |>
  mutate(lab = ifelse(grepl("coverage", metric), sprintf("%.2f", value),
                      sprintf("%.1f", value)),
         is_oracle = method == "Oracle")

p6 <- ggplot(bars, aes(x = method, y = value, fill = method, alpha = is_oracle)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = lab), vjust = -0.35, size = 3, alpha = 1) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = mcol) +
  scale_alpha_manual(values = c(`TRUE` = 0.45, `FALSE` = 0.85), guide = "none") +
  labs(title = "Curve reconstruction quality by calibration method",
       subtitle = "Faded bar = Oracle (true parameters): reference floor, not a competitor",
       x = NULL, y = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 35, hjust = 1, size = 8))

ggsave(file.path(OUT_DIR, "recon_summary_bars.png"), p6, width = 15, height = 5, dpi = 150)
cat("Saved: recon_summary_bars.png\n")

# -- 6. Excess MAE over the floor --------------------------------------------

p7 <- sim_stats |>
  filter(method %in% RANKED_METHODS) |>
  ggplot(aes(x = method, y = excess_mae, fill = method)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = mcol[["Oracle"]], linewidth = 0.7) +
  geom_violin(alpha = 0.55, scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.11, outlier.size = 0.5, alpha = 0.85) +
  scale_fill_manual(values = mcol) +
  labs(title = "Excess MAE above the stochastic floor",
       subtitle = "0 = reconstruction as good as knowing the true parameters",
       x = NULL, y = "MAE - Oracle MAE") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(OUT_DIR, "recon_excess_mae.png"), p7, width = 9, height = 5, dpi = 150)
cat("Saved: recon_excess_mae.png\n")

cat(sprintf("\nAll reconstruction outputs under:\n  %s\n", OUT_DIR))
