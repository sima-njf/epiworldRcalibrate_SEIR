# =============================================================================
#  window_all_methods_submit.R  — windowed calibration comparison
#
#  Compare ALL 5 methods when only a partial window of the epidemic curve is
#  available for calibration.
#
#  Windows: 15, 30, 60, 90 days  x  early / mid / late  = 12 windows
#
#  For each (test sim, window):
#    ABC-MCMC, ABC-SMC, Nelder-Mead, DE — calibrate on the truncated window
#    BiLSTM — predictions already computed in test_bilstm_predictions_tuned.csv
#
#  Then for each method: run SEIR with predicted params -> E->I curve ->
#    SMAPE vs full observed curve (days 2-365).
#    Timing is recorded per method.
#
#  STEP 1 — Submit:  Rscript window_all_methods_submit.R
#  STEP 2 — Collect: Rscript window_all_methods_submit.R --collect
#
#  RUNTIME: ~12 windows x 5 methods x N_SIMS workers.
#  Reduced method budgets are used to keep wall time comparable to the
#  full-calibration scripts.
# =============================================================================

PROJECT_DIR <- path.expand(
  "~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction"
)

SEIR_COMMON <- file.path(PROJECT_DIR, "seir_common.R")
source(SEIR_COMMON)
seir_set_libpath()

library(slurmR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reticulate)

SCRATCH <- file.path("/scratch/general/vast", Sys.getenv("USER"), "slurmR")
dir.create(path.expand(SCRATCH), recursive = TRUE, showWarnings = FALSE)

OUT_DIR   <- file.path(PROJECT_DIR, "results")
PLOTS_DIR <- file.path(PROJECT_DIR, "plots", "windows")
dir.create(OUT_DIR,   recursive = TRUE, showWarnings = FALSE)
dir.create(PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)

SLURM_OPTS <- list(
  account         = "vegayon-np",
  partition       = "vegayon-np",
  `cpus-per-task` = 1,
  `mem-per-cpu`   = "8G",
  time            = "06:00:00"
)

# Window definitions: 15/30/60/90 days x early/mid/late
WIN_LENGTHS <- c(15L, 30L, 60L, 90L)
T_MAX       <- SEIR_NDAYS

make_windows <- function(lengths = WIN_LENGTHS, t_max = T_MAX) {
  wins <- list()
  for (L in lengths) {
    wins[[sprintf("early_%03dd", L)]] <- list(start = 0L,            len = L, regime = "early")
    mid_s <- max(0L, (t_max - L) %/% 2L)
    wins[[sprintf("mid_%03dd",   L)]] <- list(start = mid_s,         len = L, regime = "mid")
    late_s <- max(0L, t_max - L)
    wins[[sprintf("late_%03dd",  L)]] <- list(start = late_s,        len = L, regime = "late")
  }
  wins
}

WINDOWS <- make_windows()

# Method budgets (reduced vs full-calibration scripts for tractability)
N_SAMPLE         <- 20L      # number of test sims to compare (same set as compare_all_methods.R)
TEST_SAMPLE_SEED <- 42L
EVAL_REPS        <- 300L     # SEIR ensemble runs per method for curve evaluation
NTHREADS_EVAL    <- 4L

ABC_N_SAMPLES <- 1000L
ABC_BURNIN    <- 500L
ABC_PROP_FRAC <- 0.08
ABC_EPS_FRAC  <- 0.30

SMC_N_PARTICLES   <- 100L
SMC_N_GENERATIONS <- 4L
SMC_EPS_QUANTILE  <- 0.5
SMC_BUDGET        <- 3000L

NM_MAXIT   <- 500L
NM_RELTOL  <- 1e-7
NM_RESTART <- 1L

DE_MAXITER <- 80L
DE_TOL     <- 1e-8

NDAYS        <- SEIR_NDAYS
BILSTM_WINS  <- sprintf("%s_%03dd",
                        rep(c("early", "mid", "late"), each = length(WIN_LENGTHS)),
                        rep(WIN_LENGTHS, 3))

# =============================================================================
# Load and sample test set
# =============================================================================

actual    <- read.csv(file.path(PROJECT_DIR, SEIR_TEST_PARAMS_FILE))
preds_all <- read.csv(file.path(PROJECT_DIR, "test_bilstm_predictions_tuned.csv"))

np_module <- reticulate::import("numpy")
inc_raw   <- as.matrix(
  np_module$load(path.expand(file.path(PROJECT_DIR, SEIR_TEST_INCIDENCE_NPY)))
)
stopifnot(nrow(inc_raw) == nrow(actual), ncol(inc_raw) == NDAYS)

# Keep only windows in BILSTM_WINS
bl_sub <- preds_all |>
  filter(window %in% BILSTM_WINS) |>
  select(sim_idx, window, regime,
         win_len, beta_pred, R0_pred)

# Sims that have BiLSTM predictions for all target windows
bl_by_sim <- bl_sub |>
  group_by(sim_idx) |>
  summarise(n_wins = n_distinct(window), .groups = "drop") |>
  filter(n_wins == length(WINDOWS))

common_ids <- intersect(actual$sim_idx, bl_by_sim$sim_idx)
set.seed(TEST_SAMPLE_SEED)
sampled_ids <- sort(sample(common_ids, min(N_SAMPLE, length(common_ids))))

test_params <- actual[actual$sim_idx %in% sampled_ids, ]
test_params <- test_params[order(test_params$sim_idx), ]
mat_rows    <- match(test_params$sim_idx, actual$sim_idx)
inc_matrix  <- inc_raw[mat_rows, , drop = FALSE]

cat(sprintf("Window comparison: %d sims x %d windows x 5 methods\n",
            nrow(test_params), length(WINDOWS)))

# =============================================================================
# Worker function
# =============================================================================

window_one_sim <- function(row_idx,
                           test_params, inc_matrix, bl_sub_df,
                           seir_common,
                           windows_list,
                           ndays, eval_reps, nthreads_eval,
                           abc_n, abc_burnin, abc_prop, abc_eps,
                           smc_np, smc_ng, smc_eq, smc_bud,
                           nm_maxit, nm_reltol, nm_restart,
                           de_maxiter, de_tol) {

  source(seir_common); seir_set_libpath()
  suppressPackageStartupMessages({
    library(epiworldR); library(DEoptimR); library(dplyr)
  })

  s       <- test_params[row_idx, ]
  sid     <- s$sim_idx
  obs_inc <- as.numeric(inc_matrix[row_idx, ])   # full 365-day observed curve

  known  <- list(n = s$n, prevalence = s$prevalence,
                 incub = s$incub, recov = s$recov)
  true_p <- list(beta = s$beta, recov = s$recov)
  bnd    <- seir_bounds_mat()
  d      <- length(SEIR_PARS)

  # Helper: (beta, R0) from predicted params
  beta_to_p_local <- function(beta) {
    crate <- max(as.numeric(beta), 1.0)
    ptran <- as.numeric(beta) / crate
    c(known, list(crate = crate, ptran = ptran))
  }

  # Run SEIR METHOD_REPS times, return median E->I curve
  eval_p <- function(p, seed = sid) {
    seir_run_multi_ci(p, ndays = ndays, nreps = eval_reps,
                      nthreads = nthreads_eval, seed = as.integer(seed))
  }

  smape_fn <- function(obs, pred) {
    mean(2 * abs(obs - pred) / (abs(obs) + abs(pred) + 1e-6)) * 100
  }

  # ---------------------------------------------------------------------------
  # ABC-SMC helpers (avoids mvtnorm dependency)
  # ---------------------------------------------------------------------------
  rmvn_fn <- function(mu, chol_S) as.numeric(mu + crossprod(chol_S, rnorm(length(mu))))
  dmvn_fn <- function(x, mu, S_inv, det_S) {
    dv <- as.numeric(x - mu)
    as.numeric(exp(-0.5 * dv %*% S_inv %*% dv) / sqrt((2 * pi)^d * det_S))
  }
  wq_fn <- function(x, w, probs) {
    o  <- order(x); xs <- x[o]; ws <- w[o] / sum(w); cw <- cumsum(ws)
    vapply(probs, function(p) xs[which(cw >= p)[1]], numeric(1))
  }

  out_rows <- list()

  for (win_tag in names(windows_list)) {
    w      <- windows_list[[win_tag]]
    t0_idx <- w$start + 1L           # R is 1-indexed
    t1_idx <- t0_idx + w$len - 1L
    obs_win <- obs_inc[t0_idx:t1_idx]  # calibration window
    obs_cmp <- obs_inc[-1]             # evaluation: full curve days 2-365

    # Distance on the window only
    n_evals_local <- 0L
    distance_win <- function(theta) {
      n_evals_local <<- n_evals_local + 1L
      p    <- seir_resolve(theta, known)
      pred <- seir_run_single(p, ndays = ndays)
      win_pred <- pred[t0_idx:t1_idx]
      sqrt(mean((win_pred - obs_win)^2))
    }

    set.seed(sid + as.integer(chartr("abc", "123", substr(win_tag, 1, 1))) * 1000L)

    # -- BiLSTM (from pre-computed CSV) ---------------------------------------
    bl_row <- bl_sub_df[bl_sub_df$sim_idx == sid & bl_sub_df$window == win_tag, ]
    if (nrow(bl_row) > 0) {
      bl_beta <- bl_row$beta_pred[1]
      bl_R0   <- bl_row$R0_pred[1]
      p_bl    <- beta_to_p_local(bl_beta)
      q_bl    <- eval_p(p_bl)
      pred_bl <- q_bl$med[-1]
      out_rows[[length(out_rows) + 1L]] <- data.frame(
        sim_idx    = sid,  window = win_tag,  regime = w$regime,  win_len = w$len,
        method     = "BiLSTM",
        pred_beta  = bl_beta,  pred_R0 = bl_R0,
        true_beta  = s$beta,   true_R0 = s$R0,
        smape      = smape_fn(obs_cmp, pred_bl),
        mae        = mean(abs(obs_cmp - pred_bl)),
        time_sec   = NA_real_,
        n_model_runs = eval_reps,
        stringsAsFactors = FALSE
      )
    }

    # -- ABC-MCMC on obs_win --------------------------------------------------
    n_evals_local <- 0L
    t_abc <- proc.time()
    abc_fit <- tryCatch({
      rms_obs  <- sqrt(mean(obs_win^2))
      epsilon  <- max(abc_eps * rms_obs, 1e-3)
      prop_sd  <- abc_prop * (bnd$hi - bnd$lo)
      kern     <- function(dd) exp(-(dd^2) / (2 * epsilon^2))

      cur <- (bnd$lo + bnd$hi) / 2
      if (!seir_feasible(cur, known)) {
        for (ki in 1:500) {
          cand <- bnd$lo + runif(d) * (bnd$hi - bnd$lo)
          if (seir_feasible(cand, known)) { cur <- cand; break }
        }
      }
      cur_d <- distance_win(cur); cur_k <- kern(cur_d)
      chain <- matrix(NA_real_, abc_n, d, dimnames = list(NULL, SEIR_PARS))
      n_acc <- 0L
      for (it in seq_len(abc_n)) {
        prop <- cur + rnorm(d, sd = prop_sd)
        if (seir_feasible(prop, known)) {
          pd <- distance_win(prop); pk <- kern(pd)
          acc <- if (cur_k <= 0 && pk <= 0) pd < cur_d else
            runif(1) < min(1, pk / max(cur_k, .Machine$double.xmin))
          if (acc) { cur <- prop; cur_d <- pd; cur_k <- pk; n_acc <- n_acc + 1L }
        }
        chain[it, ] <- cur
      }
      post <- chain[(abc_burnin + 1L):abc_n, , drop = FALSE]
      pred_theta <- setNames(apply(post, 2, median), SEIR_PARS)
      list(ok = TRUE, theta = pred_theta, n_sims = n_evals_local)
    }, error = function(e) list(ok = FALSE))
    tt_abc <- proc.time() - t_abc

    if (abc_fit$ok) {
      p_abc    <- seir_resolve(abc_fit$theta, known)
      q_abc    <- eval_p(p_abc)
      pred_abc <- q_abc$med[-1]
      abc_beta <- p_abc$beta
      out_rows[[length(out_rows) + 1L]] <- data.frame(
        sim_idx    = sid,  window = win_tag,  regime = w$regime,  win_len = w$len,
        method     = "ABC",
        pred_beta  = abc_beta,  pred_R0 = abc_beta / s$recov,
        true_beta  = s$beta,    true_R0 = s$R0,
        smape      = smape_fn(obs_cmp, pred_abc),
        mae        = mean(abs(obs_cmp - pred_abc)),
        time_sec   = unname(tt_abc[["elapsed"]]),
        n_model_runs = abc_fit$n_sims + eval_reps,
        stringsAsFactors = FALSE
      )
    }

    # -- ABC-SMC on obs_win ---------------------------------------------------
    n_evals_local <- 0L
    t_smc <- proc.time()
    smc_fit <- tryCatch({
      theta <- matrix(NA_real_, smc_np, d, dimnames = list(NULL, SEIR_PARS))
      filled <- 0L; tries <- 0L
      while (filled < smc_np && tries < 500L * smc_np) {
        tries <- tries + 1L
        cand  <- bnd$lo + runif(d) * (bnd$hi - bnd$lo)
        if (seir_feasible(cand, known)) { filled <- filled + 1L; theta[filled, ] <- cand }
      }
      if (filled < smc_np) stop("Could not fill prior.")
      dist <- apply(theta, 1, distance_win)
      w    <- rep(1 / smc_np, smc_np)
      eps  <- as.numeric(quantile(dist, smc_eq, names = FALSE))
      ok   <- TRUE
      for (g in 2:smc_ng) {
        mu <- colSums(w * theta); dv <- sweep(theta, 2, mu)
        S  <- 2 * (t(dv) %*% (dv * w))
        S  <- S + diag(1e-10 * pmax(diag(S), 1), d)
        chol_S <- tryCatch(chol(S), error = function(e) diag(sqrt(pmax(diag(S), 1e-12)), d))
        S_inv  <- solve(S); det_S <- det(S)
        new_theta <- matrix(NA_real_, smc_np, d, dimnames = list(NULL, SEIR_PARS))
        new_dist  <- numeric(smc_np)
        i <- 1L; atts <- 0L
        while (i <= smc_np && atts < 25L * smc_np && n_evals_local < smc_bud) {
          atts <- atts + 1L
          idx  <- sample.int(smc_np, 1L, prob = w)
          cand <- rmvn_fn(theta[idx, ], chol_S)
          if (!seir_feasible(cand, known)) next
          dd <- distance_win(cand)
          if (dd <= eps) { new_theta[i, ] <- cand; new_dist[i] <- dd; i <- i + 1L }
        }
        if (i <= smc_np) { ok <- FALSE; break }
        dens <- vapply(seq_len(smc_np), function(k)
          sum(vapply(seq_len(smc_np), function(j)
            w[j] * dmvn_fn(new_theta[k, ], theta[j, ], S_inv, det_S), numeric(1))),
          numeric(1))
        w <- 1 / pmax(dens, .Machine$double.xmin); w <- w / sum(w)
        theta <- new_theta; dist <- new_dist
        eps <- as.numeric(quantile(dist, smc_eq, names = FALSE))
      }
      pred_theta <- vapply(seq_len(ncol(theta)),
                           function(j) wq_fn(theta[, j], w, 0.5), numeric(1))
      pred_theta <- setNames(pred_theta, SEIR_PARS)
      list(ok = ok, theta = pred_theta, n_sims = n_evals_local)
    }, error = function(e) list(ok = FALSE))
    tt_smc <- proc.time() - t_smc

    if (smc_fit$ok) {
      p_smc    <- seir_resolve(smc_fit$theta, known)
      q_smc    <- eval_p(p_smc)
      pred_smc <- q_smc$med[-1]
      smc_beta <- p_smc$beta
      out_rows[[length(out_rows) + 1L]] <- data.frame(
        sim_idx    = sid,  window = win_tag,  regime = w$regime,  win_len = w$len,
        method     = "ABC-SMC",
        pred_beta  = smc_beta,  pred_R0 = smc_beta / s$recov,
        true_beta  = s$beta,    true_R0 = s$R0,
        smape      = smape_fn(obs_cmp, pred_smc),
        mae        = mean(abs(obs_cmp - pred_smc)),
        time_sec   = unname(tt_smc[["elapsed"]]),
        n_model_runs = smc_fit$n_sims + eval_reps,
        stringsAsFactors = FALSE
      )
    }

    # -- Nelder-Mead on obs_win -----------------------------------------------
    n_evals_local <- 0L
    t_nm <- proc.time()
    nm_fit <- tryCatch({
      objective <- function(theta_unc) {
        theta <- seir_from_unc(theta_unc)
        if (!seir_feasible(theta, known)) return(.Machine$double.xmax / 1e6)
        n_evals_local <<- n_evals_local + 1L
        p    <- seir_resolve(theta, known)
        pred <- seir_run_single(p, ndays = ndays)
        win_pred <- pred[t0_idx:t1_idx]
        mean((win_pred - obs_win)^2)
      }
      theta0 <- (bnd$lo + bnd$hi) / 2
      if (!seir_feasible(theta0, known)) {
        for (ki in 1:200) {
          cand <- bnd$lo + runif(d) * (bnd$hi - bnd$lo)
          if (seir_feasible(cand, known)) { theta0 <- cand; break }
        }
      }
      r <- optim(seir_to_unc(theta0), objective, method = "Nelder-Mead",
                 control = list(maxit = nm_maxit, reltol = nm_reltol))
      for (ki in seq_len(nm_restart))
        r <- optim(r$par, objective, method = "Nelder-Mead",
                   control = list(maxit = nm_maxit, reltol = nm_reltol))
      list(ok = TRUE, theta = seir_from_unc(r$par), n_sims = n_evals_local)
    }, error = function(e) list(ok = FALSE))
    tt_nm <- proc.time() - t_nm

    if (nm_fit$ok) {
      p_nm    <- seir_resolve(nm_fit$theta, known)
      q_nm    <- eval_p(p_nm)
      pred_nm <- q_nm$med[-1]
      nm_beta <- p_nm$beta
      out_rows[[length(out_rows) + 1L]] <- data.frame(
        sim_idx    = sid,  window = win_tag,  regime = w$regime,  win_len = w$len,
        method     = "NelderMead",
        pred_beta  = nm_beta,  pred_R0 = nm_beta / s$recov,
        true_beta  = s$beta,   true_R0 = s$R0,
        smape      = smape_fn(obs_cmp, pred_nm),
        mae        = mean(abs(obs_cmp - pred_nm)),
        time_sec   = unname(tt_nm[["elapsed"]]),
        n_model_runs = nm_fit$n_sims + eval_reps,
        stringsAsFactors = FALSE
      )
    }

    # -- Differential Evolution on obs_win ------------------------------------
    n_evals_local <- 0L
    t_de <- proc.time()
    de_fit <- tryCatch({
      DE_NP <- 10L * d
      objective_de <- function(theta) {
        theta <- pmin(pmax(theta, bnd$lo), bnd$hi)
        if (!seir_feasible(theta, known)) return(.Machine$double.xmax / 1e6)
        n_evals_local <<- n_evals_local + 1L
        p    <- seir_resolve(theta, known)
        pred <- seir_run_single(p, ndays = ndays)
        win_pred <- pred[t0_idx:t1_idx]
        mean((win_pred - obs_win)^2)
      }
      r <- DEoptimR::JDEoptim(lower = bnd$lo, upper = bnd$hi,
                               fn = objective_de, NP = DE_NP,
                               maxiter = de_maxiter, tol = de_tol)
      list(ok = TRUE, theta = setNames(as.numeric(r$par), SEIR_PARS),
           n_sims = n_evals_local)
    }, error = function(e) list(ok = FALSE))
    tt_de <- proc.time() - t_de

    if (de_fit$ok) {
      p_de    <- seir_resolve(de_fit$theta, known)
      q_de    <- eval_p(p_de)
      pred_de <- q_de$med[-1]
      de_beta <- p_de$beta
      out_rows[[length(out_rows) + 1L]] <- data.frame(
        sim_idx    = sid,  window = win_tag,  regime = w$regime,  win_len = w$len,
        method     = "DE",
        pred_beta  = de_beta,  pred_R0 = de_beta / s$recov,
        true_beta  = s$beta,   true_R0 = s$R0,
        smape      = smape_fn(obs_cmp, pred_de),
        mae        = mean(abs(obs_cmp - pred_de)),
        time_sec   = unname(tt_de[["elapsed"]]),
        n_model_runs = de_fit$n_sims + eval_reps,
        stringsAsFactors = FALSE
      )
    }

  } # end window loop

  if (length(out_rows) == 0) return(NULL)
  dplyr::bind_rows(out_rows)
}

# =============================================================================
# STEP 1 — Submit
# =============================================================================

args         <- commandArgs(trailingOnly = TRUE)
collect_only <- "--collect" %in% args

if (!collect_only) {
  n_sims <- nrow(test_params)
  cat(sprintf("Submitting window comparison: %d sims x %d windows x 5 methods\n",
              n_sims, length(WINDOWS)))
  cat(sprintf("  Windows: %s\n", paste(names(WINDOWS), collapse = ", ")))
  cat(sprintf("  ABC: %d samples / DE: %d iters / NM: %d maxit\n",
              ABC_N_SAMPLES, DE_MAXITER, NM_MAXIT))

  job <- Slurm_lapply(
    X   = as.list(seq_len(n_sims)),
    FUN = function(i) {
      window_one_sim(
        i, test_params, inc_matrix, bl_sub, SEIR_COMMON, WINDOWS,
        NDAYS, EVAL_REPS, NTHREADS_EVAL,
        ABC_N_SAMPLES, ABC_BURNIN, ABC_PROP_FRAC, ABC_EPS_FRAC,
        SMC_N_PARTICLES, SMC_N_GENERATIONS, SMC_EPS_QUANTILE, SMC_BUDGET,
        NM_MAXIT, NM_RELTOL, NM_RESTART,
        DE_MAXITER, DE_TOL
      )
    },
    njobs      = min(n_sims, 200),
    mc.cores   = 1,
    job_name   = "seir_window_cmp",
    plan       = "submit",
    sbatch_opt = SLURM_OPTS,
    export     = c("window_one_sim",
                   "test_params", "inc_matrix", "bl_sub", "SEIR_COMMON", "WINDOWS",
                   "NDAYS", "EVAL_REPS", "NTHREADS_EVAL",
                   "ABC_N_SAMPLES", "ABC_BURNIN", "ABC_PROP_FRAC", "ABC_EPS_FRAC",
                   "SMC_N_PARTICLES", "SMC_N_GENERATIONS", "SMC_EPS_QUANTILE", "SMC_BUDGET",
                   "NM_MAXIT", "NM_RELTOL", "NM_RESTART",
                   "DE_MAXITER", "DE_TOL"),
    tmp_path   = SCRATCH
  )

  cat("\nJobs submitted. When done:\n")
  cat("  Rscript window_all_methods_submit.R --collect\n\n")
  quit(save = "no")
}

# =============================================================================
# STEP 2 — Collect
# =============================================================================

cat("Collecting window comparison results...\n")

job_path    <- file.path(path.expand(SCRATCH), "seir_window_cmp")
results_raw <- Slurm_collect(read_slurm_job(job_path), any. = TRUE)
results_raw <- Filter(Negate(is.null), results_raw)
cat(sprintf("Non-null results: %d / %d\n", length(results_raw), nrow(test_params)))

all_df <- dplyr::bind_rows(results_raw) |>
  mutate(
    method = factor(method, levels = c("ABC", "ABC-SMC", "NelderMead", "DE", "BiLSTM")),
    regime = factor(regime, levels = c("early", "mid", "late")),
    win_len = as.integer(win_len)
  )

write.csv(all_df, file.path(OUT_DIR, "window_all_methods_results.csv"), row.names = FALSE)
cat(sprintf("Saved %d rows -> %s/window_all_methods_results.csv\n",
            nrow(all_df), OUT_DIR))

# Summary: mean SMAPE per (method, window, regime)
summary_df <- all_df |>
  group_by(method, window, regime, win_len) |>
  summarise(
    n_sims        = n(),
    mean_smape    = mean(smape,    na.rm = TRUE),
    median_smape  = median(smape,  na.rm = TRUE),
    sd_smape      = sd(smape,      na.rm = TRUE),
    mean_mae      = mean(mae,      na.rm = TRUE),
    mean_time_sec = mean(time_sec, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(win_len, regime, method)

write.csv(summary_df, file.path(OUT_DIR, "window_all_methods_summary.csv"), row.names = FALSE)

cat("\n========== Window Comparison Summary ==========\n")
print(as.data.frame(summary_df), digits = 3, row.names = FALSE)

# =============================================================================
# Plots
# =============================================================================

MC <- c(
  ABC        = "#E65100",
  `ABC-SMC`  = "#D4537E",
  NelderMead = "#6A1B9A",
  DE         = "#00695C",
  BiLSTM     = "#1565C0"
)

# -- Plot 1: sMAPE vs window length, faceted by regime -----------------------
p1 <- ggplot(summary_df,
             aes(x = win_len, y = mean_smape, color = method, group = method)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.5) +
  facet_wrap(~ regime, ncol = 3,
             labeller = labeller(regime = c(early = "Early window",
                                            mid   = "Mid window",
                                            late  = "Late window"))) +
  scale_x_continuous(breaks = WIN_LENGTHS) +
  scale_color_manual(values = MC) +
  labs(
    title    = "sMAPE vs Observation Window Length — All 5 Methods",
    subtitle = "Calibrated on first/mid/last N days | Evaluated on full 365-day E->I curve",
    x = "Observation window (days)",
    y = "Mean sMAPE (%)",
    color = "Method"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "bottom",
        plot.title       = element_text(face = "bold"))

ggsave(file.path(PLOTS_DIR, "window_smape_vs_length.png"),
       p1, width = 12, height = 5, dpi = 150)
cat("Saved: window_smape_vs_length.png\n")

# -- Plot 2: sMAPE boxplots, faceted by (regime, win_len) --------------------
p2 <- ggplot(all_df,
             aes(x = method, y = smape, fill = method)) +
  geom_boxplot(alpha = 0.65, outlier.size = 0.5) +
  facet_grid(regime ~ factor(win_len),
             labeller = labeller(
               win_len = function(x) paste0(x, "d"),
               regime  = c(early = "Early", mid = "Mid", late = "Late")
             )) +
  scale_fill_manual(values = MC) +
  labs(
    title    = "sMAPE by Window Length and Regime — 5 Methods",
    subtitle = "Calibrated on partial window | Evaluated on full 365-day E->I curve",
    x = NULL, y = "sMAPE (%)", fill = "Method"
  ) +
  theme_bw(base_size = 9) +
  theme(axis.text.x      = element_text(angle = 40, hjust = 1, size = 7),
        panel.grid.minor = element_blank(),
        legend.position  = "bottom",
        plot.title       = element_text(face = "bold"),
        strip.text       = element_text(size = 8))

ggsave(file.path(PLOTS_DIR, "window_smape_boxplot.png"),
       p2, width = 14, height = 9, dpi = 150)
cat("Saved: window_smape_boxplot.png\n")

# -- Plot 3: timing comparison (methods that have timing) --------------------
timing_sum <- all_df |>
  filter(!is.na(time_sec)) |>
  group_by(method, win_len) |>
  summarise(mean_time = mean(time_sec, na.rm = TRUE), .groups = "drop")

# Add BiLSTM timing from generate_tuned_predictions.py if available
timing_path <- file.path(PROJECT_DIR, "bilstm_timing.csv")
if (file.exists(timing_path)) {
  bl_t <- read.csv(timing_path)
  bl_t_val <- mean(bl_t$mean_time_sec, na.rm = TRUE)
  bl_rows <- data.frame(method   = "BiLSTM",
                         win_len  = WIN_LENGTHS,
                         mean_time = bl_t_val,
                         stringsAsFactors = FALSE)
  timing_sum <- bind_rows(timing_sum, bl_rows) |>
    mutate(method = factor(method, levels = c("ABC","ABC-SMC","NelderMead","DE","BiLSTM")))
}

p3 <- ggplot(timing_sum, aes(x = factor(win_len), y = mean_time,
                              fill = method, group = method)) +
  geom_col(position = "dodge", alpha = 0.80) +
  scale_y_log10() +
  scale_fill_manual(values = MC) +
  labs(
    title    = "Calibration Time vs Window Length",
    subtitle = "Log scale | BiLSTM = forward pass only (amortised training)",
    x = "Window length (days)", y = "Mean wall time (s, log scale)",
    fill = "Method"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title       = element_text(face = "bold"),
        legend.position  = "bottom")

ggsave(file.path(PLOTS_DIR, "window_timing_comparison.png"),
       p3, width = 10, height = 5, dpi = 150)
cat("Saved: window_timing_comparison.png\n")

# -- Plot 4: beta prediction error vs window length --------------------------
beta_err <- all_df |>
  mutate(err_beta_pct = abs(pred_beta - true_beta) / (true_beta + 1e-9) * 100) |>
  group_by(method, win_len, regime) |>
  summarise(mean_err = mean(err_beta_pct, na.rm = TRUE), .groups = "drop")

p4 <- ggplot(beta_err,
             aes(x = win_len, y = mean_err, color = method, group = method)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.5) +
  facet_wrap(~ regime, ncol = 3) +
  scale_x_continuous(breaks = WIN_LENGTHS) +
  scale_color_manual(values = MC) +
  labs(
    title = "Beta Prediction Error vs Window Length",
    x = "Window (days)", y = "Mean |pred beta - true beta| / true beta  (%)",
    color = "Method"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "bottom",
        plot.title       = element_text(face = "bold"))

ggsave(file.path(PLOTS_DIR, "window_beta_error.png"),
       p4, width = 12, height = 5, dpi = 150)
cat("Saved: window_beta_error.png\n")

cat(sprintf("\nAll window comparison outputs saved under:\n  %s\n", PLOTS_DIR))
