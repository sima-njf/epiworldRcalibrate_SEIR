# =============================================================================
#  window_30day_submit_tuned.R
#
#  Same as bernardo/window_30day_submit.R but uses the tuned BiLSTM model
#  (model/model_bilstm_tuned.pt) instead of the Bernardo model.
#
#  Window: obs_inc[2:31]  (days 2-31, 30 days, no seed day)
#
#  BiLSTM (Tuned): uses pre-computed "noseed_030d" predictions from
#    test_bilstm_predictions_tuned.csv
#    (run generate_tuned_predictions.py first if not up-to-date)
#
#  ABC / ABC-SMC / NM / DE: calibrate in real-time on the 30-day window.
#
#  STEP 1 — Submit:   Rscript window_30day_submit_tuned.R
#  STEP 2 — Collect:  Rscript window_30day_submit_tuned.R --collect
# =============================================================================

PROJECT_DIR <- path.expand(
  "~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction"
)
TUNED_DIR   <- file.path(PROJECT_DIR, "tuned_w30")

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
dir.create(TUNED_DIR, recursive = TRUE, showWarnings = FALSE)

SLURM_OPTS <- list(
  account         = "vegayon-np",
  partition       = "vegayon-np",
  `cpus-per-task` = 4L,
  `mem-per-cpu`   = "4G",
  time            = "04:00:00"
)

NDAYS        <- SEIR_NDAYS
N_SAMPLE     <- 2000L
SAMPLE_SEED  <- 42L
USE_ALL_SIMS <- FALSE
EVAL_REPS    <- 500L
NTHREADS     <- 4L

# Window: days 2-31 (R 1-indexed: obs_inc[2:31]), 30 days, no seed
WIN_START <- 2L
WIN_END   <- 31L
WIN_LEN   <- 30L
WIN_TAG   <- "noseed_030d"

# Reduced calibration budgets
ABC_N      <- 2000L;  ABC_BURNIN <- 1000L
ABC_PROP   <- 0.08;   ABC_EPS    <- 0.30
SMC_NP     <- 150L;   SMC_NG     <- 4L;    SMC_EQ <- 0.5;  SMC_BUD <- 5000L
NM_MAXIT   <- 800L;   NM_RELTOL  <- 1e-7;  NM_RESTART <- 2L
DE_MAXITER <- 100L;   DE_TOL     <- 1e-8

METHOD_ORDER  <- SEIR_METHOD_ORDER
METHOD_COLORS <- SEIR_METHOD_COLORS[names(SEIR_METHOD_COLORS) != "Oracle"]

PLOTS_DIR <- file.path(TUNED_DIR, "plots")
dir.create(PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Load data
# =============================================================================

actual    <- read.csv(file.path(PROJECT_DIR, SEIR_TEST_PARAMS_FILE))
preds_all <- read.csv(file.path(PROJECT_DIR, "test_bilstm_predictions_tuned.csv"))

np_module <- reticulate::import("numpy")
inc_raw   <- as.matrix(
  np_module$load(path.expand(file.path(PROJECT_DIR, SEIR_TEST_INCIDENCE_NPY)))
)
stopifnot(nrow(inc_raw) == nrow(actual), ncol(inc_raw) == NDAYS)

# Tuned BiLSTM noseed_030d predictions
bl_noseed <- preds_all |> filter(window == WIN_TAG)

if (nrow(bl_noseed) == 0)
  stop("No 'noseed_030d' predictions found. Re-run generate_tuned_predictions.py first.")

cat(sprintf("Tuned BiLSTM noseed_030d predictions: %d sims\n", nrow(bl_noseed)))

# Common sims
common_ids <- intersect(actual$sim_idx, bl_noseed$sim_idx)
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

bl_sub <- bl_noseed[bl_noseed$sim_idx %in% sampled_ids, ]

cat(sprintf("Evaluating %d sims x 5 methods on first 30 days (no seed)\n",
            nrow(test_params)))

# =============================================================================
# Worker: one sim -> 5 methods on 30-day window
# =============================================================================

w30_one_sim <- function(row_idx,
                        test_params, inc_matrix, bl_sub,
                        seir_common, method_order,
                        ndays, win_start, win_end,
                        eval_reps, nthreads,
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
  obs_inc <- as.numeric(inc_matrix[row_idx, ])

  obs_win <- obs_inc[win_start:win_end]
  obs_cmp <- obs_inc[-1]

  known  <- list(n = s$n, prevalence = s$prevalence, incub = s$incub, recov = s$recov)
  bnd    <- seir_bounds_mat()
  d      <- length(SEIR_PARS)

  beta_to_p_local <- function(beta) {
    crate <- max(as.numeric(beta), 1.0)
    ptran <- as.numeric(beta) / crate
    list(n = s$n, prevalence = s$prevalence, incub = s$incub, recov = s$recov,
         crate = crate, ptran = ptran)
  }

  eval_p <- function(p) {
    seir_run_multi_ci(p, ndays = ndays, nreps = eval_reps,
                      nthreads = nthreads, seed = as.integer(sid))
  }

  smape_fn <- function(obs, pred)
    mean(2 * abs(obs - pred) / (abs(obs) + abs(pred) + 1e-6)) * 100

  distance_win <- function(theta) {
    p    <- seir_resolve(theta, known)
    pred <- seir_run_single(p, ndays = ndays)
    win_pred <- pred[win_start:win_end]
    sqrt(mean((win_pred - obs_win)^2))
  }

  rmvn_fn <- function(mu, chol_S) as.numeric(mu + crossprod(chol_S, rnorm(length(mu))))
  dmvn_fn <- function(x, mu, S_inv, det_S) {
    dv <- as.numeric(x - mu)
    as.numeric(exp(-0.5 * dv %*% S_inv %*% dv) / sqrt((2 * pi)^d * det_S))
  }
  wq_fn <- function(x, w, probs) {
    o <- order(x); xs <- x[o]; ws <- w[o] / sum(w); cw <- cumsum(ws)
    vapply(probs, function(p) xs[which(cw >= p)[1]], numeric(1))
  }

  out_rows <- list()

  # -- BiLSTM (Tuned) ----------------------------------------------------------
  bl_row <- bl_sub[bl_sub$sim_idx == sid, ]
  if (nrow(bl_row) > 0) {
    bl_beta <- bl_row$beta_pred[1]
    bl_R0   <- bl_row$R0_pred[1]
    p_bl    <- beta_to_p_local(bl_beta)
    q_bl    <- tryCatch(eval_p(p_bl), error = function(e) NULL)
    if (!is.null(q_bl)) {
      pred_bl <- q_bl$med[-1]
      out_rows[[length(out_rows) + 1L]] <- data.frame(
        sim_idx  = sid, method = "BiLSTM",
        pred_beta = bl_beta, pred_R0 = bl_R0,
        true_beta = s$beta,  true_R0 = s$R0,
        smape     = smape_fn(obs_cmp, pred_bl),
        mae       = mean(abs(obs_cmp - pred_bl)),
        pearson_r = cor(obs_cmp, pred_bl),
        err_beta_pct = abs(bl_beta - s$beta) / s$beta * 100,
        err_R0_pct   = abs(bl_R0   - s$R0)   / s$R0   * 100,
        time_sec  = NA_real_, n_model_runs = eval_reps,
        stringsAsFactors = FALSE
      )
    }
  }

  # -- ABC-MCMC ----------------------------------------------------------------
  n_ev <- 0L
  t_abc <- proc.time()
  abc_fit <- tryCatch({
    rms_obs <- sqrt(mean(obs_win^2))
    eps     <- max(abc_eps * rms_obs, 1e-3)
    prop_sd <- abc_prop * (bnd$hi - bnd$lo)
    kern    <- function(dd) exp(-(dd^2) / (2 * eps^2))
    cur     <- (bnd$lo + bnd$hi) / 2
    if (!seir_feasible(cur, known))
      for (ki in 1:500) { cand <- bnd$lo + runif(d)*(bnd$hi-bnd$lo); if (seir_feasible(cand,known)){cur<-cand;break} }
    cur_d <- distance_win(cur); cur_k <- kern(cur_d)
    chain <- matrix(NA_real_, abc_n, d, dimnames = list(NULL, SEIR_PARS))
    for (it in seq_len(abc_n)) {
      prop <- cur + rnorm(d, sd = prop_sd)
      if (seir_feasible(prop, known)) {
        pd <- distance_win(prop); pk <- kern(pd)
        acc <- if (cur_k <= 0 && pk <= 0) pd < cur_d else
          runif(1) < min(1, pk / max(cur_k, .Machine$double.xmin))
        if (acc) { cur <- prop; cur_d <- pd; cur_k <- pk }
      }
      chain[it, ] <- cur
    }
    post       <- chain[(abc_burnin + 1L):abc_n, , drop = FALSE]
    pred_theta <- setNames(apply(post, 2, median), SEIR_PARS)
    list(ok = TRUE, theta = pred_theta)
  }, error = function(e) list(ok = FALSE))
  tt_abc <- proc.time() - t_abc

  if (abc_fit$ok) {
    p_abc    <- seir_resolve(abc_fit$theta, known)
    q_abc    <- tryCatch(eval_p(p_abc), error = function(e) NULL)
    if (!is.null(q_abc)) {
      pred_abc <- q_abc$med[-1]
      abc_beta <- p_abc$beta
      out_rows[[length(out_rows) + 1L]] <- data.frame(
        sim_idx  = sid, method = "ABC",
        pred_beta = abc_beta, pred_R0 = abc_beta / s$recov,
        true_beta = s$beta,   true_R0 = s$R0,
        smape     = smape_fn(obs_cmp, pred_abc),
        mae       = mean(abs(obs_cmp - pred_abc)),
        pearson_r = cor(obs_cmp, pred_abc),
        err_beta_pct = abs(abc_beta - s$beta) / s$beta * 100,
        err_R0_pct   = abs(abc_beta/s$recov - s$R0) / s$R0 * 100,
        time_sec  = unname(tt_abc[["elapsed"]]), n_model_runs = eval_reps,
        stringsAsFactors = FALSE
      )
    }
  }

  # -- ABC-SMC -----------------------------------------------------------------
  n_ev <- 0L
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
      S  <- 2 * (t(dv) %*% (dv * w)) + diag(1e-10 * pmax(diag(2*(t(dv)%*%(dv*w))), 1), d)
      chol_S <- tryCatch(chol(S), error = function(e) diag(sqrt(pmax(diag(S), 1e-12)), d))
      S_inv  <- solve(S); det_S <- det(S)
      new_theta <- matrix(NA_real_, smc_np, d, dimnames = list(NULL, SEIR_PARS))
      new_dist  <- numeric(smc_np)
      i <- 1L; atts <- 0L
      while (i <= smc_np && atts < 25L * smc_np && n_ev < smc_bud) {
        atts <- atts + 1L
        idx  <- sample.int(smc_np, 1L, prob = w)
        cand <- rmvn_fn(theta[idx, ], chol_S)
        if (!seir_feasible(cand, known)) next
        dd <- distance_win(cand); n_ev <- n_ev + 1L
        if (dd <= eps) { new_theta[i, ] <- cand; new_dist[i] <- dd; i <- i + 1L }
      }
      if (i <= smc_np) { ok <- FALSE; break }
      dens  <- vapply(seq_len(smc_np), function(k)
        sum(vapply(seq_len(smc_np), function(j)
          w[j] * dmvn_fn(new_theta[k, ], theta[j, ], S_inv, det_S), numeric(1))),
        numeric(1))
      w     <- 1 / pmax(dens, .Machine$double.xmin); w <- w / sum(w)
      theta <- new_theta; dist <- new_dist
      eps   <- as.numeric(quantile(dist, smc_eq, names = FALSE))
    }
    pred_theta <- setNames(vapply(seq_len(d), function(j)
      wq_fn(theta[, j], w, 0.5), numeric(1)), SEIR_PARS)
    list(ok = ok, theta = pred_theta)
  }, error = function(e) list(ok = FALSE))
  tt_smc <- proc.time() - t_smc

  if (smc_fit$ok) {
    p_smc    <- seir_resolve(smc_fit$theta, known)
    q_smc    <- tryCatch(eval_p(p_smc), error = function(e) NULL)
    if (!is.null(q_smc)) {
      pred_smc <- q_smc$med[-1]
      smc_beta <- p_smc$beta
      out_rows[[length(out_rows) + 1L]] <- data.frame(
        sim_idx  = sid, method = "ABC-SMC",
        pred_beta = smc_beta, pred_R0 = smc_beta / s$recov,
        true_beta = s$beta,   true_R0 = s$R0,
        smape     = smape_fn(obs_cmp, pred_smc),
        mae       = mean(abs(obs_cmp - pred_smc)),
        pearson_r = cor(obs_cmp, pred_smc),
        err_beta_pct = abs(smc_beta - s$beta) / s$beta * 100,
        err_R0_pct   = abs(smc_beta/s$recov - s$R0) / s$R0 * 100,
        time_sec  = unname(tt_smc[["elapsed"]]), n_model_runs = eval_reps,
        stringsAsFactors = FALSE
      )
    }
  }

  # -- Nelder-Mead -------------------------------------------------------------
  n_ev <- 0L
  t_nm <- proc.time()
  nm_fit <- tryCatch({
    obj <- function(theta_unc) {
      theta <- seir_from_unc(theta_unc)
      if (!seir_feasible(theta, known)) return(.Machine$double.xmax / 1e6)
      p    <- seir_resolve(theta, known)
      pred <- seir_run_single(p, ndays = ndays)
      mean((pred[win_start:win_end] - obs_win)^2)
    }
    theta0 <- (bnd$lo + bnd$hi) / 2
    if (!seir_feasible(theta0, known))
      for (ki in 1:200) { cand <- bnd$lo + runif(d)*(bnd$hi-bnd$lo); if(seir_feasible(cand,known)){theta0<-cand;break} }
    r <- optim(seir_to_unc(theta0), obj, method = "Nelder-Mead",
               control = list(maxit = nm_maxit, reltol = nm_reltol))
    for (ki in seq_len(nm_restart))
      r <- optim(r$par, obj, method = "Nelder-Mead",
                 control = list(maxit = nm_maxit, reltol = nm_reltol))
    list(ok = TRUE, theta = seir_from_unc(r$par))
  }, error = function(e) list(ok = FALSE))
  tt_nm <- proc.time() - t_nm

  if (nm_fit$ok) {
    p_nm    <- seir_resolve(nm_fit$theta, known)
    q_nm    <- tryCatch(eval_p(p_nm), error = function(e) NULL)
    if (!is.null(q_nm)) {
      pred_nm <- q_nm$med[-1]
      nm_beta <- p_nm$beta
      out_rows[[length(out_rows) + 1L]] <- data.frame(
        sim_idx  = sid, method = "NelderMead",
        pred_beta = nm_beta, pred_R0 = nm_beta / s$recov,
        true_beta = s$beta,  true_R0 = s$R0,
        smape     = smape_fn(obs_cmp, pred_nm),
        mae       = mean(abs(obs_cmp - pred_nm)),
        pearson_r = cor(obs_cmp, pred_nm),
        err_beta_pct = abs(nm_beta - s$beta) / s$beta * 100,
        err_R0_pct   = abs(nm_beta/s$recov - s$R0) / s$R0 * 100,
        time_sec  = unname(tt_nm[["elapsed"]]), n_model_runs = eval_reps,
        stringsAsFactors = FALSE
      )
    }
  }

  # -- Differential Evolution --------------------------------------------------
  n_ev <- 0L
  t_de <- proc.time()
  de_fit <- tryCatch({
    obj_de <- function(theta) {
      theta <- pmin(pmax(theta, bnd$lo), bnd$hi)
      if (!seir_feasible(theta, known)) return(.Machine$double.xmax / 1e6)
      p    <- seir_resolve(theta, known)
      pred <- seir_run_single(p, ndays = ndays)
      mean((pred[win_start:win_end] - obs_win)^2)
    }
    r <- DEoptimR::JDEoptim(lower = bnd$lo, upper = bnd$hi, fn = obj_de,
                             NP = 10L * d, maxiter = de_maxiter, tol = de_tol)
    list(ok = TRUE, theta = setNames(as.numeric(r$par), SEIR_PARS))
  }, error = function(e) list(ok = FALSE))
  tt_de <- proc.time() - t_de

  if (de_fit$ok) {
    p_de    <- seir_resolve(de_fit$theta, known)
    q_de    <- tryCatch(eval_p(p_de), error = function(e) NULL)
    if (!is.null(q_de)) {
      pred_de <- q_de$med[-1]
      de_beta <- p_de$beta
      out_rows[[length(out_rows) + 1L]] <- data.frame(
        sim_idx  = sid, method = "DE",
        pred_beta = de_beta, pred_R0 = de_beta / s$recov,
        true_beta = s$beta,  true_R0 = s$R0,
        smape     = smape_fn(obs_cmp, pred_de),
        mae       = mean(abs(obs_cmp - pred_de)),
        pearson_r = cor(obs_cmp, pred_de),
        err_beta_pct = abs(de_beta - s$beta) / s$beta * 100,
        err_R0_pct   = abs(de_beta/s$recov - s$R0) / s$R0 * 100,
        time_sec  = unname(tt_de[["elapsed"]]), n_model_runs = eval_reps,
        stringsAsFactors = FALSE
      )
    }
  }

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
  cat(sprintf("Submitting %d SLURM jobs (30-day window, no seed, tuned model)...\n", n_sims))

  job <- Slurm_lapply(
    X   = as.list(seq_len(n_sims)),
    FUN = function(i) {
      w30_one_sim(i, test_params, inc_matrix, bl_sub,
                  SEIR_COMMON, METHOD_ORDER,
                  NDAYS, WIN_START, WIN_END,
                  EVAL_REPS, NTHREADS,
                  ABC_N, ABC_BURNIN, ABC_PROP, ABC_EPS,
                  SMC_NP, SMC_NG, SMC_EQ, SMC_BUD,
                  NM_MAXIT, NM_RELTOL, NM_RESTART,
                  DE_MAXITER, DE_TOL)
    },
    njobs      = min(n_sims, 500),
    mc.cores   = 1,
    job_name   = "w30_tuned",
    plan       = "submit",
    sbatch_opt = SLURM_OPTS,
    export     = c("w30_one_sim",
                   "test_params", "inc_matrix", "bl_sub",
                   "SEIR_COMMON", "METHOD_ORDER",
                   "NDAYS", "WIN_START", "WIN_END",
                   "EVAL_REPS", "NTHREADS",
                   "ABC_N", "ABC_BURNIN", "ABC_PROP", "ABC_EPS",
                   "SMC_NP", "SMC_NG", "SMC_EQ", "SMC_BUD",
                   "NM_MAXIT", "NM_RELTOL", "NM_RESTART",
                   "DE_MAXITER", "DE_TOL"),
    tmp_path   = SCRATCH
  )

  cat("\nJobs submitted. When done:\n")
  cat("  Rscript window_30day_submit_tuned.R --collect\n\n")
  quit(save = "no")
}

# =============================================================================
# STEP 2 — Collect + Plots
# =============================================================================

cat("Collecting 30-day window results (tuned model)...\n")

job_path    <- file.path(path.expand(SCRATCH), "w30_tuned")
results_raw <- Slurm_collect(read_slurm_job(job_path), any. = TRUE)
results_raw <- Filter(is.data.frame, results_raw)
cat(sprintf("Non-null results: %d / %d\n", length(results_raw), nrow(test_params)))

metrics_all <- bind_rows(results_raw) |>
  mutate(method = factor(method, levels = METHOD_ORDER))

write.csv(metrics_all,
          file.path(TUNED_DIR, "w30_5method_metrics.csv"),
          row.names = FALSE)
cat("Saved: tuned_w30/w30_5method_metrics.csv\n")

summary_tbl <- metrics_all |>
  group_by(method) |>
  summarise(
    n_sims        = n(),
    mean_smape    = mean(smape,        na.rm = TRUE),
    median_smape  = median(smape,      na.rm = TRUE),
    mean_mae      = mean(mae,          na.rm = TRUE),
    mean_pearson  = mean(pearson_r,    na.rm = TRUE),
    mean_err_beta = mean(err_beta_pct, na.rm = TRUE),
    mean_err_R0   = mean(err_R0_pct,   na.rm = TRUE),
    mean_time_sec = mean(time_sec,     na.rm = TRUE),
    .groups = "drop"
  )

cat("\n========== 30-Day Window (No Seed) — Tuned BiLSTM — 5-Method Summary ==========\n")
print(as.data.frame(summary_tbl), digits = 4, row.names = FALSE)

write.csv(summary_tbl,
          file.path(TUNED_DIR, "w30_5method_summary.csv"),
          row.names = FALSE)

n_sims <- n_distinct(metrics_all$sim_idx)

# -- Plot 1: sMAPE violin + boxplot ------------------------------------------

p1 <- ggplot(metrics_all, aes(x = method, y = smape, fill = method)) +
  geom_violin(alpha = 0.55, scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.10, outlier.size = 0.7, alpha = 0.85) +
  geom_jitter(width = 0.07, size = 1.5, alpha = 0.5, color = "black") +
  scale_fill_manual(values = METHOD_COLORS) +
  labs(title    = "sMAPE: First 30 Days (No Seed) — Tuned BiLSTM",
       subtitle = sprintf("N=%d sims | calibrate on days 2-31 | evaluate on days 2-365",
                          n_sims),
       x = NULL, y = "sMAPE (%)", fill = "Method") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "none",
        plot.title = element_text(face = "bold"))

ggsave(file.path(PLOTS_DIR, "w30_5method_smape.png"),
       p1, width = 8, height = 5, dpi = 150)
cat("Saved: w30_5method_smape.png\n")

# -- Plot 2: all metrics boxplots --------------------------------------------

metrics_long <- metrics_all |>
  select(sim_idx, method, smape, mae, pearson_r, err_beta_pct, err_R0_pct) |>
  pivot_longer(-c(sim_idx, method), names_to = "metric", values_to = "value")

metric_labels <- c(
  smape        = "sMAPE (%)",
  mae          = "MAE (incidence)",
  pearson_r    = "Pearson r (shape)",
  err_beta_pct = "beta Error (%)",
  err_R0_pct   = "R0 Error (%)"
)
metrics_long$metric <- factor(metric_labels[metrics_long$metric],
                              levels = unname(metric_labels))

p2 <- ggplot(metrics_long, aes(x = method, y = value, fill = method)) +
  geom_boxplot(alpha = 0.65, outlier.size = 0.7) +
  geom_jitter(width = 0.15, size = 1.6, alpha = 0.55, color = "black") +
  facet_wrap(~ metric, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = METHOD_COLORS) +
  labs(title    = "All Metrics: First 30 Days (No Seed) — Tuned BiLSTM",
       subtitle = sprintf("N=%d sims | calibrate on days 2-31 | evaluate on days 2-365",
                          n_sims),
       x = NULL, y = "Value", fill = "Method") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold"), strip.text = element_text(size = 9),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 7))

ggsave(file.path(PLOTS_DIR, "w30_5method_all_metrics.png"),
       p2, width = 14, height = 7, dpi = 150)
cat("Saved: w30_5method_all_metrics.png\n")

# -- Plot 3: predicted vs true beta / R0 ------------------------------------

p3a <- ggplot(metrics_all,
              aes(x = true_beta, y = pred_beta, color = method, shape = method)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_color_manual(values = METHOD_COLORS) +
  labs(title = "Predicted vs True beta (30-day window, tuned)",
       x = "True beta", y = "Predicted beta", color = "Method", shape = "Method") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

p3b <- ggplot(metrics_all,
              aes(x = true_R0, y = pred_R0, color = method, shape = method)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_color_manual(values = METHOD_COLORS) +
  labs(title = "Predicted vs True R0 (30-day window, tuned)",
       x = "True R0", y = "Predicted R0", color = "Method", shape = "Method") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

ggsave(file.path(PLOTS_DIR, "w30_5method_param_scatter.png"),
       plot_grid(p3a, p3b, ncol = 2), width = 13, height = 5, dpi = 150)
cat("Saved: w30_5method_param_scatter.png\n")

# -- Plot 4: timing bar ------------------------------------------------------

timing_tbl <- summary_tbl |>
  filter(!is.na(mean_time_sec)) |>
  mutate(label = ifelse(method == "BiLSTM",
                        sprintf("%.3f s\n(amortised)", mean_time_sec),
                        sprintf("%.1f s", mean_time_sec)))

if (nrow(timing_tbl) > 0) {
  p4 <- ggplot(timing_tbl, aes(x = method, y = mean_time_sec, fill = method)) +
    geom_col(alpha = 0.80, width = 0.6) +
    geom_text(aes(label = label), vjust = -0.4, size = 3.2) +
    scale_y_log10() +
    scale_fill_manual(values = METHOD_COLORS) +
    labs(title    = "Calibration Time — 30-Day Window (No Seed, Tuned BiLSTM)",
         subtitle = "Log scale | BiLSTM = forward pass only",
         x = NULL, y = "Mean wall time (s, log scale)", fill = "Method") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "none",
          plot.title = element_text(face = "bold"))

  ggsave(file.path(PLOTS_DIR, "w30_5method_timing.png"),
         p4, width = 7, height = 5, dpi = 150)
  cat("Saved: w30_5method_timing.png\n")
}

cat(sprintf("\nAll tuned 30-day window outputs saved to:\n  %s\n", TUNED_DIR))
