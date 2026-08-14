# =============================================================================
#  abc_seir_submit.R  — ABC-MCMC over the multi-parameter space
#
#  STEP 1 — Submit:   Rscript abc_seir_submit.R
#  STEP 2 — Collect:  Rscript abc_seir_submit.R --collect
#
#  Calibrated : SEIR_PARS from seir_common.R (default crate + ptran)
#  Prior      : independent Uniform over SEIR_BOUNDS
#  Proposal   : Gaussian random walk in PARAMETER space (symmetric), rejected
#               outside the box. Working in theta space rather than a
#               transformed space keeps the Metropolis ratio free of Jacobian
#               terms, so no correction is needed.
#  Kernel     : exp(-d^2 / 2 eps^2) with d = RMSE on days 2..365
#  Epsilon    : EPS_FRAC x RMS(observed), in incidence-count units
#
#  Kept in the comparison as the naive Bayesian baseline. In more than one
#  dimension a single random-walk chain mixes noticeably worse than ABC-SMC,
#  especially along a ridge, so the accept_rate and the per-parameter
#  acceptance diagnostics matter more here than they did in the 1-D study.
#
#  Simulations are UNSEEDED on purpose (the simulator is the likelihood).
# =============================================================================

PROJECT_DIR <- path.expand(
  "~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction"
)

SEIR_COMMON <- file.path(PROJECT_DIR, "seir_common.R")
source(SEIR_COMMON)
seir_set_libpath()

library(slurmR)
library(dplyr)
library(reticulate)

SCRATCH <- file.path("/scratch/general/vast", Sys.getenv("USER"), "slurmR")
OUT_DIR <- PROJECT_DIR
dir.create(path.expand(SCRATCH), recursive = TRUE, showWarnings = FALSE)

SLURM_OPTS <- list(
  account         = "vegayon-np",
  partition       = "vegayon-np",
  `cpus-per-task` = 1,
  `mem-per-cpu`   = "8G",
  time            = "12:00:00"
)

NDAYS       <- SEIR_NDAYS
N_SAMPLES   <- 3000L    # raised for the higher-dimensional space
BURNIN      <- 1500L
PROP_FRAC   <- 0.08     # proposal sd = this fraction of each box width
EPS_FRAC    <- 0.30     # epsilon = EPS_FRAC * RMS(observed)

# =============================================================================
# Load
# =============================================================================

actual <- read.csv(file.path(PROJECT_DIR, SEIR_TEST_PARAMS_FILE))

np_module <- reticulate::import("numpy")
inc_raw   <- as.matrix(
  np_module$load(path.expand(file.path(PROJECT_DIR, SEIR_TEST_INCIDENCE_NPY)))
)
stopifnot(nrow(inc_raw) == nrow(actual), ncol(inc_raw) == NDAYS)

need <- c("beta", "R0", "recov", "incub", "n", "prevalence")
miss <- setdiff(need, names(actual))
if (length(miss) > 0) {
  stop(SEIR_TEST_PARAMS_FILE, " is missing columns: ", paste(miss, collapse = ", "))
}

test_params <- actual
inc_matrix  <- inc_raw

cat(sprintf("Test set: %d sims | calibrating: %s\n",
            nrow(test_params), paste(SEIR_PARS, collapse = ", ")))

args         <- commandArgs(trailingOnly = TRUE)
collect_only <- "--collect" %in% args

seir_wquantile <- function(x, w, probs) {
  o  <- order(x); xs <- x[o]; ws <- w[o] / sum(w); cw <- cumsum(ws)
  vapply(probs, function(p) xs[which(cw >= p)[1]], numeric(1))
}

# =============================================================================
# Worker
# =============================================================================

abc_one_sim <- function(row_idx, test_params, inc_matrix, seir_common,
                        ndays, n_samples, burnin, prop_frac, eps_frac) {

  source(seir_common)
  seir_set_libpath()
  library(epiworldR)

  s       <- test_params[row_idx, ]
  sim_idx <- s$sim_idx
  d       <- length(SEIR_PARS)
  bnd     <- seir_bounds_mat()

  known  <- list(n = s$n, prevalence = s$prevalence, incub = s$incub, recov = s$recov)
  true_p <- list(beta = s$beta, recov = s$recov)

  obs_inc <- as.numeric(inc_matrix[row_idx, ])
  obs_cmp <- obs_inc[-1]

  set.seed(sim_idx + 100L)
  n_sims_used <- 0L

  distance <- function(theta) {
    n_sims_used <<- n_sims_used + 1L
    p    <- seir_resolve(theta, known)
    pred <- seir_run_single(p, ndays = ndays)
    sqrt(mean((pred[-1] - obs_cmp)^2))
  }

  rms_obs  <- sqrt(mean(obs_cmp^2))
  epsilon  <- max(eps_frac * rms_obs, 1e-3)
  prop_sd  <- prop_frac * (bnd$hi - bnd$lo)

  kern <- function(dd) exp(-(dd^2) / (2 * epsilon^2))

  cat(sprintf("  [sim_idx %d] ABC-MCMC: eps=%.3f RMS(obs)=%.2f D=%d\n",
              sim_idx, epsilon, rms_obs, d))

  t0 <- proc.time()

  fit <- tryCatch({
    cur <- (bnd$lo + bnd$hi) / 2       # same start as the optimisers
    if (!seir_feasible(cur, known)) {
      for (k in 1:500) {
        cand <- bnd$lo + runif(d) * (bnd$hi - bnd$lo)
        if (seir_feasible(cand, known)) { cur <- cand; break }
      }
    }
    cur_d   <- distance(cur)
    cur_k   <- kern(cur_d)

    chain   <- matrix(NA_real_, n_samples, d,
                      dimnames = list(NULL, SEIR_PARS))
    dists   <- numeric(n_samples)
    n_acc   <- 0L

    for (it in seq_len(n_samples)) {
      prop <- cur + rnorm(d, sd = prop_sd)
      if (seir_feasible(prop, known)) {
        pd <- distance(prop)
        pk <- kern(pd)
        # Uniform prior + symmetric proposal => ratio is the kernel ratio.
        # Guard the 0/0 case explicitly instead of letting it propagate.
        acc <- if (cur_k <= 0 && pk <= 0) {
          pd < cur_d          # both underflow: fall back to the distance itself
        } else {
          runif(1) < min(1, pk / max(cur_k, .Machine$double.xmin))
        }
        if (acc) { cur <- prop; cur_d <- pd; cur_k <- pk; n_acc <- n_acc + 1L }
      }
      chain[it, ] <- cur
      dists[it]   <- cur_d
    }

    list(ok = TRUE, chain = chain, dists = dists,
         acc_rate = n_acc / n_samples)
  }, error = function(e) {
    cat(sprintf("  [sim_idx %d] ABC FAILED: %s\n", sim_idx, e$message))
    list(ok = FALSE)
  })

  tt <- proc.time() - t0
  converged <- isTRUE(fit$ok)

  post_cols <- list()

  if (converged) {
    post <- fit$chain[(burnin + 1L):n_samples, , drop = FALSE]
    w    <- rep(1 / nrow(post), nrow(post))

    pred_theta <- setNames(apply(post, 2, median), SEIR_PARS)
    pred_p     <- seir_resolve(pred_theta, known)

    for (j in seq_along(SEIR_PARS)) {
      nm <- SEIR_PARS[j]
      q  <- unname(quantile(post[, j], c(0.025, 0.975)))
      post_cols[[paste0(nm, "_lo_95")]]  <- q[1]
      post_cols[[paste0(nm, "_hi_95")]]  <- q[2]
      post_cols[[paste0(nm, "_covered")]] <-
        as.integer(true_p[[nm]] >= q[1] && true_p[[nm]] <= q[2])
    }

    beta_draw  <- if ("beta" %in% SEIR_PARS) post[, "beta"] else
      apply(post, 1, function(r) { p <- seir_resolve(r, known); p$crate * p$ptran })
    recov_draw <- if ("recov" %in% SEIR_PARS) post[, "recov"] else rep(known$recov, nrow(post))
    R0_draw    <- beta_draw / recov_draw

    qb <- unname(quantile(beta_draw, c(0.025, 0.975)))
    qr <- unname(quantile(R0_draw,   c(0.025, 0.975)))
    post_cols$beta_lo_95 <- qb[1]; post_cols$beta_hi_95 <- qb[2]
    post_cols$R0_lo_95   <- qr[1]; post_cols$R0_hi_95   <- qr[2]
    tb <- true_p$beta
    post_cols$beta_covered <- as.integer(tb >= qb[1] && tb <= qb[2])
    post_cols$R0_covered   <- as.integer((tb / true_p$recov) >= qr[1] &&
                                         (tb / true_p$recov) <= qr[2])

    post_cols$cv_beta <- sd(beta_draw) / mean(beta_draw)

    post_cols$accept_rate   <- fit$acc_rate
    post_cols$epsilon       <- epsilon
    post_cols$n_unique_frac <- nrow(unique(post)) / nrow(post)

    # Build pred_p for SEIR run
    pred_p_run <- seir_resolve(pred_theta, known)
    pred_inc   <- seir_run_single(pred_p_run, ndays = ndays, seed = sim_idx + 200L)
    pred_cmp <- pred_inc[-1]

    pred_p_met <- list(beta = pred_theta[["beta"]], recov = known$recov)
    summ <- seir_metrics(sim_idx, "ABC", obs_cmp, pred_cmp,
                         true_p, pred_p_met, known, converged = TRUE,
                         time_sec = unname(tt[["elapsed"]]),
                         cpu_sec  = unname(tt[["user.self"]] + tt[["sys.self"]]),
                         n_evals = n_sims_used, n_model_runs = n_sims_used + 1L)

    if (fit$acc_rate < 0.02) {
      cat(sprintf("  [sim_idx %d] *** acceptance %.3f - chain barely moved ***\n",
                  sim_idx, fit$acc_rate))
    }
  } else {
    pred_cmp   <- rep(NA_real_, length(obs_cmp))
    pred_p_met <- true_p
    summ <- seir_metrics(sim_idx, "ABC", obs_cmp, pred_cmp,
                         true_p, pred_p_met, known, converged = FALSE,
                         time_sec = unname(tt[["elapsed"]]),
                         cpu_sec  = unname(tt[["user.self"]] + tt[["sys.self"]]),
                         n_evals = n_sims_used, n_model_runs = n_sims_used)
  }

  for (nm in names(post_cols)) summ[[nm]] <- post_cols[[nm]]

  daily <- data.frame(sim_idx = sim_idx, method = "ABC", day = 2:ndays,
                      obs_inc = obs_cmp, pred_inc = pred_cmp,
                      stringsAsFactors = FALSE)

  list(summary = summ, daily = daily)
}

# =============================================================================
# STEP 1 — Submit
# =============================================================================

if (!collect_only) {
  n_sims <- nrow(test_params)
  cat(sprintf("Submitting ABC-MCMC for %d sims | pars: %s (D=%d)\n",
              n_sims, paste(SEIR_PARS, collapse = ", "), length(SEIR_PARS)))
  cat(sprintf("  Samples %d | burn-in %d | proposal %.0f%% of box | eps %.2f x RMS\n",
              N_SAMPLES, BURNIN, 100 * PROP_FRAC, EPS_FRAC))

  job <- Slurm_lapply(
    X   = as.list(seq_len(n_sims)),
    FUN = function(i) {
      abc_one_sim(i, test_params, inc_matrix, SEIR_COMMON, NDAYS,
                  N_SAMPLES, BURNIN, PROP_FRAC, EPS_FRAC)
    },
    njobs      = min(n_sims, 200),
    mc.cores   = 1,
    job_name   = "seir_abc",
    plan       = "submit",
    sbatch_opt = SLURM_OPTS,
    export     = c("abc_one_sim", "test_params", "inc_matrix", "SEIR_COMMON",
                   "NDAYS", "N_SAMPLES", "BURNIN", "PROP_FRAC", "EPS_FRAC"),
    tmp_path   = SCRATCH
  )

  cat("\nJobs submitted. When finished, run:\n")
  cat("  Rscript abc_seir_submit.R --collect\n\n")
  quit(save = "no")
}

# =============================================================================
# STEP 2 — Collect
# =============================================================================

cat("Collecting ABC-MCMC results...\n")

job_path    <- file.path(path.expand(SCRATCH), "seir_abc")
results_raw <- Slurm_collect(read_slurm_job(job_path), any. = TRUE)
results_raw <- Filter(Negate(is.null), results_raw)
cat(sprintf("Non-null results: %d / %d\n", length(results_raw), nrow(test_params)))

all_summary <- bind_rows(lapply(results_raw, `[[`, "summary"))
all_daily   <- bind_rows(lapply(results_raw, `[[`, "daily"))

write.csv(all_summary, file.path(OUT_DIR, "abc_seir_summary.csv"), row.names = FALSE)
write.csv(all_daily,   file.path(OUT_DIR, "abc_seir_daily.csv"),   row.names = FALSE)
cat("Saved: abc_seir_summary.csv / abc_seir_daily.csv\n")

ok <- all_summary[all_summary$converged %in% TRUE, ]

cat("\n========== ABC-MCMC Results ==========\n")
cat(sprintf("Attempted / converged : %d / %d (%.1f%%)\n",
            nrow(all_summary), nrow(ok),
            100 * nrow(ok) / max(nrow(all_summary), 1)))
cat(sprintf("Mean acceptance rate  : %.3f\n", mean(ok$accept_rate, na.rm = TRUE)))
cat(sprintf("Chains below 2%% accept: %d  <- these did not explore\n",
            sum(ok$accept_rate < 0.02, na.rm = TRUE)))

cat("\n--- Parameter recovery (mean %% error of the posterior median) ---\n")
for (nm in c("crate", "ptran", "recov")) {
  tag <- if (nm %in% SEIR_PARS) "calibrated" else "known/fixed"
  cat(sprintf("  %-6s : %7.1f%%   (%s)\n", nm,
              mean(ok[[paste0("err_", nm, "_pct")]], na.rm = TRUE), tag))
}
cat(sprintf("  %-6s : %7.1f%%   (derived)\n", "beta", mean(ok$err_beta_pct, na.rm = TRUE)))
cat(sprintf("  %-6s : %7.1f%%   (derived)\n", "R0",   mean(ok$err_R0_pct,   na.rm = TRUE)))

cat("\n--- 95%% credible interval coverage (nominal 95%%) ---\n")
for (nm in SEIR_PARS) {
  cat(sprintf("  %-6s : %5.1f%%\n", nm,
              100 * mean(ok[[paste0(nm, "_covered")]], na.rm = TRUE)))
}
cat(sprintf("  %-6s : %5.1f%%\n", "beta", 100 * mean(ok$beta_covered, na.rm = TRUE)))
cat(sprintf("  %-6s : %5.1f%%\n", "R0",   100 * mean(ok$R0_covered,   na.rm = TRUE)))

if (!is.null(ok$post_cor_crate_ptran)) {
  cat("\n--- Identifiability signature ---\n")
  cat(sprintf("Mean posterior cor(crate, ptran) : %+.3f\n",
              mean(ok$post_cor_crate_ptran, na.rm = TRUE)))
  cat(sprintf("Mean posterior CV: crate %.3f | ptran %.3f | beta %.3f\n",
              mean(ok$cv_crate, na.rm = TRUE), mean(ok$cv_ptran, na.rm = TRUE),
              mean(ok$cv_beta,  na.rm = TRUE)))
}

cat("\n--- Cost ---\n")
cat(sprintf("Mean wall / CPU sec   : %.2f / %.2f\n",
            mean(ok$time_sec, na.rm = TRUE), mean(ok$cpu_sec, na.rm = TRUE)))
cat(sprintf("Mean model runs / sim : %.0f\n", mean(ok$n_model_runs, na.rm = TRUE)))

cat("\nRun validate_methods.R to compare all methods.\n")
