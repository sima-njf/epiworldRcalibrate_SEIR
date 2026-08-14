# =============================================================================
#  abcsmc_seir_submit.R  — ABC-SMC over the multi-parameter space
#
#  STEP 1 — Submit:   Rscript abcsmc_seir_submit.R
#  STEP 2 — Collect:  Rscript abcsmc_seir_submit.R --collect
#
#  Calibrated : SEIR_PARS from seir_common.R (default crate + ptran)
#  Prior      : independent Uniform over SEIR_BOUNDS
#  Kernel     : multivariate Gaussian, covariance = 2 x weighted covariance of
#               the previous population (Beaumont et al. optimal scaling,
#               multivariate form -- this is what lets the particle cloud align
#               itself with the crate*ptran ridge instead of fighting it)
#  Distance   : RMSE on days 2..365
#
#  THIS IS THE METHOD THAT SHOULD SHINE HERE. crate and ptran are only weakly
#  separable, and a posterior can say so: expect wide marginals for crate and
#  ptran, a strong negative correlation between them, and a tight marginal for
#  their product. The optimisers return one confident point on that ridge with
#  no warning. The posterior correlation is saved so you can quantify it.
#
#  Simulations are deliberately UNSEEDED: in ABC the simulator is the
#  likelihood, so each draw must be an independent realisation.
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

NDAYS            <- SEIR_NDAYS
N_PARTICLES      <- 200L
N_GENERATIONS    <- 6L
EPS_QUANTILE     <- 0.5
MAX_ATTEMPT_MULT <- 25L
SIM_BUDGET       <- 8000L
SAVE_PARTICLES   <- TRUE    # keep the final cloud for posterior plots

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

# =============================================================================
# Small multivariate helpers (avoid a mvtnorm dependency on the workers)
# =============================================================================

seir_rmvn <- function(mu, chol_S) {
  as.numeric(mu + crossprod(chol_S, rnorm(length(mu))))
}

seir_dmvn <- function(x, mu, S_inv, det_S) {
  d  <- length(x)
  dx <- as.numeric(x - mu)
  as.numeric(exp(-0.5 * dx %*% S_inv %*% dx) / sqrt((2 * pi)^d * det_S))
}

seir_wquantile <- function(x, w, probs) {
  o  <- order(x); xs <- x[o]; ws <- w[o] / sum(w); cw <- cumsum(ws)
  vapply(probs, function(p) xs[which(cw >= p)[1]], numeric(1))
}

# =============================================================================
# Worker
# =============================================================================

abcsmc_one_sim <- function(row_idx, test_params, inc_matrix, seir_common,
                           ndays, n_particles, n_generations, eps_quantile,
                           max_attempt_mult, sim_budget, save_particles,
                           rmvn_fn, dmvn_fn, wq_fn) {

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

  set.seed(sim_idx + 900L)
  n_sims_used <- 0L

  distance <- function(theta) {
    n_sims_used <<- n_sims_used + 1L
    p    <- seir_resolve(theta, known)
    pred <- seir_run_single(p, ndays = ndays)
    sqrt(mean((pred[-1] - obs_cmp)^2))
  }

  t0 <- proc.time()

  fit <- tryCatch({

    # ---- Generation 1: uniform prior, accept everything -------------------
    # Rejection-sample the prior so every starting particle is feasible.
    theta <- matrix(NA_real_, n_particles, d,
                    dimnames = list(NULL, SEIR_PARS))
    filled <- 0L; tries <- 0L
    while (filled < n_particles && tries < 500L * n_particles) {
      tries <- tries + 1L
      cand  <- bnd$lo + runif(d) * (bnd$hi - bnd$lo)
      if (seir_feasible(cand, known)) {
        filled <- filled + 1L
        theta[filled, ] <- cand
      }
    }
    if (filled < n_particles) {
      stop("Could not draw a feasible prior sample; check SEIR_BOUNDS.")
    }
    dist <- apply(theta, 1, distance)
    w    <- rep(1 / n_particles, n_particles)
    eps  <- as.numeric(stats::quantile(dist, eps_quantile, names = FALSE))

    trace <- data.frame(generation = 1L, epsilon = NA_real_, next_epsilon = eps,
                        accept_rate = 1.0, n_sims = n_sims_used,
                        ess = n_particles, stringsAsFactors = FALSE)

    gens_done <- 1L
    ok <- TRUE

    for (g in 2:n_generations) {

      prev_theta <- theta
      prev_w     <- w

      # Weighted covariance, then Beaumont scaling. The multivariate form is
      # what lets the kernel tilt along the crate*ptran ridge.
      mu <- colSums(prev_w * prev_theta)
      dv <- sweep(prev_theta, 2, mu)
      S  <- 2 * (t(dv) %*% (dv * prev_w))
      S  <- S + diag(1e-10 * pmax(diag(S), 1), d)   # numerical floor

      chol_S <- tryCatch(chol(S), error = function(e) {
        diag(sqrt(pmax(diag(S), 1e-12)), d)          # fall back to diagonal
      })
      S_inv <- solve(S)
      det_S <- det(S)

      new_theta <- matrix(NA_real_, n_particles, d,
                          dimnames = list(NULL, SEIR_PARS))
      new_dist  <- numeric(n_particles)

      i <- 1L; attempts <- 0L
      max_attempts <- max_attempt_mult * n_particles
      sims_before  <- n_sims_used

      while (i <= n_particles && attempts < max_attempts &&
             n_sims_used < sim_budget) {
        attempts <- attempts + 1L
        idx  <- sample.int(n_particles, 1L, prob = prev_w)
        cand <- rmvn_fn(prev_theta[idx, ], chol_S)
        if (!seir_feasible(cand, known)) next   # box + valid implied ptran
        dd <- distance(cand)
        if (dd <= eps) {
          new_theta[i, ] <- cand
          new_dist[i]    <- dd
          i <- i + 1L
        }
      }

      if (i <= n_particles) {
        cat(sprintf("  [sim_idx %d] gen %d exhausted budget (%d/%d)\n",
                    sim_idx, g, i - 1L, n_particles))
        ok <- FALSE
        break
      }

      dens <- vapply(seq_len(n_particles), function(k) {
        sum(vapply(seq_len(n_particles), function(j) {
          prev_w[j] * dmvn_fn(new_theta[k, ], prev_theta[j, ], S_inv, det_S)
        }, numeric(1)))
      }, numeric(1))

      w     <- 1 / pmax(dens, .Machine$double.xmin)
      w     <- w / sum(w)
      theta <- new_theta
      dist  <- new_dist

      next_eps <- as.numeric(stats::quantile(dist, eps_quantile, names = FALSE))
      ess      <- 1 / sum(w^2)

      trace <- rbind(trace, data.frame(
        generation = g, epsilon = eps, next_epsilon = next_eps,
        accept_rate = n_particles / attempts,
        n_sims = n_sims_used - sims_before, ess = ess,
        stringsAsFactors = FALSE))

      cat(sprintf("  [sim_idx %d] gen %d: eps=%.3f acc=%.3f ESS=%.1f\n",
                  sim_idx, g, eps, n_particles / attempts, ess))

      eps <- next_eps
      gens_done <- g
    }

    list(ok = ok, theta = theta, w = w, trace = trace,
         gens_done = gens_done, final_eps = eps)

  }, error = function(e) {
    cat(sprintf("  [sim_idx %d] ABC-SMC FAILED: %s\n", sim_idx, e$message))
    list(ok = FALSE)
  })

  tt <- proc.time() - t0
  converged <- isTRUE(fit$ok)

  # ---- Posterior summaries --------------------------------------------------
  post_cols <- list()

  if (converged) {
    th <- fit$theta; w <- fit$w

    pred_theta <- vapply(seq_len(ncol(th)),
                         function(j) wq_fn(th[, j], w, 0.5), numeric(1))
    pred_theta <- setNames(pred_theta, SEIR_PARS)
    pred_p     <- seir_resolve(pred_theta, known)

    for (j in seq_along(SEIR_PARS)) {
      nm <- SEIR_PARS[j]
      q  <- wq_fn(th[, j], w, c(0.025, 0.975))
      post_cols[[paste0(nm, "_lo_95")]] <- q[1]
      post_cols[[paste0(nm, "_hi_95")]] <- q[2]
      post_cols[[paste0(nm, "_covered")]] <-
        as.integer(true_p[[nm]] >= q[1] && true_p[[nm]] <= q[2])
    }

    # Derived posteriors propagated particle-by-particle
    beta_draw  <- if ("beta" %in% SEIR_PARS) th[, "beta"] else
      apply(th, 1, function(r) { p <- seir_resolve(r, known); p$crate * p$ptran })
    recov_draw <- if ("recov" %in% SEIR_PARS) th[, "recov"] else rep(known$recov, nrow(th))
    R0_draw    <- beta_draw / recov_draw

    qb <- wq_fn(beta_draw, w, c(0.025, 0.5, 0.975))
    qr <- wq_fn(R0_draw,   w, c(0.025, 0.5, 0.975))
    post_cols$beta_lo_95   <- qb[1]; post_cols$beta_hi_95 <- qb[3]
    post_cols$R0_lo_95     <- qr[1]; post_cols$R0_hi_95   <- qr[3]
    post_cols$beta_covered  <- as.integer(true_p$beta >= qb[1] && true_p$beta <= qb[3])
    post_cols$R0_covered    <- as.integer(
      (true_p$beta / true_p$recov) >= qr[1] && (true_p$beta / true_p$recov) <= qr[3])
    post_cols$cv_beta <- sd(beta_draw) / mean(beta_draw)

    post_cols$ess           <- 1 / sum(w^2)
    post_cols$generations   <- fit$gens_done
    post_cols$final_epsilon <- fit$final_eps

    pred_p_run <- seir_resolve(pred_theta, known)
    pred_inc   <- seir_run_single(pred_p_run, ndays = ndays, seed = sim_idx + 500L)
    pred_cmp   <- pred_inc[-1]

    pred_p_met <- list(beta = pred_theta[["beta"]], recov = known$recov)
    summ <- seir_metrics(sim_idx, "ABC-SMC", obs_cmp, pred_cmp,
                         true_p, pred_p_met, known, converged = TRUE,
                         time_sec = unname(tt[["elapsed"]]),
                         cpu_sec  = unname(tt[["user.self"]] + tt[["sys.self"]]),
                         n_evals = n_sims_used, n_model_runs = n_sims_used + 1L)
  } else {
    pred_cmp   <- rep(NA_real_, length(obs_cmp))
    pred_p_met <- true_p
    summ <- seir_metrics(sim_idx, "ABC-SMC", obs_cmp, pred_cmp,
                         true_p, pred_p_met, known, converged = FALSE,
                         time_sec = unname(tt[["elapsed"]]),
                         cpu_sec  = unname(tt[["user.self"]] + tt[["sys.self"]]),
                         n_evals = n_sims_used, n_model_runs = n_sims_used)
  }

  for (nm in names(post_cols)) summ[[nm]] <- post_cols[[nm]]

  daily <- data.frame(sim_idx = sim_idx, method = "ABC-SMC", day = 2:ndays,
                      obs_inc = obs_cmp, pred_inc = pred_cmp,
                      stringsAsFactors = FALSE)

  particles <- NULL
  if (converged && save_particles) {
    particles <- data.frame(sim_idx = sim_idx, fit$theta, weight = fit$w,
                            stringsAsFactors = FALSE)
  }

  list(summary = summ, daily = daily,
       trace = if (converged) cbind(sim_idx = sim_idx, fit$trace) else NULL,
       particles = particles)
}

# =============================================================================
# STEP 1 — Submit
# =============================================================================

if (!collect_only) {
  n_sims <- nrow(test_params)
  cat(sprintf("Submitting ABC-SMC for %d sims | pars: %s (D=%d)\n",
              n_sims, paste(SEIR_PARS, collapse = ", "), length(SEIR_PARS)))
  cat(sprintf("  Particles %d | Generations %d | eps quantile %.2f | budget %d\n",
              N_PARTICLES, N_GENERATIONS, EPS_QUANTILE, SIM_BUDGET))

  job <- Slurm_lapply(
    X   = as.list(seq_len(n_sims)),
    FUN = function(i) {
      abcsmc_one_sim(i, test_params, inc_matrix, SEIR_COMMON, NDAYS,
                     N_PARTICLES, N_GENERATIONS, EPS_QUANTILE,
                     MAX_ATTEMPT_MULT, SIM_BUDGET, SAVE_PARTICLES,
                     seir_rmvn, seir_dmvn, seir_wquantile)
    },
    njobs      = min(n_sims, 200),
    mc.cores   = 1,
    job_name   = "seir_abcsmc",
    plan       = "submit",
    sbatch_opt = SLURM_OPTS,
    export     = c("abcsmc_one_sim", "test_params", "inc_matrix", "SEIR_COMMON",
                   "NDAYS", "N_PARTICLES", "N_GENERATIONS", "EPS_QUANTILE",
                   "MAX_ATTEMPT_MULT", "SIM_BUDGET", "SAVE_PARTICLES",
                   "seir_rmvn", "seir_dmvn", "seir_wquantile"),
    tmp_path   = SCRATCH
  )

  cat("\nJobs submitted. When finished, run:\n")
  cat("  Rscript abcsmc_seir_submit.R --collect\n\n")
  quit(save = "no")
}

# =============================================================================
# STEP 2 — Collect
# =============================================================================

cat("Collecting ABC-SMC results...\n")

job_path    <- file.path(path.expand(SCRATCH), "seir_abcsmc")
results_raw <- Slurm_collect(read_slurm_job(job_path), any. = TRUE)
results_raw <- Filter(Negate(is.null), results_raw)
cat(sprintf("Non-null results: %d / %d\n", length(results_raw), nrow(test_params)))

all_summary   <- bind_rows(lapply(results_raw, `[[`, "summary"))
all_daily     <- bind_rows(lapply(results_raw, `[[`, "daily"))
all_trace     <- bind_rows(lapply(results_raw, `[[`, "trace"))
all_particles <- bind_rows(lapply(results_raw, `[[`, "particles"))

write.csv(all_summary,   file.path(OUT_DIR, "abcsmc_seir_summary.csv"),   row.names = FALSE)
write.csv(all_daily,     file.path(OUT_DIR, "abcsmc_seir_daily.csv"),     row.names = FALSE)
write.csv(all_trace,     file.path(OUT_DIR, "abcsmc_seir_trace.csv"),     row.names = FALSE)
if (nrow(all_particles) > 0) {
  write.csv(all_particles, file.path(OUT_DIR, "abcsmc_seir_particles.csv"), row.names = FALSE)
  cat("Saved: abcsmc_seir_particles.csv  (final posterior clouds)\n")
}
cat("Saved: abcsmc_seir_summary.csv / _daily.csv / _trace.csv\n")

ok <- all_summary[all_summary$converged %in% TRUE, ]

cat("\n========== ABC-SMC Results ==========\n")
cat(sprintf("Attempted / converged : %d / %d (%.1f%%)\n",
            nrow(all_summary), nrow(ok),
            100 * nrow(ok) / max(nrow(all_summary), 1)))
cat(sprintf("Mean generations      : %.2f / %d\n",
            mean(ok$generations, na.rm = TRUE), N_GENERATIONS))
cat(sprintf("Mean ESS (of %d)      : %.1f\n", N_PARTICLES, mean(ok$ess, na.rm = TRUE)))

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
  cat("A strong negative correlation with CV(beta) << CV(crate) is the ridge:\n")
  cat("the data pin down the product but not the factors. This is the result\n")
  cat("no point estimator can express.\n")
}

cat("\n--- Cost ---\n")
cat(sprintf("Mean wall / CPU sec   : %.2f / %.2f\n",
            mean(ok$time_sec, na.rm = TRUE), mean(ok$cpu_sec, na.rm = TRUE)))
cat(sprintf("Mean model runs / sim : %.0f\n", mean(ok$n_model_runs, na.rm = TRUE)))

if (nrow(all_trace) > 0) {
  cat("\n--- Mean schedule across sims ---\n")
  print(as.data.frame(all_trace |> group_by(generation) |>
    summarise(mean_eps = mean(epsilon, na.rm = TRUE),
              mean_accept = mean(accept_rate, na.rm = TRUE),
              mean_ess = mean(ess, na.rm = TRUE),
              mean_sims = mean(n_sims, na.rm = TRUE), .groups = "drop")),
    digits = 4, row.names = FALSE)
}

cat("\nRun validate_methods.R to compare all methods.\n")
