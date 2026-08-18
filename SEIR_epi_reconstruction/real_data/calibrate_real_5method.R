# =============================================================================
#  real_data/calibrate_real_5method.R
#
#  Calibrates all 5 methods on real incidence datasets using early windows of
#  30, 60, and 90 days. For each method x window:
#    1. Estimate beta (and R0) from the observed window
#    2. Run SEIRconn 1000 times to get a 95% CI band
#    3. Plot the CI band + real data as a solid line
#
#  BiLSTM: uses Bernardo model predictions from bernardo_predict_real.py
#  ABC / ABC-SMC / NM / DE: calibrated in real-time here
#
#  Prerequisites:
#    Rscript real_data/prepare_utah_covid.R
#    Rscript real_data/prepare_measles.R
#    python  real_data/bernardo_predict_real.py
#
#  Usage:
#    Rscript real_data/calibrate_real_5method.R
# =============================================================================

suppressPackageStartupMessages({
  library(epiworldR)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(cowplot)
  library(DEoptimR)
  library(jsonlite)
})

PROJECT_DIR <- path.expand(
  "~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction"
)
REAL_DIR  <- file.path(PROJECT_DIR, "real_data")
PLOTS_DIR <- file.path(REAL_DIR, "plots")
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)

source(file.path(PROJECT_DIR, "seir_common.R"))
seir_set_libpath()

# =============================================================================
# Settings
# =============================================================================

WINDOWS    <- c(60L)                    # calibration window length (days) - single
                                         # window, mirroring the SIR vignette's
                                         # fixed 61-day analysis window
EVAL_REPS  <- 2000L                     # SEIR ensemble reps for CI (matches SIR vignette)
NTHREADS   <- 12L

# Population scales: read from the single shared config, NOT hard-coded here.
# bernardo_predict_real.py reads the same file for its own scale (n_bilstm).
# Do not duplicate these numbers in this script -- edit seir_scale_config.json
# instead, or the two scripts can silently drift out of sync again (this is
# exactly the bug pattern found 2026-08-17 on the SIR side: three scripts
# each hard-coding their own population/recovery-rate assumption with no
# cross-check). See seir_scale_config.json for the full rationale on why
# there are two different scales (BiLSTM vs classical methods).
scale_cfg <- jsonlite::fromJSON(file.path(REAL_DIR, "seir_scale_config.json"))
N_SCALED  <- as.integer(scale_cfg$n_classical)
N_BILSTM  <- as.integer(scale_cfg$n_bilstm)  # only used to sanity-check the BiLSTM CSV below

# Calibration budgets (runs locally, so keep reasonable)
ABC_N      <- 3000L; ABC_BURNIN <- 1500L
ABC_PROP   <- 0.08;  ABC_EPS    <- 0.30
# SMC_NG/SMC_EQ/SMC_BUD: diagnosed 2026-08-17 -- with the old NG=4, EQ=0.5,
# eps only shrank 163->127 over 4 generations while a good fit needs eps~24,
# so ABC-SMC always terminated early with beta still an order of magnitude
# too high (R0 in the 7-17 range vs ABC/NM/DE's ~1.5-5.7). Instrumented runs
# showed the population converges to the same beta as ABC/NM/DE (weighted
# median matching to 3 sig figs) by generation 6-7; generation 8 buys
# negligible extra accuracy for roughly double the cost of generation 7 (the
# per-generation eval cost escalates fast as eps nears the single-run
# objective's noise floor) and was the generation observed to intermittently
# stall/fail depending on ambient RNG state (see the 50x attempt-cap comment
# in calibrate_smc()). NG=7 stops right at the plateau. BUD=20000 gives
# headroom for the escalating per-generation cost near the noise floor
# (worst observed cumulative usage through gen 7 was ~5334).
# SMC_ATT_MULT (per-generation attempt cap = SMC_ATT_MULT * SMC_NP): 50 was
# tuned for wave1/current61's eval cost (~0.1-0.3s/eval at N=7000). Smaller
# populations (e.g. n=187) evaluate in ~1ms, so the same 50x cap can starve
# out a generation (observed 2026-08-17: stalled at gen 5/7, 72/150 filled,
# after the full 50x=7500 attempts) purely because eps has room to keep
# shrinking fast relative to how many candidates that eval-budget can try --
# not because the population is failing to converge. Datasets with cheap
# evals should override this (and SMC_NG/SMC_BUD) upward after sourcing.
SMC_NP     <- 150L;  SMC_NG     <- 7L;  SMC_EQ <- 0.3; SMC_BUD <- 20000L
SMC_ATT_MULT <- 50L
NM_MAXIT   <- 2000L; NM_RELTOL  <- 1e-8; NM_RESTART <- 3L
DE_MAXITER <- 150L;  DE_TOL     <- 1e-8

METHOD_ORDER  <- SEIR_METHOD_ORDER[SEIR_METHOD_ORDER != "Oracle"]
METHOD_COLORS <- SEIR_METHOD_COLORS[names(SEIR_METHOD_COLORS) != "Oracle"]

# =============================================================================
# Calibration helpers
# =============================================================================

beta_to_p <- function(beta, known) {
  crate <- max(as.numeric(beta), 1.0)
  ptran <- as.numeric(beta) / crate
  list(n = known$n, prevalence = known$prevalence,
       incub = known$incub, recov = known$recov,
       crate = crate, ptran = ptran)
}

run_ci <- function(beta, known, ndays, reps = EVAL_REPS, nthreads = NTHREADS) {
  p  <- beta_to_p(beta, known)
  tryCatch(
    seir_run_multi_ci(p, ndays = ndays, nreps = reps, nthreads = nthreads),
    error = function(e) NULL
  )
}

# REVERTED 2026-08-17: tried replacing raw pointwise MSE with a shape-based
# summary-statistic distance (total cases / peak / peak day / phase means)
# to fix the explosive-R0 preference diagnosed on the measles1861 smoke
# test. It didn't fix that dataset (the real cause turned out to be
# n_pop=188 forcing a ~100% attack rate within the observed window -- a
# structural well-mixed-model limitation documented in prepare_measles.R
# and prepare_measles1861_epiestim.R, not something any distance function
# can fix) AND it broke ABC on datasets that were working fine: current61's
# ABC went from R0~1.5-1.7 to R0=64.4, wave1's from R0~4.5-5.7 to R0=106.1
# (NM/ABC-SMC/DE stayed roughly stable, so the regression is specific to
# calibrate_abc()'s epsilon -- the empirical-quantile eps introduced
# alongside this change was evidently mis-scaled). Back to raw pointwise
# MSE, which is what's actually validated against current61/wave1.
distance_fn <- function(theta, known, obs_win, win_len, ndays) {
  p    <- seir_resolve(theta, known)
  pred <- seir_run_single(p, ndays = ndays)
  sqrt(mean((pred[seq_len(win_len)] - obs_win)^2))
}

# -- ABC-MCMC -----------------------------------------------------------------
calibrate_abc <- function(obs_win, known, ndays) {
  win_len <- length(obs_win)
  bnd     <- seir_bounds_mat()
  d       <- length(SEIR_PARS)
  prop_sd <- ABC_PROP * (bnd$hi - bnd$lo)
  rms_obs <- sqrt(mean(obs_win^2))
  eps     <- max(ABC_EPS * rms_obs, 1e-3)
  kern    <- function(dd) exp(-(dd^2) / (2 * eps^2))
  cur     <- (bnd$lo + bnd$hi) / 2
  if (!seir_feasible(cur, known))
    for (ki in 1:500) { cand <- bnd$lo + runif(d)*(bnd$hi-bnd$lo)
      if (seir_feasible(cand,known)){cur<-cand;break} }
  cur_d <- distance_fn(cur, known, obs_win, win_len, ndays)
  cur_k <- kern(cur_d)
  chain <- matrix(NA_real_, ABC_N, d, dimnames=list(NULL,SEIR_PARS))
  for (it in seq_len(ABC_N)) {
    prop <- cur + rnorm(d, sd=prop_sd)
    if (seir_feasible(prop, known)) {
      pd <- distance_fn(prop, known, obs_win, win_len, ndays); pk <- kern(pd)
      acc <- if (cur_k<=0 && pk<=0) pd<cur_d else
        runif(1) < min(1, pk/max(cur_k,.Machine$double.xmin))
      if (acc) { cur<-prop; cur_d<-pd; cur_k<-pk }
    }
    chain[it,] <- cur
  }
  post <- chain[(ABC_BURNIN+1L):ABC_N,,drop=FALSE]
  setNames(apply(post,2,median), SEIR_PARS)
}

# -- ABC-SMC ------------------------------------------------------------------
calibrate_smc <- function(obs_win, known, ndays) {
  win_len <- length(obs_win)
  bnd     <- seir_bounds_mat()
  d       <- length(SEIR_PARS)
  rmvn_fn <- function(mu,chol_S) as.numeric(mu+crossprod(chol_S,rnorm(d)))
  dmvn_fn <- function(x,mu,S_inv,det_S) {
    dv <- as.numeric(x-mu)
    as.numeric(exp(-0.5*dv%*%S_inv%*%dv)/sqrt((2*pi)^d*det_S))
  }
  wq_fn <- function(x,w,probs) {
    o<-order(x);xs<-x[o];ws<-w[o]/sum(w);cw<-cumsum(ws)
    vapply(probs,function(p) xs[which(cw>=p)[1]],numeric(1))
  }
  theta  <- matrix(NA_real_,SMC_NP,d,dimnames=list(NULL,SEIR_PARS))
  filled <- 0L; tries <- 0L
  while (filled<SMC_NP && tries<500L*SMC_NP) {
    tries<-tries+1L; cand<-bnd$lo+runif(d)*(bnd$hi-bnd$lo)
    if (seir_feasible(cand,known)) { filled<-filled+1L; theta[filled,]<-cand }
  }
  if (filled<SMC_NP) return(NULL)
  dist <- apply(theta,1,function(th) distance_fn(th,known,obs_win,win_len,ndays))
  w    <- rep(1/SMC_NP,SMC_NP)
  eps  <- as.numeric(quantile(dist,SMC_EQ,names=FALSE))
  n_ev <- 0L
  for (g in 2:SMC_NG) {
    mu<-colSums(w*theta); dv<-sweep(theta,2,mu)
    S<-2*(t(dv)%*%(dv*w))+diag(1e-10*pmax(diag(2*(t(dv)%*%(dv*w))),1),d)
    chol_S<-tryCatch(chol(S),error=function(e) diag(sqrt(pmax(diag(S),1e-12)),d))
    S_inv<-solve(S); det_S<-det(S)
    new_theta<-matrix(NA_real_,SMC_NP,d,dimnames=list(NULL,SEIR_PARS))
    new_dist<-numeric(SMC_NP)
    i<-1L; atts<-0L
    # 50x, not 25x: per-generation eval cost roughly doubles each round as
    # eps approaches the single-stochastic-run objective's noise floor (14s
    # -> 113s across generations, observed 2026-08-17), so a tighter attempt
    # cap intermittently starves out a generation depending on the ambient
    # RNG state (this call runs after calibrate_abc() has already consumed
    # a full 3000-iteration chain's worth of draws, so which trajectory it
    # gets is effectively random per run) -- confirmed on wave1 across 4
    # seeds with this cap: 0/4 stalls vs frequent stalls at 25x.
    while (i<=SMC_NP && atts<SMC_ATT_MULT*SMC_NP && n_ev<SMC_BUD) {
      atts<-atts+1L; idx<-sample.int(SMC_NP,1L,prob=w)
      cand<-rmvn_fn(theta[idx,],chol_S)
      if (!seir_feasible(cand,known)) next
      dd<-distance_fn(cand,known,obs_win,win_len,ndays); n_ev<-n_ev+1L
      if (dd<=eps) { new_theta[i,]<-cand; new_dist[i]<-dd; i<-i+1L }
    }
    # Graceful degradation (added 2026-08-17, small-N datasets): per-generation
    # eval cost can explode as eps approaches the single-realization
    # objective's noise floor -- observed on n=187 doubling each generation
    # (553 -> 1258 -> 3814 -> 11275 -> 23517 -> 38836 evals) even at a very
    # generous 300x attempt cap / 300k budget, without the population having
    # actually stopped improving (weighted median beta was still drifting
    # generation to generation, not oscillating around a fixed point). Rather
    # than discarding every completed generation's work and returning NULL
    # (the old behaviour -- this is what made ABC-SMC FAIL outright on the
    # first measles1861 run), stop and return the last generation that DID
    # complete: a slightly-less-converged answer beats no answer, especially
    # since ABC/NM/DE all reach a comparable region unaided.
    if (i<=SMC_NP) {
      cat(sprintf("    (ABC-SMC stalled at generation %d, %d/%d filled -- using generation %d's result)\n",
                  g, i-1L, SMC_NP, g-1L))
      break
    }
    dens<-vapply(seq_len(SMC_NP),function(k)
      sum(vapply(seq_len(SMC_NP),function(j)
        w[j]*dmvn_fn(new_theta[k,],theta[j,],S_inv,det_S),numeric(1))),numeric(1))
    w<-1/pmax(dens,.Machine$double.xmin); w<-w/sum(w)
    theta<-new_theta; dist<-new_dist
    eps<-as.numeric(quantile(dist,SMC_EQ,names=FALSE))
  }
  setNames(vapply(seq_len(d),function(j) wq_fn(theta[,j],w,0.5),numeric(1)),SEIR_PARS)
}

# -- Nelder-Mead --------------------------------------------------------------
calibrate_nm <- function(obs_win, known, ndays) {
  win_len <- length(obs_win)
  bnd     <- seir_bounds_mat(); d <- length(SEIR_PARS)
  obj <- function(z) {
    theta <- seir_from_unc(z)
    if (!seir_feasible(theta,known)) return(.Machine$double.xmax/1e6)
    distance_fn(theta, known, obs_win, win_len, ndays)
  }

  # SEIR_PARS currently has a single calibrated parameter (beta), so this is
  # a 1-D optimization. optim(method="Nelder-Mead") explicitly warns
  # "one-dimensional optimization by Nelder-Mead is unreliable" -- its
  # simplex degenerates with only 2 vertices and it was observed to walk
  # straight to a search-bound corner (beta=10, R0=100 on real Utah COVID
  # data; bug found 2026-08-17). Use optim(method="Brent") instead, which is
  # the method R itself recommends for 1-D problems and requires finite
  # bounds -- the unconstrained logit space saturates well within +-15, so
  # +-20 safely covers the full feasible range without truncating it.
  if (d == 1L) {
    r <- optim(seir_to_unc(theta0 <- (bnd$lo + bnd$hi) / 2), obj,
               method = "Brent", lower = -20, upper = 20,
               control = list(reltol = NM_RELTOL))
    return(seir_from_unc(r$par))
  }

  theta0 <- (bnd$lo+bnd$hi)/2
  if (!seir_feasible(theta0,known))
    for (ki in 1:200) { cand<-bnd$lo+runif(d)*(bnd$hi-bnd$lo)
      if(seir_feasible(cand,known)){theta0<-cand;break} }
  r <- optim(seir_to_unc(theta0), obj, method="Nelder-Mead",
             control=list(maxit=NM_MAXIT,reltol=NM_RELTOL))
  for (ki in seq_len(NM_RESTART))
    r <- optim(r$par, obj, method="Nelder-Mead",
               control=list(maxit=NM_MAXIT,reltol=NM_RELTOL))
  seir_from_unc(r$par)
}

# -- Differential Evolution ---------------------------------------------------
calibrate_de <- function(obs_win, known, ndays) {
  win_len <- length(obs_win)
  bnd     <- seir_bounds_mat(); d <- length(SEIR_PARS)
  obj_de <- function(theta) {
    theta <- pmin(pmax(theta,bnd$lo),bnd$hi)
    if (!seir_feasible(theta,known)) return(.Machine$double.xmax/1e6)
    distance_fn(theta, known, obs_win, win_len, ndays)
  }
  r <- tryCatch(
    DEoptimR::JDEoptim(lower=bnd$lo, upper=bnd$hi, fn=obj_de,
                       NP=10L*d, maxiter=DE_MAXITER, tol=DE_TOL),
    error=function(e) NULL
  )
  if (is.null(r)) return(NULL)
  setNames(as.numeric(r$par), SEIR_PARS)
}

# =============================================================================
# Process one dataset × one window
# =============================================================================

calibrate_one_window <- function(obs, known, win_len, bilstm_beta = NA_real_) {
  ndays   <- length(obs)
  obs_win <- obs[seq_len(win_len)]

  results <- list()
  cat(sprintf("    Window %d days...\n", win_len))

  # -- BiLSTM (Bernardo) -------------------------------------------------------
  if (!is.na(bilstm_beta)) {
    cat("      BiLSTM... ")
    ci <- run_ci(bilstm_beta, known, ndays)
    if (!is.null(ci)) {
      results[["BiLSTM"]] <- list(
        beta = bilstm_beta,
        R0   = bilstm_beta / known$recov,
        ci   = ci
      )
      cat(sprintf("beta=%.4f R0=%.3f\n", bilstm_beta, bilstm_beta/known$recov))
    } else cat("FAILED\n")
  }

  # -- ABC ---------------------------------------------------------------------
  cat("      ABC... ")
  t0 <- proc.time()[["elapsed"]]
  abc_theta <- tryCatch(calibrate_abc(obs_win, known, ndays), error=function(e) NULL)
  cat(sprintf("%.1f s", proc.time()[["elapsed"]]-t0))
  if (!is.null(abc_theta)) {
    abc_p  <- seir_resolve(abc_theta, known)
    abc_ci <- run_ci(abc_p$beta, known, ndays)
    if (!is.null(abc_ci)) {
      results[["ABC"]] <- list(beta=abc_p$beta, R0=abc_p$beta/known$recov, ci=abc_ci)
      cat(sprintf("  beta=%.4f R0=%.3f\n", abc_p$beta, abc_p$beta/known$recov))
    } else cat("  CI FAILED\n")
  } else cat("  FAILED\n")

  # -- ABC-SMC -----------------------------------------------------------------
  cat("      ABC-SMC... ")
  t0 <- proc.time()[["elapsed"]]
  smc_theta <- tryCatch(calibrate_smc(obs_win, known, ndays), error=function(e) NULL)
  cat(sprintf("%.1f s", proc.time()[["elapsed"]]-t0))
  if (!is.null(smc_theta)) {
    smc_p  <- seir_resolve(smc_theta, known)
    smc_ci <- run_ci(smc_p$beta, known, ndays)
    if (!is.null(smc_ci)) {
      results[["ABC-SMC"]] <- list(beta=smc_p$beta, R0=smc_p$beta/known$recov, ci=smc_ci)
      cat(sprintf("  beta=%.4f R0=%.3f\n", smc_p$beta, smc_p$beta/known$recov))
    } else cat("  CI FAILED\n")
  } else cat("  FAILED\n")

  # -- Nelder-Mead -------------------------------------------------------------
  cat("      NM... ")
  t0 <- proc.time()[["elapsed"]]
  nm_theta <- tryCatch(calibrate_nm(obs_win, known, ndays), error=function(e) NULL)
  cat(sprintf("%.1f s", proc.time()[["elapsed"]]-t0))
  if (!is.null(nm_theta)) {
    nm_p  <- seir_resolve(nm_theta, known)
    nm_ci <- run_ci(nm_p$beta, known, ndays)
    if (!is.null(nm_ci)) {
      results[["NelderMead"]] <- list(beta=nm_p$beta, R0=nm_p$beta/known$recov, ci=nm_ci)
      cat(sprintf("  beta=%.4f R0=%.3f\n", nm_p$beta, nm_p$beta/known$recov))
    } else cat("  CI FAILED\n")
  } else cat("  FAILED\n")

  # -- DE ----------------------------------------------------------------------
  cat("      DE... ")
  t0 <- proc.time()[["elapsed"]]
  de_theta <- tryCatch(calibrate_de(obs_win, known, ndays), error=function(e) NULL)
  cat(sprintf("%.1f s", proc.time()[["elapsed"]]-t0))
  if (!is.null(de_theta)) {
    de_p  <- seir_resolve(de_theta, known)
    de_ci <- run_ci(de_p$beta, known, ndays)
    if (!is.null(de_ci)) {
      results[["DE"]] <- list(beta=de_p$beta, R0=de_p$beta/known$recov, ci=de_ci)
      cat(sprintf("  beta=%.4f R0=%.3f\n", de_p$beta, de_p$beta/known$recov))
    } else cat("  CI FAILED\n")
  } else cat("  FAILED\n")

  results
}

# =============================================================================
# Plot: all 5 method CI ribbons + real data line
# =============================================================================

plot_ci_all_methods <- function(obs, window_results, win_len,
                                 dates = NULL, dataset_label = "",
                                 out_png = NULL) {
  ndays <- length(obs)
  xv    <- if (!is.null(dates)) dates else seq_len(ndays)

  df_obs <- data.frame(x = xv, y = obs)

  ribbon_rows <- list()
  for (meth in names(window_results)) {
    ci <- window_results[[meth]]$ci
    if (is.null(ci)) next
    ribbon_rows[[meth]] <- data.frame(
      x      = xv,
      lower  = ci$lower,
      med    = ci$med,
      upper  = ci$upper,
      method = meth
    )
  }

  if (length(ribbon_rows) == 0) {
    cat("    No CI data to plot.\n")
    return(invisible(NULL))
  }

  df_ribbons <- bind_rows(ribbon_rows) |>
    mutate(method = factor(method, levels = METHOD_ORDER))

  avail_methods <- levels(droplevels(df_ribbons$method))
  colors_use    <- METHOD_COLORS[avail_methods]

  g <- ggplot() +
    geom_ribbon(data = df_ribbons,
                aes(x = x, ymin = lower, ymax = upper, fill = method),
                alpha = 0.20) +
    geom_line(data = df_ribbons,
              aes(x = x, y = med, colour = method),
              linewidth = 0.9) +
    # vertical line marking end of calibration window
    geom_vline(xintercept = if (!is.null(dates)) dates[win_len] else win_len,
               linetype = "dashed", colour = "grey40", linewidth = 0.6) +
    annotate("text",
             x     = if (!is.null(dates)) dates[win_len] else win_len,
             y     = Inf, vjust = 1.5, hjust = -0.1, size = 3,
             colour = "grey40",
             label = sprintf("day %d", win_len)) +
    # real data on top
    geom_line(data = df_obs,
              aes(x = x, y = y),
              colour = "black", linewidth = 0.7) +
    geom_point(data = df_obs,
               aes(x = x, y = y),
               colour = "black", size = 0.6, alpha = 0.5) +
    scale_fill_manual(values   = colors_use, name = "Method") +
    scale_colour_manual(values = colors_use, name = "Method") +
    labs(
      title    = sprintf("%s — %d-day window", dataset_label, win_len),
      subtitle = sprintf(
        "Black line = real data | Ribbons = 95%% CI (%d SEIR runs) | Dashed = calibration cutoff",
        EVAL_REPS),
      x = if (!is.null(dates)) "Date" else "Day",
      y = "Daily incidence"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position    = "top",
      panel.grid.minor   = element_blank(),
      plot.title         = element_text(face = "bold")
    )

  if (!is.null(out_png)) {
    ggsave(out_png, g, width = 10, height = 5, dpi = 150)
    cat(sprintf("    Saved: %s\n", basename(out_png)))
  }
  invisible(g)
}

# =============================================================================
# Process one full dataset (all windows)
# =============================================================================

process_dataset <- function(label, inc_csv, meta_csv, bilstm_csv) {

  cat(sprintf("\n=== %s ===\n", label))
  slug <- gsub("[^A-Za-z0-9]", "_", tolower(label))

  if (!file.exists(inc_csv) || !file.exists(meta_csv)) {
    cat("  Missing CSV — skipping.\n"); return(invisible(NULL))
  }

  inc_df <- read.csv(inc_csv)
  meta   <- read.csv(meta_csv)[1L, ]
  obs    <- as.numeric(inc_df$daily_cases)
  ndays  <- length(obs)
  dates  <- if ("date" %in% names(inc_df)) as.Date(inc_df$date) else NULL

  n_true <- as.numeric(meta$n_pop)
  # Cap the modeled population at N_SCALED, but do NOT scale obs down by
  # n_use/n_true. obs is real case counts (tens/day) and N_SCALED (5000) is
  # chosen to be the same order of magnitude -- exactly like the SIR
  # vignette, which used N=5000 directly against raw daily_cases with no
  # reference to Utah's true population at all. Diluting obs by a true-pop
  # fraction (as this used to do when N_SCALED was 100000) forced the
  # discrete stochastic simulator into a 0/1/2-cases/day quantized regime
  # during calibration, and rescaling the CI back up amplified every unit of
  # that integer noise into a huge jump in the final plot -- see
  # n_classical_comment in seir_scale_config.json.
  n_use <- if (n_true > N_SCALED) N_SCALED else n_true

  # The observed series starts after cases are already being reported. Convert
  # its opening E -> I incidence into both an Exposed stock (incidence times the
  # latent duration) and an Infected stock (incidence times the infectious
  # duration). ModelSEIRCONN otherwise puts every seed in Exposed, which makes
  # every reconstructed curve start near zero and forces calibration methods
  # to compensate for a phase/initial-condition error.
  init <- seir_initial_conditions(
    obs, n_use, as.numeric(meta$incub_days), as.numeric(meta$recov_rate)
  )
  known <- list(
    n          = n_use,
    prevalence = init$prevalence,
    initial_infected_fraction = init$initial_infected_fraction,
    incub      = as.numeric(meta$incub_days),
    recov      = as.numeric(meta$recov_rate)
  )

  cat(sprintf(paste0(
      "  %d days | n_true=%.0f, n_used=%.0f, initial E=%.1f, ",
      "initial I=%.1f, incub=%.0fd, recov=%.3f\n"),
      ndays, n_true, n_use, init$initial_exposed, init$initial_infected,
      known$incub, known$recov))

  # Load Bernardo predictions if available
  bilstm_preds <- if (file.exists(bilstm_csv)) read.csv(bilstm_csv) else NULL
  if (is.null(bilstm_preds)) {
    cat("  BiLSTM predictions not found — run bernardo_predict_real.py first.\n")
  } else if ("n_bilstm_used" %in% names(bilstm_preds)) {
    # Cross-check: was this CSV generated with the scale seir_scale_config.json
    # currently specifies? Catches a stale CSV after the config changes, or
    # the two scripts drifting apart -- the exact bug pattern found
    # 2026-08-17 (three scripts silently inventing their own assumptions).
    csv_n_bilstm <- unique(bilstm_preds$n_bilstm_used)
    if (length(csv_n_bilstm) != 1L || csv_n_bilstm != N_BILSTM) {
      stop(sprintf(
        "bilstm_csv was generated with n_bilstm=%s but seir_scale_config.json ",
        paste(csv_n_bilstm, collapse = ",")), sprintf(
        "now specifies n_bilstm=%d. Re-run bernardo_predict_real.py before ",
        N_BILSTM), "using this CSV, or the BiLSTM comparison is stale.")
    }
  } else {
    cat("  WARNING: bilstm_csv has no n_bilstm_used column (stale/pre-fix ",
        "output) -- cannot verify it used the current scale config. ",
        "Re-run bernardo_predict_real.py.\n", sep = "")
  }

  all_params <- list()

  for (win_len in WINDOWS) {
    if (win_len > ndays) {
      cat(sprintf("  Skipping %d-day window (data only has %d days)\n",
                  win_len, ndays))
      next
    }

    # Look up Bernardo prediction for early_XXXd window
    win_tag     <- sprintf("early_%03dd", win_len)
    bilstm_beta <- NA_real_
    if (!is.null(bilstm_preds) && win_tag %in% bilstm_preds$window) {
      bilstm_beta <- bilstm_preds$beta_pred[bilstm_preds$window == win_tag][1L]
    }

    cat(sprintf("  --- Window: %d days (BiLSTM beta=%s) ---\n",
        win_len, if (is.na(bilstm_beta)) "NA" else sprintf("%.4f", bilstm_beta)))

    win_res <- calibrate_one_window(obs, known, win_len, bilstm_beta)

    # obs and the simulated CI are already on the same n_use scale -- no
    # rescaling needed (see the n_use/prevalence comment above).
    obs_plot <- obs

    # Save parameter table
    for (meth in names(win_res)) {
      all_params[[length(all_params)+1L]] <- data.frame(
        dataset  = label,
        win_len  = win_len,
        method   = meth,
        beta     = win_res[[meth]]$beta,
        R0       = win_res[[meth]]$R0,
        stringsAsFactors = FALSE
      )
    }

    # Plot
    out_png <- file.path(PLOTS_DIR,
      sprintf("%s_w%03d_5method_ci.png", slug, win_len))
    plot_ci_all_methods(
      obs            = obs_plot,
      window_results = win_res,
      win_len        = win_len,
      dates          = dates,
      dataset_label  = label,
      out_png        = out_png
    )
  }

  if (length(all_params) > 0) {
    params_df <- bind_rows(all_params)
    out_csv   <- file.path(REAL_DIR,
      sprintf("%s_5method_params.csv", slug))
    write.csv(params_df, out_csv, row.names = FALSE)
    cat(sprintf("  Params saved: %s\n", basename(out_csv)))

    cat("\n  ========== Parameter estimates ==========\n")
    print(params_df, row.names = FALSE, digits = 4)
  }

  invisible(NULL)
}

# =============================================================================
# Run datasets
# =============================================================================

if (Sys.getenv("SKIP_WAVE1", "0") != "1") {
  process_dataset(
    label      = "Utah COVID-19 Wave 1",
    inc_csv    = file.path(REAL_DIR, "utah_covid_wave1.csv"),
    meta_csv   = file.path(REAL_DIR, "utah_covid_meta.csv"),
    bilstm_csv = file.path(REAL_DIR, "bernardo_real_covid_predictions.csv")
  )
}

if (Sys.getenv("SKIP_MEASLES", "0") != "1") {
  process_dataset(
    label      = "Measles Hagelloch 1861",
    inc_csv    = file.path(REAL_DIR, "measles_hagelloch_incidence.csv"),
    meta_csv   = file.path(REAL_DIR, "measles_hagelloch_meta.csv"),
    bilstm_csv = file.path(REAL_DIR, "bernardo_real_measles_predictions.csv")
  )
}

cat("\nDone. Plots saved to real_data/plots/\n")
