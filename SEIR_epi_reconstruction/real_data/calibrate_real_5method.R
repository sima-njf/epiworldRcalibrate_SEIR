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

WINDOWS    <- c(30L, 60L, 90L)         # calibration window lengths (days)
EVAL_REPS  <- 400L                      # SEIR ensemble reps for CI
NTHREADS   <- 4L
N_SCALED   <- 100000L                   # scale large populations down to this for memory efficiency

# Calibration budgets (runs locally, so keep reasonable)
ABC_N      <- 3000L; ABC_BURNIN <- 1500L
ABC_PROP   <- 0.08;  ABC_EPS    <- 0.30
SMC_NP     <- 150L;  SMC_NG     <- 4L;  SMC_EQ <- 0.5; SMC_BUD <- 6000L
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
  ok   <- TRUE; n_ev <- 0L
  for (g in 2:SMC_NG) {
    mu<-colSums(w*theta); dv<-sweep(theta,2,mu)
    S<-2*(t(dv)%*%(dv*w))+diag(1e-10*pmax(diag(2*(t(dv)%*%(dv*w))),1),d)
    chol_S<-tryCatch(chol(S),error=function(e) diag(sqrt(pmax(diag(S),1e-12)),d))
    S_inv<-solve(S); det_S<-det(S)
    new_theta<-matrix(NA_real_,SMC_NP,d,dimnames=list(NULL,SEIR_PARS))
    new_dist<-numeric(SMC_NP)
    i<-1L; atts<-0L
    while (i<=SMC_NP && atts<25L*SMC_NP && n_ev<SMC_BUD) {
      atts<-atts+1L; idx<-sample.int(SMC_NP,1L,prob=w)
      cand<-rmvn_fn(theta[idx,],chol_S)
      if (!seir_feasible(cand,known)) next
      dd<-distance_fn(cand,known,obs_win,win_len,ndays); n_ev<-n_ev+1L
      if (dd<=eps) { new_theta[i,]<-cand; new_dist[i]<-dd; i<-i+1L }
    }
    if (i<=SMC_NP) { ok<-FALSE; break }
    dens<-vapply(seq_len(SMC_NP),function(k)
      sum(vapply(seq_len(SMC_NP),function(j)
        w[j]*dmvn_fn(new_theta[k,],theta[j,],S_inv,det_S),numeric(1))),numeric(1))
    w<-1/pmax(dens,.Machine$double.xmin); w<-w/sum(w)
    theta<-new_theta; dist<-new_dist
    eps<-as.numeric(quantile(dist,SMC_EQ,names=FALSE))
  }
  if (!ok) return(NULL)
  setNames(vapply(seq_len(d),function(j) wq_fn(theta[,j],w,0.5),numeric(1)),SEIR_PARS)
}

# -- Nelder-Mead --------------------------------------------------------------
calibrate_nm <- function(obs_win, known, ndays) {
  win_len <- length(obs_win)
  bnd     <- seir_bounds_mat(); d <- length(SEIR_PARS)
  obj <- function(z) {
    theta <- seir_from_unc(z)
    if (!seir_feasible(theta,known)) return(.Machine$double.xmax/1e6)
    pred <- seir_run_single(seir_resolve(theta,known), ndays=ndays)
    mean((pred[seq_len(win_len)] - obs_win)^2)
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
    pred <- seir_run_single(seir_resolve(theta,known), ndays=ndays)
    mean((pred[seq_len(win_len)] - obs_win)^2)
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
  # Scale large populations down to reduce memory usage during calibration.
  # beta is a rate so it is unaffected; incidence scales proportionally.
  scale_factor <- if (n_true > N_SCALED) N_SCALED / n_true else 1.0
  n_use        <- if (n_true > N_SCALED) N_SCALED else n_true
  obs          <- obs * scale_factor   # scaled incidence for calibration

  known <- list(
    n          = n_use,
    prevalence = as.numeric(meta$prevalence),
    incub      = as.numeric(meta$incub_days),
    recov      = as.numeric(meta$recov_rate)
  )

  cat(sprintf("  %d days | n_true=%.0f, n_used=%.0f (scale=%.4f), incub=%.0fd, recov=%.3f\n",
      ndays, n_true, n_use, scale_factor, known$incub, known$recov))

  # Load Bernardo predictions if available
  bilstm_preds <- if (file.exists(bilstm_csv)) read.csv(bilstm_csv) else NULL
  if (is.null(bilstm_preds))
    cat("  BiLSTM predictions not found — run bernardo_predict_real.py first.\n")

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

    # Scale CI back up to original population size for plotting
    if (scale_factor < 1.0) {
      for (meth in names(win_res)) {
        ci <- win_res[[meth]]$ci
        if (!is.null(ci)) {
          win_res[[meth]]$ci$lower <- ci$lower / scale_factor
          win_res[[meth]]$ci$med   <- ci$med   / scale_factor
          win_res[[meth]]$ci$upper <- ci$upper / scale_factor
        }
      }
    }

    # Restore original obs for plotting
    obs_plot <- obs / scale_factor

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

process_dataset(
  label      = "Utah COVID-19 Wave 1",
  inc_csv    = file.path(REAL_DIR, "utah_covid_wave1.csv"),
  meta_csv   = file.path(REAL_DIR, "utah_covid_meta.csv"),
  bilstm_csv = file.path(REAL_DIR, "bernardo_real_covid_predictions.csv")
)

process_dataset(
  label      = "Measles Hagelloch 1861",
  inc_csv    = file.path(REAL_DIR, "measles_hagelloch_incidence.csv"),
  meta_csv   = file.path(REAL_DIR, "measles_hagelloch_meta.csv"),
  bilstm_csv = file.path(REAL_DIR, "bernardo_real_measles_predictions.csv")
)

cat("\nDone. Plots saved to real_data/plots/\n")
