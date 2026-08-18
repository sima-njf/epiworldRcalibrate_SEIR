suppressPackageStartupMessages(library(epiworldR))
setwd("~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction")
source("seir_common.R"); seir_set_libpath()

inc_df <- read.csv("real_data/measles1861_epiestim.csv")
obs <- as.numeric(inc_df$daily_cases)
known <- list(n = 187, prevalence = obs[1]/187, incub = 10, recov = 0.125)
obs_win <- obs; win_len <- length(obs_win); ndays <- win_len

distance_fn <- function(theta, known, obs_win, win_len, ndays) {
  p    <- seir_resolve(theta, known)
  pred <- seir_run_single(p, ndays = ndays)
  sqrt(mean((pred[seq_len(win_len)] - obs_win)^2))
}

SMC_NP <- 150L; SMC_NG <- 15L; SMC_EQ <- 0.3; SMC_BUD <- 300000L; ATT_MULT <- 300L

set.seed(1)
bnd <- seir_bounds_mat(); d <- length(SEIR_PARS)
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
cat("filled:", filled, "/", SMC_NP, "\n")
dist <- apply(theta,1,function(th) distance_fn(th,known,obs_win,win_len,ndays))
w    <- rep(1/SMC_NP,SMC_NP)
eps  <- as.numeric(quantile(dist,SMC_EQ,names=FALSE))
cat(sprintf("Gen 1: eps=%.3f  weighted median beta=%.4f  dist range=[%.2f,%.2f]\n",
    eps, wq_fn(theta[,1],w,0.5), min(dist), max(dist)))

n_ev <- 0L
for (g in 2:SMC_NG) {
  mu<-colSums(w*theta); dv<-sweep(theta,2,mu)
  S<-2*(t(dv)%*%(dv*w))+diag(1e-10*pmax(diag(2*(t(dv)%*%(dv*w))),1),d)
  chol_S<-tryCatch(chol(S),error=function(e) diag(sqrt(pmax(diag(S),1e-12)),d))
  S_inv<-solve(S); det_S<-det(S)
  new_theta<-matrix(NA_real_,SMC_NP,d,dimnames=list(NULL,SEIR_PARS))
  new_dist<-numeric(SMC_NP)
  i<-1L; atts<-0L; n_ev0<-n_ev
  t0<-proc.time()[["elapsed"]]
  while (i<=SMC_NP && atts<ATT_MULT*SMC_NP && n_ev<SMC_BUD) {
    atts<-atts+1L; idx<-sample.int(SMC_NP,1L,prob=w)
    cand<-rmvn_fn(theta[idx,],chol_S)
    if (!seir_feasible(cand,known)) next
    dd<-distance_fn(cand,known,obs_win,win_len,ndays); n_ev<-n_ev+1L
    if (dd<=eps) { new_theta[i,]<-cand; new_dist[i]<-dd; i<-i+1L }
  }
  el<-proc.time()[["elapsed"]]-t0
  if (i<=SMC_NP) { cat(sprintf("Gen %d: STALLED at i=%d/150, atts=%d, n_ev_this_gen=%d (%.1fs)\n", g, i-1, atts, n_ev-n_ev0, el)); break }
  dens<-vapply(seq_len(SMC_NP),function(k)
    sum(vapply(seq_len(SMC_NP),function(j)
      w[j]*dmvn_fn(new_theta[k,],theta[j,],S_inv,det_S),numeric(1))),numeric(1))
  w<-1/pmax(dens,.Machine$double.xmin); w<-w/sum(w)
  theta<-new_theta; dist<-new_dist
  eps<-as.numeric(quantile(dist,SMC_EQ,names=FALSE))
  cat(sprintf("Gen %d: n_ev_this_gen=%d (%.1fs) eps=%.4f  theta range=[%.3f,%.3f]  weighted median beta=%.4f  ESS=%.1f\n",
      g, n_ev-n_ev0, el, eps, min(theta), max(theta), wq_fn(theta[,1],w,0.5), 1/sum(w^2)))
}

# Also: single-realization noise check at a "reasonable" beta (R0~5, matching
# textbook lower end for measles) vs an "explosive" beta (R0~30, matching
# what ABC/NM found) -- run each 15 times to see distance-score variance.
cat("\n--- Noise check: distance_fn variance at two candidate betas ---\n")
for (beta_test in c(0.625, 3.8)) {  # R0 = beta/0.125 = 5, 30.4
  ds <- replicate(15, distance_fn(beta_test, known, obs_win, win_len, ndays))
  cat(sprintf("beta=%.3f (R0=%.1f): dist mean=%.2f sd=%.2f range=[%.2f,%.2f]\n",
      beta_test, beta_test/known$recov, mean(ds), sd(ds), min(ds), max(ds)))
}
