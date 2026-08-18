suppressPackageStartupMessages(library(epiworldR))
setwd("~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction")
source("seir_common.R"); seir_set_libpath()

known <- list(n = 187, prevalence = 1/187, incub = 10, recov = 0.125)
mkp <- function(beta) {
  crate <- max(beta,1.0); ptran <- beta/crate
  list(n=known$n, prevalence=known$prevalence, incub=known$incub, recov=known$recov, crate=crate, ptran=ptran)
}
betas <- c(BiLSTM=0.2175, NelderMead=0.7891, DE=1.1621, `ABC-SMC`=2.3903, ABC=5.6365)
for (nm in names(betas)) {
  p <- mkp(betas[nm])
  totals <- vapply(1:300, function(i) sum(seir_run_single(p, ndays=48, seed=i)), numeric(1))
  extinct <- mean(totals <= 1)  # only the seed case itself, no onward spread
  cat(sprintf("%-12s beta=%.4f R0=%5.2f  frac_extinct(<=1 total case)=%.1f%%  median_total=%.0f  mean_total=%.1f\n",
      nm, betas[nm], betas[nm]/known$recov, extinct*100, median(totals), mean(totals)))
}
