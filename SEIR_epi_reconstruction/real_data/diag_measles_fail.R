suppressPackageStartupMessages(library(epiworldR))
setwd("~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction")
source("seir_common.R"); seir_set_libpath()

known <- list(n = 187, prevalence = 1/187, incub = 10, recov = 0.125)
beta <- 0.2175
crate <- max(beta, 1.0); ptran <- beta / crate
p <- list(n = known$n, prevalence = known$prevalence, incub = known$incub,
          recov = known$recov, crate = crate, ptran = ptran)
cat("p:", paste(names(p), unlist(p), sep="=", collapse="  "), "\n")

ci <- tryCatch(
  seir_run_multi_ci(p, ndays = 48, nreps = 20, nthreads = 1),
  error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(ci)) { cat("SUCCESS\n"); print(head(ci)) }
