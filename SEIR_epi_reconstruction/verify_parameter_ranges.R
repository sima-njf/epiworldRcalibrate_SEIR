# =============================================================================
#  verify_parameter_ranges.R  — RUN THIS FIRST
#
#  Multi-parameter calibration only works if the bounds in seir_common.R
#  actually contain the true values. This prints the observed ranges and tells
#  you exactly what to set SEIR_BOUNDS to.
#
#  It also checks that the columns the pipeline needs are present:
#    crate (contact rate), ptran (transmission rate), recov, incub,
#    n, prevalence, and ideally beta and R0.
#
#  Usage: Rscript verify_parameter_ranges.R
# =============================================================================

PROJECT_DIR <- path.expand(
  "~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction"
)

source(file.path(PROJECT_DIR, "seir_common.R"))
seir_set_libpath()

actual <- read.csv(file.path(PROJECT_DIR, "test_actual_parameters.csv"))

cat("========== Columns present ==========\n")
print(names(actual))

# ---- Locate the contact-rate and transmission-rate columns -----------------

find_col <- function(df, cands) {
  hit <- intersect(tolower(cands), tolower(names(df)))
  if (length(hit) == 0) return(NA_character_)
  names(df)[tolower(names(df)) == hit[1]][1]
}

c_crate <- find_col(actual, c("crate", "contact_rate", "contactrate"))
c_ptran <- find_col(actual, c("ptran", "transmission_rate", "trate", "ptransmit"))
c_recov <- find_col(actual, c("recov", "recovery_rate", "gamma"))
c_incub <- find_col(actual, c("incub", "incubation_days"))
c_beta  <- find_col(actual, c("beta"))
c_R0    <- find_col(actual, c("r0"))

cat("\n========== Column mapping ==========\n")
for (nm in c("crate", "ptran", "recov", "incub", "beta", "R0")) {
  v <- get(paste0("c_", if (nm == "R0") "R0" else nm))
  cat(sprintf("  %-6s -> %s\n", nm, ifelse(is.na(v), "*** NOT FOUND ***", v)))
}

missing <- c("crate", "ptran")[is.na(c(c_crate, c_ptran))]
if (length(missing) > 0) {
  cat("\n*** test_actual_parameters.csv does not contain ", 
      paste(missing, collapse = " and "), ". ***\n", sep = "")
  cat("*** Multi-parameter calibration is impossible without the true    ***\n")
  cat("*** contact rate and transmission rate. Regenerate the test set   ***\n")
  cat("*** carrying these columns through from theta_use_seir.csv.       ***\n")
  quit(status = 1, save = "no")
}

# ---- Ranges ----------------------------------------------------------------

cat("\n========== Observed ranges ==========\n")
rng <- function(x) sprintf("[%.5g, %.5g]  median %.5g", min(x), max(x), median(x))
cat(sprintf("  crate : %s\n", rng(actual[[c_crate]])))
cat(sprintf("  ptran : %s\n", rng(actual[[c_ptran]])))
cat(sprintf("  recov : %s\n", rng(actual[[c_recov]])))
if (!is.na(c_incub)) cat(sprintf("  incub : %s\n", rng(actual[[c_incub]])))

beta_true <- actual[[c_crate]] * actual[[c_ptran]]
cat(sprintf("  beta = crate*ptran : %s\n", rng(beta_true)))
cat(sprintf("  R0   = beta/recov  : %s\n", rng(beta_true / actual[[c_recov]])))

# ---- Consistency with any stored beta / R0 ---------------------------------

cat("\n========== Consistency checks ==========\n")
if (!is.na(c_beta)) {
  gap <- max(abs(actual[[c_beta]] - beta_true))
  cat(sprintf("max |stored beta - crate*ptran| : %.3e %s\n", gap,
              ifelse(gap < 1e-6, "OK", "*** MISMATCH ***")))
}
if (!is.na(c_R0)) {
  gap <- max(abs(actual[[c_R0]] - beta_true / actual[[c_recov]]))
  cat(sprintf("max |stored R0 - beta/recov|    : %.3e %s\n", gap,
              ifelse(gap < 1e-6, "OK", "*** MISMATCH ***")))
}

# ---- Recommended bounds ----------------------------------------------------
# 10% padding so true values are never at the boundary, where optimisers stall
# and the uniform prior gives zero support.

pad <- function(x, frac = 0.10) {
  lo <- min(x); hi <- max(x); w <- hi - lo
  c(max(lo - frac * w, 0), hi + frac * w)
}

b_crate <- pad(actual[[c_crate]])
b_ptran <- pad(actual[[c_ptran]]); b_ptran[2] <- min(b_ptran[2], 1.0)
b_recov <- pad(actual[[c_recov]]); b_recov[2] <- min(b_recov[2], 1.0)
b_R0    <- pad(beta_true / actual[[c_recov]])

cat("\n========== Paste this into seir_common.R ==========\n")
cat("SEIR_BOUNDS <- list(\n")
cat(sprintf("  crate = c(%.5g, %.5g),\n", b_crate[1], b_crate[2]))
cat(sprintf("  ptran = c(%.5g, %.5g),\n", b_ptran[1], b_ptran[2]))
cat(sprintf("  recov = c(%.5g, %.5g),\n", b_recov[1], b_recov[2]))
cat(sprintf("  R0    = c(%.5g, %.5g)\n",  b_R0[1],    b_R0[2]))
cat(")\n")

cat("\n--- Parameterisation ---\n")
cat("  SEIR_PARS <- c(\"crate\", \"ptran\")   search the simulator's arguments\n")
cat("  SEIR_PARS <- c(\"crate\", \"R0\")      search R0 directly, ptran derived\n")
cat("Both fit the same model; they differ in which quantity carries the flat\n")
cat("prior. Use the R0 version if R0 is what you report.\n")

# Under (crate, R0), some in-box pairs imply ptran > 1 and are infeasible.
gr  <- expand.grid(crate = seq(b_crate[1], b_crate[2], length.out = 60),
                   R0    = seq(b_R0[1],    b_R0[2],    length.out = 60))
med_recov <- median(actual[[c_recov]])
frac_ok   <- mean(gr$R0 * med_recov / gr$crate <= 1)
cat(sprintf("\nAt the median recov (%.4g), %.1f%% of the (crate, R0) box implies a\n",
            med_recov, 100 * frac_ok))
cat("valid ptran <= 1. If that fraction is small, tighten the R0 upper bound\n")
cat("or raise the crate lower bound, or ABC will reject most proposals.\n")

cur <- seir_bounds_mat(c("crate", "ptran"))
inside <- all(actual[[c_crate]] >= cur$lo[1] & actual[[c_crate]] <= cur$hi[1]) &&
          all(actual[[c_ptran]] >= cur$lo[2] & actual[[c_ptran]] <= cur$hi[2])
cat(sprintf("\nCurrent bounds contain every true value: %s\n",
            ifelse(inside, "YES", "NO -- update before running anything")))

# ---- Identifiability preview -----------------------------------------------

cat("\n========== Identifiability note ==========\n")
cat("Mean SEIR dynamics depend on crate and ptran only through their product.\n")
cat(sprintf("Spread of the ridge in this test set: crate varies %.1fx, ptran %.1fx,\n",
            max(actual[[c_crate]]) / min(actual[[c_crate]]),
            max(actual[[c_ptran]]) / min(actual[[c_ptran]])))
cat(sprintf("while beta varies only %.1fx.\n", max(beta_true) / min(beta_true)))
cat("Expect beta to recover far better than crate or ptran separately.\n")
