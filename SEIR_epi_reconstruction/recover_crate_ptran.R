# =============================================================================
#  recover_crate_ptran.R  — restore contact rate and transmission rate
#                           into test_actual_parameters.csv
#
#  WHY THIS IS NEEDED
#    test_actual_parameters.csv carries beta but not its factors. Since
#    beta = crate * ptran, the factors CANNOT be recovered from beta alone --
#    that is the same non-identifiability the calibration study is about.
#    They must be retrieved from the source file that generated beta.
#
#  WHERE THEY LIVE
#    The BiLSTM notebook builds the target as
#        theta_df['beta'] = theta_df['ptran'] * theta_df['crate']
#    so theta_use_seir.csv contains both columns.
#
#  WHAT THIS SCRIPT DOES
#    1. Locates theta_use_seir.csv
#    2. Works out how sim_idx maps to its rows (explicit id column, 0-based
#       row index, or 1-based row index) by testing each candidate and
#       checking whether the KNOWN columns (recov, incub, n, prevalence, beta)
#       agree after the merge. Agreement on those is the proof the mapping is
#       right -- it is not assumed.
#    3. Writes test_actual_parameters_full.csv with crate and ptran added.
#
#  Usage: Rscript recover_crate_ptran.R
# =============================================================================

PROJECT_DIR <- path.expand(
  "~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction"
)

# theta_use_seir.csv lives in the CHPC home directory per the notebook, but
# check a few plausible locations.
THETA_CANDIDATES <- c(
  file.path(PROJECT_DIR, "theta_use_seir.csv"),
  path.expand("~/theta_use_seir.csv"),
  path.expand("~/sima/theta_use_seir.csv"),
  path.expand("~/sima/epiworldRcalibrate_SEIR/theta_use_seir.csv")
)

TOL_REL <- 1e-4   # relative tolerance when checking column agreement

# =============================================================================
# Locate the files
# =============================================================================

actual_path <- file.path(PROJECT_DIR, "test_actual_parameters.csv")
if (!file.exists(actual_path)) stop("Not found: ", actual_path)

theta_path <- THETA_CANDIDATES[file.exists(THETA_CANDIDATES)][1]
if (is.na(theta_path)) {
  cat("Could not find theta_use_seir.csv in any of:\n")
  cat(paste0("  ", THETA_CANDIDATES, collapse = "\n"), "\n\n")
  cat("Locate it with:\n")
  cat("  find ~ -name 'theta_use_seir.csv' 2>/dev/null\n")
  cat("then set THETA_CANDIDATES at the top of this script.\n")
  quit(status = 1, save = "no")
}

actual <- read.csv(actual_path)
theta  <- read.csv(theta_path)

cat(sprintf("test set : %s  (%d rows)\n", basename(actual_path), nrow(actual)))
cat(sprintf("theta    : %s  (%d rows)\n", theta_path, nrow(theta)))
cat("\ntheta columns:\n"); print(names(theta))

# =============================================================================
# Find the crate / ptran columns
# =============================================================================

find_col <- function(df, cands) {
  hit <- intersect(tolower(cands), tolower(names(df)))
  if (length(hit) == 0) return(NA_character_)
  names(df)[tolower(names(df)) == hit[1]][1]
}

t_crate <- find_col(theta, c("crate", "contact_rate", "contactrate"))
t_ptran <- find_col(theta, c("ptran", "transmission_rate", "trate", "ptransmit"))

cat(sprintf("\ncrate column in theta : %s\n", ifelse(is.na(t_crate), "NOT FOUND", t_crate)))
cat(sprintf("ptran column in theta : %s\n", ifelse(is.na(t_ptran), "NOT FOUND", t_ptran)))

if (is.na(t_crate) || is.na(t_ptran)) {
  cat("\n*** theta_use_seir.csv does not contain both factors either.      ***\n")
  cat("*** Then only beta was ever stored, and crate/ptran cannot be     ***\n")
  cat("*** recovered from any existing file: beta = crate*ptran is not   ***\n")
  cat("*** invertible. You would need to re-run the data generation and  ***\n")
  cat("*** save both columns. Stick to single-parameter beta calibration ***\n")
  cat("*** until then.                                                   ***\n")
  quit(status = 1, save = "no")
}

# =============================================================================
# Work out the sim_idx -> theta row mapping, and PROVE it
# =============================================================================

# Columns present in both files that must agree if the mapping is correct
check_cols <- intersect(names(actual), names(theta))
check_cols <- setdiff(check_cols, "sim_idx")
cat(sprintf("\nColumns available to verify the mapping: %s\n",
            paste(check_cols, collapse = ", ")))
if (length(check_cols) == 0) {
  stop("No shared columns to verify a mapping against. Cannot proceed safely.")
}

agreement_score <- function(rows) {
  # Fraction of (row, column) comparisons that match within tolerance
  if (any(is.na(rows)) || any(rows < 1) || any(rows > nrow(theta))) return(-1)
  tot <- 0; hit <- 0
  for (cl in check_cols) {
    a <- as.numeric(actual[[cl]])
    b <- as.numeric(theta[[cl]][rows])
    good <- is.finite(a) & is.finite(b)
    tot <- tot + sum(good)
    hit <- hit + sum(abs(a[good] - b[good]) <= TOL_REL * (abs(a[good]) + 1e-9))
  }
  if (tot == 0) return(-1)
  hit / tot
}

candidates <- list()

t_id <- find_col(theta, c("sim_idx", "sim", "id", "index"))
if (!is.na(t_id)) {
  candidates[["explicit id column"]] <- match(actual$sim_idx, theta[[t_id]])
}
candidates[["sim_idx as 0-based row index"]] <- actual$sim_idx + 1L
candidates[["sim_idx as 1-based row index"]] <- actual$sim_idx

cat("\n========== Testing candidate mappings ==========\n")
scores <- vapply(candidates, agreement_score, numeric(1))
for (nm in names(candidates)) {
  cat(sprintf("  %-30s : %s\n", nm,
              ifelse(scores[[nm]] < 0, "invalid (out of range)",
                     sprintf("%.4f agreement", scores[[nm]]))))
}

best <- names(which.max(scores))
if (max(scores) < 0.999) {
  cat("\n*** No mapping reproduces the known columns. Do NOT proceed --   ***\n")
  cat("*** a wrong mapping would attach the wrong crate/ptran to every   ***\n")
  cat("*** simulation and silently corrupt the entire study.             ***\n")
  cat(sprintf("*** Best candidate was '%s' at %.4f agreement.\n", best, max(scores)))
  quit(status = 1, save = "no")
}

rows <- candidates[[best]]
cat(sprintf("\n=> Mapping: %s  (agreement %.6f)\n", best, max(scores)))

# =============================================================================
# Merge and validate
# =============================================================================

out <- actual
out$crate <- as.numeric(theta[[t_crate]][rows])
out$ptran <- as.numeric(theta[[t_ptran]][rows])

cat("\n========== Consistency checks ==========\n")

beta_rebuilt <- out$crate * out$ptran
gap_beta <- max(abs(out$beta - beta_rebuilt) / (abs(out$beta) + 1e-9))
cat(sprintf("max relative |stored beta - crate*ptran| : %.3e  %s\n",
            gap_beta, ifelse(gap_beta < 1e-4, "OK", "*** MISMATCH ***")))

if ("R0" %in% names(out)) {
  gap_R0 <- max(abs(out$R0 - beta_rebuilt / out$recov) / (abs(out$R0) + 1e-9))
  cat(sprintf("max relative |stored R0 - beta/recov|    : %.3e  %s\n",
              gap_R0, ifelse(gap_R0 < 1e-4, "OK", "*** MISMATCH ***")))
}

if (gap_beta >= 1e-4) {
  cat("\n*** crate*ptran does not reproduce the stored beta. The mapping   ***\n")
  cat("*** or the column choice is wrong. Not writing the output file.   ***\n")
  quit(status = 1, save = "no")
}

cat("\n========== Recovered ranges ==========\n")
rng <- function(x) sprintf("[%.5g, %.5g]  median %.5g", min(x), max(x), median(x))
cat(sprintf("  crate : %s\n", rng(out$crate)))
cat(sprintf("  ptran : %s\n", rng(out$ptran)))
cat(sprintf("  beta  : %s\n", rng(out$beta)))

# Is there anything to separate? If crate never varies, the two-parameter
# study collapses back to the one-parameter one.
if (diff(range(out$crate)) < 1e-9) {
  cat("\n*** crate is CONSTANT across the test set. There is no contact-   ***\n")
  cat("*** rate variation to recover: ptran is just beta/crate. Multi-    ***\n")
  cat("*** parameter calibration adds nothing here -- stay with beta.     ***\n")
} else {
  cat(sprintf("\ncrate varies %.1fx, ptran varies %.1fx, beta varies %.1fx.\n",
              max(out$crate) / min(out$crate),
              max(out$ptran) / min(out$ptran),
              max(out$beta)  / min(out$beta)))
  cat("If beta varies much less than its factors, the identifiability ridge\n")
  cat("is large in this test set and the two-parameter study is worthwhile.\n")
}

# ptran must be a valid probability for the simulator
bad <- sum(out$ptran <= 0 | out$ptran > 1)
cat(sprintf("\nptran outside (0, 1]: %d %s\n", bad,
            ifelse(bad == 0, "OK", "*** these sims cannot be simulated ***")))

# =============================================================================
# Write
# =============================================================================

out_path <- file.path(PROJECT_DIR, "test_actual_parameters_full.csv")
write.csv(out, out_path, row.names = FALSE)

cat(sprintf("\nSaved: %s\n", out_path))
cat("\nNext:\n")
cat("  1. Back up the original:\n")
cat("       cp test_actual_parameters.csv test_actual_parameters_beta_only.csv\n")
cat("  2. Swap the full version in:\n")
cat("       cp test_actual_parameters_full.csv test_actual_parameters.csv\n")
cat("  3. Rscript verify_parameter_ranges.R\n")
