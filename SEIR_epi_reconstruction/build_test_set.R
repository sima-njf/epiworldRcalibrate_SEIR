# =============================================================================
#  build_test_set.R  — construct a fresh, self-consistent test set
#
#  WHY
#    test_actual_parameters.csv carried beta but not crate/ptran, and the
#    home-directory theta_use_seir.csv is a DIFFERENT (older, 7-column) file
#    that does not align with it. Rather than guess a mapping, this script
#    rebuilds the test set directly from the source pair in data_construction:
#
#       theta_use_seir.csv      20000 rows: n, initial_infected, preval,
#                               crate, incub, recov, R0, ptran, beta
#       incidence_seir_365.csv  20000 rows x 365 days
#
#    Row i of one corresponds to row i of the other by construction, so no
#    mapping has to be inferred.
#
#  HELD-OUT SPLIT
#    The BiLSTM was trained on 80% of these sims
#    (train_test_split(np.arange(N_SIMS), test_size=0.2, random_state=42)).
#    Sampling from all 20000 would evaluate the network on its own training
#    data. This script therefore prefers the sim_idx list already present in
#    the existing test_actual_parameters.csv (the validation split) and only
#    falls back to the full pool with an explicit warning.
#
#  OUTPUTS (old files are backed up first)
#    test_actual_parameters.csv   sampled parameters, incl. crate and ptran
#    test_incidence_raw.npy       matching incidence rows
#    test_incidence_raw.csv       same, as a fallback
#
#  Usage: Rscript build_test_set.R
# =============================================================================

PROJECT_DIR <- path.expand(
  "~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction"
)
DATA_DIR <- file.path(PROJECT_DIR, "data_construction")

THETA_CSV     <- file.path(DATA_DIR, "theta_use_seir.csv")
INCIDENCE_CSV <- file.path(DATA_DIR, "incidence_seir_365.csv")

N_SAMPLE   <- 1000L
SAMPLE_SEED <- 2024L
NDAYS      <- 365L
N_INC_CHECK <- 6L      # sims used for the incidence-definition check
N_INC_REPS  <- 25L     # replicates per sim in that check

source(file.path(PROJECT_DIR, "seir_common.R"))
seir_set_libpath()

suppressPackageStartupMessages({
  library(data.table)
  library(reticulate)
  library(epiworldR)
})

# =============================================================================
# 1. Load
# =============================================================================

for (f in c(THETA_CSV, INCIDENCE_CSV)) {
  if (!file.exists(f)) stop("Not found: ", f)
}

theta <- as.data.frame(fread(THETA_CSV))
cat(sprintf("theta     : %d rows x %d cols\n", nrow(theta), ncol(theta)))
cat("  columns : "); cat(paste(names(theta), collapse = ", "), "\n")

inc <- as.matrix(fread(INCIDENCE_CSV))
cat(sprintf("incidence : %d rows x %d cols\n", nrow(inc), ncol(inc)))

# Drop a leading row-index column if one snuck in
if (ncol(inc) == NDAYS + 1L) {
  first <- inc[, 1]
  if (all(abs(first - seq_len(nrow(inc))) < 1e-9) ||
      all(abs(first - (seq_len(nrow(inc)) - 1L)) < 1e-9)) {
    cat("  dropping leading row-index column\n")
    inc <- inc[, -1, drop = FALSE]
  }
}

if (ncol(inc) != NDAYS) {
  stop(sprintf("incidence has %d columns, expected %d", ncol(inc), NDAYS))
}
if (nrow(inc) != nrow(theta)) {
  stop(sprintf("row mismatch: theta %d, incidence %d", nrow(theta), nrow(inc)))
}

# Column aliases
col_of <- function(cands) {
  hit <- intersect(tolower(cands), tolower(names(theta)))
  if (length(hit) == 0) return(NA_character_)
  names(theta)[tolower(names(theta)) == hit[1]][1]
}

c_n     <- col_of(c("n"))
c_prev  <- col_of(c("preval", "prevalence", "prev"))
c_crate <- col_of(c("crate", "contact_rate"))
c_ptran <- col_of(c("ptran", "transmission_rate"))
c_recov <- col_of(c("recov", "recovery_rate"))
c_incub <- col_of(c("incub", "incubation_days"))
c_beta  <- col_of(c("beta"))
c_R0    <- col_of(c("r0"))
c_init  <- col_of(c("initial_infected", "init_infected"))

need <- c(c_n, c_prev, c_crate, c_ptran, c_recov, c_incub)
if (any(is.na(need))) stop("theta is missing one of n/preval/crate/ptran/recov/incub")

# =============================================================================
# 2. Structural checks on the source data
# =============================================================================

cat("\n========== Source data checks ==========\n")

beta_calc <- theta[[c_crate]] * theta[[c_ptran]]
if (!is.na(c_beta)) {
  g <- max(abs(theta[[c_beta]] - beta_calc) / (abs(theta[[c_beta]]) + 1e-12))
  cat(sprintf("beta == crate*ptran   : max rel gap %.3e  %s\n", g,
              ifelse(g < 1e-4, "OK", "*** MISMATCH ***")))
}
if (!is.na(c_R0)) {
  R0_calc <- beta_calc / theta[[c_recov]]
  g <- max(abs(theta[[c_R0]] - R0_calc) / (abs(theta[[c_R0]]) + 1e-12))
  cat(sprintf("R0 == beta/recov      : max rel gap %.3e  %s\n", g,
              ifelse(g < 1e-4, "OK", "*** MISMATCH ***")))
  if (g >= 1e-4) {
    stop("R0 is not beta/recov in the source data. seir_derive() in ",
         "seir_common.R assumes it is. Resolve before continuing.")
  }
}
if (!is.na(c_init)) {
  g <- max(abs(theta[[c_prev]] - theta[[c_init]] / theta[[c_n]]))
  cat(sprintf("preval == init/n      : max abs gap %.3e  %s\n", g,
              ifelse(g < 1e-6, "OK", "note: prevalence defined differently")))
}
cat(sprintf("ptran within (0, 1]   : %d violations\n",
            sum(theta[[c_ptran]] <= 0 | theta[[c_ptran]] > 1)))

# =============================================================================
# 3. Choose the sampling pool (prefer the held-out split)
# =============================================================================

cat("\n========== Sampling pool ==========\n")

old_actual_path <- file.path(PROJECT_DIR, "test_actual_parameters.csv")
pool <- NULL

if (file.exists(old_actual_path)) {
  old <- read.csv(old_actual_path)
  if ("sim_idx" %in% names(old)) {
    rows0 <- old$sim_idx + 1L      # sim_idx is 0-based (Python)
    if (all(rows0 >= 1 & rows0 <= nrow(theta))) {
      # Verify against columns copied verbatim between the two files
      chk <- c(n = c_n, recov = c_recov, incub = c_incub)
      frac <- vapply(names(chk), function(k) {
        if (!k %in% names(old)) return(NA_real_)
        a <- as.numeric(old[[k]]); b <- as.numeric(theta[[chk[k]]][rows0])
        mean(abs(a - b) <= 1e-4 * (abs(a) + 1e-9))
      }, numeric(1))
      cat("Existing test set vs data_construction theta, per column:\n")
      for (k in names(frac)) cat(sprintf("  %-6s : %.4f\n", k, frac[k]))

      if (all(frac[!is.na(frac)] > 0.999)) {
        pool <- old$sim_idx
        cat(sprintf("=> Verified. Sampling from the %d held-out sims.\n", length(pool)))
      }
    }
  }
}

if (is.null(pool)) {
  cat("*** Could not verify the existing held-out split against this theta ***\n")
  cat("*** file. Falling back to the full pool of 20000 sims.              ***\n")
  cat("*** WARNING: ~80%% of those were BiLSTM TRAINING data. Any BiLSTM   ***\n")
  cat("*** number computed on this test set is optimistically biased and   ***\n")
  cat("*** must not be reported. ABC/ABC-SMC/NM/DE are unaffected: they    ***\n")
  cat("*** never saw the data at all.                                      ***\n")
  pool <- seq_len(nrow(theta)) - 1L
}

# =============================================================================
# 4. Sample
# =============================================================================

set.seed(SAMPLE_SEED)
sim_idx <- sort(sample(pool, min(N_SAMPLE, length(pool))))
rows    <- sim_idx + 1L

cat(sprintf("\nSampled %d sims (seed %d)\n", length(sim_idx), SAMPLE_SEED))

test_params <- data.frame(
  sim_idx    = sim_idx,
  n          = theta[[c_n]][rows],
  prevalence = theta[[c_prev]][rows],
  crate      = theta[[c_crate]][rows],
  ptran      = theta[[c_ptran]][rows],
  recov      = theta[[c_recov]][rows],
  incub      = theta[[c_incub]][rows],
  stringsAsFactors = FALSE
)
test_params$beta <- test_params$crate * test_params$ptran
test_params$R0   <- test_params$beta / test_params$recov

inc_sample <- inc[rows, , drop = FALSE]
storage.mode(inc_sample) <- "double"

# =============================================================================
# 5. Incidence-definition check (S->E vs E->I)
# =============================================================================

cat("\n========== Incidence definition check ==========\n")
cat("Simulating with the true parameters and comparing both conventions.\n")

both_conventions <- function(model, ndays) {
  tm <- as.data.frame(get_hist_transition_matrix(model))
  cc <- seir_transition_cols(tm)
  fr <- as.character(tm[[cc$from]]); to <- as.character(tm[[cc$to]])
  dt <- as.integer(tm[[cc$date]]);  ct <- as.numeric(tm[[cc$counts]])
  grab <- function(fp, tp) {
    k <- grepl(fp, fr, ignore.case = TRUE) & grepl(tp, to, ignore.case = TRUE) &
         dt >= 1L & dt <= ndays
    o <- numeric(ndays)
    if (any(k)) { a <- tapply(ct[k], dt[k], sum); o[as.integer(names(a))] <- as.numeric(a) }
    o
  }
  list(s2e = grab("^suscep", "^expos"), e2i = grab("^expos", "^infect"))
}

set.seed(7)
chk <- sample(seq_len(nrow(test_params)), min(N_INC_CHECK, nrow(test_params)))
err_s2e <- err_e2i <- cor_s2e <- cor_e2i <- numeric(length(chk))

for (k in seq_along(chk)) {
  s   <- test_params[chk[k], ]
  obs <- as.numeric(inc_sample[chk[k], ])
  p   <- list(n = s$n, prevalence = s$prevalence, crate = s$crate,
              ptran = s$ptran, recov = s$recov, incub = s$incub)

  m_s2e <- m_e2i <- matrix(0, N_INC_REPS, NDAYS)
  for (r in seq_len(N_INC_REPS)) {
    mo <- seir_build(p)
    run(mo, ndays = NDAYS + 1L, seed = as.integer(1000 * k + r))
    bb <- both_conventions(mo, NDAYS)
    m_s2e[r, ] <- bb$s2e; m_e2i[r, ] <- bb$e2i
  }
  med_s2e <- apply(m_s2e, 2, median); med_e2i <- apply(m_e2i, 2, median)

  err_s2e[k] <- abs(which.max(med_s2e) - which.max(obs))
  err_e2i[k] <- abs(which.max(med_e2i) - which.max(obs))
  cor_s2e[k] <- seir_cor(obs, med_s2e)
  cor_e2i[k] <- seir_cor(obs, med_e2i)

  cat(sprintf("  sim %d: obs peak %d | S->E %d | E->I %d\n",
              s$sim_idx, which.max(obs), which.max(med_s2e), which.max(med_e2i)))
}

cat(sprintf("\nMean |peak-day error|  S->E %.2f  |  E->I %.2f\n",
            mean(err_s2e), mean(err_e2i)))
cat(sprintf("Mean correlation       S->E %.4f |  E->I %.4f\n",
            mean(cor_s2e, na.rm = TRUE), mean(cor_e2i, na.rm = TRUE)))

better <- if (mean(err_e2i) <= mean(err_s2e)) "E->I" else "S->E"
cat(sprintf("=> incidence_seir_365.csv looks like %s\n", better))
if (better != "E->I") {
  cat("*** The pipeline extracts E->I. Change seir_incidence_e2i() in    ***\n")
  cat("*** seir_common.R to use Susceptible->Exposed, or every metric    ***\n")
  cat("*** will compare two different quantities.                        ***\n")
}

# =============================================================================
# 6. Write (backing up whatever is there)
# =============================================================================

cat("\n========== Writing ==========\n")

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
for (f in c("test_actual_parameters.csv", "test_incidence_raw.npy")) {
  p <- file.path(PROJECT_DIR, f)
  if (file.exists(p)) {
    bk <- file.path(PROJECT_DIR, sprintf("%s.bak_%s", f, stamp))
    file.copy(p, bk, overwrite = TRUE)
    cat(sprintf("  backed up %s -> %s\n", f, basename(bk)))
  }
}

write.csv(test_params, file.path(PROJECT_DIR, "test_actual_parameters.csv"),
          row.names = FALSE)
cat("  wrote test_actual_parameters.csv\n")

np <- reticulate::import("numpy", convert = FALSE)
np$save(file.path(PROJECT_DIR, "test_incidence_raw.npy"), inc_sample)
cat("  wrote test_incidence_raw.npy\n")

fwrite(as.data.frame(inc_sample),
       file.path(PROJECT_DIR, "test_incidence_raw.csv"), col.names = FALSE)
cat("  wrote test_incidence_raw.csv (fallback)\n")

# Read the npy back and confirm it round-trips
chk_npy <- as.matrix(np$load(file.path(PROJECT_DIR, "test_incidence_raw.npy")))
cat(sprintf("  npy round-trip: %d x %d, max diff %.3e\n",
            nrow(chk_npy), ncol(chk_npy), max(abs(chk_npy - inc_sample))))

# =============================================================================
# 7. Ranges and bounds
# =============================================================================

rng <- function(x) sprintf("[%.5g, %.5g]  median %.5g", min(x), max(x), median(x))

cat("\n========== Test-set ranges ==========\n")
for (nm in c("n", "prevalence", "crate", "ptran", "recov", "incub", "beta", "R0")) {
  cat(sprintf("  %-10s : %s\n", nm, rng(test_params[[nm]])))
}

cat(sprintf("\ncrate varies %.1fx | ptran varies %.1fx | beta varies %.1fx\n",
            max(test_params$crate)/min(test_params$crate),
            max(test_params$ptran)/min(test_params$ptran),
            max(test_params$beta)/min(test_params$beta)))

pad <- function(x, f = 0.10) {
  lo <- min(x); hi <- max(x); w <- hi - lo
  c(max(lo - f * w, 0), hi + f * w)
}
b_crate <- pad(test_params$crate)
b_ptran <- pad(test_params$ptran); b_ptran[2] <- min(b_ptran[2], 1.0)
b_recov <- pad(test_params$recov); b_recov[2] <- min(b_recov[2], 1.0)
b_R0    <- pad(test_params$R0)

cat("\n========== Paste into seir_common.R ==========\n")
cat("SEIR_BOUNDS <- list(\n")
cat(sprintf("  crate = c(%.5g, %.5g),\n", b_crate[1], b_crate[2]))
cat(sprintf("  ptran = c(%.5g, %.5g),\n", b_ptran[1], b_ptran[2]))
cat(sprintf("  recov = c(%.5g, %.5g),\n", b_recov[1], b_recov[2]))
cat(sprintf("  R0    = c(%.5g, %.5g)\n",  b_R0[1],    b_R0[2]))
cat(")\n")

cat("\nNext:\n")
cat("  1. paste the bounds above into seir_common.R\n")
cat("  2. Rscript verify_parameter_ranges.R\n")
cat("  3. Rscript abc_seir_submit.R / abcsmc_seir_submit.R / optimize_seir_submit.R\n")
cat("  4. ... --collect each, then Rscript validate_methods.R\n")
cat("\nNOTE: any existing test_bilstm_predictions_tuned.csv refers to the OLD\n")
cat("test set. Regenerate BiLSTM predictions for these sim_idx values, or\n")
cat("validate_methods.R will simply omit BiLSTM.\n")
