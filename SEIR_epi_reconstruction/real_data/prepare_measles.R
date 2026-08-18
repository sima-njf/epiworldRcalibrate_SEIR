# prepare_measles.R
# Loads the Hagelloch 1861 measles outbreak from the 'outbreaks' R package
# and aggregates individual-level prodrome onset dates into daily incidence.
#
# Dataset: measles_hagelloch_1861
#   188 observed measles cases, Hagelloch, Germany, 1861
#
# For comparison with the SEIR model:
#   Model incidence = daily E -> I transitions
#   Observed incidence = daily date_of_prodrome counts
#
# date_of_prodrome is used as the empirical proxy for transition into
# the infectious state. date_of_rash occurs later and is therefore
# not used as the incidence date.
#
# Output (real_data/):
#   measles_hagelloch_incidence.csv
#   measles_hagelloch_meta.csv
#
# Usage:
#   Rscript real_data/prepare_measles.R

PROJECT_DIR <- path.expand(
  "~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction"
)

REAL_DIR <- file.path(PROJECT_DIR, "real_data")
dir.create(REAL_DIR, showWarnings = FALSE, recursive = TRUE)

# Install 'outbreaks' if needed
if (!requireNamespace("outbreaks", quietly = TRUE)) {
  cat("Installing 'outbreaks' package...\n")

  install.packages(
    "outbreaks",
    repos = "https://cloud.r-project.org",
    lib = "/uufs/chpc.utah.edu/common/home/u1418987/R/x86_64-pc-linux-gnu-library/4.4"
  )
}

suppressPackageStartupMessages(
  library(outbreaks)
)

cat("=== Preparing Hagelloch 1861 measles outbreak data ===\n\n")

# -------------------------------------------------------------------------
# Load data
# -------------------------------------------------------------------------

data("measles_hagelloch_1861", package = "outbreaks")
dat <- measles_hagelloch_1861

cat(sprintf("Loaded %d case records\n", nrow(dat)))
cat(sprintf(
  "Columns: %s\n\n",
  paste(names(dat), collapse = ", ")
))

# -------------------------------------------------------------------------
# Construct observed incidence
# -------------------------------------------------------------------------
# epiworldR SEIR incidence represents E -> I transitions.
#
# The Hagelloch dataset does not directly observe the mechanistic E -> I
# transition. We therefore use date_of_prodrome as the closest available
# clinical proxy for the beginning of the infectious/symptomatic phase.
#
# We DO NOT use date_of_rash because rash generally occurs later.

onset <- dat$date_of_prodrome

n_missing <- sum(is.na(onset))

if (n_missing > 0) {
  cat(sprintf(
    "Note: %d cases with missing prodrome date excluded.\n",
    n_missing
  ))
}

onset <- onset[!is.na(onset)]

# Ensure Date class
onset <- as.Date(onset)

start_date <- min(onset)
end_date   <- max(onset)

all_dates <- seq(
  from = start_date,
  to   = end_date,
  by   = "day"
)

daily_counts <- vapply(
  all_dates,
  function(d) sum(onset == d),
  integer(1L)
)

df <- data.frame(
  date        = all_dates,
  day_index   = seq_along(all_dates),
  daily_cases = daily_counts
)

cat(sprintf(
  "Observed period    : %s to %s (%d days)\n",
  start_date,
  end_date,
  nrow(df)
))

cat(sprintf(
  "Total cases        : %d\n",
  sum(df$daily_cases)
))

cat(sprintf(
  "Peak incidence     : %d cases on day %d (%s)\n",
  max(df$daily_cases),
  which.max(df$daily_cases),
  df$date[which.max(df$daily_cases)]
))

# -------------------------------------------------------------------------
# Save incidence curve
# -------------------------------------------------------------------------

incidence_file <- file.path(
  REAL_DIR,
  "measles_hagelloch_incidence.csv"
)

write.csv(
  df,
  incidence_file,
  row.names = FALSE
)

# -------------------------------------------------------------------------
# SEIR model assumptions
# -------------------------------------------------------------------------
#
# UPDATE 2026-08-17: n_pop=188 is NOT just a modeling shortcut -- it is the
# historically documented susceptible population. Per Neal & Roberts (2004,
# Biostatistics 5(2):249-261) and the `surveillance` package's hagelloch
# docs, Hagelloch had 577 total inhabitants in 1861, of whom exactly 188
# children aged <=15 were susceptible to measles (adults were presumably
# already immune from earlier outbreaks) -- and all 188 of them were
# infected over the course of the outbreak (100% attack rate among
# susceptibles). So n_pop=188 correctly reflects the closed susceptible
# pool; using the full village population (577) would overstate it.
#
# CONSEQUENCE for calibration: reproducing a ~100% attack rate THIS FAST
# (within the observed window, from a single seed case) in a homogeneously-
# mixed model (ModelSEIRCONN has no household/school structure) requires an
# inflated "effective R0" -- diagnosed 2026-08-17 on the companion EpiEstim
# smoke-test dataset (measles1861_epiestim), where ABC/NM/DE/ABC-SMC
# consistently converged to R0 in the 6-45 range regardless of the distance
# metric used (raw pointwise MSE vs shape-based summary statistics -- see
# distance_fn in calibrate_real_5method.R). This is not a calibration bug:
# Neal & Roberts' own paper exists because a simple well-mixed model can't
# explain this outbreak without added household/school-class contact
# structure, which ModelSEIRCONN doesn't have. Any R0 this pipeline reports
# for this dataset should be read as "the effective R0 a homogeneously-
# mixed model needs to match the outbreak's speed", not literally the
# biological R0 of measles.
#
# Fixed inputs used here:
#
#   n_pop      = 188
#       Modeled population size.
#
#   incub_days = 10
#       Assumed mean exposed/latent duration before transition E -> I.
#
#   recov_rate = 1 / 8 = 0.125 per day
#       Assumes a mean infectious duration of 8 days.
#
#   prevalence = 1 / 188
#       Initializes the model with approximately one infectious individual.
#
# These are fixed modeling assumptions, not parameters directly estimated
# from the Hagelloch line-list data.

N_POP <- 188L
INCUB_DAYS <- 10.0
INFECTIOUS_DAYS <- 8.0
RECOV_RATE <- 1.0 / INFECTIOUS_DAYS
INITIAL_INFECTED <- 1L
PREVALENCE <- INITIAL_INFECTED / N_POP

meta <- data.frame(
  dataset = "measles_hagelloch_1861",

  incidence_definition = "date_of_prodrome",
  model_incidence = "E_to_I",

  n_cases = nrow(dat),

  n_pop = N_POP,
  incub_days = INCUB_DAYS,
  infectious_days = INFECTIOUS_DAYS,
  recov_rate = RECOV_RATE,

  initial_infected = INITIAL_INFECTED,
  prevalence = PREVALENCE,

  n_days = nrow(df),
  peak_cases = max(df$daily_cases),
  total_cases = sum(df$daily_cases),

  note = paste(
    "Observed incidence is aggregated from date_of_prodrome",
    "and used as the empirical proxy for SEIR E->I incidence.",
    "n_pop=188 is the documented susceptible population (Neal & Roberts 2004:",
    "Hagelloch had 577 total inhabitants, 188 children <=15 susceptible to",
    "measles, all 188 infected), NOT the full village population.",
    "A homogeneously-mixed SEIR needs an inflated effective R0 to reproduce",
    "this outbreak's speed since it lacks the household/school contact",
    "structure Neal & Roberts show is needed to properly explain this data;",
    "R0 estimates from this pipeline should be read accordingly."
  ),

  stringsAsFactors = FALSE
)

meta_file <- file.path(
  REAL_DIR,
  "measles_hagelloch_meta.csv"
)

write.csv(
  meta,
  meta_file,
  row.names = FALSE
)

# -------------------------------------------------------------------------
# Summary
# -------------------------------------------------------------------------

cat("\nSaved:\n")
cat(sprintf("  %s\n", incidence_file))
cat(sprintf("  %s\n", meta_file))

cat("\nSEIR assumptions:\n")

cat(sprintf(
  "  Population size      : %d\n",
  N_POP
))

cat(sprintf(
  "  Initial infected     : %d\n",
  INITIAL_INFECTED
))

cat(sprintf(
  "  Initial prevalence   : %.6f\n",
  PREVALENCE
))

cat(sprintf(
  "  Latent period        : %.1f days\n",
  INCUB_DAYS
))

cat(sprintf(
  "  Infectious period    : %.1f days\n",
  INFECTIOUS_DAYS
))

cat(sprintf(
  "  Recovery rate        : %.3f / day\n",
  RECOV_RATE
))

cat("\nIncidence mapping:\n")
cat("  Real data : date_of_prodrome\n")
cat("  SEIR model: E -> I transitions\n")
