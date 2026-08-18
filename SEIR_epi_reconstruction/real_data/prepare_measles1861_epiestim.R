# =============================================================================
#  real_data/prepare_measles1861_epiestim.R
#
#  Pipeline smoke-test dataset requested 2026-08-17: the ready-made daily
#  measles incidence bundled with the EpiEstim package (Measles1861), used
#  instead of the individual-level Hagelloch line-list (measles_hagelloch_*)
#  specifically to test that the SEIR real-data pipeline runs end-to-end
#  without the line-list aggregation step. No date cleaning / missing-day
#  imputation needed -- EpiEstim::Measles1861$incidence is already a clean
#  48-day daily vector.
#
#  NOTE: this is very likely the same underlying 1861 Hagelloch outbreak,
#  just pre-aggregated by EpiEstim from (probably) date-of-symptom-onset
#  rather than date-of-prodrome -- sum(incidence)=187 here vs 188 total
#  cases in measles_hagelloch_meta.csv. Per instructions, the E->I vs
#  symptom-onset distinction is deliberately ignored for this first test.
#
#  Usage:
#    Rscript real_data/prepare_measles1861_epiestim.R
# =============================================================================

suppressPackageStartupMessages(library(EpiEstim))

REAL_DIR <- path.expand("~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction/real_data")

data("Measles1861")
inc <- as.numeric(Measles1861$incidence)
ndays <- length(inc)

cat(sprintf("Measles1861 (EpiEstim): %d days, total cases=%d, peak=%d on day %d\n",
    ndays, sum(inc), max(inc), which.max(inc)))

# Incidence CSV -- day index only (no real calendar dates provided by EpiEstim
# for this dataset), matching process_dataset()'s optional "date" column.
inc_df <- data.frame(day_index = seq_len(ndays), daily_cases = inc)
write.csv(inc_df, file.path(REAL_DIR, "measles1861_epiestim.csv"), row.names = FALSE)

# Meta CSV -- same measles assumptions already used for measles_hagelloch_meta.csv
# (incubation ~10d, infectious period ~8d), so the two are directly comparable
# once the "scientifically correct" pass happens. n_pop = total observed cases,
# same modeling convention as the existing Hagelloch meta (a population-size
# stand-in, not a claim about the true susceptible pool).
meta <- data.frame(
  dataset      = "measles1861_epiestim",
  n_pop        = sum(inc),
  incub_days   = 10.0,
  infectious_days = 8.0,
  recov_rate   = 1 / 8,
  n_days       = ndays,
  peak_cases   = max(inc),
  total_cases  = sum(inc),
  note = paste(
    "EpiEstim::Measles1861$incidence, used as-is per instructions (no E->I",
    "vs symptom-onset correction yet). Smoke-test dataset for the SEIR",
    "real-data pipeline. n_pop = sum(incidence) is the documented",
    "susceptible population, not a proxy: per Neal & Roberts (2004,",
    "Biostatistics 5(2):249-261), Hagelloch (this is the same outbreak)",
    "had 577 total inhabitants but only 188 children <=15 susceptible to",
    "measles, all of whom were infected (100% attack rate among",
    "susceptibles). A homogeneously-mixed SEIR needs an inflated effective",
    "R0 to reproduce this outbreak's speed without household/school",
    "contact structure -- see calibrate_real_5method.R's distance_fn",
    "comment and prepare_measles.R for the full diagnosis (2026-08-17)."
  )
)
write.csv(meta, file.path(REAL_DIR, "measles1861_epiestim_meta.csv"), row.names = FALSE)

cat("Saved: measles1861_epiestim.csv, measles1861_epiestim_meta.csv\n")
