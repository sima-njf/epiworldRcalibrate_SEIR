# prepare_utah_covid.R
# Downloads Utah COVID-19 daily confirmed cases.
#
# Primary source : Utah COVID dashboard ZIP (discontinued ~2023)
# Fallback source: Johns Hopkins CSSE US time series (GitHub archive)
#
# Output (real_data/):
#   utah_covid_full.csv     — complete time series
#   utah_covid_wave1.csv    — first 365-day window starting at first confirmed case
#   utah_covid_meta.csv     — SEIR parameter priors for this dataset
#
# Usage:
#   Rscript real_data/prepare_utah_covid.R

suppressPackageStartupMessages(library(dplyr))

PROJECT_DIR <- path.expand(
  "~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction"
)
REAL_DIR <- file.path(PROJECT_DIR, "real_data")
dir.create(REAL_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== Preparing Utah COVID-19 incidence data ===\n\n")

df <- NULL

# ── Attempt A: Utah COVID dashboard ───────────────────────────────────────────
tryCatch({
  cat("Trying Utah COVID dashboard...\n")
  temp <- tempfile()
  download.file(
    "https://coronavirus-dashboard.utah.gov/Utah_COVID19_data.zip",
    temp, mode = "wb", quiet = TRUE
  )
  znames     <- unzip(temp, list = TRUE)$Name
  trend_file <- grep("Trends_Epidemic", znames,
                     ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(trend_file)) {
    raw      <- read.csv(unz(temp, trend_file))
    raw$Date <- as.Date(raw$Date)
    raw      <- raw %>% arrange(Date)
    df <- data.frame(
      date        = raw$Date,
      daily_cases = as.integer(pmax(raw$Daily.Cases, 0))
    )
    cat("  Success: Utah COVID dashboard.\n")
  }
  unlink(temp)
}, error = function(e) {
  cat("  Utah dashboard unavailable:", conditionMessage(e), "\n")
})

# ── Attempt B: Johns Hopkins CSSE archive (GitHub) ────────────────────────────
if (is.null(df)) {
  tryCatch({
    cat("Trying JHU CSSE archive (GitHub)...\n")
    url <- paste0(
      "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/",
      "csse_covid_19_data/csse_covid_19_time_series/",
      "time_series_covid19_confirmed_US.csv"
    )
    raw        <- read.csv(url, check.names = FALSE)
    utah       <- raw %>% filter(Province_State == "Utah")
    date_cols  <- names(utah)[12:ncol(utah)]
    dates      <- as.Date(date_cols, format = "%m/%d/%y")
    cumulative <- colSums(utah[, date_cols], na.rm = TRUE)
    daily      <- as.integer(pmax(c(cumulative[1], diff(cumulative)), 0))
    df <- data.frame(date = dates, daily_cases = daily)
    cat("  Success: JHU CSSE archive.\n")
  }, error = function(e) {
    stop("Both sources failed: ", conditionMessage(e), call. = FALSE)
  })
}

cat(sprintf("Full series: %d days, %s to %s\n",
    nrow(df), min(df$date), max(df$date)))

# ── Extract first-wave window (365 days from first confirmed case) ─────────────
first_case  <- min(df$date[df$daily_cases > 0], na.rm = TRUE)
wave_start  <- first_case
wave_end    <- wave_start + 364L

wave1 <- df %>%
  filter(date >= wave_start & date <= wave_end) %>%
  mutate(day_index = as.integer(date - wave_start + 1L))

# Pad to exactly 365 rows if the data ends early
if (nrow(wave1) < 365L) {
  pad_dates <- seq(max(wave1$date) + 1L, wave_start + 364L, by = "day")
  pad_df    <- data.frame(
    date        = pad_dates,
    daily_cases = 0L,
    day_index   = as.integer(pad_dates - wave_start + 1L)
  )
  wave1 <- rbind(wave1, pad_df)
}

cat(sprintf("Wave 1: %s to %s (%d days)\n", wave_start, wave_end, nrow(wave1)))
cat(sprintf("  Peak : %d cases on %s (day %d)\n",
    max(wave1$daily_cases),
    wave1$date[which.max(wave1$daily_cases)],
    which.max(wave1$daily_cases)))
cat(sprintf("  Total: %d cases\n", sum(wave1$daily_cases)))

# ── Save ──────────────────────────────────────────────────────────────────────
write.csv(df,    file.path(REAL_DIR, "utah_covid_full.csv"),  row.names = FALSE)
write.csv(wave1, file.path(REAL_DIR, "utah_covid_wave1.csv"), row.names = FALSE)

# SEIR parameter priors:
#   n          = 2020 Utah population (3.27 M)
#   incub_days = incubation period in days (COVID-19 ~ 5 days)
#   recov_rate = recovery rate per day (1/10 ~ 10-day infectious period)
#   prevalence = initial infected fraction (approx 1 seed at epidemic start)
meta <- data.frame(
  dataset      = "utah_covid_wave1",
  n_pop        = 3271616L,
  incub_days   = 5.0,
  recov_rate   = 0.1,
  prevalence   = 1.0 / 3271616,
  wave_start   = as.character(wave_start),
  wave_end     = as.character(wave_end),
  n_days       = nrow(wave1),
  peak_cases   = max(wave1$daily_cases),
  total_cases  = sum(wave1$daily_cases),
  note = paste(
    "Daily confirmed cases used as proxy for E->I incidence.",
    "Expect systematic under-count and ~5-day reporting lag."
  ),
  stringsAsFactors = FALSE
)
write.csv(meta, file.path(REAL_DIR, "utah_covid_meta.csv"), row.names = FALSE)

cat("\nSaved:\n")
cat("  real_data/utah_covid_full.csv\n")
cat("  real_data/utah_covid_wave1.csv\n")
cat("  real_data/utah_covid_meta.csv\n")
cat(sprintf("SEIR priors: n=%d, incub=%.0f days, recov=%.2f/day\n",
    meta$n_pop, meta$incub_days, meta$recov_rate))
