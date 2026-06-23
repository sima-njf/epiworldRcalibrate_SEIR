# Load required libraries
library(epiworldR)
library(data.table)
library(parallel)

# --------------------------
# Global Simulation Settings
# --------------------------
model_ndays <- 365
model_seed  <- 122
N_SIMS      <- 20000

# --------------------------
# Generate Parameter Sets
# --------------------------
set.seed(model_seed)

theta <- data.table(
  # N: population size, known
  n = sample(5000:10000, N_SIMS, replace = TRUE),

  # recovery rate, known
  recov = runif(N_SIMS, min = 0.071, max = 0.25),

  # contact rate
  crate = runif(N_SIMS, min = 5, max = 15),

  # R0
  R0 = runif(N_SIMS, min = 1, max = 5)
)

# Initial infected count: Uniform(100, 2000)
theta[, initial_infected := sample(100:2000, N_SIMS, replace = TRUE)]

# Prevalence = initial infected proportion
theta[, preval := initial_infected / n]

# Transmission rate:
# R0 = ptran * crate / recov
# therefore ptran = R0 * recov / crate
theta[, ptran := R0 * recov / crate]

# Optional: combined transmission force beta = ptran * crate
theta[, beta := ptran * crate]

# Check relationship:
# beta should equal R0 * recov
theta[, beta_check := R0 * recov]

# --------------------------
# Final dataset
# --------------------------
theta_use <- theta[, .(
  n,
  preval,
  crate,
  recov,
  R0,
  ptran
)]

# --------------------------
# Save parameter sets
# --------------------------
data.table::fwrite(theta_use, "~/sima/epiworldRcalibrate_SEIR/SIR_epi_reconstruction/theta_use_sir.csv")

# --------------------------
# Quick checks
# --------------------------
cat("Theta generated successfully.\n")
cat("Dimensions:\n")
print(dim(theta_use))

cat("\nPreview:\n")
print(head(theta_use))

cat("\nSummary:\n")
print(summary(theta_use))

cat("\nCheck max absolute difference between beta and R0 * recov:\n")
print(max(abs(theta$beta - theta$beta_check)))
