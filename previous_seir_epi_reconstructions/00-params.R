# Load required libraries
library(epiworldR)
library(data.table)
library(parallel)

# Optional plotting libraries
library(ggplot2)
library(dplyr)
library(tidyverse)
library(gridExtra)
library(cowplot)

# --------------------------
# Global Simulation Settings
# --------------------------
model_ndays <- 365    # simulation duration
model_seed  <- 122    # seed for reproducibility
N_SIMS      <- 20000  # number of simulations

# --------------------------
# Generate Parameter Sets using Theta (SEIR)
# --------------------------
set.seed(model_seed)

theta <- data.table(
  # N: population size, known
  n = sample(5000:10000, N_SIMS, replace = TRUE),

  # recovery rate, known
  recov = runif(N_SIMS, min = 0.071, max = 0.25),

  # contact rate
  crate = runif(N_SIMS, min = 5, max = 15),

  # incubation days for SEIR
  incub = runif(N_SIMS, min = 3, max = 21),

  # basic reproduction number
  R0 = runif(N_SIMS, min = 1, max = 5)
)

# Initial infected count: Uniform(100, 2000)
theta[, initial_infected := sample(100:2000, N_SIMS, replace = TRUE)]

# Prevalence = initial infected proportion
theta[, preval := initial_infected / n]

# Transmission rate:
# For SEIR/SIRCONN:
# R0 = ptran * crate / recov
# therefore ptran = R0 * recov / crate
theta[, ptran := R0 * recov / crate]

# Optional combined transmission force
theta[, beta := ptran * crate]

# Check:
# beta should equal R0 * recov
theta[, beta_check := R0 * recov]

# --------------------------
# Final dataset
# --------------------------
theta_use <- theta[, .(
  n,
  initial_infected,
  preval,
  crate,
  incub,
  recov,
  R0,
  ptran,
  beta
)]

# --------------------------
# Print summary
# --------------------------
cat("Theta generated successfully for SEIR.\n")
cat("\nDimensions:\n")
print(dim(theta_use))

cat("\nSample of parameter sets:\n")
print(head(theta_use, 10))

cat("\nSummary:\n")
print(summary(theta_use))

cat("\nCheck max absolute difference between beta and R0 * recov:\n")
print(max(abs(theta$beta - theta$beta_check)))

# --------------------------
# Save parameter sets
# --------------------------
data.table::fwrite(theta_use, "theta_use_seir.csv")
