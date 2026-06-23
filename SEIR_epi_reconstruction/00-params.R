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
data.table::fwrite(theta_use, "~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction/theta_use_seir.csv")



library(epiworldR)
library(data.table)
library(parallel)
ndays=365
run_seir_simulations <- function(N, ndays, ncores, theta, seeds,
                                 disease_name = "Disease",
                                 output_file_csv = "incidence_seir.csv") {

  cat("==============================================\n")
  cat("Running", N, "SEIR simulations on", ncores, "cores\n")
  cat("==============================================\n")

  # Run simulations in parallel
  incidence_list <- parallel::mclapply(1:N, FUN = function(i) {
    set.seed(seeds[i])
    # Create SEIR model
    m <- epiworldR::ModelSEIRCONN(
      name              = disease_name,
      n                 = theta$n[i],
      prevalence        = theta$preval[i],
      contact_rate      = theta$crate[i],
      incubation_days   = theta$incub[i],
      transmission_rate = theta$ptran[i],
      recovery_rate     = theta$recov[i]
    )

    # Turn off verbose output
    epiworldR::verbose_off(m)

    # Run the simulation
    epiworldR::run(m, ndays = ndays)

    # Extract incidence data
    ans <- list(

      incidence = epiworldR::plot_incidence(m, plot = FALSE)

    )

    # Filling
    ans <- lapply(ans, data.table::as.data.table)

    # Replacing NaN and NAs with the previous value
    # in each element in the list
    # Replace NA values with the last observed value
    ans$incidence <- ans$incidence[as.integer(rownames(ans$incidence)) <= (ndays + 1),]

    # Reference table for merging
    # ndays <- epiworldR::get_ndays(m)

    ref_table <- data.table::data.table(
      date = 0:ndays
    )

    # Replace the $ with the [[ ]] to avoid the warning in the next
    # two lines



    # Generating the data.table with necessary columns
    ans <- data.table::data.table(
      infected  = ans[["incidence"]][["Infected"]],
      recovered = ans[["incidence"]][["Recovered"]]
    )

    # Replace NA values with the last observed value for all columns
    nafill_cols <- c("infected", "recovered")

    for (col in nafill_cols) {
      ans[[col]] <- data.table::nafill(ans[[col]], type = "locf")
    }
    # Extract only infected counts (as a vector)
    infected_vector <- ans$infected[2:(ndays+1)]
    return(infected_vector)

  }, mc.cores = ncores)

  # Combine all simulations
  cat("\n==============================================\n")
  cat("Combining results from all simulations...\n")
  cat("==============================================\n")

  # Convert list to matrix (each row = one simulation, each column = one day)
  infected_matrix <- do.call(rbind, incidence_list)

  # Convert to data.table
  infected_dt <- data.table::as.data.table(infected_matrix)

  # Rename columns to V1, V2, V3, ... V(ndays+1)
  colnames(infected_dt) <- paste0("V", 2:(ndays + 1))

  # Save as CSV
  cat("\n==============================================\n")
  cat("Saving incidence data to:", output_file_csv, "\n")
  cat("==============================================\n")

  data.table::fwrite(infected_dt, file = output_file_csv)

  cat("\n==============================================\n")
  cat("Simulation complete!\n")
  cat("Total simulations:", N, "\n")
  cat("Days per simulation:", ndays + 1, "\n")
  cat("Matrix dimensions:", dim(infected_dt), "\n")
  cat("==============================================\n")

  cat("\nInfected incidence preview:\n")
  print(head(infected_dt, 10))

  return(infected_dt)
}

# --------------------------
# 1. Load theta parameters from CSV
# --------------------------
theta_use <- data.table::fread("~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction/theta_use_seir.csv")

# --------------------------
# 2. Set number of simulations
# --------------------------
N_SIMS <- nrow(theta_use)  # Use all rows from theta.csv

# --------------------------
# 3. Generate seeds
# --------------------------
set.seed(122)
seeds <- sample(1:1000000, N_SIMS)

# --------------------------
# 4. Run simulations
# --------------------------
incidence_data <- run_seir_simulations(
  N            = N_SIMS,
  ndays        = ndays,  # This will create V1 to V61 (days 0-60)
  ncores       = 11,
  theta        = theta_use,
  seeds        = seeds,
  disease_name = "General Disease",
  output_file_csv = "incidence.csv"
)

# --------------------------
# 5. View results
# --------------------------
cat("\n==============================================\n")
cat("Final Results:\n")
cat("==============================================\n")

print(head(incidence_data))
print(dim(incidence_data))




# OR using base R:
write.csv(incidence_data, "~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction/incidence_seir.csv", row.names = FALSE)

