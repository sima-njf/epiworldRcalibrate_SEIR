library(epiworldR)
library(data.table)
library(parallel)

ndays <- 365

run_sir_simulations <- function(N, ndays, ncores, theta, seeds,
                                disease_name = "Disease",
                                output_file_csv = "incidence_SIR_365.csv") {

  cat("==============================================\n")
  cat("Running", N, "SIR simulations on", ncores, "cores\n")
  cat("==============================================\n")

  incidence_list <- parallel::mclapply(1:N, FUN = function(i) {

    set.seed(seeds[i])

    # Create SIR model
    m <- epiworldR::ModelSIRCONN(
      name              = disease_name,
      n                 = theta$n[i],
      prevalence        = theta$preval[i],
      contact_rate      = theta$crate[i],
      transmission_rate = theta$ptran[i],
      recovery_rate     = theta$recov[i]
    )

    # Turn off verbose output
    epiworldR::verbose_off(m)

    # Run simulation
    epiworldR::run(m, ndays = ndays)

    # Extract incidence data
    incidence <- epiworldR::plot_incidence(m, plot = FALSE)
    incidence <- data.table::as.data.table(incidence)

    # Keep days 0 to ndays
    incidence <- incidence[as.integer(rownames(incidence)) <= (ndays + 1), ]

    # Create clean data.table
    ans <- data.table::data.table(
      infected  = incidence[["Infected"]],
      recovered = incidence[["Recovered"]]
    )

    # Fill missing values with last observed value
    nafill_cols <- c("infected", "recovered")

    for (col in nafill_cols) {
      ans[[col]] <- data.table::nafill(ans[[col]], type = "locf")
    }

    # Extract infected counts only
    # This removes day 0 and keeps days 1 to ndays
    infected_vector <- ans$infected[2:(ndays + 1)]

    return(infected_vector)

  }, mc.cores = ncores)

  cat("\n==============================================\n")
  cat("Combining results from all simulations...\n")
  cat("==============================================\n")

  # Each row = one simulation, each column = one day
  infected_matrix <- do.call(rbind, incidence_list)

  infected_dt <- data.table::as.data.table(infected_matrix)

  # Columns are days 1 to ndays
  colnames(infected_dt) <- paste0("V", 1:ndays)

  cat("\n==============================================\n")
  cat("Saving incidence data to:", output_file_csv, "\n")
  cat("==============================================\n")

  data.table::fwrite(infected_dt, file = output_file_csv)

  cat("\n==============================================\n")
  cat("Simulation complete!\n")
  cat("Total simulations:", N, "\n")
  cat("Days per simulation:", ndays, "\n")
  cat("Matrix dimensions:", dim(infected_dt), "\n")
  cat("==============================================\n")

  cat("\nInfected incidence preview:\n")
  print(head(infected_dt, 10))

  return(infected_dt)
}

# --------------------------
# 1. Load theta parameters from CSV
# --------------------------
theta_use <- data.table::fread("~/sima/epiworldRcalibrate_SEIR/SIR_epi_reconstruction/theta_use_sir.csv")

# Required columns:
# n, preval, crate, ptran, recov
# No incub column is needed for SIR.

# --------------------------
# 2. Set number of simulations
# --------------------------
N_SIMS <- nrow(theta_use)

# --------------------------
# 3. Generate seeds
# --------------------------
set.seed(122)
seeds <- sample(1:1000000, N_SIMS)

# --------------------------
# 4. Run simulations
# --------------------------
incidence_data <- run_sir_simulations(
  N               = N_SIMS,
  ndays           = ndays,
  ncores          = 11,
  theta           = theta_use,
  seeds           = seeds,
  disease_name    = "General Disease",
  output_file_csv = "~/sima/epiworldRcalibrate_SEIR/SIR_epi_reconstruction/incidence_SIR_365.csv"
)

# --------------------------
# 5. View results
# --------------------------
cat("\n==============================================\n")
cat("Final Results:\n")
cat("==============================================\n")

print(head(incidence_data))
print(dim(incidence_data))
