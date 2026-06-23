library(epiworldR)
library(data.table)
library(parallel)
ndays=365
run_seir_simulations <- function(N, ndays, ncores, theta, seeds,
                                 disease_name = "Disease",
                                 output_file_csv = "incidence.csv") {

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
theta_use <- data.table::fread("theta_use_seir.csv")

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



# Save the parameter sets
data.table::fwrite(incidence_data, "incidence_365.csv")

# OR using base R:
write.csv(incidence_data, "incidence_365.csv", row.names = FALSE)
