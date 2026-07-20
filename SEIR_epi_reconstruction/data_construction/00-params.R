#  Purpose:
#    1. Generate SEIR parameter sets.
#    2. Save the parameter table.
#    3. Run all SEIR simulations in parallel.
#    4. Save the daily incidence matrix.
#
#  Expected location:
#    ~/DeepIMC/data_construction
#
#  Outputs:
#    ~/DeepIMC/data_construction/theta_use_seir.csv
#    ~/DeepIMC/data_construction/incidence_seir.csv
# =============================================================================

# =============================================================================
# Required libraries
# =============================================================================

required_packages <- c(
  "epiworldR",
  "data.table",
  "parallel"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Please install the following packages before running this script:\n",
      paste0("  - ", missing_packages, collapse = "\n")
    ),
    call. = FALSE
  )
}

library(epiworldR)
library(data.table)
library(parallel)

# =============================================================================
# Portable project paths
# =============================================================================

# Every user should place or clone the repository directly under their home:
#   ~/DeepIMC
PROJECT_DIR <- path.expand("~/DeepIMC")

# All files created by this script are stored here:
DATA_CONSTRUCTION_DIR <- file.path(
  PROJECT_DIR,
  "data_construction"
)

dir.create(
  DATA_CONSTRUCTION_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

THETA_FILE <- file.path(
  DATA_CONSTRUCTION_DIR,
  "theta_use_seir.csv"
)

INCIDENCE_FILE <- file.path(
  DATA_CONSTRUCTION_DIR,
  "incidence_seir.csv"
)

# =============================================================================
# Global simulation settings
# =============================================================================

MODEL_NDAYS <- 365L
MODEL_SEED  <- 122L
N_SIMS      <- 20000L
N_CORES     <- 11L

# Avoid requesting more cores than are available.
available_cores <- parallel::detectCores(logical = TRUE)

if (!is.na(available_cores)) {
  N_CORES <- min(
    N_CORES,
    max(1L, available_cores - 1L)
  )
}

cat("Project directory:\n  ", PROJECT_DIR, "\n", sep = "")
cat("Data-construction directory:\n  ", DATA_CONSTRUCTION_DIR, "\n", sep = "")
cat("Number of simulations: ", N_SIMS, "\n", sep = "")
cat("Number of days: ", MODEL_NDAYS, "\n", sep = "")
cat("Number of cores: ", N_CORES, "\n\n", sep = "")

# =============================================================================
# Part 1: Generate SEIR parameter sets
# =============================================================================

set.seed(MODEL_SEED)

theta <- data.table(
  # Population size: known input
  n = sample(
    5000:10000,
    N_SIMS,
    replace = TRUE
  ),

  # Recovery rate: known input
  recov = runif(
    N_SIMS,
    min = 0.071,
    max = 0.25
  ),

  # Contact rate
  crate = runif(
    N_SIMS,
    min = 5,
    max = 15
  ),

  # Incubation period in days: known input
  incub = runif(
    N_SIMS,
    min = 3,
    max = 21
  ),

  # Basic reproduction number
  R0 = runif(
    N_SIMS,
    min = 1,
    max = 5
  )
)

# Initial infected count: Uniform(100, 2000)
theta[
  ,
  initial_infected := sample(
    100:2000,
    N_SIMS,
    replace = TRUE
  )
]

# Initial prevalence
theta[
  ,
  preval := initial_infected / n
]

# For ModelSEIRCONN:
#   R0 = ptran * crate / recov
# Therefore:
#   ptran = R0 * recov / crate
theta[
  ,
  ptran := R0 * recov / crate
]

# Combined transmission force
theta[
  ,
  beta := ptran * crate
]

# Numerical consistency check
theta[
  ,
  beta_check := R0 * recov
]

theta_use <- theta[
  ,
  .(
    n,
    initial_infected,
    preval,
    crate,
    incub,
    recov,
    R0,
    ptran,
    beta
  )
]

cat("Theta generated successfully for SEIR.\n\n")

cat("Dimensions:\n")
print(dim(theta_use))

cat("\nFirst 10 parameter sets:\n")
print(head(theta_use, 10))

cat("\nParameter summary:\n")
print(summary(theta_use))

max_beta_difference <- max(
  abs(theta$beta - theta$beta_check)
)

cat(
  "\nMaximum absolute difference between beta and R0 * recov: ",
  max_beta_difference,
  "\n",
  sep = ""
)

data.table::fwrite(
  theta_use,
  file = THETA_FILE
)

cat("\nSaved parameter sets to:\n  ", THETA_FILE, "\n", sep = "")

# =============================================================================
# Helper: Run one SEIR simulation
# =============================================================================

run_one_seir_simulation <- function(
    i,
    theta,
    seeds,
    ndays,
    disease_name
) {

  n_i <- max(
    as.integer(round(theta$n[i])),
    10L
  )

  prevalence_i <- max(
    min(theta$preval[i], 1),
    1 / n_i
  )

  model <- epiworldR::ModelSEIRCONN(
    name              = disease_name,
    n                 = n_i,
    prevalence        = prevalence_i,
    contact_rate      = theta$crate[i],
    incubation_days   = theta$incub[i],
    transmission_rate = theta$ptran[i],
    recovery_rate     = theta$recov[i]
  )

  epiworldR::verbose_off(model)

  epiworldR::run(
    model,
    ndays = ndays + 1L,
    seed = as.integer(seeds[i])
  )

  incidence_table <- data.table::as.data.table(
    epiworldR::plot_incidence(
      model,
      plot = FALSE
    )
  )

  if (!"Infected" %in% names(incidence_table)) {
    stop(
      paste0(
        "The incidence output for simulation ",
        i,
        " does not contain an 'Infected' column."
      ),
      call. = FALSE
    )
  }

  infected <- as.numeric(
    incidence_table[["Infected"]]
  )

  # Drop the day-0 initialization value.
  if (length(infected) >= ndays + 1L) {
    infected <- infected[2:(ndays + 1L)]
  } else {
    infected <- infected[-1]

    if (length(infected) < ndays) {
      infected <- c(
        infected,
        rep(0, ndays - length(infected))
      )
    }
  }

  as.numeric(infected[seq_len(ndays)])
}

# =============================================================================
# Helper: Run all simulations
# =============================================================================

run_seir_simulations <- function(
    N,
    ndays,
    ncores,
    theta,
    seeds,
    disease_name = "Disease",
    output_file_csv
) {

  cat("\n==============================================\n")
  cat(
    "Running ",
    N,
    " SEIR simulations on ",
    ncores,
    " cores\n",
    sep = ""
  )
  cat("==============================================\n")

  simulation_indices <- seq_len(N)

  if (ncores <= 1L) {

    incidence_list <- lapply(
      simulation_indices,
      run_one_seir_simulation,
      theta = theta,
      seeds = seeds,
      ndays = ndays,
      disease_name = disease_name
    )

  } else if (.Platform$OS.type == "windows") {

    # PSOCK clusters work on Windows.
    cluster <- parallel::makeCluster(ncores)

    on.exit(
      parallel::stopCluster(cluster),
      add = TRUE
    )

    parallel::clusterEvalQ(
      cluster,
      {
        library(epiworldR)
        library(data.table)
        NULL
      }
    )

    parallel::clusterExport(
      cluster,
      varlist = c("run_one_seir_simulation"),
      envir = environment()
    )

    incidence_list <- parallel::parLapply(
      cluster,
      simulation_indices,
      run_one_seir_simulation,
      theta = theta,
      seeds = seeds,
      ndays = ndays,
      disease_name = disease_name
    )

  } else {

    # Forked processing is efficient on Linux and macOS.
    incidence_list <- parallel::mclapply(
      simulation_indices,
      FUN = run_one_seir_simulation,
      theta = theta,
      seeds = seeds,
      ndays = ndays,
      disease_name = disease_name,
      mc.cores = ncores
    )
  }

  cat("\nCombining results from all simulations...\n")

  infected_matrix <- do.call(
    rbind,
    incidence_list
  )

  infected_dt <- data.table::as.data.table(
    infected_matrix
  )

  data.table::setnames(
    infected_dt,
    paste0("day_", seq_len(ndays))
  )

  cat("Saving incidence data to:\n  ", output_file_csv, "\n", sep = "")

  data.table::fwrite(
    infected_dt,
    file = output_file_csv
  )

  cat("\n==============================================\n")
  cat("Simulation complete.\n")
  cat("Total simulations: ", N, "\n", sep = "")
  cat("Days per simulation: ", ndays, "\n", sep = "")
  cat(
    "Matrix dimensions: ",
    nrow(infected_dt),
    " x ",
    ncol(infected_dt),
    "\n",
    sep = ""
  )
  cat("==============================================\n")

  cat("\nIncidence preview:\n")
  print(head(infected_dt, 10))

  infected_dt
}

# =============================================================================
# Part 2: Read parameter table and generate seeds
# =============================================================================

theta_use <- data.table::fread(
  THETA_FILE
)

N_SIMS <- nrow(theta_use)

set.seed(MODEL_SEED)

simulation_seeds <- sample.int(
  n = 1000000000L,
  size = N_SIMS,
  replace = FALSE
)

# =============================================================================
# Part 3: Run simulations and save incidence data
# =============================================================================

incidence_data <- run_seir_simulations(
  N               = N_SIMS,
  ndays           = MODEL_NDAYS,
  ncores          = N_CORES,
  theta           = theta_use,
  seeds           = simulation_seeds,
  disease_name    = "General Disease",
  output_file_csv = INCIDENCE_FILE
)

# =============================================================================
# Final checks
# =============================================================================

cat("\n==============================================\n")
cat("Final results\n")
cat("==============================================\n")

print(head(incidence_data))
print(dim(incidence_data))

cat("\nCreated files:\n")
cat("  ", THETA_FILE, "\n", sep = "")
cat("  ", INCIDENCE_FILE, "\n", sep = "")
cat("\nFinished successfully.\n")
